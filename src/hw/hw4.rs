use crate::util::read;
use std::{
    cell::RefCell,
    cmp::max,
    collections::{HashMap, HashSet, LinkedList},
    io::Error,
    rc::Rc,
};

const NODE_KEY: &str = "V";
const START_KEY: &str = "START";
const END_KEY: &str = "END";
const BASE_KEYS: [char; 5] = ['A', 'C', 'G', 'T', 'N'];

#[derive(Eq, Hash, PartialEq, Debug, Clone, Copy)]
enum Base {
    A,
    C,
    G,
    T,
    N,
}

impl Base {
    fn from_char(b: char) -> Base {
        match b {
            'A' => Base::A,
            'C' => Base::C,
            'G' => Base::G,
            'T' => Base::T,
            _ => Base::N,
        }
    }
}

#[derive(Debug)]
struct WeightedDirectedAcyclicGraph {
    nodes: HashMap<String, Rc<RefCell<Node>>>,
}

impl WeightedDirectedAcyclicGraph {
    fn new(nodes: HashMap<String, Rc<RefCell<Node>>>) -> Self {
        return WeightedDirectedAcyclicGraph { nodes };
    }

    fn get_path_scores(&self, min_score: isize, constraints: HashSet<String>) -> HashMap<String, TraceScore> {
        let mut scores: HashMap<String, TraceScore> = HashMap::new();
        let mut queue = LinkedList::new();

        self.nodes
            .iter()
            .filter(|(n, _)| constraints.is_empty() || constraints.contains(*n))
            .for_each(|(name, node)| {
                if node.borrow_mut().parents.is_empty() {
                    scores.insert(
                        name.to_string(),
                        TraceScore {
                            score: 0,
                            node_name: name.to_string(),
                            edge_name: "".to_string(),
                            parent: None,
                        },
                    );
                    queue.push_back(Rc::clone(&node));
                }
            });

        while queue.len() > 0 {
            let node = queue.pop_front().unwrap();
            let node_mut = node.borrow_mut();
            for (child_name, edge) in &node_mut.children {
                queue.push_back(Rc::clone(&edge.to));
                let parent = scores.get(&node_mut.name).unwrap();
                let score = max(min_score, edge.weight + parent.score);
                if let Some(current_weight) = scores.get(child_name) {
                    if current_weight.score >= score {
                        continue;
                    }
                }
                scores.insert(
                    child_name.to_string(),
                    TraceScore {
                        score,
                        node_name: child_name.to_string(),
                        edge_name: edge.name.to_string(),
                        parent: Some(Box::from(parent.clone())),
                    },
                );
            }
        }

        return scores;
    }

    fn print_best_path(&self) {
        let scores = self.get_path_scores(0, HashSet::new());
        let max_score = scores.values().max_by_key(|t| t.score).unwrap();
        let (begin, path) = max_score.trace_back(0, None);

        println!("Score: {}", max_score.score);
        println!("Begin: {}", begin);
        println!("End: {}", max_score.node_name);
        println!("Path: {}", path);
    }

    fn print_best_path_nodes(&self, start: String, end: String) {
        let constraints = HashSet::from([start.clone()]);
        let scores = self.get_path_scores(isize::MIN, constraints);
        let end_score = scores.get(&end).unwrap();
        let (begin, path) = end_score.trace_back(isize::MIN, Some(start));

        println!("Score: {}", end_score.score);
        println!("Begin: {}", begin);
        println!("End: {}", end_score.node_name);
        println!("Path: {}", path);
    }
}

#[derive(Debug, Clone)]
struct TraceScore {
    score: isize,
    node_name: String,
    edge_name: String,
    parent: Option<Box<TraceScore>>,
}

impl TraceScore {
    fn trace_back(&self, min_score: isize, stop_node: Option<String>) -> (String, String) {
        return self.trace_back_inner(min_score, stop_node, "".to_string());
    }

    fn trace_back_inner(&self, min_score: isize, stop_node: Option<String>, current_path: String) -> (String, String) {
        if let Some(parent) = &self.parent {
            let should_stop = stop_node.clone().map(|n| n == self.node_name).unwrap_or_default();
            if parent.score > min_score && !should_stop {
                return parent.trace_back_inner(min_score, stop_node, format!("{}{}", self.edge_name, current_path));
            } else {
                return (parent.node_name.clone(), format!("{}{}", self.edge_name, current_path));
            }
        } else {
            return (self.node_name.clone(), current_path);
        }
    }
}

#[derive(Debug)]
struct Node {
    name: String,
    parents: HashMap<String, Edge>,
    children: HashMap<String, Edge>,
}

impl Node {
    fn new(name: String) -> Self {
        return Node {
            name,
            parents: HashMap::new(),
            children: HashMap::new(),
        };
    }

    fn add_parent(&mut self, name: String, edge: Edge) {
        self.parents.insert(name, edge);
    }

    fn add_child(&mut self, name: String, edge: Edge) {
        self.children.insert(name, edge);
    }
}

#[derive(Debug, Clone)]
struct Edge {
    name: String,
    weight: isize,
    to: Rc<RefCell<Node>>,
}

pub fn run(file_path1: &str, file_path2: &str) -> Result<(), Error> {
    let (dag, start_node, end_node) = parse_dag(file_path1)?;

    println!("Part 1");
    dag.print_best_path();

    println!("\nPart 2");
    dag.print_best_path_nodes(start_node.unwrap(), end_node.unwrap());

    println!("\nPart 3");
    score_genome(file_path2)?;

    Ok(())
}

fn parse_dag(file_path: &str) -> Result<(WeightedDirectedAcyclicGraph, Option<String>, Option<String>), Error> {
    let mut start_node: Option<String> = None;
    let mut end_node: Option<String> = None;
    let mut nodes: HashMap<String, Rc<RefCell<Node>>> = HashMap::new();

    let lines = read::lines(file_path)?;
    for line in lines {
        if let Ok(ip) = line {
            let parts: Vec<&str> = ip.split_whitespace().collect();
            let graph_type = parts.get(0).unwrap();
            if *graph_type == NODE_KEY {
                let name = parts.get(1).unwrap().to_string();
                nodes.insert(name.clone(), Rc::from(RefCell::from(Node::new(name.clone()))));
                if let Some(k) = parts.get(2) {
                    match *k {
                        START_KEY => start_node = Some(name),
                        END_KEY => end_node = Some(name),
                        _ => panic!("Uknown key type '{k}' encountered"),
                    }
                }
            } else {
                // Edge
                let name = parts.get(1).unwrap().to_string();
                let from = nodes.get(&parts.get(2).unwrap().to_string()).unwrap();
                let to = nodes.get(&parts.get(3).unwrap().to_string()).unwrap();
                let weight: isize = parts.get(4).unwrap().trim().parse().unwrap();

                let edge = Edge {
                    name,
                    weight,
                    to: Rc::clone(to),
                };
                let mut to_mut = to.borrow_mut();
                let mut from_mut = from.borrow_mut();
                to_mut.add_parent(from_mut.name.clone(), edge.clone());
                from_mut.add_child(to_mut.name.clone(), edge.clone());
            }
        }
    }

    Ok((WeightedDirectedAcyclicGraph::new(nodes), start_node, end_node))
}

fn score_base(base: Base) -> f64 {
    match base {
        Base::A | Base::T => -1.49,
        Base::C | Base::G => 0.74,
        _ => 0.0,
    }
}

fn score_genome(file_path: &str) -> Result<(), Error> {
    let mut sequence = String::with_capacity(read::file_size(file_path) as usize);
    let mut header: String = String::from("NO HEADER");
    let mut base_counts: HashMap<Base, usize> = HashMap::new();
    let mut non_alpha_count: i32 = 0;
    let mut i = 0;
    let mut score: f64 = 0.0;
    let mut high_score: f64 = 0.0;
    let mut start = 0;
    let mut best_start = 0;
    let mut best_end = 0;

    let lines = read::lines(file_path)?;
    for line in lines {
        if let Ok(ip) = line {
            if ip.starts_with('>') {
                header = ip;
                continue;
            }
            for c in ip.to_uppercase().chars() {
                if BASE_KEYS.contains(&c) {
                    let base = Base::from_char(c);
                    score += score_base(base);
                    if score <= 0.0 {
                        score = 0.0;
                        start = i + 1
                    } else if score > high_score {
                        high_score = score;
                        best_start = start;
                        best_end = i + 1;
                    }
                    sequence.push(c);
                    *base_counts.entry(base).or_default() += 1;
                } else if c != ' ' {
                    // Don't count spaces for some reason?
                    non_alpha_count += 1;
                }
                i += 1;
            }
        }
    }

    println!("Fasta: {}", read::file_name_from_path(&file_path));
    println!("Non-alphabetic characters: {}", non_alpha_count);
    println!("{}", header);
    println!("*={}", base_counts.iter().fold(0, |t, (_, b)| t + b));
    for key in BASE_KEYS {
        println!("{}={}", key, base_counts.get(&Base::from_char(key)).unwrap_or(&0))
    }

    println!("\nScore: {:.2}", high_score);
    println!("Begin: {best_start}");
    println!("End: {best_end}");
    println!("Path: {}", sequence[best_start..best_end].to_string());
    println!("Description: TODO");

    Ok(())
}
