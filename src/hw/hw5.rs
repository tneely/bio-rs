use crate::util::read;
use bio::scores::blosum62;
use itertools::Itertools;
use std::cell::RefCell;
use std::cmp::max;
use std::collections::{HashMap, LinkedList};
use std::io::Error;
use std::ops::Range;
use std::rc::Rc;

const GAP_PENALTY: isize = -6;
const GAP_CHAR: char = '-';
const AA_KEYS: [char; 23] = [
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X',
];

#[derive(Debug)]
struct WeightedDirectedAcyclicGraph {
    nodes: HashMap<String, Rc<RefCell<Node>>>,
}

impl WeightedDirectedAcyclicGraph {
    fn new() -> Self {
        return WeightedDirectedAcyclicGraph { nodes: HashMap::new() };
    }

    fn add_node(&mut self, node: &Rc<RefCell<Node>>) {
        self.nodes.insert(node.borrow().name.clone(), Rc::clone(node));
    }

    fn get_node(&self, name: &str) -> Option<&Rc<RefCell<Node>>> {
        return self.nodes.get(name);
    }

    fn get_path_scores(&self, start: &str) -> HashMap<String, TraceScore> {
        let mut scores: HashMap<String, TraceScore> = HashMap::new();
        let mut queue = LinkedList::new();

        let start_node = self.nodes.get(start).unwrap();
        scores.insert(
            start.to_string(),
            TraceScore {
                score: 0,
                edge_name: "".to_string(),
                parent: None,
            },
        );
        queue.push_back(Rc::clone(start_node));

        while queue.len() > 0 {
            let node = queue.pop_front().unwrap();
            let node_mut = node.borrow_mut();
            for (child_name, edge) in &node_mut.children {
                if !scores.contains_key(child_name) {
                    queue.push_back(Rc::clone(&edge.to));
                }
                let parent = scores.get(&node_mut.name).unwrap();
                let score = max(0, edge.weight + parent.score);
                if let Some(current_weight) = scores.get(child_name) {
                    if current_weight.score >= score {
                        continue;
                    }
                }
                scores.insert(
                    child_name.to_string(),
                    TraceScore {
                        score,
                        edge_name: edge.name.to_string(),
                        parent: Some(Box::from(parent.clone())),
                    },
                );
            }
        }

        return scores;
    }

    fn print_graph(&self) {
        self.nodes.iter().for_each(|(name, _)| println!("V {name}"));
        self.nodes.iter().for_each(|(from, node)| {
            node.borrow()
                .children
                .iter()
                .for_each(|(to, e)| println!("E {} {from} {to}, {}", e.name, e.weight))
        });
    }

    fn print_edges(&self) {
        let mut edge_weights = HashMap::new();
        let mut edge_counts: HashMap<String, usize> = HashMap::new();
        for (_, node) in &self.nodes {
            for (_, edge) in &node.borrow().children {
                edge_weights.insert(edge.name.clone(), edge.weight);
                *edge_counts.entry(edge.name.clone()).or_default() += 1;
            }
        }
        println!("\nEdge Weights:");
        edge_weights.iter().sorted().for_each(|(name, weight)| println!("{name} = {weight}"));
        println!("\nEdge Counts:");
        edge_counts.iter().sorted().for_each(|(name, count)| println!("{name} = {count}"));
    }
}

#[derive(Debug, Clone)]
struct TraceScore {
    score: isize,
    edge_name: String,
    parent: Option<Box<TraceScore>>,
}

impl TraceScore {
    fn trace_back(&self) {
        println!("{}", self.trace_back_inner("".to_string()));
    }

    fn trace_back_inner(&self, current_path: String) -> String {
        return if let Some(parent) = &self.parent {
            if parent.score > 0 {
                parent.trace_back_inner(format!("{}\n{}", self.edge_name, current_path))
            } else {
                format!("{}\n{}", self.edge_name, current_path)
            }
        } else {
            current_path
        };
    }
}

#[derive(Debug)]
struct Node {
    name: String,
    children: HashMap<String, Edge>,
}

impl Node {
    fn new(name: String) -> Self {
        return Node {
            name,
            children: HashMap::new(),
        };
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

pub fn run(file_path1: &str, file_path2: &str, file_path3: &str) -> Result<(), Error> {
    let seq1 = load_sequence(file_path1)?;
    let seq2 = load_sequence(file_path2)?;
    let seq3 = load_sequence(file_path3)?;

    let dag = create_dag(seq1, seq2, seq3);
    let scores = dag.get_path_scores("(0,0,0)");
    let max_score = scores.values().max_by_key(|t| t.score).unwrap();
    println!("Score: {}", max_score.score);
    dag.print_edges();
    println!("\nLocal Alignment:");
    max_score.trace_back();

    Ok(())
}

fn load_sequence(file_path: &str) -> Result<String, Error> {
    let mut sequence = String::with_capacity(read::file_size(file_path) as usize);

    let lines = read::lines(file_path)?;
    for line in lines {
        if let Ok(ip) = line {
            if ip.starts_with('>') {
                continue;
            }
            for c in ip.to_uppercase().chars() {
                if AA_KEYS.contains(&c) {
                    sequence.push(c);
                }
            }
        }
    }

    Ok(sequence)
}

fn create_dag(seq1: String, seq2: String, seq3: String) -> WeightedDirectedAcyclicGraph {
    let mut dag = WeightedDirectedAcyclicGraph::new();

    for i in 0..seq1.len() + 1 {
        for j in 0..seq2.len() + 1 {
            for k in 0..seq3.len() + 1 {
                let name = format!("({i},{j},{k})");
                if dag.nodes.contains_key(&name) {
                    continue;
                }
                let node = Rc::from(RefCell::from(Node::new(name.clone())));
                dag.add_node(&node);
                score_edges(&node, &seq1, &seq2, &seq3, i, j, k).iter().for_each(|(p_name, edge)| {
                    if let Some(p_node) = dag.get_node(p_name) {
                        p_node.borrow_mut().add_child(name.clone(), edge.clone());
                    } else {
                        let p_node = Rc::from(RefCell::from(Node::new(p_name.clone())));
                        dag.add_node(&p_node);
                        p_node.borrow_mut().add_child(name.clone(), edge.clone());
                    }
                });
            }
        }
    }

    return dag;
}

fn get_range(i: usize) -> Range<usize> {
    return if i == 0 { 0..i + 1 } else { i - 1..i + 1 };
}

fn score_edges(node: &Rc<RefCell<Node>>, seq1: &str, seq2: &str, seq3: &str, i: usize, j: usize, k: usize) -> Vec<(String, Edge)> {
    let mut edges: Vec<(String, Edge)> = Vec::new();

    for i2 in get_range(i) {
        for j2 in get_range(j) {
            for k2 in get_range(k) {
                if i2 == i && j2 == j && k2 == k {
                    continue;
                }
                let node_name = format!("({i2},{j2},{k2})");
                let r1 = if i2 == i { GAP_CHAR } else { seq1.as_bytes()[i2] as char };
                let r2 = if j2 == j { GAP_CHAR } else { seq2.as_bytes()[j2] as char };
                let r3 = if k2 == k { GAP_CHAR } else { seq3.as_bytes()[k2] as char };
                let edge_name = format!("{r1}{r2}{r3}");
                let score = score_edge(&edge_name);
                edges.push((
                    node_name,
                    Edge {
                        name: edge_name,
                        weight: score,
                        to: Rc::clone(node),
                    },
                ));
            }
        }
    }

    return edges;
}

fn score_edge(edge_name: &str) -> isize {
    let mut score: isize = 0;

    let i = edge_name.as_bytes()[0];
    let j = edge_name.as_bytes()[1];
    let k = edge_name.as_bytes()[2];

    score += score_pair(i, j);
    score += score_pair(i, k);
    score += score_pair(j, k);

    return score;
}

fn score_pair(a: u8, b: u8) -> isize {
    return if a == GAP_CHAR as u8 && b == GAP_CHAR as u8 {
        0
    } else if a == GAP_CHAR as u8 || b == GAP_CHAR as u8 {
        GAP_PENALTY
    } else {
        blosum62(a, b) as isize
    };
}
