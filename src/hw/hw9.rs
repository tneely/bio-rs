use itertools::Itertools;

use crate::util::read;
use std::{collections::HashMap, io::Error};

const NEUTRAL_FREQ_FILE: &str = "/Users/tneely/dev/bio-rs/data/hw/hw9/STATE1_anc_rep_counts.txt";
const CONSERVED_FREQ_FILE: &str = "/Users/tneely/dev/bio-rs/data/hw/hw9/STATE2_codon1_2_counts.txt";

pub fn run(file_path: &str) -> Result<(), Error> {
    let neutral_state = load_freqs(NEUTRAL_FREQ_FILE)?;
    let conserved_state = load_freqs(CONSERVED_FREQ_FILE)?;
    let states = [&neutral_state, &conserved_state];
    let start_probs = [0.95_f64.ln(), 0.05_f64.ln()];
    let transition_probs = [[0.95_f64.ln(), 0.05_f64.ln()], [0.10_f64.ln(), 0.90_f64.ln()]];

    let mut state_segments: [Vec<(usize, usize)>; 2] = [Vec::new(), Vec::new()];
    let mut lines = read::lines(file_path)?;
    while let Some(line) = lines.next() {
        if let Ok(ip) = line {
            if ip.starts_with("#") {
                let (start, _) = parse_range(&ip);
                let aln = parse_sequence(&lines.next().unwrap()?, &lines.next().unwrap()?, &lines.next().unwrap()?);
                let mut current_state: Option<usize> = None;
                let (mut seg_start, mut seg_end) = (0, 0);
                for (i, em) in aln.iter().enumerate() {
                    // TODO: Compute forward (see hw8), then backtrack to find highest probability path?
                    if let Some(cs) = current_state {
                        let p1 = transition_probs[cs][0] + states[0].get_prob(em);
                        let p2 = transition_probs[cs][1] + states[1].get_prob(em);
                        let next_state = if p1 > p2 { 0 } else { 1 };
                        if cs == next_state {
                            seg_end += 1
                        } else {
                            state_segments[cs].push((seg_start, seg_end));
                            (seg_start, seg_end) = (start + i, start + i);
                        }
                    } else {
                        let p1 = start_probs[0] + states[0].get_prob(em);
                        let p2 = start_probs[1] + states[1].get_prob(em);
                        current_state = if p1 > p2 { Some(0) } else { Some(1) };
                        (seg_start, seg_end) = (start + i, start + i);
                    }
                }
            }
        }
    }

    println!("\nState Histogram:");
    println!("TODO");

    println!("\nSegment Histogram:");
    println!("TODO");

    println!("\nInitial State Probabilities:");
    (0..states.len()).for_each(|i| println!("{}={:.5}", i + 1, start_probs[i].exp()));

    println!("\nTransition Probabilities:");
    (0..states.len()).for_each(|i| (0..states.len()).for_each(|j| println!("{},{}={:.5}", i + 1, j + 1, transition_probs[i][j].exp())));

    println!("\nEmission Probabilities:");
    // neutral_state.print_probs("1");
    // conserved_state.print_probs("2");

    println!("\nLongest Segment List:");
    println!("{state_segments:?}");
    let mut top_conserved: Vec<&(usize, usize)> = state_segments[1].iter().sorted_by_key(|(s, e)| s - e).collect();
    top_conserved.truncate(10);
    top_conserved.iter().for_each(|(s, e)| println!("{s} {e}"));

    println!("\nAnnotations:");
    println!("TODO");

    Ok(())
}

struct StateProbability {
    emission_probabilities: HashMap<String, f64>,
}

impl StateProbability {
    fn from_counts(emission_counts: HashMap<String, usize>) -> Self {
        let mut emission_probabilities = HashMap::new();
        let total_count = emission_counts.iter().fold(0.0, |prev, (_, next)| prev + *next as f64);
        for (em, c) in emission_counts {
            emission_probabilities.insert(em, c as f64 / total_count);
        }
        return StateProbability { emission_probabilities };
    }

    fn get_prob(&self, em: &String) -> f64 {
        return self.emission_probabilities[em];
    }

    fn print_probs(&self, id: &str) {
        for (em, c) in self.emission_probabilities.iter().sorted_by_key(|(em, _)| *em) {
            println!("{id},{em}={c:.5}")
        }
    }
}

fn load_freqs(file_path: &str) -> Result<StateProbability, Error> {
    let mut emission_counts = HashMap::new();
    let lines = read::lines(file_path)?;
    for line in lines {
        if let Ok(ip) = line {
            let mut line_parts = ip.split("\t");
            let emission = line_parts.next().unwrap();
            let count = line_parts.next().unwrap().parse::<usize>().unwrap();
            emission_counts.insert(emission.to_string(), count);
        }
    }

    let state_prob = StateProbability::from_counts(emission_counts);

    Ok(state_prob)
}

// Expected format: # chrX:152767491-152767698
fn parse_range(string: &str) -> (usize, usize) {
    let range_str = string.split(":").nth(1).unwrap();
    let mut range_parts = range_str.split("-");
    let start: usize = range_parts.next().unwrap().parse().unwrap();
    let end: usize = range_parts.next().unwrap().parse().unwrap();

    return (start, end);
}

// Expected format: hg18\tATAAAA
fn parse_sequence(l1: &str, l2: &str, l3: &str) -> Vec<String> {
    let seq1 = l1.split("\t").nth(1).unwrap();
    let seq2 = l2.split("\t").nth(1).unwrap();
    let seq3 = l3.split("\t").nth(1).unwrap();

    let mut aln = Vec::new();
    for i in 0..seq1.len() {
        aln.push(format!(
            "{}{}{}",
            seq1.as_bytes()[i] as char,
            seq2.as_bytes()[i] as char,
            seq3.as_bytes()[i] as char
        ))
    }

    return aln;
}
