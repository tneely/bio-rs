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

    let (aln, start, end) = load_alignment(file_path)?;

    // Compute scores
    let mut score_tracker = Vec::from([vec![0.0; aln.len()], vec![0.0; aln.len()]]);
    let mut state_tracker = Vec::from([vec![0; aln.len()], vec![0_usize; aln.len()]]);
    for (i, em) in aln.iter().enumerate() {
        if i == 0 {
            for s in 0..states.len() {
                score_tracker[s][i] = start_probs[s] + states[s].get_prob(em)
            }
        } else {
            for s in 0..states.len() {
                let p0 = score_tracker[0][i - 1] + transition_probs[0][s] + states[s].get_prob(em);
                let p1 = score_tracker[1][i - 1] + transition_probs[1][s] + states[s].get_prob(em);
                if p0 > p1 {
                    score_tracker[s][i] = p0;
                    state_tracker[s][i] = 0;
                } else {
                    score_tracker[s][i] = p1;
                    state_tracker[s][i] = 1;
                }
            }
        }
    }

    // Backtrack on best path
    let mut state_segments: [Vec<(isize, isize)>; 2] = [Vec::new(), Vec::new()];
    let mut current_state: Option<usize> = None;
    let mut best_path = vec![0; aln.len()];
    let mut seg_end = 0_isize;
    for i in (0..aln.len()).rev() {
        if let Some(cs) = current_state {
            best_path[i] = cs;
            let next_state = state_tracker[cs][i];
            if cs != next_state || i == 0 {
                let current_pos = start + i as isize;
                state_segments[cs].push((current_pos, seg_end));
                seg_end = current_pos - 1;
            }
            current_state = Some(next_state)
        } else {
            let p0 = score_tracker[0][i];
            let p1 = score_tracker[1][i];
            current_state = if p0 > p1 { Some(0) } else { Some(1) };
            seg_end = end;
        }
    }

    println!("\nState Histogram:");
    (0..states.len()).for_each(|i| println!("{}={:.5}", i + 1, state_segments[i].iter().fold(0, |t, (s, e)| t + e - s)));

    println!("\nSegment Histogram:");
    (0..states.len()).for_each(|i| println!("{}={:.5}", i + 1, state_segments[i].len()));

    println!("\nInitial State Probabilities:");
    (0..states.len()).for_each(|i| println!("{}={:.5}", i + 1, start_probs[i].exp()));

    println!("\nTransition Probabilities:");
    (0..states.len()).for_each(|i| (0..states.len()).for_each(|j| println!("{},{}={:.5}", i + 1, j + 1, transition_probs[i][j].exp())));

    println!("\nEmission Probabilities:");
    neutral_state.print_probs("1");
    conserved_state.print_probs("2");

    println!("\nLongest Segment List:");
    let mut top_conserved: Vec<&(isize, isize)> = state_segments[1].iter().sorted_by_key(|(s, e)| s - e).collect();
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
            emission_probabilities.insert(em, (c as f64 / total_count).ln());
        }
        return StateProbability { emission_probabilities };
    }

    fn get_prob(&self, em: &String) -> f64 {
        return self.emission_probabilities[em];
    }

    fn print_probs(&self, id: &str) {
        for (em, c) in self.emission_probabilities.iter().sorted_by_key(|(em, _)| *em) {
            println!("{id},{em}={:.5}", c.exp())
        }
    }
}

fn load_alignment(file_path: &str) -> Result<(Vec<String>, isize, isize), Error> {
    let mut lines = read::lines(file_path)?;
    let mut aln = Vec::new();
    let mut start: Option<isize> = None;
    let mut end: Option<isize> = None;
    while let Some(line) = lines.next() {
        if let Ok(ip) = line {
            if ip.starts_with("#") {
                let (s, e) = parse_range(&ip);
                if start == None {
                    start = Some(s);
                } else {
                    end = Some(e);
                }

                // Expected format: hg18\tATAAAA
                let l1 = lines.next().unwrap()?;
                let seq1 = l1.split("\t").nth(1).unwrap();
                let l2 = lines.next().unwrap()?;
                let seq2 = l2.split("\t").nth(1).unwrap();
                let l3 = lines.next().unwrap()?;
                let seq3 = l3.split("\t").nth(1).unwrap();

                for i in 0..seq1.len() {
                    aln.push(format!(
                        "{}{}{}",
                        seq1.as_bytes()[i] as char,
                        seq2.as_bytes()[i] as char,
                        seq3.as_bytes()[i] as char
                    ))
                }
            }
        }
    }

    return Ok((aln, start.unwrap(), end.unwrap()));
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
fn parse_range(string: &str) -> (isize, isize) {
    let range_str = string.split(":").nth(1).unwrap();
    let mut range_parts = range_str.split("-");
    let start: isize = range_parts.next().unwrap().parse().unwrap();
    let end: isize = range_parts.next().unwrap().parse().unwrap();

    return (start, end);
}
