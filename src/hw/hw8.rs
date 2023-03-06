use itertools::Itertools;

use crate::util::read;
use std::{collections::HashMap, io::Error};

pub fn run(file_path: &str) -> Result<(), Error> {
    let mut start_probs = [0.996_f64.ln(), 0.004_f64.ln()];
    let mut transition_probs = [[0.999_f64.ln(), 0.001_f64.ln()], [0.01_f64.ln(), 0.99_f64.ln()]];
    let emissions_1 = HashMap::from([('A', 0.3_f64.ln()), ('T', 0.3_f64.ln()), ('G', 0.2_f64.ln()), ('C', 0.2_f64.ln())]);
    let emissions_2 = HashMap::from([('A', 0.15_f64.ln()), ('T', 0.15_f64.ln()), ('G', 0.35_f64.ln()), ('C', 0.35_f64.ln())]);
    let mut emission_probs = [emissions_1, emissions_2];

    let sequence = load_fasta(file_path)?;
    let n = sequence.len();

    let mut iterations = 0;
    let mut prev_ll: f64 = -1.0;
    let mut current_ll: f64 = 0.0;
    while (prev_ll - current_ll).abs() > 0.1 {
        println!(
            "Diff: {}\nLL: {current_ll}\nSP: {start_probs:?}\nTP: {transition_probs:?}\nEP: {emission_probs:?}",
            prev_ll - current_ll
        );

        iterations += 1;
        prev_ll = current_ll;

        // Compute forward and backward probabilities
        let forward = compute_scores(sequence.chars().into_iter(), &emission_probs, &start_probs, &transition_probs);
        // current_ll = -forward.iter().fold(0.0, |t, c| t + sum_log_prob(c[0], c[1]));
        let last = forward.last().unwrap();
        current_ll = sum_log_prob(last[0], last[1]);
        let backward = compute_scores(sequence.chars().rev(), &emission_probs, &start_probs, &transition_probs);

        let mut fb = Vec::from([0.0, 0.0]);
        let mut fwb = Vec::from([0.0, 0.0]);
        for t in 0..n {
            for s_i in 0..transition_probs.len() {
                fb[s_i] = if t == 0 {
                    forward[t][s_i] + backward[n - t - 1][s_i] - current_ll
                } else {
                    sum_log_prob(fb[s_i], forward[t][s_i] + backward[n - t - 1][s_i] - current_ll)
                };
                if t < n - 1 {
                    for s_j in 0..transition_probs.len() {
                        fwb[s_i] = if t == 0 {
                            forward[t][s_i] + backward[n - t - 2][s_j] + transition_probs[s_i][s_j] - current_ll
                        } else {
                            sum_log_prob(fwb[s_i], forward[t][s_i] + backward[n - t - 2][s_j] + transition_probs[s_i][s_j] - current_ll)
                        };
                    }
                }
            }
        }

        // New start probs
        start_probs[0] = forward[0][0] + backward[n - 2][0] - current_ll;
        start_probs[1] = forward[0][1] + backward[n - 2][1] - current_ll;

        // Zero emission probs
        for s_i in 0..emission_probs.len() {
            *emission_probs[s_i].get_mut(&'A').unwrap() = 0.0;
            *emission_probs[s_i].get_mut(&'C').unwrap() = 0.0;
            *emission_probs[s_i].get_mut(&'G').unwrap() = 0.0;
            *emission_probs[s_i].get_mut(&'T').unwrap() = 0.0;
        }

        let mut new_t = [[0.0, 0.0], [0.0, 0.0]];
        for t in 0..n {
            for s_i in 0..emission_probs.len() {
                let res = sequence.as_bytes()[t] as char;
                let old = emission_probs[s_i][&res];
                *emission_probs[s_i].get_mut(&res).unwrap() = sum_log_prob(old, forward[t][s_i] + backward[n - t - 1][s_i]);
                if t < n - 1 {
                    for s_j in 0..transition_probs.len() {
                        new_t[s_i][s_j] = sum_log_prob(new_t[s_i][s_j], forward[t][s_j] + backward[n - t - 2][s_i] + transition_probs[s_i][s_j]);
                    }
                }
            }
        }
        transition_probs = new_t;

        // Divide transition probs
        for s_i in 0..transition_probs.len() {
            for s_j in 0..transition_probs.len() {
                transition_probs[s_i][s_j] -= fwb[s_i] - current_ll;
                // transition_probs[s_i][s_j] -= current_ll;
            }
        }

        // Divide emission probs
        for s_i in 0..emission_probs.len() {
            *emission_probs[s_i].get_mut(&'A').unwrap() -= fb[s_i] - current_ll;
            *emission_probs[s_i].get_mut(&'C').unwrap() -= fb[s_i] - current_ll;
            *emission_probs[s_i].get_mut(&'G').unwrap() -= fb[s_i] - current_ll;
            *emission_probs[s_i].get_mut(&'T').unwrap() -= fb[s_i] - current_ll;

            // *emission_probs[s_i].get_mut(&'A').unwrap() -= current_ll;
            // *emission_probs[s_i].get_mut(&'C').unwrap() -= current_ll;
            // *emission_probs[s_i].get_mut(&'G').unwrap() -= current_ll;
            // *emission_probs[s_i].get_mut(&'T').unwrap() -= current_ll;
        }
    }

    println!("\nIterations for Convergence:\n{iterations}");

    println!("\nLog Likelihood:\n{:.3}", current_ll);

    println!("\nInitial State Probabilities:");
    (0..start_probs.len()).for_each(|i| println!("{i}={:.3e}", start_probs[i].exp()));

    println!("\nTransition Probabilities:");
    (0..transition_probs.len()).for_each(|i| (0..transition_probs.len()).for_each(|j| println!("{},{}={:.3e}", i + 1, j + 1, transition_probs[i][j].exp())));

    println!("\nEmission Probabilities:");
    (0..emission_probs.len()).for_each(|i| {
        emission_probs[i]
            .iter()
            .sorted_by_key(|(c, _)| *c)
            .for_each(|(c, p)| println!("{},{c}={:.3e}", i + 1, p.exp()))
    });

    Ok(())
}

fn load_fasta(file_path: &str) -> Result<String, Error> {
    let mut sequence = String::with_capacity(read::file_size(file_path) as usize);
    let lines = read::lines(file_path)?;
    for line in lines {
        if let Ok(ip) = line {
            if ip.starts_with('>') {
                println!("Fasta: {}", read::file_name_from_path(&file_path));
                println!("{}", ip);
                continue;
            }
            for c in ip.to_uppercase().chars() {
                sequence.push(c);
            }
        }
    }

    Ok(sequence)
}

fn compute_scores<I>(mut char_iter: I, emission_probs: &[HashMap<char, f64>; 2], start_probs: &[f64; 2], transition_probs: &[[f64; 2]; 2]) -> Vec<Vec<f64>>
where
    I: Iterator<Item = char>,
{
    let mut scores: Vec<Vec<f64>> = Vec::new();
    let first_c = char_iter.next().unwrap();
    let values: Vec<f64> = (0..emission_probs.len()).map(|i| (start_probs[i] + emission_probs[i][&first_c])).collect();
    scores.push(values);

    let mut prev = 0;
    for c in char_iter {
        let values = (0..emission_probs.len())
            .map(|i| {
                let first = scores[prev][0] + emission_probs[i][&c] + transition_probs[0][i];
                sum_log_prob(first, scores[prev][1] + emission_probs[i][&c] + transition_probs[1][i])
            })
            .collect();
        scores.push(values);
        prev += 1;
    }

    return scores;
}

fn sum_log_prob(a: f64, b: f64) -> f64 {
    return if a > b { a + logp1exp(b - a) } else { b + logp1exp(a - b) };
}

fn logp1exp(x: f64) -> f64 {
    return if x < -709.089565713 { 0.0 } else { x.exp().ln_1p() };
}
