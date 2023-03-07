use itertools::Itertools;

use crate::util::read;
use std::{collections::HashMap, io::Error};

pub fn run(file_path: &str) -> Result<(), Error> {
    let states = 2;
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
        iterations += 1;
        prev_ll = current_ll;

        // Compute forward and backward probabilities
        let forward = compute_forward_scores(&sequence, &emission_probs, &start_probs, &transition_probs);
        let last = forward.last().unwrap();
        current_ll = sum_log_prob(last[0], last[1]);
        let backward = compute_backward_scores(&sequence, &emission_probs, &transition_probs);

        // New start probs
        for s_i in 0..states {
            start_probs[s_i] = forward[0][s_i] + backward[0][s_i] - current_ll;
        }

        // Calculate gamma and xi for 1..T-1
        let mut gamma = Vec::from([0.0, 0.0]);
        let mut xi = Vec::from([Vec::from([0.0, 0.0]), Vec::from([0.0, 0.0])]);
        for t in 0..n - 1 {
            for s_i in 0..states {
                gamma[s_i] = if t == 0 {
                    forward[t][s_i] + backward[t][s_i] - current_ll
                } else {
                    sum_log_prob(gamma[s_i], forward[t][s_i] + backward[t][s_i] - current_ll)
                };
                let res = sequence.as_bytes()[t + 1] as char;
                for s_j in 0..states {
                    xi[s_i][s_j] = if t == 0 {
                        forward[t][s_i] + backward[t + 1][s_j] + transition_probs[s_i][s_j] + emission_probs[s_j][&res] - current_ll
                    } else {
                        sum_log_prob(
                            xi[s_i][s_j],
                            forward[t][s_i] + backward[t + 1][s_j] + transition_probs[s_i][s_j] + emission_probs[s_j][&res] - current_ll,
                        )
                    };
                }
            }
        }

        // Calculate transition probs
        for s_i in 0..states {
            for s_j in 0..states {
                transition_probs[s_i][s_j] = xi[s_i][s_j] - gamma[s_i];
            }
        }

        // Add in last element to gamma for 1..T
        for s_i in 0..states {
            gamma[s_i] = sum_log_prob(gamma[s_i], forward[n - 1][s_i] + backward[n - 1][s_i] - current_ll)
        }

        // Calculate emission probs
        for s_i in 0..states {
            for t in 0..n {
                let res = sequence.as_bytes()[t] as char;
                *emission_probs[s_i].get_mut(&res).unwrap() = if t == 0 {
                    forward[t][s_i] + backward[t][s_i] - current_ll
                } else {
                    let cum = emission_probs[s_i][&res];
                    sum_log_prob(cum, forward[t][s_i] + backward[t][s_i] - current_ll)
                }
            }

            *emission_probs[s_i].get_mut(&'A').unwrap() -= gamma[s_i];
            *emission_probs[s_i].get_mut(&'C').unwrap() -= gamma[s_i];
            *emission_probs[s_i].get_mut(&'G').unwrap() -= gamma[s_i];
            *emission_probs[s_i].get_mut(&'T').unwrap() -= gamma[s_i];
        }
    }

    println!("\nIterations for Convergence:\n{iterations}");

    println!("\nLog Likelihood:\n{:.3}", current_ll);

    println!("\nInitial State Probabilities:");
    (0..states).for_each(|i| println!("{}={:.3e}", i + 1, start_probs[i].exp()));

    println!("\nTransition Probabilities:");
    (0..states).for_each(|i| (0..states).for_each(|j| println!("{},{}={:.3e}", i + 1, j + 1, transition_probs[i][j].exp())));

    println!("\nEmission Probabilities:");
    (0..states).for_each(|i| {
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

fn compute_forward_scores(seq: &str, emission_probs: &[HashMap<char, f64>; 2], start_probs: &[f64; 2], transition_probs: &[[f64; 2]; 2]) -> Vec<Vec<f64>> {
    let states = start_probs.len();
    let mut scores: Vec<Vec<f64>> = Vec::new();
    let mut char_iter = seq.chars().into_iter();
    let first_c = char_iter.next().unwrap();
    let values: Vec<f64> = (0..states).map(|i| (start_probs[i] + emission_probs[i][&first_c])).collect();
    scores.push(values);

    let mut prev = 0;
    for c in char_iter {
        let values = (0..states)
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

fn compute_backward_scores(seq: &str, emission_probs: &[HashMap<char, f64>; 2], transition_probs: &[[f64; 2]; 2]) -> Vec<Vec<f64>> {
    let states = emission_probs.len();
    let mut scores: Vec<Vec<f64>> = vec![Vec::from([0.0, 0.0]); seq.len()];
    let pos_iter = (0..seq.len() - 1).rev();

    for t in pos_iter {
        let res = seq.as_bytes()[t + 1] as char;
        let values = (0..states)
            .map(|i| {
                let first = scores[t + 1][0] + emission_probs[0][&res] + transition_probs[i][0];
                sum_log_prob(first, scores[t + 1][1] + emission_probs[1][&res] + transition_probs[i][1])
            })
            .collect();
        scores[t] = values;
    }

    return scores;
}

fn sum_log_prob(a: f64, b: f64) -> f64 {
    return if a > b { a + logp1exp(b - a) } else { b + logp1exp(a - b) };
}

fn logp1exp(x: f64) -> f64 {
    return if x < -709.089565713 { 0.0 } else { x.exp().ln_1p() };
}
