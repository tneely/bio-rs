use crate::util::read;
use std::{collections::HashMap, io::Error};

pub fn run(file_path: &str) -> Result<(), Error> {
    let mut start_probs = [0.996, 0.004];
    let mut transition_probs = [[0.999, 0.001], [0.01, 0.99]];
    let mut state_1 = HashMap::from([('A', 0.3), ('T', 0.3), ('G', 0.2), ('C', 0.2)]);
    let mut state_2 = HashMap::from([('A', 0.15), ('T', 0.15), ('G', 0.35), ('C', 0.35)]);
    let states = [state_1, state_2];

    let sequence = load_fasta(file_path)?;

    let mut iterations = 0;
    let mut prev_ll: f64 = -1.0;
    let mut current_ll: f64 = 0.0;
    while (prev_ll - current_ll).abs() > 0.1 {
        iterations += 1;
        prev_ll = current_ll;
        current_ll = 0.0;

        // Compute forward and reverse probabilities
        let forward = compute_scores(sequence.chars().into_iter(), &states, &start_probs, &transition_probs);
        let reverse = compute_scores(sequence.chars().rev(), &states, &start_probs, &transition_probs);

        // TODO: Everything else
    }

    println!("\nIterations for Convergence:\n{iterations}");

    println!("\nLog Likelihood:\n{:.3}", current_ll);

    println!("\nInitial State Probabilities:");
    (0..states.len()).for_each(|i| println!("{i}={:.3e}", start_probs[i]));

    println!("\nTransition Probabilities:");
    (0..states.len()).for_each(|i| (0..states.len()).for_each(|j| println!("{},{}={:.3e}", i + 1, j + 1, transition_probs[i][j])));

    println!("\nEmission Probabilities:");
    (0..states.len()).for_each(|i| states[i].iter().for_each(|(c, p)| println!("{},{c}={p:.3e}", i + 1)));

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

fn compute_scores<I>(mut char_iter: I, states: &[HashMap<char, f64>; 2], start_probs: &[f64; 2], transition_probs: &[[f64; 2]; 2]) -> Vec<Vec<f64>>
where
    I: Iterator<Item = char>,
{
    let mut scores: Vec<Vec<f64>> = Vec::new();
    let first_c = char_iter.next().unwrap();
    let values: Vec<f64> = (0..states.len()).map(|i| start_probs[i] * states[i][&first_c]).collect();
    scores.push(values);

    let mut prev = 0;
    for c in char_iter {
        let values = (0..states.len())
            .map(|i| (0..states.len()).fold(0.0, |t, j| t + scores[prev][j] * states[i][&c] * transition_probs[j][i]))
            .collect();
        scores.push(values);
        prev += 1;
    }

    return scores;
}
