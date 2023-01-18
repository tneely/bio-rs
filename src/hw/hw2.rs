use crate::util::read;
use rand::Rng;
use std::collections::HashMap;
use std::fs::File;
use std::io::{Error, Write};
use std::path::Path;

const BASE_KEYS: [char; 5] = ['A', 'C', 'G', 'T', 'N'];
const DATA_PATH: &str = "./data/hw/hw2/";

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

struct FrequencyDistribution {
    counts: HashMap<Base, HashMap<Option<Base>, usize>>,
    freqs: HashMap<Base, HashMap<Option<Base>, f64>>,
    cond_freqs: HashMap<Base, HashMap<Option<Base>, f64>>,
    base_count: usize,
}

impl FrequencyDistribution {
    fn new(counts: HashMap<Base, HashMap<Option<Base>, usize>>) -> FrequencyDistribution {
        let base_count = counts.iter().fold(0, |t1, (_, pair)| t1 + pair.get(&None).unwrap_or(&0));
        let pair_count = counts.iter().fold(0, |t1, (_, pair)| t1 + pair.iter().fold(0, |t2, (_, count)| t2 + count)) - base_count;

        let mut freqs: HashMap<Base, HashMap<Option<Base>, f64>> = HashMap::new();
        counts.iter().for_each(|(base1, pair)| {
            pair.iter().for_each(|(base2, count)| {
                let total_count = if let Some(_) = base2 { pair_count } else { base_count };
                let freq = *count as f64 / total_count as f64;
                freqs.entry(*base1).or_default().insert(*base2, freq);
            });
        });

        let mut cond_freqs: HashMap<Base, HashMap<Option<Base>, f64>> = HashMap::new();
        for base1 in [Base::A, Base::C, Base::G, Base::T] {
            let base_freq = freqs
                .get(&base1)
                .unwrap_or(&HashMap::new())
                .iter()
                .fold(0.0, |t, (pair, freq)| if let Some(_) = pair { t + freq } else { t });
            for base2 in [Base::A, Base::C, Base::G, Base::T] {
                let pair_freq = freqs.get(&base1).unwrap_or(&HashMap::new()).get(&Some(base2)).unwrap_or(&0.0) / base_freq;
                cond_freqs.entry(base1).or_default().insert(Some(base2), pair_freq);
            }
        }

        return FrequencyDistribution {
            counts,
            freqs,
            cond_freqs,
            base_count,
        };
    }

    fn predict_next_base(&self, prev_base: Option<Base>) -> Base {
        let mut rng = rand::thread_rng();
        let mut score: f64 = rng.gen();

        if let Some(base1) = prev_base {
            for base2 in [Base::A, Base::C, Base::G, Base::T] {
                score -= self.get_conditional_freq(base1, base2);
                if score < 0.0 {
                    return base2;
                }
            }
        } else {
            for base in [Base::A, Base::C, Base::G, Base::T] {
                score -= self.get_base_freq(base);
                if score < 0.0 {
                    return base;
                }
            }
        }

        panic!("It should be impossible to get here!");
    }

    fn get_base_count(&self, base: Base) -> usize {
        return *self.counts.get(&base).unwrap_or(&HashMap::new()).get(&None).unwrap_or(&0);
    }

    fn get_base_freq(&self, base: Base) -> f64 {
        return *self.freqs.get(&base).unwrap_or(&HashMap::new()).get(&None).unwrap_or(&0.0);
    }

    fn get_pair_count(&self, prev_base: Base, curr_base: Base) -> usize {
        return *self.counts.get(&prev_base).unwrap_or(&HashMap::new()).get(&Some(curr_base)).unwrap_or(&0);
    }

    fn get_pair_freq(&self, prev_base: Base, curr_base: Base) -> f64 {
        return *self.freqs.get(&prev_base).unwrap_or(&HashMap::new()).get(&Some(curr_base)).unwrap_or(&0.0);
    }

    fn get_conditional_freq(&self, prev_base: Base, curr_base: Base) -> f64 {
        return *self.cond_freqs.get(&prev_base).unwrap_or(&HashMap::new()).get(&Some(curr_base)).unwrap_or(&0.0);
    }

    fn print_base_count(&self) {
        println!("{}={}", "*", self.base_count);
        for base in [Base::A, Base::C, Base::G, Base::T, Base::N] {
            println!("{:?}={}", base, self.get_base_count(base));
        }
    }

    fn print_base_freq(&self) {
        println!("\nNucleotide Frequencies:");
        for base in [Base::A, Base::C, Base::G, Base::T] {
            println!("{:?}={:.4}", base, self.get_base_freq(base));
        }
    }

    fn print_pair_count(&self) {
        println!("\nDinucleotide Count Matrix:");
        for base1 in [Base::A, Base::C, Base::G, Base::T] {
            print!("{:?}=", base1);
            for base2 in [Base::A, Base::C, Base::G, Base::T] {
                print!("{} ", self.get_pair_count(base1, base2));
            }
            println!();
        }
    }

    fn print_pair_freq(&self) {
        println!("\nDinucleotide Frequency Matrix:");
        for base1 in [Base::A, Base::C, Base::G, Base::T] {
            print!("{:?}=", base1);
            for base2 in [Base::A, Base::C, Base::G, Base::T] {
                print!("{:.4} ", self.get_pair_freq(base1, base2));
            }
            println!();
        }
    }

    fn print_pair_conditional_freq(&self) {
        println!("\nConditional Frequency Matrix:");
        for base1 in [Base::A, Base::C, Base::G, Base::T] {
            print!("{:?}=", base1);
            for base2 in [Base::A, Base::C, Base::G, Base::T] {
                print!("{:.4} ", self.get_conditional_freq(base1, base2));
            }
            println!();
        }
    }
}

pub fn run(file_path1: &str) -> Result<(), Error> {
    println!("Fasta 1: {}", file_name_from_path(file_path1));

    let file_freq_dist = count_bases(file_path1)?;

    let equal_base_pairs = [Base::A, Base::C, Base::G, Base::T].map(|b| (b, HashMap::from([(None, 1)])));
    let equal_freq_dist = FrequencyDistribution::new(HashMap::from(equal_base_pairs));

    let file_path2 = format!("{}{}", DATA_PATH, "simulated_equal_freq.fa");
    println!("\nFasta 2: {}", file_name_from_path(file_path2.as_str()));
    gen_sequence(file_path2.as_str(), &equal_freq_dist, file_freq_dist.base_count, false)?;
    count_bases(file_path2.as_str())?;

    let file_path3 = format!("{}{}", DATA_PATH, "simulated_markov_0.fa");
    println!("\nFasta 3: {}", file_name_from_path(file_path3.as_str()));
    gen_sequence(file_path3.as_str(), &file_freq_dist, file_freq_dist.base_count, false)?;
    count_bases(file_path3.as_str())?;

    let file_path4 = format!("{}{}", DATA_PATH, "simulated_markov_1.fa");
    println!("\nFasta 4: {}", file_name_from_path(file_path4.as_str()));
    gen_sequence(file_path4.as_str(), &file_freq_dist, file_freq_dist.base_count, true)?;
    count_bases(file_path4.as_str())?;

    Ok(())
}

fn file_name_from_path(file_path: &str) -> &str {
    return Path::new(file_path).file_name().unwrap().to_str().unwrap();
}

fn count_bases(file_path: &str) -> Result<FrequencyDistribution, Error> {
    let mut header: String = String::from("NO HEADER");
    let mut base_counts: HashMap<Base, HashMap<Option<Base>, usize>> = HashMap::new();
    let mut non_alpha_count: i32 = 0;
    let mut prev_base: Option<Base> = None;

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
                    *base_counts.entry(base).or_default().entry(None).or_default() += 1;
                    if let Some(prev) = prev_base {
                        *base_counts.entry(prev).or_default().entry(Some(base)).or_default() += 1;
                    }
                    prev_base = Some(base);
                } else if c != ' ' {
                    // Don't count spaces for some reason?
                    non_alpha_count += 1;
                }
            }
        }
    }

    println!("Non-alphabetic characters: {}", non_alpha_count);
    println!("{}", header);

    let freq_dist = FrequencyDistribution::new(base_counts);

    freq_dist.print_base_count();
    freq_dist.print_base_freq();
    freq_dist.print_pair_count();
    freq_dist.print_pair_freq();
    freq_dist.print_pair_conditional_freq();

    Ok(freq_dist)
}

fn gen_sequence(file_path: &str, freq_dist: &FrequencyDistribution, len: usize, use_prev: bool) -> Result<(), Error> {
    let mut file = File::create(file_path)?;
    let mut prev_base: Option<Base> = None;

    for _ in 0..len {
        let base = freq_dist.predict_next_base(prev_base);
        file.write(format!("{:?}", base).as_ref())?;
        if use_prev {
            prev_base = Some(base);
        }
    }

    Ok(())
}
