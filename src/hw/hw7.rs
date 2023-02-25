use crate::util::read;
use itertools::{Itertools, Position};
use rand::Rng;
use std::{collections::HashMap, io::Error};

const D_SCORE_1: f64 = -20.0;
const D_SCORE_2: f64 = -5.0;
const BACKGROUND_N: f64 = 8_422_401.0;

pub fn run(file_path: &str) -> Result<(), Error> {
    let default_scores = HashMap::from([(0, -0.1077), (1, 0.47720), (2, 1.0622), (3, 1.6748)]);
    let rh = parse_sequence(file_path, D_SCORE_1, &default_scores)?;

    rh.print_background_freqs();
    rh.print_elevated_freqs();
    rh.print_scoring_scheme();

    let rh_custom = parse_sequence(file_path, D_SCORE_2, &rh.get_scoring_scheme())?;
    println!("\nReal data:");
    rh_custom.print_score_histogram();

    let rh_simulated = rh.simulate_new(D_SCORE_2);
    println!("\nSimulated data:");
    rh_simulated.print_score_histogram();
    println!("\nRatios of simulated data:");
    rh_simulated.print_score_ratios();

    Ok(())
}

struct ReadHistogram {
    segs: Vec<(isize, isize, f64)>,
    non_elevated_copies: HashMap<isize, isize>,
    elevated_copies: HashMap<isize, isize>,
}

impl ReadHistogram {
    fn new() -> Self {
        return ReadHistogram {
            segs: Vec::new(),
            non_elevated_copies: HashMap::new(),
            elevated_copies: HashMap::new(),
        };
    }

    fn simulate_new(&self, d_score: f64) -> Self {
        let background_freqs = self.get_background_freq();
        let scoring_scheme = self.get_scoring_scheme();
        let total = self.get_total_background() as isize;
        let mut rh = ReadHistogram::new();
        let mut cum: f64 = 0.0;
        let mut max: f64 = 0.0;
        let mut start: isize = 1;
        let mut end: isize = 1;
        let mut rng = rand::thread_rng();

        for i in 0..total {
            let rnd = rng.gen::<f64>();
            let cnt: isize = if rnd < background_freqs[&0] {
                0
            } else if rnd < background_freqs[&0] + background_freqs[&1] {
                1
            } else if rnd < background_freqs[&0] + background_freqs[&1] + background_freqs[&2] {
                2
            } else {
                3
            };

            cum += get_read_score(cnt, &scoring_scheme);
            if cum >= max {
                max = cum;
                end = i;
            }

            if cum <= 0.0 || cum <= max + d_score || i == total - 1 {
                if max >= -d_score {
                    rh.segs.push((start, end, max));
                }
                max = 0.0;
                cum = 0.0;
                start = i + 1;
                end = i + 1;
            }
        }

        return rh;
    }

    fn get_total_background(&self) -> f64 {
        let total_elevated = self.elevated_copies.values().fold(0.0, |t, v| t + *v as f64);
        let total_non_elevated = self.non_elevated_copies.values().fold(0.0, |t, v| t + *v as f64);
        return total_elevated + total_non_elevated - BACKGROUND_N;
    }

    fn get_background_freq(&self) -> HashMap<isize, f64> {
        let total = self.get_total_background();
        let mut freqs = HashMap::new();
        for (key, value) in self.elevated_copies.iter() {
            let mut count = (value + self.non_elevated_copies[key]) as f64;
            if *key == 0 {
                count -= BACKGROUND_N;
            }
            freqs.insert(*key, count / total);
        }

        return freqs;
    }

    fn get_evelated_freq(&self) -> HashMap<isize, f64> {
        let total = self.elevated_copies.values().fold(0.0, |t, v| t + *v as f64);
        let mut freqs = HashMap::new();
        for (key, value) in self.elevated_copies.iter() {
            freqs.insert(*key, *value as f64 / total);
        }

        return freqs;
    }

    fn get_scoring_scheme(&self) -> HashMap<isize, f64> {
        let mut scores = HashMap::new();
        let background_freqs = self.get_background_freq();
        for (key, value) in self.get_evelated_freq().iter() {
            scores.insert(*key, value.log2() - background_freqs[key].log2());
        }
        return scores;
    }

    fn print_background_freqs(&self) {
        println!("\nBackground frequencies:");
        for (cnt, freq) in self.get_background_freq().iter().sorted_by_key(|(k, _)| *k) {
            if *cnt == 3 {
                println!(">={cnt}={freq:.4}")
            } else {
                println!("{cnt}={freq:.4}")
            }
        }
    }

    fn print_elevated_freqs(&self) {
        println!("\nTarget frequencies:");
        for (cnt, freq) in self.get_evelated_freq().iter().sorted_by_key(|(k, _)| *k) {
            if *cnt == 3 {
                println!(">={cnt}={freq:.4}")
            } else {
                println!("{cnt}={freq:.4}")
            }
        }
    }

    fn print_scoring_scheme(&self) {
        println!("\nScoring scheme:");
        for (cnt, score) in self.get_scoring_scheme().iter().sorted_by_key(|(k, _)| *k) {
            if *cnt == 3 {
                println!(">={cnt}={score:.4}")
            } else {
                println!("{cnt}={score:.4}")
            }
        }
    }

    fn print_score_histogram(&self) {
        for i in 5..31 {
            let count = self.segs.iter().fold(0, |t, (_, _, score)| if *score >= i as f64 { t + 1 } else { t });
            println!("{i} {count}");
        }
    }

    fn print_score_ratios(&self) {
        let mut prev_count = -1.0;
        for i in 5..31 {
            let count = self.segs.iter().fold(0, |t, (_, _, score)| if *score >= i as f64 { t + 1 } else { t }) as f64;
            if prev_count >= 0.0 {
                let ratio = if count > 0.0 { prev_count / count } else { -1.0 };
                println!("N_seg({})/N_seg({i}) {ratio:.2}", i - 1);
            }
            prev_count = count;
        }
    }
}

#[inline]
fn get_read_score(reads: isize, scoring_scheme: &HashMap<isize, f64>) -> f64 {
    return if reads >= 3 { return scoring_scheme[&3] } else { scoring_scheme[&reads] };
}

fn parse_sequence(file_path: &str, d_score: f64, scoring_scheme: &HashMap<isize, f64>) -> Result<ReadHistogram, Error> {
    let mut rh = ReadHistogram::new();
    let mut cum: f64 = 0.0;
    let mut max: f64 = 0.0;
    let mut start: isize = 1;
    let mut end: isize = 1;

    let (mut c_0, mut c_1, mut c_2, mut c_3) = (0, 0, 0, 0);
    let (mut t_0, mut t_1, mut t_2, mut t_3) = (0, 0, 0, 0);
    let (mut m_0, mut m_1, mut m_2, mut m_3) = (0, 0, 0, 0);

    let lines = read::lines(file_path)?;
    for line in lines.enumerate().with_position() {
        let (is_last, line) = match line {
            Position::Middle((_, res)) => (false, res),
            Position::Last((_, res)) => (true, res),
            Position::First((_, res)) => (false, res),
            Position::Only((_, res)) => (true, res),
        };
        if let Ok(ip) = line {
            let mut iter = ip.split_whitespace();
            let _chr = iter.next().unwrap();
            let pos: isize = iter.next().unwrap().parse().unwrap();
            let cnt: isize = iter.next().unwrap().parse().unwrap();

            match cnt {
                0 => {
                    c_0 += 1;
                    t_0 += 1;
                }
                1 => {
                    c_1 += 1;
                    t_1 += 1;
                }
                2 => {
                    c_2 += 1;
                    t_2 += 1;
                }
                _ => {
                    // >=3
                    c_3 += 1;
                    t_3 += 1;
                }
            }

            cum += get_read_score(cnt, scoring_scheme);
            if cum >= max {
                max = cum;
                end = pos;
                (m_0, m_1, m_2, m_3) = (c_0, c_1, c_2, c_3)
            }

            if cum <= 0.0 || cum <= max + d_score || is_last {
                if max >= -d_score {
                    rh.segs.push((start, end, max));
                    *rh.elevated_copies.entry(0).or_default() += m_0;
                    *rh.elevated_copies.entry(1).or_default() += m_1;
                    *rh.elevated_copies.entry(2).or_default() += m_2;
                    *rh.elevated_copies.entry(3).or_default() += m_3;
                }
                (c_0, c_1, c_2, c_3) = (0, 0, 0, 0);
                max = 0.0;
                cum = 0.0;
                start = pos + 1;
                end = pos + 1;
            }
        }
    }

    *rh.non_elevated_copies.entry(0).or_default() += t_0 - rh.elevated_copies[&0];
    *rh.non_elevated_copies.entry(1).or_default() += t_1 - rh.elevated_copies[&1];
    *rh.non_elevated_copies.entry(2).or_default() += t_2 - rh.elevated_copies[&2];
    *rh.non_elevated_copies.entry(3).or_default() += t_3 - rh.elevated_copies[&3];

    Ok(rh)
}
