use gb_io::feature_kind;
use gb_io::reader::SeqReader;
use gb_io::seq::Location;
use itertools::Itertools;
use std::cmp::{max, min};
use std::collections::HashMap;
use std::fs::File;
use std::io::Error;

const BASE_OFFSET: i64 = 10;

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
            'a' | 'A' => Base::A,
            'c' | 'C' => Base::C,
            'g' | 'G' => Base::G,
            't' | 'T' => Base::T,
            _ => Base::N,
        }
    }
}

struct PositionalDistribution {
    counts: HashMap<i64, HashMap<Base, usize>>,
    cds_locs: Vec<Location>,
    freqs: HashMap<i64, HashMap<Base, f64>>,
    weights: HashMap<i64, HashMap<Base, f64>>,
    background: BackgroundDistribution,
}

impl PositionalDistribution {
    fn new(counts: HashMap<i64, HashMap<Base, usize>>, cds_locs: Vec<Location>, background: BackgroundDistribution) -> PositionalDistribution {
        let mut freqs: HashMap<i64, HashMap<Base, f64>> = HashMap::new();
        counts.iter().for_each(|(p, m)| {
            let pos_count = m.iter().fold(0, |t, (b, c)| if *b != Base::N { t + c } else { t });
            m.iter().for_each(|(b, count)| {
                let freq = *count as f64 / pos_count as f64;
                freqs.entry(*p).or_default().insert(*b, freq);
            });
        });

        let mut weights: HashMap<i64, HashMap<Base, f64>> = HashMap::new();
        freqs.iter().for_each(|(p, m)| {
            m.iter().for_each(|(b, p_site)| {
                let weight = if *p_site == 0.0 {
                    -99.00
                } else {
                    p_site.log2() - background.get_base_freq(*b).log2()
                };
                weights.entry(*p).or_default().insert(*b, weight);
            });
        });

        return PositionalDistribution {
            counts,
            cds_locs,
            freqs,
            weights,
            background,
        };
    }

    fn get_pos_count(&self, pos: i64, base: Base) -> usize {
        return *self.counts.get(&pos).unwrap_or(&HashMap::new()).get(&base).unwrap_or(&0);
    }

    fn get_pos_freq(&self, pos: i64, base: Base) -> f64 {
        return *self.freqs.get(&pos).unwrap_or(&HashMap::new()).get(&base).unwrap_or(&0.0);
    }

    fn get_pos_weight(&self, pos: i64, base: Base) -> f64 {
        return *self.weights.get(&pos).unwrap_or(&HashMap::new()).get(&base).unwrap_or(&0.0);
    }

    fn get_max_score(&self) -> f64 {
        return self
            .weights
            .iter()
            .fold(0.0, |t, (_, m)| t + *m.iter().max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap()).unwrap().1);
    }

    fn in_cds_range(&self, pos: i64) -> bool {
        return self.cds_locs.iter().any(|l| {
            let (start, end) = l.find_bounds().unwrap();
            return start <= pos && pos <= end;
        });
    }

    fn print_pos_count(&self) {
        println!("\nCount Matrix:");
        for (p, _) in self.counts.iter().sorted_by_key(|w| w.0) {
            print!("{} ", p);
            for b in [Base::A, Base::C, Base::G, Base::T] {
                print!("{} ", self.get_pos_count(*p, b));
            }
            println!();
        }
    }

    fn print_pos_freq(&self) {
        println!("\nFrequency Matrix:");
        for (p, _) in self.freqs.iter().sorted_by_key(|w| w.0) {
            print!("{} ", p);
            for b in [Base::A, Base::C, Base::G, Base::T] {
                print!("{:.4} ", self.get_pos_freq(*p, b));
            }
            println!();
        }
    }

    fn print_pos_weight(&self) {
        println!("\nWeight Matrix:");
        for (p, _) in self.weights.iter().sorted_by_key(|w| w.0) {
            print!("{} ", p);
            for b in [Base::A, Base::C, Base::G, Base::T] {
                print!("{:.4} ", self.get_pos_weight(*p, b));
            }
            println!();
        }
    }
}

struct BackgroundDistribution {
    base_counts: HashMap<Base, usize>,
    forward_counts: HashMap<Base, usize>,
    known_count: usize,
}

impl BackgroundDistribution {
    fn new(base_counts: HashMap<Base, usize>, forward_counts: HashMap<Base, usize>) -> BackgroundDistribution {
        let total_count: usize = base_counts.iter().fold(0, |t, (b, c)| if *b != Base::N { t + c } else { t });
        return BackgroundDistribution {
            base_counts,
            forward_counts,
            known_count: total_count,
        };
    }

    fn get_forward_count(&self, base: Base) -> usize {
        return *self.forward_counts.get(&base).unwrap_or(&0);
    }

    fn get_base_count(&self, base: Base) -> usize {
        return *self.base_counts.get(&base).unwrap_or(&0);
    }

    fn get_base_freq(&self, base: Base) -> f64 {
        return self.get_base_count(base) as f64 / self.known_count as f64;
    }

    fn print_base_count(&self) {
        println!("\nNucleotide Histogram:");
        for base in [Base::A, Base::C, Base::G, Base::T, Base::N] {
            println!("{:?}={}", base, self.get_forward_count(base));
        }
    }

    fn print_base_freq(&self) {
        println!("\nBackground Frequency:");
        for base in [Base::A, Base::C, Base::G, Base::T] {
            println!("{:?}={:.4}", base, self.get_base_freq(base));
        }
    }
}

pub fn run(file_path: &str) -> Result<(), Error> {
    let pos_dist = count_positions(file_path)?;

    pos_dist.background.print_base_count();
    pos_dist.background.print_base_freq();
    pos_dist.print_pos_count();
    pos_dist.print_pos_freq();
    pos_dist.print_pos_weight();
    println!("\nMaximum Score: {:.10}", pos_dist.get_max_score());

    score_positions(file_path, &pos_dist)?;

    Ok(())
}

fn count_positions(file_path: &str) -> Result<PositionalDistribution, Error> {
    let file = File::open(file_path).unwrap();
    let mut base_counts: HashMap<Base, usize> = HashMap::new();
    let mut forward_counts: HashMap<Base, usize> = HashMap::new();
    let mut pos_counts: HashMap<i64, HashMap<Base, usize>> = HashMap::new();
    let mut cds_locs: Vec<Location> = Vec::new();

    for seq in SeqReader::new(file) {
        let seq = seq.unwrap();
        seq.seq.iter().for_each(|c| {
            *base_counts.entry(Base::from_char(*c as char)).or_default() += 1;
            *forward_counts.entry(Base::from_char(*c as char)).or_default() += 1;
        });
        seq.revcomp()
            .seq
            .iter()
            .for_each(|c| *base_counts.entry(Base::from_char(*c as char)).or_default() += 1);
        seq.features
            .iter()
            .filter(|f| f.kind == feature_kind!("CDS"))
            .map(|f| get_start_location(&f.location))
            .for_each(|l| {
                let mut p = -BASE_OFFSET;
                for c in seq.extract_location(&l).unwrap() {
                    *pos_counts.entry(p).or_default().entry(Base::from_char(c as char)).or_default() += 1;
                    p += 1;
                }
                cds_locs.push(l);
            });
    }
    Ok(PositionalDistribution::new(
        pos_counts,
        cds_locs,
        BackgroundDistribution::new(base_counts, forward_counts),
    ))
}

fn score_positions(file_path: &str, pos_dist: &PositionalDistribution) -> Result<(), Error> {
    let file = File::open(file_path).unwrap();
    let window_size = (BASE_OFFSET * 2 + 1) as usize;
    let mut all_score: HashMap<isize, usize> = HashMap::new();
    let mut cds_score: HashMap<isize, usize> = HashMap::new();
    let mut outliers: Vec<(Location, f64)> = Vec::new();

    for seq in SeqReader::new(file) {
        let seq = seq.unwrap();
        seq.seq.windows(window_size).enumerate().for_each(|(i, w)| {
            let mut p = -BASE_OFFSET;
            let mut score = 0.0;
            for c in w {
                score += pos_dist.get_pos_weight(p, Base::from_char(*c as char));
                p += 1;
            }
            *all_score.entry(bin_score(score)).or_default() += 1;
            if score >= 10.0 {
                let s = i as i64 + BASE_OFFSET + 1;
                if !pos_dist.in_cds_range(s) {
                    outliers.push((Location::single(s), score));
                }
            }
        });
        seq.revcomp().seq.windows(window_size).enumerate().for_each(|(i, w)| {
            let mut p = -BASE_OFFSET;
            let mut score = 0.0;
            for c in w {
                score += pos_dist.get_pos_weight(p, Base::from_char(*c as char));
                p += 1;
            }
            *all_score.entry(bin_score(score)).or_default() += 1;
            if score >= 10.0 {
                let s = i as i64 + BASE_OFFSET + 1;
                if !pos_dist.in_cds_range(s) {
                    outliers.push((Location::Complement(Box::from(Location::single(s))), score));
                }
            }
        });
        pos_dist.cds_locs.iter().for_each(|l| {
            let mut p = -BASE_OFFSET;
            let mut score = 0.0;
            for c in seq.extract_location(&l).unwrap() {
                score += pos_dist.get_pos_weight(p, Base::from_char(c as char));
                p += 1;
            }
            *cds_score.entry(bin_score(score)).or_default() += 1;
        });
    }

    println!("\nScore Histogram CDS:");
    for (p, c) in cds_score.iter().sorted_by_key(|w| w.0) {
        println!("{} {}", p, c);
    }

    println!("\nScore Histogram All:");
    for (p, c) in all_score.iter().sorted_by_key(|w| w.0) {
        println!("{} {}", p, c);
    }

    println!("\nPosition List:");
    for (l, s) in outliers.iter().sorted_by_key(|(l, _)| l.find_bounds().unwrap().0) {
        let strand = match l {
            Location::Complement(_) => 1,
            _ => 0,
        };
        println!("{} {} {:.4}", l.find_bounds().unwrap().0, strand, s);
    }

    Ok(())
}

fn bin_score(score: f64) -> isize {
    return min(max(-51, score.floor() as isize), 51);
}

fn get_start_location(l: &Location) -> Location {
    return if let Location::Complement(lc) = l {
        match lc.as_ref() {
            Location::Range(_, (end, _)) => Location::Complement(Box::from(Location::simple_range(end - BASE_OFFSET - 1, end + BASE_OFFSET))),
            Location::Join(ls) => {
                let mut sls = Vec::new();
                let mut rem = BASE_OFFSET;
                for (i, li) in ls.iter().enumerate().rev() {
                    let (s, e) = li.find_bounds().unwrap();

                    let end = if i == ls.len() - 1 { e + BASE_OFFSET } else { e };
                    let start = if rem - (e - s) <= 0 { e - rem - 1 } else { s };

                    sls.push(Location::simple_range(start, end));
                    rem -= e - s;

                    if rem <= 0 {
                        break;
                    }
                }
                sls.reverse();
                Location::Complement(Box::from(Location::Join(sls)))
            }
            _ => panic!("Don't know how to convert location '{}'", l),
        }
    } else {
        match l {
            Location::Range((start, _), _) => Location::simple_range(start - BASE_OFFSET, start + BASE_OFFSET + 1),
            Location::Join(ls) => {
                let mut sls = Vec::new();
                let mut rem = BASE_OFFSET;
                for (i, li) in ls.iter().enumerate() {
                    let (s, e) = li.find_bounds().unwrap();

                    let start = if i == 0 { s - BASE_OFFSET } else { s };
                    let end = if rem - (e - s) <= 0 { s + rem + 1 } else { e };

                    sls.push(Location::simple_range(start, end));
                    rem -= e - s;

                    if rem <= 0 {
                        break;
                    }
                }
                Location::Join(sls)
            }
            _ => panic!("Don't know how to convert location '{}'", l),
        }
    };
}
