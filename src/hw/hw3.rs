use gb_io::feature_kind;
use gb_io::reader::SeqReader;
use gb_io::seq::Location;
use std::collections::HashMap;
use std::fs::File;
use std::io::Error;
use std::ops::ControlFlow;

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
    freqs: HashMap<i64, HashMap<Base, f64>>,
    weights: HashMap<i64, HashMap<Base, f64>>,
    background: BackgroundDistribution,
}

impl PositionalDistribution {
    fn new(counts: HashMap<i64, HashMap<Base, usize>>, background: BackgroundDistribution) -> PositionalDistribution {
        let mut freqs: HashMap<i64, HashMap<Base, f64>> = HashMap::new();

        let mut weights: HashMap<i64, HashMap<Base, f64>> = HashMap::new();

        return PositionalDistribution {
            counts,
            freqs,
            weights,
            background,
        };
    }

    fn get_pos_count(&self, pos: i64, base: Base) -> usize {
        return *self.counts.get(&pos).unwrap_or(&HashMap::new()).get(&base).unwrap_or(&0);
    }

    fn print_pos_count(&self) {
        println!("\nCount Matrix:");
        for p in -BASE_OFFSET..BASE_OFFSET + 1 {
            print!("{} ", p);
            for b in [Base::A, Base::C, Base::G, Base::T] {
                print!("{} ", self.get_pos_count(p, b));
            }
            println!();
        }
    }
}

struct BackgroundDistribution {
    counts: HashMap<Base, usize>,
    known_count: usize,
}

impl BackgroundDistribution {
    fn new(counts: HashMap<Base, usize>) -> BackgroundDistribution {
        let total_count: usize = counts.iter().fold(0, |t, (b, c)| if *b != Base::N { t + c } else { t });
        return BackgroundDistribution {
            counts,
            known_count: total_count,
        };
    }

    fn get_base_count(&self, base: Base) -> usize {
        return *self.counts.get(&base).unwrap_or(&0);
    }

    fn print_base_count(&self) {
        println!("\nNucleotide Histogram:");
        for base in [Base::A, Base::C, Base::G, Base::T, Base::N] {
            println!("{:?}={}", base, self.get_base_count(base));
        }
    }

    fn print_base_freq(&self) {
        println!("\nBackground Frequencies:");
        for base in [Base::A, Base::C, Base::G, Base::T] {
            println!("{:?}={:.4}", base, self.get_base_count(base) as f64 / self.known_count as f64);
        }
    }
}

pub fn run(file_path: &str) -> Result<(), Error> {
    let pos_dist = count_positions(file_path)?;

    pos_dist.background.print_base_count();
    pos_dist.background.print_base_freq();
    pos_dist.print_pos_count();

    Ok(())
}

fn count_positions(file_path: &str) -> Result<PositionalDistribution, Error> {
    let file = File::open(file_path).unwrap();
    let mut base_counts: HashMap<Base, usize> = HashMap::new();
    let mut pos_counts: HashMap<i64, HashMap<Base, usize>> = HashMap::new();

    for seq in SeqReader::new(file) {
        let seq = seq.unwrap();
        seq.seq.iter().for_each(|c| *base_counts.entry(Base::from_char(*c as char)).or_default() += 1);
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
            });
    }
    Ok(PositionalDistribution::new(pos_counts, BackgroundDistribution::new(base_counts)))
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

                    let end = if i == 0 { e + BASE_OFFSET } else { e };
                    let start = if rem - (e - s) <= 0 { e - rem - 1 } else { s };

                    sls.push(Location::simple_range(start, end));
                    rem -= e - s;

                    if rem <= 0 {
                        break;
                    }
                }
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
