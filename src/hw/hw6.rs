use itertools::{Itertools, Position};

use crate::util::read;
use std::{collections::HashMap, io::Error};

const D_SCORE: f64 = -33.219;
const S_SCORE: f64 = -D_SCORE;

pub fn run(file_path1: &str) -> Result<(), Error> {
    let rh = parse_sequence(file_path1)?;

    rh.print_seg_list();
    println!("\n Annotations:\nTODO");
    rh.print_non_elevated();
    rh.print_elevated();

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

    fn print_seg_list(&self) {
        let elevated = self.segs.len();
        let non_elevated = elevated + 1;

        println!("Segment Histogram:");
        println!("Non-Elevated CN Segments={non_elevated}");
        println!("Elevated CN Segments={elevated}");

        println!("\nSegment List:");
        for (p_start, p_end, s) in self.segs.iter() {
            if *s >= S_SCORE {
                println!("{p_start} {p_end} {s:.2}");
            }
        }
    }

    fn print_non_elevated(&self) {
        println!("\nRead start histogram for non-elevated copy-number segments:");
        for (cnt, tot) in self.non_elevated_copies.iter().sorted() {
            if *cnt == 3 {
                println!(">={cnt}={tot}")
            } else {
                println!("{cnt}={tot}")
            }
        }
    }

    fn print_elevated(&self) {
        println!("\nRead start histogram for elevated copy-number segments:");
        for (cnt, tot) in self.elevated_copies.iter().sorted() {
            if *cnt == 3 {
                println!(">={cnt}={tot}")
            } else {
                println!("{cnt}={tot}")
            }
        }
    }
}

fn get_read_score(reads: isize) -> f64 {
    return match reads {
        0 => -0.3464,
        1 => 0.2488,
        2 => 0.8439,
        _ => 1.5337, // >= 3
    };
}

fn parse_sequence(file_path: &str) -> Result<ReadHistogram, Error> {
    let mut rh = ReadHistogram::new();
    let mut cum: f64 = 0.0;
    let mut max: f64 = 0.0;
    let mut start: isize = 1;
    let mut end: isize = 1;
    let mut copies: HashMap<isize, isize> = HashMap::new();
    let mut copies_total: HashMap<isize, isize> = HashMap::new();
    let mut copies_at_max: HashMap<isize, isize> = HashMap::new();

    let lines = read::lines(file_path)?;
    for line in lines.enumerate().with_position() {
        let (is_last, line) = match line {
            Position::Only((_, res)) => (true, res),
            Position::Last((_, res)) => (true, res),
            Position::First((_, res)) => (false, res),
            Position::Middle((_, res)) => (false, res),
        };
        if let Ok(ip) = line {
            let mut iter = ip.split_whitespace();
            let _chr = iter.next().unwrap();
            let pos: isize = iter.next().unwrap().parse().unwrap();
            let cnt: isize = iter.next().unwrap().parse().unwrap();

            let key = if cnt < 3 { cnt } else { 3 };
            *copies.entry(key).or_default() += 1;
            *copies_total.entry(key).or_default() += 1;

            cum += get_read_score(cnt);
            if cum >= max {
                max = cum;
                end = pos;
                copies_at_max = copies.clone();
            }

            if cum <= 0.0 || cum <= max + D_SCORE || is_last {
                if max >= S_SCORE {
                    rh.segs.push((start, end, max));
                    for (k, v) in copies_at_max.iter() {
                        *rh.elevated_copies.entry(*k).or_default() += v;
                    }
                }
                copies.clear();
                max = 0.0;
                cum = 0.0;
                start = pos + 1;
                end = pos + 1;
            }
        }
    }

    for (k, v) in copies_total.iter() {
        *rh.non_elevated_copies.entry(*k).or_default() += v - rh.elevated_copies[k];
    }

    Ok(rh)
}
