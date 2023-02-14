use itertools::{Itertools, Position};

use crate::util::read;
use std::{collections::HashMap, io::Error};

const D_SCORE: f64 = -20.0;
const S_SCORE: f64 = -D_SCORE;

pub fn run(file_path1: &str) -> Result<(), Error> {
    let rh = parse_sequence(file_path1)?;

    rh.print_seg_list();
    rh.print_annotations();
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

    fn print_annotations(&self) {
        println!("\nAnnotations:");
        self.segs
            .iter()
            .sorted_unstable_by(|(_, _, score1), (_, _, score2)| score2.partial_cmp(score1).unwrap())
            .take(3)
            .for_each(|(start, end, _)| {
                println!("\nStart: {start}");
                println!("End: {end}");
                println!("Description: TODO");
            })
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

#[inline]
fn get_read_score(reads: isize) -> f64 {
    return match reads {
        0 => -0.1077,
        1 => 0.4772,
        2 => 1.0622,
        _ => 1.6748, // >= 3
    };
}

fn parse_sequence(file_path: &str) -> Result<ReadHistogram, Error> {
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

            cum += get_read_score(cnt);
            if cum >= max {
                max = cum;
                end = pos;
                (m_0, m_1, m_2, m_3) = (c_0, c_1, c_2, c_3)
            }

            if cum <= 0.0 || cum <= max + D_SCORE || is_last {
                if max >= S_SCORE {
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
