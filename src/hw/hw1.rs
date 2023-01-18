use crate::util::read;
use itertools::Itertools;
use std::cmp::min;
use std::collections::{HashMap, HashSet};
use std::io::Error;
use std::path::Path;

const ALL_ALPHA_KEY: char = '*';

const ALPHA_CHARS: [char; 5] = ['A', 'C', 'G', 'T', 'N'];

#[derive(Debug, Hash, Eq, PartialEq)]
struct SuffixPointer<'a> {
    start_index: usize,
    string: &'a str,
}

impl<'a> SuffixPointer<'a> {
    fn new(start_index: usize, string: &'a str) -> SuffixPointer<'a> {
        return SuffixPointer { start_index, string };
    }

    fn subsequence(&self, offset: usize) -> &'a str {
        return &self.string[self.start_index..self.start_index + offset];
    }

    fn suffix(&self) -> &'a str {
        return &self.string[self.start_index..];
    }
}

pub fn run(file_path1: &str, file_path2: &str) -> Result<(), Error> {
    let file_name1 = file_name_from_path(file_path1);
    println!("Fasta 1: {}", file_name1);
    let seq1 = load_sequence(file_path1)?;

    let file_name2 = file_name_from_path(file_path2);
    println!("\nFasta 2: {}", file_name2);
    let seq2 = load_sequence(file_path2)?;
    let seq2_rev = reverse_complement(&seq2);

    let mut suffix_array = build_suffix_array(Vec::from([seq1.as_str(), seq2.as_str(), seq2_rev.as_str()]));
    suffix_array.sort_unstable_by_key(|s| &s.string[s.start_index..]);

    let mut len_histogram: HashMap<usize, i32> = HashMap::new();
    let mut longest_len = 0;
    let mut longest_matches: HashSet<&SuffixPointer> = HashSet::new();
    let mut match_string = "";
    for (i, s1) in suffix_array.iter().enumerate() {
        if std::ptr::eq(seq1.as_str(), s1.string) {
            let mut s2_idx = 0;
            let mut max_len = 0;
            for j in (0..i).rev() {
                if !std::ptr::eq(seq1.as_str(), suffix_array[j].string) {
                    max_len = count_common_prefix(&s1.suffix(), &suffix_array[j].suffix());
                    s2_idx = j;
                    break;
                }
            }

            for j in i + 1..suffix_array.len() {
                if !std::ptr::eq(seq1.as_str(), suffix_array[j].string) {
                    let len = count_common_prefix(&s1.suffix(), &suffix_array[j].suffix());
                    if len > max_len {
                        max_len = len;
                        s2_idx = j;
                    }
                    break;
                }
            }

            if max_len > longest_len {
                longest_matches.clear();
                longest_len = max_len;
                match_string = s1.subsequence(longest_len);
                longest_matches.insert(&s1);
                longest_matches.insert(&suffix_array[s2_idx]);
            } else if max_len == longest_len {
                longest_matches.insert(&s1);
                longest_matches.insert(&suffix_array[s2_idx]);
            }
            *len_histogram.entry(max_len).or_default() += 1
        }
    }

    println!("\nMatch Length Histogram:");
    for key in len_histogram.keys().sorted() {
        println!("{} {}", key, len_histogram[key]);
    }

    println!("\nThe longest match length: {}", longest_len);
    let unique_seqs: HashSet<&str> = HashSet::from_iter(longest_matches.clone().into_iter().map(|s| s.subsequence(longest_len)));
    println!("Number of match strings: {}", unique_seqs.len());

    println!("\nMatch string: {}", match_string);
    println!("Description: This sequence comes from [look up entry in .gbff annotation file using the position information below]");

    for s in longest_matches {
        if std::ptr::eq(seq1.as_str(), s.string) {
            println!("\nFasta: {}\nPosition: {}\nStrand: forward", file_name1, s.start_index + 1)
        } else if std::ptr::eq(seq2.as_str(), s.string) {
            println!("\nFasta: {}\nPosition: {}\nStrand: forward", file_name2, s.start_index + 1)
        } else {
            println!("\nFasta: {}\nPosition: {}\nStrand: reverse", file_name2, s.start_index + 1)
        }
    }

    Ok(())
}

fn file_name_from_path(file_path: &str) -> &str {
    return Path::new(file_path).file_name().unwrap().to_str().unwrap();
}

fn load_sequence(file_name: &str) -> Result<String, Error> {
    let mut sequence = String::with_capacity(read::file_size(file_name) as usize);
    let mut header: String = String::from("NO HEADER");
    let mut base_counts = HashMap::new();
    let mut alpha_count: u32 = 0;
    let mut non_alpha_count: u32 = 0;

    let lines = read::lines(file_name)?;
    for line in lines {
        if let Ok(ip) = line {
            if ip.starts_with('>') {
                header = ip;
                continue;
            }
            for c in ip.to_uppercase().chars() {
                if ALPHA_CHARS.contains(&c) {
                    sequence.push(c);
                    *base_counts.entry(c).or_default() += 1;
                    alpha_count += 1;
                } else if c != ' ' {
                    // Don't count spaces for some reason?
                    non_alpha_count += 1;
                }
            }
        }
    }

    println!("Non-alphabetic characters: {}", non_alpha_count);
    println!("{}", header);
    println!("{}={}", ALL_ALPHA_KEY, alpha_count);
    for key in ALPHA_CHARS {
        println!("{}={}", key, base_counts.get(&key).unwrap_or(&0))
    }

    Ok(sequence)
}

fn build_suffix_array(strings: Vec<&str>) -> Vec<SuffixPointer> {
    let total_size = strings.iter().fold(0, |sum, s| sum + s.len());
    let mut suffix_array = Vec::with_capacity(total_size);
    for string in strings {
        for (i, _) in string.char_indices() {
            suffix_array.push(SuffixPointer::new(i, string));
        }
    }
    return suffix_array;
}

fn reverse_complement(sequence: &str) -> String {
    let mut rev_complement = String::with_capacity(sequence.len());

    for c in sequence.chars().rev() {
        rev_complement.push(match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => 'N',
        })
    }

    return rev_complement;
}

fn count_common_prefix(s1: &str, s2: &str) -> usize {
    let min_len = min(s1.len(), s2.len());
    let mut common_count = 0;

    for i in 0..min_len {
        if s1.chars().nth(i) == s2.chars().nth(i) {
            common_count += 1;
        } else {
            break;
        }
    }

    return common_count;
}
