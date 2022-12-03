use std::collections::HashMap;
use std::io::Error;
use crate::util::{read};

const ALL_CHAR: char = '*';

pub fn run(file_name: &str) -> Result<(), Error>{
    let base_counts = count_bases(file_name)?;
    for key in [ALL_CHAR, 'A', 'C', 'G', 'T', 'N'] {
        println!("{}={}", key, base_counts.get(&key).unwrap_or(&0))
    }

    Ok(())
}

fn count_bases(file_name: &str) -> Result<HashMap<char, u32>, Error> {
    let lines = read::lines(file_name)?;
    let mut base_counts = HashMap::new();
    let mut total_count: u32 = 0;
    for line in lines {
        if let Ok(ip) = line {
            if ip.starts_with('>') {
                continue;
            }
            for c in ip.to_uppercase().chars() {
                *base_counts.entry(c).or_default() += 1;
                total_count += 1;
            }
        }
    }
    base_counts.insert(ALL_CHAR, total_count);
    Ok(base_counts)
}
