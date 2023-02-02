use std::fs::{self, File};
use std::io::{self, BufRead};
use std::path::Path;

pub fn file(file_path: &str) -> String {
    return fs::read_to_string(file_path).expect(&format!("Unable to read file '{}'", file_path));
}

pub fn file_size(file_path: &str) -> u64 {
    return fs::metadata(file_path).expect(&format!("Unable to read file '{}'", file_path)).len();
}

// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
pub fn lines<P>(file_path: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(file_path)?;
    Ok(io::BufReader::new(file).lines())
}

pub fn file_name_from_path(file_path: &str) -> &str {
    return Path::new(file_path).file_name().unwrap().to_str().unwrap();
}
