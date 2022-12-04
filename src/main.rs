extern crate core;

use std::time::Instant;

use clap::{ArgGroup, Parser};

mod aoc;
mod hw;
mod util;

// Simple program to run assignments
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
#[command(group(
ArgGroup::new("vers")
.required(true)
.args(["hw", "aoc"]),
))]
struct Args {
    /// Genome 540 homework assignment to run
    #[arg(long)]
    hw: Option<u8>,

    /// Advent of code day to run
    #[arg(long)]
    aoc: Option<u8>,
}

fn main() {
    let args = Args::parse();

    if let Some(hw) = args.hw {
        println!("Running homework assignment '{}':", hw);
        let now = Instant::now();
        match hw {
            0 => hw::hw0::run("./data/hw/hw0/test.fna").expect("Homework 0 failed!"),
            _ => panic!("This assignment hasn't been completed!"),
        }
        println!(
            "Homework '{}' completed in '{}' seconds",
            hw,
            now.elapsed().as_secs()
        );
    } else if let Some(aoc) = args.aoc {
        println!("Running advent of code day '{}':", aoc);
        let now = Instant::now();
        match aoc {
            1 => aoc::day1::run("./data/aoc/day1/input.txt").expect("Day 1 failed!"),
            2 => aoc::day2::run("./data/aoc/day2/input.txt").expect("Day 2 failed!"),
            3 => aoc::day3::run("./data/aoc/day3/input.txt").expect("Day 3 failed!"),
            _ => panic!("This day hasn't been completed!"),
        }
        println!("Day '{}' completed in '{:#?}'", aoc, now.elapsed());
    } else {
        panic!("How'd you get here?!")
    }
}
