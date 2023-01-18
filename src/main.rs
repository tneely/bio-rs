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
            1 => hw::hw1::run("./data/hw/hw1/CP001872.fna", "./data/hw/hw1/CP003913.fna").expect("Homework 1 should pass!"),
            2 => hw::hw2::run("./data/hw/hw1/CP003913.fna").expect("Homework 2 should pass!"),
            _ => panic!("This assignment hasn't been completed!"),
        }
        println!("Homework '{}' completed in '{}' seconds", hw, now.elapsed().as_secs());
    } else if let Some(aoc) = args.aoc {
        println!("Running advent of code day '{}':", aoc);
        let now = Instant::now();
        match aoc {
            1 => aoc::day1::run("./data/aoc/day1/input.txt").expect("Day 1 failed!"),
            2 => aoc::day2::run("./data/aoc/day2/input.txt").expect("Day 2 failed!"),
            3 => aoc::day3::run("./data/aoc/day3/input.txt").expect("Day 3 failed!"),
            4 => aoc::day4::run("./data/aoc/day4/input.txt").expect("Day 4 failed!"),
            5 => aoc::day5::run("./data/aoc/day5/input.txt").expect("Day 5 failed!"),
            6 => aoc::day6::run("./data/aoc/day6/input.txt").expect("Day 6 failed!"),
            7 => aoc::day7::run("./data/aoc/day7/input.txt").expect("Day 7 failed!"),
            8 => aoc::day8::run("./data/aoc/day8/example.txt").expect("Day 8 failed!"),
            _ => panic!("This day hasn't been completed!"),
        }
        println!("Day '{}' completed in '{:#?}'", aoc, now.elapsed());
    } else {
        panic!("How'd you get here?!")
    }
}
