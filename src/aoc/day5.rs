/***
--- Day 5: Supply Stacks ---

The expedition can depart as soon as the final supplies have been unloaded from the ships.
Supplies are stored in stacks of marked crates, but because the needed supplies are buried under many other crates,
the crates need to be rearranged.

The ship has a giant cargo crane capable of moving crates between stacks. To ensure none of the crates get crushed or fall over,
the crane operator will rearrange them in a series of carefully-planned steps. After the crates are rearranged,
the desired crates will be at the top of each stack.

The Elves don't want to interrupt the crane operator during this delicate procedure, but they forgot to ask her which crate will end up where,
and they want to be ready to unload them as soon as possible so they can embark.

They do, however, have a drawing of the starting stacks of crates and the rearrangement procedure (your puzzle input). For example:

    [D]
[N] [C]
[Z] [M] [P]
 1   2   3

move 1 from 2 to 1
move 3 from 1 to 3
move 2 from 2 to 1
move 1 from 1 to 2

In this example, there are three stacks of crates. Stack 1 contains two crates: crate Z is on the bottom, and crate N is on top.
Stack 2 contains three crates; from bottom to top, they are crates M, C, and D. Finally, stack 3 contains a single crate, P.

Then, the rearrangement procedure is given. In each step of the procedure, a quantity of crates is moved from one stack to a different stack.
In the first step of the above rearrangement procedure, one crate is moved from stack 2 to stack 1, resulting in this configuration:

[D]
[N] [C]
[Z] [M] [P]
 1   2   3

In the second step, three crates are moved from stack 1 to stack 3. Crates are moved one at a time,
so the first crate to be moved (D) ends up below the second and third crates:

        [Z]
        [N]
    [C] [D]
    [M] [P]
 1   2   3

Then, both crates are moved from stack 2 to stack 1. Again, because crates are moved one at a time, crate C ends up below crate M:

        [Z]
        [N]
[M]     [D]
[C]     [P]
 1   2   3

Finally, one crate is moved from stack 1 to stack 2:

        [Z]
        [N]
        [D]
[C] [M] [P]
 1   2   3

The Elves just need to know which crate will end up on top of each stack; in this example,
the top crates are C in stack 1, M in stack 2, and Z in stack 3, so you should combine these together and give the Elves the message CMZ.

After the rearrangement procedure completes, what crate ends up on top of each stack?

--- Part Two ---

As you watch the crane operator expertly rearrange the crates, you notice the process isn't following your prediction.

Some mud was covering the writing on the side of the crane, and you quickly wipe it away.
The crane isn't a CrateMover 9000 - it's a CrateMover 9001.

The CrateMover 9001 is notable for many new and exciting features: air conditioning, leather seats,
an extra cup holder, and the ability to pick up and move multiple crates at once.

Again considering the example above, the crates begin in the same configuration:

    [D]
[N] [C]
[Z] [M] [P]
 1   2   3

Moving a single crate from stack 2 to stack 1 behaves the same as before:

[D]
[N] [C]
[Z] [M] [P]
 1   2   3

However, the action of moving three crates from stack 1 to stack 3 means that those three moved crates stay in the same order,
resulting in this new configuration:

        [D]
        [N]
    [C] [Z]
    [M] [P]
 1   2   3

Next, as both crates are moved from stack 2 to stack 1, they retain their order as well:

        [D]
        [N]
[C]     [Z]
[M]     [P]
 1   2   3

Finally, a single crate is still moved from stack 1 to stack 2, but now it's crate C that gets moved:

        [D]
        [N]
        [Z]
[M] [C] [P]
 1   2   3

In this example, the CrateMover 9001 has put the crates in a totally different order: MCD.

Before the rearrangement process finishes, update your simulation so that the Elves know where they should stand to be ready to unload the final supplies.
After the rearrangement procedure completes, what crate ends up on top of each stack?

https://adventofcode.com/2022/day/5
 */

use regex::Regex;
use std::error::Error;

use crate::util::read;

const STACK_REGEX_PATTERN: &str = r"[\[\]A-Z]{3}| {4}";
const CRATE_REGEX_PATTERN: &str = r"\[[A-Z]\]";
const ACTION_REGEX_PATTERN: &str = r"move ([0-9]+) from ([0-9]+) to ([0-9]+)";

pub fn run(file_name: &str) -> Result<(), Box<dyn Error>> {
    let lines = read::lines(file_name)?;

    let mut stacks: Vec<Vec<String>> = Vec::new();
    let mut loading_phase = true;

    // Load Crates
    for line in lines {
        if let Ok(ip) = line {
            if loading_phase && ip.contains('1') {
                for stack in &mut stacks {
                    stack.reverse();
                }
                println!("Loaded stacks: {:#?}", stacks);
                loading_phase = false;
            } else if ip.is_empty() {
                continue;
            } else if loading_phase {
                load_crates(ip, &mut stacks);
            } else {
                move_crates_2(ip, &mut stacks);
            }
        }
    }
    println!("Moved stacks: {:#?}", stacks);
    println!("Top stacks:");
    for stack in stacks {
        print!(
            "{}",
            stack[stack.len() - 1].replace("[", "").replace("]", "")
        );
    }
    println!();
    Ok(())
}

fn load_crates(line: String, stacks: &mut Vec<Vec<String>>) {
    let stack_regex = Regex::new(STACK_REGEX_PATTERN).unwrap();
    let crate_regex = Regex::new(CRATE_REGEX_PATTERN).unwrap();
    let mut stack_num = 0;

    for stack_match in stack_regex.find_iter(&line) {
        if stacks.len() <= stack_num {
            stacks.push(Vec::new());
        }
        let stack = stack_match.as_str().trim();
        if let Some(crate_match) = crate_regex.find(stack) {
            stacks[stack_num].push(crate_match.as_str().parse().unwrap());
        }
        stack_num += 1;
    }
}

fn move_crates(line: String, stacks: &mut Vec<Vec<String>>) {
    let action_regex = Regex::new(ACTION_REGEX_PATTERN).unwrap();
    let caps = action_regex.captures(&line).unwrap();
    let n = caps.get(1).unwrap().as_str().parse::<usize>().unwrap();
    let start_stack = caps.get(2).unwrap().as_str().parse::<usize>().unwrap() - 1;
    let end_stack = caps.get(3).unwrap().as_str().parse::<usize>().unwrap() - 1;
    for _i in 0..n {
        let moved_crate = stacks[start_stack].pop().unwrap();
        stacks[end_stack].push(moved_crate);
    }
}

fn move_crates_2(line: String, stacks: &mut Vec<Vec<String>>) {
    let action_regex = Regex::new(ACTION_REGEX_PATTERN).unwrap();
    let caps = action_regex.captures(&line).unwrap();
    let n = caps.get(1).unwrap().as_str().parse::<usize>().unwrap();
    let start_stack = caps.get(2).unwrap().as_str().parse::<usize>().unwrap() - 1;
    let end_stack = caps.get(3).unwrap().as_str().parse::<usize>().unwrap() - 1;

    let end_i = stacks[start_stack].len() - n;
    let moved_crates = stacks[start_stack].split_off(end_i);
    stacks[end_stack].extend(moved_crates);
}
