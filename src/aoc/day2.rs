/***
--- Day 2: Rock Paper Scissors ---

The Elves begin to set up camp on the beach. To decide whose tent gets to be closest to the snack storage,
a giant Rock Paper Scissors tournament is already in progress.

Rock Paper Scissors is a game between two players. Each game contains many rounds; in each round,
the players each simultaneously choose one of Rock, Paper, or Scissors using a hand shape.
Then, a winner for that round is selected: Rock defeats Scissors, Scissors defeats Paper, and Paper defeats Rock.
If both players choose the same shape, the round instead ends in a draw.

Appreciative of your help yesterday, one Elf gives you an encrypted strategy guide (your puzzle input) that they say will be sure to help you win.
"The first column is what your opponent is going to play: A for Rock, B for Paper, and C for Scissors.
The second column--" Suddenly, the Elf is called away to help with someone's tent.

The second column, you reason, must be what you should play in response: X for Rock, Y for Paper, and Z for Scissors.
Winning every time would be suspicious, so the responses must have been carefully chosen.

The winner of the whole tournament is the player with the highest score. Your total score is the sum of your scores for each round.
The score for a single round is the score for the shape you selected (1 for Rock, 2 for Paper, and 3 for Scissors) plus the score for
the outcome of the round (0 if you lost, 3 if the round was a draw, and 6 if you won).

Since you can't be sure if the Elf is trying to help you or trick you, you should calculate the score you
would get if you were to follow the strategy guide.

For example, suppose you were given the following strategy guide:

A Y
B X
C Z

This strategy guide predicts and recommends the following:

    In the first round, your opponent will choose Rock (A), and you should choose Paper (Y). This ends in a win for you with a score of 8 (2 because you chose Paper + 6 because you won).
    In the second round, your opponent will choose Paper (B), and you should choose Rock (X). This ends in a loss for you with a score of 1 (1 + 0).
    The third round is a draw with both players choosing Scissors, giving you a score of 3 + 3 = 6.

In this example, if you were to follow the strategy guide, you would get a total score of 15 (8 + 1 + 6).

What would your total score be if everything goes exactly according to your strategy guide?

https://adventofcode.com/2022/day/2
 */

use std::error::Error;

use crate::util::read;

const ELF_ROCK: &str = "A";
const ELF_PAPER: &str = "B";
const ELF_SCISSORS: &str = "C";

const SELF_ROCK: &str = "X";
const SELF_PAPER: &str = "Y";
const SELF_SCISSORS: &str = "Z";

#[derive(Clone, Copy, PartialEq)]
enum Move {
    Rock,
    Paper,
    Scissors
}

pub fn run(file_name: &str) -> Result<(), Box<dyn Error>> {
    let lines = read::lines(file_name)?;
    let mut score = 0;
    for line in lines {
        if let Ok(ip) = line {
            let mut moves = ip.split(" ");
            let elf_move = translate_elf_move(moves.next().unwrap());
            let self_move = translate_self_move(moves.next().unwrap());
            score += calc_score(elf_move, self_move);
        }
    }
    println!("Got a score of '{}'", score);
    Ok(())
}

fn translate_elf_move(elf_move: &str) -> Move {
    match elf_move {
        ELF_ROCK=>Move::Rock,
        ELF_PAPER=>Move::Paper,
        ELF_SCISSORS=>Move::Scissors,
        _=>panic!("'{}' is not a legal elf move!", elf_move)
    }
}

fn translate_self_move(self_move: &str) -> Move {
    match self_move {
        SELF_ROCK=>Move::Rock,
        SELF_PAPER=>Move::Paper,
        SELF_SCISSORS=>Move::Scissors,
        _=>panic!("'{}' is not a legal self move!", self_move)
    }
}

fn calc_score(elf_move: Move, self_move: Move) -> i32 {
    return move_score(self_move) + game_score(elf_move, self_move)
}

fn move_score(self_move: Move) -> i32 {
    match self_move {
        Move::Rock=>1,
        Move::Paper=>2,
        Move::Scissors=>3,
    }
}

fn game_score(elf_move: Move, self_move: Move) -> i32 {
    return if elf_move == self_move {
        return 3 // Draw
    } else if (elf_move == Move::Rock && self_move == Move::Paper) ||
        (elf_move == Move::Paper && self_move == Move::Scissors) ||
        (elf_move == Move::Scissors && self_move == Move::Rock) {
        6 // Win
    } else {
        0 // Lose
    }
}
