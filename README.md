## bio-rs

This repository serves as a workspace for learning Rust.

It was originally created to follow along with
UW's [Genome 540 Introduction to Computational Molecular Biology:
Genome and Protein Sequence Analysis](http://bozeman.mbt.washington.edu/compbio/mbt599/) course, but has recently been co-opted
for [Advent of Code](https://adventofcode.com/).

### Usage

To run a specific homework assignment, run:

```shell
cargo run -- --hw=0
```

For an advent of code day, run:

```shell
cargo run -- --aoc=1
```

Alternatively, you can build the project and run the binary directly:

```shell
cargo build
./target/debug/bio-rs --hw=0
```

If you want to run an optimized build, use the `--release` flag and run from `./target/release/bio-rs`.

### Development

To add a new package, you can run:

```shell
cargo add <<package-name>>
```
