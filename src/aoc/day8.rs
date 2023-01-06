/***
--- Day 8: Treetop Tree House ---

The expedition comes across a peculiar patch of tall trees all planted carefully in a grid.
The Elves explain that a previous expedition planted these trees as a reforestation effort.
Now, they're curious if this would be a good location for a tree house.

First, determine whether there is enough tree cover here to keep a tree house hidden. To do this,
you need to count the number of trees that are visible from outside the grid when looking directly
along a row or column.

The Elves have already launched a quadcopter to generate a map with the height of each tree (your puzzle input).
For example:

30373
25512
65332
33549
35390

Each tree is represented as a single digit whose value is its height, where 0 is the shortest and 9 is the tallest.

A tree is visible if all of the other trees between it and an edge of the grid are shorter than it.
Only consider trees in the same row or column; that is, only look up, down, left, or right from any given tree.

All of the trees around the edge of the grid are visible - since they are already on the edge,
there are no trees to block the view. In this example, that only leaves the interior nine trees to consider:

    The top-left 5 is visible from the left and top. (It isn't visible from the right or bottom since other trees of height 5 are in the way.)
    The top-middle 5 is visible from the top and right.
    The top-right 1 is not visible from any direction; for it to be visible, there would need to only be trees of height 0 between it and an edge.
    The left-middle 5 is visible, but only from the right.
    The center 3 is not visible from any direction; for it to be visible, there would need to be only trees of at most height 2 between it and an edge.
    The right-middle 3 is visible from the right.
    In the bottom row, the middle 5 is visible, but the 3 and 4 are not.

With 16 trees visible on the edge and another 5 visible in the interior,
a total of 21 trees are visible in this arrangement.

Consider your map; how many trees are visible from outside the grid?

--- Part Two ---

Content with the amount of tree cover available, the Elves just need to know the best spot to build their tree house:
they would like to be able to see a lot of trees.

To measure the viewing distance from a given tree, look up, down, left, and right from that tree;
stop if you reach an edge or at the first tree that is the same height or taller than the tree under consideration.
(If a tree is right on the edge, at least one of its viewing distances will be zero.)

The Elves don't care about distant trees taller than those found by the rules above;
the proposed tree house has large eaves to keep it dry, so they wouldn't be able to see higher than
the tree house anyway.

In the example above, consider the middle 5 in the second row:

30373
25512
65332
33549
35390

    Looking up, its view is not blocked; it can see 1 tree (of height 3).
    Looking left, its view is blocked immediately; it can see only 1 tree (of height 5, right next to it).
    Looking right, its view is not blocked; it can see 2 trees.
    Looking down, its view is blocked eventually; it can see 2 trees (one of height 3, then the tree of height 5 that blocks its view).

A tree's scenic score is found by multiplying together its viewing distance in each of the four directions.
For this tree, this is 4 (found by multiplying 1 * 1 * 2 * 2).

However, you can do even better: consider the tree of height 5 in the middle of the fourth row:

30373
25512
65332
33549
35390

    Looking up, its view is blocked at 2 trees (by another tree with a height of 5).
    Looking left, its view is not blocked; it can see 2 trees.
    Looking down, its view is also not blocked; it can see 1 tree.
    Looking right, its view is blocked at 2 trees (by a massive tree of height 9).

This tree's scenic score is 8 (2 * 2 * 1 * 2); this is the ideal spot for the tree house.

Consider each tree on your map. What is the highest scenic score possible for any tree?

https://adventofcode.com/2022/day/8
 */

use std::cell::RefCell;
use std::error::Error;
use std::rc::Rc;

use crate::util::read;

#[derive(Debug)]
struct Tree {
    height: i32,
    is_visible: bool,
    tallest_up: i32,
    tallest_left: i32,
    tallest_down: i32,
    tallest_right: i32,
}

pub fn run(file_name: &str) -> Result<(), Box<dyn Error>> {
    let lines = read::lines(file_name)?;

    let mut trees = Vec::new();
    let mut num_visible = 0;
    lines.enumerate().for_each(|(_i, line)| {
        if let Ok(ip) = line {
            let mut row = Vec::new();
            ip.chars().enumerate().for_each(|(_j, c)| {
                let height = c.to_digit(10).unwrap() as i32;
                let tree = Rc::new(RefCell::new(Tree {
                    height,
                    is_visible: false,
                    tallest_up: -1,
                    tallest_left: -1,
                    tallest_down: -1,
                    tallest_right: -1,
                }));
                row.push(tree);
            });
            trees.push(row);
        }
    });

    for i in 0..trees.len() {
        for j in 0..trees[i].len() {
            let current_tree = trees[i][j].clone();

            let up_tree = get_tree(&trees, i as i32 - 1, j as i32);
            let (visible_up, tallest_up) = if let Some(tree) = up_tree {
                if current_tree.borrow().height > tree.borrow().tallest_up {
                    (true, current_tree.borrow().height)
                } else {
                    (false, tree.borrow().tallest_up)
                }
            } else {
                (true, current_tree.borrow().height)
            };

            if !current_tree.borrow().is_visible && visible_up {
                num_visible += 1;
            }

            current_tree.borrow_mut().tallest_up = tallest_up;
            current_tree.borrow_mut().is_visible = visible_up;

            let left_tree = get_tree(&trees, i as i32, j as i32 - 1);
            let (visible_left, tallest_left) = if let Some(tree) = left_tree {
                if current_tree.borrow().height > tree.borrow().tallest_left {
                    (true, current_tree.borrow().height)
                } else {
                    (false, tree.borrow().tallest_left)
                }
            } else {
                (true, current_tree.borrow().height)
            };

            if !current_tree.borrow().is_visible && visible_left {
                num_visible += 1;
            }

            current_tree.borrow_mut().tallest_left = tallest_left;
            current_tree.borrow_mut().is_visible = visible_up || visible_left;
        }
    }

    for i in (0..trees.len()).rev() {
        for j in (0..trees[i].len()).rev() {
            let current_tree = trees[i][j].clone();

            let down_tree = get_tree(&trees, (i + 1) as i32, j as i32);
            let (visible_down, tallest_down) = if let Some(tree) = down_tree {
                if current_tree.borrow().height > tree.borrow().tallest_down {
                    (true, current_tree.borrow().height)
                } else {
                    (false, tree.borrow().tallest_down)
                }
            } else {
                (true, current_tree.borrow().height)
            };

            current_tree.borrow_mut().tallest_down = tallest_down;

            let right_tree = get_tree(&trees, i as i32, (j + 1) as i32);
            let (visible_right, tallest_right) = if let Some(tree) = right_tree {
                if current_tree.borrow().height > tree.borrow().tallest_right {
                    (true, current_tree.borrow().height)
                } else {
                    (false, tree.borrow().tallest_right)
                }
            } else {
                (true, current_tree.borrow().height)
            };

            let is_visible = current_tree.borrow().is_visible;
            if !is_visible && (visible_right || visible_down) {
                num_visible += 1;
            }

            current_tree.borrow_mut().tallest_right = tallest_right;
            current_tree.borrow_mut().is_visible = is_visible || visible_right || visible_down;
        }
    }

    for row in trees.iter() {
        for tree in row.iter() {
            print!("{}", tree.borrow().is_visible as i32)
        }
        println!()
    }
    println!("There are '{}' visible trees!", num_visible);
    Ok(())
}

fn get_tree(trees: &Vec<Vec<Rc<RefCell<Tree>>>>, i: i32, j: i32) -> Option<Rc<RefCell<Tree>>> {
    if i < 0 || j < 0 {
        return None;
    } else if let Some(row) = trees.get(i as usize) {
        if let Some(tree) = row.get(j as usize) {
            return Some(tree.clone());
        }
    }
    return None;
}
