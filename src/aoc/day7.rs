/***
--- Day 7: No Space Left On Device ---

You can hear birds chirping and raindrops hitting leaves as the expedition proceeds.
Occasionally, you can even hear much louder sounds in the distance; how big do the animals get out here, anyway?

The device the Elves gave you has problems with more than just its communication system.
You try to run a system update:

$ system-update --please --pretty-please-with-sugar-on-top
Error: No space left on device

Perhaps you can delete some files to make space for the update?

You browse around the filesystem to assess the situation and save the resulting terminal output (your puzzle input).
For example:

$ cd /
$ ls
dir a
14848514 b.txt
8504156 c.dat
dir d
$ cd a
$ ls
dir e
29116 f
2557 g
62596 h.lst
$ cd e
$ ls
584 i
$ cd ..
$ cd ..
$ cd d
$ ls
4060174 j
8033020 d.log
5626152 d.ext
7214296 k

The filesystem consists of a tree of files (plain data) and directories (which can contain other directories or files).
The outermost directory is called /. You can navigate around the filesystem, moving into or out of directories and
listing the contents of the directory you're currently in.

Within the terminal output, lines that begin with $ are commands you executed, very much like some modern computers:

    cd means change directory. This changes which directory is the current directory, but the specific result depends on the argument:
        cd x moves in one level: it looks in the current directory for the directory named x and makes it the current directory.
        cd .. moves out one level: it finds the directory that contains the current directory, then makes that directory the current directory.
        cd / switches the current directory to the outermost directory, /.
    ls means list. It prints out all of the files and directories immediately contained by the current directory:
        123 abc means that the current directory contains a file named abc with size 123.
        dir xyz means that the current directory contains a directory named xyz.

Given the commands and output in the example above, you can determine that the filesystem looks visually like this:

- / (dir)
  - a (dir)
    - e (dir)
      - i (file, size=584)
    - f (file, size=29116)
    - g (file, size=2557)
    - h.lst (file, size=62596)
  - b.txt (file, size=14848514)
  - c.dat (file, size=8504156)
  - d (dir)
    - j (file, size=4060174)
    - d.log (file, size=8033020)
    - d.ext (file, size=5626152)
    - k (file, size=7214296)

Here, there are four directories: / (the outermost directory), a and d (which are in /), and e (which is in a).
These directories also contain files of various sizes.

Since the disk is full, your first step should probably be to find directories that are good candidates for deletion.
To do this, you need to determine the total size of each directory. The total size of a directory is the sum of the sizes of the files it contains,
directly or indirectly. (Directories themselves do not count as having any intrinsic size.)

The total sizes of the directories above can be found as follows:

    The total size of directory e is 584 because it contains a single file i of size 584 and no other directories.
    The directory a has total size 94853 because it contains files f (size 29116), g (size 2557), and h.lst (size 62596), plus file i indirectly (a contains e which contains i).
    Directory d has total size 24933642.
    As the outermost directory, / contains every file. Its total size is 48381165, the sum of the size of every file.

To begin, find all of the directories with a total size of at most 100000, then calculate the sum of their total sizes.
In the example above, these directories are a and e; the sum of their total sizes is 95437 (94853 + 584).
(As in this example, this process can count files more than once!)

Find all of the directories with a total size of at most 100000.
What is the sum of the total sizes of those directories?

--- Part Two ---

Now, you're ready to choose a directory to delete.

The total disk space available to the filesystem is 70000000. To run the update, you need unused space of at least 30000000. You need to find a directory you can delete that will free up enough space to run the update.

In the example above, the total size of the outermost directory (and thus the total amount of used space) is 48381165; this means that the size of the unused space must currently be 21618835, which isn't quite the 30000000 required by the update. Therefore, the update still requires a directory with total size of at least 8381165 to be deleted before it can run.

To achieve this, you have the following options:

    Delete directory e, which would increase unused space by 584.
    Delete directory a, which would increase unused space by 94853.
    Delete directory d, which would increase unused space by 24933642.
    Delete directory /, which would increase unused space by 48381165.

Directories e and a are both too small; deleting them would not free up enough space. However, directories d and / are both big enough! Between these, choose the smallest: d, increasing unused space by 24933642.

Find the smallest directory that, if deleted, would free up enough space on the filesystem to run the update. What is the total size of that directory?

https://adventofcode.com/2022/day/7
 */

use regex::Regex;
use std::cell::RefCell;
use std::error::Error;
use std::iter::repeat;
use std::rc::Rc;

use crate::util::read;

const CD_REGEX_PATTERN: &str = r"\$ cd (.*)";
const FILE_REGEX_PATTERN: &str = r"([0-9]+) (.*)";
const ROOT_DIR: &str = "/";
const UP_DIR: &str = "..";

#[derive(Debug)]
struct Directory {
    name: String,
    size: usize,
    child_size: usize,
    children: Vec<Rc<RefCell<Directory>>>,
    parent: Option<Rc<RefCell<Directory>>>,
}

impl Directory {
    fn new(name: String, parent: Option<Rc<RefCell<Directory>>>) -> Directory {
        return Directory {
            name,
            size: 0,
            child_size: 0,
            parent,
            children: Vec::new(),
        };
    }

    fn add_child(&mut self, child: Rc<RefCell<Directory>>) {
        self.children.push(child);
    }

    fn add_file(&mut self, file_size: usize) {
        self.size += file_size;
    }

    fn update_size(&mut self) -> usize {
        self.child_size += self
            .children
            .iter()
            .map(|c| c.borrow_mut().update_size())
            .reduce(|a, b| a + b)
            .unwrap_or(0);
        return self.child_size + self.size;
    }

    fn print(&self, depth: usize) {
        let total_size = self.size + self.child_size;
        let spaces = repeat("  ").take(depth * 2).collect::<String>();
        println!(
            "{}- {} (size={}, files={})",
            spaces, self.name, total_size, self.size
        );
        self.children
            .iter()
            .for_each(|child| child.borrow().print(depth + 1));
    }

    fn sum_subdirs_of_max_size(&self, max_size: usize) -> usize {
        return filter_val(self.size + self.child_size, max_size)
            + self
                .children
                .iter()
                .map(|c| c.borrow_mut().sum_subdirs_of_max_size(max_size))
                .reduce(|a, b| a + b)
                .unwrap_or(0);
    }

    fn smallest_dir_larger_than(&self, min_size: usize, smallest: Option<usize>) -> Option<usize> {
        let mut smaller = smallest.clone();
        let total_size = self.size + self.child_size;

        if total_size < min_size {
            println!(
                "Dir {} is too small to delete ({})",
                self.name,
                self.size + self.child_size
            );
            return smallest;
        } else if let Some(smol) = smallest.clone() {
            if total_size < smol {
                println!(
                    "Found smaller dir {} ({})",
                    self.name,
                    self.size + self.child_size
                );
                smaller = Some(total_size);
            }
        } else {
            println!(
                "Found smaller dir {} ({})",
                self.name,
                self.size + self.child_size
            );
            smaller = Some(total_size);
        }

        for child in self.children.iter() {
            if let Some(smol1) = child.borrow().smallest_dir_larger_than(min_size, smaller) {
                if let Some(smol2) = smaller {
                    if smol1 < smol2 {
                        println!(
                            "Found smaller dir {} ({})",
                            child.borrow().name,
                            child.borrow().size + child.borrow().child_size
                        );
                        smaller = Some(smol1)
                    }
                } else {
                    println!(
                        "Found smaller dir {} ({})",
                        child.borrow().name,
                        child.borrow().size + child.borrow().child_size
                    );
                    smaller = Some(smol1)
                }
            }
        }

        return smaller;
    }
}

fn filter_val(val: usize, max_val: usize) -> usize {
    return if val <= max_val { val } else { 0 };
}

enum Command {
    Cd,
    Ls,
}

struct CommandAction {
    command: Command,
    target: Option<String>,
}

pub fn run(file_name: &str) -> Result<(), Box<dyn Error>> {
    let lines = read::lines(file_name)?;
    let root_dir = Rc::new(RefCell::new(Directory::new("/".to_string(), None)));
    lines.fold(Rc::clone(&root_dir), |cwd, line| {
        if let Ok(ip) = line {
            if let Some(cmd) = parse_command(&ip) {
                match cmd.command {
                    Command::Ls => (), // Do nothing for LS and let contents be picked up in next loop
                    Command::Cd => {
                        let target = cmd.target.unwrap();
                        return navigate_directory(Rc::clone(&cwd), target);
                    }
                }
            } else {
                parse_ls(Rc::clone(&cwd), &ip)
            }
        }
        return cwd;
    });

    let total_size = root_dir.borrow_mut().update_size();
    root_dir.borrow_mut().print(0);

    let small_dir_size = root_dir.borrow_mut().sum_subdirs_of_max_size(100000);
    println!("There are '{}' bytes in small directories", small_dir_size);

    let size_to_free = total_size - 40000000;
    let smallest_dir_to_free = root_dir
        .borrow_mut()
        .smallest_dir_larger_than(size_to_free, None)
        .unwrap();
    println!("Deleting directory '?' will free up '{}' bytes, which is more than the space needed to free '{}'", 
              smallest_dir_to_free, size_to_free);

    Ok(())
}

fn parse_command(command: &str) -> Option<CommandAction> {
    let cd_regex = Regex::new(CD_REGEX_PATTERN).unwrap();
    return if command.starts_with("$") {
        if let Some(caps) = cd_regex.captures(command) {
            Some(CommandAction {
                command: Command::Cd,
                target: Some(caps.get(1).unwrap().as_str().to_string()),
            })
        } else {
            // We only have to support two commands here so if it's not `cd` it must be `ls`
            Some(CommandAction {
                command: Command::Ls,
                target: None,
            })
        }
    } else {
        None
    };
}

fn navigate_directory(cwd: Rc<RefCell<Directory>>, cd_target: String) -> Rc<RefCell<Directory>> {
    return if cd_target.eq(UP_DIR) {
        match cwd.borrow().parent {
            None => cwd.clone(),
            Some(ref parent) => parent.clone(),
        }
    } else if cd_target.eq(ROOT_DIR) {
        // We should never get here after the first `cd` command,
        // but if we do this implementation won't work!
        cwd.clone()
    } else {
        let child = Rc::new(RefCell::new(Directory::new(cd_target, Some(cwd.clone()))));
        cwd.borrow_mut().add_child(child.clone());
        child
    };
}

fn parse_ls(cwd: Rc<RefCell<Directory>>, line: &str) {
    let file_regex = Regex::new(FILE_REGEX_PATTERN).unwrap();
    if let Some(file_caps) = file_regex.captures(line) {
        let size = file_caps.get(1).unwrap().as_str().parse::<usize>().unwrap();
        cwd.borrow_mut().add_file(size)
    }
    // Ignore directories
}
