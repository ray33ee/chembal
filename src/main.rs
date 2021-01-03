
mod solve;
mod parser;

//extern crate num;

use crate::solve::matrices::Augmented;
use crate::parser::equation_parser::TokenIterator;
use crate::parser::equation_parser::TokenType;

use num_rational::Ratio;
use std::collections::HashMap;

use num_traits::identities::One;
use num_traits::identities::Zero;

use clap::{Arg, App};

use std::time::{Instant};

fn parse_group<'a>(group: & 'a [u8], map: & mut HashMap<& 'a [u8], Ratio<i32>>, scalar: Ratio<i32>) -> Option<String> {

    for token in TokenIterator::new(group) {
        match token {
            TokenType::Symbol(element, _, quantity) => {

                let num_ref_option = map.get_mut(element);

                match num_ref_option {
                    Some(num_ref) => {
                        *num_ref += quantity * scalar;
                    },
                    None => {
                        map.insert(element, quantity * scalar);
                    }
                }
            },
            TokenType::Group(group, _, quantity) => {
                let mut err = parse_group(group, map, scalar * quantity);
                if err.is_some() {
                    return err.take();
                }
            },
            TokenType::Invalid(c) => {
                return Some(format!("Invalid token '{}'", unsafe { std::str::from_utf8_unchecked(c) }));
            },
            TokenType::Error(slice, error) => {
                return Some(format!("{} ({})", error, slice));
            },
            TokenType::Separator(sep) => {
                return Some(format!("Invalid symbol ({}) found within parenthesis", sep))
                //Separator within group, error
            }
        }
    }

    None

}

fn send_column<'a>(table: & mut HashMap<&'a[u8], Ratio<i32>>,col: & mut Vec<Ratio<i32>>,
               mat: & mut Augmented,master: HashMap<&'a [u8], Ratio<i32>>,verbose: bool) {
    if verbose {
        for (element, quantity) in table.iter() {
            print!("({} {}) ", unsafe { std::str::from_utf8_unchecked(element) }, quantity);
        }
        println!("\n");
    }

    for quantity in table.values() {
        col.push(*quantity);
    }

    mat.add_column(&col);

    col.clear();

    *table = master.clone();
}

fn solve_equation(equation_asbytes: &[u8], verbose: bool) -> Result<String, String> {

    let mut master_table = HashMap::<&[u8], Ratio<i32>>::new();

    let mut equals_count = 0;
    let mut equals_index = 0;

    for ch in equation_asbytes {
        if *ch > 127 {
            return Err(String::from("Invalid character detected. Please only use valid Ascii characters (Unicode is not supported)"));
        }

        //Get the total number of equals signs
        if *ch == 61 {
            equals_count += 1;
        }

        //Get the index of the first occurence of an equals sign
        if equals_count == 0 {
            equals_index += 1;
        }
    }

    if equals_count != 1 {
        return Err(String::from("Formula must have exactly one equals"));
    }

    if equals_index == 0 || equals_index == equation_asbytes.len() - 1 {
        return Err(String::from("Formula must have at least one reactant and one product"));
    }

    //Perform the initial run looking for ions and symbols
    for token in TokenIterator::new(equation_asbytes) {
        match token {
            TokenType::Symbol(s, _, _) => {
                master_table.insert(s, Ratio::zero());
            },
            TokenType::Group(group, _, _) => {
                let err = parse_group(group, & mut master_table, Ratio::zero());
                if err.is_some() {
                    return Err(err.unwrap());
                };
            },
            TokenType::Invalid(c) => {

                return Err(format!("Invalid token '{}'", unsafe { std::str::from_utf8_unchecked(c) }));

            },
            TokenType::Error(slice, error) => {
                return Err(format!("{} ({})", error, slice));
            },
            TokenType::Separator(_) => {
            }
        }
    }

    if equals_count != 1 {
        return Err(String::from("Equation must contain exactly one equals sign ('=')"));
    }

    //if equation_asbytes.len() < {}

    if verbose {
        println!("Symbol table");
        for element in master_table.keys() {
            print!("{} ", unsafe { std::str::from_utf8_unchecked(element) });
        }
        println!("\n");
    }

    let mut symbol_table = master_table.clone();

    let mut sign = Ratio::<i32>::one();

    let mut matrix = Augmented::new(symbol_table.len() );

    let mut column = Vec::with_capacity(symbol_table.len() );

    if verbose {
        println!("Molecule matrix");
    }

    //Iterate over each token, adding columns to the matrix as we find new molecules
    for token in TokenIterator::new(equation_asbytes) {


        match token {
            TokenType::Symbol(element, _, quantity) => {
                let num_ref = symbol_table.get_mut(element).unwrap();
                *num_ref += quantity * sign;
            },
            TokenType::Group(group, _, quantity) => {
                parse_group(group, & mut symbol_table, quantity * sign);
            },
            TokenType::Invalid(_) => {
            },
            TokenType::Separator(sep) => {

                send_column(& mut symbol_table, & mut column, & mut matrix, master_table.clone(),verbose);

                //If the separator is an equals, flip the sign
                if sep == 61 {
                    sign = -Ratio::<i32>::one();
                }
            },
            TokenType::Error(_,_) => {
            }
        }

    }

    send_column(& mut symbol_table, & mut column, & mut matrix, master_table.clone(),verbose);

    matrix.augment();

    if verbose {
        println!("Augmented matrix");
        matrix.print();
    }

    matrix.row_reduce();

    if verbose {
        println!("Row reduced matrix");
        matrix.print();
    }

    match matrix.solve() {
        Ok(solution) => {
            if verbose {

                println!("Solved matrix");
                matrix.print();

                println!("Solution vector");
                for coeff in solution.iter() {
                    print!("{} ", coeff);
                }
                println!("\n");
            }

            let mut result = String::with_capacity(equation_asbytes.len());

            if solution[0] != 1 {
                result = format!("{}{}", result, solution[0])
            }

            let mut solution_index = 1;

            for token in TokenIterator::new(equation_asbytes) {



                match token {
                    TokenType::Symbol(symbol, string, _) => {

                        result.push_str(string);

                        if unsafe { std::str::from_utf8_unchecked(symbol) } == "charge" && string != "e" {
                            result.push_str("}");
                        }

                    },
                    TokenType::Group(_, string, _) => {

                        result.push_str(string);

                    },
                    TokenType::Separator(sep) => {

                        result.push_str(
                            if sep == 43 {
                                "+"
                            } else if sep == 61 {
                                "="
                            } else {
                                ""
                            }
                        );

                        if solution[solution_index] != 1 {
                            result = format!("{}{}", result, solution[solution_index])
                        }

                        solution_index += 1;

                    }
                    _ => {

                    }
                }

            }

            Ok(result)
        },
        Err(e) => {
            Err(e)
        }
    }

}

fn remove_whitespace(string: & str) -> String {
    let mut result = String::new();

    for ch in string.as_bytes() {
        if *ch != 9 && *ch != 10 && *ch != 13 && *ch != 32 {
            result.push(char::from(*ch));
        }
    }

    result
}

fn main() {

    let matches = App::new("Chemical Equation Balancer")
        .version("0.2.3")
        .author("Will Cooper")
        .about("Command line tool to balance chemical equations")
        .arg(Arg::with_name("equation")
            .short("e")
            .long("equation")
            .takes_value(true)
            .help("Chemical equation to balance")
            .required(true))
        .arg(Arg::with_name("verbose")
            .short("v")
            .long("verbose")
            .takes_value(false)
            .help("Displays intermediate steps"))
        .arg(Arg::with_name("duration")
            .short("d")
            .long("duration")
            .takes_value(false)
            .help("Displays computation time"))
        .get_matches();

    let start = Instant::now();

    let verbose = matches.is_present("verbose");

    let cleaned_equation = remove_whitespace(matches.value_of("equation").unwrap());

    let balance = solve_equation(
        cleaned_equation.as_bytes(),
        verbose);


    match balance {
        Ok(s) => {


            if verbose {
                println!("Equation: {}", s);
            }
            else {
                println!("{}", s);
            }

            if matches.is_present("duration") {
                println!("\nElapsed: {:?}", start.elapsed());
            }


        }
        ,
        Err(e) => println!("Cannot solve equation. {}", e)
    };


}