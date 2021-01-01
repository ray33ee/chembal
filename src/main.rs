
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

fn parse_group<'a>(group: & 'a [u8], map: & mut HashMap<& 'a [u8], Ratio<i32>>, scalar: Ratio<i32>) -> Option<& 'a [u8]> {

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
                return Some(c)
            },
            _ => {

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

fn better_solve_equation(equation: &str, verbose: bool) -> Result<String, String> {

    let mut master_table = HashMap::<&[u8], Ratio<i32>>::new();

    let equation_asbytes = equation.as_bytes();

    let mut equals_count = 0;
    let mut equals_index = 0;

    for ch in equation_asbytes {
        if *ch > 127 {
            return Err(format!("Invalid character detected. Please only use valid Ascii characters (Unicode is not supported)"));
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
        return Err(format!("Formula must have exactly one equals"));
    }

    if equals_index == 0 || equals_index == equation_asbytes.len() - 1 {
        return Err(format!("Formula must have at least one reactant and one product"));
    }

    //Perform the initial run looking for ions and symbols
    for token in TokenIterator::new(equation_asbytes) {
        match token {
            TokenType::Symbol(s, _, _) => {
                master_table.insert(s, Ratio::zero());
            },
            TokenType::Group(group, _, _) => {
                parse_group(group, & mut master_table, Ratio::zero());
            }
            _ => {}
        }
    }

    if equals_count != 1 {
        return Err(format!("Equation must contain exactly one equals sign ('=')"));
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

    for token in TokenIterator::new(equation_asbytes) {


        match token {
            TokenType::Symbol(element, _, quantity) => {

                let num_ref = symbol_table.get_mut(element).unwrap();
                *num_ref += quantity * sign;

            },
            TokenType::Group(group, _, quantity) => {
                let err = parse_group(group, & mut symbol_table, quantity * sign);
                if err.is_some() {
                    return Err(format!("Invalid token '{}' in formula {}", unsafe { std::str::from_utf8_unchecked(err.unwrap()) }, equation));
                }
            },
            TokenType::Invalid(c) => {

                return Err(format!("Invalid token '{}' in formula {}", unsafe { std::str::from_utf8_unchecked(c) }, equation));

            },
            TokenType::Separator(sep) => {

                send_column(& mut symbol_table, & mut column, & mut matrix, master_table.clone(),verbose);

                //If the separator is an equals, flip the sign
                if sep == 61 {
                    sign = -Ratio::<i32>::one();
                }
            },
            TokenType::Whitespace => {
                // println!("Whitespace");
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
                for quantity in solution.iter() {
                    print!("{} ", quantity);
                }
                println!("\n");
            }

            let mut result = String::with_capacity(equation.len());

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

fn main() {

    let matches = App::new("Chemical Equation Balancer")
        .version("0.2.1")
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
            .help("Shows calculation duration"))
        .get_matches();

    let start = Instant::now();

    let verbose = matches.is_present("verbose");

    let balance = better_solve_equation(
        matches.value_of("equation").unwrap(),
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