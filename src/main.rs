
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

/*fn float_string_to_ratio(st: &String) -> Ratio<i32> {
    st.parse().expect("float to string parse error")
}

fn get_table_simple(formula: &[u8], map: & mut HashMap<String, Ratio<i32>>, scalar: Ratio<i32>)  {

    lazy_static! {

        //static ref quantity_reg: Regex = Regex::new(r"[A-Z][a-z]*[0-9]*(\.[0-9]+)?").unwrap();
        //static ref symbol_reg: Regex = Regex::new(r"[A-Z][a-z]*").unwrap();
    };

    let quantity_reg: Regex = Regex::new(r"[A-Z][a-z]*[0-9]*(\.[0-9]+)?").unwrap();
    let symbol_reg: Regex = Regex::new(r"[A-Z][a-z]*").unwrap();

    if formula.len() != 0 {

        for cap in quantity_reg.find_iter(String::from_utf8_lossy(formula).as_ref()) {
            //print!("{} ", cap.as_str());

            let symbol_number = symbol_reg.find(cap.as_str()).unwrap();

            let number = &cap.as_str().as_bytes()[symbol_number.end()..];

            let number = if number.is_empty() {
                Ratio::one()
            }
            else {
                float_string_to_ratio(&String::from_utf8_lossy(number).to_string())
            } * scalar;

            let symbol = String::from(symbol_number.as_str());

            let num = map.get_mut(symbol.as_str()).unwrap();

            *num = *num + number;


        }
        //println!()
    }

}

fn get_table(formula: &[u8], map: & mut HashMap<String, Ratio<i32>>, scalar: Ratio<i32>) {

    let mut first = 0;
    let mut last_splice = 0;

    let mut balance = 0;

    let mut ship: Option<&[u8]> = Option::None;

    let mut sub = String::new();

    for i in 0..formula.len() { //Ca(OH(Ch2Ah3)3Fh)2Ah(Fr)3
        if ship.is_some()  {
            if formula[i] >= 48 && formula[i] <= 57 || formula[i] == 46 {
                sub.push(char::from(formula[i]));
                last_splice = last_splice + 1;
            }
            else {
                get_table(ship.unwrap(), map, float_string_to_ratio(&sub) * scalar);
                sub = String::new();
                ship = None;
            }
        }

        if formula[i] == 40 {
            if balance == 0 { //If we're not already within a nested parenthesis, we are now!
                get_table_simple(&formula[last_splice..i], map, scalar);
                first = i;
            }
            balance = balance + 1;
        }
        else if formula[i] == 41 {
            balance = balance - 1;
            if balance == 0 { //If balance is zero, we have found the correct closing parenthesis
                last_splice = i+1;
                ship = Some(&formula[first+1..i]);
            }
        }

    }

    //If there is still a value to ship, i.e. a closing bracket followed by a number at the end of the string, ship it
    if ship.is_some()  {
        get_table(ship.unwrap(), map, float_string_to_ratio(&sub) * scalar);
    }

    get_table_simple(&formula[last_splice..], map, scalar);

}

fn balance_equation(equation: &String, verbose: bool) -> Result<String, String>  {
    lazy_static! {
        //static ref molecule_pattern: Regex = Regex::new(r"([A-Z][a-z]*[0-9]*(\.[0-9]+)?|\(|\)[0-9]*(\.[0-9]+)?)+").unwrap();
        //static ref term_pattern: Regex = Regex::new(r"([A-Z][a-z]*[0-9]*(\.[0-9]+)?|\(|\)[0-9]*(\.[0-9]+)?)+[\+=]?").unwrap();
        //static ref symbol_reg: Regex = Regex::new(r"[A-Z][a-z]*").unwrap();
    };

    let start = Instant::now();

    let molecule_pattern = Regex::new(r"([A-Z][a-z]*[0-9]*(\.[0-9]+)?|\(|\)[0-9]*(\.[0-9]+)?)+").unwrap();
    let term_pattern = Regex::new(r"([A-Z][a-z]*[0-9]*(\.[0-9]+)?|\(|\)[0-9]*(\.[0-9]+)?)+[\+=]?").unwrap();
    let symbol_reg: Regex = Regex::new(r"[A-Z][a-z]*").unwrap();

    println!("Regex: {:?}", start.elapsed());

    let mut symbol_table: HashMap<String, Ratio<i32>> = HashMap::new();

    // Perform an initial parse to get a set of all symbols used
    for cap in symbol_reg.find_iter(equation) {
        symbol_table.insert(String::from(cap.as_str()), Ratio::zero());
    }

    if verbose {
        print!("Symbols: ");
        for element in symbol_table.keys() {
            print!("{} ", element);
        }
        println!("\n");
    }

    //Construct an augmented matrix
    let mut matrix = Augmented::new(symbol_table.len() );

    let equals_position = match equation.find("=")
    {
        Some(s) => s,
        None => {
            return Err(String::from("Invalid chemical equation - no equals sign found"));
        }
    };



    let is_first_molecule = true;


    for mat in molecule_pattern.find_iter(equation){
        if verbose {
            print!("Molecule: {}\n\t", mat.as_str());
        }
        let mut table = symbol_table.clone();

        get_table(mat.as_str().as_bytes(), &mut table,
      if mat.end() <= equals_position {
                Ratio::one()
            }
            else {
                -Ratio::one()
            }
        );

        if verbose {
            for (key, value) in table.iter() {
                print!("( {} {} )", key, value);
            }
            println!("\n");
        }

        let mut column: Vec<Ratio<i32>> = Vec::new();

        for value in table.values() {
            column.push(*value)
        }

        if is_first_molecule {
            column.push(Ratio::one());
        }
        else {
            column.push(Ratio::zero());
        }

        

        matrix.add_column(&column);


    }

    matrix.augment();

    println!("Construct augmented matrix: {:?}", start.elapsed());


    //matrix.convert();

    if verbose {
        println!("Augmented matrix:");
        matrix.print();
    }

    matrix.row_reduce();


    println!("Row reduction: {:?}", start.elapsed());

    if verbose {
        println!("Row reduced matrix:");
        matrix.print();
    }

    match matrix.solve() {
        Ok(r) => {


            println!("Matrix solved: {:?}", start.elapsed());

            if verbose {
                println!("Solved matrix:");
                matrix.print();
            }

            if verbose {
                print!("Coefficients: ");
                for ratio in r.iter() {
                    print!("{} ", ratio);
                }
                println!("\n");
            }
            let mut output = String::new();

            let mut ind = 0;

            for mat in term_pattern.find_iter(equation) {
                let replacement = r[ind];

                let st = if replacement == 1 {
                    String::from(mat.as_str())
                }
                else {
                    format!("{}{}", replacement, mat.as_str())
                };

                output.push_str(st.as_str());


                ind += 1;
            }

            Ok(output)
        }
        Err(e) => Err(e)
    }
}
*/
fn parse_group<'a>(group: & 'a [u8], map: & mut HashMap<& 'a [u8], Ratio<i32>>, scalar: Ratio<i32>) -> Option<u8> {

    for token in TokenIterator::new(group) {
        match token {
            TokenType::Symbol(element, quantity) => {

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
            TokenType::Group(group, quantity) => {
                let mut err = parse_group(group, map, scalar * quantity);
                if err.is_some() {
                    return err.take();
                }
            },
            TokenType::Garbage(c) => {
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
            print!("({} {}) ", std::str::from_utf8(element).unwrap(), quantity);
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

    //Perform the initial run looking for ions and symbols
    for token in TokenIterator::new(equation_asbytes) {
        match token {
            TokenType::Symbol(s, _) => {
                master_table.insert(s, Ratio::zero());
            },
            TokenType::Group(group, _) => {
                parse_group(group, & mut master_table, Ratio::zero());
            }
            _ => {}
        }
    }

    if verbose {
        println!("Symbol table");
        for element in master_table.keys() {
            print!("{} ", std::str::from_utf8(element).unwrap());
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
            TokenType::Symbol(element, quantity) => {

                println!("Element: {} {}", std::str::from_utf8(element).unwrap(), quantity);
                let num_ref = symbol_table.get_mut(element).unwrap();
                *num_ref += quantity * sign;

            },
            TokenType::Group(group, quantity) => {
                let err = parse_group(group, & mut symbol_table, quantity * sign);
                if err.is_some() {
                    return Err(format!("Invalid token '{}' in formula {}", std::str::from_utf8(&[err.unwrap()]).unwrap(), equation));
                }
            },
            TokenType::Garbage(c) => {
                return Err(format!("Invalid token '{}' in formula {}", std::str::from_utf8(&[c]).unwrap(), equation));
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
        Ok(r) => {
            if verbose {

                println!("Solved matrix");
                matrix.print();

                println!("Solution vector");
                for quantity in r.iter() {
                    print!("{} ", quantity);
                }
                println!("\n");
            }

            let mut result = String::with_capacity(equation.len());
            /*let equation = String::from(equation);
            let mut start = true;
            let mut coeff_index = 0;


            for ch in equation.chars() {

                if ch == '\n' || ch == ' ' || ch == '\t' || ch == '\r' {
                    continue;
                }

                if start {
                    if r[coeff_index] != 1 {
                        result = format!("{}{}", result, r[coeff_index]);
                    }

                    result.push(ch);

                    start = false;

                    coeff_index += 1;
                }
                else {
                    result.push(ch);

                    if ch == '+' || ch == '=' {
                        if r[coeff_index] != 1 {
                            result = format!("{}{}", result, r[coeff_index]);
                        }

                        coeff_index += 1;
                    }
                }


            }*/

            Ok(result)
        },
        Err(e) => {
            Err(e)
        }
    }

}

fn main() {

    let matches = App::new("Chemical Equation Balancer")
        .version("0.2.0")
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

    let balance = better_solve_equation(
        matches.value_of("equation").unwrap(),
        matches.is_present("verbose"));


    match balance {
        Ok(s) => {
            println!("Equation: {}", s);
            if matches.is_present("duration") {
                println!("Elapsed: {:?}", start.elapsed());
            }
        }
        ,
        Err(e) => println!("Cannot solve equation. {}", e)
    };


}