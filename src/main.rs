mod solve;

//extern crate num;

use crate::solve::matrices::Augmented;
use num_rational::Ratio;
use std::collections::HashMap;

use num_traits::identities::One;
use num_traits::identities::Zero;
use regex::Regex;

use clap::{Arg, App};

#[macro_use] extern crate lazy_static;

fn float_string_to_ratio(st: &String) -> Ratio<i32> {
    st.parse().expect("float to string parse error")
}

fn get_table_simple(formula: &[u8], map: & mut HashMap<String, Ratio<i32>>, scalar: Ratio<i32>)  {

    lazy_static! {

        static ref quantity_reg: Regex = Regex::new(r"[A-Z][a-z]*[0-9]*(\.[0-9]+)?").unwrap();
        static ref symbol_reg: Regex = Regex::new(r"[A-Z][a-z]*").unwrap();
    };

    //let quantity_reg: Regex = Regex::new(r"[A-Z][a-z]*[0-9]*(\.[0-9]+)?").unwrap();
    //let symbol_reg: Regex = Regex::new(r"[A-Z][a-z]*").unwrap();

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

    get_table_simple(&formula[last_splice..], map, scalar);

}

fn balance_equation(equation: &String, verbose: bool) -> Result<String, String>  {
    lazy_static! {
        //static ref molecule_pattern: Regex = Regex::new(r"([A-Z][a-z]*[0-9]*(\.[0-9]+)?|\(|\)[0-9]*(\.[0-9]+)?)+").unwrap();
        //static ref term_pattern: Regex = Regex::new(r"([A-Z][a-z]*[0-9]*(\.[0-9]+)?|\(|\)[0-9]*(\.[0-9]+)?)+[\+=]?").unwrap();
        //static ref symbol_reg: Regex = Regex::new(r"[A-Z][a-z]*").unwrap();
    };

    let molecule_pattern = Regex::new(r"([A-Z][a-z]*[0-9]*(\.[0-9]+)?|\(|\)[0-9]*(\.[0-9]+)?)+").unwrap();
    let term_pattern = Regex::new(r"([A-Z][a-z]*[0-9]*(\.[0-9]+)?|\(|\)[0-9]*(\.[0-9]+)?)+[\+=]?").unwrap();
    let symbol_reg: Regex = Regex::new(r"[A-Z][a-z]*").unwrap();

    let mut symbol_table: HashMap<String, Ratio<i32>> = HashMap::new();

    // Perform an initial parse to get a set of all symbols used
    for cap in symbol_reg.find_iter(equation) {
        symbol_table.insert(String::from(cap.as_str()), Ratio::zero());
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
            println!("Molecule: {}", mat.as_str());
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
            println!();
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




    if verbose {
        println!("Augmented matrix:");
        matrix.print();
    }

    matrix.row_reduce();

    if verbose {
        println!("Row reduced matrix:");
        matrix.print();
    }

    match matrix.solve() {
        Ok(r) => {

            if verbose {
                println!("Solved matrix:");
                matrix.print();
            }

            if verbose {
                print!("Coefficients: ");
                for ratio in r.iter() {
                    print!("{} ", ratio);
                }
                println!();
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

fn main() {

    //let text = "H2+O2=H2O";
    //let text = "H2SO4F3D2+NaOH=Na2SO4F2D2+H2O";
    //let text = "H2SO4+NaOH=Na2SO4+H2O";
    //let text = "P4O10+H2O=H3PO4";
    //let text = "H2O=HO2";
    //let text = "H2+O5=H3+O7";
    //let text = "H2=H3";
    //let text = "HNO3+Ag=AgNO3+NO2+H2O";

    let matches = App::new("Chemical Equation Balancer")
        .version("0.1.0")
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
            .help("Displays intermediate steps")
        )
        .get_matches();


    match balance_equation(&String::from(
        matches.value_of("equation").unwrap()),
    matches.is_present("verbose")
    ) {
        Ok(s) => println!("Equation: {}", s),
        Err(e) => println!("Cannot solve equation. {}", e)
    };

}