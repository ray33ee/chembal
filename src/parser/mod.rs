


pub mod equation_parser {

    use std::iter::Iterator;

    use num_rational::{Ratio, Rational32};

    use num_traits::One;

    //Converts a string floating point number into a rational number
    fn float_string_to_ratio(st: &[u8]) -> Result<Rational32, String> {
        if st.len() == 0 {
            Ok(Ratio::one())
        } else {

            //We try and convert using the native parse() function, this will only work if the number is an integer
            match String::from_utf8_lossy(st).parse() {
                Ok(number) => {
                    Ok(number)
                },
                Err(_) => {

                    match String::from_utf8_lossy(st).parse::<f64>() {
                        Ok(decimal) => {

                            let mut decimal_point_index = 0;

                            // Get the index of the floating point
                            for ch in st {
                                if *ch == '.' as u8 {
                                    break;
                                }
                                decimal_point_index += 1;
                            }

                            let exponent = st.len() - decimal_point_index - 1;

                            let denom = 10i32.pow(exponent as u32);

                            let numer: i32 = (denom as f64 * decimal) as i32;

                            //println!("Numer: {}, Demon: {}", numer, denom);

                            Ok(Ratio::new(numer, denom))
                        },
                        Err(_) => {
                            //Invalid number
                            Err(format!("'{}' is not a valid number", unsafe { std::str::from_utf8_unchecked(&st) }))
                        }
                    }
                }
            }
        }
    }

    //Get the number part of the string, with start and end index too
    fn get_num_index(tok_str: & [u8]) -> Result<(usize, usize, Ratio<i32>), String> {
        let mut start = 0;
        let mut end = 0;

        //Locate first character in number
        for ch in tok_str {

            if (*ch >= 48 && *ch <= 57) || *ch == '.' as u8 {
                break;
            }

            start += 1;
        }

        //Locate final character in number
        for ch in &tok_str[start..] {

            if (*ch < 48 || *ch > 57) && *ch != '.' as u8 {
                break;
            }

            end += 1;
        }

        match float_string_to_ratio(&tok_str[start..start+end]) {
            Ok(ratio) => Ok((start, end, ratio)),
            Err(error) => Err(error)
        }


    }

    // Possible tokens returned when parsing a molecule
    pub enum TokenType<'a> {
        Symbol(& 'a [u8], & 'a str, Ratio<i32>), // A symbol, followed by an optional quantity, i.e. H, Na2, Mg3, etc.
        Group(& 'a [u8], & 'a str, Ratio<i32>), // A series of tokens within brackets, with optional quantity, i.e. (OH)2, (CH3), (SO4)2, etc.
        Separator(u8), //Molecule separator, either a plus or an equals
        Whitespace, //Any amount of contiguous whitespace
        Invalid(& 'a [u8]),
        Error(& 'a str, String)
    }

    //Contains information for identifying and parsing tokens
    struct TokenComponent<'a> {
        _start_condition: fn(&u8) -> bool,
        _end_condition: fn(&u8) -> bool,
        _parse: fn(& 'a [u8]) -> TokenType<'a>,
        _ignore_end: bool //tells the parser whether or not to ignore the terminating character
    }

    pub struct TokenIterator<'a> {
        _formula: & 'a [u8],
        _index: usize
    }

    impl<'a> TokenIterator<'a> {
        pub fn new(st: & 'a [u8]) -> Self {
            TokenIterator {
                _formula: st,
                _index: 0
            }
        }
    }

    impl<'a> Iterator for TokenIterator<'a> {

        type Item = TokenType<'a>;

        fn next(&mut self) -> Option<Self::Item> {

            if self._index == self._formula.len() {
                return None;
            }

            #[allow(non_snake_case)]
            let COMPONENT_LIST: [TokenComponent; 5] = [
                TokenComponent { //Symbol
                    _start_condition: |ch| *ch >= 65 && *ch <= 90,
                    _end_condition: |ch| (*ch < 48 || *ch > 57) && (*ch < 97 || *ch > 122) && *ch != '.' as u8,
                    _parse: |tok_str: & 'a [u8]| {

                        match get_num_index(tok_str) {
                            Ok((num_start, _, ratio)) =>
                                TokenType::Symbol(&tok_str[..num_start], unsafe { std::str::from_utf8_unchecked(&tok_str) },ratio),
                            Err(error) =>
                                TokenType::Error(unsafe { std::str::from_utf8_unchecked(&tok_str) }, error)
                        }
                    },
                    _ignore_end: false
                },
                TokenComponent { //Whitespace
                    _start_condition: |ch| *ch == 9 || *ch == 10 || *ch == 13 || *ch == 32,
                    _end_condition: |ch| *ch != 9 && *ch != 10 && *ch != 13 && *ch != 32,
                    _parse: |_tok_str: & 'a [u8]| {
                        TokenType::Whitespace
                    },
                    _ignore_end: false
                },
                TokenComponent { //Separator
                    _start_condition: |ch| *ch == 43 || *ch == 61,
                    _end_condition: |ch| *ch != 43 && *ch != 61,
                    _parse: |tok_str: & 'a [u8]| {
                        TokenType::Separator(tok_str[0])
                    },
                    _ignore_end: false
                },
                TokenComponent { //Hydrate
                    _start_condition: |ch| *ch == 42,
                    _end_condition: |ch| *ch == 43 || *ch == 61,
                    _parse: |tok_str: & 'a [u8]| {

                        match get_num_index(&tok_str[1..]) {
                            Ok((_, num_end, ratio)) =>
                                TokenType::Group(&tok_str[num_end+1..], unsafe { std::str::from_utf8_unchecked(&tok_str) },ratio),
                            Err(error) =>
                                TokenType::Error(unsafe { std::str::from_utf8_unchecked(&tok_str) }, error)
                        }


                    },
                    _ignore_end: false
                },
                TokenComponent { //Charge
                    _start_condition: |ch| *ch == 123,
                    _end_condition: |ch| *ch == 125,
                    _parse: |tok_str: & 'a [u8]| {

                        let truncated = &tok_str[1..];

                        match get_num_index(truncated) {
                            Ok((start, end, mut ratio)) => {

                                //println!("start: {}, end: {}, len: {}, slice: {}", start, end, truncated.len(), unsafe { std::str::from_utf8_unchecked(truncated) });

                                if start != 0 && end != truncated.len() -1 {
                                    //Number is in the middle of string, error
                                    return TokenType::Invalid(truncated);
                                }

                                let sign_range = if start == 0 {
                                    end..truncated.len()
                                }
                                else {
                                    0..start
                                };

                                for ch in &truncated[sign_range] {
                                    if *ch == 45 {
                                        ratio *= -1;
                                    }
                                        else if *ch != 43 {
                                            //Invalid character in charge, error
                                            return TokenType::Invalid(truncated);
                                        }

                                }

                                //println!("charge: {}", unsafe { std::str::from_utf8_unchecked(&tok_str) });


                                TokenType::Symbol("charge".as_bytes(), unsafe { std::str::from_utf8_unchecked(&tok_str) },ratio)
                            }
                            Err(error) =>
                                TokenType::Error(unsafe { std::str::from_utf8_unchecked(&tok_str) }, error)
                        }




                    },
                    _ignore_end: true
                }
            ];

            let first_char = self._formula[self._index];

            let mut index = self._index+1;

            if first_char == 40 {
                //If the first character is an open parenthesis, we search for the closing parenthesis and the optional quantity
                let mut matching = 1;

                //Search for the closing parenthesis
                for ch in &self._formula[index..] {
                    if *ch == 40 {
                        matching += 1;
                    }
                    else if *ch == 41 {
                        matching -= 1;
                    }
                    if matching == 0 {
                        break
                    }
                    index += 1;
                }

                //println!("Par index: {}", index);

                match get_num_index(&self._formula[index+1..]) {
                    Ok((_, last, ratio)) => {
                        let slice = &self._formula[self._index+1..index];

                        let str_slice = unsafe { std::str::from_utf8_unchecked(&self._formula[self._index..index+last+1]) };

                        //Advance the iterator forward
                        self._index = index+last+1;

                        Some(TokenType::Group(slice, str_slice,ratio))
                    },
                    Err(error) => {
                        Some(TokenType::Error(unsafe { std::str::from_utf8_unchecked(&self._formula[index+1..]) }, error))
                    }
                }


            }
            else {
                // If its not, we look to COMPONENT_LIST to identify the token
                let mut result = None;

                let mut matched = false;

                for component in COMPONENT_LIST.iter() {
                    if (component._start_condition)(&first_char) { // If the component matches the first character
                        matched = true;
                        //Keep searching until the stop condition is satisfied
                        for ch in &self._formula[index..] {
                            if (component._end_condition)(ch) {
                                break
                            }
                            index += 1
                        }

                        //println!("Index: {}", std::str::from_utf8(&self._formula[self._index..index]).unwrap());

                        result = Some((component._parse)(&self._formula[self._index..index]));

                        //If the component has ignore end set, increment index to push past ignroed character
                        if component._ignore_end {
                            index += 1;
                        }

                        break;
                    }
                }

                //println!("Token: {}", std::str::from_utf8(&self._formula[self._index..index]).unwrap());

                //If we haven't identified the token, return it as an invalid token
                let result = if matched {
                    result
                } else {
                    Some(TokenType::Invalid(&self._formula[self._index..self._index+1]))
                };

                self._index = index;

                result

            }






        }

    }

}