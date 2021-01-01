


pub mod equation_parser {


    use std::iter::Iterator;

    use num_rational::Ratio;

    use num_traits::One;

    //Converts a slice floating point number into a rational number
    fn float_string_to_ratio(st: &[u8]) -> Ratio<i32> {
        if st.len() == 0 {
            Ratio::one()
        } else {
            String::from_utf8_lossy(st).parse().expect("float to string parse error")
        }
    }

    //Get the number part of the string, with start and end index too
    fn get_num_index(tok_str: & [u8]) -> (usize, usize, Ratio<i32>) {
        let mut start = 0;
        let mut end = 0;

        //Locate first character in number
        for ch in tok_str {

            if (*ch >= 48 && *ch <= 57) || *ch == 46 {
                break;
            }

            start += 1;
        }

        //Locate final character in number
        for ch in &tok_str[start..] {

            if (*ch < 48 || *ch > 57) && *ch != 46 {
                break;
            }

            end += 1;
        }

        (start, end, float_string_to_ratio(&tok_str[start..start+end]))
    }

    // Possible tokens returned when parsing a molecule
    pub enum TokenType<'a> {
        Symbol(& 'a [u8], & 'a str, Ratio<i32>), // A symbol, followed by an optional quantity, i.e. H, Na2, Mg3, etc.
        Group(& 'a [u8], & 'a str, Ratio<i32>), // A series of tokens within brackets, with optional quantity, i.e. (OH)2, (CH3), (SO4)2, etc.
        Separator(u8), //Molecule separator, either a plus or an equals
        Whitespace, //Any amount of contiguous whitespace
        Invalid(& 'a [u8])
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

            let COMPONENT_LIST: [TokenComponent; 5] = [
                TokenComponent { //Symbol
                    _start_condition: |ch| *ch >= 65 && *ch <= 90,
                    _end_condition: |ch| (*ch < 48 || *ch > 57) && (*ch < 97 || *ch > 122) && *ch != 46,
                    _parse: |tok_str: & 'a [u8]| {
                        let (num_start, _, ratio) =  get_num_index(tok_str);

                        TokenType::Symbol(&tok_str[..num_start], unsafe { std::str::from_utf8_unchecked(&tok_str) },ratio)
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

                        let (_, num_end, num) = get_num_index(&tok_str[1..]);

                        TokenType::Group(&tok_str[num_end+1..], unsafe { std::str::from_utf8_unchecked(&tok_str) },num)
                    },
                    _ignore_end: false
                }/*,
                TokenComponent { //Electron
                    _start_condition: |ch| *ch == 101,
                    _end_condition: |ch| *ch != 101,
                    _parse: |tok_str: & 'a [u8]| {
                        println!("electron: {}", unsafe { std::str::from_utf8_unchecked(&tok_str) });
                        TokenType::Symbol("charge".as_bytes(), unsafe { std::str::from_utf8_unchecked(&tok_str) },-Ratio::one())
                    },
                    _ignore_end: false
                }*/,
                TokenComponent { //Charge
                    _start_condition: |ch| *ch == 123,
                    _end_condition: |ch| *ch == 125,
                    _parse: |tok_str: & 'a [u8]| {

                        let truncated = &tok_str[1..];

                        let (start, end, mut num) = get_num_index(truncated);

                        /*println!("start: {}, end: {}, len: {}, slice: {}", start, end, truncated.len(), unsafe { std::str::from_utf8_unchecked(truncated) });

                        if start != 0 && end != truncated.len() -1 {
                            //Number is in the middle of string, error
                            return TokenType::Invalid(truncated);
                        }*/

                        let sign_range = if start == 0 {
                            end..truncated.len()
                        }
                        else {
                            0..start
                        };

                        for ch in &truncated[sign_range] {
                            if *ch == 45 {
                                num *= -1;
                            }
                            else if *ch != 43 {
                                //Invalid character in charge, error
                                return TokenType::Invalid(truncated);
                            }

                        }

                        //println!("charge: {}", unsafe { std::str::from_utf8_unchecked(&tok_str) });


                        TokenType::Symbol("charge".as_bytes(), unsafe { std::str::from_utf8_unchecked(&tok_str) },num)
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

                //Search for a number after
                let mut num_index = index;

                for ch in &self._formula[index+1..] {
                    if (*ch < 48 || *ch > 57) && *ch != 46 {
                        break;
                    }
                    num_index += 1;
                }

                let ratio = float_string_to_ratio(&self._formula[index+1..num_index+1]);
                let slice = &self._formula[self._index+1..index];

                let str_slice = unsafe { std::str::from_utf8_unchecked(&self._formula[self._index..num_index+1]) };

                //Advance the iterator forward
                self._index = num_index+1;

                Some(TokenType::Group(slice, str_slice,ratio))
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