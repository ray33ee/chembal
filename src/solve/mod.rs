pub mod matrices
{
    use num_rational::Ratio;
    use num_traits::identities::Zero;
    use num_traits::identities::One;

    use std::collections::HashSet;

    // Represents an nxn augmented matrix
    pub struct Augmented {
        matrix: Vec<Vec<Ratio<i32>>>
    }

    impl Augmented {
        pub fn new(n: usize) -> Augmented {
            let mut mat: Vec<Vec<Ratio<i32>>> = Vec::new();

            for _i in 0..n {
                mat.push(Vec::new())
            }

            Augmented {
                matrix: mat
            }
        }

        // Add column to matrix
        pub fn add_column(&mut self, column: &Vec<Ratio<i32>>) {
            for i in 0..self.matrix.len() {
                self.matrix[i].push(column[i])
            }
        }

        // Multiply a row by scalar (use multiply trait instead)
        fn scalar(&mut self, scale: Ratio<i32>, index: usize) {
            for element in self.matrix[index].iter_mut() {
                *element = *element * scale;
            }
        }

        fn addmultiple(&mut self, destination: usize, source: usize, scalar: Ratio<i32>) {
            for i in 0..self.matrix[0].len() {
                self.matrix[destination][i] = self.matrix[destination][i] + self.matrix[source][i] * scalar;
            }
        }

        fn swap(&mut self, a: usize, b: usize) {
            if a > b {
                self.swap(b, a);
            } else {
                let row = self.matrix.remove(a);

                self.matrix.insert(b, row);

                let row = self.matrix.remove(b - 1);

                self.matrix.insert(a, row);
            }
        }

        // Add augmented column
        pub fn augment(&mut self) {
            for i in 0..self.matrix.len() {
                self.matrix[i].push(Ratio::new(0, 1));
            }
        }

        // Use the display trait instead of this function
        pub fn print(&self) {
            for row in self.matrix.iter() {
                for element in row {
                    print!("{} ", *element);
                }
                println!();
            }
            println!();
        }

        // Gaussian elimination using elementary row operations
        pub fn row_reduce(&mut self) {
            let row_count = if self.matrix.len() < self.matrix[0].len() {
                self.matrix.len()
            }
            else {
                self.matrix[0].len()-1
            };

            // Iterate over each row
            for i in 0..row_count {
                let mut pivot = self.matrix[i][i];

                // If pivot = 0, search for a row such that pivot != 0 and swap
                if pivot == Ratio::zero() {
                    let mut found = false;
                    for j in i + 1..self.matrix.len() {
                        if self.matrix[j][i] != Ratio::zero() {
                            self.swap(i, j);
                            found = true;
                            break;
                        }
                    }
                    pivot = self.matrix[i][i];
                    //println!("({})", found);
                    if found == false {
                        continue;
                    }
                }

                //self.print();
                self.scalar(Ratio::new(1, 1) / pivot, i);

                for j in 0..self.matrix.len() {
                    if i != j {
                        self.addmultiple(j, i, self.matrix[j][i] * -Ratio::one());
                    }
                }
                //self.print();
            }
        }

        pub fn solve(& mut self) -> Result<Vec<i32>, String> {
            //Assumes the matrix is row reduced

            //Basically here we compare the number of rows and columns to determine if the system of equations has
            //one unique solution, or infinitely many.


            // Now we need to identify all the independent variables, substitute dummy values (x_n = 1)
            // And substitute values in for the dependent variables

            // Iterate over each leading one, then over the row to determine independent variables

            //Iterate over each row
            //  Find the first one
            //  If none exists, continue
            // Every column containing a non-zero entry represents an independent variable that needs to be substituted

            let mut independent_set: HashSet<usize> = HashSet::new();
            let mut remove_set : Vec<usize> = Vec::new();

            for i in 0..self.matrix.len() {
                let mut leading_found = false; //Bool representing whether the leading one has been found

                for j in 0..self.matrix[0].len() {

                    let element = self.matrix[i][j];

                    if leading_found == true && element != Ratio::zero() {
                        independent_set.insert(j);
                    }

                    if element == Ratio::one() {
                        leading_found = true;
                    }
                }

                // If no leading one is found, the row is empty and should be marked for removal
                if leading_found == false {
                    remove_set.insert(0, i);
                }
            }

            // Remove empty rows
            for kill in remove_set.iter() {
                self.matrix.remove(*kill);
            }

            //At this point, if the matrix without the augment is square, we have a unique solution
            if self.matrix.len() == self.matrix[0].len() - 1 {
                Err(String::from("Trivial solution detected, impossible chemical equation"))
            }
            else {

                for index in independent_set.iter() {
                    let mut row = Vec::new();

                    for i in 0..self.matrix[0].len() {
                        if i == *index || i == self.matrix[0].len()-1 {
                            row.push(Ratio::one())
                        }
                        else {
                            row.push(Ratio::zero())
                        }
                    }

                    self.matrix.push(row);

                }

                self.row_reduce();

                let mut result : Vec<i32> = Vec::new();

                let mut lm = 1;

                for i in 0..self.matrix.len() {
                    let ratio = self.matrix[i][self.matrix[0].len() - 1];
                    lm = num::integer::lcm(lm, *ratio.denom());
                }

                for i in 0..self.matrix.len() {
                    let ratio = self.matrix[i][self.matrix[0].len() - 1] * lm;

                    result.push(*(ratio.numer()));
                }

                Ok(result)


            }



        }
    }
}