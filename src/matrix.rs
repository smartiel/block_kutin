// Masks for row and column operations
use std::fmt;
const COL_MASK_2: [u64; 2] = [3, 12];
const COL_MASK_4: [u64; 4] = [15, 240, 3840, 61440];
const COL_MASK_6: [u64; 6] = [63, 4032, 258048, 16515072, 1056964608, 67645734912];
const COL_MASK_8: [u64; 8] = [
    255,
    65280,
    16711680,
    4278190080,
    1095216660480,
    280375465082880,
    71776119061217280,
    18374686479671623680,
];
const ROW_MASK_2: [u64; 2] = [5, 10];
const ROW_MASK_4: [u64; 4] = [4369, 8738, 17476, 34952];
const ROW_MASK_6: [u64; 6] = [
    1090785345,
    2181570690,
    4363141380,
    8726282760,
    17452565520,
    34905131040,
];
const ROW_MASK_8: [u64; 8] = [
    72340172838076673,
    144680345676153346,
    289360691352306692,
    578721382704613384,
    1157442765409226768,
    2314885530818453536,
    4629771061636907072,
    9259542123273814144,
];
#[derive(Clone)]
pub struct MatrixF2<'a> {
    pub n: &'a usize,
    pub m: &'a usize,
    // This can store a full matrix description up to 8x8 which is enough for us
    // data is store in column major (for perf reasons)
    pub data: u64,
}

impl<'a> MatrixF2<'a> {
    /// Reconstruct a matrix from a cannonical representation
    pub fn from_representation(n: &'a usize, m: &'a usize, data: &u64) -> MatrixF2<'a> {
        return Self {
            n,
            m,
            data: data.clone(),
        };
    }
    /// Checks if the matrix satisfy the target property for the step 2 of the algorithm
    /// (i.e. its bottom right block is 0)
    pub fn check_step_2(&self) -> bool {
        let p = self.n / 2;
        for i in p..*self.n {
            for j in p..*self.n {
                if self.get(i, j) {
                    return false;
                }
            }
        }
        return true;
    }
    /// Construct a fresh identity matrix
    pub fn identity(n: &'a usize, m: &'a usize) -> MatrixF2<'a> {
        assert!([2, 4, 6, 8].contains(n));
        assert!([2, 4, 6, 8].contains(m));
        let mut data: u64 = 0;
        for i in 0..std::cmp::min(*n, *m) {
            data |= 1 << (i * n + i);
        }
        return Self { n, m, data };
    }
    /// Construct a fresh identity matrix
    pub fn identity_half(n: &'a usize, m: &'a usize) -> MatrixF2<'a> {
        assert!([2, 4, 6, 8].contains(n));
        assert!([2, 4, 6, 8].contains(m));
        let mut data: u64 = 0;
        for i in 0..(std::cmp::min(*n, *m) / 2) {
            data |= 1 << (i * n + i);
        }
        return Self { n, m, data };
    }
    /// Build an all-0 matrix
    pub fn zero(n: &'a usize, m: &'a usize) -> MatrixF2<'a> {
        assert!([2, 4, 6, 8].contains(n));
        assert!([2, 4, 6, 8].contains(m));
        return Self { n, m, data: 0 };
    }
    /// Performs a row operation
    pub fn row_op(&mut self, i: usize, j: usize) {
        match self.n {
            2 => {
                if i < j {
                    self.data ^= (self.data & ROW_MASK_2[i]) << (j - i);
                } else {
                    self.data ^= (self.data & ROW_MASK_2[i]) >> (i - j);
                }
            }
            4 => {
                if i < j {
                    self.data ^= (self.data & ROW_MASK_4[i]) << (j - i);
                } else {
                    self.data ^= (self.data & ROW_MASK_4[i]) >> (i - j);
                }
            }
            6 => {
                if i < j {
                    self.data ^= (self.data & ROW_MASK_6[i]) << (j - i);
                } else {
                    self.data ^= (self.data & ROW_MASK_6[i]) >> (i - j);
                }
            }
            8 => {
                if i < j {
                    self.data ^= (self.data & ROW_MASK_8[i]) << (j - i);
                } else {
                    self.data ^= (self.data & ROW_MASK_8[i]) >> (i - j);
                }
            }
            _ => panic!(),
        }
    }
    /// Performs a column operation
    pub fn col_op(&mut self, i: usize, j: usize) {
        match self.m {
            2 => {
                if i < j {
                    self.data ^= (self.data & COL_MASK_2[i]) << (self.n * (j - i));
                } else {
                    self.data ^= (self.data & COL_MASK_2[i]) >> (self.n * (i - j));
                }
            }
            4 => {
                if i < j {
                    self.data ^= (self.data & COL_MASK_4[i]) << (self.n * (j - i));
                } else {
                    self.data ^= (self.data & COL_MASK_4[i]) >> (self.n * (i - j));
                }
            }
            6 => {
                if i < j {
                    self.data ^= (self.data & COL_MASK_6[i]) << (self.n * (j - i));
                } else {
                    self.data ^= (self.data & COL_MASK_6[i]) >> (self.n * (i - j));
                }
            }
            8 => {
                if i < j {
                    self.data ^= (self.data & COL_MASK_8[i]) << (self.n * (j - i));
                } else {
                    self.data ^= (self.data & COL_MASK_8[i]) >> (self.n * (i - j));
                }
            }
            _ => panic!("Size {} is not supported!", self.m),
        }
    }
    /// Updates the matrix by left application of a list of CNOT gates (possibly in reverse order)
    pub fn apply_operator(&mut self, operator: &Vec<(usize, usize)>) {
        for (a, b) in operator {
            self.row_op(*a, *b);
        }
    }
    /// Returns a single entry of the matrix
    pub fn get(&self, i: usize, j: usize) -> bool {
        return (self.data >> (j * self.n + i)) & 1 != 0;
    }
    pub fn set(&mut self, i: usize, j: usize) {
        self.data |= 1 << (j * self.n + i);
    }
    /// Clones the matrix and shift its data by m/2 columns (basically fetches the second half of the matrix)
    pub fn clone_shift(&self) -> MatrixF2<'a> {
        return Self {
            n: self.n,
            m: self.m,
            data: self.data >> ((self.m / 2) * self.n),
        };
    }
    /// Computes a cannonical representation of the matrix as a u64 (for step 2)
    /// This might modify the matrix.
    pub fn cannonical_representation(&mut self) -> u64 {
        self.half_column_echelon_form();

        let mut shifted_self = self.clone_shift();

        shifted_self.half_column_echelon_form();

        let new_data = (self.data & ((1 << ((self.m / 2) * self.n)) - 1))
            | (shifted_self.data << ((self.m / 2) * self.n));

        self.data = new_data;
        if *self.n < 8 || *self.m < 8 {
            return new_data % (1 << (self.n * self.m));
        }
        return new_data;
    }
    /// Computes a cannonical representation of the matrix as a u64 (for step 1)
    /// This might modify the matrix.
    pub fn cannonical_representation_step1(&mut self) -> u64 {
        self.half_column_echelon_form();
        return self.data;
    }
    /// Computes a cannonical representation of the matrix as a u64 (for step 3)
    /// In that case the representant is trivial (we return the matrix data)
    pub fn cannonical_representation_step3(&self) -> u64 {
        return self.data;
    }

    /// Put the first n/2 columns of the matrix in column echelon form
    fn half_column_echelon_form(&mut self) {
        let p = self.n / 2;
        let mut h = 0;
        let mut k = 0;
        while h < *self.n && k < p {
            let mut pivot: Option<usize> = None;
            for i in k..p {
                if self.get(h, i) {
                    pivot = Some(i);
                    break;
                }
            }

            if let Some(pivot) = pivot {
                if pivot != k {
                    self.col_op(pivot, k);
                }
                for i in 0..(self.m / 2) {
                    if self.get(h, i) && i != k {
                        self.col_op(k, i);
                    }
                }
                k += 1;
            }
            h += 1;
        }
    }
}

impl fmt::Display for MatrixF2<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in 0..*self.n {
            for j in 0..*self.m {
                write!(f, "{} ", (self.data >> (j * self.n + i)) & 1)?;
            }
            write!(f, "\n")?;
        }
        return fmt::Result::Ok(());
    }
}

#[cfg(test)]
mod test_matrix {
    use super::MatrixF2;

    #[test]
    fn test_id() {
        let mat = MatrixF2::identity(&4, &4);
        for i in 0..4 {
            assert!(mat.get(i, i));
        }
        let mat = MatrixF2::identity(&4, &8);
        println!("{mat}");
        for i in 0..4 {
            assert!(mat.get(i, i));
        }
    }
    #[test]
    fn test_row_op() {
        let mut mat = MatrixF2::identity(&6, &6);
        mat.row_op(0, 1);
        assert!(mat.get(1, 0));
        mat.row_op(0, 2);
        assert!(mat.get(2, 0));
        mat.row_op(1, 2);
        assert!(!mat.get(2, 0));
        assert!(mat.get(2, 1));
    }
    #[test]
    fn test_col_op() {
        let mut mat = MatrixF2::identity(&6, &6);
        mat.col_op(0, 1);
        assert!(mat.get(0, 1));
        mat.col_op(0, 2);
        assert!(mat.get(0, 2));
        mat.col_op(1, 2);
        assert!(!mat.get(0, 2));
        assert!(mat.get(1, 2));
    }
    #[test]
    fn test_repr() {
        let mut mat = MatrixF2::identity(&8, &8);
        mat.row_op(0, 1);
        mat.row_op(0, 2);
        println!("{mat}");
        println!("{}", mat.cannonical_representation());
        println!("{mat}");
    }
    #[test]
    fn test_repr2() {
        let mut mat = MatrixF2::identity(&6, &6);
        mat.row_op(0, 4);
        mat.row_op(0, 5);
        mat.row_op(3, 5);
        println!("{mat}");
        println!("{}", mat.cannonical_representation());
        println!("{mat}");
        println!(
            "{}",
            MatrixF2::from_representation(&6, &6, &mat.cannonical_representation())
        );
    }
    #[test]
    fn test_col_op_exhaust() {
        for size in [4, 6, 8] {
            for i in 0..size {
                for j in 0..size {
                    if i == j {
                        continue;
                    }
                    for k in 0..size {
                        let mut mat = MatrixF2::zero(&size, &size);
                        mat.set(k, i);
                        mat.col_op(i, j);
                        assert!(mat.get(k, j));
                        mat.col_op(j, i);
                        assert!(!mat.get(k, i));
                    }
                }
            }
        }
    }
}
