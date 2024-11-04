//! Implementations of common linear algebra operations. For this scope of the
//! assignment only contains operations optimized for modulo-2 arithmetic,
//! implemented with the [ZMod2] struct.

use crate::basic_types::ZMod2;
use crate::basic_types::*;

use na::base::{Matrix, DMatrix};

/// Computes the reduced row echelon form modulo 2 in place for a [Matrix] of
/// which the elements are of type [ZMod2].
/// 
/// The reduced row echelon form of the matrix is computed by an implementation
/// of the Gauss-Jordan elimination algorithm adapted to modulo 2 arithmetic.
/// 
/// # Arguments
/// * `m` - Used as the input matrix and the resulting output matrix as the
///         computation is done in place.
/// 
/// # Returns
/// A vector containing the indices of the columns, in the resulting matrix, that
/// do not contain a pivot.
/// 
/// # Arguments
/// * `m` - Used as the input matrix and the resulting output matrix as the
///         computation is done in place.
/// # Note
/// A more memory and computationallity efficient could be acchieved as follows:
/// * Group the elements of each row into groups of size W, where W is the number
///   of bits in a word for the specific compute arcitecture.
/// * Map each element inside a group to the corresponding bit in a word.
/// * This way we can represent each row by words storing the values of each
///   group. We additionally need to keep track on how many bits are used in the
///   right most or left most (depending on your choice of implementation) word
///   for each row.
///       
/// It is easy to see that this way the memory complexity is reduced by a factor W.
/// Since the computational complexity of this function is linear in the
/// computational complexity of elementary row operations, and the above modification
/// reduces this complexity by a factor W, the overall computational complexity of
/// the algorithm is also reduced by a factor W. Since on most modern CPU
/// architectures W = 64, such a modification will cause a significant performance
/// gain. To keep things simple for the implementation made for this assignment, this
/// modification is not applied in the code implented here.
pub fn to_reduced_row_echelon <R,C,S> (m: &mut Matrix<ZMod2,R,C,S>) -> Vec<usize>
where R: na::Dim,
      C: na::Dim,
      S: na::RawStorageMut<ZMod2, R, C>
{
    let mut columns_wo_pivot = Vec::new();
    let mut columns_with_pivot = Vec::new();
    let mut c = 0;

    ////// Gaussian elimination: //////
    // Invariants:
    // 1. The submatrix of m consisting of the first columns_with_pivot.len()
    //    rows of m are in echelon form.
    // 2. The submatrix of m consisting of the first c columns of m are in
    //    echelon form.
    // 3. For the submatrix of m consisting of the first c columns of m,
    //    columns_with_pivot contains all indices of columns that contain a
    //    pivot entry. More specifically for any valid i,
    //    columns_with_pivot[i] denotes the column index of the i'th column
    //    with a pivot entry.
    // 4. For the submatrix of m consisting of the first c columns of m,
    //    columns_wo_pivot contains all indices of columns that not contain
    //    a pivot entry. More specifically for any valid i,
    //    columns_wo_pivot[i] denotes the column index of the i'th column
    //    with no pivot entry.
    while c < m.ncols() && columns_with_pivot.len() < m.nrows() {
        let mut r = columns_with_pivot.len();

        while r < m.nrows() && m[(r,c)].is_zero(){
            r += 1;
        }
        if r < m.nrows() { // if this column contains a non-zero entry
            if r != columns_with_pivot.len() {
                m.swap_rows(r, columns_with_pivot.len());
            }
            // Note: no need to divide the row by the leading entry, because it is mod 2
            for r in columns_with_pivot.len()+1 .. m.nrows() {
                if !m[(r,c)].is_zero() {
                    // Subtract the top row from the r'th row (unfortunately I cannot
                    // come up with a nicer expression, due to the borrowing rules).
                    // Note: also here no need to multiply the top row by m[(r,c)]
                    //       before subtracting because of mod 2
                    for i in c .. m.ncols() {
                        m[(r,i)] = m[(r,i)] - m[(columns_with_pivot.len(),i)];
                    }
                }
            }
            columns_with_pivot.push(c);
        }
        else {
            columns_wo_pivot.push(c);
        }

        c += 1;
    }
    columns_wo_pivot.extend(c .. m.ncols());
    // Now all 4 of the invariants, where c is substituted by m.ncols(), hold

    ////// The "Jordan" part of Gauss-Jordan elimination: //////
    for col_with_pivot_idx in (0 .. columns_with_pivot.len()).rev() {
        // Note: columns_with_pivot[col_with_pivot_idx] is the col_with_pivot_idx'th
        //       column of m that has a pivot entry. Observe that col_with_pivot_idx
        //       is also the index of the row of m corresponding to that same pivot
        //       entry
        c = columns_with_pivot[col_with_pivot_idx];
        for r2 in 0..col_with_pivot_idx {
            if !m[(r2,c)].is_zero() {
                for c2 in c .. m.ncols() {
                    m[(r2,c2)] = m[(r2,c2)] - m[(col_with_pivot_idx,c2)];
                }
            }
        }
    }
    columns_wo_pivot
}


/// Computes the modulo-2 null space of a [Matrix] of which the elements are of type
/// [ZMod2].
/// 
/// # Arguments
/// * `m` - The matrix to compute the null space of. After calling this function this
///         matrix will be modified.
///
/// # Returns
/// A matrix of which the columns form basis for the null space of matrix `m`.
/// As a consequence the number of rows of the resulting matrix equals the number of
/// columns of `m`, and the number of columns equals the nullity of m.
pub fn null_space<R,C,S> (m: &mut Matrix<ZMod2,R,C,S>) -> DMatrix<ZMod2>
where R: na::Dim,
      C: na::Dim,
      S: na::RawStorageMut<ZMod2, R, C>
{
    let columns_wo_pivot = to_reduced_row_echelon(m);
    let mut null_space_matrix = DMatrix::<ZMod2>::zeros(m.ncols(), columns_wo_pivot.len());

    for c in 0..columns_wo_pivot.len() {
        let mut next_col_wo_pivot = 0;

        for r in 0..m.ncols() {
            if next_col_wo_pivot < columns_wo_pivot.len() && r == columns_wo_pivot[next_col_wo_pivot] {
                next_col_wo_pivot += 1;
                if r == columns_wo_pivot[c] {
                    null_space_matrix[(r,c)] = ZMod2::one();
                }
            }
            else if r - next_col_wo_pivot < m.nrows() {
                // Note: negation not necessary since we are doing mod 2 arithmetic 
                null_space_matrix[(r,c)] = m[(r - next_col_wo_pivot, columns_wo_pivot[c])]; 
            }
        }
    }
    null_space_matrix
}



mod tests {
    use super::*;

    #[test]
    fn test_null_space() {
        let mut m = na::Matrix3x5::<ZMod2>::new(ZMod2::one(), ZMod2::one(), ZMod2::zero(), ZMod2::one(), ZMod2::one(),
                                                                                      ZMod2::zero(), ZMod2::zero(), ZMod2::one(), ZMod2::one(), ZMod2::one(),
                                                                                      ZMod2::zero(), ZMod2::zero(),ZMod2::zero(),ZMod2::zero(),ZMod2::zero(),);
        let null_space = null_space(&mut m);
        assert!((m*null_space).iter().all(|&x| x.is_zero()))
    }
}