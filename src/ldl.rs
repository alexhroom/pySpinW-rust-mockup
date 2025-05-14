use nalgebra::{Const, DMatrix, DVector};
use num_complex::Complex;

type C64 = Complex<f64>;

/// Computes the LDL^T factorization.
///
/// The input matrix `p` is assumed to be Hermitian and this decomposition will only read
/// the lower-triangular part of `p`.
pub fn ldl(p: DMatrix<C64>) -> (DMatrix<C64>, DVector<C64>) {
    let n = p.ncols();

    let n_dim = p.shape_generic().1;

    let mut d = DVector::<C64>::zeros_generic(n_dim, Const::<1>);
    let mut l = DMatrix::<C64>::zeros_generic(n_dim, n_dim);

    for j in 0..n {
        let mut d_j = p[(j, j)].clone();

        if j > 0 {
            for k in 0..j {
                d_j -= d[k].clone() * l[(j, k)].clone().powi(2);
            }
        }

        d[j] = d_j;

        for i in j..n {
            let mut l_ij = p[(j, i)].clone();

            for k in 0..j {
                l_ij -= d[k].clone() * l[(j, k)].clone() * l[(i, k)].clone();
            }

            if d[j] == Complex::from(0.) {
                l[(i, j)] = Complex::from(0.) 
            } else {
                l[(i, j)] = l_ij / d[j].clone();
            }
        }
    }

    ( l, d )
}

