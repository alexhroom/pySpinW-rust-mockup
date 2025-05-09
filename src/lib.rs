use pyo3::prelude::*;
use nalgebra::*;

/// Temporary description of the coupling between atoms.
#[pyclass]
struct Coupling {
    index1: usize,
    index2: usize,
    matrix: Matrix3x3<f64>,
    inter_site_vector: Vec<f64> 
}

/// Run the main calculation step for a spinwave calculation.
#[pyfunction]
fn pyo3_spinwave_calculation(rotations: , magnitudes:, q_vectors:, couplings: ) -> PyResult<Vec<f64>> {
    todo!() 
}

#[allow(non_snake_case)]
fn spinwave_calculation(rotations: Vec<Matrix3<f64>>, magnitudes: Vec<f64>, q_vectors: Vec<Vector3<f64>>, couplings: Vec<Coupling>) -> Vec<Vec<f64>> {
    let n_sites = rotations.len();

    // in the notation of Petit (2011)
    // eta[i] is the direction of the i'th moment in Cartesian coordinates
    let z = Complex {re: _get_rotation_component(&rotations, 0), im: _get_rotation_component(&rotations, 1)};
    let eta = _get_rotation_component(&rotations, 2);

    // make spin coefficients array
    // so spin_coefficients[i, j] = sqrt(S_i S_j) / 2
    let root_mags = DVector::<f64>::from_iterator(n_sites, magnitudes.into_iter().map(|x| (0.5*x).sqrt()));
    let spin_coefficients = (root_mags * root_mags.transpose()).transpose();

    let mut C = DMatrix::<f64>::zeros(n_sites, n_sites);
    for coupling in couplings {
        let j = coupling.index2;
        for l in 0..n_sites {
            C[(j, j)] += spin_coefficients[(l, l)] * eta[]
        }
    };

    q_vectors.into_iter().map(|q| _spinwave_single_q(q, &C, &z, &spin_coefficients, &couplings)).collect() 
} 

#[allow(non_snake_case)]
fn _spinwave_single_q(q: Vector3<f64>, C: &DMatrix<f64>, z: &Complex<MatrixXx3<f64>>, spin_coefficients: &DMatrix<f64>, couplings: &Vec<Coupling>) -> Vec<f64>{
    todo!()
}

fn _get_rotation_component(rotations: &Vec<Matrix3<f64>>, index: usize) -> MatrixXx3<f64> {

}

/// A Python module implemented in Rust.
#[pymodule]
fn pyspinw(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(pyo3_spinwave_calculation, m)?)?;
    Ok(())
}
