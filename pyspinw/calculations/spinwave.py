"""Spinwave Calculations"""

from dataclasses import dataclass
from enum import Enum

import numpy as np
from scipy.linalg import ldl

from pyspinw.checks import check_sizes


# Disable linting for bad variable names, because they should match the docs
# ruff: noqa: E741

@dataclass
class Coupling:
    """Temporary description of the coupling between atoms"""

    index1: int
    index2: int
    matrix: np.ndarray
    inter_site_vector: np.ndarray

class CalculationMethod(Enum):
    """Type of method used (for debugging purposes)"""

    CHOLESKY = 0
    LDL = 1

@dataclass
class SpinwaveResult:
    """Results from a spinwave calculation"""

    q_vectors: np.ndarray
    raw_energies: list[np.ndarray]
    method: list[CalculationMethod]


@check_sizes(rotations=('n_sites',3,3),
             magnitudes=('n_sites',),
             q_vectors=('n_q', 3))
def spinwave_calculation(rotations: np.ndarray,
                         magnitudes: np.ndarray,
                         q_vectors: np.ndarray,
                         couplings: list[Coupling]):
    """Main calculation step

    Unlike the main interface it takes indexed arrays, the meaning of the arrays is set elsewhere

    """
    n_sites = rotations.shape[0]

    z = rotations[:,:,0] + 1j*rotations[:,:,1] # n-by-3, complex
    z_conj = z.conj()
    eta = rotations[:,:,2] # n-by-3 real

    # Create the A, B, and C matrices

    # Create the matrix sqrt(S_i S_j)/2
    root_mags = np.sqrt(0.5*magnitudes) # S_i / sqrt(2)
    spin_coefficients = root_mags.reshape(-1, 1) * root_mags.reshape(1, -1)

    energies = []
    methods = []
    for q in q_vectors:

        A = np.zeros((n_sites, n_sites), dtype=complex)
        B = np.zeros((n_sites, n_sites), dtype=complex)
        C = np.zeros((n_sites, n_sites), dtype=complex)

        # Add the terms up to the total spin coefficient
        for coupling in couplings:

            phase_factor = np.exp((2j*np.pi)*np.dot(q, coupling.inter_site_vector))

            i = coupling.index1
            j = coupling.index2

            A[i, j] += z[i, :] @ coupling.matrix @ z_conj[j, :] * phase_factor
            B[i, j] += z[i, :] @ coupling.matrix @ z[j, :] * phase_factor


            for l in range(n_sites):
                C[j, j] += spin_coefficients[l,l] * eta[j, :].T @ coupling.matrix @ eta[l, :]

        A *= spin_coefficients
        B *= spin_coefficients
        C *= 2
        #
        # print("A", A)
        # print("B", B)
        # print("C", C)


        # hamiltonian_matrix = np.block([[A - C, B], [B.conj().T, A.conj().T - C]])
        hamiltonian_matrix = np.block([[A - C, B], [B.conj().T, A.conj().T - C]])

        #
        # We need to enforce the bosonic commutation properties, we do this
        # by finding the 'square root' of the matrix (i.e. finding K such that KK^dagger = H)
        # and then negating the second half.
        #
        # In matrix form we do
        #
        #     M = K^dagger C K
        #
        # where C is a diagonal matrix of length 2n, with the first n entries being 1, and the
        # remaining entries being -1. Note the adjoint is on the other side to the definition in
        # of the decomposition
        #
        # We can also do this via an LDL decomposition, but the method is very slightly different


        try:
            sqrt_hamiltonian = np.linalg.cholesky(hamiltonian_matrix)

            sqrt_hamiltonian_with_commutation = sqrt_hamiltonian.copy()
            sqrt_hamiltonian_with_commutation[:, n_sites:] *= -1

            to_diagonalise = np.conj(sqrt_hamiltonian_with_commutation).T @ sqrt_hamiltonian

            methods.append(CalculationMethod.CHOLESKY)

        except np.linalg.LinAlgError: # Catch postive definiteness errors


            # l, d, perm = ldl(hamiltonian_matrix) # To LDL^\dagger (i.e. adjoint on right)
            l, d, _ = ldl(hamiltonian_matrix) # To LDL^\dagger (i.e. adjoint on right)
            sqrt_hamiltonian = l @ np.sqrt(d)

            # TODO: Check for actual diagonal (could potentially contain non-diagonal 2x2 blocks)

            sqrt_hamiltonian_with_commutation = sqrt_hamiltonian.copy()
            sqrt_hamiltonian_with_commutation[:, n_sites:] *= -1

            to_diagonalise = np.conj(sqrt_hamiltonian_with_commutation).T @ sqrt_hamiltonian

            methods.append(CalculationMethod.LDL)

        eig_res = np.linalg.eig(to_diagonalise)

        energies.append(eig_res.eigenvalues) # These are currently the square energies

    # pylint: enable=C0103

    return SpinwaveResult(
                q_vectors=q_vectors,
                raw_energies=energies,
                method=methods)
