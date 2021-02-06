#include "export_correction.hpp"
#include <losc/correction.hpp>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

void export_correction(py::module &m)
{
    // ao_hamiltonian_correction
    m.def("ao_hamiltonian_correction", &losc::ao_hamiltonian_correction,
          R"pddoc(
    Calculate LOSC effective Hamiltonian under AO basis.

    The LOSC effective Hamiltonian is constructed with LOs fixed. The
    expression of the effective Hamiltonian is shown as Eq. S25 in the
    supporting information of the original LOSC paper
    (https://doi.org/10.1093/nsr/nwx111). This effective Hamiltonian is exact
    in the developed version of SCF-LOSC (self-consistent LOSC). See reference
    (J. Phys. Chem. Lett. 2020, 11, 23, 10269-10277) for more details about
    how to perform reliable SCF-LOSC calculations.

    Parameters
    ----------
    S: np.ndarray [nbasis, nbasis]
        AO overlap matrix.
    C_lo: np.ndarray [nbasis, nbasis]
        LO coefficient matrix under AO basis.
    Curvature: np.ndarray [nbasis, nbasis]
        LOSC curvature matrix with dimension [nlo, nlo].
    LocalOcc: np.ndarray [nbasis, nbasis]
        LOSC local occupation matrix.

    Returns
    -------
    out: np.ndarray [nbasis, nbasis]
        The LOSC effective Hamiltonian under AOs.

    See Also
    --------
    CurvatureV1.kappa(), CurvatureV2.kappa(), ...: Return the curvature matrix.
    local_occupation(): Return the local occupation matrix.
    )pddoc");

    // energy_correction
    m.def("energy_correction", &losc::energy_correction, R"pddoc(
    Calculate the total energy correction from LOSC.

    This is just the energy correction from LOSC, NOT the total energy of
    LOSC-DFA. Total energy of LOSC-DFA is: E_losc_dfa = E_dfa + E_losc.

    Parameters
    ----------
    Curvature: np.ndarray [nlo, nlo]
        The LOSC curvature matrix.
    LocalOcc: np.ndarray [nlo, nlo]
        The LOSC local occupation matrix.

    Returns
    -------
    out: float
        The correction from LOSC to the total energy.

    See Also
    --------
    CurvatureV1.kappa(), CurvatureV2.kappa(), ...: Return the curvature matrix.
    local_occupation(): Return the local occupation matrix.
    )pddoc");

    // orbital_energy_post_scf
    m.def("orbital_energy_post_scf", &losc::orbital_energy_post_scf,
          R"pddoc(
    Calculate corrected orbital energy from LOSC in a post-SCF approach.

    This function gives the final orbital energies WITH the correction
    from LOSC. Note the difference to the function energy_correction() that only
    calculates the energy correction. The corrected orbital energies are the
    expectation values of converged DFA's COs on the LOSC-DFA Hamiltonian,
    that is,
    $\epsilon_i = \langle \psi_i | H_{\rm{dfa}} + H_{\rm{losc}} | \psi_i \rangle.$

    Parameters
    ----------
    H_dfa: np.ndarray [nbasis, nbasis]
        The DFA Hamiltonian under AOs.
    H_losc: np.ndarray [nbasis, nbasis]
        The LOSC effective Hamiltonian under AOs.
    C_co: np.ndarray [nbasis, n], n <= basis, which is the number of COs.
        The coefficient matrix of converged COs under AOs from DFA.

    Returns
    -------
    out: list
        The corrected orbital energies from LOSC with size of n. The order of
        orbital energies match the order of input COs (order of columns in
        C_co matrix).

    See Also
    --------
    ao_hamiltonian_correction(): obtain the LOSC effective Hamiltonian under AOs.

    Notes
    -----
    This function is just one of the ways to construct the LOSC corrected orbital
    energy in a post-SCF LOSC calculation. It is the way we used to produce
    results in the published paper for the post-SCF LOSC calculations. Besides
    this way, there are another two ways to calculate corrected orbital energies:
    (1) diagonalize the corrected LOSC-DFA Hamiltonian; (2) Follow Eq. 11 to
    calculate the corrections to orbital energies. These three ways usually
    prouce very similar results. Here, we only provide the commonly used
    approach.
    )pddoc");
}
