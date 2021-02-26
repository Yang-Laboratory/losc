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

    Parameters
    ----------
    S : numpy.ndarray
        AO overlap matrix with dimension [nbasis, nbasis].
    C_lo : numpy.array
        LO coefficient matrix under AO basis with dimension of [nbasis, nbasis].
        The coefficients of i-th LO is the i-th column of `C_lo`.
    Curvature : numpy.array
        LOSC curvature matrix with dimension [nlo, nlo].
    LocalOcc : numpy.array
        LOSC local occupation matrix with dimension [nlo, nlo].

    Returns
    -------
    numpy.array
        The LOSC effective Hamiltonian under AOs with dimension [nbasis, nbasis].

    See Also
    --------
    CurvatureV1.kappa, CurvatureV2.kappa, local_occupation

    Notes
    -----
    The LOSC effective Hamiltonian is constructed with LOs fixed. The
    expression of the effective Hamiltonian is shown as Eq. S25 in the
    supporting information of the `original LOSC paper
    <https://doi.org/10.1093/nsr/nwx111>`_. This effective Hamiltonian is exact
    in the `developed version of SCF-LOSC
    <https://doi.org/10.1021/acs.jpclett.0c03133>`_
    )pddoc");

    // energy_correction
    m.def("energy_correction", &losc::energy_correction, R"pddoc(
    Calculate the total energy correction from LOSC.

    Parameters
    ----------
    Curvature : numpy.ndarray
        The LOSC curvature matrix with dimension [nlo, nlo].
    LocalOcc: numpy.ndarray
        The LOSC local occupation matrix with dimension [nlo, nlo].

    Returns
    -------
    float
        The correction from LOSC to the total energy.

    See Also
    --------
    CurvatureV1.kappa, CurvatureV2.kappa, local_occupation

    Notes
    -----
    This is just the energy correction from LOSC, NOT the total energy of
    LOSC-DFA. Total energy of LOSC-DFA is ``E_losc_dfa = E_dfa + E_losc``.
    )pddoc");

    // orbital_energy_post_scf
    m.def("orbital_energy_post_scf", &losc::orbital_energy_post_scf,
          R"pddoc(
    Calculate corrected orbital energy from LOSC in a post-SCF approach.



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
    ao_hamiltonian_correction

    Notes
    -----
    This function gives the final corrected orbital energies from LOSC.
    Note the difference to the function `energy_correction`.

    The corrected orbital energies are the expectation values of converged
    DFA's COs on the LOSC-DFA Hamiltonian, that is,

    .. math:: \epsilon_i = \langle \psi_i | H_{\rm{dfa}} + H_{\rm{losc}}
       | \psi_i \rangle.


    This is just one of the ways to calculate the LOSC corrected orbital
    energy in a post-SCF LOSC calculation. It is the way usually used to produce
    results in the `published paper <https://doi.org/10.1093/nsr/nwx111>`_
    for the post-SCF LOSC calculations. Besides this way, there are another two
    ways to calculate corrected orbital energies: (1) diagonalize the corrected
    LOSC-DFA Hamiltonian; (2) Follow
    `Eq. 11 <https://doi.org/10.1093/nsr/nwx111>`_ to calculate the corrections
    to orbital energies. These three ways usually produce similar results.
    )pddoc");
}
