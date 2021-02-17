"""
Import this module will update function `psi4.proc.scf_wavefunction_factory()`
defined in psi4 to be an extended version. Making such update is to enable
`psi4.energy()` to be compatible with LOSC calculation.
"""

import psi4
import psi4_losc
import psi4_losc.wfn
from psi4 import core
from psi4.driver.p4util.exceptions import ValidationError

# Function `psi4.proc.scf_wavefunction_factory()` will be extended in this
# module. Here, we save the original one first.
_psi4_scf_wavefunction_factory = psi4.proc.scf_wavefunction_factory


def _scf_wavefunction_factory_losc(name, ref_wfn, reference, **kwargs):
    """
    Build an SCF wavefunction for SCF-LOSC. This function is modified on the
    basis of `psi4.proc.scf_wavefunction_factory()`.
    """
    losc_data = kwargs.get('losc_data', {})
    # Figure out functional and dispersion
    superfunc, _disp_functor = psi4.proc.build_disp_functor(
        name, restricted=(reference in ["RKS", "RHF"]), **kwargs)

    # Build the wavefunction
    psi4.core.prepare_options_for_module("SCF")
    # Do SCF-LOSC: build LOSC wavefunction.
    if not losc_data:
        raise Exception('need LOSC data to build LOSC class.')
    curvature = losc_data.get('curvature', [])
    C_lo = losc_data.get('C_lo', [])
    if not curvature:
        raise Exception('need LOSC curvature matrix to build LOSC class.')
    if not C_lo:
        raise Exception('need LOSC localized orbital to build LOSC class.')

    # restricted
    if reference in ["RHF", "RKS"]:
        if len(curvature) != 1 or len(C_lo) != 1:
            raise Exception(
                'restricted case: size of curvature or C_lo is not 1.')
        wfn = psi4_losc.wfn.RLOSC(ref_wfn, superfunc, curvature, C_lo)
    # unrestricted
    elif reference in ['UHF', 'UKS']:
        if len(curvature) != 2 or len(C_lo) != 2:
            raise Exception(
                'unrestricted case: size of curvature or C_lo is not 2.')
        wfn = psi4_losc.wfn.ULOSC(ref_wfn, superfunc, curvature, C_lo)
    # error reference
    else:
        raise ValidationError(
            f"reference ({reference}) not supported for LOSC wavefunction.")

    if _disp_functor and _disp_functor.engine != 'nl':
        wfn._disp_functor = _disp_functor

    # Set the DF basis sets
    if (("DF" in core.get_global_option("SCF_TYPE")) or
            (core.get_option("SCF", "DF_SCF_GUESS") and (core.get_global_option("SCF_TYPE") == "DIRECT"))):
        aux_basis = core.BasisSet.build(wfn.molecule(), "DF_BASIS_SCF",
                                        core.get_option("SCF", "DF_BASIS_SCF"),
                                        "JKFIT", core.get_global_option(
                                            'BASIS'),
                                        puream=wfn.basisset().has_puream())
        wfn.set_basisset("DF_BASIS_SCF", aux_basis)
    else:
        wfn.set_basisset("DF_BASIS_SCF", core.BasisSet.zero_ao_basis_set())

    # Set the relativistic basis sets
    if core.get_global_option("RELATIVISTIC") in ["X2C", "DKH"]:
        decon_basis = core.BasisSet.build(wfn.molecule(), "BASIS_RELATIVISTIC",
                                          core.get_option(
                                              "SCF", "BASIS_RELATIVISTIC"),
                                          "DECON", core.get_global_option(
                                              'BASIS'),
                                          puream=wfn.basisset().has_puream())
        wfn.set_basisset("BASIS_RELATIVISTIC", decon_basis)

    # Set the multitude of SAD basis sets
    if (core.get_option("SCF", "GUESS") in ["SAD", "SADNO", "HUCKEL"]):
        sad_basis_list = core.BasisSet.build(wfn.molecule(), "ORBITAL",
                                             core.get_global_option("BASIS"),
                                             puream=wfn.basisset().has_puream(),
                                             return_atomlist=True)
        wfn.set_sad_basissets(sad_basis_list)

        if ("DF" in core.get_option("SCF", "SAD_SCF_TYPE")):
            # We need to force this to spherical regardless of any user or other demands.
            optstash = psi4.p4util.OptionsState(['PUREAM'])
            core.set_global_option('PUREAM', True)
            sad_fitting_list = core.BasisSet.build(wfn.molecule(), "DF_BASIS_SAD",
                                                   core.get_option(
                                                       "SCF", "DF_BASIS_SAD"),
                                                   puream=True,
                                                   return_atomlist=True)
            wfn.set_sad_fitting_basissets(sad_fitting_list)
            optstash.restore()

    # Deal with the EXTERN issues
    if hasattr(core, "EXTERN"):
        wfn.set_external_potential(core.EXTERN)

    return wfn


def _scf_wavefunction_factory_extended_version(name, ref_wfn, reference, **kwargs):
    """
    Build an SCF wavefunction. This is the extended version for
    `psi4.proc.scf_wavefunction_factory()`.
    """
    reference = reference.upper()
    losc_data = kwargs.get('losc_data', {})
    # psi4 default wavefunction: integer and aufbau system
    if not losc_data:
        wfn = _psi4_scf_wavefunction_factory(
            name, ref_wfn, reference, **kwargs)
    # LOSC wavefunction or wavefunction with customized occupation numbers.
    else:
        wfn = _scf_wavefunction_factory_losc(
            name, ref_wfn, reference, **kwargs)
    return wfn


# Here, we update `psi4.proc.scf_wavefunction_factory()` with the extended
# version.
psi4.proc.scf_wavefunction_factory = _scf_wavefunction_factory_extended_version
