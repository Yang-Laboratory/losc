psi4_losc
=========

Introduction
------------
``psi4_losc`` is a Python module that enables doing LOSC calculations in
`psi4 <https://psicode.org/>`_, an open-source package to perform quantum
chemistry calculations. The ``psi4_losc`` module interacts with the ``psi4``
package via its Python interface to achieve LOSC calculation. Using this
``psi4_losc`` module along with ``psi4`` (the `python model
<https://psicode.org/psi4manual/master/psiapi.html>`_ of psi4), you can
perform post-SCF-LOSC and SCF-LOSC calculations for aufbau/non-aufbau systems
with integer/fractional number of electrons.

Installing psi4
---------------
``psi4_losc`` module is developed based on psi4 version ``1.3.2``.
To have stable performance, it would be better to install this version of psi4.
If you never use psi4 before, you are suggested to just install the conda binary
package of psi4:

.. code-block:: bash

    >>> conda install psi4=1.3.2 -c psi4

If you want to set up psi4 from the scratch, you can get more guidance for
the installation of psi4 from its `documentation
<https://psicode.org/psi4manual/master/external.html>`_.

First Step
----------
LOSC can significantly eliminate the delocalization error exisiting in
conventional density functional approximations (DFAs), such as local density
approximations (LDA), generalized gradient approximations (GGAs) and hybrid
functionals like B3LYP. With minimal delocalization error from LOSC-DFA
calculations, you can expect good description of physical properties from
LOSC-DFA, such as the ionization potentials (IPs), electron affinities
(EAs), the photoemision spectra.

For most cases, you would be interested in performing LOSC calculations for
aufbau systems (ground-state) with integer number of electrons. There are two
ways to apply LOSC to an associated DFA for these systems.
One is the post-SCF-LOSC approach [#losc1]_ [#losc2]_, which means directly
applying LOSC correction to the SCF-DFA calculations. In this approach,
post-SCF-LOSC only corrects the total energy and orbital energies for the
associated DFA. The electron density is not touched for LOSC to the DFA.
The other approach is the SCF-LOSC approach [#scf-losc]_, which means using the
LOSC effective Hamiltonian to correct the associated DFA Hamiltonian and
perform the SCF calculation for LOSC-DFA. The total energy, orbital energies
and electron densities are all corrected from this approach.

No matter doing SCF-LOSC or post-SCF-LOSC by using this ``psi4_losc`` module,
an SCF-DFA calculation should always be performed in advance to generate
the corresponding converged psi4 wavefunction object. Although psi4
supports calculations with point group symmetry, LOSC does not support
symmetry because the usage of localized orbitals. So you need to turn off
the symmetry in psi4 calculation.

In the following text, following two lines are assumed to be executed to
import ``psi4`` and ``psi4_losc`` modules.

.. code-block:: python

    import psi4
    import psi4_losc

Now, taking a simple molecule, a stretched H2+ molecule with only one
electron, as an example. You first calculate the SCF-DFA in psi4.

.. code-block:: python

    # A stretched H2 molecule with 10A bond length with symmetry turned off.
    mol = psi4.geometry("""
        0 1
        H 0 0 0
        H 10 0 0
        symmetry c1
        """)

    # some basic setting in psi4.
    psi4.set_options({'basis': '6-31g',
                      'guess': 'core',
                      'reference': 'uks'})

    # Here we do a b3lyp calculation and let psi4 return the wfn object.
    E_dfa, dfa_wfn = psi4.energy('b3lyp', return_wfn=True)
    print(E_dfa)

    # You would see E_dfa = -0.5765567776548441

post-SCF-LOSC for ground state integer system
*********************************************
Following the aforementioned H2 example, let's do a post-SCF-LOSC calculation.
You need to call the ``psi4_losc.post_scf_losc()`` function to perform the
calculation. There are two key arguments you need to specify: the first one
is the DFA functional type, and the second one the corresponding converged
DFA wavefunction object. ``psi4_losc.post_scf_losc()`` will return two
variables: the corrected total energy and orbital energies.

.. code-block:: python

    # do post-SCF-LOSC-B3LYP calculation
    E_losc, Orb_losc = psi4_losc.post_scf_losc(psi4_losc.B3LYP, dfa_wfn)

    print(E_losc)
    print(E_losc - E_dfa)

    # You would see total eneriges:
    # E_losc = -0.758073589993662
    # E_losc - E_dfa = 0.13428778258589236

    print(Orb_losc)

    # You would see all LOSC corrected orbital energies in a.u. for
    # both alpha and beta spin:
    # [array([-6.98665404, -5.54668941, 24.04349262, 24.04349294]),
    #  array([-6.98665404, -5.54668941, 24.04349262, 24.04349294])]



SCF-LOSC for ground state integer system
****************************************
Following the aforementioned H2 example, let's do an SCF-LOSC calculation.
You need to call the ``psi4_losc.scf_losc()`` function to perform the
calculation. The input arguments are very similar to
``psi4_losc.post_scf_losc()`` function. But the return is a psi4 wavefunction
object.

.. code-block:: python

    # do SCF-LOSC-B3LYP calculation
    losc_wfn = psi4_losc.scf_losc(psi4_losc.B3LYP, dfa_wfn)

    E_losc = losc_wfn.energy()
    print(E_losc)
    print(E_losc - E_dfa)

    # You would see total eneriges:
    # E_losc = xxx
    # E_losc - E_dfa = xxx

    import numpy as np
    Orb_losc = np.asarray(losc_wfn.epsilon_a())
    print(Orb_losc)

    # You would see all LOSC corrected orbital energies in a.u. for alpha spin.
    #  xxx


Advanced
--------
TODO: fractional calculations.


Python interface for psi4
-------------------------
Below is the complete list of Python interfaces (functions) of ``psi4``
that are used in ``psi4_losc`` module:

- `psi4.driver <>` module

    - `driver.p4util.OptionsState <https://psicode.org/psi4manual/master/optionshandling.html#handling-options-in-driver>`_

        - OptionsState.restore()
    - `driver.scf_wavefunction_factory()`

- `psi4.core.clean() <>`_
- `psi4.core.get_active_molecule() <>`_
- `psi4.core.get_option() <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.get_option>`_
- `psi4.core.set_local_option() <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.set_local_option>`_
- `psi4.core.BasisSet <https://psicode.org/psi4manual/master/api/psi4.core.basisset>`_

    - `BasisSet.nbf() <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.BasisSet.nbf>`_
    - `BasisSet.build() <>`_

- `psi4.core.MintsHelper <https://psicode.org/psi4manual/master/api/psi4.core.mintshelper>`_

    - `psi4.core.MintsHelper.ao_dipole() <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.MintsHelper.ao_dipole>`_
    - `psi4.core.MintsHelper.ao_eri() <>'_

- `psi4.core.Wavefunction <https://psicode.org/psi4manual/master/api/psi4.core.wavefunction>`_

    - `Wavefunction.molecule() <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.molecule>`_
    - Wavefunction.Ca(), Wavefunction.Cb()
    - Wavefunction.Fa(), Wavefunction.Fb()
    - Wavefunction.Da(), Wavefunction.Db()
    - Wavefunction.nalphe(), Wavefunction.nbeta()
    - Wavefunction.epsilon_a(), Wavefunction.epsilon_b()
    - Wavefunction.energy()
    - Wavefunction.basisset()
    - Wavefunction.same_a_b_orbs()
    - Wavefunction.same_a_b_dens()
    - Wavefunction.build()
    - Wavefunction.form_Shalf()
    - Wavefunction.guess()

- `psi4.core.Matrix`

- `psi4.core.JK`

    - JK.build()
    - JK.initialize()
    - JK.compute()
    - JK.C_clear()
    - JK.C_left_add()
    - JK.J()
    - JK.K()

- `psi4.core.HF <https://psicode.org/psi4manual/master/api/psi4.core.hf>`_

    - HF.functional()
    - HF.V_potential()

- `psi4.core.SuperFunctional`

    - `SuperFunctional.needs_xc()`
    - `SuperFunctional.x_alpha()`

- psi4.core.Molecule

    - Molecule.schoenflies_symbol()
    - Molecule.print_out()
    - Molecule.nuclear_repulsion_energy()

- psi4.core.VBase

    - VBase.set_D()
    - VBase.compute_V()
    - VBase.quadrature_values()



References
----------
**psi4_losc module**
********************
.. automodule:: psi4_losc.psi4_losc
    :members:
    :inherited-members:

**psi4_losc.scf module**
********************
.. automodule:: psi4_losc.scf
    :members:
    :inherited-members:


.. [#losc1] `Li, Chen, et al. "Localized orbital scaling correction for
    systematic elimination of delocalization error in density functional
    approximations." National Science Review 5.2 (2018): 203-215.
    <https://doi.org/10.1093/nsr/nwx111>`_
.. [#losc2] `Su, Neil Qiang, Aaron Mahler, and Weitao Yang.
    "Preserving symmetry and degeneracy in the localized orbital scaling
    correction approach." The journal of physical chemistry letters 11.4
    (2020): 1528-1535.
    <https://doi.org/10.1021/acs.jpclett.9b03888>`_
.. [#scf-losc] `Mei, Yuncai, Zehua Chen, and Weitao Yang.
    "Self-Consistent Calculation of the Localized Orbital Scaling
    Correction for Correct Electron Densities and Energy-Level Alignments
    in Density Functional Theory."
    The Journal of Physical Chemistry Letters 11.23 (2020): 10269-10277.
    <https://doi.org/10.1021/acs.jpclett.0c03133>`_
