=========
psi4_losc
=========

------------
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

---------------------

---------------
Installing psi4
---------------

``psi4_losc`` module is developed based on psi4 version ``1.3.2``.
To have stable performance, it would be better to install this version of psi4.
If you never use psi4 before, you are suggested to just install the conda binary
package of psi4:

.. code-block:: bash

    conda install psi4=1.3.2 -c psi4

If you want to set up psi4 from the scratch, you can get more guidance for
the installation of psi4 from its `documentation
<https://psicode.org/psi4manual/master/external.html>`_.

---------------------

----------
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
---------------------------------------------

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

    # You would see total energies:
    # E_losc = -0.758073589993662
    # E_losc - E_dfa = 0.13428778258589236

    print(Orb_losc)

    # You would see all LOSC corrected orbital energies in a.u. for
    # both alpha and beta spin:
    # [array([-6.98665404, -5.54668941, 24.04349262, 24.04349294]),
    #  array([-6.98665404, -5.54668941, 24.04349262, 24.04349294])]


SCF-LOSC for ground state integer system
----------------------------------------

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

    # You would see total energies:
    # E_losc = xxx
    # E_losc - E_dfa = xxx

    import numpy as np
    Orb_losc = np.asarray(losc_wfn.epsilon_a())
    print(Orb_losc)

    # You would see all LOSC corrected orbital energies in a.u. for alpha spin.
    #  xxx

---------------------

--------
Advanced
--------

TODO: fractional calculations.

---------------------

----------
References
----------

**psi4_losc module**
--------------------

**psi4_losc: used psi4 Python API**
***********************************

Below is the complete list of psi4 Python interfaces that are used in
``psi4_losc`` module.

- psi4.driver

    - `driver.p4util.OptionsState
      <https://psicode.org/psi4manual/master/optionshandling.html#handling-options-in-driver>`_

        - OptionsState.restore()

    - driver.scf_wavefunction_factory()

- `psi4.core.clean()
  <https://psicode.org/psi4manual/master/api/psi4.core.clean.html>`_
- `psi4.core.get_active_molecule()
  <https://psicode.org/psi4manual/master/api/psi4.core.get_active_molecule.html>`_
- `psi4.core.get_option()
  <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.get_option>`_
- `psi4.core.set_local_option()
  <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.set_local_option>`_

- `psi4.core.BasisSet
  <https://psicode.org/psi4manual/master/api/psi4.core.basisset>`_

    - `nbf()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.BasisSet.nbf>`_
    - `build()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.BasisSet.build>`_

- `psi4.core.MintsHelper
  <https://psicode.org/psi4manual/master/api/psi4.core.mintshelper>`_

    - `ao_dipole()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.MintsHelper.ao_dipole>`_

- `psi4.core.Wavefunction
  <https://psicode.org/psi4manual/master/api/psi4.core.wavefunction>`_

    - `molecule()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.molecule>`_
    - `S()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.S>`_
    - `Ca()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.Ca>`_
    - `Cb()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.Cb>`_
    - `Fa()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.Fa>`_
    - `Fb()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.Fb>`_
    - `Da()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.Da>`_
    - `Db()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.Db>`_
    - `nalphe()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.nalpha>`_
    - `nbeta()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.nbeta>`_
    - `epsilon_a()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.epsilon_a>`_
    - `epsilon_b()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.epsilon_b>`_
    - `energy()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.energy>`_
    - `basisset()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.basisset>`_
    - `same_a_b_orbs()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.same_a_b_orbs>`_
    - `same_a_b_dens()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.same_a_b_dens>`_
    - `build()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.build>`_

- `psi4.core.HF
  <https://psicode.org/psi4manual/master/api/psi4.core.hf>`_

    - `functional()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.HF.functional>`_

- `psi4.core.SuperFunctional
  <https://psicode.org/psi4manual/master/api/psi4.core.superfunctional>`_

    - `is_x_lrc()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.SuperFunctional.is_x_lrc>`_
    - `is_c_hybrid()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.SuperFunctional.is_c_hybrid>`_
    - `is_meta()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.SuperFunctional.is_meta>`_

- `psi4.core.Molecule
  <https://psicode.org/psi4manual/master/api/psi4.core.molecule>`_

    - `schoenflies_symbol()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Molecule.schoenflies_symbol>`_
    - `print_out()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Molecule.print_out>`_
    - `nuclear_repulsion_energy()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Molecule.nuclear_repulsion_energy>`_


**psi4_losc: API**
******************

.. autodata:: psi4_losc.B3LYP
.. autodata:: psi4_losc.BLYP
.. autodata:: psi4_losc.PBE
.. autodata:: psi4_losc.PBE0
.. autodata:: psi4_losc.GGA
.. autofunction:: psi4_losc.post_scf_losc
.. autofunction:: psi4_losc.scf_losc


**psi4_losc.scf module**
------------------------

**psi4_losc.scf: used psi4 Python API**
***************************************

Below is the complete list of psi4 Python interfaces that are used in
``psi4_losc.scf`` module.

- psi4.driver

    - `driver.p4util.OptionsState
      <https://psicode.org/psi4manual/master/optionshandling.html#handling-options-in-driver>`_

        - OptionsState.restore()

    - driver.scf_wavefunction_factory()

- `psi4.core.clean()
  <https://psicode.org/psi4manual/master/api/psi4.core.clean.html>`_
- `psi4.core.get_active_molecule()
  <https://psicode.org/psi4manual/master/api/psi4.core.get_active_molecule.html>`_
- `psi4.core.get_option()
  <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.get_option>`_
- `psi4.core.set_local_option()
  <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.set_local_option>`_

- `psi4.core.BasisSet
  <https://psicode.org/psi4manual/master/api/psi4.core.basisset>`_

    - `nbf()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.BasisSet.nbf>`_
    - `build()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.BasisSet.build>`_

- `psi4.core.MintsHelper
  <https://psicode.org/psi4manual/master/api/psi4.core.mintshelper>`_

    - `ao_dipole()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.MintsHelper.ao_dipole>`_
    - `ao_eri()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.MintsHelper.ao_eri>`_

- `psi4.core.Wavefunction
  <https://psicode.org/psi4manual/master/api/psi4.core.wavefunction>`_

    - `molecule()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.molecule>`_
    - `Ca()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.Ca>`_
    - `Cb()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.Cb>`_
    - `Fa()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.Fa>`_
    - `Fb()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.Fb>`_
    - `Da()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.Da>`_
    - `Db()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.Db>`_
    - `nalphe()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.nalpha>`_
    - `nbeta()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.nbeta>`_
    - `epsilon_a()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.epsilon_a>`_
    - `epsilon_b()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.epsilon_b>`_
    - `energy()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.energy>`_
    - `basisset()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.basisset>`_
    - `same_a_b_orbs()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.same_a_b_orbs>`_
    - `same_a_b_dens()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.same_a_b_dens>`_
    - `build()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Wavefunction.build>`_

- `psi4.core.Matrix
  <https://psicode.org/psi4manual/master/api/psi4.core.matrix>`_

- `psi4.core.JK
  <https://psicode.org/psi4manual/master/api/psi4.core.jk>`_

    - `build()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.JK.build>`_
    - `initialize()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.JK.initialize>`_
    - `compute()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.JK.compute>`_
    - `C_clear()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.JK.C_clear>`_
    - `C_left_add()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.JK.C_left_add>`_
    - `J()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.JK.J>`_
    - `K()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.JK.K>`_

- `psi4.core.HF
  <https://psicode.org/psi4manual/master/api/psi4.core.hf>`_

    - `form_Shalf()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.HF.form_Shalf>`_
    - `guess()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.HF.guess>`_
    - `functional()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.HF.functional>`_
    - `V_potential()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.HF.V_potential>`_

- `psi4.core.SuperFunctional
  <https://psicode.org/psi4manual/master/api/psi4.core.superfunctional>`_

    - `needs_xc()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.SuperFunctional.needs_xc>`_
    - `x_alpha()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.SuperFunctional.x_alpha>`_

- `psi4.core.Molecule
  <https://psicode.org/psi4manual/master/api/psi4.core.molecule>`_

    - `schoenflies_symbol()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Molecule.schoenflies_symbol>`_
    - `print_out()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Molecule.print_out>`_
    - `nuclear_repulsion_energy()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.Molecule.nuclear_repulsion_energy>`_

- `psi4.core.VBase
  <https://psicode.org/psi4manual/master/api/psi4.core.vbase>`_

    - `set_D()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.VBase.set_D>`_
    - `compute_V()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.VBase.compute_V>`_
    - `quadrature_values()
      <https://psicode.org/psi4manual/master/psi4api.html#psi4.core.VBase.quadrature_values>`_


**psi4_losc.scf: API**
**********************

.. automodule:: psi4_losc.scf
    :members:
    :inherited-members:

**Literature**
--------------

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
