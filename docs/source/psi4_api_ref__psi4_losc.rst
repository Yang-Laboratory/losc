**Full list of psi4 Python API used in psi4_losc module**
*********************************************************

.. Warning:: The list may not be complete. The list is for the package
   developers to show the interaction between ``psi4_losc`` module and
   ``psi4`` Python API.

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
