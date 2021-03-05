**Full list of psi4 Python API used in psi4_losc.scf module**
*************************************************************

.. Warning:: The list may not be complete. The list is for the package
   developers to show the interaction between ``psi4_losc.scf`` module and
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
