==========================
Use LOSC library in psi4
==========================

------------
Introduction
------------

``psi4_losc`` is a Python module that extends the `psi4 <https://psicode.org/>`_,
an open-source package for quantum chemistry calculations, to perform the
calculations of LOSC method. The ``psi4_losc`` module interacts with the ``psi4``
package via its Python interface. Using this ``psi4_losc`` module along with
``psi4`` (the `python model
<https://psicode.org/psi4manual/master/psiapi.html>`_ of psi4), you can
perform post-SCF-LOSC [#losc1]_ [#losc2]_ and SCF-LOSC [#scf-losc]_
calculations for aufbau/non-aufbau systems
with integer/fractional numbers of electrons.

---------------------

---------------
Installing psi4
---------------

.. Warning:: ``psi4_losc`` module is developed based on psi4 version ``1.3.2``.
   To have stable performance, it would be better to install this version of
   psi4. Other versions of ``psi4`` are not tested.

To perform LOSC calculation with ``psi4_losc`` in ``psi4``, you need to
install ``psi4`` package first. If never use psi4 before, you are suggested to
just install the conda binary package of psi4:

.. code-block:: bash

    conda install psi4=1.3.2 -c psi4

If you want to set up psi4 from the scratch, you can get more guidance for
the installation of psi4 from its `documentation
<https://psicode.org/psi4manual/master/external.html>`_.

---------------------

-----------
Basic Guide
-----------

LOSC can significantly eliminate the delocalization error exisiting in
conventional density functional approximations (DFAs), such as local density
approximations (LDA), generalized gradient approximations (GGAs) and hybrid
functionals like B3LYP. With minimal delocalization error from LOSC-DFA
calculations, you can expect good description of physical properties,
such as the ionization potentials (IPs), electron affinities
(EAs), the photoemision spectra.

For most cases, you would be interested in performing LOSC calculations for
aufbau systems (ground-state) with integer number of electrons. There are two
ways to apply LOSC to an associated DFA:

1. the post-SCF-LOSC approach [#losc1]_ [#losc2]_

    Directly applying LOSC correction to the SCF-DFA calculations.
    In this approach, the total energy and orbital energies from the associated
    DFA are corrected. The electron density is not corrected from LOSC.

2. the SCF-LOSC approach [#scf-losc]_

    Using the LOSC effective Hamiltonian to correct the associated DFA
    Hamiltonian and perform the SCF calculation for LOSC-DFA.
    In this approach, the total energy, orbital energies and electron
    densities are all corrected.

To perform the calculation of LOSC via ``psi4_losc``, always import following:

.. code-block:: python

    import psi4
    import psi4_losc

No matter doing SCF-LOSC or post-SCF-LOSC,
an SCF calculation for the associated DFA should be performed in advance to
generate the corresponding converged psi4 wavefunction object.

.. note:: Although psi4 supports calculations with point group symmetry, LOSC
   does not support symmetry because the usage of localized orbitals. So you
   need to turn off the symmetry in psi4 calculation (set ``symettr c1`` for
   the system).

Now, taking a simple molecule, a stretched H2 molecule, as an example.
You first calculate the SCF-DFA in psi4.

.. code-block:: python
    :emphasize-lines: 6

    # A stretched H2 molecule with 10A bond length with symmetry turned off.
    mol = psi4.geometry("""
        0 1
        H 0 0 0
        H 10 0 0
        symmetry c1  # turn off symmetry
        """)

    # some basic setting in psi4.
    psi4.set_options({'basis': '6-31g',
                      'guess': 'core',
                      'reference': 'uks'})

    # Here we do a b3lyp calculation and let psi4 return the wfn object.
    E_dfa, dfa_wfn = psi4.energy('b3lyp', return_wfn=True)
    print(E_dfa)

    # You would see E_dfa = -0.5765567776548441

post-SCF-LOSC for integer system with aufbau occupation
-------------------------------------------------------

Following the aforementioned H2 example, you can do a post-SCF-LOSC calculation
via calling ``psi4_losc.post_scf_losc()`` function.
There are two key arguments you need to specify: the first one
is the DFA functional type, and the second one the corresponding converged
DFA wavefunction object. ``psi4_losc.post_scf_losc()`` will return two
variables: the corrected total energy and orbital energies.

.. code-block:: python
   :emphasize-lines: 1-3

    # do post-SCF-LOSC-B3LYP calculation to obtain corrected total energy
    # and the orbital energies.
    E_losc, Orb_losc = psi4_losc.post_scf_losc(psi4_losc.B3LYP, dfa_wfn)

    # You would see total energies:
    # E_losc = -0.758073589993662
    # E_losc - E_dfa = 0.13428778258589236
    print(E_losc)
    print(E_losc - E_dfa)

    # You would see all LOSC corrected orbital energies in a.u. for
    # both alpha and beta spin:
    # [array([-6.98665404, -5.54668941, 24.04349262, 24.04349294]),
    #  array([-6.98665404, -5.54668941, 24.04349262, 24.04349294])]
    print(Orb_losc)


SCF-LOSC for integer system with aufbau occupation
--------------------------------------------------

Following the aforementioned H2 example, you can do an SCF-LOSC calculation via
calling ``psi4_losc.scf_losc()`` function. The input arguments are very similar
to ``psi4_losc.post_scf_losc()`` function. But the return type is a psi4
wavefunction object.

.. code-block:: python
   :emphasize-lines: 1-2

    # do SCF-LOSC-B3LYP calculation to obtain a LOSC-B3LYP wavefunction object.
    losc_wfn = psi4_losc.scf_losc(psi4_losc.B3LYP, dfa_wfn)

    # You would see total energies:
    # E_losc = xxx
    # E_losc - E_dfa = xxx
    E_losc = losc_wfn.energy()
    print(E_losc)
    print(E_losc - E_dfa)

    # You would see all LOSC corrected orbital energies in a.u. for alpha spin.
    #  xxx
    import numpy as np
    Orb_losc = np.asarray(losc_wfn.epsilon_a())
    print(Orb_losc)


---------------------

--------
Advanced
--------

Configure LOSC calculations with options
----------------------------------------

You can configure the options for the LOSC calculations in ``psi4_losc`` module.
These configurations include option settings related to LOSC curvatures,
localizations and so on. See details in :ref:`this section
<psi4_losc__options_label>`.


Setting energy window for LOSC calculation
------------------------------------------

The LOSC corrections are constructed from a set of localized orbitals (LOs).
These LOs are obtained by a unitary transformation from a set of canonical
orbitals (COs) from the associated DFA. In default (for
``psi4_losc.post_scf_losc`` and ``psi4_losc.scf_losc``), all the COs (the same
number of basis sets) are involved
to perform the localization. As results, LOSC corrects all the orbital energies.
To simplify the calculation, you can only select a subset of COs to perform the
LOSC calculation, and you should expect this simplified approach to produce
similar results from the calculation will all COs localized.

``psi4_losc`` module accepts the selection of COs with an energy window
(**in eV**) setting, that is selecting all the COs whose orbital energies are
inside the window. To enable the window, you use the `window` key argument.

.. code-block:: python

    # do post-SCF-LOSC-B3LYP calculation with setting a window of -30 - 10 eV.
    E_losc, Orb_losc = psi4_losc.post_scf_losc(psi4_losc.B3LYP, dfa_wfn,
                                               window=[-30, 10])

    # do SCF-LOSC-B3LYP calculation with setting a window of -30 - 10 eV.
    losc_wfn = psi4_losc.scf_losc(psi4_losc.B3LYP, dfa_wfn, window=[-30, 10])

.. note::
   Using window setting in the calculation of LOSC reduces the space of LOs and
   makes the calculation much faster. However, keep in mind that doing so
   would only produce corrections to orbital energies of the selected COs.
   **The orbital energies of non-selected COs will not be touched.**

   **To use the window setting, you should select most of the valence orbitals.**
   The usual energy window is -30 - 10 eV. Using this window
   should give you very similar results in total energy, orbital energies and
   electron density to the ones calculated from all COs localized.
   Excluding the core orbitals to the localization is a reasonble choice.
   This is because: (1) core orbitals usually are already localized, thus,
   they do not contribute the total energy correction because of the integer
   occupation number (see the LOSC energy correction formular [#losc1]_);
   (2) the core orbital energies are not as much interesting as the valence
   orbital energies.

---------------------


Systems with fractional numbers of electrons or non-aufbau occupation
---------------------------------------------------------------------
Psi4 mainly supports the calculations of systems with integer number of
electrons and aufbau occupations. However, it would be interesting to calculate
systems with fractional number of electrons (such as study of the delocalization
error) or non-aufbau occupation (such as :math:`\Delta`-SCF calculations).

``psi4_losc.scf`` module provides the extended SCF procedure to enable
calculations for systems with fractional numbers of electrons or non-aufbau
occupation. See :ref:`this <psi4_losc__scf_label>` for more details.

----------
References
----------

**psi4_losc Module**
--------------------

.. This includes the full list of psi4 API used in psi4_losc module.

.. include:: psi4_api_ref__psi4_losc.rst

**API of psi4_losc module**
***************************

.. autodata:: psi4_losc.B3LYP
.. autodata:: psi4_losc.BLYP
.. autodata:: psi4_losc.PBE
.. autodata:: psi4_losc.PBE0
.. autodata:: psi4_losc.GGA
.. autofunction:: psi4_losc.post_scf_losc
.. autofunction:: psi4_losc.scf_losc

------------


**psi4_losc.scf Module**
------------------------

.. This includes the full list of psi4 API used in psi4_losc.scf module.

.. include:: psi4_api_ref__psi4_losc__scf.rst

.. _psi4_losc__scf_label:

**API of psi4_losc.scf module**
*******************************

.. Warning:: This module is not fully tested. Use it with caution.

.. automodule:: psi4_losc.scf
    :members:
    :inherited-members:

--------------------

.. _psi4_losc__options_label:

**Configure LOSC options in psi4_losc**
---------------------------------------

The configuration of LOSC calculation for `psi4_losc` module is controlled
by `psi4_losc.options` variable. These configurations include the adjustments
to the LOSC curvature, and localizations. The details are shown below.

.. autodata:: psi4_losc.options
.. autoclass:: psi4_losc.losc_options.Options
    :members:

**Literature**
--------------
.. [#losc1] Li, C.; Zheng, X.; Su, N. Q.; Yang, W. Localized Orbital Scaling
   Correction for System- atic Elimination of Delocalization Error in
   Density Functional Approximations.
   `Natl. Sci. Rev. 2018, 5, 203−215. 203-215.
   <https://doi.org/10.1093/nsr/nwx111>`_

.. [#losc2] Su, N. Q.; Mahler, A.; Yang, W. Preserving Symmetry and
   Degeneracy in the Localized Orbital Scaling Correction Approach. J.
   `Phys. Chem. Lett. 2020, 11, 1528−1535.
   <https://doi.org/10.1021/acs.jpclett.9b03888>`_

.. [#scf-losc] Mei, Y.; Chen, Z.; Yang, W.
   Self-Consistent Calculation of the Localized Orbital Scaling Correction
   for Correct Electron Densities and Energy-Level Alignments in Density
   Functional Theory.
   `J. Phys. Chem. Lett. 2020, 11, 23, 10269–10277.
   <https://doi.org/10.1021/acs.jpclett.0c03133>`_
