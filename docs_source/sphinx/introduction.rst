============
Introduction
============

------------------------
What is the LOSC method?
------------------------

The localized orbital scaling correction (LOSC) is a newly developed method in
density functional theory (DFT) to eliminate the delocalization error (DE)
of many conventional density functional approximations (DFAs). Bascally,
the LOSC method involves a set of localized orbitals (LOs) to construct
the local occupation number, the LOSC curvature matrix and the LOSC corrections.
The references of LOSC that describes the detailed methodology are listed
:ref:`here <introduction:references of losc>`.

-----------------------------
What is the LOSC library for?
-----------------------------
The LOSC library is developed for two goals:

- It provides sub-libraries with compatible interfaces to several popular
  programming languages in quantum chemistry, including
  :ref:`C <c_losc:losc c library>`, :ref:`C++ <losc:losc c++ library>` and
  :ref:`Python <py_losc:losc python library>`, for the developers who would
  be interested to implement LOSC method in their favorite quantum chemistry
  packages. These sub-libraries provide the functionalities to perform the
  essential calculations of the LOSC method, such as constructing the LOs,
  LOSC curvature matrix, and LOSC corrections.

- It provides the implementation of LOSC method in an open-source quantum
  chemistry package, `psi4 <https://psicode.org/>`_, with the Python interface.
  If you are a user of pis4, you can directly use this library with psi4
  to perform LOSC calculations.
  See :ref:`this section <psi4_losc:use LOSC library in psi4>`.

------------------
References of LOSC
------------------

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
