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
The references for the detailed methodology of LOSC are listed
:ref:`here <introduction:references of losc>`.

-----------------------------
What is the LOSC library for?
-----------------------------
The LOSC library is for two things:

- It provides sub-libraries with compatibilities to several programming
  languages, including :ref:`C <c_losc:losc c library>`,
  :ref:`C++ <losc:losc c++ library>` and
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

.. [#losc1] `Li, Chen, et al. "Localized orbital scaling correction for
   systematic elimination of delocalization error in density functional
   approximations." National Science Review 5.2 (2018): 203-215.
   <https://doi.org/10.1093/nsr/nwx111>`_

.. [#losc2] `Su, Neil Qiang, Aaron Mahler, and Weitao Yang.
   "Preserving symmetry and degeneracy in the localized orbital scaling
   correction approach."
   The journal of physical chemistry letters 11.4
   (2020): 1528-1535.
   <https://doi.org/10.1021/acs.jpclett.9b03888>`_

.. [#scf-losc] `Mei, Yuncai, Zehua Chen, and Weitao Yang.
   "Self-Consistent Calculation of the Localized Orbital Scaling
   Correction for Correct Electron Densities and Energy-Level Alignments
   in Density Functional Theory."
   The Journal of Physical Chemistry Letters 11.23 (2020): 10269-10277.
   <https://doi.org/10.1021/acs.jpclett.0c03133>`_
