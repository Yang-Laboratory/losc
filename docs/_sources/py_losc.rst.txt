===================
LOSC Python library
===================

------------
Introduction
------------
``py_losc`` is the Python module that provides the Python interface to do
calculations of LOSC. It wraps the
:ref:`LOSC C++ library <losc:losc c++ library>` with the help of
`pybind11 library <https://pybind11.readthedocs.io/en/latest/>`_.
The functionalities provided in ``py_losc`` library are similar to the
:ref:`LOSC C++ library <losc:losc c++ library>`.

------------
Basic Guide
------------

To use Python interface for the LOSC library, all you need to do is
``import py_losc``.

- To construct LOSC curvature matrix,
  see :ref:`this section <py_losc:losc curvature>`.
- To perform LOSC localization,
  see :ref:`this section <py_losc:LOSC Localization>`.
- To calculate LOSC corrections,
  see :ref:`this section <py_losc:LOSC corrections>`.

To implement post-SCF-LOSC and SCF-LOSC calculations in Python with
the ``py_losc`` module, please refer to :ref:`this section
<psi4_losc:use losc library in psi4>`, which demonstrates a real
example of using ``py_losc`` in psi4 package.

-------------------------------
Detailed References for the API
-------------------------------

.. automodule:: py_losc

DFA Representation in LOSC Library
**********************************
.. autoclass:: py_losc.DFAInfo
    :members:
    :inherited-members:

LOSC Curvature
**************
.. autoclass:: py_losc.CurvatureV1
    :members:
    :inherited-members:
.. autoclass:: py_losc.CurvatureV2
    :members:
    :inherited-members:


LOSC Localization
*****************
.. autoclass:: py_losc.LocalizerV2
    :members:
    :inherited-members:

LOSC Local Occupation
*********************
.. autofunction:: py_losc.local_occupation

LOSC Corrections
****************
.. autofunction:: py_losc.ao_hamiltonian_correction
.. autofunction:: py_losc.energy_correction
.. autofunction:: py_losc.orbital_energy_post_scf
