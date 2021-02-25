py_losc
=======

Introduction
------------
``py_losc`` is the Python module that provides the Python interface for the LOSC
library. It wraps the ``losc`` C++ library with the help of `pybind11
<https://pybind11.readthedocs.io/en/latest/>`_.


References
----------

**LOSC Curvature**
******************
.. autoclass:: py_losc.DFAInfo
    :members:
    :inherited-members:
.. autoclass:: py_losc.CurvatureV1
    :members:
    :inherited-members:
.. autoclass:: py_losc.CurvatureV2
    :members:
    :inherited-members:


**LOSC Localization**
*********************
.. autofunction:: py_losc.local_occupation
.. autoclass:: py_losc.LocalizerV2
    :members:
    :inherited-members:

**LOSC Corrections**
********************
.. autofunction:: py_losc.ao_hamiltonian_correction
.. autofunction:: py_losc.energy_correction
.. autofunction:: py_losc.orbital_energy_post_scf
