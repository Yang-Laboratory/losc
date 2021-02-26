"""Python interface for localized orbital scaling correction (LOSC) library.

Notes
-----
All the matrices are represented and stored in the 2-dimensional
`numpy.array` object. The storage order of matrices is required to be
C-style (row major). If the input matrix is not C-style, it will be converted
interanally.

Following notations are used in the docstrings:

    - LO: the localized orbital.
    - AO: the atomic orbitals.
    - CO: the canonical orbital.
    - nlo: the number of LOs.
    - nbasis: the number of AOs.
    - nfitbasis: the number of fitting basis for density fitting.
    - npts: the number of grids.
"""

# ==> Interface for LOSC corrections <==
from py_losc.py_losc_core import ao_hamiltonian_correction
from py_losc.py_losc_core import orbital_energy_post_scf
from py_losc.py_losc_core import energy_correction

# ==> Interface for LOSC local occupation <==
from py_losc.py_losc_core import local_occupation


# ==> Interface for LOSC localization <==
from py_losc.py_losc import DFAInfo
from py_losc.py_losc import CurvatureV1
from py_losc.py_losc import CurvatureV2

# ==> Interface for LOSC curvature matrix <==
from py_losc.py_losc import LocalizerV2
