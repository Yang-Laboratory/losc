"""
Python interface for localized orbital scaling correction (LOSC) library.

Notations
---------
LO refers to the localized orbital.
AO refers to the atomic orbitals.
CO refers to the canonical orbital.
`nlo` refers to the number of LOs.
`nbasis` refers to the number of AOs.
`nfitbasis` refers to the number of fitting basis for density fitting.
`npts` refers to the number of grids.

Notes
-----
1. All the matrices are represented and stored in 2-dimensional
`np.ndarray` object. The storage order of matrices is required to be
C-style (row major). If the input matrix is not C-style, it will be converted
interanally.
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
