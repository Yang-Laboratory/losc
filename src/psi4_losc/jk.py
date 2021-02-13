from numpy.lib.function_base import unwrap
import psi4
import numpy as np

from psi4.driver.driver import optimize


class JK_psi4_jk:
    """
    JK matrix constructor that is wrapped from psi4.core.JK class.

    Notes
    -----
    1. This class can only construct the JK matrices for aufbau and integer
    systems. This is limited by the psi4.core.JK class.

    2. The psi4 global `SCF_TYPE` settings sets the algorithm to build the JK
    matrices.
    """

    def __init__(self, wfn, Cocc):
        """
        Parameters
        ----------
        wfn: psi4.core.wavefunction
            The wavefunction object.
        Cocc: [psi4.core.Matrix, ...], list of psi4.core.Matrix
            The occupied CO coefficient matrices.
        """
        jk = psi4.core.JK.build(wfn.basisset())
        jk.initialize()
        self._Cocc = Cocc
        self._jk = jk
        self._nspin = len(Cocc)

    def compute(self):
        """
        Compute the JK matrices based on current coefficient matrix.
        """
        for Ci in self._Cocc:
            self._jk.C_left_add(Ci)
        self._jk.compute()
        self._jk.C_clear()

    def J(self):
        """
        Returns
        -------
        out: [np.array, ...], list of np.array
            The J matrices. The size of list is the same as the input Cocc.
        """
        return [np.asarray(self._jk.J()[i]) for i in range(self._nspin)]

    def K(self):
        """
        Returns
        -------
        out: [np.array, ...], list of np.array
            The K matrices. The size of list is the same as the input Cocc.
        """
        return [np.asarray(self._jk.K()[i]) for i in range(self._nspin)]


class JK_psi4_mints:
    """
    JK matrix constructor that is based on psi4.core.libmints library.

    Notes
    -----
    1. This class can construct the JK matrices for general systems, including
    non-aufbau or fractional systems.

    2. The psi4 global `SCF_TYPE` settings sets the algorithm to build the JK
    matrices. Currently, it only supports `SCF_TYPE DIRECT` algorithm.
    """

    def __init__(self, wfn, occ_val, Cocc):
        """
        Parameters
        ----------
        wfn: psi4.core.wavefunction
            The wavefunction object.
        occ_val: list or np.array
            The occupation number corrresponding to `Cocc`.
        Cocc: [np.array, ...], list of np.array with dimension [nbasis, nocc].
            The occupied CO coefficient matrices.
        """
        # same sanity-check
        # check spin.
        if not occ_val or not Cocc:
            raise Exception(
                "Empty input is not allowed to construct JK_psi4_mints")
        if len(occ_val) != len(Cocc):
            raise Exception("Mismatched length for occ_val and Cocc.")
        self._nspin = len(occ_val)
        self._nocc = [len(x) for x in occ_val]
        self._nbf = Cocc[0].shape[0]
        # check dimension of Cocc.
        for s in range(self._nspin):
            if Cocc[s].shape != (self._nbf, self._nocc[s]):
                raise Exception(f"Detect wrong dimension for Cocc[{s}]")
        # check algorithm
        self._support_algo = ['DIRECT', 'DF']
        self._algo = psi4.core.get_global_option('SCF_TYPE')
        if self._algo not in self._support_algo:
            raise Exception(
                f"Algorithm used for building JK is not supported: {self._algo}")

        # create psi4 MintsHelper object.
        self._wfn = wfn
        self._mints = psi4.core.MintsHelper(wfn.basisset())
        self._occ_val = occ_val
        self._Cocc = Cocc

        # build AO intergral
        self._uvst = None  # AO ERI
        self._df_pmn = None  # DF AO ERI
        self._df_Vpq_inv = None  # DF inverse
        self._D = [None] * self._nspin
        if self._algo == 'DIRECT':
            self._compute_ao_eri_direct()
        elif self._algo == 'DF':
            self._compute_ao_eri_df()
        else:
            raise Exception(f"SCF type not support: {self._algo}")


    def _compute_ao_eri_direct(self):
        self._uvst = self._mints.ao_eri()

    def _compute_ao_eri_df(self):
        basis = self._wfn.basisset()
        zero_bas = psi4.core.BasisSet.zero_ao_basis_set()
        aux_bas_name = psi4.core.get_global_option('DF_BASIS_SCF').lower()
        aux_bas = psi4.core.BasisSet.build(self._wfn.molecule(), "ORBITAL",
                                           aux_bas_name)

        # build three-center integral <fitbasis|lo, lo>
        self._df_pmn = np.asarray(self._mints.ao_eri(aux_bas, zero_bas,
                                                     basis, basis))
        self._df_pmn = np.squeeze(self._df_pmn)

        # build density fitting Vpq inverse
        df_Vpq = np.asarray(self._mints.ao_eri(
            aux_bas, zero_bas, aux_bas, zero_bas))
        df_Vpq = np.squeeze(df_Vpq)
        self._df_Vpq_inv = np.linalg.inv(df_Vpq)

    def compute(self):
        """
        Compute the density matrix based on current Cocc.
        """
        self._D = [np.einsum('i,ui,vi->uv', np.asarray(self._occ_val[s]),
                             self._Cocc[s], self._Cocc[s], optimize=True)
                   for s in range(self._nspin)]

    def J(self):
        """
        Returns
        -------
        out: [np.array, ...], list of np.array
            The J matrices. The size of list is the same as the input Cocc.
        """
        if self._algo == 'DIRECT':
            J = [np.einsum('pqrs,rs->pq', self._uvst, self._D[s], optimize=True)
                 for s in range(self._nspin)]
        elif self._algo == 'DF':
            J = [np.einsum('pqm,mn,nrs,rs->pq', self._df_pmn, self._df_Vpq_inv,
                           self._df_pmn, self._D[s], optimize=True)
                 for s in range(self._nspin)]
        else:
            raise Exception(f"SCF type not support: {self._algo}")
        return J

    def K(self):
        """
        Returns
        -------
        out: [np.array, ...], list of np.array
            The K matrices. The size of list is the same as the input Cocc.
        """
        if self._algo == 'DIRECT':
            K = [np.einsum('prqs,rs->pq', self._uvst, self._D[s], optimize=True)
                 for s in range(self._nspin)]
        elif self._algo == 'DF':
            K = [np.einsum('prm,mn,nqs,rs->pq', self._df_pmn, self._df_Vpq_inv,
                           self._df_pmn, self._D[s], optimize=True)
                 for s in range(self._nspin)]
        else:
            raise Exception(f"SCF type not support: {self._algo}")
        return K
