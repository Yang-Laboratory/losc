#include <Eigen/Dense>
#include <losc/localization.hpp>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // to bind std::vector to python list.
#include <vector>

namespace py = pybind11;
using namespace losc;
using namespace pybind11::literals;
using namespace Eigen;

template <class Base> class PyBase : public Base {
    // Here, Base class is the LocalizerBase class.
  public:
    using Base::Base; // Inherit constructors.

    // We should list override ALL the virtual functions.
    // Here, all the virtual functions means virtual functions in LocalizerBase
    // class.
    // It looks like virtual deconstructor can be ignored.
    vector<LOSCMatrix> lo_U(const string &guess = "identity") override
    {
        PYBIND11_OVERRIDE_PURE(vector<LOSCMatrix>, /* Return type */
                               Base,               /* Parent class */
                               lo_U,               /* Name of function */
                               guess               /* Arguments */
        );
    }
    vector<LOSCMatrix> lo_U(ConstRefMat &U_guess,
                            double threshold = 1e-8) override
    {
        PYBIND11_OVERRIDE_PURE(vector<LOSCMatrix>, /* Return type */
                               Base,               /* Parent class */
                               lo_U,               /* Name of function */
                               U_guess,            /* Arguments */
                               threshold           /* Arguments */
        );
    }
    void C_API_lo_U(RefMat L, RefMat U) override
    {
        PYBIND11_OVERRIDE_PURE(void,       /* Return type */
                               Base,       /* Parent class */
                               C_API_lo_U, /* Name of function */
                               L,          /* Arguments */
                               U           /* Arguments */
        );
    }
    LOSCMatrix lo(const string &guess = "identity") override
    {
        PYBIND11_OVERRIDE(LOSCMatrix, /* Return type */
                          Base,       /* Parent class */
                          lo,         /* Name of function */
                          guess       /* Arguments */
        );
    }
    LOSCMatrix lo(ConstRefMat &U_guess, double threshold = 1e-8) override
    {
        PYBIND11_OVERRIDE(LOSCMatrix, /* Return type */
                          Base,       /* Parent class */
                          lo,         /* Name of function */
                          U_guess     /* Arguments */
        );
    }
    double cost_func(ConstRefMat &lo) const override
    {
        PYBIND11_OVERRIDE_PURE(double,    /* Return type */
                               Base,      /* Parent class */
                               cost_func, /* Name of function */
                               lo         /* Arguments */
        );
    }
};

template <class Derive_L1> class PyDerive_L1 : public PyBase<Derive_L1> {
    // Derive_L1 is the level 1 derived class from base.
    // Level 1 derived class includes: LocalizerV2 ...
  public:
    using PyBase<Derive_L1>::PyBase; // Inherit constructors.

    // We should list override all the virtual functions.
    // Here, all the virtual functions means virtual functions in
    // Derive_L1 class.
    //
    // Here is the trick. PyDerive_L1 is publicly inherited from
    // PyBase<Derive_L1>. PyBase<Derive_L1> class generates all the virtual
    // functions in Derived_L1 class, which are inherited from base.
    // Such public inheriting from PyBase<Derived_L1> helps us to avoid
    // duplicated codes to override these virtual functions that are already
    // overode in PyBase<Derive_L1> class.
    // Now, the only thing left is that we can just look at the Derive_L1
    // class, and override ALL the virtual functions showing in Derive_L1 class.
    vector<LOSCMatrix> lo_U(const string &guess = "identity") override
    {
        PYBIND11_OVERRIDE(vector<LOSCMatrix>, /* Return type */
                          Derive_L1,          /* Parent class */
                          lo_U,               /* Name of function */
                          guess               /* Arguments */
        );
    }
    vector<LOSCMatrix> lo_U(ConstRefMat &U_guess,
                            double threshold = 1e-8) override
    {
        PYBIND11_OVERRIDE(vector<LOSCMatrix>, /* Return type */
                          Derive_L1,          /* Parent class */
                          lo_U,               /* Name of function */
                          U_guess,            /* Arguments */
                          threshold           /* Arguments */
        );
    }
    void C_API_lo_U(RefMat L, RefMat U) override
    {
        PYBIND11_OVERRIDE(void,       /* Return type */
                          Derive_L1,  /* Parent class */
                          C_API_lo_U, /* Name of function */
                          L,          /* Arguments */
                          U           /* Arguments */
        );
    }
    double cost_func(ConstRefMat &lo) const override
    {
        PYBIND11_OVERRIDE(double,    /* Return type */
                          Derive_L1, /* Parent class */
                          cost_func, /* Name of function */
                          lo         /* Arguments */
        );
    }
};

void export_localization_base(py::module &m)
{
    /* losc::LocalizerBase */
    py::class_<LocalizerBase, PyBase<LocalizerBase>>(
        m, "LocalizerBase", "LOSC localization base", py::dynamic_attr())
        // constructor
        .def(py::init<ConstRefMat & // C_lo_basis
                      >(),
             "C_lo_basis"_a.noconvert())
        // set_max_iter
        .def("set_max_iter", &LocalizerBase::set_max_iter, R"pddoc(
        Set the maximum iteration number for localization.

        Parameters
        ----------
        max_iter : int
            The maximum number of iterations.

        Returns
        -------
        None
        )pddoc")
        // set_convergence
        .def("set_convergence", &LocalizerBase::set_convergence, R"pddoc(
        Set the convergence tolerance for localization.

        Parameters
        ----------
        tol : float
            The convergence tolerance.

        Returns
        -------
        None
        )pddoc")
        // set_random_permutation
        .def("set_random_permutation", &LocalizerBase::set_random_permutation,
             R"pddoc(
        Set the flag to perform random permutation or not in the localization
        with Jacobi-Sweep algorithm.

        Parameters
        ----------
        flag : bool

        Returns
        -------
        None
        )pddoc")
        // lo_U: overload 1
        .def("lo_U",
             static_cast<vector<LOSCMatrix> (LocalizerBase::*)(const string &)>(
                 &LocalizerBase::lo_U),
             "guess"_a = "identity",
             R"pddoc(
        Calculate the LOs and the unitary transformation matrix.

        Parameters
        ----------
        guess : {'identity', 'random', 'random_fixed_seed'}, default='identity'
            Initial guess of the unitary matrix to do localization.

            - 'identity': initial U matrix is set as an identity matrix.
            - 'random': initial U matrix is set as a random unitary matrix.
            - 'random_fixed_seed': initial U matrix is set as a random unitary
              matrix with fixed random seed.

        Returns
        -------
        C_lo : numpy.array
            The LO coefficient matrix. ``C_lo[:, i]`` is the coefficients of
            the i-th LO represented on AO.
        U : numpy.array
            The corresponding unitary transformation matrix for `C_lo`.
            ``U[:, i]`` corresponds to the transformation of i-th LO.
        )pddoc")
        // lo_U: overload 2
        .def("lo_U",
             static_cast<vector<LOSCMatrix> (LocalizerBase::*)(
                 ConstRefMat &, double)>(&LocalizerBase::lo_U),
             "U_guess"_a, "threshold"_a = 1.e-8,
             R"pddoc(
        Calculate the LOs and the unitary transformation matrix with a given
        U matrix as the initial guess.

        Parameters
        ----------
        U_guess : numpy.array
            The initial guess of U matrix. Its data will be copied for
            localization. Its unitarity will be verified and throw an exception
            if the validation fails.
        threshold : float, default=1.0e-8
            The threshold used to check the unitarity.

        Returns
        -------
        C_lo : numpy.array
            The LO coefficient matrix. ``C_lo[:, i]`` is the coefficients of
            the i-th LO represented on AO.
        U : numpy.array
            The corresponding unitary transformation matrix for `C_lo`.
            ``U[:, i]`` corresponds to the transformation of i-th LO.
        )pddoc")
        // lo: overload 1
        .def("lo",
             static_cast<LOSCMatrix (LocalizerBase::*)(const string &)>(
                 &LocalizerBase::lo),
             "guess"_a = "identity",
             R"pddoc(
        Calculate the LOs' coefficient matrix under AO.

        Parameters
        ----------
        guess : {'identity', 'random', 'random_fixed_seed'}, default='identity'
            Initial guess of the unitary matrix to do localization.

            - 'identity': initial U matrix is set as an identity matrix.
            - 'random': initial U matrix is set as a random unitary matrix.
            - 'random_fixed_seed': initial U matrix is set as a random unitary
              matrix with fixed random seed.

        Returns
        -------
        C_lo : numpy.array
            The LO coefficient matrix. ``C_lo[:, i]`` is the coefficients of
            the i-th LO represented on AO.

        See Also
        --------
        lo_U
        )pddoc")
        // lo: overload 2
        .def("lo",
             static_cast<LOSCMatrix (LocalizerBase::*)(ConstRefMat &, double)>(
                 &LocalizerBase::lo),
             "U_guess"_a, "threshold"_a = 1.e-8,
             R"pddoc(
        Calculate the LOs' coefficient matrix under AO with a given U matrix as
        the initial guess.

        Parameters
        ----------
        U_guess : numpy.array
            The initial guess of U matrix. Its data will be copied for
            localization. Its unitarity will be verified and throw an exception
            if the validation fails.
        threshold : float, default=1.0e-8
            The threshold used to check the unitarity.

        Returns
        -------
        C_lo : numpy.array
            The LO coefficient matrix. ``C_lo[:, i]`` is the coefficients of
            the i-th LO represented on AO.

        See Also
        --------
        lo_U
        )pddoc")
        // cost_func
        .def("cost_func", &LocalizerBase::cost_func, "lo"_a,
             R"pddoc(
        Calculate the localization cost function value for the given LOs.

        Parameters
        ----------
        lo : numpy.array
            The LOs coefficient matrix under AOs with dimension [nbasis, nlo].

        Returns
        -------
        float
            The cost function value in hartree.
        )pddoc")
        // nsteps
        .def("steps", &LocalizerBase::steps,
             R"pddoc(
        Return the number of iteration steps for the most recent localization
        performed by calling the `lo_U` and `lo` functions.

        Returns
        -------
        int
            The number of iterations.
        )pddoc")

        // is_converged
        .def("is_converged", &LocalizerBase::is_converged,
             R"pddoc(
        Return the convergence of the most recent localization performed by
        calling the `lo_U()` and `lo()` functions.

        Returns
        -------
        bool
            The convergence of localization.
        )pddoc");
}

void export_localization_v2(py::module &m)
{
    /* losc::LocalizerV2 */
    py::class_<LocalizerV2, LocalizerBase, PyDerive_L1<LocalizerV2>>(
        m, "LocalizerV2", "LOSC localization version 2", py::dynamic_attr())
        // constructor
        .def(py::init<ConstRefMat &,              // C_lo_basis
                      ConstRefMat &,              // H_ao
                      const vector<RefConstMat> & // D_ao
                      >(),
             "C_lo_basis"_a.noconvert(), "H_ao"_a.noconvert(),
             "D_ao"_a.noconvert(),
             R"pddoc(
        Constructor of LOSC localizerV2.

        Parameters
        ----------
        C_lo_basis : numpy.array
            The coefficient matrix of LO basis that is expanded under AO with
            dimension [nbasis, nlo]. ``C_lo_basis[:, i]`` is the coefficient of
            the i-th LO basis. The LO basis refers to a set of MOs that
            will be unitary transformed to be the LOs. The basis is usually
            being the COs from the associated DFA.

        H_ao : numpy.array
            The Hamiltonian matrix under AO representation for the associated
            DFA. The dimension is [nbasis, nbasis].

        D_ao : list of numpy.array
            The dipole matrix under AO representation for x, y and z directions.
        )pddoc")
        // set_gamma
        .def("set_gamma", &LocalizerV2::set_gamma, "gamma"_a,
             R"pddoc(
        Set the localization parameter gamma.

        Returns
        -------
        None
        )pddoc")
        // set_c
        .def("set_c", &LocalizerV2::set_c, "c"_a,
             R"pddoc(
        Set the localization parameter c.

        Returns
        -------
        None
        )pddoc");
}
