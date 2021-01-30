#include "../../src/localization.h"
#include <Eigen/Dense>
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
    vector<MatrixXd> lo_U(const string &guess = "identity") const override
    {
        PYBIND11_OVERRIDE_PURE(vector<MatrixXd>, /* Return type */
                               Base,             /* Parent class */
                               lo_U,             /* Name of function */
                               guess             /* Arguments */
        );
    }
    vector<MatrixXd> lo_U(ConstRefMat &U_guess,
                          double threshold = 1e-8) const override
    {
        PYBIND11_OVERRIDE_PURE(vector<MatrixXd>, /* Return type */
                               Base,             /* Parent class */
                               lo_U,             /* Name of function */
                               U_guess,          /* Arguments */
                               threshold         /* Arguments */
        );
    }
    void lo_U(RefMat L, RefMat U) const override
    {
        PYBIND11_OVERRIDE_PURE(void, /* Return type */
                               Base, /* Parent class */
                               lo_U, /* Name of function */
                               L,    /* Arguments */
                               U     /* Arguments */
        );
    }
    MatrixXd lo(const string &guess = "identity") const override
    {
        PYBIND11_OVERRIDE(MatrixXd, /* Return type */
                          Base,     /* Parent class */
                          lo,       /* Name of function */
                          guess     /* Arguments */
        );
    }
    MatrixXd lo(ConstRefMat &U_guess, double threshold = 1e-8) const override
    {
        PYBIND11_OVERRIDE(MatrixXd, /* Return type */
                          Base,     /* Parent class */
                          lo,       /* Name of function */
                          U_guess   /* Arguments */
        );
    }
};

template <class Derive_L1> class PyDerive_L1 : public PyBase<Derive_L1> {
    // Derive_L1 is the level 1 derived class from base.
    // Level 1 derived class includes: LoscLocalizerV2 ...
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
    vector<MatrixXd> lo_U(const string &guess = "identity") const override
    {
        PYBIND11_OVERRIDE(vector<MatrixXd>, /* Return type */
                          Derive_L1,        /* Parent class */
                          lo_U,             /* Name of function */
                          guess             /* Arguments */
        );
    }
    vector<MatrixXd> lo_U(ConstRefMat &U_guess,
                          double threshold = 1e-8) const override
    {
        PYBIND11_OVERRIDE(vector<MatrixXd>, /* Return type */
                          Derive_L1,        /* Parent class */
                          lo_U,             /* Name of function */
                          U_guess,          /* Arguments */
                          threshold         /* Arguments */
        );
    }
    void lo_U(RefMat L, RefMat U) const override
    {
        PYBIND11_OVERRIDE(void,      /* Return type */
                          Derive_L1, /* Parent class */
                          lo_U,      /* Name of function */
                          L,         /* Arguments */
                          U          /* Arguments */
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
        max_iter: int
            The max number of iterations.

        Returns
        -------
        out: None
        )pddoc")
        // set_convergence
        .def("set_convergence", &LocalizerBase::set_convergence, R"pddoc(
        Set the convergence tolerance for localization.

        Parameters
        ----------
        tol: double
            The convergence tolerance.

        Returns
        -------
        out: None
        )pddoc")
        // set_random_permutation
        .def("set_random_permutation", &LocalizerBase::set_random_permutation,
             R"pddoc(
        Set the flag for performing random permutation or not in localization
        with Jacobi-Sweep algorithm.

        Parameters
        ----------
        flag: bool

        Returns
        -------
        out: None
        )pddoc")
        // lo_U: overload 1
        .def("lo_U",
             static_cast<vector<MatrixXd> (LocalizerBase::*)(const string &)
                             const>(&LocalizerBase::lo_U),
             R"pddoc(
        Calculate the LOs and the unitary transformation matrix.

        Parameters
        ----------
        guess: str
            Initial guess of the unitary matrix to do localization. The choices
            are ['identity', 'random', 'random_fixed_seed']. Default to 'identity'.
            'identity': initial U matrix is set as an identity matrix.
            'random': initial U matrix is set as a random unitary matrix.
            'random_fixed_seed': initial U matrix is set as a random unitary
                matrix with fixed random seed.

        Returns
        -------
        out: (np.ndarray, np.ndarray)
            The first one is the LO coefficient matrix, and the second one is the
            corresponding U matrix.
        )pddoc")
        // lo_U: overload 2
        .def("lo_U",
             static_cast<vector<MatrixXd> (LocalizerBase::*)(
                 ConstRefMat &, double) const>(&LocalizerBase::lo_U),
             R"pddoc(
        Calculate the LOs and the unitary transformation matrix with a given
        U matrix as the initial guess.

        Parameters
        ----------
        U_guess: np.ndarray
            The initial guess of U matrix. Its data will be copied for
            localization. Its unitarity will be verified and throw an exception
            if the validation fails.
        threshold: float
            The threshold used to check the unitarity. Default to 1e-8.

        Returns
        -------
        out: (np.ndarray, np.ndarray)
            The first one is the LO coefficient matrix, and the second one is the
            corresponding U matrix.
        )pddoc")
        // lo: overload 1
        .def("lo",
             static_cast<MatrixXd (LocalizerBase::*)(const string &) const>(
                 &LocalizerBase::lo),
             R"pddoc(
        Calculate the LOs' coefficient matrix under AO.

        Returns
        -------
        out: np.ndarray
            The the LO coefficient matrix.

        Notes
        -----
        Same interface of arguments as lo_U.

        See Also
        --------
        lo_U: Return both LOs and the U matrix.
        )pddoc")
        // lo: overload 2
        .def("lo",
             static_cast<MatrixXd (LocalizerBase::*)(ConstRefMat &, double)
                             const>(&LocalizerBase::lo),
             R"pddoc(
        Calculate the LOs' coefficient matrix under AO with a given U matrix as
        the initial guess.

        Returns
        -------
        out: np.ndarray
            The the LO coefficient matrix.

        Notes
        -----
        Same interface of arguments as lo_U.

        See Also
        --------
        lo_U: return both LOs and the U matrix.
        )pddoc");
}

void export_localization_v2(py::module &m)
{
    /* losc::LoscLocalizerV2 */
    py::class_<LoscLocalizerV2, LocalizerBase, PyDerive_L1<LoscLocalizerV2>>(
        m, "LoscLocalizerV2", "LOSC localization version 2", py::dynamic_attr())
        // constructor
        .def(py::init<ConstRefMat &,              // C_lo_basis
                      ConstRefMat &,              // H_ao
                      const vector<RefConstMat> & // D_ao
                      >(),
             "C_lo_basis"_a.noconvert(), "H_ao"_a.noconvert(),
             "D_ao"_a.noconvert());
    // LoscLocalizerV2 class has no new functions compared to LocalizerBase.
    // So no more functions can be exported here.
}
