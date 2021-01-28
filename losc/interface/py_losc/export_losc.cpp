#include "../../src/curvature.h"
#include <Eigen/Dense>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace losc;
using namespace pybind11::literals;
using namespace Eigen;

class PyCurvatureBase : public CurvatureBase {
  public:
    /* Inherit the constructors */
    using CurvatureBase::CurvatureBase;

    /* Trampoline (need one for each virtual function) */
    MatrixXd kappa() const override
    {
        PYBIND11_OVERRIDE_PURE(
            MatrixXd,      /* Return type */
            CurvatureBase, /* Parent class */
            kappa, /* Name of function in C++ (must match Python name) */
        );
    }
    void kappa(RefMat K) const override
    {
        PYBIND11_OVERRIDE_PURE(
            void,          /* Return type */
            CurvatureBase, /* Parent class */
            kappa, /* Name of function in C++ (must match Python name) */
        );
    }
};

PYBIND11_MODULE(py_losc_core, m)
{
    m.doc() = "Localizer Orbital Scaling Correction (LOSC) Library";

    /* losc::DFAInfo */
    py::class_<DFAInfo>(m, "DFAInfo",
                        "density functional approximation information",
                        py::dynamic_attr())
        .def(py::init<double, double, const string &>(), "gga_x_wt"_a,
             "hf_x_wt"_a, "name"_a = "",
             R"pddoc(
             Constructor of DFAInfo class which represents a DFA.

             ---------
             Parameters:
             name: str
                The description of the DFA. Default to an empty string.
             hf_x_wt: float
                The weight of HF exchange in the DFA.
             gga_x_wt: float
                The weight of GGA type exchange in the DFA.
             )pddoc");

    /* losc::CurvatureBase */
    py::class_<CurvatureBase, PyCurvatureBase>(
        m, "CurvatureBase", "Curvature base class", py::dynamic_attr())
        .def(py::init<const DFAInfo &,
                      ConstRefMat &, // C_lo
                      ConstRefMat &, // df_pii
                      ConstRefMat &, // df_vpq_inv
                      ConstRefMat &, // grid_basis_val
                      ConstRefVec &  // grid_weight
                      >())
        .def("kappa", static_cast<MatrixXd (CurvatureBase::*)() const>(
                          &CurvatureBase::kappa))
        .def("nlo", &CurvatureBase::nlo)
        .def("nbasis", &CurvatureBase::nbasis)
        .def("nfitbasis", &CurvatureBase::nfitbasis)
        .def("npts", &CurvatureBase::npts);

    /* losc::CurvatureV1*/
    py::class_<CurvatureV1, CurvatureBase>(
        m, "CurvatureV1", "Curvature version 1", py::dynamic_attr())
        .def(py::init<const DFAInfo &,
                      ConstRefMat &, // C_lo
                      ConstRefMat &, // df_pii
                      ConstRefMat &, // df_vpq_inv
                      ConstRefMat &, // grid_basis_val
                      ConstRefVec &  // grid_weight
                      >())
        .def(
            "kappa",
            static_cast<MatrixXd (CurvatureV1::*)() const>(&CurvatureV1::kappa),
            "Compute the LOSC curvature version 1 matrix.");

    /* losc::CurvatureV2*/
    py::class_<CurvatureV2, CurvatureBase>(
        m, "CurvatureV2", "Curvature version 2", py::dynamic_attr())
        .def(py::init<const DFAInfo &,
                      ConstRefMat &, // C_lo
                      ConstRefMat &, // df_pii
                      ConstRefMat &, // df_vpq_inv
                      ConstRefMat &, // grid_basis_val
                      ConstRefVec &  // grid_weight
                      >())
        .def(
            "kappa",
            static_cast<MatrixXd (CurvatureV2::*)() const>(&CurvatureV2::kappa),
            "Compute the LOSC curvature version 2 matrix.");
}
