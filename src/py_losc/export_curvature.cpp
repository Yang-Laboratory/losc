#include <Eigen/Dense>
#include <losc/curvature.hpp>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace losc;
using namespace pybind11::literals;
using namespace Eigen;

template <class Base> class PyBase : public Base {
    // Here, Base class is the CurvatureBase class.
  public:
    using Base::Base; // Inherit constructors.

    void C_API_kappa(RefMat K) const override
    {
        PYBIND11_OVERRIDE_PURE(void,        /* Return type */
                               Base,        /* Parent class */
                               C_API_kappa, /* Name of function */
                               K            /* Arguments */
        );
    }
    MatrixXd kappa() const override
    {
        PYBIND11_OVERRIDE(MatrixXd, /* Return type */
                          Base,     /* Parent class */
                          kappa,    /* Name of function */
                                    /* This function has no argument. We need
                                     * to keep the comma as suggested by
                                     * pybind11.*/
        );
    }
};

template <class Derive_L1> class PyDerive_L1 : public PyBase<Derive_L1> {
    // Derive_L1 is the level 1 derived class from base.
    // Level 1 derived class includes: CurvatureV1, CurvatureV2, ...
  public:
    using PyBase<Derive_L1>::PyBase; // Inherit constructors.
    using PyBase<Derive_L1>::kappa; // bring other overloaded kappa functions to
                                    // avoid hiding.

    void C_API_kappa(RefMat K) const override
    {
        PYBIND11_OVERRIDE(void,        /* Return type */
                          Derive_L1,   /* Parent class */
                          C_API_kappa, /* Name of function */
                          K            /* Arguments */
        );
    }
};

void export_curvature_base(py::module &m)
{
    /* losc::DFAInfo */
    py::class_<DFAInfo>(m, "DFAInfo",
                        "density functional approximation information",
                        py::dynamic_attr())
        .def(py::init<double, double, const string &>(), "gga_x"_a, "hf_x"_a,
             "name"_a = "")
        .def("name", &DFAInfo::name,
             R"pddoc(
             Returns
             -------
             out: str
                The description of the DFA.
             )pddoc")
        .def("gga_x", &DFAInfo::gga_x,
             R"pddoc(
             Returns
             -------
             out: float
                The weight of the GGA type exchange of the DFA.
             )pddoc")
        .def("hf_x", &DFAInfo::hf_x,
             R"pddoc(
             Returns
             -------
             out: float
                The weight of HF exchange of the DFA.
             )pddoc");

    /* losc::CurvatureBase */
    py::class_<CurvatureBase, PyBase<CurvatureBase>>(
        m, "CurvatureBase", "LOSC curvature base", py::dynamic_attr())
        // constructor
        .def(py::init<const DFAInfo &, // dfa_info
                      ConstRefMat &,   // df_pii
                      ConstRefMat &,   // df_Vpq_inv
                      ConstRefMat &,   // grid_lo
                      ConstRefVec &    // grid_weight
                      >(),
             "dfa_info"_a, "df_pii"_a.noconvert(), "df_Vpq_inv"_a.noconvert(),
             "grid_lo"_a.noconvert(), "grid_weight"_a.noconvert())
        // nlo
        .def("nlo", &CurvatureBase::nlo, R"pddoc(
            Returns
            -------
            out: int
                the number of LOs.
            )pddoc")
        // nfitbasis
        .def("nfitbasis", &CurvatureBase::nfitbasis, R"pddoc(
            Returns
            -------
            out: int
                the number of fitting basis in density fitting.
            )pddoc")
        // npts
        .def("npts", &CurvatureBase::npts, R"pddoc(
            Returns
            -------
            out: int
                the number of grid points.
            )pddoc")
        // kappa: only export the version of creating a kappa matrix.
        .def("kappa",
             static_cast<MatrixXd (CurvatureBase::*)() const>(
                 &CurvatureBase::kappa),
             R"pddoc(
            Returns
            -------
            out: np.ndarray [nlo, nlo]
                The LOSC curvature matrix.
            )pddoc");
}

void export_curvature_v1(py::module &m)
{
    /* losc::CurvatureV1 */
    py::class_<CurvatureV1, CurvatureBase, PyDerive_L1<CurvatureV1>>(
        m, "CurvatureV1", "LOSC curvature version 1", py::dynamic_attr())
        // constructor
        .def(py::init<const DFAInfo &, // dfa_info
                      ConstRefMat &,   // df_pii
                      ConstRefMat &,   // df_Vpq_inv
                      ConstRefMat &,   // grid_lo
                      ConstRefVec &    // grid_weight
                      >(),
             "dfa_info"_a, "df_pii"_a.noconvert(), "df_Vpq_inv"_a.noconvert(),
             "grid_lo"_a.noconvert(), "grid_weight"_a.noconvert());
    // CurvatureV1 class has no new functions compared to CurvatureBase.
    // So no more functions need to be exported here.
}

void export_curvature_v2(py::module &m)
{
    /* losc::CurvatureV2 */
    py::class_<CurvatureV2, CurvatureBase, PyDerive_L1<CurvatureV2>>(
        m, "CurvatureV2", "LOSC curvature version 2", py::dynamic_attr())
        // constructor
        .def(py::init<const DFAInfo &, // dfa_info
                      ConstRefMat &,   // df_pii
                      ConstRefMat &,   // df_Vpq_inv
                      ConstRefMat &,   // grid_lo
                      ConstRefVec &    // grid_weight
                      >(),
             "dfa_info"_a, "df_pii"_a.noconvert(), "df_Vpq_inv"_a.noconvert(),
             "grid_lo"_a.noconvert(), "grid_weight"_a.noconvert());
    // CurvatureV2 class has no new functions compared to CurvatureBase.
    // So no more functions need to be exported here.
}
