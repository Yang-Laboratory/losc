#ifndef _LOSC_SRC_PY_LOSC_EXPORT_CURVATURE_HPP_
#define _LOSC_SRC_PY_LOSC_EXPORT_CURVATURE_HPP_

#include <pybind11/pybind11.h>
namespace py = pybind11;

void export_curvature_base(py::module &m);
void export_curvature_v1(py::module &m);
void export_curvature_v2(py::module &m);

#endif
