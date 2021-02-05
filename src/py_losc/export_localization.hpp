#ifndef _LOSC_SRC_PY_LOSC_EXPORT_LOCALIZATION_HPP_
#define _LOSC_SRC_PY_LOSC_EXPORT_LOCALIZATION_HPP_

#include <pybind11/pybind11.h>
namespace py = pybind11;

void export_localization_base(py::module &m);
void export_localization_v2(py::module &m);

#endif
