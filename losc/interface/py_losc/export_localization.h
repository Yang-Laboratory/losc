#ifndef __LOSC_INTERFACE_PY_LOSC_EXPORT_LOCALIZATION_H__
#define __LOSC_INTERFACE_PY_LOSC_EXPORT_LOCALIZATION_H__

#include <pybind11/pybind11.h>
namespace py = pybind11;

void export_localization_base(py::module &m);
void export_localization_v2(py::module &m);

#endif
