/**
 * @file export_py_losc.cpp
 * @brief Binding C++ LOSC library with Python with using pybind11 library.
 *
 * This is the main file in which the binding code is implemented. It merges
 * all the branches of binding code (for curvature, localization, ...) into one.
 * For details of binding code, see other `export_xxx.cpp` file. At the end,
 * these `export_xxx.cpp` files will be compiled to give a binary library
 * which serves as a customized C++ extension module to be used in Python.
 *
 * @note
 * 1. We don't want to the Python user to directly import the binary C++
 * extension module, although it can be successfully imported. Instead,
 * we provide a Python module that further wraps the C++ extension module.
 * See `py_losc.py` for details. These python wrapping modules provide
 * better and much readable documentations. In addition, they also provide
 * some utility functions, which enriches the core classes and functions
 * are defined in C++ extension, to have better and easier usage in Python
 * for LOSC calculations. As a whole, the Python users only need to
 * import `py_losc.py` in their python code to use the LOSC library.
 *
 * For developers:
 * 1. Functions in `export_xxx.cpp` files that have documentation strings is the
 * final version of the exporting.
 * 2. Functions in `export_xxx.cpp` files that have NO documentation strings
 * may be further wrapped in python module to produce better implementation.
 * For example, CurvatureV2 class (its constructor) is further wrapped in
 * python module to enable the implicit conversion of storage order between
 * np.ndarray and Eigen::MatrixXd.
 */

#include "export_py_losc.hpp"
#include "export_curvature.hpp"
#include "export_localization.hpp"
#include "export_correction.hpp"
#include "export_local_occupation.hpp"

PYBIND11_MODULE(PY_LOSC_MODULE_NAME, m)
{
    m.doc() = "Localizer Orbital Scaling Correction (LOSC) Library";

    // curvature
    export_curvature_base(m);
    export_curvature_v1(m);
    export_curvature_v2(m);

    // localization
    export_localization_base(m);
    export_localization_v2(m);

    // correction
    export_correction(m);

    // local_occ
    export_local_occupation(m);
}
