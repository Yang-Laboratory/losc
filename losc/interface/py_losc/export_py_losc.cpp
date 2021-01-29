#include "export_curvature.h"
#include "export_localization.h"
#include "export_py_losc.h"

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
}
