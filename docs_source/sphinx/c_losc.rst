==============
LOSC C library
==============

Introduction
------------

The LOSC library provides a C interface to perform calculations of LOSC,
including the construction of LOSC curvature matrix, LOSC localization,
local occupation number and LOSC corrections. The C interface is developed
by wrapping the :ref:`C++ interface <losc:losc c++ library>` in C-style.

.. Warning:: Not sure about how to let C code handle the exceptions that may
   be thrown in the C++ code.

The provided C interface (C header files) includes

    ================================  ======================================
    Header                              Description
    ================================  ======================================
    ``<losc/matrix.h>``               provide the C interface for the matrix
                                      objects that are used in the LOSC library.
    ``<losc/curvature.h>``            for the LOSC curvature matrix.
    ``<losc/localization.h>``         for the LOSC localization.
    ``<losc/local_occupation.h>``     for the LOSC local occupation.
    ``<losc/correction.h>``           for the LOSC corrections.
    ``<losc/losc.h>``                 A convenient header that includes
                                      everything you need.
    ================================  ======================================

In order to use the C interface of the LOSC library, following the
:ref:`steps of installation <installation:installation>`. Then all you
need to do is to include the corresponding headers in your project
to use the LOSC C library.


Data Structure
--------------
The LOSC C library defines matrix object with ``struct losc_matrix`` that is
declared in ``<losc/matrix.h>`` header. You can call the constructor function
``losc_matrix_create`` to create a ``struct losc_matrix`` with any size,
or you can pass a C pointer ``double *`` that points to an array to generate a
``struct losc_matrix``, which does not own the data. Remember the storage of
the matrix is ALWAYS interpreted in row-wise (C-style) major. See
``<losc/matrix.h>`` for more details.

Detailed References for the API
-------------------------------

Please refer `this section <./doxygen/c_losc/html/index.html>`_ for the full
documentation of C interface of the LOSC library.
