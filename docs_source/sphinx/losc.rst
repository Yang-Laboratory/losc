================
LOSC C++ library
================

Introduction
------------

The LOSC library provides a C++ interface to perform calculations of LOSC,
including the construction of LOSC curvature matrix, LOSC localization,
local occupation number and LOSC corrections.

The provided C++ interface (C++ header files) includes

    ==================================  ======================================
    Header                              Description
    ==================================  ======================================
    ``<losc/curvature.hpp>``            for the LOSC curvature matrix.
    ``<losc/localization.hpp>``         for the LOSC localization.
    ``<losc/local_occupation.hpp>``     for the LOSC local occupation.
    ``<losc/correction.hpp>``           for the LOSC corrections.
    ``<losc/exception.hpp>``            for the exceptions that may be thrown
                                        in LOSC library.
    ``<losc/eigen_def.hpp>``            for the declarations of matrix used
                                        in LOSC library.
    ``<losc/losc.hpp>``                 A convenient header that includes
                                        everything you need.
    ==================================  ======================================

In order to use the C++ interface of the LOSC library, following the
:ref:`steps of installation <installation:installation>`. Then all you
need to do is to include the corresponding headers in your project
to use the LOSC C++ library.

Data Structure
--------------
The main data structure in the LOSC C++ library is the representation of
matrix/vector objects. In the LOSC C++ library, we use the
popular `Eigen library <https://eigen.tuxfamily.org/dox/>`_,
which lets us achieve the manipulation of matrices and vectors easily.
Particularly, we use ``Eigen::Matrix`` to represent matrices, and
``Eigen::Vector`` to represent vectors.

Being different to the default behavior in ``Eigen`` that the storage of matrix
is column-wise (Fortran-style), the storage of matrix in the LOSC C++ library
is row-wise (C-style). To avoid declaring the row-wise ``Eigen::Matrix``
in a cumbersome way being like

.. code-block:: C++

   Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

an alias ``LOSCMatrix`` is defined in ``<losc/eigen_def.hpp>`` as

.. code-block:: C++

   using LOSCMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

There is a similar alias ``LOSCVector`` for the vector objects. For other
alias of ``Eigen`` types that are used the LOSC C++ library, please
refer to ``<losc/eigen_def.hpp>``.

Detailed References for the API
-------------------------------

Please refer `this section <./doxygen/losc/html/index.html>`_ for the full
documentation of C++ interface of the LOSC library.
