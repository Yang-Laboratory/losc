#ifndef _LOSC_SRC_C_LOSC_MATRIX_IMPL_HPP_
#define _LOSC_SRC_C_LOSC_MATRIX_IMPL_HPP_

#include <Eigen/Dense>
#include <losc/exception.hpp>
#include <losc/eigen_def.hpp>

/**
 * A simple matrix struct that is used to communicate matrix data between
 * the users and LOSC library.
 */
class losc_matrix {
  public:
    double *data_;
    size_t row_;
    size_t col_;
    bool own_data_;

    losc_matrix(size_t row, size_t col) : row_{row}, col_{col}
    {
        data_ = new double[row_ * col_];
        own_data_ = true;
    }

    losc_matrix(size_t row, size_t col, double *data)
        : row_{row}, col_{col}, data_{data}, own_data_{false}
    {
        if (row_ != 0 && col_ != 0 && data_ == nullptr) {
            throw losc::exception::LoscException(
                "losc_matrix: detect null pointer to an non-empty matrix.");
        }
    }

    ~losc_matrix()
    {
        if (own_data_)
            delete[] data_;
    }
};

/**
 * Map losc_matrix object to a const Eigen::Ref object.
 */
inline losc::RefConstMat
losc_matrix_to_eigen_const(const losc_matrix *m)
{
    return losc::MapConstMat(m->data_, m->row_, m->col_);
}

/**
 * Map losc_matrix object to a Eigen::Ref object.
 */
inline losc::RefMat losc_matrix_to_eigen(const losc_matrix *m)
{
    return losc::MapMat(m->data_, m->row_, m->col_);
}

#endif
