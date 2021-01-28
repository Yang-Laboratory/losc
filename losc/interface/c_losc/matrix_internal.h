#ifndef __LOSC_INTERFACE_C_LOSC_MATRIX_INTERNAL_H__
#define __LOSC_INTERFACE_C_LOSC_MATRIX_INTERNAL_H__

#include "../../src/exception.h"
#include <Eigen/Dense>

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
Eigen::Ref<const Eigen::MatrixXd> losc_matrix_to_eigen_const(const losc_matrix *m)
{
    return Eigen::Map<const Eigen::MatrixXd>(m->data_, m->row_, m->col_);
}

/**
 * Map losc_matrix object to a Eigen::Ref object.
 */
Eigen::Ref<Eigen::MatrixXd> losc_matrix_to_eigen(const losc_matrix *m)
{
    return Eigen::Map<Eigen::MatrixXd>(m->data_, m->row_, m->col_);
}

#endif
