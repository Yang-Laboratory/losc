/**
 * @file
 * @brief A collection of alias to ``Eigen`` types that used in LOSC.
 */

#ifndef _LOSC_INCLUDE_LOSC_EIGEN_DEF_HPP_
#define _LOSC_INCLUDE_LOSC_EIGEN_DEF_HPP_

#include <Eigen/Dense>

namespace losc {

using Eigen::Map;
using Eigen::Ref;

/**
 * @brief Matrix object used in LOSC.
 * @note The storage order is row-wise (C-style), which is opposite to the
 * ``Eigen`` default behavior.
 */
using LOSCMatrix =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

/**
 * @brief Vector object used in LOSC.
 */
using LOSCVector = Eigen::VectorXd;

/**
 * @brief Constatnt reference to a constant LOSCMatrix object
 */
using ConstRefMat = const Ref<const LOSCMatrix>;

/**
 * @brief Constatnt reference to a constant LOSCVector object
 */
using ConstRefVec = const Ref<const LOSCVector>;

/**
 * @brief Reference to a constant LOSCMatrix object
 */
using RefConstMat = Ref<const LOSCMatrix>;

/**
 * @brief Reference to a constant LOSCVector object
 */
using RefConstVec = Ref<const LOSCMatrix>;

/**
 * @brief Reference to a LOSCMatrix object
 */
using RefMat = Ref<LOSCMatrix>;

/**
 * @brief Reference to a LOSCVector object
 */
using RefVec = Ref<LOSCVector>;

/**
 * @brief Constatnt map to a constant LOSCMatrix object
 */
using ConstMapMat = const Map<const LOSCMatrix>;

/**
 * @brief Constatnt map to a constant LOSCVetor object
 */
using ConstMapVec = const Map<const LOSCVector>;

/**
 * @brief Map to a constant LOSCMatrix object
 */
using MapConstMat = Map<const LOSCMatrix>;

/**
 * @brief Map to a constant LOSCVetor object
 */
using MapConstVec = Map<const LOSCVector>;

/**
 * @brief Map to a LOSCMatrix object
 */
using MapMat = Map<LOSCMatrix>;

/**
 * @brief Map to a LOSCVetor object
 */
using MapVec = Map<LOSCVector>;

} // namespace losc

#endif
