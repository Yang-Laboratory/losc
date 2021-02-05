#ifndef _LOSC_INCLUDE_LOSC_EIGEN_DEF_HPP_
#define _LOSC_INCLUDE_LOSC_EIGEN_DEF_HPP_

#include <Eigen/Dense>

namespace losc {

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::Ref;
using Eigen::VectorXd;
using ConstRefMat = const Ref<const MatrixXd>;
using ConstRefVec = const Ref<const VectorXd>;
using RefConstMat = Ref<const MatrixXd>;
using RefConstVec = Ref<const VectorXd>;
using RefMat = Ref<MatrixXd>;
using RefVec = Ref<VectorXd>;

} // namespace losc

#endif
