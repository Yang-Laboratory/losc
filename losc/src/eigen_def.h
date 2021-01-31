#ifndef __LOSC_SRC_EIGEN_DEF_H__
#define __LOSC_SRC_EIGEN_DEF_H__

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
