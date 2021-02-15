#ifndef _LOSC_INCLUDE_LOSC_EIGEN_DEF_HPP_
#define _LOSC_INCLUDE_LOSC_EIGEN_DEF_HPP_

#include <Eigen/Dense>

namespace losc {

using Eigen::Map;
using Eigen::Ref;
using LOSCMatrix =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using LOSCVector = Eigen::VectorXd;

using ConstRefMat = const Ref<const LOSCMatrix>;
using ConstRefVec = const Ref<const LOSCVector>;
using RefConstMat = Ref<const LOSCMatrix>;
using RefConstVec = Ref<const LOSCMatrix>;
using RefMat = Ref<LOSCMatrix>;
using RefVec = Ref<LOSCVector>;

using ConstMapMat = const Map<const LOSCMatrix>;
using ConstMapVec = const Map<const LOSCVector>;
using MapConstMat = Map<const LOSCMatrix>;
using MapConstVec = Map<const LOSCVector>;
using MapMat = Map<LOSCMatrix>;
using MapVec = Map<LOSCVector>;

} // namespace losc

#endif
