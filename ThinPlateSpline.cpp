#include "ThinPlateSpline.hpp"
#include <Eigen/QR>

ThinPlateSpline::ThinPlateSpline(const std::vector<Eigen::Vector3d> &src,
                                 const std::vector<Eigen::Vector3d> &dst)
    : mSrcPoints(src), mDstPoints(dst) {}

void ThinPlateSpline::solve() {
  if (mSrcPoints.size() != mDstPoints.size())
    return;

  const int num(int(mSrcPoints.size()));
  const int rows(num + 3 + 1);

  // Create L Matrix
  mL = Eigen::MatrixXd::Zero(rows, rows);

  for (int i(0); i < num; ++i) {

    int j(i + 1);

    for (; j < num; ++j)
      mL(i, j) = mL(j, i) = radialBasis(
          (mSrcPoints[std::size_t(i)] - mSrcPoints[std::size_t(j)]).norm());

    mL(j, i) = mL(i, j) = 1.0;
    ++j;

    for (int posElm(0); j < rows; ++posElm, ++j)
      mL(j, i) = mL(i, j) = mSrcPoints[std::size_t(i)][posElm];
  }

  // Create Y Matrix
  Eigen::MatrixXd Y = Eigen::MatrixXd::Zero(rows, 3);

  for (int i(0); i < num; ++i)
    Y.row(i) = mDstPoints[std::size_t(i)];

  // Solve L W^T = Y as W^T = L^-1 Y
  mW = mL.colPivHouseholderQr().solve(Y);
}

Eigen::Vector3d ThinPlateSpline::interpolate(const Eigen::Vector3d &p) const {
  Eigen::Vector3d res = Eigen::Vector3d::Zero();
  int i(0);

  for (; i < mW.rows() - (3 + 1); ++i) {

    double rb = radialBasis((mSrcPoints[std::size_t(i)] - p).norm());

    res += mW.row(i) * rb;
  }

  res += mW.row(i);
  i++;

  for (int j(0); j < 3; ++j, ++i)
    res += mW.row(i) * p[j];

  return res;
}
