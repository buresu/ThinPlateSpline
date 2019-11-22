#pragma once

#include <Eigen/Core>
#include <Eigen/StdVector>

class ThinPlateSpline {
public:
  typedef std::vector<Eigen::Vector3d,
                      Eigen::aligned_allocator<Eigen::Vector3d>>
      PointList;

  ThinPlateSpline() {}
  ThinPlateSpline(const PointList &src, const PointList &dst);
  ~ThinPlateSpline() {}

  /* Solve */
  void solve();

  /* Interpolate */
  Eigen::Vector3d interpolate(const Eigen::Vector3d &p) const;

  /* Source Points */
  const PointList &srcPoints() const { return mSrcPoints; }

  /* Set Source Points */
  void setSrcPoints(const PointList &points) { mSrcPoints = points; }

  /* Destination Points */
  const PointList &dstPoints() const { return mDstPoints; }

  /* Set Destination Points */
  void setDstPoints(const PointList &points) { mDstPoints = points; }

protected:
  /* Radial Basis Function */
  static inline double radialBasis(double r) {
    return r == 0.0 ? r : r * r * log(r);
  }

  /* Data */
  PointList mSrcPoints;
  PointList mDstPoints;
  Eigen::MatrixXd mW;
  Eigen::MatrixXd mL;
};
