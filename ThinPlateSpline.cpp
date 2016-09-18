#include "ThinPlateSpline.h"
#include <Eigen/QR>

ThinPlateSpline::ThinPlateSpline(const std::vector<Eigen::Vector3d> &src,
                                 const std::vector<Eigen::Vector3d> &dst)
    : mSrcPoints(src)
    , mDstPoints(dst)
{}

void ThinPlateSpline::solve()
{
    if (mSrcPoints.size() != mDstPoints.size())
        return;

    const int num(int(mSrcPoints.size()));
    const int rows(num + 3 + 1);

    // Create L Matrix
    mL = Eigen::MatrixXd::Zero(rows, rows);

    for (int i(0); i < num; ++i) {

        int j(i + 1);

        for (; j < num; ++j)
            mL(i, j) = mL(j, i) = ThinPlateSpline::radialBasis((mSrcPoints[std::size_t(i)]
                                                               - mSrcPoints[std::size_t(j)]).norm());

        mL(j, i) = mL(i, j) = 1.0;

        ++j;

        for (int posElm(0); j < rows; ++posElm, ++j)
            mL(j, i) = mL(i, j) = mSrcPoints[std::size_t(i)][posElm];
    }

    // Create Y Matrix
    Eigen::MatrixXd Y = Eigen::MatrixXd::Zero(rows, 3);

    for (int i(0); i < num; ++i)
        for (int j(0); j < 3; ++j)
            Y(i, j) = mDstPoints[std::size_t(i)][j];

    // Solve L W^T = Y as W^T = L^-1 Y
    mW = mL.colPivHouseholderQr().solve(Y);
}

Eigen::Vector3d ThinPlateSpline::interpolate(const Eigen::Vector3d &p) const
{
    Eigen::Vector3d res = Eigen::Vector3d::Zero();
    int i(0);

    for (; i < mW.rows() - (3 + 1); ++i) {

        double rb = ThinPlateSpline::radialBasis((mSrcPoints[std::size_t(i)] - p).norm());

        for (int j(0); j < 3; ++j)
            res[j] += mW(i, j) * rb;
    }

    for (int j(0); j < 3; ++j)
        res[j] += mW(i, j);

    i++;

    for (int j(0); j < 3; ++j, ++i)
        for (int k(0); k < 3; ++k)
            res[k] += mW(i, k) * p[j];

    return res;
}
