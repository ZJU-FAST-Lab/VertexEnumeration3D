#ifndef GEOUTILS_HPP
#define GEOUTILS_HPP

#include "quickhull.hpp"
#include "sdlp.hpp"

#include <Eigen/Eigen>

#include <cfloat>
#include <cstdint>
#include <set>
#include <chrono>

namespace geoutils
{

    // Each row of hPoly is defined by h0, h1, h2, h3 as
    // h0*x + h1*y + h2*z + h3 <= 0
    inline bool findInterior(const Eigen::MatrixX4d &hPoly,
                             Eigen::Vector3d &interior)
    {
        const int m = hPoly.rows();

        Eigen::MatrixX4d A(m, 4);
        Eigen::VectorXd b(m);
        Eigen::Vector4d c, x;
        const Eigen::ArrayXd hNorm = hPoly.leftCols<3>().rowwise().norm();
        A.leftCols<3>() = hPoly.leftCols<3>().array().colwise() / hNorm;
        A.rightCols<1>().setConstant(1.0);
        b = -hPoly.rightCols<1>().array() / hNorm;
        c.setZero();
        c(3) = -1.0;

        double minmaxsd = sdlp::linprog<4>(c, A, b, x);
        interior = x.head<3>();

        return minmaxsd < 0.0 && !std::isinf(minmaxsd);
    }

    struct filterLess
    {
        inline bool operator()(const Eigen::Vector3d &l,
                               const Eigen::Vector3d &r)
        {
            return l(0) < r(0) ||
                   (l(0) == r(0) &&
                    (l(1) < r(1) ||
                     (l(1) == r(1) &&
                      l(2) < r(2))));
        }
    };

    inline void filterVs(const Eigen::Matrix3Xd &rV,
                         const double &epsilon,
                         Eigen::Matrix3Xd &fV)
    {
        double mag = std::max(fabs(rV.maxCoeff()), fabs(rV.minCoeff()));
        double res = mag * std::max(fabs(epsilon) / mag, DBL_EPSILON);
        std::set<Eigen::Vector3d, filterLess> filter;
        fV = rV;
        int offset = 0;
        Eigen::Vector3d quanti;
        for (int i = 0; i < rV.cols(); i++)
        {
            quanti = (rV.col(i) / res).array().round();
            if (filter.find(quanti) == filter.end())
            {
                filter.insert(quanti);
                fV.col(offset) = rV.col(i);
                offset++;
            }
        }
        fV = fV.leftCols(offset).eval();
        return;
    }

    // Each row of hPoly is defined by h0, h1, h2, h3 as
    // h0*x + h1*y + h2*z + h3 <= 0
    // proposed epsilon is 1.0e-6
    inline void enumerateVs(const Eigen::MatrixX4d &hPoly,
                            const Eigen::Vector3d &inner,
                            Eigen::Matrix3Xd &vPoly,
                            const double epsilon = 1.0e-6)
    {
        Eigen::VectorXd b = -hPoly.rightCols<1>() - hPoly.leftCols<3>() * inner;
        Eigen::Matrix<double, 3, -1, Eigen::ColMajor> A =
            (hPoly.leftCols<3>().array().colwise() / b.array()).transpose();

        quickhull::QuickHull<double> qh;
        double qhullEps = std::min(epsilon, quickhull::defaultEps<double>());
        // CCW is false because the normal in quickhull towards interior
        const auto cvxHull = qh.getConvexHull(A.data(), A.cols(), false, true, qhullEps);
        const auto &idBuffer = cvxHull.getIndexBuffer();
        int hNum = idBuffer.size() / 3;
        Eigen::Matrix3Xd rV(3, hNum);
        Eigen::Vector3d normal, point, edge0, edge1;
        for (int i = 0; i < hNum; i++)
        {
            point = A.col(idBuffer[3 * i + 1]);
            edge0 = point - A.col(idBuffer[3 * i]);
            edge1 = A.col(idBuffer[3 * i + 2]) - point;
            normal = edge0.cross(edge1); //cross in CW gives an outter normal
            rV.col(i) = normal / normal.dot(point);
        }
        filterVs(rV, epsilon, vPoly);
        vPoly = (vPoly.array().colwise() + inner.array()).eval();
        return;
    }

    // Each row of hPoly is defined by h0, h1, h2, h3 as
    // h0*x + h1*y + h2*z + h3 <= 0
    // proposed epsilon is 1.0e-6
    inline bool enumerateVs(const Eigen::MatrixX4d &hPoly,
                            Eigen::Matrix3Xd &vPoly,
                            const double epsilon = 1.0e-6)
    {
        Eigen::Vector3d inner;
        if (findInterior(hPoly, inner))
        {
            enumerateVs(hPoly, inner, vPoly, epsilon);
            return true;
        }
        else
        {
            return false;
        }
    }

} // namespace geoutils

#endif
