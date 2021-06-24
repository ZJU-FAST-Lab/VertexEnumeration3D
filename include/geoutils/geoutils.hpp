#ifndef GEOUTILS_HPP
#define GEOUTILS_HPP

#include "quickhull.hpp"
#include "lbfgs.hpp"

#include <Eigen/Eigen>

#include <cfloat>
#include <cstdint>
#include <set>
#include <chrono>

namespace geoutils
{

inline double objectiveCVXLP(void *data,
                             const double *x,
                             double *grad,
                             const int n)
{
    const int *pM = (int *)data;
    const double *pScale = (double *)(pM + 1);
    const double *pAb = pScale + 1;

    const int M = *pM;
    const double scale = *pScale;
    Eigen::Map<const Eigen::MatrixXd> Ab(pAb, M, 4);
    Eigen::Map<const Eigen::Vector3d> p(x);
    Eigen::Map<Eigen::Vector3d> g(grad);

    double cost = 0;
    g.setZero();
    double d, denSqrt;
    Eigen::VectorXd consViola = (Ab.leftCols<3>() * p - Ab.col(3)) * scale;
    for (int i = 0; i < M; i++)
    {
        d = consViola(i);

        if (d > 0.0)
        {
            cost += ((0.5 * d + 1.0) * d + 1.0);
            g += (d + 1.0) * scale * Ab.block<1, 3>(i, 0).transpose();
        }
        else
        {
            cost += 1.0 / ((0.5 * d - 1.0) * d + 1.0);
            denSqrt = (0.5 * d - 1.0) * d + 1.0;
            g += (1.0 - d) / (denSqrt * denSqrt) * scale * Ab.block<1, 3>(i, 0).transpose();
        }
    }

    return cost;
}

// Each col of hPoly denotes a facet (outter_normal^T,point^T)^T
// The outter_normal is assumed to be NORMALIZED
// The proposed resolution is 0.001
inline bool findInterior(const Eigen::MatrixXd &hPoly,
                         Eigen::Vector3d &interior,
                         const double resolution = 0.001)
{

    uint8_t *data = new uint8_t[sizeof(int) + (1 + 4 * hPoly.cols()) * sizeof(double)];
    int *pM = (int *)data;
    double *pScale = (double *)(pM + 1);
    double *pAb = pScale + 1;

    *pM = hPoly.cols();
    // The scale should be strictly larger than (sqrt(2M-1)-1)/resolution
    *pScale = sqrt((*pM > 2 ? *pM : 3) * 2.0) / (resolution < 1.0 ? resolution : 1.0);
    Eigen::Map<Eigen::MatrixXd> Ab(pAb, *pM, 4);
    Ab.leftCols<3>() = hPoly.topRows<3>().transpose();
    Ab.col(3) = hPoly.topRows<3>().cwiseProduct(hPoly.bottomRows<3>()).colwise().sum().transpose();

    Eigen::Vector3d x(0.0, 0.0, 0.0);
    double minCost;
    lbfgs::lbfgs_parameter_t cvxlp_params;
    lbfgs::lbfgs_load_default_parameters(&cvxlp_params);
    cvxlp_params.g_epsilon = FLT_EPSILON;
    cvxlp_params.mem_size = 9;
    cvxlp_params.max_iterations = 64;

    lbfgs::lbfgs_optimize(3,
                          x.data(),
                          &minCost,
                          &objectiveCVXLP,
                          nullptr,
                          nullptr,
                          data,
                          &cvxlp_params);

    interior = x;
    bool interiorFound = (Ab.leftCols<3>() * interior - Ab.col(3)).maxCoeff() < 0;

    delete[] data;

    return interiorFound;
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

inline void filterVs(const Eigen::MatrixXd &rV,
                     const double &epsilon,
                     Eigen::MatrixXd &fV)
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

// Each col of hPoly denotes a facet (outter_normal^T,point^T)^T
// The outter_normal is assumed to be NORMALIZED
// proposed epsilon is 1.0e-6
inline void enumerateVs(const Eigen::MatrixXd &hPoly,
                        const Eigen::Vector3d &inner,
                        Eigen::MatrixXd &vPoly,
                        const double epsilon = 1.0e-6)
{
    Eigen::RowVectorXd b = hPoly.topRows<3>().cwiseProduct(hPoly.bottomRows<3>()).colwise().sum() -
                           inner.transpose() * hPoly.topRows<3>();
    Eigen::MatrixXd A = hPoly.topRows<3>().array().rowwise() / b.array();

    quickhull::QuickHull<double> qh;
    double qhullEps = std::min(epsilon, quickhull::defaultEps<double>());
    // CCW is false because the normal in quickhull towards interior
    const auto cvxHull = qh.getConvexHull(A.data(), A.cols(), false, true, qhullEps);
    const auto &idBuffer = cvxHull.getIndexBuffer();
    int hNum = idBuffer.size() / 3;
    Eigen::MatrixXd rV(3, hNum);
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

// Each col of hPoly denotes a facet (outter_normal^T,point^T)^T
// The outter_normal is assumed to be NORMALIZED
// proposed resolution is 0.001 for iterior finding
// proposed epsilon is 1.0e-6
inline bool enumerateVs(const Eigen::MatrixXd &hPoly,
                        Eigen::MatrixXd &vPoly,
                        const double resolution = 0.001,
                        const double epsilon = 1.0e-6)
{
    Eigen::Vector3d inner;
    if (findInterior(hPoly, inner, resolution))
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