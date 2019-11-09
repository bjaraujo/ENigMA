// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include <Eigen/Sparse>

namespace ENigMA {

namespace sle {

    enum EMatrixType {
        MT_UNKNOWN = -1,
        MT_SPARSE_SYMMETRIC,
        MT_SPARSE,
        MT_DENSE
    };

    template <typename Real>
    class CSleSystem {
    public:
        CSleSystem();
        ~CSleSystem();

        EMatrixType matrixType;

        Eigen::SparseMatrix<Real> matrixA;
        Eigen::Matrix<Real, Eigen::Dynamic, 1> vectorB;

        Eigen::Matrix<Real, Eigen::Dynamic, 1> solve();

        CSleSystem<Real>& operator=(const Real right);
        CSleSystem<Real>& operator=(Eigen::Matrix<Real, Eigen::Dynamic, 1> right);
    };

    template <typename Real>
    CSleSystem<Real> operator-(const CSleSystem<Real>& left, const Real right);

    template <typename Real>
    CSleSystem<Real> operator+(const CSleSystem<Real>& left, const Real right);

    template <typename Real>
    CSleSystem<Real> operator+(const CSleSystem<Real>& left, const CSleSystem<Real>& right);

    template <typename Real>
    CSleSystem<Real> operator-(const CSleSystem<Real>& left, const CSleSystem<Real>& right);

    template <typename Real>
    CSleSystem<Real> operator*(const Real left, const CSleSystem<Real>& right);

    template <typename Real>
    CSleSystem<Real> operator*(const Eigen::Matrix<Real, Eigen::Dynamic, 1>& left, const CSleSystem<Real>& right);
}
}

#include "SleSystem_Imp.hpp"
