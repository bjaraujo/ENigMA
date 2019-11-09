// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

#ifdef USE_VIENNACL

#define VIENNACL_HAVE_EIGEN 1
//#define VIENNACL_WITH_OPENCL 1

#include "viennacl/compressed_matrix.hpp"
#include "viennacl/linalg/bicgstab.hpp"
#include "viennacl/linalg/cg.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/vector.hpp"

#endif

namespace ENigMA {

namespace sle {

    template <typename Real>
    CSleSystem<Real>::CSleSystem()
    {

        matrixType = MT_UNKNOWN;
    }

    template <typename Real>
    CSleSystem<Real>::~CSleSystem()
    {
    }

    template <typename Real>
    Eigen::Matrix<Real, Eigen::Dynamic, 1> CSleSystem<Real>::solve()
    {

        Eigen::Matrix<Real, Eigen::Dynamic, 1> x;

        if (matrixType == MT_SPARSE_SYMMETRIC) {

#ifdef USE_VIENNACL

            x = viennacl::linalg::solve(matrixA, vectorB, viennacl::linalg::cg_tag());

#else

            Eigen::ConjugateGradient<Eigen::SparseMatrix<Real>> solver;
            solver.compute(matrixA);
            x = solver.solve(vectorB);

#endif

        } else if (matrixType == MT_SPARSE) {

#ifdef USE_VIENNACL

            x = viennacl::linalg::solve(matrixA, vectorB, viennacl::linalg::bicgstab_tag());

#else

            Eigen::BiCGSTAB<Eigen::SparseMatrix<Real>> solver;
            solver.compute(matrixA);
            x = solver.solve(vectorB);

#endif

        } else if (matrixType == MT_DENSE) {

            Eigen::SparseLU<Eigen::SparseMatrix<Real, Eigen::ColMajor>, Eigen::COLAMDOrdering<Integer>> solver;

            solver.analyzePattern(matrixA);
            solver.factorize(matrixA);
            x = solver.solve(vectorB);
        }

        return x;
    }

    template <typename Real>
    CSleSystem<Real>& CSleSystem<Real>::operator=(const Real right)
    {

        vectorB.array() += right;
        return *this;
    }

    template <typename Real>
    CSleSystem<Real>& CSleSystem<Real>::operator=(Eigen::Matrix<Real, Eigen::Dynamic, 1> right)
    {

        vectorB += right;
        return *this;
    }

    template <typename Real>
    CSleSystem<Real> operator-(const CSleSystem<Real>& left, const Real right)
    {

        CSleSystem<Real> aSystem;

        aSystem.matrixType = left.matrixType;

        aSystem.matrixA = left.matrixA;
        for (int k = 0; k < aSystem.matrixA.outerSize(); ++k)
            aSystem.matrixA.coeffRef(k, k) -= right;

        aSystem.vectorB = left.vectorB;
        for (int k = 0; k < aSystem.vectorB.size(); ++k)
            aSystem.vectorB[k] -= right;

        return aSystem;
    }

    template <typename Real>
    CSleSystem<Real> operator+(const CSleSystem<Real>& left, const Real right)
    {

        CSleSystem<Real> aSystem;

        aSystem.matrixType = left.matrixType;

        aSystem.matrixA = left.matrixA;
        for (int k = 0; k < aSystem.matrixA.outerSize(); ++k)
            aSystem.matrixA.coeffRef(k, k) += right;

        aSystem.vectorB = left.vectorB;
        for (int k = 0; k < aSystem.vectorB.size(); ++k)
            aSystem.vectorB[k] += right;

        return aSystem;
    }

    template <typename Real>
    CSleSystem<Real> operator+(const CSleSystem<Real>& left, const CSleSystem<Real>& right)
    {

        CSleSystem<Real> aSystem;

        if (left.matrixType > right.matrixType)
            aSystem.matrixType = left.matrixType;
        else
            aSystem.matrixType = right.matrixType;

        aSystem.matrixA = left.matrixA + right.matrixA;
        aSystem.vectorB = left.vectorB + right.vectorB;

        return aSystem;
    }

    template <typename Real>
    CSleSystem<Real> operator-(const CSleSystem<Real>& left, const CSleSystem<Real>& right)
    {

        CSleSystem<Real> aSystem;

        if (left.matrixType > right.matrixType)
            aSystem.matrixType = left.matrixType;
        else
            aSystem.matrixType = right.matrixType;

        aSystem.matrixA = left.matrixA - right.matrixA;
        aSystem.vectorB = left.vectorB - right.vectorB;

        return aSystem;
    }

    template <typename Real>
    CSleSystem<Real> operator*(const Real left, const CSleSystem<Real>& right)
    {

        CSleSystem<Real> aSystem;

        aSystem.matrixType = right.matrixType;

        aSystem.matrixA = left * right.matrixA;
        aSystem.vectorB = left * right.vectorB;

        return aSystem;
    }

    template <typename Real>
    CSleSystem<Real> operator*(const Eigen::Matrix<Real, Eigen::Dynamic, 1>& left, const CSleSystem<Real>& right)
    {

        CSleSystem<Real> aSystem;

        aSystem.matrixType = right.matrixType;

        aSystem.matrixA = right.matrixA;
        for (int k = 0; k < aSystem.matrixA.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(aSystem.matrixA, k); it; ++it) {
                aSystem.matrixA.coeffRef(it.row(), k) *= left[it.row()];
            }
        }

        aSystem.vectorB = right.vectorB;
        for (int k = 0; k < aSystem.vectorB.size(); ++k)
            aSystem.vectorB[k] *= left[k];

        return aSystem;
    }
}
}
