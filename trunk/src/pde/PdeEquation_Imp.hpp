// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "bem/BemOperators.hpp"
#include "fdm/FdmOperators.hpp"
#include "fem/FemOperators.hpp"
#include "fvm/FvmOperators.hpp"

namespace ENigMA {

namespace pde {

    template <typename Real>
    CPdeEquation<Real>::CPdeEquation(const ENigMA::sle::CSleSystem<Real>& aSystem)
    {

        m_system = aSystem;

        m_eliminationMethod = false;
    }

    template <typename Real>
    CPdeEquation<Real>::~CPdeEquation()
    {
    }

    template <typename Real>
    ENigMA::sle::CSleSystem<Real>& CPdeEquation<Real>::system()
    {

        return m_system;
    }

    template <typename Real>
    void CPdeEquation<Real>::setPenaltyFactor(ENigMA::pde::CPdeField<Real>& aField, Real aPenaltyFactor)
    {

        if (aField.discretMethod() == DM_FEM || aField.discretMethod() == DM_FDM) {

            // Penalty
            for (typename std::map<Integer, Real>::const_iterator itr = aField.uFixed.begin(); itr != aField.uFixed.end(); ++itr) {

                Integer anIndex = itr->first;
                Real aValue = itr->second;

                m_system.matrixA.coeffRef(anIndex, anIndex) = aPenaltyFactor;

                m_system.vectorB(anIndex) = aPenaltyFactor * aValue;
            }
        }
    }

    template <typename Real>
    void CPdeEquation<Real>::setElimination(ENigMA::pde::CPdeField<Real>& aField)
    {

        if (aField.discretMethod() == DM_FEM) {

            m_bDeleteIndex.resize(m_system.vectorB.size(), false);

            // Elimination
            for (typename std::map<Integer, Real>::const_iterator itr = aField.uFixed.begin(); itr != aField.uFixed.end(); ++itr) {

                Integer anIndex = itr->first;
                Real aValue = itr->second;

                m_system.vectorB -= m_system.matrixA.col(anIndex) * aValue;

                //Mark row/column for deletion
                m_bDeleteIndex[anIndex] = true;
            }

            Integer newIndex = 0;
            std::vector<Integer> oldToNewIndex;

            for (int i = 0; i < m_system.vectorB.size(); ++i) {

                oldToNewIndex.push_back(newIndex);

                if (!m_bDeleteIndex[i])
                    newIndex++;
            }

            ENigMA::sle::CSleSystem<Real> aNewSystem;

            aNewSystem.matrixType = m_system.matrixType;

            aNewSystem.matrixA.resize(newIndex, newIndex);
            aNewSystem.matrixA.reserve(newIndex);
            aNewSystem.vectorB.resize(newIndex);

            for (int i = 0; i < m_system.vectorB.size(); ++i) {

                if (!m_bDeleteIndex[i]) {
                    aNewSystem.vectorB(oldToNewIndex[i]) = m_system.vectorB(i);
                }
            }

            for (int k = 0; k < m_system.matrixA.outerSize(); ++k) {
                for (typename Eigen::SparseMatrix<Real>::InnerIterator it(m_system.matrixA, k); it; ++it) {
                    if (!m_bDeleteIndex[it.row()] && !m_bDeleteIndex[it.col()]) {

                        int r = oldToNewIndex[it.row()];
                        int c = oldToNewIndex[it.col()];

                        aNewSystem.matrixA.coeffRef(r, c) = it.value();
                    }
                }
            }

            m_system = aNewSystem;

            m_eliminationMethod = true;
        }
    }

    template <typename Real>
    void CPdeEquation<Real>::setSources(ENigMA::pde::CPdeField<Real>& aField)
    {

        for (typename std::map<Integer, Real>::const_iterator itr = aField.uSource.begin(); itr != aField.uSource.end(); ++itr) {

            Integer anIndex = itr->first;
            Real aValue = itr->second;

            m_system.vectorB(anIndex) += aValue;
        }
    }

    template <typename Real>
    void CPdeEquation<Real>::solve(ENigMA::pde::CPdeField<Real>& aField)
    {

        if (m_eliminationMethod) {

            Integer newIndex = 0;
            std::map<Integer, Integer> newToOldIndex;

            for (Integer i = 0; i < static_cast<Integer>(m_bDeleteIndex.size()); ++i) {

                if (!m_bDeleteIndex[i]) {
                    newToOldIndex[newIndex] = i;
                    newIndex++;
                }
            }

            aField.u.resize(m_bDeleteIndex.size());

            Eigen::Matrix<Real, Eigen::Dynamic, 1> x = m_system.solve();

            for (int i = 0; i < x.size(); ++i) {

                aField.u(newToOldIndex[i]) = x(i);
            }

            for (typename std::map<Integer, Real>::const_iterator itr = aField.uFixed.begin(); itr != aField.uFixed.end(); ++itr) {

                Integer anIndex = itr->first;
                Real aValue = itr->second;

                aField.u(anIndex) = aValue;
            }

        } else {

            aField.u = m_system.solve();

            // Set fixed
            for (typename std::map<Integer, Real>::const_iterator itr = aField.uFixed.begin(); itr != aField.uFixed.end(); ++itr) {

                Integer anIndex = itr->first;
                Real aValue = itr->second;

                aField.u(anIndex) = aValue;
            }
        }
    }

    // Mathematical operators
    template <typename Real>
    ENigMA::sle::CSleSystem<Real> operator*(const Real left, const ENigMA::pde::CPdeField<Real>& right)
    {

        ENigMA::sle::CSleSystem<Real> aSystem;

        aSystem.matrixA.resize(right.u.size(), right.u.size());
        aSystem.matrixA.reserve(right.u.size());

        for (int k = 0; k < right.u.size(); ++k)
            aSystem.matrixA.coeffRef(k, k) = left;

        aSystem.vectorB.resize(right.u.size());

        aSystem.vectorB.setZero();

        return aSystem;
    }

    // Pde operators
    template <typename Real>
    ENigMA::sle::CSleSystem<Real> ddt(ENigMA::pde::CPdeField<Real>& aField)
    {

        ENigMA::sle::CSleSystem<Real> aSystem;

        if (aField.discretMethod() == DM_BEM)
            bem::ddt(aSystem, aField);
        else if (aField.discretMethod() == DM_FDM)
            fdm::ddt(aSystem, aField);
        else if (aField.discretMethod() == DM_FEM)
            fem::ddt(aSystem, aField);
        else if (aField.discretMethod() == DM_FVM)
            fvm::ddt(aSystem, aField);

        return aSystem;
    }

    template <typename Real>
    ENigMA::sle::CSleSystem<Real> laplacian(ENigMA::pde::CPdeField<Real>& aField)
    {

        ENigMA::sle::CSleSystem<Real> aSystem;

        if (aField.discretMethod() == DM_BEM)
            bem::laplacian(aSystem, aField);
        else if (aField.discretMethod() == DM_FDM)
            fdm::laplacian(aSystem, aField);
        else if (aField.discretMethod() == DM_FEM)
            fem::laplacian(aSystem, aField);
        else if (aField.discretMethod() == DM_FVM)
            fvm::laplacian(aSystem, aField);

        return aSystem;
    }

    template <typename Real>
    ENigMA::sle::CSleSystem<Real> divergence(ENigMA::pde::CPdeField<Real>& aField)
    {

        ENigMA::sle::CSleSystem<Real> aSystem;

        if (aField.discretMethod() == DM_BEM)
            bem::divergence(aSystem, aField);
        else if (aField.discretMethod() == DM_FDM)
            fdm::divergence(aSystem, aField);
        else if (aField.discretMethod() == DM_FEM)
            fem::divergence(aSystem, aField);
        else if (aField.discretMethod() == DM_FVM)
            fvm::divergence(aSystem, aField);

        return aSystem;
    }

    template <typename Real>
    ENigMA::sle::CSleSystem<Real> divergence(ENigMA::pde::CPdeField<Real>& aField1, ENigMA::pde::CPdeField<Real>& aField2, Real dt)
    {

        ENigMA::sle::CSleSystem<Real> aSystem;

        if (aField1.discretMethod() == DM_FEM)
            fem::divergence(aSystem, aField1, aField2, dt);

        return aSystem;
    }

    template <typename Real>
    ENigMA::sle::CSleSystem<Real> divergence(ENigMA::pde::CPdeField<Real>& aField1, ENigMA::pde::CPdeField<Real>& aField2, ENigMA::pde::CPdeField<Real>& aField3, Real dt)
    {

        ENigMA::sle::CSleSystem<Real> aSystem;

        if (aField1.discretMethod() == DM_FEM)
            fem::divergence(aSystem, aField1, aField2, aField3, dt);

        return aSystem;
    }

    template <typename Real>
    ENigMA::sle::CSleSystem<Real> gradient(ENigMA::pde::CPdeField<Real>& aField, const EComponent aComponent)
    {

        ENigMA::sle::CSleSystem<Real> aSystem;

        if (aField.discretMethod() == DM_FEM)
            fem::gradient(aSystem, aField, aComponent);

        return aSystem;
    }

    template <typename Real>
    Eigen::Matrix<Real, Eigen::Dynamic, 1> source(ENigMA::pde::CPdeField<Real>& aField, const Real aSource)
    {

        Eigen::Matrix<Real, Eigen::Dynamic, 1> aVectorB;

        if (aField.discretMethod() == DM_BEM)
            bem::source(aVectorB, aField, aSource);
        else if (aField.discretMethod() == DM_FDM)
            fdm::source(aVectorB, aField, aSource);
        else if (aField.discretMethod() == DM_FEM)
            fem::source(aVectorB, aField, aSource);
        else if (aField.discretMethod() == DM_FVM)
            fvm::source(aVectorB, aField, aSource);

        return aVectorB;
    }
}
}
