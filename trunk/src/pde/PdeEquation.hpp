// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "PdeField.hpp"
#include "SleSystem.hpp"

namespace ENigMA {

namespace pde {

    enum EComponent {
        CP_X = 0,
        CP_Y = 1,
        CP_Z = 2
    };

    template <typename Real>
    class CPdeEquation {
    private:
        ENigMA::sle::CSleSystem<Real> m_system;
        std::vector<bool> m_bDeleteIndex;
        bool m_eliminationMethod;

    public:
        explicit CPdeEquation(const ENigMA::sle::CSleSystem<Real>& aSystem);
        virtual ~CPdeEquation();

        ENigMA::sle::CSleSystem<Real>& system();

        void setPenaltyFactor(ENigMA::pde::CPdeField<Real>& aField, Real aPenaltyFactor);
        void setElimination(ENigMA::pde::CPdeField<Real>& aField);

        void setSources(ENigMA::pde::CPdeField<Real>& aField);

        void solve(ENigMA::pde::CPdeField<Real>& aField);
    };

    // Mathematical operators
    template <typename Real>
    ENigMA::sle::CSleSystem<Real> operator*(const Real left, const ENigMA::pde::CPdeField<Real>& right);

    // Pde operators
    template <typename Real>
    ENigMA::sle::CSleSystem<Real> ddt(ENigMA::pde::CPdeField<Real>& aField);

    template <typename Real>
    ENigMA::sle::CSleSystem<Real> laplacian(ENigMA::pde::CPdeField<Real>& aField);

    template <typename Real>
    ENigMA::sle::CSleSystem<Real> divergence(ENigMA::pde::CPdeField<Real>& aField);

    template <typename Real>
    ENigMA::sle::CSleSystem<Real> divergence(ENigMA::pde::CPdeField<Real>& aField1, ENigMA::pde::CPdeField<Real>& aField2, Real dt);

    template <typename Real>
    ENigMA::sle::CSleSystem<Real> divergence(ENigMA::pde::CPdeField<Real>& aField1, ENigMA::pde::CPdeField<Real>& aField2, ENigMA::pde::CPdeField<Real>& aField3, Real dt);

    template <typename Real>
    ENigMA::sle::CSleSystem<Real> gradient(ENigMA::pde::CPdeField<Real>& aField, const EComponent aComponent);

    template <typename Real>
    Eigen::Matrix<Real, Eigen::Dynamic, 1> source(ENigMA::pde::CPdeField<Real>& aField, const Real aSource);
}
}

#include "PdeEquation_Imp.hpp"
