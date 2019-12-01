// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "FvmTemperatureSolver.hpp"

namespace ENigMA {
namespace fvm {
    template <typename Real>
    class CFvmVofSolver : public CFvmPisoSolver<Real> {
    protected:
	typedef std::map<Integer, Real> varMap;

        Real m_dens0;
        Real m_visc0;

        Real m_dens1;
        Real m_visc1;

        varMap m_s0; // Cell center gamma - previous time step
        varMap m_s; // Cell center gamma
        varMap m_sf; // Face center gamma

        varMap m_betaf; // Face center interpolation coefficient
        varMap m_Co; // Cell center Courant

        virtual void storePreviousQuantities();

        virtual void updateProperties();

        virtual void calculateGammaField();

        virtual Real calculateCourant(double dt, bool bInterface);
        virtual void predictBeta(const Real aTolerance = 0.0);
        virtual void correctBeta(double dt, const Real aTolerance = 0.0);

    public:
        explicit CFvmVofSolver(CFvmMesh<Real>& aFvmMesh);
        virtual ~CFvmVofSolver();

        virtual void setMaterialProperties(const Real aDensity0, const Real aViscosity0, const Real aDensity1, const Real aViscosity1);

        void setInitialGamma(const Integer aControlVolumeId, const Real s);

        void setBoundaryGamma(const Integer aFaceId, const EBoundaryType sFaceType, const Real s);
        void setBoundaryGamma(const std::vector<Integer>& sFaceIds, const EBoundaryType sFaceType, const Real s);

        virtual void iterate(const Real dt, const bool bInit = false);
        virtual void residual(Real& ru, Real& rv, Real& rw, Real& rp, Real& rs);

        Real s(const Integer aControlVolumeId);

        Real sf(const Integer aFaceId);

        CGeoVector<Real> grads(const Integer aControlVolumeId);
    };
}
}

#include "FvmVofSolver_Imp.hpp"
