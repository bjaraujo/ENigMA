// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "FvmPisoSolver.hpp"

namespace ENigMA {
namespace fvm {
    template <typename Real>
    class CFvmTemperatureSolver : public CFvmPisoSolver<Real> {
    protected:
	typedef std::map<Integer, Real> varMap;

        bool m_calcT;

        Real m_bthcond; // Boundary thermal conductivity

        varMap m_thcond; // Cell thermal conductivity
        varMap m_spheat; // Cell specific heat

        varMap m_T0; // Cell center temperature - previous time step
        varMap m_T; // Cell center temperature
        varMap m_Tf; // Face center temperature

        void calculateTemperatureField();

    public:
        explicit CFvmTemperatureSolver(CFvmMesh<Real>& aFvmMesh);
        virtual ~CFvmTemperatureSolver();

        virtual void storePreviousQuantities();

        virtual void setMaterialProperties(const Real aDensity, const Real aViscosity, const Real aThermalConductivity, const Real aSpecificHeat);

        void setBoundaryTemperature(const std::vector<Integer>& sFaceIds, const EBoundaryType sFaceType, const Real T);

        virtual void iterate(const Real dt, const bool bInit = false);
        virtual void residual(Real& ru, Real& rv, Real& rw, Real& rp, Real& rT);

        Real T(const Integer aControlVolumeId);

        Real Tf(const Integer aFaceId);

        CGeoVector<Real> gradT(const Integer aControlVolumeId);
    };
}
}

#include "FvmTemperatureSolver_Imp.hpp"
