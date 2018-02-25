// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include <map>
#include <vector>

#include "GeoHashGrid.hpp"

#include "FvmMesh.hpp"

namespace ENigMA
{

    namespace fvm
    {

        template <typename Real>
        class CFvmPisoSolver
        {
        protected:

            typedef std::map<Integer, Real> varMap;
            typedef std::map<Integer, CGeoVector<Real> > vecMap;

            std::map<Integer, Integer> m_mapIdToIndex;
            std::map<Integer, Integer> m_mapIndexToId;

            bool m_calcu, m_calcv, m_calcw;
            bool m_calcp;

            Real m_dt;

            Real m_gx, m_gy, m_gz;      // Gravity

            varMap m_dens;              // Cell density
            varMap m_visc;              // Cell viscosity

            varMap m_flux;              // Face flux

            varMap m_u0, m_v0, m_w0;    // Cell center velocity - previous time step
            varMap m_u, m_v, m_w;       // Cell center velocity
            varMap m_uf, m_vf, m_wf;    // Face center velocity

            varMap m_p0;                // Cell center pressure - previous time step
            varMap m_p;                 // Cell center pressure
            varMap m_pf;                // Face center pressure

            varMap m_ap;
            varMap m_bu, m_bv, m_bw;
            varMap m_Hu, m_Hv, m_Hw;

            varMap m_massError;         // Mas error

            CFvmMesh<Real> m_fvmMesh;

            CGeoVector<Real> gradient(varMap& var, varMap& varf, const Integer aControlVolumeId);

            virtual void setTimeInterval(const Real dt);
            
            virtual void storePreviousQuantities();

            virtual void calculateVelocityField();
            virtual void calculatePressureField();

            virtual void correctFlux();
            virtual void correctVelocityField();
            virtual void correctPressureField();

            Real checkMassConservationCell(const Integer aControlVolumeId);

        public:

            CFvmPisoSolver(CFvmMesh<Real>& aFvmMesh);
            ~CFvmPisoSolver();

            void setGravity(const Real gx, const Real gy, const Real gz);

            virtual void setMaterialProperties(const Real aDensity, const Real aViscosity);

            void setBoundaryVelocity(const std::vector<Integer>& sFaceIds, EBoundaryType sFaceType, const Real u, const Real v, const Real w);
            void setBoundaryPressure(const std::vector<Integer>& sFaceIds, EBoundaryType sFaceType, const Real p);

            virtual void iterate(const Real dt, const bool bInit = false);
            virtual void checkMassConservation(Real& aMassError);
            virtual void residual(Real& ru, Real& rv, Real& rw, Real& rp);

            Real u(const Integer aControlVolumeId);
            Real v(const Integer aControlVolumeId);
            Real w(const Integer aControlVolumeId);
            Real p(const Integer aControlVolumeId);

            Real uf(const Integer aFaceId);
            Real vf(const Integer aFaceId);
            Real wf(const Integer aFaceId);
            Real pf(const Integer aFaceId);

            Real flux(const Integer aFaceId);

            CGeoVector<Real> gradu(const Integer aControlVolumeId);
            CGeoVector<Real> gradv(const Integer aControlVolumeId);
            CGeoVector<Real> gradw(const Integer aControlVolumeId);
            CGeoVector<Real> gradp(const Integer aControlVolumeId);

            Real massError(const Integer aControlVolumeId);

        };

    }

}

#include "FvmPisoSolver_Imp.hpp"

