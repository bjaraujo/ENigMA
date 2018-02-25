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
#include "PdeField.hpp"

namespace ENigMA
{

    namespace fvm
    {

        template <typename Real>
        class CFemCbsSolver
        {
        protected:

            Real m_dt;

            Real m_dens;
            Real m_visc;
            
            CPdeField<Real> m_u, m_v, m_w, m_p;
            
            SparseMatrix<Real> m_G1, m_G2, m_G3;
            
            virtual void calculateVelocityField();
            virtual void calculatePressureField();
            virtual void correctVelocityField();
            
        public:

            CFemCbsSolver(CMesh<Real>& aMesh);
            ~CFemCbsSolver();

            void setGravity(const Real gx, const Real gy, const Real gz);
            
            virtual void setMaterialProperties(const Real aDensity, const Real aViscosity);
            
            virtual void setTimeInterval(const Real dt);
            
            virtual void iterate(const Real dt, const bool bInit = false);

            Real u(const Integer aNodeIndex);
            Real v(const Integer aNodeIndex);
            Real w(const Integer aNodeIndex);
            Real p(const Integer aNodeIndex);

        };

    }

}

#include "FemCBSSolver_Imp.hpp"

