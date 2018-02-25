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
        class CFemCBSSolver
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

            CFemCBSSolver();
            ~CFemCBSSolver();

            virtual void setTimeInterval(const Real dt);
            
            virtual void iterate(const Real dt, const bool bInit = false);

            Real u(const Integer aNodeId);
            Real v(const Integer aNodeId);
            Real w(const Integer aNodeId);
            Real p(const Integer aNodeId);

        };

    }

}

#include "FemCBSSolver_Imp.hpp"

