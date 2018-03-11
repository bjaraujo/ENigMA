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

#include "GeoHashGrid.hpp"
#include "MshMesh.hpp"
#include "PdeField.hpp"

using namespace ENigMA::mesh;
using namespace ENigMA::pde;

namespace ENigMA
{

    namespace fem
    {

        template <typename Real, Integer Dof>
        class CFemCbsSolver
        {
        };
    
        // Two-dimensional    
        template <typename Real>
        class CFemCbsSolver<Real, 2>
        {
        protected:

            CMshMesh<Real> m_mesh;

            Real m_dt;

            Real m_dens;
            Real m_visc;
            
            Real m_gx, m_gy;

            CPdeField<Real> m_u, m_v;
            CPdeField<Real> m_p;
            
            Eigen::SparseMatrix<Real> m_G1, m_G2;
            
            virtual void calculateVelocityField();
            virtual void calculatePressureField();
            virtual void correctVelocityField();
            
        public:

            explicit CFemCbsSolver(CMshMesh<Real>& aMesh);
            ~CFemCbsSolver();

            virtual void setGravity(const Real gx, const Real gy);
            
            virtual void setMaterialProperties(const Real aDensity, const Real aViscosity);
            
            virtual void setTimeInterval(const Real dt);
            
            virtual void iterate(const Real dt, const bool bInit = false);

            CPdeField<Real>& u();
            CPdeField<Real>& v();

            CPdeField<Real>& p();

        };
    
        // Three-dimensional    
        template <typename Real>
        class CFemCbsSolver<Real, 3> : public CFemCbsSolver<Real, 2>
        {
        protected:

            Real m_gz;

            CPdeField<Real> m_w;
            
            Eigen::SparseMatrix<Real> m_G3;

            virtual void calculateVelocityField() override;
            virtual void calculatePressureField() override;
            virtual void correctVelocityField() override;
            
        public:

            explicit CFemCbsSolver(CMshMesh<Real>& aMesh);
            ~CFemCbsSolver();

            void setGravity(const Real gx, const Real gy) = delete;
            void setGravity(const Real gx, const Real gy, const Real gz);

            CPdeField<Real>& w();

        };

    }

}

#include "FemCBSSolver_Imp.hpp"

