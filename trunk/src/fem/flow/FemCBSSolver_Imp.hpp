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

namespace ENigMA
{

    namespace fem
    {

        template <typename Real>
        CFemCbsSolver<Real, 2>::CFemCbsSolver(CMesh<Real>& aMesh)
        {

            m_G1 = gradient<Real>(u, CP_X).matrixA;
            m_G2 = gradient<Real>(v, CP_Y).matrixA;

            m_u.setMesh(aMesh);
            m_u.setDiscretMethod(DM_FEM);
            m_u.setDiscretOrder(DO_LINEAR);
            m_u.setDiscretLocation(DL_NODE);
            m_u.setSimulationType(ST_FLOW);
            m_u.setNbDofs(1);

            m_v.setMesh(aMesh);
            m_v.setDiscretMethod(DM_FEM);
            m_v.setDiscretOrder(DO_LINEAR);
            m_v.setDiscretLocation(DL_NODE);
            m_v.setSimulationType(ST_FLOW);
            m_v.setNbDofs(1);
            
            m_p.setMesh(aMesh);
            m_p.setDiscretMethod(DM_FEM);
            m_p.setDiscretOrder(DO_LINEAR);
            m_p.setDiscretLocation(DL_NODE);
            m_p.setSimulationType(ST_FLOW);
            m_p.setNbDofs(1);
            
        }

        template <typename Real>
        CFemCbsSolver<Real, 2>::~CFemCbsSolver()
        {

        }

        template <typename Real>
        void CFemCbsSolver<Real, 2>::setGravity(const Real gx, const Real gy)
        {

            m_gx = gx;
            m_gy = gy;
                    
        }
            
        template <typename Real>
        void CFemCbsSolver<Real, 2>::setMaterialProperties(const Real aDensity, const Real aViscosity)
        {

            m_dens = aDensity;
            m_visc = aViscosity;
            
        }
            
        template <typename Real>
        void CFemCbsSolver<Real, 2>::setTimeInterval(const Real dt)
        {

            m_dt = dt;

        }

        template <typename Real>
        void CFemCbsSolver<Real, 2>::calculateVelocityField()
        {

            Real nu = m_visc / m_dens;

            CPdeEquation<Real> aPdeEquationU(1.0 / m_dt * ddt<Real>(m_u) - nu * laplacian<Real>(m_u) + divergence<Real>(m_u, m_v, m_dt) = m_gx);
            aPdeEquationU.setElimination(m_u);
            aPdeEquationU.solve(m_u);

            CPdeEquation<Real> aPdeEquationV(1.0 / m_dt * ddt<Real>(m_v) - nu * laplacian<Real>(m_v) + divergence<Real>(m_u, m_v, m_dt) = m_gy);
            aPdeEquationV.setElimination(m_v);
            aPdeEquationV.solve(m_v);

        }

        template <typename Real>
        void CFemCbsSolver<Real, 2>::calculatePressureField()
        {

            // Pressure calculation
            CPdeEquation<Real> aPdeEquationP(laplacian<Real>(m_p) = m_dens / m_dt * (m_G1 * m_u.u + m_G2 * m_v.u));
            aPdeEquationP.setElimination(m_p);
            aPdeEquationP.solve(m_p);

        }

        template <typename Real>
        void CFemCbsSolver<Real, 2>::correctVelocityField()
        {

            CPdeEquation<Real> aPdeEquationUCor(ddt<Real>(m_u) = -m_dt / m_dens * m_G1 * m_p.u);
            CPdeEquation<Real> aPdeEquationVCor(ddt<Real>(m_v) = -m_dt / m_dens * m_G2 * m_p.u);

            aPdeEquationUCor.setElimination(m_u);
            aPdeEquationVCor.setElimination(m_v);

            aPdeEquationUCor.solve(m_u);
            aPdeEquationVCor.solve(m_v);

        }

        template <typename Real>
        void CFemCbsSolver<Real, 2>::iterate(const Real dt, const bool bInit)
        {

            this->setTimeInterval(dt);

            this->calculateVelocityField();
            this->calculatePressureField();

            this->correctVelocityField();
            
        }

        template <typename Real>
        Real CFemCbsSolver<Real, 2>::u(const Integer aNodeId)
        {

            return m_u.value(aControlVolumeId);

        }

        template <typename Real>
        Real CFemCbsSolver<Real, 2>::v(const Integer aNodeIndex)
        {

            return m_v.value(aNodeIndex);

        }

        template <typename Real>
        Real CFemCbsSolver<Real, 2>::p(const Integer aNodeIndex)
        {

            return m_p.value(aNodeIndex);

        }

        // Three-dimensional
        template <typename Real>
        CFemCbsSolver<Real, 3>::CFemCbsSolver(CMesh<Real>& aMesh) : CFemCbsSolver<Real, 2>(aMesh)
        {

            m_G3 = gradient<Real>(w, CP_Z).matrixA;

            m_w.setMesh(aMesh);
            m_w.setDiscretMethod(DM_FEM);
            m_w.setDiscretOrder(DO_LINEAR);
            m_w.setDiscretLocation(DL_NODE);
            m_w.setSimulationType(ST_FLOW);
            m_w.setNbDofs(1);
                        
        }

        template <typename Real>
        CFemCbsSolver<Real, 3>::~CFemCbsSolver()
        {

        }

        template <typename Real>
        void CFemCbsSolver<Real, 3>::setGravity(const Real gx, const Real gy, const Real gz)
        {

            CFemCbsSolver<Real, 2>::setGravity(gx, gy);
        
            m_gz = gz;
                    
        }
            
        template <typename Real>
        void CFemCbsSolver<Real, 3>::setMaterialProperties(const Real aDensity, const Real aViscosity)
        {

            m_dens = aDensity;
            m_visc = aViscosity;
            
        }

        template <typename Real>
        void CFemCbsSolver<Real, 3>::calculateVelocityField()
        {

            Real nu = m_visc / m_dens;

            CPdeEquation<Real> aPdeEquationU(1.0 / m_dt * ddt<Real>(m_u) - nu * laplacian<Real>(m_u) + divergence<Real>(m_u, m_v, m_w, m_dt) = m_gx);
            aPdeEquationU.setElimination(m_u);
            aPdeEquationU.solve(m_u);

            CPdeEquation<Real> aPdeEquationV(1.0 / m_dt * ddt<Real>(m_v) - nu * laplacian<Real>(m_v) + divergence<Real>(m_u, m_v, m_w, m_dt) = m_gy);
            aPdeEquationV.setElimination(m_v);
            aPdeEquationV.solve(m_v);

            CPdeEquation<Real> aPdeEquationW(1.0 / m_dt * ddt<Real>(m_w) - nu * laplacian<Real>(m_w) + divergence<Real>(m_u, m_v, m_w, m_dt) = m_gz);
            aPdeEquationW.setElimination(m_w);
            aPdeEquationW.solve(m_w);

        }

        template <typename Real>
        void CFemCbsSolver<Real, 3>::calculatePressureField()
        {

            // Pressure calculation
            CPdeEquation<Real> aPdeEquationP(laplacian<Real>(m_p) = m_dens / m_dt * (m_G1 * m_u.u + m_G2 * m_v.u + m_G3 * m_w.u));
            aPdeEquationP.setElimination(m_p);
            aPdeEquationP.solve(m_p);

        }

        template <typename Real>
        void CFemCbsSolver<Real, 3>::correctVelocityField()
        {

            CFemCbsSolver<Real, 2>::correctVelocityField();
        
            CPdeEquation<Real> aPdeEquationWCor(ddt<Real>(m_w) = -m_dt / m_dens * m_G3 * m_p.u);

            aPdeEquationWCor.setElimination(m_w);

            aPdeEquationWCor.solve(m_w);

        }

        template <typename Real>
        Real CFemCbsSolver<Real, 3>::w(const Integer aNodeIndex)
        {

            return m_w.value(aNodeIndex);

        }

    }

}


