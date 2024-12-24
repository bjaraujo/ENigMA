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

#include "PdeEquation.hpp"

namespace ENigMA
{
    namespace fem
    {
        template <typename Real>
        CFemCbsSolver<Real, 2>::CFemCbsSolver(CMshMesh<Real>& aMesh)
        {
            m_mesh = aMesh;

            m_u.setMesh(aMesh);
            m_u.setDiscretMethod(DM_FEM);
            m_u.setDiscretOrder(DO_LINEAR);
            m_u.setDiscretLocation(DL_NODE);
            m_u.setSimulationType(ST_FLOW);
            m_u.setNbDofs(1);
            m_u.u.resize(aMesh.nbNodes());

            m_v.setMesh(aMesh);
            m_v.setDiscretMethod(DM_FEM);
            m_v.setDiscretOrder(DO_LINEAR);
            m_v.setDiscretLocation(DL_NODE);
            m_v.setSimulationType(ST_FLOW);
            m_v.setNbDofs(1);
            m_v.u.resize(aMesh.nbNodes());

            m_p.setMesh(aMesh);
            m_p.setDiscretMethod(DM_FEM);
            m_p.setDiscretOrder(DO_LINEAR);
            m_p.setDiscretLocation(DL_NODE);
            m_p.setSimulationType(ST_FLOW);
            m_p.setNbDofs(1);
            m_p.u.resize(aMesh.nbNodes());

            m_G1 = gradient<Real>(m_u, CP_X).matrixA;
            m_G2 = gradient<Real>(m_v, CP_Y).matrixA;
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
            CPdeEquation<Real> aPdeEquationP(laplacian<Real>(m_p) = 1.0 / m_dt * (m_G1 * m_u.u + m_G2 * m_v.u));
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
        CPdeField<Real>& CFemCbsSolver<Real, 2>::u()
        {
            return m_u;
        }

        template <typename Real>
        CPdeField<Real>& CFemCbsSolver<Real, 2>::v()
        {
            return m_v;
        }

        template <typename Real>
        CPdeField<Real>& CFemCbsSolver<Real, 2>::p()
        {
            return m_p;
        }

        // Three-dimensional
        template <typename Real>
        CFemCbsSolver<Real, 3>::CFemCbsSolver(CMshMesh<Real>& aMesh)
        {
            m_mesh = aMesh;

            m_u.setMesh(aMesh);
            m_u.setDiscretMethod(DM_FEM);
            m_u.setDiscretOrder(DO_LINEAR);
            m_u.setDiscretLocation(DL_NODE);
            m_u.setSimulationType(ST_FLOW);
            m_u.setNbDofs(1);
            m_u.u.resize(aMesh.nbNodes());

            m_v.setMesh(aMesh);
            m_v.setDiscretMethod(DM_FEM);
            m_v.setDiscretOrder(DO_LINEAR);
            m_v.setDiscretLocation(DL_NODE);
            m_v.setSimulationType(ST_FLOW);
            m_v.setNbDofs(1);
            m_v.u.resize(aMesh.nbNodes());

            m_w.setMesh(aMesh);
            m_w.setDiscretMethod(DM_FEM);
            m_w.setDiscretOrder(DO_LINEAR);
            m_w.setDiscretLocation(DL_NODE);
            m_w.setSimulationType(ST_FLOW);
            m_w.setNbDofs(1);
            m_w.u.resize(aMesh.nbNodes());

            m_p.setMesh(aMesh);
            m_p.setDiscretMethod(DM_FEM);
            m_p.setDiscretOrder(DO_LINEAR);
            m_p.setDiscretLocation(DL_NODE);
            m_p.setSimulationType(ST_FLOW);
            m_p.setNbDofs(1);
            m_p.u.resize(aMesh.nbNodes());

            m_G1 = gradient<Real>(m_u, CP_X).matrixA;
            m_G2 = gradient<Real>(m_v, CP_Y).matrixA;
            m_G3 = gradient<Real>(m_w, CP_Z).matrixA;
        }

        template <typename Real>
        CFemCbsSolver<Real, 3>::~CFemCbsSolver()
        {
        }

        template <typename Real>
        void CFemCbsSolver<Real, 3>::setGravity(const Real gx, const Real gy, const Real gz)
        {
            m_gx = gx;
            m_gy = gy;
            m_gz = gz;
        }

        template <typename Real>
        void CFemCbsSolver<Real, 3>::setMaterialProperties(const Real aDensity, const Real aViscosity)
        {
            m_dens = aDensity;
            m_visc = aViscosity;
        }

        template <typename Real>
        void CFemCbsSolver<Real, 3>::setTimeInterval(const Real dt)
        {
            m_dt = dt;
        }

        template <typename Real>
        void CFemCbsSolver<Real, 3>::calculateVelocityField()
        {
            Real nu = m_visc / m_dens;

            CPdeEquation<Real> aPdeEquationU(1.0 / m_dt * ddt<Real>(m_u) - nu * laplacian<Real>(m_u) + divergence<Real>(m_u, m_v, m_w, m_dt) = m_gx / m_dens);
            aPdeEquationU.setElimination(m_u);
            aPdeEquationU.solve(m_u);

            CPdeEquation<Real> aPdeEquationV(1.0 / m_dt * ddt<Real>(m_v) - nu * laplacian<Real>(m_v) + divergence<Real>(m_u, m_v, m_w, m_dt) = m_gy / m_dens);
            aPdeEquationV.setElimination(m_v);
            aPdeEquationV.solve(m_v);

            CPdeEquation<Real> aPdeEquationW(1.0 / m_dt * ddt<Real>(m_w) - nu * laplacian<Real>(m_w) + divergence<Real>(m_u, m_v, m_w, m_dt) = m_gz / m_dens);
            aPdeEquationW.setElimination(m_w);
            aPdeEquationW.solve(m_w);
        }

        template <typename Real>
        void CFemCbsSolver<Real, 3>::calculatePressureField()
        {
            // Pressure calculation
            CPdeEquation<Real> aPdeEquationP(laplacian<Real>(m_p) = 1.0 / m_dt * (m_G1 * m_u.u + m_G2 * m_v.u + m_G3 * m_w.u));
            aPdeEquationP.setElimination(m_p);
            aPdeEquationP.solve(m_p);
        }

        template <typename Real>
        void CFemCbsSolver<Real, 3>::correctVelocityField()
        {
            CPdeEquation<Real> aPdeEquationUCor(ddt<Real>(m_u) = -m_dt / m_dens * m_G1 * m_p.u);
            CPdeEquation<Real> aPdeEquationVCor(ddt<Real>(m_v) = -m_dt / m_dens * m_G2 * m_p.u);
            CPdeEquation<Real> aPdeEquationWCor(ddt<Real>(m_w) = -m_dt / m_dens * m_G3 * m_p.u);

            aPdeEquationUCor.setElimination(m_u);
            aPdeEquationVCor.setElimination(m_v);
            aPdeEquationWCor.setElimination(m_w);

            aPdeEquationUCor.solve(m_u);
            aPdeEquationVCor.solve(m_v);
            aPdeEquationWCor.solve(m_w);
        }

        template <typename Real>
        void CFemCbsSolver<Real, 3>::iterate(const Real dt, const bool bInit)
        {
            this->setTimeInterval(dt);

            this->calculateVelocityField();
            this->calculatePressureField();

            this->correctVelocityField();
        }

        template <typename Real>
        CPdeField<Real>& CFemCbsSolver<Real, 3>::u()
        {
            return m_u;
        }

        template <typename Real>
        CPdeField<Real>& CFemCbsSolver<Real, 3>::v()
        {
            return m_v;
        }

        template <typename Real>
        CPdeField<Real>& CFemCbsSolver<Real, 3>::w()
        {
            return m_w;
        }

        template <typename Real>
        CPdeField<Real>& CFemCbsSolver<Real, 3>::p()
        {
            return m_p;
        }
    }
}
