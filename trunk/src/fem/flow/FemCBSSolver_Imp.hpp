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

    namespace fvm
    {

        template <typename Real>
        CFemCBSSolver<Real>::CFemCBSSolver()
        {

            m_G1 = gradient<Real>(u, CP_X).matrixA;
            m_G2 = gradient<Real>(v, CP_Y).matrixA;
            m_G3 = gradient<Real>(w, CP_Z).matrixA;

        }

        template <typename Real>
        CFemCBSSolver<Real>::~CFemCBSSolver()
        {

        }

        template <typename Real>
        void CFemCBSSolver<Real>::setTimeInterval(const Real dt)
        {

            m_dt = dt;

        }

        template <typename Real>
        void CFemCBSSolver<Real>::calculateVelocityField()
        {

            Real nu = m_visc / m_dens;

            CPdeEquation<Real> aPdeEquationU(1.0 / m_dt * ddt<Real>(m_u) - nu * laplacian<Real>(m_u) + divergence<Real>(m_u, m_v, m_w, m_dt) = 0);
            aPdeEquationU.setElimination(m_u);
            aPdeEquationU.solve(m_u);

            CPdeEquation<Real> aPdeEquationV(1.0 / m_dt * ddt<Real>(m_v) - nu * laplacian<Real>(m_v) + divergence<Real>(m_u, m_v, m_w, m_dt) = 0);
            aPdeEquationV.setElimination(m_v);
            aPdeEquationV.solve(m_v);

            CPdeEquation<Real> aPdeEquationW(1.0 / m_dt * ddt<Real>(m_w) - nu * laplacian<Real>(m_w) + divergence<Real>(m_u, m_v, m_w, m_dt) = 0);
            aPdeEquationW.setElimination(m_w);
            aPdeEquationW.solve(m_w);

        }

        template <typename Real>
        void CFemCBSSolver<Real>::calculatePressureField()
        {

            // Pressure calculation
            CPdeEquation<Real> aPdeEquationP(laplacian<Real>(m_p) = m_dens / m_dt * (m_G1 * m_u.u + m_G2 * m_v.u + m_G3 * m_w.u));
            aPdeEquationP.setElimination(m_p);
            aPdeEquationP.solve(m_p);

        }

        template <typename Real>
        void CFemCBSSolver<Real>::correctVelocityField()
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
        void CFemCBSSolver<Real>::iterate(const Real dt, const bool bInit)
        {

            this->setTimeInterval(dt);

            this->calculateVelocityField();
            this->calculatePressureField();

            this->correctVelocityField();
            
        }

        template <typename Real>
        Real CFemCBSSolver<Real>::u(const Integer aControlVolumeId)
        {

            return m_u[aControlVolumeId];

        }

        template <typename Real>
        Real CFemCBSSolver<Real>::v(const Integer aControlVolumeId)
        {

            return m_v[aControlVolumeId];

        }

        template <typename Real>
        Real CFemCBSSolver<Real>::w(const Integer aControlVolumeId)
        {

            return m_w[aControlVolumeId];

        }

        template <typename Real>
        Real CFemCBSSolver<Real>::p(const Integer aControlVolumeId)
        {

            return m_p[aControlVolumeId];

        }

    }

}


