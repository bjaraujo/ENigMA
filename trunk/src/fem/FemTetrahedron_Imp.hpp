// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// <Copyright> Copyright (c) 2012, All Rights Reserved </Copyright>
// *****************************************************************************

#pragma once

using namespace ENigMA::geometry;

namespace ENigMA
{

    namespace fem
    {

        template <typename Real>
        CFemTetrahedron<Real, 4, 1, 1>::CFemTetrahedron()
        {
            
            this->source.resize(4);
            this->source.setZero();

        }

        template <typename Real>
        CFemTetrahedron<Real, 4, 1, 1>::~CFemTetrahedron()
        {

        }

        template <typename Real>
        void CFemTetrahedron<Real, 4, 1, 1>::rebuild()
        {

            CGeoCoordinate<Real> n1 = this->m_vertices[0];
            CGeoCoordinate<Real> n2 = this->m_vertices[1];
            CGeoCoordinate<Real> n3 = this->m_vertices[2];
            CGeoCoordinate<Real> n4 = this->m_vertices[3];

            m_x1 = n1.x();
            m_x2 = n2.x();
            m_x3 = n3.x();
            m_x4 = n4.x();

            m_y1 = n1.y();
            m_y2 = n2.y();
            m_y3 = n3.y();
            m_y4 = n4.y();

            m_z1 = n1.z();
            m_z2 = n2.z();
            m_z3 = n3.z();
            m_z4 = n4.z();

            m_x14 = m_x1 - m_x4;
            m_y14 = m_y1 - m_y4;
            m_z14 = m_z1 - m_z4;

            m_x24 = m_x2 - m_x4;
            m_y24 = m_y2 - m_y4;
            m_z24 = m_z2 - m_z4;

            m_x34 = m_x3 - m_x4;
            m_y34 = m_y3 - m_y4;
            m_z34 = m_z3 - m_z4;

        }

        template <typename Real>
        void CFemTetrahedron<Real, 4, 1, 1>::setTransientTerm()
        {

            this->ddt.resize(4, 4);

            this->ddt << 0.10, 0.05, 0.05, 0.05,
                                      0.05, 0.10, 0.05, 0.05,
                                      0.05, 0.05, 0.10, 0.05,
                                      0.05, 0.05, 0.05, 0.10;

            this->ddt *= this->m_volume;

            if (this->m_dt > 0.0) 
                this->ddt /= this->m_dt;

        }

        template <typename Real>
        void CFemTetrahedron<Real, 4, 1, 1>::setDiffusionTerm()
        {

            this->laplacian.resize(4, 4);

            Real bi = m_y24 * m_z34 - m_y34 * m_z24;
            Real bj = m_y34 * m_z14 - m_y14 * m_z34;
            Real bk = m_y14 * m_z24 - m_y24 * m_z14;
            Real bl = -(bi + bj + bk);

            Real ci = m_z24 * m_x34 - m_z34 * m_x24;
            Real cj = m_z34 * m_x14 - m_z14 * m_x34;
            Real ck = m_z14 * m_x24 - m_z24 * m_x14;
            Real cl = -(ci + cj + ck);

            Real di = m_x24 * m_y34 - m_x34 * m_y24;
            Real dj = m_x34 * m_y14 - m_x14 * m_y34;
            Real dk = m_x14 * m_y24 - m_x24 * m_y14;
            Real dl = -(di + dj + dk);

            Eigen::Matrix<Real, 3, 4> B;

            B << bi, bj, bk, bl,
                 ci, cj, ck, cl,
                 di, dj, dk, dl;
            
            if (this->m_volume > 0.0)
                B *= 1.0/(6.0 * this->m_volume);

            this->laplacian = -B.transpose() * B;

            this->laplacian *= this->m_volume;

        }

        template <typename Real>
        void CFemTetrahedron<Real, 4, 1, 1>::setConvectiveTerm()
        {

        }

        template <typename Real>
        void CFemTetrahedron<Real, 4, 1, 1>::update()
        {

            // Calculate geometrical properties
            rebuild();

            // Diffusion term
            setDiffusionTerm();

            // Transient term
            if (this->m_transient)
                setTransientTerm();

        }

        template <typename Real>
        void CFemTetrahedron<Real, 4, 1, 1>::setSourceOnNode(const Integer aNodeIndex, const Real aValue)
        {

            this->source(aNodeIndex) += aValue;

        }

        template <typename Real>
        void CFemTetrahedron<Real, 4, 1, 1>::setSourceOnEdge(const Integer anEdgeIndex, const Real aValue)
        {

            this->source((anEdgeIndex + 0) % 4) += aValue * 0.5;
            this->source((anEdgeIndex + 1) % 4) += aValue * 0.5;

        }

        template <typename Real>
        void CFemTetrahedron<Real, 4, 1, 1>::setSourceOnFace(const Integer aFaceIndex, const Real aValue)
        {

        }

        template <typename Real>
        void CFemTetrahedron<Real, 4, 1, 1>::calculateB(Eigen::Matrix<Real, 3, 4>& B)
        {

            Real bi = m_y24 * m_z34 - m_y34 * m_z24;
            Real bj = m_y34 * m_z14 - m_y14 * m_z34;
            Real bk = m_y14 * m_z24 - m_y24 * m_z14;
            Real bl = -(bi + bj + bk);

            Real ci = m_z24 * m_x34 - m_z34 * m_x24;
            Real cj = m_z34 * m_x14 - m_z14 * m_x34;
            Real ck = m_z14 * m_x24 - m_z24 * m_x14;
            Real cl = -(ci + cj + ck);

            Real di = m_x24 * m_y34 - m_x34 * m_y24;
            Real dj = m_x34 * m_y14 - m_x14 * m_y34;
            Real dk = m_x14 * m_y24 - m_x24 * m_y14;
            Real dl = -(di + dj + dk);

            B << bi, bj, bk, bl,
                 ci, cj, ck, cl,
                 di, dj, dk, dl;

        }

        template <typename Real>
        void CFemTetrahedron<Real, 4, 3, 1>::calculateB(Eigen::Matrix<Real, 6, 12>& B)
        {

            Real x14 = this->m_x14;
            Real x24 = this->m_x24;
            Real x34 = this->m_x34;

            Real y14 = this->m_y14;
            Real y24 = this->m_y24;
            Real y34 = this->m_y34;

            Real z14 = this->m_z14;
            Real z24 = this->m_z24;
            Real z34 = this->m_z34;

            Eigen::Matrix<Real, 3, 3> A;

            A << y24*z34-y34*z24, y34*z14-y14*z34, y14*z24-y24*z14,
                 z24*x34-z34*x24, z34*x14-z14*x34, z14*x24-z24*x14,
                 x24*y34-x34*y24, x34*y14-x14*y34, x14*y24-x24*y14;

            A /= (6.0 * this->m_volume);

            Real A11 = A(0, 0);
            Real A12 = A(0, 1);
            Real A13 = A(0, 2);

            Real A21 = A(1, 0);
            Real A22 = A(1, 1);
            Real A23 = A(1, 2);

            Real A31 = A(2, 0);
            Real A32 = A(2, 1);
            Real A33 = A(2, 2);

            Real A1t = A11 + A12 + A13;
            Real A2t = A21 + A22 + A23;
            Real A3t = A31 + A32 + A33;

            B << A11, 0.0, 0.0, A12, 0.0, 0.0, A13, 0.0, 0.0, -A1t, 0.0, 0.0,   
                 0.0, A21, 0.0, 0.0, A22, 0.0, 0.0, A23, 0.0, 0.0, -A2t, 0.0,
                 0.0, 0.0, A31, 0.0, 0.0, A32, 0.0, 0.0, A33, 0.0, 0.0, -A3t,
                 0.0, A31, A21, 0.0, A32, A22, 0.0, A33, A23, 0.0, -A3t, -A2t,
                 A31, 0.0, A11, A32, 0.0, A12, A33, 0.0, A13, -A3t, 0.0, -A1t,
                 A21, A11, 0.0, A22, A12, 0.0, A23, A13, 0.0, -A2t, -A1t, 0.0;

        }

    }

}
