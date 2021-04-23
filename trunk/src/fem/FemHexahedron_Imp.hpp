// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

using namespace ENigMA::geometry;

namespace ENigMA
{
    namespace fem
    {
        template <typename Real>
        CFemHexahedron<Real, 8, 1, 1>::CFemHexahedron()
        {
            this->source.resize(8);
            this->source.setZero();

            CIntGaussIntegration<Real>::m_integPoints = 8;
        }

        template <typename Real>
        CFemHexahedron<Real, 8, 1, 1>::~CFemHexahedron()
        {
        }

        template <typename Real>
        void CFemHexahedron<Real, 8, 1, 1>::rebuild()
        {
            CGeoCoordinate<Real> n1 = this->m_vertices[0];
            CGeoCoordinate<Real> n2 = this->m_vertices[1];
            CGeoCoordinate<Real> n3 = this->m_vertices[2];
            CGeoCoordinate<Real> n4 = this->m_vertices[3];
            CGeoCoordinate<Real> n5 = this->m_vertices[4];
            CGeoCoordinate<Real> n6 = this->m_vertices[5];
            CGeoCoordinate<Real> n7 = this->m_vertices[6];
            CGeoCoordinate<Real> n8 = this->m_vertices[7];

            m_x1 = n1.x();
            m_x2 = n2.x();
            m_x3 = n3.x();
            m_x4 = n4.x();
            m_x5 = n5.x();
            m_x6 = n6.x();
            m_x7 = n7.x();
            m_x8 = n8.x();

            m_y1 = n1.y();
            m_y2 = n2.y();
            m_y3 = n3.y();
            m_y4 = n4.y();
            m_y5 = n5.y();
            m_y6 = n6.y();
            m_y7 = n7.y();
            m_y8 = n8.y();

            m_z1 = n1.z();
            m_z2 = n2.z();
            m_z3 = n3.z();
            m_z4 = n4.z();
            m_z5 = n5.z();
            m_z6 = n6.z();
            m_z7 = n7.z();
            m_z8 = n8.z();
        }

        template <typename Real>
        void CFemHexahedron<Real, 8, 1, 1>::setTransientTerm()
        {
            Real N[8];

            Eigen::Matrix<Real, 3, 3> J;
            Eigen::Matrix<Real, 8, 8> Ji;

            this->ddt.resize(8, 8);
            this->ddt.setZero();

            for (Integer p = 0; p < CIntGaussIntegration<Real>::m_integPoints; p++)
            {
                Real xi = this->m_xi[p];
                Real eta = this->m_eta[p];
                Real zeta = this->m_zeta[p];

                Real wxi = this->m_wxi[p];
                Real weta = this->m_weta[p];
                Real wzeta = this->m_wzeta[p];

                Real xim = (1 - xi);
                Real xip = (1 + xi);

                Real etam = (1 - eta);
                Real etap = (1 + eta);

                Real zetam = (1 - zeta);
                Real zetap = (1 + zeta);

                J << -etam * zetam * m_x1 + etam * zetam * m_x2 - etap * zetam * m_x3 + etap * zetam * m_x4 - etam * zetap * m_x5 + etam * zetap * m_x6 - etap * zetap * m_x7 + etap * zetap * m_x8,
                    -etam * zetam * m_y1 + etam * zetam * m_y2 - etap * zetam * m_y3 + etap * zetam * m_y4 - etam * zetap * m_y5 + etam * zetap * m_y6 - etap * zetap * m_y7 + etap * zetap * m_y8,
                    -etam * zetam * m_z1 + etam * zetam * m_z2 - etap * zetam * m_z3 + etap * zetam * m_z4 - etam * zetap * m_z5 + etam * zetap * m_z6 - etap * zetap * m_z7 + etap * zetap * m_z8,
                    -xim * zetam * m_x1 - xip * zetam * m_x2 + xim * zetam * m_x3 + xip * zetam * m_x4 - xim * zetap * m_x5 - xip * zetap * m_x6 + xim * zetap * m_x7 + xip * zetap * m_x8,
                    -xim * zetam * m_y1 - xip * zetam * m_y2 + xim * zetam * m_y3 + xip * zetam * m_y4 - xim * zetap * m_y5 - xip * zetap * m_y6 + xim * zetap * m_y7 + xip * zetap * m_y8,
                    -xim * zetam * m_z1 - xip * zetam * m_z2 + xim * zetam * m_z3 + xip * zetam * m_z4 - xim * zetap * m_z5 - xip * zetap * m_z6 + xim * zetap * m_z7 + xip * zetap * m_z8,
                    -etam * xim * m_x1 - etam * xip * m_x2 - etap * xim * m_x3 - etap * xip * m_x4 + etam * xim * m_x5 + etam * xip * m_x6 + etap * xim * m_x7 + etap * xip * m_x8,
                    -etam * xim * m_y1 - etam * xip * m_y2 - etap * xim * m_y3 - etap * xip * m_y4 + etam * xim * m_y5 + etam * xip * m_y6 + etap * xim * m_y7 + etap * xip * m_y8,
                    -etam * xim * m_z1 - etam * xip * m_z2 - etap * xim * m_z3 - etap * xip * m_z4 + etam * xim * m_z5 + etam * xip * m_z6 + etap * xim * m_z7 + etap * xip * m_z8;

                J *= 0.125;

                Real detJ = J.determinant();

                N[0] = 0.125 * xim * etam * zetam;
                N[1] = 0.125 * xip * etam * zetam;
                N[2] = 0.125 * xim * etap * zetam;
                N[3] = 0.125 * xip * etap * zetam;
                N[4] = 0.125 * xim * etam * zetap;
                N[5] = 0.125 * xip * etam * zetap;
                N[6] = 0.125 * xim * etap * zetap;
                N[7] = 0.125 * xip * etap * zetap;

                Ji.setZero();

                for (Integer i = 0; i < 8; ++i)
                {
                    for (Integer j = 0; j < 8; ++j)
                    {
                        Ji(i, j) = N[i] * N[j];
                    }
                }

                Ji *= wxi * weta * wzeta * fabs(detJ);

                this->ddt += Ji;
            }

            if (this->m_dt > 0.0)
                this->ddt /= this->m_dt;
        }

        template <typename Real>
        void CFemHexahedron<Real, 8, 1, 1>::setDiffusionTerm()
        {
            Eigen::Matrix<Real, 3, 3> J;
            Eigen::Matrix<Real, 3, 8> G;
            Eigen::Matrix<Real, 3, 8> B;
            Eigen::Matrix<Real, 8, 8> Ki;

            this->laplacian.resize(8, 8);
            this->laplacian.setZero();

            for (Integer p = 0; p < CIntGaussIntegration<Real>::m_integPoints; p++)
            {
                Real xi = this->m_xi[p];
                Real eta = this->m_eta[p];
                Real zeta = this->m_zeta[p];

                Real wxi = this->m_wxi[p];
                Real weta = this->m_weta[p];
                Real wzeta = this->m_wzeta[p];

                Real xim = (1 - xi);
                Real xip = (1 + xi);

                Real etam = (1 - eta);
                Real etap = (1 + eta);

                Real zetam = (1 - zeta);
                Real zetap = (1 + zeta);

                J << -etam * zetam * m_x1 + etam * zetam * m_x2 - etap * zetam * m_x3 + etap * zetam * m_x4 - etam * zetap * m_x5 + etam * zetap * m_x6 - etap * zetap * m_x7 + etap * zetap * m_x8,
                    -etam * zetam * m_y1 + etam * zetam * m_y2 - etap * zetam * m_y3 + etap * zetam * m_y4 - etam * zetap * m_y5 + etam * zetap * m_y6 - etap * zetap * m_y7 + etap * zetap * m_y8,
                    -etam * zetam * m_z1 + etam * zetam * m_z2 - etap * zetam * m_z3 + etap * zetam * m_z4 - etam * zetap * m_z5 + etam * zetap * m_z6 - etap * zetap * m_z7 + etap * zetap * m_z8,
                    -xim * zetam * m_x1 - xip * zetam * m_x2 + xim * zetam * m_x3 + xip * zetam * m_x4 - xim * zetap * m_x5 - xip * zetap * m_x6 + xim * zetap * m_x7 + xip * zetap * m_x8,
                    -xim * zetam * m_y1 - xip * zetam * m_y2 + xim * zetam * m_y3 + xip * zetam * m_y4 - xim * zetap * m_y5 - xip * zetap * m_y6 + xim * zetap * m_y7 + xip * zetap * m_y8,
                    -xim * zetam * m_z1 - xip * zetam * m_z2 + xim * zetam * m_z3 + xip * zetam * m_z4 - xim * zetap * m_z5 - xip * zetap * m_z6 + xim * zetap * m_z7 + xip * zetap * m_z8,
                    -etam * xim * m_x1 - etam * xip * m_x2 - etap * xim * m_x3 - etap * xip * m_x4 + etam * xim * m_x5 + etam * xip * m_x6 + etap * xim * m_x7 + etap * xip * m_x8,
                    -etam * xim * m_y1 - etam * xip * m_y2 - etap * xim * m_y3 - etap * xip * m_y4 + etam * xim * m_y5 + etam * xip * m_y6 + etap * xim * m_y7 + etap * xip * m_y8,
                    -etam * xim * m_z1 - etam * xip * m_z2 - etap * xim * m_z3 - etap * xip * m_z4 + etam * xim * m_z5 + etam * xip * m_z6 + etap * xim * m_z7 + etap * xip * m_z8;

                J *= 0.125;

                Real detJ = J.determinant();

                G << -etam * zetam, +etam * zetam, -etap * zetam, +etap * zetam, -etam * zetap, +etam * zetap, -etap * zetap, +etap * zetap,
                    -xim * zetam, -xip * zetam, +xim * zetam, +xip * zetam, -xim * zetap, -xip * zetap, +xim * zetap, +xip * zetap,
                    -etam * xim, -etam * xip, -etap * xim, -etap * xip, +etam * xim, +etam * xip, +etap * xim, +etap * xip;

                G *= 0.125;

                B = J.inverse() * G;

                Ki = -B.transpose() * B;

                Ki *= wxi * weta * wzeta * fabs(detJ);

                this->laplacian += Ki;
            }
        }

        template <typename Real>
        void CFemHexahedron<Real, 8, 1, 1>::setConvectiveTerm()
        {
        }

        template <typename Real>
        void CFemHexahedron<Real, 8, 1, 1>::update()
        {
            // Calculate geometrical properties
            rebuild();

            // Integration
            this->setGaussPoints();

            // Diffusion term
            setDiffusionTerm();

            // Transient term
            if (this->m_transient)
                setTransientTerm();
        }

        template <typename Real>
        void CFemHexahedron<Real, 8, 1, 1>::setSourceOnNode(Integer aNodeIndex, Real aValue)
        {
            this->source(aNodeIndex) += aValue;
        }

        template <typename Real>
        void CFemHexahedron<Real, 8, 1, 1>::setSourceOnEdge(const Integer anEdgeIndex, const Real aValue)
        {
            this->source((anEdgeIndex + 0) % 8) += aValue * 0.5;
            this->source((anEdgeIndex + 1) % 8) += aValue * 0.5;
        }

        template <typename Real>
        void CFemHexahedron<Real, 8, 1, 1>::setSourceOnFace(const Integer aFaceIndex, const Real aValue)
        {
            //this->m_source += (aValue * 0.5);
        }
    }
}
