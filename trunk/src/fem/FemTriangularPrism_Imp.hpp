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

namespace ENigMA {

namespace fem {

    template <typename Real>
    CFemTriangularPrism<Real, 6, 1, 1>::CFemTriangularPrism()
    {

        this->source.resize(6);
        this->source.setZero();

        CIntGaussIntegration<Real>::m_integPoints = 6;
    }

    template <typename Real>
    CFemTriangularPrism<Real, 6, 1, 1>::~CFemTriangularPrism()
    {
    }

    template <typename Real>
    void CFemTriangularPrism<Real, 6, 1, 1>::rebuild()
    {

        CGeoCoordinate<Real> n1 = this->m_vertices[0];
        CGeoCoordinate<Real> n2 = this->m_vertices[1];
        CGeoCoordinate<Real> n3 = this->m_vertices[2];
        CGeoCoordinate<Real> n4 = this->m_vertices[3];
        CGeoCoordinate<Real> n5 = this->m_vertices[4];
        CGeoCoordinate<Real> n6 = this->m_vertices[5];

        m_x1 = n1.x();
        m_x2 = n2.x();
        m_x3 = n3.x();
        m_x4 = n4.x();
        m_x5 = n5.x();
        m_x6 = n6.x();

        m_y1 = n1.y();
        m_y2 = n2.y();
        m_y3 = n3.y();
        m_y4 = n4.y();
        m_y5 = n5.y();
        m_y6 = n6.y();

        m_z1 = n1.z();
        m_z2 = n2.z();
        m_z3 = n3.z();
        m_z4 = n4.z();
        m_z5 = n5.z();
        m_z6 = n6.z();
    }

    template <typename Real>
    void CFemTriangularPrism<Real, 6, 1, 1>::setTransientTerm()
    {

        Real N[6];

        Eigen::Matrix<Real, 3, 3> J;
        Eigen::Matrix<Real, 6, 6> Ji;

        this->ddt.resize(6, 6);
        this->ddt.setZero();

        for (Integer p = 0; p < CIntGaussIntegration<Real>::m_integPoints; p++) {

            Real xi = this->m_xi[p];
            Real eta = this->m_eta[p];
            Real zeta = this->m_zeta[p];

            Real xieta = -xi - eta + 1;
            Real zetam = (1 - zeta);
            Real zetap = (1 + zeta);

            J << m_x5 * zetap - m_x4 * zetap + m_x2 * zetam - m_x1 * zetam,
                m_y5 * zetap - m_y4 * zetap + m_y2 * zetam - m_y1 * zetam,
                m_z5 * zetap - m_z4 * zetap + m_z2 * zetam - m_z1 * zetam,
                m_x6 * zetap - m_x4 * zetap + m_x3 * zetam - m_x1 * zetam,
                m_y6 * zetap - m_y4 * zetap + m_y3 * zetam - m_y1 * zetam,
                m_z6 * zetap - m_z4 * zetap + m_z3 * zetam - m_z1 * zetam,
                m_x5 * xi - m_x2 * xi + m_x4 * xieta - m_x1 * xieta + eta * m_x6 - eta * m_x3,
                eta * m_y6 + xi * m_y5 + xieta * m_y4 - eta * m_y3 - xi * m_y2 - xieta * m_y1,
                eta * m_z6 + xi * m_z5 + xieta * m_z4 - eta * m_z3 - xi * m_z2 - xieta * m_z1;

            J *= 0.5;

            Real detJ = J.determinant();

            N[0] = 0.5 * (1 - xi - eta) * (1 - zeta);
            N[1] = 0.5 * xi * (1 - zeta);
            N[2] = 0.5 * eta * (1 - zeta);
            N[3] = 0.5 * (1 - xi - eta) * (1 + zeta);
            N[4] = 0.5 * xi * (1 + zeta);
            N[5] = 0.5 * eta * (1 + zeta);

            Ji.setZero();

            for (Integer i = 0; i < 6; ++i) {

                for (Integer j = 0; j < 6; ++j) {

                    Ji(i, j) = N[i] * N[j];
                }
            }

            Ji *= this->m_wxi[p] * this->m_weta[p] * this->m_wzeta[p] * fabs(detJ);

            this->ddt += Ji;
        }

        if (this->m_dt > 0.0)
            this->ddt /= this->m_dt;
    }

    template <typename Real>
    void CFemTriangularPrism<Real, 6, 1, 1>::setDiffusionTerm()
    {

        Eigen::Matrix<Real, 3, 3> J;
        Eigen::Matrix<Real, 3, 6> G;
        Eigen::Matrix<Real, 3, 6> B;
        Eigen::Matrix<Real, 6, 6> Ki;

        this->laplacian.resize(6, 6);
        this->laplacian.setZero();

        for (Integer p = 0; p < this->m_integPoints; p++) {

            Real xi = this->m_xi[p];
            Real eta = this->m_eta[p];
            Real zeta = this->m_zeta[p];

            Real wxi = this->m_wxi[p];
            Real weta = this->m_weta[p];
            Real wzeta = this->m_wzeta[p];

            Real xieta = -xi - eta + 1;
            Real zetam = (1 - zeta);
            Real zetap = (1 + zeta);

            J << m_x5 * zetap - m_x4 * zetap + m_x2 * zetam - m_x1 * zetam,
                m_y5 * zetap - m_y4 * zetap + m_y2 * zetam - m_y1 * zetam,
                m_z5 * zetap - m_z4 * zetap + m_z2 * zetam - m_z1 * zetam,
                m_x6 * zetap - m_x4 * zetap + m_x3 * zetam - m_x1 * zetam,
                m_y6 * zetap - m_y4 * zetap + m_y3 * zetam - m_y1 * zetam,
                m_z6 * zetap - m_z4 * zetap + m_z3 * zetam - m_z1 * zetam,
                m_x5 * xi - m_x2 * xi + m_x4 * xieta - m_x1 * xieta + eta * m_x6 - eta * m_x3,
                eta * m_y6 + xi * m_y5 + xieta * m_y4 - eta * m_y3 - xi * m_y2 - xieta * m_y1,
                eta * m_z6 + xi * m_z5 + xieta * m_z4 - eta * m_z3 - xi * m_z2 - xieta * m_z1;

            J *= 0.5;

            Real detJ = J.determinant();

            G << -zetam, +zetam, 0.0, -zetap, +zetap, 0.0,
                -zetam, 0.0, +zetam, -zetap, 0.0, +zetap,
                -xieta, -xi, -eta, +xieta, +xi, +eta;

            G *= 0.5;

            B = J.inverse() * G;

            Ki = -B.transpose() * B;

            Ki *= wxi * weta * wzeta * fabs(detJ);

            this->laplacian += Ki;
        }
    }

    template <typename Real>
    void CFemTriangularPrism<Real, 6, 1, 1>::setConvectiveTerm()
    {
    }

    template <typename Real>
    void CFemTriangularPrism<Real, 6, 1, 1>::update()
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
    void CFemTriangularPrism<Real, 6, 1, 1>::setSourceOnNode(const Integer aNodeIndex, const Real aValue)
    {

        this->source(aNodeIndex) += aValue;
    }

    template <typename Real>
    void CFemTriangularPrism<Real, 6, 1, 1>::setSourceOnEdge(const Integer anEdgeIndex, const Real aValue)
    {

        this->source((anEdgeIndex + 0) % 6) += aValue * 0.5;
        this->source((anEdgeIndex + 1) % 6) += aValue * 0.5;
    }

    template <typename Real>
    void CFemTriangularPrism<Real, 6, 1, 1>::setSourceOnFace(const Integer aFaceIndex, const Real aValue)
    {
    }
}
}
