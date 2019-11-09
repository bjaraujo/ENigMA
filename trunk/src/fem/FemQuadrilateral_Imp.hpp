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
    CFemQuadrilateral<Real, 4, 1, 1>::CFemQuadrilateral()
    {

        this->source.resize(4);
        this->source.setZero();

        this->m_thickness = 1.0;

        this->m_integPoints = 4;
    }

    template <typename Real>
    CFemQuadrilateral<Real, 4, 1, 1>::~CFemQuadrilateral()
    {
    }

    template <typename Real>
    void CFemQuadrilateral<Real, 4, 1, 1>::rebuild()
    {

        CGeoCoordinate<Real> n1 = this->m_vertices[0];
        CGeoCoordinate<Real> n2 = this->m_vertices[1];
        CGeoCoordinate<Real> n3 = this->m_vertices[2];
        CGeoCoordinate<Real> n4 = this->m_vertices[3];

        if (n1.z() != 0.0 || n2.z() != 0.0 || n3.z() != 0.0) {

            CGeoVector<Real> v1 = n2 - n1;
            CGeoVector<Real> v2 = n3 - n1;

            CGeoVector<Real> vx = v1;
            vx.normalize();

            CGeoVector<Real> vz = vx.cross(v2);
            vz.normalize();

            CGeoVector<Real> vy = vz.cross(vx);
            vy.normalize();

            CGeoCoordinateSystem<Real> aTransform;
            aTransform.col(0) << vx;
            aTransform.col(1) << vy;
            aTransform.col(2) << vz;

            aTransform.transposeInPlace();

            n1.transform(aTransform);
            n2.transform(aTransform);
            n3.transform(aTransform);
            n4.transform(aTransform);
        }

        m_x1 = n1.x();
        m_x2 = n2.x();
        m_x3 = n3.x();
        m_x4 = n4.x();

        m_y1 = n1.y();
        m_y2 = n2.y();
        m_y3 = n3.y();
        m_y4 = n4.y();

        m_x21 = m_x2 - m_x1;
        m_y21 = m_y2 - m_y1;

        m_x32 = m_x3 - m_x2;
        m_y32 = m_y3 - m_y2;

        m_x43 = m_x4 - m_x3;
        m_y43 = m_y4 - m_y3;

        m_x14 = m_x1 - m_x4;
        m_y14 = m_y1 - m_y4;
    }

    template <typename Real>
    void CFemQuadrilateral<Real, 4, 1, 1>::setTransientTerm()
    {

        Real N[4];

        Eigen::Matrix<Real, 2, 2> J;
        Eigen::Matrix<Real, 4, 4> Ji;

        this->ddt.resize(4, 4);
        this->ddt.setZero();

        for (Integer p = 0; p < this->m_integPoints; p++) {

            Real xi = this->m_xi[p];
            Real eta = this->m_eta[p];

            Real wxi = this->m_wxi[p];
            Real weta = this->m_weta[p];

            Real xim = (1 - xi);
            Real xip = (1 + xi);
            Real etam = (1 - eta);
            Real etap = (1 + eta);

            J.resize(2, 2);

            J << -xim * m_x1 + xim * m_x2 + xip * m_x3 - xip * m_x4, -xim * m_y1 + xim * m_y2 + xip * m_y3 - xip * m_y4,
                -etam * m_x1 - etap * m_x2 + etap * m_x3 + etam * m_x4, -etam * m_y1 - etap * m_y2 + etap * m_y3 + etam * m_y4;

            J *= 0.25;

            Real detJ = J.determinant();

            N[0] = 0.25 * (1 - xi) * (1 - eta);
            N[1] = 0.25 * (1 + xi) * (1 - eta);
            N[2] = 0.25 * (1 + xi) * (1 + eta);
            N[3] = 0.25 * (1 - xi) * (1 + eta);

            Ji.setZero();

            for (Integer i = 0; i < 4; ++i) {
                for (Integer j = 0; j < 4; ++j) {

                    Ji(i, j) = N[i] * N[j];
                }
            }

            Ji *= wxi * weta * fabs(detJ);

            this->ddt += Ji;
        }

        this->ddt *= this->m_thickness;

        if (this->m_dt > 0.0)
            this->ddt /= this->m_dt;
    }

    template <typename Real>
    void CFemQuadrilateral<Real, 4, 1, 1>::setDiffusionTerm()
    {

        Eigen::Matrix<Real, 2, 2> J;
        Eigen::Matrix<Real, 2, 4> G;
        Eigen::Matrix<Real, 2, 4> B;
        Eigen::Matrix<Real, 4, 4> Ki;

        this->laplacian.resize(4, 4);
        this->laplacian.setZero();

        for (Integer p = 0; p < this->m_integPoints; p++) {

            Real xi = this->m_xi[p];
            Real eta = this->m_eta[p];

            Real wxi = this->m_wxi[p];
            Real weta = this->m_weta[p];

            Real xim = (1 - xi);
            Real xip = (1 + xi);
            Real etam = (1 - eta);
            Real etap = (1 + eta);

            J.resize(2, 2);

            J << -xim * m_x1 + xim * m_x2 + xip * m_x3 - xip * m_x4, -xim * m_y1 + xim * m_y2 + xip * m_y3 - xip * m_y4,
                -etam * m_x1 - etap * m_x2 + etap * m_x3 + etam * m_x4, -etam * m_y1 - etap * m_y2 + etap * m_y3 + etam * m_y4;

            J *= 0.25;

            Real detJ = J.determinant();

            G << -xim, +xim, +xip, -xip,
                -etam, -etap, +etap, +etam;

            G *= 0.25;

            B = J.inverse() * G;

            Ki = -B.transpose() * B;

            Ki *= wxi * weta * fabs(detJ);

            this->laplacian += Ki;
        }

        this->laplacian *= this->m_thickness;
    }

    template <typename Real>
    void CFemQuadrilateral<Real, 4, 1, 1>::setConvectiveTerm()
    {
    }

    template <typename Real>
    void CFemQuadrilateral<Real, 4, 1, 1>::update()
    {

        // Calculate geometrical properties
        this->rebuild();

        // Integration
        this->setGaussPoints();

        // Diffusion term
        this->setDiffusionTerm();

        // Transient term
        if (this->m_transient)
            this->setTransientTerm();
    }

    template <typename Real>
    void CFemQuadrilateral<Real, 4, 1, 1>::setSourceOnNode(const Integer aNodeIndex, const Real aValue)
    {

        this->source(aNodeIndex) += aValue;
    }

    template <typename Real>
    void CFemQuadrilateral<Real, 4, 1, 1>::setSourceOnEdge(const Integer anEdgeIndex, const Real aValue)
    {

        this->source((anEdgeIndex + 0) % 4) += aValue * 0.5;
        this->source((anEdgeIndex + 1) % 4) += aValue * 0.5;
    }
}
}
