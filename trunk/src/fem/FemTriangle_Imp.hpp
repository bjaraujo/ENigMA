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
    CFemTriangle<Real, 3, 1, 1>::CFemTriangle()
    {

        this->source.resize(3);
        this->source.setZero();

        this->m_thickness = 1.0;
    }

    template <typename Real>
    CFemTriangle<Real, 3, 1, 1>::~CFemTriangle()
    {
    }

    template <typename Real>
    void CFemTriangle<Real, 3, 1, 1>::rebuild()
    {

        CGeoCoordinate<Real> n1 = this->m_vertices[0];
        CGeoCoordinate<Real> n2 = this->m_vertices[1];
        CGeoCoordinate<Real> n3 = this->m_vertices[2];

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
        }

        m_x1 = n1.x();
        m_x2 = n2.x();
        m_x3 = n3.x();

        m_y1 = n1.y();
        m_y2 = n2.y();
        m_y3 = n3.y();

        m_x21 = m_x2 - m_x1;
        m_y21 = m_y2 - m_y1;

        m_x32 = m_x3 - m_x2;
        m_y32 = m_y3 - m_y2;

        m_x13 = m_x1 - m_x3;
        m_y13 = m_y1 - m_y3;
    }

    template <typename Real>
    void CFemTriangle<Real, 3, 1, 1>::setTransientTerm()
    {

        this->ddt.resize(3, 3);

        this->ddt << 1.0 / 6.0, 1.0 / 12.0, 1.0 / 12.0,
            1.0 / 12.0, 1.0 / 6.0, 1.0 / 12.0,
            1.0 / 12.0, 1.0 / 12.0, 1.0 / 6.0;

        this->ddt *= CGeoTriangle<Real>::area();

        if (this->m_dt > 0.0)
            this->ddt /= this->m_dt;
    }

    template <typename Real>
    void CFemTriangle<Real, 3, 1, 1>::calculateB(Eigen::Matrix<Real, 2, 3>& B)
    {

        B << -m_y32, -m_y13, -m_y21,
            +m_x32, +m_x13, +m_x21;
    }

    template <typename Real>
    void CFemTriangle<Real, 3, 1, 1>::setDiffusionTerm()
    {

        this->laplacian.resize(3, 3);

        Eigen::Matrix<Real, 2, 3> B;

        this->calculateB(B);

        if (CGeoTriangle<Real>::area() > 0.0)
            B *= 1.0 / (2.0 * CGeoTriangle<Real>::area());

        CFemElement<Real>::laplacian = -B.transpose() * B;
        CFemElement<Real>::laplacian *= CGeoTriangle<Real>::area() * this->m_thickness;
    }

    template <typename Real>
    void CFemTriangle<Real, 3, 1, 1>::setConvectiveTerm()
    {
    }

    template <typename Real>
    void CFemTriangle<Real, 3, 1, 1>::update()
    {

        // Calculate geometrical properties
        this->rebuild();

        // Diffusion term
        this->setDiffusionTerm();

        // Transient term
        if (this->m_transient)
            this->setTransientTerm();
    }

    template <typename Real>
    void CFemTriangle<Real, 3, 1, 1>::setSourceOnNode(Integer aNodeIndex, Real aValue)
    {

        this->source(aNodeIndex) += aValue;
    }

    template <typename Real>
    void CFemTriangle<Real, 3, 1, 1>::setSourceOnEdge(const Integer anEdgeIndex, const Real aValue)
    {

        this->source((anEdgeIndex + 0) % 3) += aValue * 0.5;
        this->source((anEdgeIndex + 1) % 3) += aValue * 0.5;
    }

    template <typename Real>
    void CFemTriangle<Real, 3, 2, 1>::calculateB(Eigen::Matrix<Real, 3, 6>& B)
    {

        Real x21 = this->m_x21;
        Real x32 = this->m_x32;
        Real x13 = this->m_x13;

        Real y21 = this->m_y21;
        Real y32 = this->m_y32;
        Real y13 = this->m_y13;

        B << -y32, 0.0, -y13, 0.0, -y21, 0.0,
            0.0, +x32, 0.0, +x13, 0.0, +x21,
            +x32, -y32, +x13, -y13, +x21, -y21;
    }
}
}
