// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "GeoTriangle.hpp"

namespace ENigMA {
namespace geometry {
    template <typename Real>
    CGeoTetrahedron<Real>::CGeoTetrahedron()
    {
    }

    template <typename Real>
    CGeoTetrahedron<Real>::~CGeoTetrahedron()
    {
    }

    template <typename Real>
    void CGeoTetrahedron<Real>::reset()
    {
        // TODO:
    }

    template <typename Real>
    void CGeoTetrahedron<Real>::calculateCentroid(bool bReCalculate)
    {
        if (!this->m_bCentroid || bReCalculate) {
            this->m_centroid = this->m_vertices[0];
            this->m_centroid += this->m_vertices[1];
            this->m_centroid += this->m_vertices[2];
            this->m_centroid += this->m_vertices[3];
            this->m_centroid /= 4.0;

            this->m_bCentroid = true;
        }
    }

    template <typename Real>
    void CGeoTetrahedron<Real>::calculateSurfaceArea(bool bReCalculate)
    {
        if (!this->m_bSurfaceArea || bReCalculate) {
            CGeoTriangle<Real> aTriangle;

            aTriangle.addVertex(this->m_vertices[0]);
            aTriangle.addVertex(this->m_vertices[1]);
            aTriangle.addVertex(this->m_vertices[2]);

            aTriangle.calculateArea();

            this->m_surfaceArea = aTriangle.area();

            for (Integer i = 0; i < 3; ++i) {
                aTriangle.reset();
                aTriangle.addVertex(this->m_vertices[(i + 0) % 3]);
                aTriangle.addVertex(this->m_vertices[(i + 1) % 3]);
                aTriangle.addVertex(this->m_vertices[3]);

                aTriangle.calculateArea();

                this->m_surfaceArea += aTriangle.area();
            }

            this->m_bSurfaceArea = true;
        }
    }

    template <typename Real>
    void CGeoTetrahedron<Real>::calculateVolume(bool bReCalculate)
    {
        if (!this->m_bVolume || bReCalculate) {
            Eigen::Matrix<Real, 4, 4> matrix;

            matrix.col(0) << this->m_vertices[0].x(), this->m_vertices[0].y(), this->m_vertices[0].z(), 1.0;
            matrix.col(1) << this->m_vertices[1].x(), this->m_vertices[1].y(), this->m_vertices[1].z(), 1.0;
            matrix.col(2) << this->m_vertices[2].x(), this->m_vertices[2].y(), this->m_vertices[2].z(), 1.0;
            matrix.col(3) << this->m_vertices[3].x(), this->m_vertices[3].y(), this->m_vertices[3].z(), 1.0;

            this->m_volume = matrix.determinant() / 6.0;

            this->m_bVolume = true;
        }
    }

    template <typename Real>
    void CGeoTetrahedron<Real>::calculateBoundingBox(bool bReCalculate)
    {
        if (!this->m_bBoundingBox || bReCalculate) {
            this->m_boundingBox.reset();

            this->m_boundingBox.addCoordinate(this->m_vertices[0]);
            this->m_boundingBox.addCoordinate(this->m_vertices[1]);
            this->m_boundingBox.addCoordinate(this->m_vertices[2]);
            this->m_boundingBox.addCoordinate(this->m_vertices[3]);

            this->m_bBoundingBox = true;
        }
    }

    template <typename Real>
    bool CGeoTetrahedron<Real>::contains(const CGeoCoordinate<Real>& aPoint, const Real aTolerance)
    {
        this->calculateBoundingBox();

        if (!this->boundingBox().contains(aPoint, aTolerance))
            return false;

        Eigen::Matrix<Real, 4, 4> d0, d1, d2, d3, d4;

        d0.col(0) << this->m_vertices[0].x(), this->m_vertices[0].y(), this->m_vertices[0].z(), 1.0;
        d0.col(1) << this->m_vertices[1].x(), this->m_vertices[1].y(), this->m_vertices[1].z(), 1.0;
        d0.col(2) << this->m_vertices[2].x(), this->m_vertices[2].y(), this->m_vertices[2].z(), 1.0;
        d0.col(3) << this->m_vertices[3].x(), this->m_vertices[3].y(), this->m_vertices[3].z(), 1.0;

        Real det0 = d0.determinant();

        if (det0 <= 0.0)
            return false;

        d1.col(0) << aPoint.x(), aPoint.y(), aPoint.z(), 1.0;
        d1.col(1) << this->m_vertices[1].x(), this->m_vertices[1].y(), this->m_vertices[1].z(), 1.0;
        d1.col(2) << this->m_vertices[2].x(), this->m_vertices[2].y(), this->m_vertices[2].z(), 1.0;
        d1.col(3) << this->m_vertices[3].x(), this->m_vertices[3].y(), this->m_vertices[3].z(), 1.0;

        Real det1 = d1.determinant();

        Real s = det1 / det0;

        if (s < 0.0)
            return false;

        d2.col(0) << this->m_vertices[0].x(), this->m_vertices[0].y(), this->m_vertices[0].z(), 1.0;
        d2.col(1) << aPoint.x(), aPoint.y(), aPoint.z(), 1.0;
        d2.col(2) << this->m_vertices[2].x(), this->m_vertices[2].y(), this->m_vertices[2].z(), 1.0;
        d2.col(3) << this->m_vertices[3].x(), this->m_vertices[3].y(), this->m_vertices[3].z(), 1.0;

        Real det2 = d2.determinant();

        Real t = det2 / det0;

        if (t < 0.0)
            return false;

        d3.col(0) << this->m_vertices[0].x(), this->m_vertices[0].y(), this->m_vertices[0].z(), 1.0;
        d3.col(1) << this->m_vertices[1].x(), this->m_vertices[1].y(), this->m_vertices[1].z(), 1.0;
        d3.col(2) << aPoint.x(), aPoint.y(), aPoint.z(), 1.0;
        d3.col(3) << this->m_vertices[3].x(), this->m_vertices[3].y(), this->m_vertices[3].z(), 1.0;

        Real det3 = d3.determinant();

        Real u = det3 / det0;

        if (u < 0.0)
            return false;

        d4.col(0) << this->m_vertices[0].x(), this->m_vertices[0].y(), this->m_vertices[0].z(), 1.0;
        d4.col(1) << this->m_vertices[1].x(), this->m_vertices[1].y(), this->m_vertices[1].z(), 1.0;
        d4.col(2) << this->m_vertices[2].x(), this->m_vertices[2].y(), this->m_vertices[2].z(), 1.0;
        d4.col(3) << aPoint.x(), aPoint.y(), aPoint.z(), 1.0;

        Real det4 = d4.determinant();

        Real v = det4 / det0;

        if (v < 0.0)
            return false;

        Real x = s + t + u + v;

        if (x < -aTolerance || x > 1.0 + aTolerance)
            return false;

        return true;
    }

    template <typename Real>
    void CGeoTetrahedron<Real>::invert()
    {
        // TODO:
    }
}
}
