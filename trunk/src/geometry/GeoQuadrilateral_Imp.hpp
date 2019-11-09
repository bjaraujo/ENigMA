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
#include "GeoVector.hpp"

namespace ENigMA {

namespace geometry {

    template <typename Real>
    CGeoQuadrilateral<Real>::CGeoQuadrilateral()
    {
    }

    template <typename Real>
    CGeoQuadrilateral<Real>::~CGeoQuadrilateral()
    {
    }

    template <typename Real>
    void CGeoQuadrilateral<Real>::reset()
    {

        CGeoVertexList<Real>::reset();

        this->m_bCentroid = false;
        this->m_bNormal = false;
        this->m_bArea = false;
        this->m_bBoundingBox = false;
    }

    template <typename Real>
    void CGeoQuadrilateral<Real>::calculateCentroid(bool bReCalculate)
    {

        if (!this->m_bCentroid || bReCalculate) {

            this->m_centroid = this->m_vertices[0];
            this->m_centroid += this->m_vertices[1];
            this->m_centroid += this->m_vertices[2];
            this->m_centroid += this->m_vertices[3];

            this->m_centroid *= 0.25;

            this->m_bCentroid = true;
        }
    }

    template <typename Real>
    void CGeoQuadrilateral<Real>::calculateNormal(bool bReCalculate)
    {

        if (!this->m_bNormal || bReCalculate) {

            CGeoVector<Real> v0 = this->m_vertices[2] - this->m_vertices[0];
            CGeoVector<Real> v1 = this->m_vertices[3] - this->m_vertices[1];

            this->m_normal = v0.cross(v1);

            this->m_area = this->m_normal.norm() * 0.5;

            if (this->m_area > 0.0)
                this->m_normal.normalize();

            this->m_bNormal = true;
            this->m_bArea = true;
        }
    }

    template <typename Real>
    void CGeoQuadrilateral<Real>::calculateArea(bool bReCalculate)
    {

        this->calculateNormal(bReCalculate);
    }

    template <typename Real>
    void CGeoQuadrilateral<Real>::calculateBoundingBox(bool bReCalculate)
    {

        if (!this->m_bBoundingBox || bReCalculate) {

            this->m_boundingBox.reset();

            for (Integer k = 0; k < 4; ++k)
                this->m_boundingBox.addCoordinate(this->m_vertices[k]);

            this->m_bBoundingBox = true;
        }
    }

    template <typename Real>
    bool CGeoQuadrilateral<Real>::contains(const CGeoCoordinate<Real>& aPoint, CGeoIntersectionType& anIntersectionType, const Real aTolerance)
    {

        CGeoTriangle<Real> aTriangle1;

        aTriangle1.addVertex(this->m_vertices[0]);
        aTriangle1.addVertex(this->m_vertices[1]);
        aTriangle1.addVertex(this->m_vertices[2]);

        if (aTriangle1.contains(aPoint, anIntersectionType, aTolerance))
            return true;

        CGeoTriangle<Real> aTriangle2;

        aTriangle2.addVertex(this->m_vertices[0]);
        aTriangle2.addVertex(this->m_vertices[2]);
        aTriangle2.addVertex(this->m_vertices[3]);

        if (aTriangle2.contains(aPoint, anIntersectionType, aTolerance))
            return true;

        return false;
    }

    template <typename Real>
    void CGeoQuadrilateral<Real>::decimate(std::vector<CGeoTriangle<Real>>& sTriangles)
    {

        CGeoTriangle<Real> aTriangle;

        aTriangle.reset();
        aTriangle.addVertex(this->m_vertices[0]);
        aTriangle.addVertex(this->m_vertices[1]);
        aTriangle.addVertex(this->m_vertices[2]);

        sTriangles.push_back(aTriangle);

        aTriangle.reset();
        aTriangle.addVertex(this->m_vertices[0]);
        aTriangle.addVertex(this->m_vertices[2]);
        aTriangle.addVertex(this->m_vertices[3]);

        sTriangles.push_back(aTriangle);
    }
}
}
