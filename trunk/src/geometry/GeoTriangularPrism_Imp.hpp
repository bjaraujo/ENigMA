// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "GeoVector.hpp"

namespace ENigMA {
namespace geometry {
    template <typename Real>
    CGeoTriangularPrism<Real>::CGeoTriangularPrism()
    {
    }

    template <typename Real>
    CGeoTriangularPrism<Real>::~CGeoTriangularPrism()
    {
    }

    template <typename Real>
    void CGeoTriangularPrism<Real>::reset()
    {
        // TODO:
    }

    template <typename Real>
    void CGeoTriangularPrism<Real>::calculateCentroid(bool bReCalculate)
    {
        // TODO
    }

    template <typename Real>
    void CGeoTriangularPrism<Real>::calculateSurfaceArea(bool bReCalculate)
    {
        // TODO
    }

    template <typename Real>
    void CGeoTriangularPrism<Real>::calculateVolume(bool bReCalculate)
    {
        if (!this->m_bVolume || bReCalculate) {
            CGeoVector<Real> v1, v2, v3;

            v1 = this->m_vertices[1] - this->m_vertices[0];
            v2 = this->m_vertices[2] - this->m_vertices[0];
            v3 = this->m_vertices[3] - this->m_vertices[0];

            m_volume = 0.5 * v1.dot(v2.cross(v3));

            this->m_bVolume = true;
        }
    }

    template <typename Real>
    void CGeoTriangularPrism<Real>::calculateBoundingBox(bool bReCalculate)
    {
        // TODO:
    }
}
}
