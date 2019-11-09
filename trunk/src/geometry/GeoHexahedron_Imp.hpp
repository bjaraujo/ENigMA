// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

// <ProjectName> $ProjectName$ </ProjectName>

#pragma once

#include "GeoVector.hpp"

namespace ENigMA {

namespace geometry {

    template <typename Real>
    CGeoHexahedron<Real>::CGeoHexahedron()
    {
    }

    template <typename Real>
    CGeoHexahedron<Real>::~CGeoHexahedron()
    {
    }

    template <typename Real>
    void CGeoHexahedron<Real>::reset()
    {

        // TODO
    }

    template <typename Real>
    void CGeoHexahedron<Real>::calculateCentroid(bool bReCalculate)
    {

        // TODO
    }

    template <typename Real>
    void CGeoHexahedron<Real>::calculateSurfaceArea(bool bReCalculate)
    {

        // TODO
    }

    template <typename Real>
    void CGeoHexahedron<Real>::calculateBoundingBox(bool bReCalculate)
    {

        // TODO
    }

    template <typename Real>
    void CGeoHexahedron<Real>::calculateVolume(bool bReCalculate)
    {

        if (!this->m_bVolume || bReCalculate) {
            CGeoVector<Real> v1, v2, v3;

            v1 = this->m_vertices[1] - this->m_vertices[0];
            v2 = this->m_vertices[3] - this->m_vertices[0];
            v3 = this->m_vertices[4] - this->m_vertices[0];

            CGeoVolume<Real>::volume() = v1.dot(v2.cross(v3));

            this->m_bVolume = true;
        }
    }

    template <typename Real>
    void CGeoHexahedron<Real>::decimate(std::vector<CGeoTetrahedron<Real>>& sTetrahedrons)
    {

        // Tetra 1
        CGeoTetrahedron<Real> aTetrahedron;

        aTetrahedron.reset();
        aTetrahedron.addVertex(this->m_vertices[0]);
        aTetrahedron.addVertex(this->m_vertices[1]);
        aTetrahedron.addVertex(this->m_vertices[3]);
        aTetrahedron.addVertex(this->m_vertices[7]);

        sTetrahedrons.push_back(aTetrahedron);

        // Tetra 2
        aTetrahedron.reset();
        aTetrahedron.addVertex(this->m_vertices[0]);
        aTetrahedron.addVertex(this->m_vertices[4]);
        aTetrahedron.addVertex(this->m_vertices[1]);
        aTetrahedron.addVertex(this->m_vertices[7]);

        sTetrahedrons.push_back(aTetrahedron);

        // Tetra 3
        aTetrahedron.reset();
        aTetrahedron.addVertex(this->m_vertices[1]);
        aTetrahedron.addVertex(this->m_vertices[4]);
        aTetrahedron.addVertex(this->m_vertices[5]);
        aTetrahedron.addVertex(this->m_vertices[7]);

        sTetrahedrons.push_back(aTetrahedron);

        // Tetra 4
        aTetrahedron.reset();
        aTetrahedron.addVertex(this->m_vertices[1]);
        aTetrahedron.addVertex(this->m_vertices[2]);
        aTetrahedron.addVertex(this->m_vertices[3]);
        aTetrahedron.addVertex(this->m_vertices[7]);

        sTetrahedrons.push_back(aTetrahedron);

        // Tetra 5
        aTetrahedron.reset();
        aTetrahedron.addVertex(this->m_vertices[1]);
        aTetrahedron.addVertex(this->m_vertices[6]);
        aTetrahedron.addVertex(this->m_vertices[2]);
        aTetrahedron.addVertex(this->m_vertices[7]);

        sTetrahedrons.push_back(aTetrahedron);

        // Tetra 6
        aTetrahedron.reset();
        aTetrahedron.addVertex(this->m_vertices[1]);
        aTetrahedron.addVertex(this->m_vertices[5]);
        aTetrahedron.addVertex(this->m_vertices[6]);
        aTetrahedron.addVertex(this->m_vertices[7]);

        sTetrahedrons.push_back(aTetrahedron);
    }
}
}
