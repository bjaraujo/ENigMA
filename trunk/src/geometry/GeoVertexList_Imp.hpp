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
    CGeoVertexList<Real>::CGeoVertexList()
    {
    }

    template <typename Real>
    CGeoVertexList<Real>::~CGeoVertexList()
    {
        this->reset();
    }

    template <typename Real>
    void CGeoVertexList<Real>::reset()
    {
        m_vertices.clear();
    }

    template <typename Real>
    Integer CGeoVertexList<Real>::nbVertices()
    {
        return static_cast<Integer>(m_vertices.size());
    }

    template <typename Real>
    void CGeoVertexList<Real>::addVertex(CGeoCoordinate<Real>& aVertex)
    {
        m_vertices.push_back(aVertex);
    }

    template <typename Real>
    void CGeoVertexList<Real>::insertVertex(const Integer aVertexIndex, CGeoCoordinate<Real>& aVertex)
    {
        m_vertices.insert(m_vertices.begin() + aVertexIndex, aVertex);
    }

    template <typename Real>
    void CGeoVertexList<Real>::removeVertex(const Integer aVertexIndex)
    {
        m_vertices.erase(m_vertices.begin() + aVertexIndex);
    }

    template <typename Real>
    CGeoCoordinate<Real>& CGeoVertexList<Real>::vertex(const Integer aVertexIndex)
    {
        if (m_vertices.size() == 0)
        {
            throw std::runtime_error("Vertex size = 0");
        }

        return m_vertices[aVertexIndex % m_vertices.size()];
    }

    template <typename Real>
    void CGeoVertexList<Real>::removeDuplicates()
    {
        m_vertices.erase(std::unique(m_vertices.begin(), m_vertices.end()), m_vertices.end());
    }

    template <typename Real>
    void CGeoVertexList<Real>::removeCollinear(const Real aTolerance)
    {
        if (m_vertices.size() == 0)
        {
            return;
        }
        
        Integer j = 1;

        for (Integer i = 1; i < static_cast<Integer>(m_vertices.size()) - 1; ++i) {

            CGeoVector<Real> v1, v2;

            v1 = m_vertices[(i + 1) % m_vertices.size()] - m_vertices[(i + 0) % m_vertices.size()];
            v2 = m_vertices[(i + 0) % m_vertices.size()] - m_vertices[(i - 1) % m_vertices.size()];

            Real angle = v1.angle(v2);

            if (angle > aTolerance) {
                // Keep
                m_vertices[j++] = m_vertices[i];
            }
        }

        m_vertices[j++] = m_vertices[m_vertices.size() - 1];

        m_vertices.resize(j);
    }

    template <typename Real>
    void CGeoVertexList<Real>::invert()
    {
        std::reverse(m_vertices.begin(), m_vertices.end());
    }
}
}
