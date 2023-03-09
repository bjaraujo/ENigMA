// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

namespace ENigMA
{
    namespace stl
    {
        template <typename Real>
        CStlEdge<Real>::CStlEdge()
            : m_neighbor(-1)
            , m_whichVertexNot(-1)
            , m_outline(false)
            , m_naked(false)
        {
        }

        template <typename Real>
        CStlEdge<Real>::~CStlEdge()
        {
        }

        template <typename Real>
        void CStlEdge<Real>::setOutline(const bool aValue)
        {
            m_outline = aValue;
        }

        template <typename Real>
        void CStlEdge<Real>::setNaked(const bool aValue)
        {
            m_naked = aValue;
        }

        template <typename Real>
        bool CStlEdge<Real>::outline()
        {
            return m_outline;
        }

        template <typename Real>
        bool CStlEdge<Real>::naked()
        {
            return m_naked;
        }

        template <typename Real>
        void CStlEdge<Real>::setNeighbor(const Integer aFacetId)
        {
            m_neighbor = aFacetId;
        }

        template <typename Real>
        void CStlEdge<Real>::setWhichVertexNot(const Integer anIndex)
        {
            m_whichVertexNot = anIndex;
        }

        template <typename Real>
        int CStlEdge<Real>::neighbor()
        {
            return m_neighbor;
        }

        template <typename Real>
        int CStlEdge<Real>::whichVertexNot()
        {
            return m_whichVertexNot;
        }

        template <typename Real>
        CStlFacet<Real>::CStlFacet()
            : extra({ 0, 0 })
        {
        }

        template <typename Real>
        CStlFacet<Real>::~CStlFacet()
        {
        }

        template <typename Real>
        CStlEdge<Real>& CStlFacet<Real>::edge(const Integer anEdgeIndex)
        {
            return m_edge[anEdgeIndex];
        }

        template <typename Real>
        Integer CStlFile<Real>::nbFacets() const
        {
            return static_cast<Integer>(m_facetIds.size());
        }

        template <typename Real>
        Integer CStlFile<Real>::facetId(const Integer aFacetIndex) const
        {
            return m_facetIds.at(aFacetIndex);
        }

        template <typename Real>
        CStlFacet<Real>& CStlFile<Real>::facet(const Integer aFacetId)
        {
            return m_facets.at(aFacetId);
        }

        template <typename Real>
        void CStlFile<Real>::addFacet(const Integer aFacetId, CStlFacet<Real>& aFacet)
        {
            m_facets[aFacetId] = aFacet;
            m_facetIds.emplace_back(aFacetId);
        }

        template <typename Real>
        void CStlFile<Real>::removeFacet(const Integer aFacetId)
        {
            m_facets.erase(aFacetId);
            m_facetIds.erase(std::find(m_facetIds.begin(), m_facetIds.end(), aFacetId));
        }

        template <typename Real>
        Integer CStlFile<Real>::nextFacetId() const
        {
            if (m_facets.size() > 0)
                return m_facets.rbegin()->first + 1;
            else
                return 0;
        }

        template <typename Real>
        void CStlFile<Real>::reset()
        {
            m_facetIds.clear();
            m_facets.clear();
        }
    }
}
