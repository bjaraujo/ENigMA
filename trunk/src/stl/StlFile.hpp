// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include <array>
#include <map>

#include "GeoLine.hpp"
#include "GeoTriangle.hpp"

namespace ENigMA
{
    namespace stl
    {
        enum EStlFileType
        {
            FT_BINARY = 0,
            FT_ASCII
        };

        template <typename Real>
        class CStlEdge : public ENigMA::geometry::CGeoLine<Real>
        {
        private:
            Integer m_neighbor;
            Integer m_whichVertexNot;

            bool m_outline;
            bool m_naked;

        public:
            CStlEdge();
            ~CStlEdge();

            void setOutline(const bool aValue);
            void setNaked(const bool aValue);

            bool outline();
            bool naked();

            void setNeighbor(const Integer aFacetId);
            void setWhichVertexNot(const Integer anIndex);

            Integer neighbor();
            Integer whichVertexNot();
        };

        template <typename Real>
        class CStlFacet : public ENigMA::geometry::CGeoTriangle<Real>
        {
        private:
            CStlEdge<Real> m_edge[3];

        public:
            CStlFacet();
            ~CStlFacet();

            std::array<char, 2> extra;

            CStlEdge<Real>& edge(const Integer anEdgeIndex);
        };

        template <typename Real>
        class CStlStats
        {
        public:
            std::streamoff fileSize;
            char header[81];
            EStlFileType type;

            ENigMA::geometry::CGeoVector<Real> max;
            ENigMA::geometry::CGeoVector<Real> min;
            ENigMA::geometry::CGeoVector<Real> size;
        };

        template <typename Real>
        class CStlFile
        {
        private:
            std::vector<Integer> m_facetIds;

            typedef std::map<Integer, CStlFacet<Real>> mapFacets;
            mapFacets m_facets;

        public:
            std::string fileName;

            void reset();

            Integer nbFacets() const;

            Integer facetId(const Integer aFacetIndex) const;

            void addFacet(const Integer aFacetId, CStlFacet<Real>& aFacet);
            void removeFacet(const Integer aFacetId);

            CStlFacet<Real>& facet(const Integer aFacetId);

            Integer nextFacetId() const;

            CStlStats<Real> stats;
        };
    }
}

#include "StlFile_Imp.hpp"
