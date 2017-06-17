// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include <map>

#include "GeoLine.hpp"
#include "GeoTriangle.hpp"

using namespace ENigMA::geometry;

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
        class CStlEdge : public CGeoLine<Real>
        {
        private:
            bool m_outline;
            bool m_naked;

            Integer m_neighbor;
            Integer m_whichVertexNot;

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
        class CStlFacet : public CGeoTriangle<Real>
        {
        private:
            CStlEdge<Real> m_edge[3];

        public:
            CStlFacet();
            ~CStlFacet();

            char extra[2];

            CStlEdge<Real>& edge(const Integer anEdgeIndex);
        };

        template <typename Real>
        class CStlStats
        {
        public:
            std::streamoff fileSize;
            char           header[81];
            EStlFileType   type;

            CGeoVector<Real> max;
            CGeoVector<Real> min;
            CGeoVector<Real> size;

        };

        template <typename Real>
        class CStlFile
        {
        private:

            std::vector<Integer> m_facetIds;

            typedef std::map<Integer, CStlFacet<Real> > mapFacets;
            mapFacets m_facets;

        public:
            std::string fileName;

            void reset();

            Integer nbFacets();

            Integer facetId(const Integer aFacetIndex);

            void addFacet(const Integer aFacetId, CStlFacet<Real>& aFacet);
            void removeFacet(const Integer aFacetId);

            CStlFacet<Real>& facet(const Integer aFacetId);

            Integer nextFacetId();

            CStlStats<Real> stats;

        };

    }

}

#include "StlFile_Imp.hpp"
