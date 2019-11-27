// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "GeoContainer.hpp"

namespace ENigMA {

namespace geometry {

    template <typename Real>
    class CGeoOctreeNode {
    private:
        std::vector<Integer> m_coordinateList;
        bool m_isLeaf;

    public:
        CGeoOctreeNode();
        virtual ~CGeoOctreeNode();

        Real xmin, xmax, ymin, ymax, zmin, zmax;
        Real xmid, ymid, zmid;

        std::vector<CGeoOctreeNode<Real>> nodes;

        Integer nbCoordinates() const;
        void addCoordinate(Integer aCoordinateIndex);
        Integer coordinate(Integer aCoordinateIndex);

        void setIsLeaf(bool isLeaf);
        bool isLeaf();

        void reset();
    };
}
}

#include "GeoOctreeNode_Imp.hpp"
