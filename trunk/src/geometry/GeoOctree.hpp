// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "GeoBoundingBox.hpp"
#include "GeoContainer.hpp"
#include "GeoCoordinate.hpp"
#include "GeoOctreeNode.hpp"

namespace ENigMA {

namespace geometry {

#define MAX_NB_ELEMENTS_PER_NODE 1000

    enum EOctreeLocation {
        OL_BACK_BOTTOM_LEFT = 0,
        OL_BACK_BOTTOM_RIGHT = 1,
        OL_FRONT_BOTTOM_LEFT = 2,
        OL_FRONT_BOTTOM_RIGHT = 3,
        OL_BACK_TOP_LEFT = 4,
        OL_BACK_TOP_RIGHT = 5,
        OL_FRONT_TOP_LEFT = 6,
        OL_FRONT_TOP_RIGHT = 7,
    };

    template <typename Real>
    class CGeoOctree : public CGeoContainer<CGeoCoordinate<Real>, Real> {
    private:
        CGeoOctreeNode<Real> m_root;

        EOctreeLocation dispatch(CGeoOctreeNode<Real>& anOctreeNode, Real x, Real y, Real z);
        void createNode(CGeoOctreeNode<Real>& anOctreeNode, const CGeoOctreeNode<Real>& anOctreeParent, Integer i);
        void createRecursive(CGeoOctreeNode<Real>& anOctreeNode);
        CGeoOctreeNode<Real>& findRecursive(CGeoOctreeNode<Real>& anOctreeNode, CGeoCoordinate<Real>& aCoordinate);

        CGeoBoundingBox<Real> m_boundingBox;

    public:
        CGeoOctree();
        virtual ~CGeoOctree();

        void reset();

        void build();

        void find(std::vector<Integer>& coordinateIds, CGeoCoordinate<Real>& aCoordinate, const Real aTolerance = 0.0);
    };
}
}

#include "GeoOctree_Imp.hpp"
