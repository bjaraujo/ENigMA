// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include <vector>

#include "GeoBoundingBox.hpp"
#include "GeoContainer.hpp"
#include "GeoCoordinate.hpp"

namespace ENigMA {

namespace geometry {

    template <typename Real>
    class CGeoAdtreeNode6 {
    public:
        CGeoAdtreeNode6<Real>* leftNode;
        CGeoAdtreeNode6<Real>* rightNode;
        CGeoAdtreeNode6<Real>* fatherNode;

        Integer id;

        Real sep;
        Real data[6];

        Integer nbChildren;

        CGeoAdtreeNode6();

        void deleteChildren();
    };

    template <typename Real>
    class CGeoAdtree6 {
    protected:
        CGeoAdtreeNode6<Real>* m_root;

        Real m_cmin[6], m_cmax[6];

        std::vector<CGeoAdtreeNode6<Real>*> m_nodes;

    public:
        CGeoAdtree6();
        ~CGeoAdtree6();

        void reset();
        void set(const Real acmin[6], const Real acmax[6]);
        void insert(const Real p[6], const Integer anId);
        void remove(Integer anId);

        void getIntersecting(const Real bmin[6], const Real bmax[6], std::vector<Integer>& sIds);
    };

    template <typename Real>
    class CGeoAdtree : public CGeoContainer<CGeoBoundingBox<Real>, Real> {
    private:
        CGeoAdtree6<Real> m_tree;
        CGeoBoundingBox<Real> m_boundingBox;

    public:
        CGeoAdtree();
        ~CGeoAdtree();

        void set(CGeoBoundingBox<Real>& aBoundingBox);

        void reset();

        void build();

        void find(std::vector<Integer>& boundingBoxIds, CGeoBoundingBox<Real>& aBoundingBox, const Real aTolerance = 0.0);

        void addGeometricObject(const Integer aBoundingBoxId, CGeoBoundingBox<Real>& aBoundingBox);

        void removeGeometricObject(const Integer aBoundingBoxId);
    };
}
}

#include "GeoAdtree_Imp.hpp"
