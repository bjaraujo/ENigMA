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
#include "RTree.h"

namespace ENigMA
{
    namespace geometry
    {
        template <typename Real>
        class CGeoRtree : public CGeoContainer<CGeoBoundingBox<Real>, Real>
        {
        private:
            RTree<Integer, Real, 3, Real> m_tree;

        public:
            CGeoRtree();
            virtual ~CGeoRtree();

            void reset();

            void build();

            void find(std::vector<Integer>& boundingBoxIds, const CGeoBoundingBox<Real>& aBoundingBox, const Real aTolerance = 0.0) override;

            void addGeometricObject(const Integer aBoundingBoxId, const CGeoBoundingBox<Real>& aBoundingBox) override;
            void removeGeometricObject(const Integer aBoundingBoxId, const CGeoBoundingBox<Real>& aBoundingBox) override;
        };
    }
}

#include "GeoRtree_Imp.hpp"
