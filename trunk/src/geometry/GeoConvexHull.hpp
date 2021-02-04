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

#include "GeoCoordinate.hpp"
#include "GeoVector.hpp"

namespace ENigMA {
namespace geometry {
    template <typename Real>
    struct SVertex {
        CGeoCoordinate<Real> v2D;
        CGeoCoordinate<Real> v3D;

        bool operator<(const SVertex& p) const
        {
            return this->v2D.x() < p.v2D.x() || (this->v2D.x() == p.v2D.x() && this->v2D.y() < p.v2D.y());
        }
    };

    // Class to compute planar convex hull of a series of vertices

    template <typename Real>
    class CGeoConvexHull {
    private:
        std::vector<CGeoCoordinate<Real>> m_vertices;

    public:
        explicit CGeoConvexHull(const std::vector<CGeoCoordinate<Real>>& sVertices);
        CGeoConvexHull();
        virtual ~CGeoConvexHull();

        void set(const std::vector<CGeoCoordinate<Real>>& sVertices);

        Integer nbVertices() const;

        void reset();

        void addVertex(CGeoCoordinate<Real>& aVertex);

        CGeoCoordinate<Real>& vertex(const Integer aVertexIndex);
    };
}
}

#include "GeoConvexHull_Imp.hpp"
