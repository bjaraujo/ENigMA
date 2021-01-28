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

namespace ENigMA {
namespace geometry {
    template <class T, typename Real>
    class CGeoContainer {
    protected:
        typedef std::map<Integer, Integer> mapGeometricObject;
        mapGeometricObject m_geometricObjectIds;

        std::vector<T> m_geometricObjects;

    public:
        CGeoContainer();
        virtual ~CGeoContainer();

        virtual void reset();

        virtual void build() = 0;

        virtual void find(std::vector<Integer>& sGeomtericObjectIds, const T& aGeometricObject, const Real aTolerance = 0.0) = 0;

        virtual void addGeometricObject(const Integer aGeomtericObjectId, const T& aGeometricObject);

        virtual void removeGeometricObject(const Integer aGeomtericObjectId);
        virtual void removeGeometricObject(const Integer aGeomtericObjectId, const T& aGeometricObject);
    };
}
}

#include "GeoContainer_Imp.hpp"
