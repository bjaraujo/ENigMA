// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

namespace ENigMA {
namespace geometry {
    template <class T, typename Real>
    CGeoContainer<T, Real>::CGeoContainer()
    {
    }

    template <class T, typename Real>
    CGeoContainer<T, Real>::~CGeoContainer()
    {
    }

    template <class T, typename Real>
    void CGeoContainer<T, Real>::reset()
    {
        m_geometricObjectIds.clear();
        m_geometricObjects.clear();
    }

    template <class T, typename Real>
    void CGeoContainer<T, Real>::addGeometricObject(const Integer aGeomtericObjectId, const T& aGeometricObject)
    {
        m_geometricObjectIds[static_cast<Integer>(m_geometricObjects.size())] = aGeomtericObjectId;
        m_geometricObjects.push_back(aGeometricObject);
    }

    template <class T, typename Real>
    void CGeoContainer<T, Real>::removeGeometricObject(const Integer aGeomtericObjectId)
    {
    }

    template <class T, typename Real>
    void CGeoContainer<T, Real>::removeGeometricObject(const Integer aGeomtericObjectId, const T& aGeometricObject)
    {
    }
}
}
