// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

using namespace ENigMA::geometry;

namespace ENigMA {

namespace fvm {

    template <typename Real>
    CFvmFace<Real>::CFvmFace()
        : m_boundaryType(BT_NONE)
        , m_controlVolumeId(0)
        , m_neighborId(0)
    {
    }

    template <typename Real>
    CFvmFace<Real>::CFvmFace(CGeoPolygon<Real>& aPolygon)
        : m_boundaryType(BT_NONE)
        , m_controlVolumeId(0)
        , m_neighborId(0)
    {

        this->set(aPolygon);
    }

    template <typename Real>
    CFvmFace<Real>::~CFvmFace()
    {
    }

    template <typename Real>
    void CFvmFace<Real>::set(CGeoPolygon<Real>& aPolygon)
    {

        CGeoPolygon<Real>::polyline().reset();

        for (Integer i = 0; i < aPolygon.polyline().nbVertices(); ++i) {

            CGeoPolygon<Real>::polyline().addVertex(aPolygon.polyline().vertex(i));
        }

        CGeoPolygon<Real>::polyline().close();
    }

    template <typename Real>
    void CFvmFace<Real>::addNode(const CFvmNode<Real>& aNode)
    {

        CGeoCoordinate<Real> aVertex(aNode.x(), aNode.y(), aNode.z());

        CGeoPolygon<Real>::polyline().addVertex(aVertex);
    }

    template <typename Real>
    void CFvmFace<Real>::setControlVolumeId(const Integer aControlVolumeId)
    {

        m_controlVolumeId = aControlVolumeId;
    }

    template <typename Real>
    Integer CFvmFace<Real>::controlVolumeId()
    {

        return m_controlVolumeId;
    }

    template <typename Real>
    Integer CFvmFace<Real>::controlVolumeId(const Integer aControlVolumeId)
    {

        if (m_controlVolumeId == aControlVolumeId)
            return m_controlVolumeId;
        else
            return m_neighborId;
    }

    template <typename Real>
    void CFvmFace<Real>::setNeighborId(const Integer aNeighborId)
    {

        m_neighborId = aNeighborId;
    }

    template <typename Real>
    Integer CFvmFace<Real>::neighborId(const Integer aControlVolumeId)
    {

        if (m_controlVolumeId == aControlVolumeId)
            return m_neighborId;
        else
            return m_controlVolumeId;
    }

    template <typename Real>
    void CFvmFace<Real>::close()
    {

        CGeoPolygon<Real>::polyline().close();
    }

    template <typename Real>
    void CFvmFace<Real>::reset()
    {

        CGeoPolygon<Real>::reset();
    }

    template <typename Real>
    CFvmBoundaryType& CFvmFace<Real>::boundaryType()
    {

        return m_boundaryType;
    }

    template <typename Real>
    void CFvmFace<Real>::setBoundaryType(CFvmBoundaryType aBoundaryType)
    {

        m_boundaryType = aBoundaryType;
    }
}
}
