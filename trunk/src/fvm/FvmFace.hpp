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

#include "FvmNode.hpp"
#include "GeoPolygon.hpp"
#include "MshFace.hpp"

namespace ENigMA {

namespace fvm {

    enum EBoundaryType {
        BT_NONE = 0,
        BT_INTERIOR_FACE,
        BT_INTERFACE,
        BT_WALL_NO_SLIP,
        BT_INLET_FLOW,
        BT_INLET_PRESSURE,
        BT_OUTLET
    };

    class CFvmBoundaryType {
    private:
        EBoundaryType m_boundaryType;

    public:
        CFvmBoundaryType()
            : m_boundaryType(BT_NONE){};
        CFvmBoundaryType(const EBoundaryType aBoundaryType)
            : m_boundaryType(aBoundaryType){};

        CFvmBoundaryType operator=(const EBoundaryType& aBoundaryType)
        {

            this->m_boundaryType = aBoundaryType;
            return *this;
        }

        bool operator==(const EBoundaryType& aBoundaryType) const
        {

            return this->m_boundaryType == aBoundaryType;
        }

        bool operator!=(const EBoundaryType& aBoundaryType) const
        {

            return this->m_boundaryType != aBoundaryType;
        }
    };

    template <typename Real>
    class CFvmFace : public ENigMA::mesh::CMshFace<Real>, public ENigMA::geometry::CGeoPolygon<Real> {
    private:
        Integer m_controlVolumeId;
        Integer m_neighborId;

        CFvmBoundaryType m_boundaryType;

    public:
        CFvmFace();
        explicit CFvmFace(ENigMA::geometry::CGeoPolygon<Real>& aPolygon);
        virtual ~CFvmFace();

        void set(ENigMA::geometry::CGeoPolygon<Real>& aPolygon);

        void addNode(const CFvmNode<Real>& aNode);

        void setControlVolumeId(const Integer aControlVolumeId);
        Integer controlVolumeId();
        Integer controlVolumeId(const Integer aControlVolumeId);

        void setNeighborId(const Integer aNeighborId);
        Integer neighborId(const Integer aControlVolumeId);

        void close();

        void reset();

        CFvmBoundaryType& boundaryType();
        void setBoundaryType(CFvmBoundaryType aBoundaryType);
    };
}
}

#include "FvmFace_Imp.hpp"
