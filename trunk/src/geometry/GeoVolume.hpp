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
#include "GeoCentroid.hpp"
#include "GeoCoordinate.hpp"

namespace ENigMA {
namespace geometry {
    template <typename Real>
    class CGeoVolume : public CGeoCentroid<Real> {
    protected:
        Real m_surfaceArea;
        Real m_volume;
        CGeoBoundingBox<Real> m_boundingBox;

        bool m_bSurfaceArea;
        bool m_bVolume;
        bool m_bBoundingBox;

    public:
        CGeoVolume();
        virtual ~CGeoVolume();

        virtual void calculateSurfaceArea(bool bReCalculate = false) = 0;
        virtual void calculateVolume(bool bReCalculate = false) = 0;
        virtual void calculateBoundingBox(bool bReCalculate = false) = 0;

        Real surfaceArea() const;        
        Real volume() const;

        CGeoBoundingBox<Real>& boundingBox();
    };
}
}

#include "GeoVolume_Imp.hpp"
