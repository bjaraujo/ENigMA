// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// <Copyright> Copyright (c) 2012, All Rights Reserved </Copyright>
// *****************************************************************************

#pragma once

#include "GeoCoordinate.hpp"
#include "GeoCentroid.hpp"
#include "GeoBoundingBox.hpp"

namespace ENigMA
{

    namespace geometry
    {

        template <typename Real>
        class CGeoVolume : public CGeoCentroid<Real>
        {
        protected:

            Real m_surfaceArea;
            Real m_volume;
            CGeoBoundingBox<Real> m_boundingBox;

            bool m_bSurfaceArea;
            bool m_bVolume;
            bool m_bBoundingBox;

        public:

            CGeoVolume();
            ~CGeoVolume();

            virtual void calculateSurfaceArea(bool bReCalculate = false) = 0;
            virtual void calculateVolume(bool bReCalculate = false) = 0;
            virtual void calculateBoundingBox(bool bReCalculate = false) = 0;

            Real& surfaceArea();
            Real& volume();
            CGeoBoundingBox<Real>& boundingBox();

        };

    }

}

#include "GeoVolume_Imp.hpp"
