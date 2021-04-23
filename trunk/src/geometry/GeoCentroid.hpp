// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "GeoCoordinate.hpp"

namespace ENigMA
{
    namespace geometry
    {
        template <typename Real>
        class CGeoCentroid
        {
        protected:
            CGeoCoordinate<Real> m_centroid;
            bool m_bCentroid;

        public:
            CGeoCentroid();
            virtual ~CGeoCentroid();

            virtual void calculateCentroid(bool bReCalculate = false) = 0;

            CGeoCoordinate<Real> centroid() const;
        };
    }
}

#include "GeoCentroid_Imp.hpp"
