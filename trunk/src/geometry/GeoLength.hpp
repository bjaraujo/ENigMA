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

namespace ENigMA
{
    namespace geometry
    {
        template <typename Real>
        class CGeoLength
        {
        protected:
            Real m_length;
            CGeoBoundingBox<Real> m_boundingBox;

            bool m_bLength;
            bool m_bBoundingBox;

        public:
            CGeoLength();
            virtual ~CGeoLength();

            virtual void calculateLength(bool bReCalculate = false) = 0;
            virtual void calculateBoundingBox(bool bReCalculate = false) = 0;

            Real length();
            CGeoBoundingBox<Real>& boundingBox();
        };
    }
}

#include "GeoLength_Imp.hpp"
