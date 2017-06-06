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

#include "GeoNormal.hpp"
#include "GeoCentroid.hpp"
#include "GeoBoundingBox.hpp"

namespace ENigMA
{

    namespace geometry
    {

        template <typename Real>
        class CGeoArea : public CGeoCentroid<Real>
        {
        protected:

            Real m_area;
            CGeoNormal<Real> m_normal;
            CGeoBoundingBox<Real> m_boundingBox;

            bool m_bNormal;
            bool m_bArea;
            bool m_bBoundingBox;

        public:

            CGeoArea();
            ~CGeoArea();

            virtual void calculateNormal(bool bReCalculate = false) = 0;
            virtual void calculateArea(bool bReCalculate = false) = 0;
            virtual void calculateBoundingBox(bool bReCalculate = false) = 0;

            CGeoNormal<Real>& normal();
            Real& area();
            CGeoBoundingBox<Real>& boundingBox();

        };

    }

}

#include "GeoArea_Imp.hpp"
