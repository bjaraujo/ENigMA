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

#include "GeoVolume.hpp"
#include "GeoVertexList.hpp"

namespace ENigMA
{

    namespace geometry
    {

        template <typename Real>
        class CGeoTetrahedron : public CGeoVolume<Real>, public CGeoVertexList<Real>
        {
        public:
            CGeoTetrahedron();
            ~CGeoTetrahedron();

            void reset();

            inline void calculateCentroid(bool bReCalculate = false);
            inline void calculateSurfaceArea(bool bReCalculate = false);
            inline void calculateVolume(bool bReCalculate = false);
            inline void calculateBoundingBox(bool bReCalculate = false);

            inline bool contains(const CGeoCoordinate<Real>& aPoint, const Real aTolerance = 0.0);

            void invert();

        };

    }
}

#include "GeoTetrahedron_Imp.hpp"

