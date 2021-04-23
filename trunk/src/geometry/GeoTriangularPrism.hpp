// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "GeoVertexList.hpp"
#include "GeoVolume.hpp"

namespace ENigMA
{
    namespace geometry
    {
        template <typename Real>
        class CGeoTriangularPrism : public CGeoVolume<Real>, public CGeoVertexList<Real>
        {
        public:
            CGeoTriangularPrism();
            virtual ~CGeoTriangularPrism();

            void reset();

            inline void calculateCentroid(bool bReCalculate = false);
            inline void calculateSurfaceArea(bool bReCalculate = false);
            inline void calculateVolume(bool bReCalculate = false);
            inline void calculateBoundingBox(bool bReCalculate = false);
        };
    }
}

#include "GeoTriangularPrism_Imp.hpp"
