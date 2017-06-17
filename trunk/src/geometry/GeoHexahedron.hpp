// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include <vector>

#include "GeoVolume.hpp"
#include "GeoVertexList.hpp"
#include "GeoTetrahedron.hpp"

namespace ENigMA
{

    namespace geometry
    {

        template <typename Real>
        class CGeoHexahedron : public CGeoVolume<Real>, public CGeoVertexList<Real>
        {
        public:
            CGeoHexahedron();
            ~CGeoHexahedron();

            void reset();

            inline void calculateCentroid(bool bReCalculate = false);
            inline void calculateSurfaceArea(bool bReCalculate = false);
            inline void calculateVolume(bool bReCalculate = false);
            inline void calculateBoundingBox(bool bReCalculate = false);

            void decimate(std::vector<CGeoTetrahedron<Real> >& sTetrahedrons);

        };

    }
}

#include "GeoHexahedron_Imp.hpp"

