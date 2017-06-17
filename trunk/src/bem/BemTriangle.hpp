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
#include "GeoTriangle.hpp"
#include "IntTriangle.hpp"

using namespace ENigMA::geometry;
using namespace ENigMA::integration;

namespace ENigMA
{

    namespace bem
    {

        template <typename Real>
        class CBemTriangle : public CIntTriangle<Real>, public CGeoTriangle<Real>
        {
        private:

        public:

            bool getIntegrationPoints(std::vector<CGeoCoordinate<Real> >& sPoints, std::vector<Real>& sWeights);
            void laplacianCoeff(const Integer i, const Integer j, CBemTriangle<Real>& aBemTriangle, Real& h, Real& g);

            CBemTriangle();
            ~CBemTriangle();
        
        };

    }

}

#include "BemTriangle_Imp.hpp"


