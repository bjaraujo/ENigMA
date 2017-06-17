// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "AnaFunction.hpp"

namespace ENigMA
{

    namespace analytical
    {

        template <typename Real>
        class CAnaTemperature : public CAnaFunction<Real>
        {
        public:

            CAnaTemperature();
            ~CAnaTemperature();

            void steadyStateHeatConduction1D(Real x, Real& T);
            void steadyStateHeatConduction2D(Real x, Real y, Real& T);
            void steadyStateHeatConduction3D(Real x, Real y, Real z, Real& T);

            void steadyStateHeatConduction1D(Real x, Real Ta, Real Tb, Real h, Real Tinf, Real k, Real perimeter, Real sectionArea, Real length, Real& T);
            void steadyStateHeatConvectionRadiation1D(Real x, Real Tb, Real h, Real e, Real k, Real perimeter, Real sectionArea, Real& T);

            void transientHeatConduction1D(Real x, Real t, Real alpha, Real &T);
            void transientHeatConduction2D(Real x, Real y, Real t, Real alpha, Real& T);
            void transientHeatConduction3D(Real x, Real y, Real z, Real t, Real alpha, Real& T);

        };

    }

}

#include "AnaTemperature_Imp.hpp"

