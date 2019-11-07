// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

using namespace ENigMA::geometry;

namespace ENigMA {

namespace fem {

    namespace thermal {

        template <typename Real>
        CFemLinearTemperatureTriangularPrism<Real, 6, 1, 1>::CFemLinearTemperatureTriangularPrism()
        {
        }

        template <typename Real>
        CFemLinearTemperatureTriangularPrism<Real, 6, 1, 1>::~CFemLinearTemperatureTriangularPrism()
        {
        }

        template <typename Real>
        void CFemLinearTemperatureTriangularPrism<Real, 6, 1, 1>::setConvectionOnFace(const Integer aFaceIndex, const Real h, const Real Tinf)
        {

            // TODO:
        }

        template <typename Real>
        void CFemLinearTemperatureTriangularPrism<Real, 6, 1, 1>::setConvectionOnFace(const Integer aFaceIndex, const Real e, const Real teta, const Real Tinf)
        {

            // Stefan-boltzmann constant
            Real sigma = 5.6704E-8;

            CFemThermalElement<Real>::setSourceOnEdge(aFaceIndex, sigma * e * teta, Tinf);
        }
    }
}
}
