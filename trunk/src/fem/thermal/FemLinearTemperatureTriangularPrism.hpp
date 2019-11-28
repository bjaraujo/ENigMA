// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "FemThermalElement.hpp"
#include "FemTriangularPrism.hpp"

namespace ENigMA {
namespace fem {
    namespace thermal {
        template <typename Real, Integer NbNodes, Integer Dof, Integer Order>
        class CFemLinearTemperatureTriangularPrism : public CFemThermalElement<Real>, public CFemTriangularPrism<Real, NbNodes, Dof, Order> {
        };

        template <typename Real>
        class CFemLinearTemperatureTriangularPrism<Real, 6, 1, 1> : public CFemThermalElement<Real>, public CFemTriangularPrism<Real, 6, 1, 1> {
        public:
            CFemLinearTemperatureTriangularPrism();
            ~CFemLinearTemperatureTriangularPrism();

            void setConvectionOnFace(const Integer anFaceIndex, const Real h, const Real Tinf);
            void setConvectionOnFace(const Integer anFaceIndex, const Real e, const Real teta, const Real Tinf);
        };
    }
}
}

#include "FemLinearTemperatureTriangularPrism_Imp.hpp"
