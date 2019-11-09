// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "FemQuadrilateral.hpp"
#include "FemThermalElement.hpp"

namespace ENigMA {

namespace fem {

    namespace thermal {

        template <typename Real, Integer NbNodes, Integer Dof, Integer Order>
        class CFemLinearTemperatureQuadrilateral : public CFemThermalElement<Real>, public CFemQuadrilateral<Real, NbNodes, Dof, Order> {
        };

        template <typename Real>
        class CFemLinearTemperatureQuadrilateral<Real, 4, 1, 1> : public CFemThermalElement<Real>, public CFemQuadrilateral<Real, 4, 1, 1> {
        public:
            CFemLinearTemperatureQuadrilateral();
            ~CFemLinearTemperatureQuadrilateral();

            void setConvectionOnEdge(const Integer anEdgeIndex, const Real h, const Real Tinf);
            void setConvectionOnEdge(const Integer anEdgeIndex, const Real e, const Real teta, const Real Tinf);
        };
    }
}
}

#include "FemLinearTemperatureQuadrilateral_Imp.hpp"
