// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "../../src/fem/FemTriangle.hpp"
#include "FemFlowElement.hpp"

namespace ENigMA
{

    namespace fem
    {

        namespace flow
        {

            template <typename Real, Integer NbNodes, Integer Dof, Integer Order>
            class CFemFlowTriangle : public CFemFlowElement<Real>, public CFemTriangle<Real, NbNodes, Dof, Order>
            {
            };

            template <typename Real>
            class CFemFlowTriangle<Real, 3, 1, 1> : public CFemFlowElement<Real>, public CFemTriangle<Real, 3, 1, 1>
            {
            private:

                void setDiffusionTerm();
                void setConvectiveTerm(Real u[3], Real v[3]);
                void setGradientTerm(const Integer aComponent);

            public:
                
                CFemFlowTriangle();
                ~CFemFlowTriangle();

                void calculateTransientTerm();
                void calculateDiffusiveTerm();
                void calculateConvectiveTerm(Real u[3], Real v[3]);
                void calculateGradientTerm(const Integer aComponent);

            };

        }

    }

}

#include "FemFlowTriangle_Imp.hpp"


