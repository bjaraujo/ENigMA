// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "FemTetrahedron.hpp"
#include "FemFlowElement.hpp"

namespace ENigMA
{

    namespace fem
    {

        namespace flow
        {

            template <typename Real, Integer NbNodes, Integer Dof, Integer Order>
            class CFemFlowTetrahedron : public CFemFlowElement<Real>, public CFemTetrahedron<Real, NbNodes, Dof, Order>
            {
            };

            template <typename Real>
            class CFemFlowTetrahedron<Real, 4, 1, 1> : public CFemFlowElement<Real>, public CFemTetrahedron<Real, 4, 1, 1>
            {
            private:

                void setDiffusionTerm();
                void setConvectiveTerm(Real u[3], Real v[3], Real w[3]);
                void setGradientTerm(const Integer aComponent);

            public:
                
                CFemFlowTetrahedron();
                ~CFemFlowTetrahedron();

                void calculateTransientTerm();
                void calculateDiffusiveTerm();
                void calculateConvectiveTerm(Real u[3], Real v[3], Real w[3]);
                void calculateGradientTerm(const Integer aComponent);

            };

        }

    }

}

#include "FemFlowTetrahedron_Imp.hpp"


