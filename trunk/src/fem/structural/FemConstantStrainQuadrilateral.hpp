// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "FemStructuralElement.hpp"
#include "FemQuadrilateral.hpp"

namespace ENigMA
{
    namespace fem
    {
        namespace structural
        {
            template <typename Real, Integer NbNodes, Integer Dof, Integer Order>
            class CFemConstantStrainQuadrilateral : public CFemStructuralElement<Real>, public CFemQuadrilateral<Real, 4, 2, 1>
            {
            };

            template <typename Real>
            class CFemConstantStrainQuadrilateral<Real, 4, 2, 1> : public CFemStructuralElement<Real>, public CFemQuadrilateral<Real, 4, 2, 1>
            {
            protected:
                void setTransientTerm();
                void setDiffusionTerm();
                void setConvectiveTerm();

            public:
                CFemConstantStrainQuadrilateral();
                ~CFemConstantStrainQuadrilateral();

                void update();
                void reCalcD();
            };
        }
    }
}

#include "FemConstantStrainQuadrilateral_Imp.hpp"
