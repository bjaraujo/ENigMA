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
#include "FemTriangularPrism.hpp"

namespace ENigMA
{
    namespace fem
    {
        namespace structural
        {
            template <typename Real, Integer NbNodes, Integer Dof, Integer Order>
            class CFemConstantStrainTriangularPrism : public CFemStructuralElement<Real>, public CFemTriangularPrism<Real, 6, 3, 1>
            {
            };

            template <typename Real>
            class CFemConstantStrainTriangularPrism<Real, 6, 3, 1> : public CFemStructuralElement<Real>, public CFemTriangularPrism<Real, 6, 3, 1>
            {
            protected:
                void setTransientTerm();
                void setDiffusionTerm();
                void setConvectiveTerm();

            public:
                CFemConstantStrainTriangularPrism();
                ~CFemConstantStrainTriangularPrism();

                void update();
                void reCalcD();
            };
        }
    }
}

#include "FemConstantStrainTriangularPrism_Imp.hpp"
