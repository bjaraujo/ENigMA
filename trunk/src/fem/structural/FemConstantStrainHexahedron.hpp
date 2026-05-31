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
#include "FemHexahedron.hpp"

namespace ENigMA
{
    namespace fem
    {
        namespace structural
        {
            template <typename Real, Integer NbNodes, Integer Dof, Integer Order>
            class CFemConstantStrainHexahedron : public CFemStructuralElement<Real>, public CFemHexahedron<Real, 8, 3, 1>
            {
            };

            template <typename Real>
            class CFemConstantStrainHexahedron<Real, 8, 3, 1> : public CFemStructuralElement<Real>, public CFemHexahedron<Real, 8, 3, 1>
            {
            protected:
                void setTransientTerm();
                void setDiffusionTerm();
                void setConvectiveTerm();

            public:
                CFemConstantStrainHexahedron();
                ~CFemConstantStrainHexahedron();

                void update();
                void reCalcD();
            };
        }
    }
}

#include "FemConstantStrainHexahedron_Imp.hpp"
