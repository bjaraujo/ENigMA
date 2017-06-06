// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// <Copyright> Copyright (c) 2012, All Rights Reserved </Copyright>
// *****************************************************************************

#pragma once

#include "../../src/fem/FemTriangle.hpp"
#include "FemStructuralElement.hpp"

namespace ENigMA
{

    namespace fem
    {

        namespace structural
        {

            template <typename Real, Integer NbNodes, Integer Dof, Integer Order>
            class CFemConstantStrainTriangle : public CFemStructuralElement<Real>, public CFemTriangle<Real, 3, 1, 1>
            {
            };

            template <typename Real>
            class CFemConstantStrainTriangle<Real, 3, 2, 1> : public CFemStructuralElement<Real>, public CFemTriangle<Real, 3, 2, 1>
            {
            protected:

                void setTransientTerm();
                void setDiffusionTerm();
                void setConvectiveTerm();

            public:

                CFemConstantStrainTriangle();
                ~CFemConstantStrainTriangle();

                void update();
                void reCalcD();
            
            };

        }

    }

}

#include "FemConstantStrainTriangle_Imp.hpp"


