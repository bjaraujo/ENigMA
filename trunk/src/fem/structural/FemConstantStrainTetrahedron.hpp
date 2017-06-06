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

#include "../../src/fem/FemTetrahedron.hpp"
#include "FemStructuralElement.hpp"

namespace ENigMA
{

    namespace fem
    {

        namespace structural
        {

            template <typename Real, Integer NbNodes, Integer Dof, Integer Order>
            class CFemConstantStrainTetrahedron : public CFemStructuralElement<Real>, public CFemTetrahedron<Real, 4, 1, 1>
            {
            };

            template <typename Real>
            class CFemConstantStrainTetrahedron<Real, 4, 3, 1> : public CFemStructuralElement<Real>, public CFemTetrahedron<Real, 4, 3, 1>
            {
            protected:

                void setTransientTerm();
                void setDiffusionTerm();
                void setConvectiveTerm();

            public:

                CFemConstantStrainTetrahedron();
                ~CFemConstantStrainTetrahedron();

                void update();
                void reCalcD();
            
            };

        }

    }

}

#include "FemConstantStrainTetrahedron_Imp.hpp"


