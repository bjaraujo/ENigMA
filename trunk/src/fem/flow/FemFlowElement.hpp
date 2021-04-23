// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

namespace ENigMA
{
    namespace fem
    {
        namespace flow
        {
            template <typename Real>
            class CFemFlowElement
            {
            private:
                Real m_density;
                Real m_viscosity;

            public:
                CFemFlowElement();
                ~CFemFlowElement();

                void setDensity(const Real aValue);
                Real density() const;

                void setViscosity(const Real aValue);
                Real viscosity() const;
            };
        }
    }
}

#include "FemFlowElement_Imp.hpp"
