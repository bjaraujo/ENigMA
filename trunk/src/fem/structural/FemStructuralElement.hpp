// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

namespace ENigMA {

namespace fem {

    namespace structural {

        template <typename Real>
        class CFemStructuralElement {
        private:
            Real m_coeffPoisson;
            Real m_elasticModulus;

        public:
            CFemStructuralElement();
            ~CFemStructuralElement();

            void setCoeffPoisson(const Real aValue);
            Real coeffPoisson() const;

            void setElasticModulus(const Real aValue);
            Real elasticModulus() const;
        };
    }
}
}

#include "FemStructuralElement_Imp.hpp"
