// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "IntGaussIntegration.hpp"

namespace ENigMA
{

    namespace integration
    {

        template <typename Real>
        class CIntBeam : public CIntGaussIntegration<Real>
        {
        protected:
            std::vector<Real> m_xi;
            std::vector<Real> m_wxi;

            void setGaussPoints();

        public:
            CIntBeam();
            virtual ~CIntBeam();
        };
    }
}

#include "IntBeam_Imp.hpp"
