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

    namespace integration
    {

        template <typename Real>
        class CIntGaussIntegration
        {
        protected:
            Integer m_integPoints;
            virtual void setGaussPoints() = 0;
        };
    }
}
