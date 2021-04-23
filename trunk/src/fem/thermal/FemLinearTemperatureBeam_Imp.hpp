// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

using namespace ENigMA::geometry;

namespace ENigMA
{
    namespace fem
    {
        namespace thermal
        {
            template <typename Real>
            CFemLinearTemperatureBeam<Real, 2, 1, 1>::CFemLinearTemperatureBeam()
            {
            }

            template <typename Real>
            CFemLinearTemperatureBeam<Real, 2, 1, 1>::~CFemLinearTemperatureBeam()
            {
            }

            template <typename Real>
            void CFemLinearTemperatureBeam<Real, 2, 1, 1>::setConvectionOnEdge(const Real h, const Real Tinf)
            {
                CFemElement<Real>::source(0) += h * Tinf * this->m_perimeter * this->m_length;
                CFemElement<Real>::source(1) += h * Tinf * this->m_perimeter * this->m_length;

                CFemElement<Real>::laplacian(0, 0) += h * this->m_perimeter * this->m_length;
                CFemElement<Real>::laplacian(1, 1) += h * this->m_perimeter * this->m_length;
            }

            template <typename Real>
            void CFemLinearTemperatureBeam<Real, 2, 1, 1>::setConvectionOnEdge(const Real e, const Real teta, const Real Tinf)
            {
                this->setConvectionOnEdge(e, teta, Tinf);
            }
        }
    }
}
