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
    namespace thermal {
        template <typename Real>
        CFemThermalElement<Real>::CFemThermalElement()
            : m_specificHeat(1.0)
            , m_density(1.0)
            , m_thermalConductivity(1.0)
        {
        }

        template <typename Real>
        CFemThermalElement<Real>::~CFemThermalElement()
        {
        }

        template <typename Real>
        Real CFemThermalElement<Real>::specificHeat()
        {
            return m_specificHeat;
        }

        template <typename Real>
        Real CFemThermalElement<Real>::density()
        {
            return m_density;
        }

        template <typename Real>
        Real CFemThermalElement<Real>::thermalConductivity()
        {
            return m_thermalConductivity;
        }
    }
}
}
