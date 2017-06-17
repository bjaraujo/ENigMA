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

        template <typename Real>
        CFemElement<Real>::CFemElement()
        {

            m_dt = 1.0;
            m_diffusionCoefficient = 1.0;
            m_convectionCoefficient = 1.0;

            m_transient = false;

        }

        template <typename Real>
        CFemElement<Real>::~CFemElement()
        {

        }

        template <typename Real>
        Real& CFemElement<Real>::dt()
        {
        
            return m_dt;
            
        }

        template <typename Real>
        Real& CFemElement<Real>::diffusionCoefficient()
        {
        
            return m_diffusionCoefficient;
            
        }

        template <typename Real>
        Real& CFemElement<Real>::convectionCoefficient()
        {
        
            return m_convectionCoefficient;
            
        }

        template <typename Real>
        bool& CFemElement<Real>::transient()
        {
        
            return m_transient;
            
        }

    }        
}

