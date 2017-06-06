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

namespace ENigMA
{

    namespace fem
    {

        namespace flow
        {

            template <typename Real>
            CFemFlowElement<Real>::CFemFlowElement()
            {

                m_density = 1.0;
                m_viscosity = 1.0;

            }

            template <typename Real>
            CFemFlowElement<Real>::~CFemFlowElement()
            {

            }
            
            template <typename Real>
            Real& CFemFlowElement<Real>::density()
            {
            
                return m_density;
                
            }

            template <typename Real>
            Real& CFemFlowElement<Real>::viscosity()
            {
            
                return m_viscosity;
                
            }

        }

    }        
}

