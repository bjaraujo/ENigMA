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

        namespace structural
        {

            template <typename Real>
            CFemStructuralElement<Real>::CFemStructuralElement()
            {

                m_coeffPoisson = 1.0;
                m_elasticModulus = 1.0;

            }

            template <typename Real>
            CFemStructuralElement<Real>::~CFemStructuralElement()
            {

            }
            
            template <typename Real>
            Real& CFemStructuralElement<Real>::coeffPoisson()
            {
            
                return m_coeffPoisson;
                
            }

            template <typename Real>
            Real& CFemStructuralElement<Real>::elasticModulus()
            {
            
                return m_elasticModulus;
                
            }

        }

    }        
}

