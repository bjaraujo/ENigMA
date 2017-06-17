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
        CFemFace<Real>::CFemFace()
        {

        }

        template <typename Real>
        CFemFace<Real>::~CFemFace()
        {

        }

        template <typename Real>
        Real& CFemFace<Real>::thickness()
        {
        
            return m_thickness;
            
        }
        
    }        
}

