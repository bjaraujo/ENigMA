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

    template <typename Real>
    CFemEdge<Real>::CFemEdge()
    {
    }

    template <typename Real>
    CFemEdge<Real>::~CFemEdge()
    {
    }

    template <typename Real>
    Real& CFemEdge<Real>::sectionArea()
    {

        return m_sectionArea;
    }

    template <typename Real>
    Real& CFemEdge<Real>::perimeter()
    {

        return m_perimeter;
    }
}
}
