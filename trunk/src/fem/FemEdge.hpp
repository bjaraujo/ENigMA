// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "CmnTypes.hpp"

namespace ENigMA {
namespace fem {
    template <typename Real>
    class CFemEdge {
    protected:
        Real m_sectionArea;
        Real m_perimeter;

    public:
        CFemEdge();
        virtual ~CFemEdge();

        Real sectionArea() const;
        Real perimeter() const;

        virtual void setSourceOnNode(const Integer aNodeIndex, const Real aValue) = 0;
    };
}
}

#include "FemEdge_Imp.hpp"
