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
    class CFemFace {
    protected:
        Real m_thickness;

    public:
        CFemFace();
        virtual ~CFemFace();

        void setThickness(const Real aValue);
        Real thickness() const;

        virtual void setSourceOnNode(const Integer aNodeIndex, const Real aValue) = 0;
        virtual void setSourceOnEdge(const Integer anEdgeIndex, const Real aValue) = 0;
    };
}
}

#include "FemFace_Imp.hpp"
