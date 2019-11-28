// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "FemEdge.hpp"
#include "FemElement.hpp"
#include "GeoLine.hpp"

using namespace ENigMA::geometry;

namespace ENigMA {
namespace fem {
    /*
           Edge/line element.
           0                 1
           *-----------------*-----> xi
        */

    template <typename Real, Integer NbNodes, Integer Dof, Integer Order>
    class CFemBeam : public CFemElement<Real>, public CFemEdge<Real>, public CGeoLine<Real> {
    };

    template <typename Real>
    class CFemBeam<Real, 2, 1, 1> : public CFemElement<Real>, public CFemEdge<Real>, public CGeoLine<Real> {
    protected:
        void setTransientTerm();
        void setDiffusionTerm();
        void setConvectiveTerm();

    public:
        CFemBeam();
        virtual ~CFemBeam();

        void rebuild();
        void update();

        void setSourceOnNode(const Integer aNodeIndex, const Real aValue);
    };
}
}

#include "FemBeam_Imp.hpp"
