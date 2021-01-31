// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "GeoQuadrilateral.hpp"

namespace ENigMA {
namespace mesh {
    template <typename Real>
    class CMshQuadrilateral : public ENigMA::geometry::CGeoQuadrilateral<Real> {
    private:
        Real m_quality;

    public:
        CMshQuadrilateral();
        virtual ~CMshQuadrilateral();

        void calculateQuality();
        Real quality() const;
    };
}
}

#include "MshQuadrilateral_Imp.hpp"
