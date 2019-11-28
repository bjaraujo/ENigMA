// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "GeoTetrahedron.hpp"

namespace ENigMA {
namespace mesh {
    template <typename Real>
    class CMshTetrahedron : public ENigMA::geometry::CGeoTetrahedron<Real> {
    private:
        Real m_quality;

    public:
        CMshTetrahedron();
        virtual ~CMshTetrahedron();

        void calculateQuality();
        Real quality();
    };
}
}

#include "MshTetrahedron_Imp.hpp"
