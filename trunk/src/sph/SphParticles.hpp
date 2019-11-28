// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include <vector>

#include "GeoBoundingBox.hpp"
#include "GeoHashGrid.hpp"
#include "GeoVector.hpp"
#include "SphConvex.hpp"
#include "SphCubicSpline.hpp"
#include "SphGaussian.hpp"
#include "SphQuintic.hpp"
#include "SphSpiky.hpp"

namespace ENigMA {
namespace sph {
    template <typename Real>
    struct CSphParticle {
        ENigMA::geometry::CGeoVector<Real> velocity;
        Real conductivity;
        Real mass;
        Real density;
    };

    template <typename Real>
    class CSphParticles {
    private:
        CSphKernel<Real>* m_kernel;

        std::vector<CSphParticle<Real>> m_particles;
        ENigMA::geometry::CGeoBoundingBox<Real> m_boundary;

        bool m_bCyclic;

        Real m_h;
        Real m_dt;

        void buildHashGrid(CPdeField<Real>& aField, ENigMA::geometry::CGeoHashGrid<Real>& aHashGrid);

        void advectParticles(CPdeField<Real>& aField);
        void addDiffusion(CPdeField<Real>& aField, ENigMA::geometry::CGeoHashGrid<Real>& aHashGrid);

    public:
        explicit CSphParticles(CSphKernel<Real>& kernel);
        virtual ~CSphParticles();

        void setBoundary(const ENigMA::geometry::CGeoBoundingBox<Real>& aBoundary);
        void setInitialVelocity(CPdeField<Real>& aField, const ENigMA::geometry::CGeoVector<Real>& aVelocity);

        void init(CPdeField<Real>& aField, Real mass, Real rho, Real diff, Real h, Real dt, bool bCyclic = false);
        void solve(CPdeField<Real>& aField);
    };
}
}

#include "SphParticles_Imp.hpp"
