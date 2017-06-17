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

#include "GeoVector.hpp"
#include "GeoBoundingBox.hpp"
#include "GeoHashGrid.hpp"
#include "SphConvex.hpp"
#include "SphGaussian.hpp"
#include "SphSpiky.hpp"
#include "SphCubicSpline.hpp"
#include "SphQuintic.hpp"

using namespace ENigMA::geometry;

namespace ENigMA
{

    namespace sph
    {

        template <typename Real>
        struct CSphParticle
        {

            CGeoVector<Real> velocity;
            Real conductivity;
            Real mass;
            Real density;

        };

        template <typename Real>
        class CSphParticles
        {
        private:
            CSphKernel<Real>* m_kernel;

            std::vector<CSphParticle<Real> > m_particles;
            CGeoBoundingBox<Real> m_boundary;

            bool m_bCyclic;

            Real m_h;
            Real m_dt;

            void buildHashGrid(CPdeField<Real>& aField, CGeoHashGrid<Real>& aHashGrid);

            void advectParticles(CPdeField<Real>& aField);
            void addDiffusion(CPdeField<Real>& aField, CGeoHashGrid<Real>& aHashGrid);

        public:
            CSphParticles(CSphKernel<Real>& kernel);
            ~CSphParticles();

            void setBoundary(const CGeoBoundingBox<Real>& aBoundary);
            void setInitialVelocity(CPdeField<Real>& aField, const CGeoVector<Real>& aVelocity);

            void init(CPdeField<Real>& aField, Real mass, Real rho, Real diff, Real h, Real dt, bool bCyclic = false);
            void solve(CPdeField<Real>& aField);

        };

    }

}

#include "SphParticles_Imp.hpp"

