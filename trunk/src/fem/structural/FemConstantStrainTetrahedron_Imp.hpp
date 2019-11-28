// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

using namespace ENigMA::geometry;

namespace ENigMA {
namespace fem {
    namespace structural {
        template <typename Real>
        CFemConstantStrainTetrahedron<Real, 4, 3, 1>::CFemConstantStrainTetrahedron()
        {
            this->source.resize(12);
            this->source.setZero();
        }

        template <typename Real>
        CFemConstantStrainTetrahedron<Real, 4, 3, 1>::~CFemConstantStrainTetrahedron()
        {
        }

        template <typename Real>
        void CFemConstantStrainTetrahedron<Real, 4, 3, 1>::setTransientTerm()
        {
        }

        template <typename Real>
        void CFemConstantStrainTetrahedron<Real, 4, 3, 1>::setDiffusionTerm()
        {
            this->laplacian.resize(12, 12);

            Eigen::Matrix<Real, 6, 12> B;

            CFemTetrahedron<Real, 4, 3, 1>::calculateB(B);

            Eigen::Matrix<Real, 6, 6> D;

            Real v = CFemStructuralElement<Real>::coeffPoisson();

            D << 1.0 - v, v, 0.0, 0.0, 0.0, 0.0,
                v, 1.0 - v, v, 0.0, 0.0, 0.0,
                v, v, 1.0 - v, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.5 - v, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.5 - v, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.5 - v;

            Real E = CFemStructuralElement<Real>::elasticModulus();

            D *= E / ((1.0 + v) * (1.0 - 2.0 * v));

            CFemElement<Real>::laplacian = B.transpose() * D * B;
            CFemElement<Real>::laplacian *= CGeoTetrahedron<Real>::volume();
        }

        template <typename Real>
        void CFemConstantStrainTetrahedron<Real, 4, 3, 1>::setConvectiveTerm()
        {
        }

        template <typename Real>
        void CFemConstantStrainTetrahedron<Real, 4, 3, 1>::update()
        {
            // Calculate geometrical properties
            CFemTetrahedron<Real, 4, 3, 1>::rebuild();

            // Diffusion term
            setDiffusionTerm();

            // Transient term
            if (CFemElement<Real>::transient())
                setTransientTerm();
        }

        template <typename Real>
        void CFemConstantStrainTetrahedron<Real, 4, 3, 1>::reCalcD()
        {
        }
    }
}
}
