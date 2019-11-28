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
        CFemConstantStrainTriangle<Real, 3, 2, 1>::CFemConstantStrainTriangle()
        {
            this->source.resize(3);
            this->source.setZero();

            this->m_thickness = 1.0;
        }

        template <typename Real>
        CFemConstantStrainTriangle<Real, 3, 2, 1>::~CFemConstantStrainTriangle()
        {
        }

        template <typename Real>
        void CFemConstantStrainTriangle<Real, 3, 2, 1>::setTransientTerm()
        {
        }

        template <typename Real>
        void CFemConstantStrainTriangle<Real, 3, 2, 1>::setDiffusionTerm()
        {
            this->laplacian.resize(6, 6);

            Eigen::Matrix<Real, 3, 6> B;

            CFemTriangle<Real, 3, 2, 1>::calculateB(B);

            if (CGeoTriangle<Real>::area() > 0.0)
                B *= 1.0 / (2.0 * CGeoTriangle<Real>::area());

            Eigen::Matrix<Real, 3, 3> D;

            Real v = CFemStructuralElement<Real>::coeffPoisson();

            D << 1.0, v, 0.0,
                v, 1.0, 0.0,
                0.0, 0.0, (1.0 - v) * 0.5;

            Real E = CFemStructuralElement<Real>::elasticModulus();

            D *= E / (1 - v * v);

            CFemElement<Real>::laplacian = B.transpose() * D * B;
            CFemElement<Real>::laplacian *= CGeoTriangle<Real>::area() * this->m_thickness;
        }

        template <typename Real>
        void CFemConstantStrainTriangle<Real, 3, 2, 1>::setConvectiveTerm()
        {
        }

        template <typename Real>
        void CFemConstantStrainTriangle<Real, 3, 2, 1>::update()
        {
            // Calculate geometrical properties
            CFemTriangle<Real, 3, 2, 1>::rebuild();

            // Diffusion term
            setDiffusionTerm();

            // Transient term
            if (CFemElement<Real>::transient())
                setTransientTerm();
        }

        template <typename Real>
        void CFemConstantStrainTriangle<Real, 3, 2, 1>::reCalcD()
        {
        }
    }
}
}
