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

namespace ENigMA
{
    namespace fem
    {
        namespace structural
        {
            template <typename Real>
            CFemConstantStrainQuadrilateral<Real, 4, 2, 1>::CFemConstantStrainQuadrilateral()
            {
                this->source.resize(4);
                this->source.setZero();

                this->m_thickness = 1.0;
            }

            template <typename Real>
            CFemConstantStrainQuadrilateral<Real, 4, 2, 1>::~CFemConstantStrainQuadrilateral()
            {
            }

            template <typename Real>
            void CFemConstantStrainQuadrilateral<Real, 4, 2, 1>::setTransientTerm()
            {
            }

            template <typename Real>
            void CFemConstantStrainQuadrilateral<Real, 4, 2, 1>::setDiffusionTerm()
            {
                this->laplacian.resize(8, 8);

                Eigen::Matrix<Real, 3, 8> B;

                CFemQuadrilateral<Real, 4, 2, 1>::calculateB(B);

                if (CGeoQuadrilateral<Real>::area() > 0.0)
                    B *= 1.0 / (2.0 * CGeoQuadrilateral<Real>::area());

                Eigen::Matrix<Real, 3, 3> D;

                Real v = CFemStructuralElement<Real>::coeffPoisson();

                D << 1.0, v, 0.0,
                    v, 1.0, 0.0,
                    0.0, 0.0, (1.0 - v) * 0.5;

                Real E = CFemStructuralElement<Real>::elasticModulus();

                D *= E / (1 - v * v);

                CFemElement<Real>::laplacian = B.transpose() * D * B;
                CFemElement<Real>::laplacian *= CGeoQuadrilateral<Real>::area() * this->m_thickness;
            }

            template <typename Real>
            void CFemConstantStrainQuadrilateral<Real, 4, 2, 1>::setConvectiveTerm()
            {
            }

            template <typename Real>
            void CFemConstantStrainQuadrilateral<Real, 4, 2, 1>::update()
            {
                // Calculate geometrical properties
                CFemQuadrilateral<Real, 4, 2, 1>::rebuild();

                // Diffusion term
                setDiffusionTerm();

                // Transient term
                if (CFemElement<Real>::transient())
                    setTransientTerm();
            }

            template <typename Real>
            void CFemConstantStrainQuadrilateral<Real, 4, 2, 1>::reCalcD()
            {
            }
        }
    }
}
