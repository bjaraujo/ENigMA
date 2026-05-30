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
                this->laplacian.setZero();

                Real v = CFemStructuralElement<Real>::coeffPoisson();
                Real E = CFemStructuralElement<Real>::elasticModulus();

                Eigen::Matrix<Real, 3, 3> D;
                D << 1.0, v, 0.0,
                     v, 1.0, 0.0,
                     0.0, 0.0, (1.0 - v) * 0.5;
                D *= E / (1.0 - v * v);

                for (Integer p = 0; p < this->m_integPoints; p++)
                {
                    Real xi  = this->m_xi[p];
                    Real eta = this->m_eta[p];

                    Real dN1dxi  = -(1.0 - eta) * 0.25;
                    Real dN2dxi  =  (1.0 - eta) * 0.25;
                    Real dN3dxi  =  (1.0 + eta) * 0.25;
                    Real dN4dxi  = -(1.0 + eta) * 0.25;

                    Real dN1deta = -(1.0 - xi) * 0.25;
                    Real dN2deta = -(1.0 + xi) * 0.25;
                    Real dN3deta =  (1.0 + xi) * 0.25;
                    Real dN4deta =  (1.0 - xi) * 0.25;

                    Eigen::Matrix<Real, 2, 2> J;
                    J(0, 0) = dN1dxi * this->m_x1 + dN2dxi * this->m_x2 + dN3dxi * this->m_x3 + dN4dxi * this->m_x4;
                    J(0, 1) = dN1dxi * this->m_y1 + dN2dxi * this->m_y2 + dN3dxi * this->m_y3 + dN4dxi * this->m_y4;
                    J(1, 0) = dN1deta * this->m_x1 + dN2deta * this->m_x2 + dN3deta * this->m_x3 + dN4deta * this->m_x4;
                    J(1, 1) = dN1deta * this->m_y1 + dN2deta * this->m_y2 + dN3deta * this->m_y3 + dN4deta * this->m_y4;

                    Real detJ = J.determinant();

                    if (std::fabs(detJ) < 1E-15)
                        continue;

                    Eigen::Matrix<Real, 2, 2> Jinv = J.inverse();

                    Real b1x = Jinv(0,0) * dN1dxi + Jinv(0,1) * dN1deta;
                    Real b1y = Jinv(1,0) * dN1dxi + Jinv(1,1) * dN1deta;
                    Real b2x = Jinv(0,0) * dN2dxi + Jinv(0,1) * dN2deta;
                    Real b2y = Jinv(1,0) * dN2dxi + Jinv(1,1) * dN2deta;
                    Real b3x = Jinv(0,0) * dN3dxi + Jinv(0,1) * dN3deta;
                    Real b3y = Jinv(1,0) * dN3dxi + Jinv(1,1) * dN3deta;
                    Real b4x = Jinv(0,0) * dN4dxi + Jinv(0,1) * dN4deta;
                    Real b4y = Jinv(1,0) * dN4dxi + Jinv(1,1) * dN4deta;

                    Eigen::Matrix<Real, 3, 8> B;
                    B << b1x, 0.0, b2x, 0.0, b3x, 0.0, b4x, 0.0,
                         0.0, b1y, 0.0, b2y, 0.0, b3y, 0.0, b4y,
                         b1y, b1x, b2y, b2x, b3y, b3x, b4y, b4x;

                    this->laplacian += B.transpose() * D * B * std::fabs(detJ) * this->m_wxi[p] * this->m_weta[p];
                }

                this->laplacian *= this->m_thickness;
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

                // Integration points
                this->setGaussPoints();

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
