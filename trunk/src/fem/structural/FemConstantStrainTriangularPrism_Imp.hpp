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
            CFemConstantStrainTriangularPrism<Real, 6, 3, 1>::CFemConstantStrainTriangularPrism()
            {
                this->source.resize(6);
                this->source.setZero();
                CIntGaussIntegration<Real>::m_integPoints = 6;
            }

            template <typename Real>
            CFemConstantStrainTriangularPrism<Real, 6, 3, 1>::~CFemConstantStrainTriangularPrism()
            {
            }

            template <typename Real>
            void CFemConstantStrainTriangularPrism<Real, 6, 3, 1>::setTransientTerm()
            {
            }

            template <typename Real>
            void CFemConstantStrainTriangularPrism<Real, 6, 3, 1>::setDiffusionTerm()
            {
                this->laplacian.resize(18, 18);
                this->laplacian.setZero();

                Real v = CFemStructuralElement<Real>::coeffPoisson();
                Real E = CFemStructuralElement<Real>::elasticModulus();

                Eigen::Matrix<Real, 6, 6> D;
                D << 1.0 - v,  v,        v,        0.0,       0.0,       0.0,
                     v,        1.0 - v,  v,        0.0,       0.0,       0.0,
                     v,        v,        1.0 - v,  0.0,       0.0,       0.0,
                     0.0,      0.0,      0.0,      0.5 - v,   0.0,       0.0,
                     0.0,      0.0,      0.0,      0.0,       0.5 - v,   0.0,
                     0.0,      0.0,      0.0,      0.0,       0.0,       0.5 - v;
                D *= E / ((1.0 + v) * (1.0 - 2.0 * v));

                for (Integer p = 0; p < CIntGaussIntegration<Real>::m_integPoints; p++)
                {
                    Real xi   = this->m_xi[p];
                    Real eta  = this->m_eta[p];
                    Real zeta = this->m_zeta[p];

                    Real wxi   = this->m_wxi[p];
                    Real weta  = this->m_weta[p];
                    Real wzeta = this->m_wzeta[p];

                    Real xieta = 1.0 - xi - eta;
                    Real zetam = 1.0 - zeta;
                    Real zetap = 1.0 + zeta;

                    // Jacobian (same as thermal, scaled by 0.5)
                    Eigen::Matrix<Real, 3, 3> J;
                    J(0, 0) = this->m_x2 * zetam - this->m_x1 * zetam + this->m_x5 * zetap - this->m_x4 * zetap;
                    J(0, 1) = this->m_y2 * zetam - this->m_y1 * zetam + this->m_y5 * zetap - this->m_y4 * zetap;
                    J(0, 2) = this->m_z2 * zetam - this->m_z1 * zetam + this->m_z5 * zetap - this->m_z4 * zetap;

                    J(1, 0) = this->m_x3 * zetam - this->m_x1 * zetam + this->m_x6 * zetap - this->m_x4 * zetap;
                    J(1, 1) = this->m_y3 * zetam - this->m_y1 * zetam + this->m_y6 * zetap - this->m_y4 * zetap;
                    J(1, 2) = this->m_z3 * zetam - this->m_z1 * zetam + this->m_z6 * zetap - this->m_z4 * zetap;

                    J(2, 0) = this->m_x5 * xi - this->m_x2 * xi + this->m_x4 * xieta - this->m_x1 * xieta + eta * this->m_x6 - eta * this->m_x3;
                    J(2, 1) = eta * this->m_y6 + xi * this->m_y5 + xieta * this->m_y4 - eta * this->m_y3 - xi * this->m_y2 - xieta * this->m_y1;
                    J(2, 2) = eta * this->m_z6 + xi * this->m_z5 + xieta * this->m_z4 - eta * this->m_z3 - xi * this->m_z2 - xieta * this->m_z1;

                    J *= 0.5;

                    Real detJ = J.determinant();
                    if (std::fabs(detJ) < 1E-15)
                        continue;

                    Eigen::Matrix<Real, 3, 3> Jinv = J.inverse();

                    // Natural derivatives of shape functions (×0.5 factor included)
                    // Rows: dN/dxi, dN/deta, dN/dzeta
                    // Cols: N0..N5
                    Eigen::Matrix<Real, 3, 6> G;
                    G << -zetam, +zetam, 0.0, -zetap, +zetap, 0.0,
                         -zetam, 0.0,   +zetam, -zetap, 0.0,  +zetap,
                         -xieta, -xi,   -eta,   +xieta, +xi,  +eta;
                    G *= 0.5;

                    // Physical derivatives: B_phys = Jinv * G (3x6)
                    Eigen::Matrix<Real, 3, 6> B_phys = Jinv * G;

                    // Strain-displacement matrix (6x18): DOF order [u,v,w] per node
                    Eigen::Matrix<Real, 6, 18> B;
                    B.setZero();
                    for (Integer i = 0; i < 6; ++i)
                    {
                        Real bx = B_phys(0, i);
                        Real by = B_phys(1, i);
                        Real bz = B_phys(2, i);

                        B(0, 3*i)   = bx;
                        B(1, 3*i+1) = by;
                        B(2, 3*i+2) = bz;
                        B(3, 3*i)   = by;
                        B(3, 3*i+1) = bx;
                        B(4, 3*i+1) = bz;
                        B(4, 3*i+2) = by;
                        B(5, 3*i)   = bz;
                        B(5, 3*i+2) = bx;
                    }

                    this->laplacian += B.transpose() * D * B * std::fabs(detJ) * wxi * weta * wzeta;
                }
            }

            template <typename Real>
            void CFemConstantStrainTriangularPrism<Real, 6, 3, 1>::setConvectiveTerm()
            {
            }

            template <typename Real>
            void CFemConstantStrainTriangularPrism<Real, 6, 3, 1>::update()
            {
                CFemTriangularPrism<Real, 6, 3, 1>::rebuild();

                this->setGaussPoints();

                setDiffusionTerm();

                if (CFemElement<Real>::transient())
                    setTransientTerm();
            }

            template <typename Real>
            void CFemConstantStrainTriangularPrism<Real, 6, 3, 1>::reCalcD()
            {
            }
        }
    }
}
