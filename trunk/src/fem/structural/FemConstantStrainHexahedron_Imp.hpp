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
            CFemConstantStrainHexahedron<Real, 8, 3, 1>::CFemConstantStrainHexahedron()
            {
                this->source.resize(8);
                this->source.setZero();
            }

            template <typename Real>
            CFemConstantStrainHexahedron<Real, 8, 3, 1>::~CFemConstantStrainHexahedron()
            {
            }

            template <typename Real>
            void CFemConstantStrainHexahedron<Real, 8, 3, 1>::setTransientTerm()
            {
            }

            template <typename Real>
            void CFemConstantStrainHexahedron<Real, 8, 3, 1>::setDiffusionTerm()
            {
                // Node ordering matches CMshBasicMesher hex output:
                // Node 0 (m_x1): xi=-1, eta=-1, zeta=-1
                // Node 1 (m_x2): xi=+1, eta=-1, zeta=-1
                // Node 2 (m_x3): xi=+1, eta=+1, zeta=-1
                // Node 3 (m_x4): xi=-1, eta=+1, zeta=-1
                // Node 4 (m_x5): xi=-1, eta=-1, zeta=+1
                // Node 5 (m_x6): xi=+1, eta=-1, zeta=+1
                // Node 6 (m_x7): xi=+1, eta=+1, zeta=+1
                // Node 7 (m_x8): xi=-1, eta=+1, zeta=+1

                this->laplacian.resize(24, 24);
                this->laplacian.setZero();

                Real v = CFemStructuralElement<Real>::coeffPoisson();
                Real E = CFemStructuralElement<Real>::elasticModulus();

                Eigen::Matrix<Real, 6, 6> D;
                D << 1.0 - v,  v,        v,        0.0,              0.0,              0.0,
                     v,        1.0 - v,  v,        0.0,              0.0,              0.0,
                     v,        v,        1.0 - v,  0.0,              0.0,              0.0,
                     0.0,      0.0,      0.0,      0.5 - v,          0.0,              0.0,
                     0.0,      0.0,      0.0,      0.0,              0.5 - v,          0.0,
                     0.0,      0.0,      0.0,      0.0,              0.0,              0.5 - v;
                D *= E / ((1.0 + v) * (1.0 - 2.0 * v));

                for (Integer p = 0; p < CIntGaussIntegration<Real>::m_integPoints; p++)
                {
                    Real xi   = this->m_xi[p];
                    Real eta  = this->m_eta[p];
                    Real zeta = this->m_zeta[p];

                    Real xim = 1.0 - xi,  xip = 1.0 + xi;
                    Real etam = 1.0 - eta, etap = 1.0 + eta;
                    Real zetam = 1.0 - zeta, zetap = 1.0 + zeta;

                    // Shape function derivatives in natural coordinates (×8, scaled later)
                    // dNi/dxi: [-etam*zetam, +etam*zetam, +etap*zetam, -etap*zetam,
                    //           -etam*zetap, +etam*zetap, +etap*zetap, -etap*zetap]
                    Real dN0dxi = -etam * zetam;
                    Real dN1dxi = +etam * zetam;
                    Real dN2dxi = +etap * zetam;
                    Real dN3dxi = -etap * zetam;
                    Real dN4dxi = -etam * zetap;
                    Real dN5dxi = +etam * zetap;
                    Real dN6dxi = +etap * zetap;
                    Real dN7dxi = -etap * zetap;

                    // dNi/deta: [-xim*zetam, -xip*zetam, +xip*zetam, +xim*zetam,
                    //            -xim*zetap, -xip*zetap, +xip*zetap, +xim*zetap]
                    Real dN0deta = -xim * zetam;
                    Real dN1deta = -xip * zetam;
                    Real dN2deta = +xip * zetam;
                    Real dN3deta = +xim * zetam;
                    Real dN4deta = -xim * zetap;
                    Real dN5deta = -xip * zetap;
                    Real dN6deta = +xip * zetap;
                    Real dN7deta = +xim * zetap;

                    // dNi/dzeta: [-xim*etam, -xip*etam, -xip*etap, -xim*etap,
                    //             +xim*etam, +xip*etam, +xip*etap, +xim*etap]
                    Real dN0dzeta = -xim * etam;
                    Real dN1dzeta = -xip * etam;
                    Real dN2dzeta = -xip * etap;
                    Real dN3dzeta = -xim * etap;
                    Real dN4dzeta = +xim * etam;
                    Real dN5dzeta = +xip * etam;
                    Real dN6dzeta = +xip * etap;
                    Real dN7dzeta = +xim * etap;

                    // Jacobian (3x3): J(row, col) = sum_i(dNi/d_row * coord_i_col) / 8
                    Eigen::Matrix<Real, 3, 3> J;
                    J(0, 0) = (dN0dxi*this->m_x1 + dN1dxi*this->m_x2 + dN2dxi*this->m_x3 + dN3dxi*this->m_x4 + dN4dxi*this->m_x5 + dN5dxi*this->m_x6 + dN6dxi*this->m_x7 + dN7dxi*this->m_x8) * 0.125;
                    J(0, 1) = (dN0dxi*this->m_y1 + dN1dxi*this->m_y2 + dN2dxi*this->m_y3 + dN3dxi*this->m_y4 + dN4dxi*this->m_y5 + dN5dxi*this->m_y6 + dN6dxi*this->m_y7 + dN7dxi*this->m_y8) * 0.125;
                    J(0, 2) = (dN0dxi*this->m_z1 + dN1dxi*this->m_z2 + dN2dxi*this->m_z3 + dN3dxi*this->m_z4 + dN4dxi*this->m_z5 + dN5dxi*this->m_z6 + dN6dxi*this->m_z7 + dN7dxi*this->m_z8) * 0.125;

                    J(1, 0) = (dN0deta*this->m_x1 + dN1deta*this->m_x2 + dN2deta*this->m_x3 + dN3deta*this->m_x4 + dN4deta*this->m_x5 + dN5deta*this->m_x6 + dN6deta*this->m_x7 + dN7deta*this->m_x8) * 0.125;
                    J(1, 1) = (dN0deta*this->m_y1 + dN1deta*this->m_y2 + dN2deta*this->m_y3 + dN3deta*this->m_y4 + dN4deta*this->m_y5 + dN5deta*this->m_y6 + dN6deta*this->m_y7 + dN7deta*this->m_y8) * 0.125;
                    J(1, 2) = (dN0deta*this->m_z1 + dN1deta*this->m_z2 + dN2deta*this->m_z3 + dN3deta*this->m_z4 + dN4deta*this->m_z5 + dN5deta*this->m_z6 + dN6deta*this->m_z7 + dN7deta*this->m_z8) * 0.125;

                    J(2, 0) = (dN0dzeta*this->m_x1 + dN1dzeta*this->m_x2 + dN2dzeta*this->m_x3 + dN3dzeta*this->m_x4 + dN4dzeta*this->m_x5 + dN5dzeta*this->m_x6 + dN6dzeta*this->m_x7 + dN7dzeta*this->m_x8) * 0.125;
                    J(2, 1) = (dN0dzeta*this->m_y1 + dN1dzeta*this->m_y2 + dN2dzeta*this->m_y3 + dN3dzeta*this->m_y4 + dN4dzeta*this->m_y5 + dN5dzeta*this->m_y6 + dN6dzeta*this->m_y7 + dN7dzeta*this->m_y8) * 0.125;
                    J(2, 2) = (dN0dzeta*this->m_z1 + dN1dzeta*this->m_z2 + dN2dzeta*this->m_z3 + dN3dzeta*this->m_z4 + dN4dzeta*this->m_z5 + dN5dzeta*this->m_z6 + dN6dzeta*this->m_z7 + dN7dzeta*this->m_z8) * 0.125;

                    Real detJ = J.determinant();
                    if (std::fabs(detJ) < 1E-15)
                        continue;

                    Eigen::Matrix<Real, 3, 3> Jinv = J.inverse();

                    // Physical derivatives: [bix, biy, biz]^T = Jinv * [dNi/dxi, dNi/deta, dNi/dzeta]^T / 8
                    auto phys = [&](Real dxi, Real deta, Real dzeta, Real& bx, Real& by, Real& bz) {
                        dxi *= 0.125; deta *= 0.125; dzeta *= 0.125;
                        bx = Jinv(0,0)*dxi + Jinv(0,1)*deta + Jinv(0,2)*dzeta;
                        by = Jinv(1,0)*dxi + Jinv(1,1)*deta + Jinv(1,2)*dzeta;
                        bz = Jinv(2,0)*dxi + Jinv(2,1)*deta + Jinv(2,2)*dzeta;
                    };

                    Real bx[8], by[8], bz[8];
                    phys(dN0dxi, dN0deta, dN0dzeta, bx[0], by[0], bz[0]);
                    phys(dN1dxi, dN1deta, dN1dzeta, bx[1], by[1], bz[1]);
                    phys(dN2dxi, dN2deta, dN2dzeta, bx[2], by[2], bz[2]);
                    phys(dN3dxi, dN3deta, dN3dzeta, bx[3], by[3], bz[3]);
                    phys(dN4dxi, dN4deta, dN4dzeta, bx[4], by[4], bz[4]);
                    phys(dN5dxi, dN5deta, dN5dzeta, bx[5], by[5], bz[5]);
                    phys(dN6dxi, dN6deta, dN6dzeta, bx[6], by[6], bz[6]);
                    phys(dN7dxi, dN7deta, dN7dzeta, bx[7], by[7], bz[7]);

                    // Strain-displacement matrix B (6x24)
                    // DOF order per node: [u, v, w] -> cols 3i, 3i+1, 3i+2
                    // Rows: [eps_xx, eps_yy, eps_zz, gam_xy, gam_yz, gam_xz]
                    Eigen::Matrix<Real, 6, 24> B;
                    B.setZero();
                    for (Integer i = 0; i < 8; ++i)
                    {
                        B(0, 3*i)   = bx[i];
                        B(1, 3*i+1) = by[i];
                        B(2, 3*i+2) = bz[i];
                        B(3, 3*i)   = by[i];
                        B(3, 3*i+1) = bx[i];
                        B(4, 3*i+1) = bz[i];
                        B(4, 3*i+2) = by[i];
                        B(5, 3*i)   = bz[i];
                        B(5, 3*i+2) = bx[i];
                    }

                    this->laplacian += B.transpose() * D * B * std::fabs(detJ) * this->m_wxi[p] * this->m_weta[p] * this->m_wzeta[p];
                }
            }

            template <typename Real>
            void CFemConstantStrainHexahedron<Real, 8, 3, 1>::setConvectiveTerm()
            {
            }

            template <typename Real>
            void CFemConstantStrainHexahedron<Real, 8, 3, 1>::update()
            {
                // Calculate geometrical properties
                CFemHexahedron<Real, 8, 3, 1>::rebuild();

                // Integration points
                this->setGaussPoints();

                // Diffusion term
                setDiffusionTerm();

                // Transient term
                if (CFemElement<Real>::transient())
                    setTransientTerm();
            }

            template <typename Real>
            void CFemConstantStrainHexahedron<Real, 8, 3, 1>::reCalcD()
            {
            }
        }
    }
}
