// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// <Copyright> Copyright (c) 2012, All Rights Reserved </Copyright>
// *****************************************************************************

#pragma once

using namespace ENigMA::geometry;

namespace ENigMA
{

    namespace fem
    {

        namespace flow
        {

            template <typename Real>
            CFemFlowTetrahedron<Real, 4, 1, 1>::CFemFlowTetrahedron()
            {
                
            }

            template <typename Real>
            CFemFlowTetrahedron<Real, 4, 1, 1>::~CFemFlowTetrahedron()
            {

            }

            template <typename Real>
            void CFemFlowTetrahedron<Real, 4, 1, 1>::setDiffusionTerm()
            {

                CFemTetrahedron<Real, 4, 1, 1>::setDiffusionTerm();

            }

            template <typename Real>
            void CFemFlowTetrahedron<Real, 4, 1, 1>::setConvectiveTerm(Real u[4], Real v[4], Real w[4])
            {

                Eigen::Matrix<Real, 3, 4> B;

                CFemTetrahedron<Real, 4, 1, 1>::calculateB(B);

                Real bi = B(0,0);
                Real bj = B(0,1);
                Real bk = B(0,2);
                Real bl = B(0,3);

                Real ci = B(1,0);
                Real cj = B(1,1);
                Real ck = B(1,2);
                Real cl = B(1,3);

                Real di = B(2,0);
                Real dj = B(2,1);
                Real dk = B(2,2);
                Real dl = B(2,3);

                Eigen::Matrix<Real, 4, 4> C;

                Real us = u[0] + u[1] + u[2] + u[3];
                Real vs = v[0] + v[1] + v[2] + v[3];
                Real ws = w[0] + w[1] + w[2] + w[3];

                C << (us + u[0]) * bi + (vs + v[0]) * ci + (ws + w[0]) * di, (us + u[0]) * bj + (vs + v[0]) * cj + (ws + w[0]) * dj, (us + u[0]) * bk + (vs + v[0]) * ck + (ws + w[0]) * dk, (us + u[0]) * bl + (vs + v[0]) * cl + (ws + w[0]) * dl,
                     (us + u[1]) * bi + (vs + v[1]) * ci + (ws + w[1]) * di, (us + u[1]) * bj + (vs + v[1]) * cj + (ws + w[1]) * dj, (us + u[1]) * bk + (vs + v[1]) * ck + (ws + w[1]) * dk, (us + u[1]) * bl + (vs + v[1]) * cl + (ws + w[1]) * dl,
                     (us + u[2]) * bi + (vs + v[2]) * ci + (ws + w[2]) * di, (us + u[2]) * bj + (vs + v[2]) * cj + (ws + w[2]) * dj, (us + u[2]) * bk + (vs + v[2]) * ck + (ws + w[2]) * dk, (us + u[2]) * bl + (vs + v[2]) * cl + (ws + w[2]) * dl,
                     (us + u[3]) * bi + (vs + v[3]) * ci + (ws + w[3]) * di, (us + u[3]) * bj + (vs + v[3]) * cj + (ws + w[3]) * dj, (us + u[2]) * bk + (vs + v[3]) * ck + (ws + w[3]) * dk, (us + u[3]) * bl + (vs + v[3]) * cl + (ws + w[3]) * dl;
                
                C *= 1.0 / 120.0;

                Eigen::Matrix<Real, 4, 4> Ks1, Ks2, Ks3, Ks4, Ks5, Ks6, Ks7, Ks8, Ks9;

                Real ua = us / 4.0;
                Real va = vs / 4.0;
                Real wa = ws / 4.0;

                Ks1 << bi * bi, bi * bj, bi * bk, bi * bl,
                       bj * bi, bj * bj, bj * bk, bj * bl,
                       bk * bi, bk * bj, bk * bk, bk * bl,
                       bl * bi, bl * bj, bl * bk, bl * bl;

                Ks1 *= ua * us / (144.0 * CGeoTetrahedron<Real>::volume());

                Ks2 << bi * ci, bi * cj, bi * ck, bi * cl,
                       bj * ci, bj * cj, bj * ck, bj * cl,
                       bk * ci, bk * cj, bk * ck, bk * cl,
                       bl * ci, bl * cj, bl * ck, bl * cl;

                Ks2 *= ua * vs / (144.0 * CGeoTetrahedron<Real>::volume());

                Ks3 << bi * di, bi * dj, bi * dk, bi * dl,
                       bj * di, bj * dj, bj * dk, bj * dl,
                       bk * di, bk * dj, bk * dk, bk * dl,
                       bl * di, bl * dj, bl * dk, bl * dl;

                Ks3 *= ua * ws / (144.0 * CGeoTetrahedron<Real>::volume());

                Ks4 << ci * bi, ci * bj, ci * bk, ci * bl,
                       cj * bi, cj * bj, cj * bk, cj * bl,
                       ck * bi, ck * bj, ck * bk, ck * bl,
                       cl * bi, cl * bj, cl * bk, cl * bl;

                Ks4 *= va * us / (144.0 * CGeoTetrahedron<Real>::volume());

                Ks5 << ci * ci, ci * cj, ci * ck, ci * cl,
                       cj * ci, cj * cj, cj * ck, cj * cl,
                       ck * ci, ck * cj, ck * ck, ck * cl,
                       cl * ci, cl * cj, cl * ck, cl * cl;

                Ks5 *= va * vs / (144.0 * CGeoTetrahedron<Real>::volume());

                Ks6 << ci * di, ci * dj, ci * dk, ci * dl,
                       cj * di, cj * dj, cj * dk, cj * dl,
                       ck * di, ck * dj, ck * dk, ck * dl,
                       cl * di, cl * dj, cl * dk, cl * dl;

                Ks6 *= va * ws / (144.0 * CGeoTetrahedron<Real>::volume());

                Ks7 << di * bi, di * bj, di * bk, di * bl,
                       dj * bi, dj * bj, dj * bk, dj * bl,
                       dk * bi, dk * bj, dk * bk, dk * bl,
                       dl * bi, dl * bj, dl * bk, dl * bl;

                Ks7 *= wa * us / (144.0 * CGeoTetrahedron<Real>::volume());

                Ks8 << di * ci, di * cj, di * ck, di * cl,
                       dj * ci, dj * cj, dj * ck, dj * cl,
                       dk * ci, dk * cj, dk * ck, dk * cl,
                       dl * ci, dl * cj, dl * ck, dl * cl;

                Ks8 *= wa * vs / (144.0 * CGeoTetrahedron<Real>::volume());

                Ks9 << di * di, di * dj, di * dk, di * dl,
                       dj * di, dj * dj, dj * dk, dj * dl,
                       dk * di, dk * dj, dk * dk, dk * dl,
                       dl * di, dl * dj, dl * dk, dl * dl;

                Ks9 *= wa * ws / (144.0 * CGeoTetrahedron<Real>::volume());

                CFemElement<Real>::divergence = C - CFemElement<Real>::m_dt * (Ks1 + Ks2 + Ks3 + Ks4 + Ks5 + Ks6 + Ks7 + Ks8 + Ks9);

            }

            template <typename Real>
            void CFemFlowTetrahedron<Real, 4, 1, 1>::setGradientTerm(const Integer aComponent)
            {

                Eigen::Matrix<Real, 3, 4> B;

                CFemTetrahedron<Real, 4, 1, 1>::calculateB(B);

                Real bi = B(0,0);
                Real bj = B(0,1);
                Real bk = B(0,2);
                Real bl = B(0,3);

                Real ci = B(1,0);
                Real cj = B(1,1);
                Real ck = B(1,2);
                Real cl = B(1,3);

                Real di = B(2,0);
                Real dj = B(2,1);
                Real dk = B(2,2);
                Real dl = B(2,3);

                if (aComponent == 0)
                {

                    Eigen::Matrix<Real, 4, 4> G1;

                    G1 << bi, bj, bk, bl,
                          bi, bj, bk, bl,
                          bi, bj, bk, bl,
                          bi, bj, bk, bl;

                    G1 *= 1.0 / 24.0;

                    CFemElement<Real>::gradient = G1;

                }
                else if (aComponent == 1)
                {

                    Eigen::Matrix<Real, 4, 4> G2;

                    G2 << ci, cj, ck, cl,
                          ci, cj, ck, cl,
                          ci, cj, ck, cl,
                          ci, cj, ck, cl;

                    G2 *= 1.0 / 24.0;

                    CFemElement<Real>::gradient = G2;

                }
                else if (aComponent == 2)
                {

                    Eigen::Matrix<Real, 4, 4> G3;

                    G3 << di, dj, dk, dl,
                          di, dj, dk, dl,
                          di, dj, dk, dl,
                          di, dj, dk, dl;

                    G3 *= 1.0 / 24.0;

                    CFemElement<Real>::gradient = G3;

                }

            }

            template <typename Real>
            void CFemFlowTetrahedron<Real, 4, 1, 1>::calculateTransientTerm()
            {

                // Calculate geometrical properties
                CFemTetrahedron<Real, 4, 1, 1>::rebuild();

                // Transient term
                CFemTetrahedron<Real, 4, 1, 1>::setTransientTerm();

            }

            template <typename Real>
            void CFemFlowTetrahedron<Real, 4, 1, 1>::calculateDiffusiveTerm()
            {

                // Calculate geometrical properties
                CFemTetrahedron<Real, 4, 1, 1>::rebuild();

                // Diffusion term
                setDiffusionTerm();

            }

            template <typename Real>
            void CFemFlowTetrahedron<Real, 4, 1, 1>::calculateConvectiveTerm(Real u[3], Real v[3], Real w[3])
            {

                // Calculate geometrical properties
                CFemTetrahedron<Real, 4, 1, 1>::rebuild();

                // Convective term
                setConvectiveTerm(u, v, w);

            }

            template <typename Real>
            void CFemFlowTetrahedron<Real, 4, 1, 1>::calculateGradientTerm(const Integer aComponent)
            {

                // Calculate geometrical properties
                CFemTetrahedron<Real, 4, 1, 1>::rebuild();

                // Convective term
                setGradientTerm(aComponent);

            }

        }

    }

}
