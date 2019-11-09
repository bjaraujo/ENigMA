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

    namespace flow {

        template <typename Real>
        CFemFlowTriangle<Real, 3, 1, 1>::CFemFlowTriangle()
        {
        }

        template <typename Real>
        CFemFlowTriangle<Real, 3, 1, 1>::~CFemFlowTriangle()
        {
        }

        template <typename Real>
        void CFemFlowTriangle<Real, 3, 1, 1>::setDiffusionTerm()
        {

            CFemTriangle<Real, 3, 1, 1>::setDiffusionTerm();
        }

        template <typename Real>
        void CFemFlowTriangle<Real, 3, 1, 1>::setConvectiveTerm(Real u[3], Real v[3])
        {

            Eigen::Matrix<Real, 2, 3> B;

            CFemTriangle<Real, 3, 1, 1>::calculateB(B);

            Real bi = B(0, 0);
            Real bj = B(0, 1);
            Real bk = B(0, 2);

            Real ci = B(1, 0);
            Real cj = B(1, 1);
            Real ck = B(1, 2);

            Eigen::Matrix<Real, 3, 3> C;

            Real us = u[0] + u[1] + u[2];
            Real vs = v[0] + v[1] + v[2];

            C << (us + u[0]) * bi + (vs + v[0]) * ci, (us + u[0]) * bj + (vs + v[0]) * cj, (us + u[0]) * bk + (vs + v[0]) * ck,
                (us + u[1]) * bi + (vs + v[1]) * ci, (us + u[1]) * bj + (vs + v[1]) * cj, (us + u[1]) * bk + (vs + v[1]) * ck,
                (us + u[2]) * bi + (vs + v[2]) * ci, (us + u[2]) * bj + (vs + v[2]) * cj, (us + u[2]) * bk + (vs + v[2]) * ck;

            C *= 1.0 / 24.0;

            Eigen::Matrix<Real, 3, 3> Ks1, Ks2, Ks3, Ks4;

            Real ua = us / 3.0;
            Real va = vs / 3.0;

            Ks1 << bi * bi, bi * bj, bi * bk,
                bj * bi, bj * bj, bj * bk,
                bk * bi, bk * bj, bk * bk;

            Ks1 *= ua * us / (12.0 * CGeoTriangle<Real>::area());

            Ks2 << bi * ci, bi * cj, bi * ck,
                bj * ci, bj * cj, bj * ck,
                bk * ci, bk * cj, bk * ck;

            Ks2 *= ua * vs / (12.0 * CGeoTriangle<Real>::area());

            Ks3 << ci * bi, ci * bj, ci * bk,
                cj * bi, cj * bj, cj * bk,
                ck * bi, ck * bj, ck * bk;

            Ks3 *= va * us / (12.0 * CGeoTriangle<Real>::area());

            Ks4 << ci * ci, ci * cj, ci * ck,
                cj * ci, cj * cj, cj * ck,
                ck * ci, ck * cj, ck * ck;

            Ks4 *= va * vs / (12.0 * CGeoTriangle<Real>::area());

            CFemElement<Real>::divergence = C - CFemElement<Real>::m_dt * (Ks1 + Ks2 + Ks3 + Ks4);
        }

        template <typename Real>
        void CFemFlowTriangle<Real, 3, 1, 1>::setGradientTerm(const Integer aComponent)
        {

            Eigen::Matrix<Real, 2, 3> B;

            CFemTriangle<Real, 3, 1, 1>::calculateB(B);

            Real bi = B(0, 0);
            Real bj = B(0, 1);
            Real bk = B(0, 2);

            Real ci = B(1, 0);
            Real cj = B(1, 1);
            Real ck = B(1, 2);

            if (aComponent == 0) {

                Eigen::Matrix<Real, 3, 3> G1;

                G1 << bi, bj, bk,
                    bi, bj, bk,
                    bi, bj, bk;

                G1 *= 1.0 / 6.0;

                CFemElement<Real>::gradient = G1;

            } else if (aComponent == 1) {

                Eigen::Matrix<Real, 3, 3> G2;

                G2 << ci, cj, ck,
                    ci, cj, ck,
                    ci, cj, ck;

                G2 *= 1.0 / 6.0;

                CFemElement<Real>::gradient = G2;
            }
        }

        template <typename Real>
        void CFemFlowTriangle<Real, 3, 1, 1>::calculateTransientTerm()
        {

            // Calculate geometrical properties
            CFemTriangle<Real, 3, 1, 1>::rebuild();

            // Transient term
            CFemTriangle<Real, 3, 1, 1>::setTransientTerm();
        }

        template <typename Real>
        void CFemFlowTriangle<Real, 3, 1, 1>::calculateDiffusiveTerm()
        {

            // Calculate geometrical properties
            CFemTriangle<Real, 3, 1, 1>::rebuild();

            // Diffusion term
            setDiffusionTerm();
        }

        template <typename Real>
        void CFemFlowTriangle<Real, 3, 1, 1>::calculateConvectiveTerm(Real u[3], Real v[3])
        {

            // Calculate geometrical properties
            CFemTriangle<Real, 3, 1, 1>::rebuild();

            // Convective term
            setConvectiveTerm(u, v);
        }

        template <typename Real>
        void CFemFlowTriangle<Real, 3, 1, 1>::calculateGradientTerm(const Integer aComponent)
        {

            // Calculate geometrical properties
            CFemTriangle<Real, 3, 1, 1>::rebuild();

            // Convective term
            setGradientTerm(aComponent);
        }
    }
}
}
