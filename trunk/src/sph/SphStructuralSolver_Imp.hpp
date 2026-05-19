// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include <iostream>
#include <Eigen/SparseLU>

namespace ENigMA
{
    namespace sph
    {
        template <typename Real>
        CSphStructuralSolver<Real>::CSphStructuralSolver()
            : m_kernel(nullptr), m_E(1.0), m_nu(0.3), m_h(1.0), 
              m_lambda(0.0), m_mu(0.0), m_nParticles(0)
        {
        }

        template <typename Real>
        CSphStructuralSolver<Real>::~CSphStructuralSolver()
        {
        }

        template <typename Real>
        void CSphStructuralSolver<Real>::setKernel(CSphKernel<Real>* kernel)
        {
            m_kernel = kernel;
        }

        template <typename Real>
        void CSphStructuralSolver<Real>::setSmoothingLength(const Real h)
        {
            m_h = h;
        }

        template <typename Real>
        void CSphStructuralSolver<Real>::setMaterialProperties(const Real E, const Real nu)
        {
            m_E = E;
            m_nu = nu;

            // Calculate Lamé parameters
            // λ = E·ν / ((1+ν)·(1-2ν))
            // μ = E / (2·(1+ν))
            m_lambda = m_E * m_nu / ((1.0 + m_nu) * (1.0 - 2.0 * m_nu));
            m_mu = m_E / (2.0 * (1.0 + m_nu));
        }

        template <typename Real>
        void CSphStructuralSolver<Real>::addParticle(const ENigMA::geometry::CGeoCoordinate<Real>& p)
        {
            m_points.push_back(p);
            m_displacement.push_back(ENigMA::geometry::CGeoVector<Real>(0.0, 0.0, 0.0));
            m_strain.push_back(Eigen::Matrix<Real, 3, 3>::Zero());
            m_stress.push_back(Eigen::Matrix<Real, 3, 3>::Zero());
            m_externalForce.push_back(Eigen::Matrix<Real, 3, 1>::Zero());
            m_isFixed.push_back(false);
            m_fixedDisplacement.push_back(ENigMA::geometry::CGeoVector<Real>(0.0, 0.0, 0.0));

            m_nParticles = m_points.size();
        }

        template <typename Real>
        void CSphStructuralSolver<Real>::addParticles(const std::vector<ENigMA::geometry::CGeoCoordinate<Real>>& points)
        {
            for (const auto& p : points)
                addParticle(p);
        }

        template <typename Real>
        Integer CSphStructuralSolver<Real>::nbParticles() const
        {
            return m_nParticles;
        }

        template <typename Real>
        void CSphStructuralSolver<Real>::setFixedBoundaryCondition(
            const Integer ptId, const ENigMA::geometry::CGeoVector<Real>& u)
        {
            if (ptId < 0 || ptId >= m_nParticles)
                return;

            m_isFixed[ptId] = true;
            m_fixedDisplacement[ptId] = u;
            m_displacement[ptId] = u;
        }

        template <typename Real>
        void CSphStructuralSolver<Real>::setForce(const Integer ptId, const ENigMA::geometry::CGeoVector<Real>& f)
        {
            if (ptId < 0 || ptId >= m_nParticles)
                return;

            m_externalForce[ptId] = f;
        }

        template <typename Real>
        void CSphStructuralSolver<Real>::calculateStrain()
        {
            if (m_kernel == nullptr)
                return;

            const Real tolerance = m_h * 1e-6;

            for (Integer i = 0; i < m_nParticles; ++i)
            {
                Eigen::Matrix<Real, 3, 3> strainRate = Eigen::Matrix<Real, 3, 3>::Zero();

                // Find neighbors using hash grid
                std::vector<Integer> neighbors;
                m_hashGrid.find(neighbors, m_points[i], m_h);

                for (Integer nIdx : neighbors)
                {
                    if (nIdx == i)
                        continue;

                    ENigMA::geometry::CGeoVector<Real> r = m_points[nIdx] - m_points[i];
                    Real dist = r.norm();

                    if (dist > m_h || dist < tolerance)
                        continue;

                    // Gradient of kernel: ∇W_ij
                    ENigMA::geometry::CGeoVector<Real> gradW = m_kernel->gradientW(r, m_h, tolerance);

                    // Displacement difference: u_j - u_i
                    Eigen::Matrix<Real, 3, 1> du(
                        m_displacement[nIdx].x() - m_displacement[i].x(),
                        m_displacement[nIdx].y() - m_displacement[i].y(),
                        m_displacement[nIdx].z() - m_displacement[i].z()
                    );

                    // Volume per particle (assume unit volume for now)
                    Real V_j = 1.0;

                    // Contribution to strain rate: u_j ⊗ ∇W_ij
                    // ∂u_d/∂x_k ≈ Σⱼ (u_d,j - u_d,i) ∂W_ij/∂x_k V_j
                    for (int d = 0; d < 3; ++d)
                    {
                        strainRate(d, 0) += du(d) * gradW.x() * V_j;
                        strainRate(d, 1) += du(d) * gradW.y() * V_j;
                        strainRate(d, 2) += du(d) * gradW.z() * V_j;
                    }
                }

                // Symmetrize: ε = 1/2(∇u + ∇u^T)
                m_strain[i] = 0.5 * (strainRate + strainRate.transpose());
            }
        }

        template <typename Real>
        void CSphStructuralSolver<Real>::calculateStress()
        {
            // σ = λ tr(ε) I + 2μ ε
            for (Integer i = 0; i < m_nParticles; ++i)
            {
                Real traceStrain = m_strain[i].trace();

                m_stress[i] = m_lambda * traceStrain * Eigen::Matrix<Real, 3, 3>::Identity() +
                              2.0 * m_mu * m_strain[i];
            }
        }

        template <typename Real>
        Eigen::Matrix<Real, Eigen::Dynamic, 1> CSphStructuralSolver<Real>::assembleResidual()
        {
            if (m_kernel == nullptr)
                return Eigen::Matrix<Real, Eigen::Dynamic, 1>::Zero(m_nParticles * 3);

            Integer nDofs = m_nParticles * 3;
            Eigen::Matrix<Real, Eigen::Dynamic, 1> residual = Eigen::Matrix<Real, Eigen::Dynamic, 1>::Zero(nDofs);

            const Real tolerance = m_h * 1e-6;

            for (Integer i = 0; i < m_nParticles; ++i)
            {
                // Force: F_i = -∫ ∇·σ_i dV + f_ext,i
                Eigen::Matrix<Real, 3, 1> force = Eigen::Matrix<Real, 3, 1>::Zero();

                // Internal force from stress divergence
                std::vector<Integer> neighbors;
                m_hashGrid.find(neighbors, m_points[i], m_h);

                for (Integer nIdx : neighbors)
                {
                    if (nIdx == i)
                        continue;

                    ENigMA::geometry::CGeoVector<Real> r = m_points[nIdx] - m_points[i];
                    Real dist = r.norm();

                    if (dist > m_h || dist < tolerance)
                        continue;

                    ENigMA::geometry::CGeoVector<Real> gradW = m_kernel->gradientW(r, m_h, tolerance);

                    // Symmetric formulation for internal forces:
                    // F_i = -Σⱼ (σ_i + σ_j)·∇W_ij V_j / 2
                    // (stabilized SPH formulation)
                    Real V_j = 1.0;

                    Eigen::Matrix<Real, 3, 1> sig_i_gradW = m_stress[i] * Eigen::Matrix<Real, 3, 1>(
                        gradW.x(), gradW.y(), gradW.z());
                    Eigen::Matrix<Real, 3, 1> sig_j_gradW = m_stress[nIdx] * Eigen::Matrix<Real, 3, 1>(
                        gradW.x(), gradW.y(), gradW.z());

                    force -= 0.5 * (sig_i_gradW + sig_j_gradW) * V_j;
                }

                // Add external forces
                force(0) += m_externalForce[i](0);
                force(1) += m_externalForce[i](1);
                force(2) += m_externalForce[i](2);

                residual(i * 3 + 0) = force(0);
                residual(i * 3 + 1) = force(1);
                residual(i * 3 + 2) = force(2);
            }

            return residual;
        }

        template <typename Real>
        void CSphStructuralSolver<Real>::buildStiffnessMatrix()
        {
            // Build approximate stiffness matrix using perturbation method
            // K_ij ≈ ∂F_i/∂u_j
            // For linear elasticity with small deformations, we use constant stiffness

            Integer nDofs = m_nParticles * 3;
            std::vector<Eigen::Triplet<Real>> triplets;

            const Real tolerance = m_h * 1e-6;
            const Real delta = 1e-8;

            // For efficiency, use approximate stiffness based on material properties
            // This is a simplified approach - for better accuracy, compute full Jacobian
            for (Integer i = 0; i < m_nParticles; ++i)
            {
                std::vector<Integer> neighbors;
                m_hashGrid.find(neighbors, m_points[i], m_h);

                for (Integer j : neighbors)
                {
                    if (i == j)
                        continue;

                    ENigMA::geometry::CGeoVector<Real> r = m_points[j] - m_points[i];
                    Real dist = r.norm();

                    if (dist > m_h || dist < tolerance)
                        continue;

                    ENigMA::geometry::CGeoVector<Real> gradW = m_kernel->gradientW(r, m_h, tolerance);
                    Real V_j = 1.0;

                    // Approximate stiffness contribution
                    // K_ij ≈ (2μ + λ) * ∇W_ij^T * ∇W_ij * V_j
                    Real diag_coeff = (2.0 * m_mu + m_lambda) * V_j;
                    Real off_diag_coeff = m_lambda * V_j;

                    for (int d = 0; d < 3; ++d)
                    {
                        // Diagonal terms
                        triplets.push_back(Eigen::Triplet<Real>(
                            i * 3 + d, i * 3 + d,
                            diag_coeff * (gradW.x() * gradW.x() +
                                          gradW.y() * gradW.y() +
                                          gradW.z() * gradW.z())));

                        // Off-diagonal coupling
                        for (int e = d + 1; e < 3; ++e)
                        {
                            Real coupling = 0.0;
                            if (d == 0 && e == 1) coupling = gradW.x() * gradW.y();
                            if (d == 0 && e == 2) coupling = gradW.x() * gradW.z();
                            if (d == 1 && e == 2) coupling = gradW.y() * gradW.z();
                            coupling *= off_diag_coeff;

                            triplets.push_back(Eigen::Triplet<Real>(i * 3 + d, i * 3 + e, coupling));
                            triplets.push_back(Eigen::Triplet<Real>(i * 3 + e, i * 3 + d, coupling));
                        }
                    }
                }
            }

            m_stiffnessMatrix.resize(nDofs, nDofs);
            m_stiffnessMatrix.setFromTriplets(triplets.begin(), triplets.end());
            m_stiffnessMatrix.makeCompressed();
        }

        template <typename Real>
        void CSphStructuralSolver<Real>::applyBoundaryConditions(
            Eigen::SparseMatrix<Real, Eigen::RowMajor>& K,
            Eigen::Matrix<Real, Eigen::Dynamic, 1>& F)
        {
            // Apply boundary conditions using penalty method
            const Real penalty = 1e10 * m_E;  // Large penalty parameter

            Integer nDofs = m_nParticles * 3;

            // Create modified stiffness matrix with penalty terms
            std::vector<Eigen::Triplet<Real>> triplets;

            // First, copy existing matrix entries
            for (Integer k = 0; k < K.outerSize(); ++k)
            {
                for (Eigen::SparseMatrix<Real, Eigen::RowMajor>::InnerIterator it(K, k); it; ++it)
                {
                    Integer row = it.row();
                    Integer col = it.col();
                    Integer ptId = row / 3;

                    // Only add entry if not a fixed DOF, or if it's a diagonal term for a fixed DOF
                    if (m_isFixed[ptId] && row == col)
                    {
                        // Add penalty to diagonal
                        triplets.push_back(Eigen::Triplet<Real>(row, col, it.value() + penalty));
                    }
                    else if (!m_isFixed[ptId])
                    {
                        triplets.push_back(Eigen::Triplet<Real>(row, col, it.value()));
                    }
                }
            }

            // Add penalty terms for fixed DOFs that don't have entries yet
            for (Integer i = 0; i < m_nParticles; ++i)
            {
                if (m_isFixed[i])
                {
                    for (int d = 0; d < 3; ++d)
                    {
                        Integer dofId = i * 3 + d;
                        bool hasEntry = false;

                        // Check if this DOF already has a diagonal entry
                        for (const auto& t : triplets)
                        {
                            if (t.row() == dofId && t.col() == dofId)
                            {
                                hasEntry = true;
                                break;
                            }
                        }

                        if (!hasEntry)
                            triplets.push_back(Eigen::Triplet<Real>(dofId, dofId, penalty));
                    }
                }
            }

            K.setFromTriplets(triplets.begin(), triplets.end());

            // Modify force vector for fixed DOFs
            for (Integer i = 0; i < m_nParticles; ++i)
            {
                if (m_isFixed[i])
                {
                    F(i * 3 + 0) = penalty * m_fixedDisplacement[i].x();
                    F(i * 3 + 1) = penalty * m_fixedDisplacement[i].y();
                    F(i * 3 + 2) = penalty * m_fixedDisplacement[i].z();
                }
            }
        }

        template <typename Real>
        void CSphStructuralSolver<Real>::solve(const Integer maxIterations, const Real tolerance)
        {
            if (m_kernel == nullptr || m_nParticles == 0)
                return;

            // Build spatial hash grid for neighbor search
            m_hashGrid.reset();
            for (Integer i = 0; i < m_nParticles; ++i)
                m_hashGrid.addGeometricObject(i, m_points[i]);
            m_hashGrid.build();

            Integer nDofs = m_nParticles * 3;

            for (Integer iter = 0; iter < maxIterations; ++iter)
            {
                // Calculate strain and stress
                calculateStrain();
                calculateStress();

                // Assemble residual
                Eigen::Matrix<Real, Eigen::Dynamic, 1> residual = assembleResidual();
                Real res_norm = residual.norm();

                std::cout << "SPH Solver - Iteration " << iter << ": residual = "
                          << std::scientific << res_norm << std::endl;

                if (res_norm < tolerance)
                {
                    std::cout << "Convergence achieved at iteration " << iter << std::endl;
                    break;
                }

                // Build stiffness matrix
                buildStiffnessMatrix();

                // Apply boundary conditions
                Eigen::SparseMatrix<Real, Eigen::RowMajor> K = m_stiffnessMatrix;
                Eigen::Matrix<Real, Eigen::Dynamic, 1> F = -residual;

                applyBoundaryConditions(K, F);

                // Solve linear system: K * du = F
                Eigen::SparseLU<Eigen::SparseMatrix<Real>> solver;
                solver.compute(K);

                if (solver.info() != Eigen::Success)
                {
                    std::cout << "Sparse LU decomposition failed!" << std::endl;
                    break;
                }

                Eigen::Matrix<Real, Eigen::Dynamic, 1> du = solver.solve(F);

                if (solver.info() != Eigen::Success)
                {
                    std::cout << "Solving linear system failed!" << std::endl;
                    break;
                }

                // Update displacement
                for (Integer i = 0; i < m_nParticles; ++i)
                {
                    if (!m_isFixed[i])
                    {
                        m_displacement[i] += ENigMA::geometry::CGeoVector<Real>(
                            du(i * 3 + 0), du(i * 3 + 1), du(i * 3 + 2));
                    }
                }

                // Line search with simple relaxation
                Real relaxation = 0.7;
                for (Integer i = 0; i < m_nParticles; ++i)
                {
                    if (!m_isFixed[i])
                    {
                        m_displacement[i] = m_displacement[i] * relaxation;
                    }
                }
            }
        }

        template <typename Real>
        ENigMA::geometry::CGeoVector<Real> CSphStructuralSolver<Real>::getDisplacement(const Integer ptId) const
        {
            if (ptId < 0 || ptId >= m_nParticles)
                return ENigMA::geometry::CGeoVector<Real>(0.0, 0.0, 0.0);

            return m_displacement[ptId];
        }

        template <typename Real>
        Eigen::Matrix<Real, 3, 3> CSphStructuralSolver<Real>::getStrain(const Integer ptId) const
        {
            if (ptId < 0 || ptId >= m_nParticles)
                return Eigen::Matrix<Real, 3, 3>::Zero();

            return m_strain[ptId];
        }

        template <typename Real>
        Eigen::Matrix<Real, 3, 3> CSphStructuralSolver<Real>::getStress(const Integer ptId) const
        {
            if (ptId < 0 || ptId >= m_nParticles)
                return Eigen::Matrix<Real, 3, 3>::Zero();

            return m_stress[ptId];
        }

        template <typename Real>
        void CSphStructuralSolver<Real>::reset()
        {
            m_points.clear();
            m_displacement.clear();
            m_strain.clear();
            m_stress.clear();
            m_externalForce.clear();
            m_isFixed.clear();
            m_fixedDisplacement.clear();
            m_nParticles = 0;
        }
    }
}
