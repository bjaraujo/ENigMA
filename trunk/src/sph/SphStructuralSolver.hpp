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
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "CmnTypes.hpp"
#include "GeoCoordinate.hpp"
#include "GeoVector.hpp"
#include "GeoHashGrid.hpp"
#include "SphKernel.hpp"

namespace ENigMA
{
    namespace sph
    {
        template <typename Real>
        class CSphStructuralSolver
        {
        private:
            std::vector<ENigMA::geometry::CGeoCoordinate<Real>> m_points;           // Material points
            std::vector<ENigMA::geometry::CGeoVector<Real>> m_displacement;         // u displacement vector
            std::vector<Eigen::Matrix<Real, 3, 3>> m_strain;                        // ε strain tensor
            std::vector<Eigen::Matrix<Real, 3, 3>> m_stress;                        // σ stress tensor
            std::vector<Eigen::Matrix<Real, 3, 1>> m_externalForce;                 // f external force
            std::vector<bool> m_isFixed;                                             // Whether point is fixed
            std::vector<ENigMA::geometry::CGeoVector<Real>> m_fixedDisplacement;     // Fixed displacement BC

            ENigMA::geometry::CGeoHashGrid<Real> m_hashGrid;
            std::vector<Integer> m_hashGridData;

            CSphKernel<Real>* m_kernel;

            Real m_E;          // Young's modulus
            Real m_nu;         // Poisson's ratio
            Real m_h;          // Support radius / smoothing length
            Real m_lambda;     // Lamé parameter: λ
            Real m_mu;         // Lamé parameter: μ (shear modulus)

            Integer m_nParticles;
            Eigen::SparseMatrix<Real, Eigen::RowMajor> m_stiffnessMatrix;
            Eigen::Matrix<Real, Eigen::Dynamic, 1> m_residual;

            void calculateStrain();
            void calculateStress();
            Eigen::Matrix<Real, Eigen::Dynamic, 1> assembleResidual();
            void buildStiffnessMatrix();
            void applyBoundaryConditions(Eigen::SparseMatrix<Real, Eigen::RowMajor>& K,
                                          Eigen::Matrix<Real, Eigen::Dynamic, 1>& F);

        public:
            CSphStructuralSolver();
            ~CSphStructuralSolver();

            // Setup methods
            void setKernel(CSphKernel<Real>* kernel);
            void setSmoothingLength(const Real h);
            void setMaterialProperties(const Real E, const Real nu);

            // Point management
            void addParticle(const ENigMA::geometry::CGeoCoordinate<Real>& p);
            void addParticles(const std::vector<ENigMA::geometry::CGeoCoordinate<Real>>& points);
            Integer nbParticles() const;

            // Boundary conditions
            void setFixedBoundaryCondition(const Integer ptId, const ENigMA::geometry::CGeoVector<Real>& u);
            void setForce(const Integer ptId, const ENigMA::geometry::CGeoVector<Real>& f);

            // Solver
            void solve(const Integer maxIterations = 20, const Real tolerance = 1e-6);

            // Query results
            ENigMA::geometry::CGeoVector<Real> getDisplacement(const Integer ptId) const;
            Eigen::Matrix<Real, 3, 3> getStrain(const Integer ptId) const;
            Eigen::Matrix<Real, 3, 3> getStress(const Integer ptId) const;

            // Utilities
            void reset();
        };
    }
}

#include "SphStructuralSolver_Imp.hpp"
