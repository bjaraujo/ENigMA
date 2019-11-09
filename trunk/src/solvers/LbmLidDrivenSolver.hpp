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

#include "CmnTypes.hpp"

namespace ENigMA {

namespace lbm {

    template <typename Real, Integer nDim, Integer nSpeeds>
    class CLbmLidDrivenSolver {
    };

    template <typename Real>
    class CLbmLidDrivenSolver<Real, 2, 9> {
    private:
        Integer m_nx, m_ny, m_nz;

        std::vector<std::vector<Integer>> m_e;
        std::vector<Real> m_w;
        std::vector<Real> m_op;

        std::vector<std::vector<Integer>> m_B;
        std::vector<std::vector<Real>> m_rho;

        std::vector<std::vector<std::vector<Real>>> m_u, m_u0;
        std::vector<std::vector<std::vector<Real>>> m_f;
        std::vector<std::vector<std::vector<Real>>> m_F;

        Real feq2(Integer k, Real rho, std::vector<Real>& u);

        void resizeMatrix(std::vector<std::vector<Integer>>& mat);
        void resizeMatrix(std::vector<std::vector<Real>>& mat);
        void resizeMatrix(std::vector<std::vector<std::vector<Real>>>& mat, Integer n);

    public:
        CLbmLidDrivenSolver(Integer nx, Integer ny, Integer nz);
        virtual ~CLbmLidDrivenSolver();

        void setDensity(Integer i, Integer j, Real aValue);

        void setVelocity(Integer i, Integer j, Integer d, Real aValue);
        Real getVelocity(Integer i, Integer j, Integer d);

        void setBoundary(Integer i, Integer j, Integer aValue);
        Integer getBoundary(Integer i, Integer j);

        void init();
        void evolve(Real tau_f);
    };
}
}

#include "LbmLidDrivenSolver_Imp.hpp"
