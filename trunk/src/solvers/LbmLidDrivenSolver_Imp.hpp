// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

namespace ENigMA {

namespace lbm {

    template <typename Real>
    CLbmLidDrivenSolver<Real, 2, 9>::CLbmLidDrivenSolver(Integer nx, Integer ny, Integer nz)
        : m_nx(nx)
        , m_ny(ny)
        , m_nz(nz)
    {

        m_e.resize(9);
        for (Integer i = 0; i < 9; ++i)
            m_e[i].resize(2);

        m_w.resize(9);
        m_op.resize(9);

        m_op[0] = 0;
        m_op[1] = 3;
        m_op[2] = 4;
        m_op[3] = 1;
        m_op[4] = 2;
        m_op[5] = 7;
        m_op[6] = 8;
        m_op[7] = 5;
        m_op[8] = 6;

        m_e[0][0] = +0;
        m_e[0][1] = +0;
        m_w[0] = 4.0 / 9.0;

        m_e[1][0] = +1;
        m_e[1][1] = +0;
        m_w[1] = 1.0 / 9.0;

        m_e[2][0] = +0;
        m_e[2][1] = +1;
        m_w[2] = 1.0 / 9.0;

        m_e[3][0] = -1;
        m_e[3][1] = +0;
        m_w[3] = 1.0 / 9.0;

        m_e[4][0] = +0;
        m_e[4][1] = -1;
        m_w[4] = 1.0 / 9.0;

        m_e[5][0] = +1;
        m_e[5][1] = +1;
        m_w[5] = 1.0 / 36.0;

        m_e[6][0] = -1;
        m_e[6][1] = +1;
        m_w[6] = 1.0 / 36.0;

        m_e[7][0] = -1;
        m_e[7][1] = -1;
        m_w[7] = 1.0 / 36.0;

        m_e[8][0] = +1;
        m_e[8][1] = -1;
        m_w[8] = 1.0 / 36.0;

        this->resizeMatrix(m_u, 2);
        this->resizeMatrix(m_u0, 2);
        this->resizeMatrix(m_rho);
        this->resizeMatrix(m_B);
        this->resizeMatrix(m_f, 9);
        this->resizeMatrix(m_F, 9);
    }

    template <typename Real>
    CLbmLidDrivenSolver<Real, 2, 9>::~CLbmLidDrivenSolver()
    {
    }

    template <typename Real>
    void CLbmLidDrivenSolver<Real, 2, 9>::resizeMatrix(std::vector<std::vector<Integer>>& mat)
    {

        mat.resize(m_nx + 1);

        for (Integer i = 0; i <= m_nx; ++i) {

            mat[i].resize(m_ny + 1);
        }
    }

    template <typename Real>
    void CLbmLidDrivenSolver<Real, 2, 9>::resizeMatrix(std::vector<std::vector<Real>>& mat)
    {

        mat.resize(m_nx + 1);

        for (Integer i = 0; i <= m_nx; ++i) {

            mat[i].resize(m_ny + 1);
        }
    }

    template <typename Real>
    void CLbmLidDrivenSolver<Real, 2, 9>::resizeMatrix(std::vector<std::vector<std::vector<Real>>>& mat, Integer n)
    {

        mat.resize(m_nx + 1);

        for (Integer i = 0; i <= m_nx; ++i) {

            mat[i].resize(m_ny + 1);

            for (Integer j = 0; j <= m_ny; ++j)
                mat[i][j].resize(n);
        }
    }

    template <typename Real>
    Real CLbmLidDrivenSolver<Real, 2, 9>::feq2(Integer k, Real rho, std::vector<Real>& u)
    {

        Real eu = m_e[k][0] * u[0] + m_e[k][1] * u[1];
        Real uv = u[0] * u[0] + u[1] * u[1];

        return m_w[k] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * uv);
    }

    template <typename Real>
    void CLbmLidDrivenSolver<Real, 2, 9>::evolve(Real tau_f)
    {

        // Evolution: streaming and collision
        for (Integer i = 1; i < m_nx; ++i) {
            for (Integer j = 1; j < m_ny; ++j) {
                for (Integer k = 0; k < 9; ++k) {

                    int ip = i - m_e[k][0];
                    int jp = j - m_e[k][1];
                    m_F[i][j][k] = m_f[ip][jp][k] + (this->feq2(k, m_rho[ip][jp], m_u[ip][jp]) - m_f[ip][jp][k]) / tau_f;
                }
            }
        }

        std::vector<Real> Ft;
        Ft.resize(9);

        for (Integer i = 0; i <= m_nx; ++i) {
            for (Integer j = 0; j <= m_ny; ++j) {

                if (m_B[i][j] == 1) {

                    for (Integer k = 0; k < 9; ++k)
                        Ft[k] = m_F[i][j][k];

                    for (Integer k = 0; k < 9; ++k)
                        m_F[i][j][k] = Ft[m_op[k]];
                }
            }
        }

        // Calculation of macroscopic quantities
        for (Integer i = 1; i < m_nx; ++i) {
            for (Integer j = 1; j < m_ny; ++j) {

                m_u0[i][j][0] = m_u[i][j][0];
                m_u0[i][j][1] = m_u[i][j][1];

                m_rho[i][j] = 0;
                m_u[i][j][0] = 0;
                m_u[i][j][1] = 0;

                for (Integer k = 0; k < 9; ++k) {
                    m_f[i][j][k] = m_F[i][j][k];
                    m_rho[i][j] += m_f[i][j][k];
                    m_u[i][j][0] += m_e[k][0] * m_f[i][j][k];
                    m_u[i][j][1] += m_e[k][1] * m_f[i][j][k];
                }

                m_u[i][j][0] /= m_rho[i][j];
                m_u[i][j][1] /= m_rho[i][j];
            }
        }

        for (Integer i = 0; i <= m_nx; ++i) // Upper boundary
        {
            for (Integer k = 0; k < 9; ++k) {
                m_rho[i][m_ny] = m_rho[i][m_ny - 1];
                m_f[i][m_ny][k] = this->feq2(k, m_rho[i][m_ny], m_u[i][m_ny]) + m_f[i][m_ny - 1][k] - this->feq2(k, m_rho[i][m_ny - 1], m_u[i][m_ny - 1]);
            }
        }
    }

    template <typename Real>
    void CLbmLidDrivenSolver<Real, 2, 9>::setDensity(Integer i, Integer j, Real aValue)
    {

        m_rho[i][j] = aValue;
    }

    template <typename Real>
    void CLbmLidDrivenSolver<Real, 2, 9>::setVelocity(Integer i, Integer j, Integer d, Real aValue)
    {

        m_u[i][j][d] = aValue;
    }

    template <typename Real>
    Real CLbmLidDrivenSolver<Real, 2, 9>::getVelocity(Integer i, Integer j, Integer d)
    {

        return m_u[i][j][d];
    }

    template <typename Real>
    void CLbmLidDrivenSolver<Real, 2, 9>::setBoundary(Integer i, Integer j, Integer aValue)
    {

        m_B[i][j] = aValue;
    }

    template <typename Real>
    Integer CLbmLidDrivenSolver<Real, 2, 9>::getBoundary(Integer i, Integer j)
    {

        return m_B[i][j];
    }

    template <typename Real>
    void CLbmLidDrivenSolver<Real, 2, 9>::init()
    {

        for (Integer i = 0; i <= m_nx; ++i) {

            for (Integer j = 0; j <= m_ny; ++j) {

                m_u[i][j][0] = 0;
                m_u[i][j][1] = 0;
                m_rho[i][j] = 1.0;

                for (int k = 0; k < 9; ++k)
                    m_f[i][j][k] = this->feq2(k, m_rho[i][j], m_u[i][j]);
            }
        }
    }
}
}
