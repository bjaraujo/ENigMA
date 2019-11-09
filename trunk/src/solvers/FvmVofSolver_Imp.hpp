// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include <Eigen/Sparse>

namespace ENigMA {

namespace fvm {

    template <typename Real>
    CFvmVofSolver<Real>::CFvmVofSolver(CFvmMesh<Real>& aFvmMesh)
        : CFvmPisoSolver(aFvmMesh)
    {

        for (Integer i = 0; i < m_fvmMesh.nbControlVolumes(); ++i) {

            Integer aControlVolumeId = m_fvmMesh.controlVolumeId(i);

            m_s0[aControlVolumeId] = m_s[aControlVolumeId] = 0.0;
        }

        for (Integer i = 0; i < m_fvmMesh.nbFaces(); ++i) {

            Integer aFaceId = m_fvmMesh.faceId(i);

            m_sf[aFaceId] = 0.0;
        }
    }

    template <typename Real>
    CFvmVofSolver<Real>::~CFvmVofSolver()
    {
    }

    template <typename Real>
    void CFvmVofSolver<Real>::storePreviousQuantities()
    {

        CFvmPisoSolver::storePreviousQuantities();

        m_s0 = m_s;
    }

    template <typename Real>
    void CFvmVofSolver<Real>::setMaterialProperties(const Real aDensity0, const Real aViscosity0, const Real aDensity1, const Real aViscosity1)
    {

        m_dens0 = aDensity0;
        m_visc0 = aViscosity0;

        m_dens1 = aDensity1;
        m_visc1 = aViscosity1;

        this->updateProperties();
    }

    template <typename Real>
    void CFvmVofSolver<Real>::setInitialGamma(const Integer aControlVolumeId, const Real s)
    {

        m_s[aControlVolumeId] = s;
    }

    template <typename Real>
    void CFvmVofSolver<Real>::setBoundaryGamma(const Integer aFaceId, const EBoundaryType sFaceType, const Real s)
    {

        m_sf[aFaceId] = s;
        m_fvmMesh.face(aFaceId).setBoundaryType(sFaceType);
    }

    template <typename Real>
    void CFvmVofSolver<Real>::setBoundaryGamma(const std::vector<Integer>& sFaceIds, const EBoundaryType sFaceType, const Real s)
    {

        for (Integer i = 0; i < static_cast<Integer>(sFaceIds.size()); ++i) {

            Integer aFaceId = sFaceIds[i];

            m_sf[aFaceId] = s;

            m_fvmMesh.face(aFaceId).setBoundaryType(sFaceType);
        }
    }

    template <typename Real>
    Real CFvmVofSolver<Real>::calculateCourant(double dt, bool bInterface)
    {

        Real maxCp = 0.0;

        for (Integer i = 0; i < m_fvmMesh.nbControlVolumes(); ++i) {

            Integer aControlVolumeId = m_fvmMesh.controlVolumeId(i);

            Real volume = m_fvmMesh.controlVolume(aControlVolumeId).originalVolume();

            Real cs = 1.0;

            if (bInterface) {
                Real s = std::min(std::max(m_s.at(aControlVolumeId), 0.0), 1.0);
                cs = (1.0 - s) * (1.0 - s) * s * s * 16.0;
            }

            Real Cpp = 0.0;

            for (Integer j = 0; j < m_fvmMesh.controlVolume(aControlVolumeId).nbFaces(); ++j) {

                Integer aFaceId = m_fvmMesh.controlVolume(aControlVolumeId).faceId(j);

                Real flux;

                if (m_fvmMesh.face(aFaceId).controlVolumeId() == aControlVolumeId)
                    flux = +m_flux.at(aFaceId);
                else
                    flux = -m_flux.at(aFaceId);

                Real area = m_fvmMesh.controlVolume(aControlVolumeId).faceArea(aFaceId);

                Real Cj = std::max(-flux * area * dt / volume, 0.0);
                Cpp += cs * Cj;
            }

            maxCp = std::max(maxCp, Cpp);

            m_Co[aControlVolumeId] = Cpp;
        }

        return maxCp;
    }

    template <typename Real>
    void CFvmVofSolver<Real>::predictBeta(const Real aTolerance)
    {

        Integer anAcceptor, aDonor;
        CGeoVector<Real> grads;
        Real dot;
        Real l1, l2;

        for (Integer i = 0; i < m_fvmMesh.nbFaces(); ++i) {

            Real betaj = 1.0;

            Integer aFaceId = m_fvmMesh.faceId(i);

            Integer aControlVolumeId = m_fvmMesh.face(aFaceId).controlVolumeId();

            if (m_fvmMesh.face(aFaceId).hasPair()) {

                Integer aNeighborId = m_fvmMesh.face(aFaceId).neighborId(aControlVolumeId);

                Real area = m_fvmMesh.controlVolume(aControlVolumeId).faceArea(aFaceId);
                Real dist = m_fvmMesh.controlVolume(aControlVolumeId).faceDist(aFaceId);

                Real flux;

                if (m_fvmMesh.face(aFaceId).controlVolumeId() == aControlVolumeId)
                    flux = +m_flux.at(aFaceId);
                else
                    flux = -m_flux.at(aFaceId);

                if (flux < 0.0) {

                    anAcceptor = aNeighborId;
                    aDonor = aControlVolumeId;

                    grads = this->gradient(m_s, m_sf, aDonor);

                    dot = grads.dot(m_fvmMesh.face(aFaceId).normal()) * dist;

                    l1 = grads.norm();
                    l2 = dist;

                } else {

                    anAcceptor = aControlVolumeId;
                    aDonor = aNeighborId;

                    grads = this->gradient(m_s, m_sf, aDonor);

                    dot = grads.dot(m_fvmMesh.face(aFaceId).normal()) * dist;

                    l1 = grads.norm();
                    l2 = dist;
                }

                Real su = std::min(std::max(m_s[anAcceptor] - 2 * dot, 0.0), 1.0);

                Real Cod = std::min(m_Co[aDonor], 1.0);

                if (fabs(m_s[anAcceptor] - su) > aTolerance) {

                    Real sdn = (m_s[aDonor] - su) / (m_s[anAcceptor] - su);

                    Real sjnCBC = 0.0;

                    if (sdn >= 0.0 && sdn <= 1.0 && fabs(Cod) > aTolerance)
                        sjnCBC = std::min(1.0, sdn / Cod);
                    else
                        sjnCBC = sdn;

                    Real sjnUQ = 0.0;

                    if (sdn >= 0.0 && sdn <= 1.0)
                        sjnUQ = std::min((8.0 * Cod * sdn + (1.0 - Cod) * (6.0 * sdn + 3.0)) / 8.0, sjnCBC);
                    else
                        sjnUQ = sdn;

                    Real ang = 0;

                    if (fabs(l1 * l2) > aTolerance)
                        ang = fabs(dot / (l1 * l2));
                    else
                        ang = fabs(dot / aTolerance);

                    if (ang > 1.0)
                        ang = 1.0;

                    Real tetaj = acos(ang);

                    Real kq = 2.0;

                    Real qj = std::min(kq * 0.5 * (cos(2 * tetaj) + 1.0), 1.0);

                    Real sjn = qj * sjnCBC + (1 - qj) * sjnUQ;

                    if (fabs(1.0 - sdn) > aTolerance)
                        betaj = std::min(std::max((sjn - sdn) / (1.0 - sdn), 0.0), 1.0);
                }
            }

            m_betaf[aFaceId] = betaj;
        }
    }

    template <typename Real>
    void CFvmVofSolver<Real>::correctBeta(double dt, const Real aTolerance)
    {

        Integer anAcceptor, aDonor;

        for (Integer i = 0; i < m_fvmMesh.nbFaces(); ++i) {

            Integer aFaceId = m_fvmMesh.faceId(i);

            Integer aControlVolumeId = m_fvmMesh.face(aFaceId).controlVolumeId();

            Real volume = m_fvmMesh.controlVolume(aControlVolumeId).originalVolume();

            Real betaj = m_betaf.at(aFaceId);

            if (betaj < 1E-2)
                continue;

            Real cbetaj = 0.0;

            if (m_fvmMesh.face(aFaceId).hasPair()) {

                Integer aNeighborId = m_fvmMesh.face(aFaceId).neighborId(aControlVolumeId);

                Real area = m_fvmMesh.controlVolume(aControlVolumeId).faceArea(aFaceId);

                Real flux;

                if (m_fvmMesh.face(aFaceId).controlVolumeId() == aControlVolumeId)
                    flux = +m_flux.at(aFaceId);
                else
                    flux = -m_flux.at(aFaceId);

                if (flux != 0.0) {

                    if (flux < 0.0) {
                        anAcceptor = aNeighborId;
                        aDonor = aControlVolumeId;
                    } else {
                        anAcceptor = aControlVolumeId;
                        aDonor = aNeighborId;
                    }

                    Real Cj = std::min(std::max(-flux * area * dt / volume, 0.0), 1.0);

                    Real ds = 0.5 * (m_s0[anAcceptor] + m_s[anAcceptor]) - 0.5 * (m_s0[aDonor] + m_s[aDonor]);

                    if (m_s[aDonor] < 0.0) {
                        Real Em = std::max(-m_s[aDonor], 0.0);

                        // Donor value < 0.0 Ex: sd = -0.1 -> Em = +0.1
                        if (Em > aTolerance && Cj > aTolerance) {
                            if (ds > Em) {
                                cbetaj = Em * (2 + Cj - 2 * Cj * betaj) / (2 * Cj * (ds - Em));

                                cbetaj = std::min(cbetaj, betaj);
                            }
                        }
                    }

                    if (m_s[aDonor] > 1.0) {

                        Real Ep = std::max(m_s[aDonor] - 1.0, 0.0);

                        // Donor value > 1.0 Ex: sd = 1.1 -> Ep = +0.1
                        if (Ep > std::numeric_limits<Real>::epsilon() && Cj > std::numeric_limits<Real>::epsilon()) {
                            if (ds < -Ep) {
                                cbetaj = Ep * (2 + Cj - 2 * Cj * betaj) / (2 * Cj * (-ds - Ep));

                                cbetaj = std::min(cbetaj, betaj);
                            }
                        }
                    }
                }
            }

            betaj -= cbetaj;

            betaj = std::max(betaj, 0.0);

            m_betaf[aFaceId] = betaj;
        }
    }

    template <typename Real>
    void CFvmVofSolver<Real>::updateProperties()
    {

        for (Integer i = 0; i < m_fvmMesh.nbControlVolumes(); ++i) {

            Integer aControlVolumeId = m_fvmMesh.controlVolumeId(i);

            m_dens[aControlVolumeId] = m_dens0 * (1.0 - m_s.at(aControlVolumeId)) + m_dens1 * m_s.at(aControlVolumeId);
            m_visc[aControlVolumeId] = m_visc0 * (1.0 - m_s.at(aControlVolumeId)) + m_visc1 * m_s.at(aControlVolumeId);
        }
    }

    template <typename Real>
    void CFvmVofSolver<Real>::calculateGammaField()
    {

        Eigen::SparseMatrix<Real> A;
        Eigen::Matrix<Real, Eigen::Dynamic, 1> b;

        A.resize(m_fvmMesh.nbControlVolumes(), m_fvmMesh.nbControlVolumes());
        A.reserve(m_fvmMesh.nbControlVolumes());

        b.resize(m_fvmMesh.nbControlVolumes());
        b.setZero();

        // Assemble gamma matrix
        for (Integer i = 0; i < m_fvmMesh.nbControlVolumes(); ++i) {

            Integer aControlVolumeId = m_fvmMesh.controlVolumeId(i);

            Real volume = m_fvmMesh.controlVolume(aControlVolumeId).originalVolume();

            Integer anIndexP = m_mapIdToIndex.at(aControlVolumeId);

            for (Integer j = 0; j < m_fvmMesh.controlVolume(aControlVolumeId).nbFaces(); ++j) {

                Integer aFaceId = m_fvmMesh.controlVolume(aControlVolumeId).faceId(j);

                Real area = m_fvmMesh.controlVolume(aControlVolumeId).faceArea(aFaceId);
                Real dist = m_fvmMesh.controlVolume(aControlVolumeId).faceDist(aFaceId);

                Real flux;

                if (m_fvmMesh.face(aFaceId).controlVolumeId() == aControlVolumeId)
                    flux = +m_flux.at(aFaceId);
                else
                    flux = -m_flux.at(aFaceId);

                if (m_fvmMesh.face(aFaceId).hasPair()) {

                    Integer aNeighborId = m_fvmMesh.face(aFaceId).neighborId(aControlVolumeId);

                    Integer anIndexN = m_mapIdToIndex.at(aNeighborId);

                    if (m_fvmMesh.controlVolume(aNeighborId).containsFace(aFaceId)) {

                        Real beta = 0.0;

                        if (flux < 0.0)
                            beta = m_betaf.at(aFaceId);
                        else
                            beta = 1.0 - m_betaf.at(aFaceId);

                        // Convection
                        A.coeffRef(anIndexP, anIndexP) += 0.5 * (1.0 - beta) * flux;
                        A.coeffRef(anIndexP, anIndexN) += 0.5 * beta * flux;

                        b.at(anIndexP) += -0.5 * (1.0 - beta) * flux * m_s.at(aControlVolumeId);
                        b.at(anIndexP) += -0.5 * beta * flux * m_s.at(aNeighborId);
                    }

                } else {

                    // Convection
                    b.at(anIndexP) += -1.0 * flux * m_sf.at(aFaceId);
                }
            }

            if (m_dt > 0) {

                // Unsteady term - Euler
                A.coeffRef(anIndexP, anIndexP) += volume / m_dt;
                b.at(anIndexP) += volume / m_dt * m_s0.at(aControlVolumeId);
            }
        }

        A.finalize();

        Eigen::BiCGSTAB<Eigen::SparseMatrix<Real>> solver;
        solver.compute(A);

        Eigen::Matrix<Real, Eigen::Dynamic, 1> s = solver.solve(b);

        for (int i = 0; i < s.rows(); ++i) {

            Integer aControlVolumeId = m_mapIndexToId.at(i);

            m_s[aControlVolumeId] = std::min(std::max(s.at(i), 0.0), 1.0);
        }
    }

    template <typename Real>
    void CFvmVofSolver<Real>::iterate(const Real dt, const bool bInit)
    {

        this->updateProperties();

        CFvmPisoSolver::iterate(dt, bInit);

        Real maxCp = this->calculateCourant(dt, true);

        Integer n = std::min((Integer)(maxCp * 2) + 1, 100);

        for (Integer j = 0; j < n; ++j) {

            // Currant number
            //maxCp = this->calculateCourant(dt / n, false);

            this->predictBeta(1E-16);

            // CICSAM corrections
            Integer m = 2;

            for (Integer l = 0; l < m; ++l) {

                this->calculateCourant(dt / m, false);

                this->calculateGammaField();

                this->correctBeta(dt / m, 1E-16);
            }
        }
    }

    template <typename Real>
    void CFvmVofSolver<Real>::residual(Real& ru, Real& rv, Real& rw, Real& rp, Real& rs)
    {

        ru = 0;
        rv = 0;
        rw = 0;
        rp = 0;
        rs = 0;

        for (Integer i = 0; i < m_fvmMesh.nbControlVolumes(); ++i) {

            Integer aControlVolumeId = m_fvmMesh.controlVolumeId(i);

            ru += (m_u.at(aControlVolumeId) - m_u0.at(aControlVolumeId)) * (m_u.at(aControlVolumeId) - m_u0.at(aControlVolumeId));
            rv += (m_v.at(aControlVolumeId) - m_v0.at(aControlVolumeId)) * (m_v.at(aControlVolumeId) - m_v0.at(aControlVolumeId));
            rw += (m_w.at(aControlVolumeId) - m_w0.at(aControlVolumeId)) * (m_w.at(aControlVolumeId) - m_w0.at(aControlVolumeId));
            rp += (m_p.at(aControlVolumeId) - m_p0.at(aControlVolumeId)) * (m_p.at(aControlVolumeId) - m_p0.at(aControlVolumeId));
            rs += (m_s.at(aControlVolumeId) - m_s0.at(aControlVolumeId)) * (m_s.at(aControlVolumeId) - m_s0.at(aControlVolumeId));
        }
    }

    template <typename Real>
    Real CFvmVofSolver<Real>::s(const Integer aControlVolumeId)
    {

        return m_s.at(aControlVolumeId);
    }

    template <typename Real>
    Real CFvmVofSolver<Real>::sf(const Integer aFaceId)
    {

        return m_sf.at(aFaceId);
    }

    template <typename Real>
    CGeoVector<Real> CFvmVofSolver<Real>::grads(const Integer aControlVolumeId)
    {

        return this->gradient(m_s, m_sf, aControlVolumeId);
    }
}
}
