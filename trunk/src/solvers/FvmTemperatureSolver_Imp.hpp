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
    CFvmTemperatureSolver<Real>::CFvmTemperatureSolver(CFvmMesh<Real>& aFvmMesh)
        : CFvmPisoSolver(aFvmMesh), 
        m_bthcond(1.0)
    {
        for (Integer i = 0; i < m_fvmMesh.nbControlVolumes(); ++i) {
            Integer aControlVolumeId = m_fvmMesh.controlVolumeId(i);

            m_T0[aControlVolumeId] = m_T[aControlVolumeId] = 0.0;
        }

        for (Integer i = 0; i < m_fvmMesh.nbFaces(); ++i) {
            Integer aFaceId = m_fvmMesh.faceId(i);

            m_Tf[aFaceId] = 0.0;
        }
    }

    template <typename Real>
    CFvmTemperatureSolver<Real>::~CFvmTemperatureSolver()
    {
    }

    template <typename Real>
    void CFvmTemperatureSolver<Real>::storePreviousQuantities()
    {
        CFvmPisoSolver::storePreviousQuantities();

        m_T0 = m_T;
    }

    template <typename Real>
    void CFvmTemperatureSolver<Real>::setMaterialProperties(const Real aDensity, const Real aViscosity, const Real aThermalConductivity, const Real aSpecificHeat)
    {
        for (Integer i = 0; i < m_fvmMesh.nbControlVolumes(); ++i) {
            Integer aControlVolumeId = m_fvmMesh.controlVolumeId(i);

            m_dens[aControlVolumeId] = aDensity;
            m_visc[aControlVolumeId] = aViscosity;
            m_thcond[aControlVolumeId] = aThermalConductivity;
            m_spheat[aControlVolumeId] = aSpecificHeat;
        }
    }

    template <typename Real>
    void CFvmTemperatureSolver<Real>::setBoundaryTemperature(const std::vector<Integer>& sFaceIds, const EBoundaryType sFaceType, const Real T)
    {
        for (Integer i = 0; i < static_cast<Integer>(sFaceIds.size()); ++i) {
            Integer aFaceId = sFaceIds[i];

            m_Tf[aFaceId] = T;

            m_fvmMesh.face(aFaceId).setBoundaryType(CFvmBoundaryType(sFaceType));
        }
    }

    template <typename Real>
    void CFvmTemperatureSolver<Real>::calculateTemperatureField()
    {
        Eigen::SparseMatrix<Real> A;
        Eigen::Matrix<Real, Eigen::Dynamic, 1> b;

        A.resize(m_fvmMesh.nbControlVolumes(), m_fvmMesh.nbControlVolumes());
        A.reserve(m_fvmMesh.nbControlVolumes());

        b.resize(m_fvmMesh.nbControlVolumes());
        b.setZero();

        // Assemble temperature matrix
        for (Integer i = 0; i < m_fvmMesh.nbControlVolumes(); ++i) {
            Integer aControlVolumeId = m_fvmMesh.controlVolumeId(i);

            Real volume = m_fvmMesh.controlVolume(aControlVolumeId).originalVolume();

            Integer anIndexP = m_mapIdToIndex.at(aControlVolumeId);

            Real dens = m_dens.at(aControlVolumeId);
            Real visc = m_visc.at(aControlVolumeId);
            Real thcond = m_thcond.at(aControlVolumeId);
            Real spheat = m_spheat.at(aControlVolumeId);

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
                        area += m_fvmMesh.controlVolume(aNeighborId).faceArea(aFaceId);
                        area *= 0.5;
                        dist += m_fvmMesh.controlVolume(aNeighborId).faceDist(aFaceId);

                        dens += m_dens[aNeighborId];
                        dens *= 0.5;
                        thcond += m_thcond[aNeighborId];
                        thcond *= 0.5;
                        spheat += m_spheat[aNeighborId];
                        spheat *= 0.5;

                        // Conduction
                        A.coeffRef(anIndexP, anIndexP) += thcond * area / dist;
                        A.coeffRef(anIndexP, anIndexN) += -thcond * area / dist;

                        // Convection
                        Real xsi;
                        if (flux > 0.0)
                            xsi = 0.0;
                        else
                            xsi = 1.0; // UPWIND

                        A.coeffRef(anIndexP, anIndexP) += (1.0 - xsi) * dens * spheat * flux;
                        A.coeffRef(anIndexP, anIndexN) += xsi * dens * spheat * flux;
                    }

                } else {
                    thcond += m_bthcond;
                    thcond *= 0.5;

                    // Conduction
                    A.coeffRef(anIndexP, anIndexP) += thcond * area / dist;
                    b[anIndexP] += thcond * area / dist / volume * m_Tf.at(aFaceId);

                    // Convection
                    Real xsi;
                    if (flux > 0.0)
                        xsi = 0.0;
                    else
                        xsi = 1.0; // UPWIND

                    A.coeffRef(anIndexP, anIndexP) += (1.0 - xsi) * dens * spheat * flux;
                    b[anIndexP] += -xsi * dens * spheat * flux * m_Tf.at(aFaceId);
                }
            }

            if (m_dt > 0) {
                // Unsteady term - Euler
                A.coeffRef(anIndexP, anIndexP) += volume * spheat * dens / m_dt;
                b[anIndexP] += volume * spheat * dens / m_dt * m_T0.at(aControlVolumeId);
            }
        }

        A.finalize();

        Eigen::BiCGSTAB<Eigen::SparseMatrix<Real>> solver;
        solver.compute(A);

        Eigen::Matrix<Real, Eigen::Dynamic, 1> T = solver.solve(b);

        for (int i = 0; i < T.rows(); ++i) {
            Integer aControlVolumeId = m_mapIndexToId[i];

            m_T[aControlVolumeId] = T[i];
        }
    }

    template <typename Real>
    void CFvmTemperatureSolver<Real>::iterate(const Real dt, const bool bInit)
    {
        CFvmPisoSolver::iterate(dt, bInit);

        this->calculateTemperatureField();
    }

    template <typename Real>
    void CFvmTemperatureSolver<Real>::residual(Real& ru, Real& rv, Real& rw, Real& rp, Real& rT)
    {
        ru = 0;
        rv = 0;
        rw = 0;
        rp = 0;
        rT = 0;

        for (Integer i = 0; i < m_fvmMesh.nbControlVolumes(); ++i) {
            Integer aControlVolumeId = m_fvmMesh.controlVolumeId(i);

            ru += (m_u.at(aControlVolumeId) - m_u0.at(aControlVolumeId)) * (m_u.at(aControlVolumeId) - m_u0.at(aControlVolumeId));
            rv += (m_v.at(aControlVolumeId) - m_v0.at(aControlVolumeId)) * (m_v.at(aControlVolumeId) - m_v0.at(aControlVolumeId));
            rw += (m_w.at(aControlVolumeId) - m_w0.at(aControlVolumeId)) * (m_w.at(aControlVolumeId) - m_w0.at(aControlVolumeId));
            rp += (m_p.at(aControlVolumeId) - m_p0.at(aControlVolumeId)) * (m_p.at(aControlVolumeId) - m_p0.at(aControlVolumeId));
            rT += (m_T.at(aControlVolumeId) - m_T0.at(aControlVolumeId)) * (m_T.at(aControlVolumeId) - m_T0.at(aControlVolumeId));
        }
    }

    template <typename Real>
    Real CFvmTemperatureSolver<Real>::T(const Integer aControlVolumeId)
    {
        return m_T.at(aControlVolumeId);
    }

    template <typename Real>
    Real CFvmTemperatureSolver<Real>::Tf(const Integer aFaceId)
    {
        return m_Tf.at(aFaceId);
    }

    template <typename Real>
    CGeoVector<Real> CFvmTemperatureSolver<Real>::gradT(const Integer aControlVolumeId)
    {
        return this->gradient(m_T, m_Tf, aControlVolumeId);
    }
}
}
