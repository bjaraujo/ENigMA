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

namespace ENigMA
{
    namespace fvm
    {
        template <typename Real>
        CFvmTemperatureSolver<Real>::CFvmTemperatureSolver(CFvmMesh<Real>& aFvmMesh)
            : CFvmPisoSolver<Real>::CFvmPisoSolver(aFvmMesh)
            , m_bthcond(1.0)
        {
            for (Integer i = 0; i < CFvmPisoSolver<Real>::m_fvmMesh.nbControlVolumes(); ++i)
            {
                Integer aControlVolumeId = CFvmPisoSolver<Real>::m_fvmMesh.controlVolumeId(i);

                m_T0[aControlVolumeId] = m_T[aControlVolumeId] = 0.0;
            }

            for (Integer i = 0; i < CFvmPisoSolver<Real>::m_fvmMesh.nbFaces(); ++i)
            {
                Integer aFaceId = CFvmPisoSolver<Real>::m_fvmMesh.faceId(i);

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
            CFvmPisoSolver<Real>::storePreviousQuantities();

            m_T0 = m_T;
        }

        template <typename Real>
        void CFvmTemperatureSolver<Real>::setMaterialProperties(const Real aDensity, const Real aViscosity, const Real aThermalConductivity, const Real aSpecificHeat)
        {
            for (Integer i = 0; i < CFvmPisoSolver<Real>::m_fvmMesh.nbControlVolumes(); ++i)
            {
                Integer aControlVolumeId = CFvmPisoSolver<Real>::m_fvmMesh.controlVolumeId(i);

                CFvmPisoSolver<Real>::m_dens[aControlVolumeId] = aDensity;
                CFvmPisoSolver<Real>::m_visc[aControlVolumeId] = aViscosity;
                m_thcond[aControlVolumeId] = aThermalConductivity;
                m_spheat[aControlVolumeId] = aSpecificHeat;
            }
        }

        template <typename Real>
        void CFvmTemperatureSolver<Real>::setBoundaryTemperature(const std::vector<Integer>& sFaceIds, const EBoundaryType sFaceType, const Real T)
        {
            for (Integer i = 0; i < static_cast<Integer>(sFaceIds.size()); ++i)
            {
                Integer aFaceId = sFaceIds[i];

                m_Tf[aFaceId] = T;

                CFvmPisoSolver<Real>::m_fvmMesh.face(aFaceId).setBoundaryType(CFvmBoundaryType(sFaceType));
            }
        }

        template <typename Real>
        void CFvmTemperatureSolver<Real>::calculateTemperatureField()
        {
            Eigen::SparseMatrix<Real> A;
            Eigen::Matrix<Real, Eigen::Dynamic, 1> b;

            A.resize(CFvmPisoSolver<Real>::m_fvmMesh.nbControlVolumes(), CFvmPisoSolver<Real>::m_fvmMesh.nbControlVolumes());
            A.reserve(CFvmPisoSolver<Real>::m_fvmMesh.nbControlVolumes());

            b.resize(CFvmPisoSolver<Real>::m_fvmMesh.nbControlVolumes());
            b.setZero();

            // Assemble temperature matrix
            for (Integer i = 0; i < CFvmPisoSolver<Real>::m_fvmMesh.nbControlVolumes(); ++i)
            {
                Integer aControlVolumeId = CFvmPisoSolver<Real>::m_fvmMesh.controlVolumeId(i);

                Real volume = CFvmPisoSolver<Real>::m_fvmMesh.controlVolume(aControlVolumeId).originalVolume();

                Integer anIndexP = CFvmPisoSolver<Real>::m_mapIdToIndex.at(aControlVolumeId);

                Real dens = CFvmPisoSolver<Real>::m_dens.at(aControlVolumeId);
                Real visc = CFvmPisoSolver<Real>::m_visc.at(aControlVolumeId);
                Real thcond = m_thcond.at(aControlVolumeId);
                Real spheat = m_spheat.at(aControlVolumeId);

                for (Integer j = 0; j < CFvmPisoSolver<Real>::m_fvmMesh.controlVolume(aControlVolumeId).nbFaces(); ++j)
                {
                    Integer aFaceId = CFvmPisoSolver<Real>::m_fvmMesh.controlVolume(aControlVolumeId).faceId(j);

                    Real area = CFvmPisoSolver<Real>::m_fvmMesh.controlVolume(aControlVolumeId).faceArea(aFaceId);
                    Real dist = CFvmPisoSolver<Real>::m_fvmMesh.controlVolume(aControlVolumeId).faceDist(aFaceId);

                    Real flux;

                    if (CFvmPisoSolver<Real>::m_fvmMesh.face(aFaceId).controlVolumeId() == aControlVolumeId)
                        flux = +CFvmPisoSolver<Real>::m_flux.at(aFaceId);
                    else
                        flux = -CFvmPisoSolver<Real>::m_flux.at(aFaceId);

                    if (CFvmPisoSolver<Real>::m_fvmMesh.face(aFaceId).hasPair())
                    {
                        Integer aNeighborId = CFvmPisoSolver<Real>::m_fvmMesh.face(aFaceId).neighborId(aControlVolumeId);

                        Integer anIndexN = CFvmPisoSolver<Real>::m_mapIdToIndex.at(aNeighborId);

                        if (CFvmPisoSolver<Real>::m_fvmMesh.controlVolume(aNeighborId).containsFace(aFaceId))
                        {
                            area += CFvmPisoSolver<Real>::m_fvmMesh.controlVolume(aNeighborId).faceArea(aFaceId);
                            area *= 0.5;
                            dist += CFvmPisoSolver<Real>::m_fvmMesh.controlVolume(aNeighborId).faceDist(aFaceId);

                            dens += CFvmPisoSolver<Real>::m_dens[aNeighborId];
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
                    }
                    else
                    {
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

                if (CFvmPisoSolver<Real>::m_dt > 0)
                {
                    // Unsteady term - Euler
                    A.coeffRef(anIndexP, anIndexP) += volume * spheat * dens / CFvmPisoSolver<Real>::m_dt;
                    b[anIndexP] += volume * spheat * dens / CFvmPisoSolver<Real>::m_dt * m_T0.at(aControlVolumeId);
                }
            }

            A.finalize();

            Eigen::BiCGSTAB<Eigen::SparseMatrix<Real>> solver;
            solver.compute(A);

            Eigen::Matrix<Real, Eigen::Dynamic, 1> T = solver.solve(b);

            for (int i = 0; i < T.rows(); ++i)
            {
                Integer aControlVolumeId = CFvmPisoSolver<Real>::m_mapIndexToId[i];

                m_T[aControlVolumeId] = T[i];
            }
        }

        template <typename Real>
        void CFvmTemperatureSolver<Real>::iterate(const Real dt, const bool bInit)
        {
            CFvmPisoSolver<Real>::iterate(dt, bInit);

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

            for (Integer i = 0; i < CFvmPisoSolver<Real>::m_fvmMesh.nbControlVolumes(); ++i)
            {
                Integer aControlVolumeId = CFvmPisoSolver<Real>::m_fvmMesh.controlVolumeId(i);

                ru += (CFvmPisoSolver<Real>::m_u.at(aControlVolumeId) - CFvmPisoSolver<Real>::m_u0.at(aControlVolumeId)) * (CFvmPisoSolver<Real>::m_u.at(aControlVolumeId) - CFvmPisoSolver<Real>::m_u0.at(aControlVolumeId));
                rv += (CFvmPisoSolver<Real>::m_v.at(aControlVolumeId) - CFvmPisoSolver<Real>::m_v0.at(aControlVolumeId)) * (CFvmPisoSolver<Real>::m_v.at(aControlVolumeId) - CFvmPisoSolver<Real>::m_v0.at(aControlVolumeId));
                rw += (CFvmPisoSolver<Real>::m_w.at(aControlVolumeId) - CFvmPisoSolver<Real>::m_w0.at(aControlVolumeId)) * (CFvmPisoSolver<Real>::m_w.at(aControlVolumeId) - CFvmPisoSolver<Real>::m_w0.at(aControlVolumeId));
                rp += (CFvmPisoSolver<Real>::m_p.at(aControlVolumeId) - CFvmPisoSolver<Real>::m_p0.at(aControlVolumeId)) * (CFvmPisoSolver<Real>::m_p.at(aControlVolumeId) - CFvmPisoSolver<Real>::m_p0.at(aControlVolumeId));
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
