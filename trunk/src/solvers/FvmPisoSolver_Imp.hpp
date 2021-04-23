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
        CFvmPisoSolver<Real>::CFvmPisoSolver(CFvmMesh<Real>& aFvmMesh)
            : m_dt(0.0)
            , m_gx(0.0)
            , m_gy(0.0)
            , m_gz(0.0)
        {

            m_fvmMesh = aFvmMesh;

            for (Integer i = 0; i < m_fvmMesh.nbControlVolumes(); ++i)
            {

                Integer aControlVolumeId = m_fvmMesh.controlVolumeId(i);

                m_fvmMesh.controlVolume(aControlVolumeId).calculateVolume();
                m_fvmMesh.controlVolume(aControlVolumeId).calculateOriginalVolume();
                m_fvmMesh.controlVolume(aControlVolumeId).calculateCentroid();

                for (Integer j = 0; j < m_fvmMesh.controlVolume(aControlVolumeId).nbFaces(); ++j)
                {

                    Integer aFaceId = m_fvmMesh.controlVolume(aControlVolumeId).faceId(j);

                    m_fvmMesh.controlVolume(aControlVolumeId).calculateFaceArea(aFaceId);
                    m_fvmMesh.controlVolume(aControlVolumeId).calculateFaceCentroid(aFaceId);
                }

                m_u0[aControlVolumeId] = m_u[aControlVolumeId] = 0.0;
                m_v0[aControlVolumeId] = m_v[aControlVolumeId] = 0.0;
                m_w0[aControlVolumeId] = m_w[aControlVolumeId] = 0.0;

                m_p0[aControlVolumeId] = m_p[aControlVolumeId] = 0.0;

                m_mapIdToIndex[aControlVolumeId] = i;
                m_mapIndexToId[i] = aControlVolumeId;
            }

            for (Integer i = 0; i < m_fvmMesh.nbFaces(); ++i)
            {

                Integer aFaceId = m_fvmMesh.faceId(i);

                m_fvmMesh.face(aFaceId).calculateArea();
                m_fvmMesh.face(aFaceId).calculateNormal();

                m_flux[aFaceId] = 0.0;

                m_uf[aFaceId] = 0.0;
                m_vf[aFaceId] = 0.0;
                m_wf[aFaceId] = 0.0;

                m_pf[aFaceId] = 0.0;
            }

            m_calcu = m_calcv = m_calcw = false;
            m_calcp = false;
        }

        template <typename Real>
        CFvmPisoSolver<Real>::~CFvmPisoSolver()
        {
        }

        template <typename Real>
        CGeoVector<Real> CFvmPisoSolver<Real>::gradient(varMap& var, varMap& varf, const Integer aControlVolumeId)
        {

            CGeoVector<Real> aVector(0.0, 0.0, 0.0);

            for (Integer j = 0; j < m_fvmMesh.controlVolume(aControlVolumeId).nbFaces(); ++j)
            {

                Integer aFaceId = m_fvmMesh.controlVolume(aControlVolumeId).faceId(j);

                Real v = 0;
                CGeoNormal<Real>& aNormal = m_fvmMesh.controlVolume(aControlVolumeId).faceNormal(aFaceId);
                Real area = m_fvmMesh.controlVolume(aControlVolumeId).faceArea(aFaceId);

                if (m_fvmMesh.face(aFaceId).hasPair())
                {

                    Integer aNeighborId = m_fvmMesh.face(aFaceId).neighborId(aControlVolumeId);

                    v = var.at(aControlVolumeId) + var[aNeighborId];
                    v *= 0.5;
                    aVector += v * aNormal * area;
                }
                else
                {

                    v = varf.at(aFaceId);
                    aVector += v * aNormal * area;
                }
            }

            Real aVolume = m_fvmMesh.controlVolume(aControlVolumeId).originalVolume();

            if (aVolume > 0.0)
                aVector /= aVolume;

            return aVector;
        }

        template <typename Real>
        void CFvmPisoSolver<Real>::setTimeInterval(Real dt)
        {

            m_dt = dt;
        }

        template <typename Real>
        void CFvmPisoSolver<Real>::storePreviousQuantities()
        {

            m_u0 = m_u;
            m_v0 = m_v;
            m_w0 = m_w;

            m_p0 = m_p;
        }

        template <typename Real>
        void CFvmPisoSolver<Real>::setGravity(const Real gx, const Real gy, const Real gz)
        {

            m_gx = gx;
            m_gy = gy;
            m_gz = gz;
        }

        template <typename Real>
        void CFvmPisoSolver<Real>::setMaterialProperties(const Real aDensity, const Real aViscosity)
        {

            for (Integer i = 0; i < m_fvmMesh.nbControlVolumes(); ++i)
            {

                Integer aControlVolumeId = m_fvmMesh.controlVolumeId(i);

                m_dens[aControlVolumeId] = aDensity;
                m_visc[aControlVolumeId] = aViscosity;
            }
        }

        template <typename Real>
        void CFvmPisoSolver<Real>::setBoundaryVelocity(const std::vector<Integer>& sFaceIds, EBoundaryType sFaceType, const Real u, const Real v, const Real w)
        {

            for (Integer i = 0; i < static_cast<Integer>(sFaceIds.size()); ++i)
            {

                Integer aFaceId = sFaceIds[i];

                m_uf[aFaceId] = u;
                m_vf[aFaceId] = v;
                m_wf[aFaceId] = w;

                m_fvmMesh.face(aFaceId).setBoundaryType(CFvmBoundaryType(sFaceType));
            }
        }

        template <typename Real>
        void CFvmPisoSolver<Real>::setBoundaryPressure(const std::vector<Integer>& sFaceIds, EBoundaryType sFaceType, const Real p)
        {

            for (Integer i = 0; i < static_cast<Integer>(sFaceIds.size()); ++i)
            {

                Integer aFaceId = sFaceIds[i];

                m_pf[aFaceId] = p;

                m_fvmMesh.face(aFaceId).setBoundaryType(CFvmBoundaryType(sFaceType));
            }
        }

        template <typename Real>
        void CFvmPisoSolver<Real>::calculateVelocityField()
        {

            Eigen::SparseMatrix<Real> A;
            Eigen::Matrix<Real, Eigen::Dynamic, 1> bu, bv, bw;

            A.resize(m_fvmMesh.nbControlVolumes(), m_fvmMesh.nbControlVolumes());
            A.reserve(m_fvmMesh.nbControlVolumes());

            bu.resize(m_fvmMesh.nbControlVolumes());
            bu.setZero();

            bv.resize(m_fvmMesh.nbControlVolumes());
            bv.setZero();

            bw.resize(m_fvmMesh.nbControlVolumes());
            bw.setZero();

            // Assemble momentum matrix
            for (Integer i = 0; i < m_fvmMesh.nbControlVolumes(); ++i)
            {

                Integer aControlVolumeId = m_fvmMesh.controlVolumeId(i);

                Integer anIndexP = m_mapIdToIndex.at(aControlVolumeId);

                Real volume = m_fvmMesh.controlVolume(aControlVolumeId).originalVolume();

                if (volume <= 0.0)
                {
                    std::cout << "Warning: element " << aControlVolumeId << " has volume = 0!" << std::endl;
                    continue;
                }

                for (Integer j = 0; j < m_fvmMesh.controlVolume(aControlVolumeId).nbFaces(); ++j)
                {

                    Integer aFaceId = m_fvmMesh.controlVolume(aControlVolumeId).faceId(j);

                    Real dens = m_dens.at(aControlVolumeId);
                    Real visc = m_visc.at(aControlVolumeId);

                    CGeoNormal<Real>& aNormal = m_fvmMesh.controlVolume(aControlVolumeId).faceNormal(aFaceId);

                    Real area = m_fvmMesh.controlVolume(aControlVolumeId).faceArea(aFaceId);
                    Real dist = m_fvmMesh.controlVolume(aControlVolumeId).faceDist(aFaceId);

                    Real flux;

                    if (m_fvmMesh.face(aFaceId).controlVolumeId() == aControlVolumeId)
                        flux = +m_flux.at(aFaceId);
                    else
                        flux = -m_flux.at(aFaceId);

                    if (m_fvmMesh.face(aFaceId).hasPair())
                    {

                        Integer aNeighborId = m_fvmMesh.face(aFaceId).neighborId(aControlVolumeId);

                        Integer anIndexN = m_mapIdToIndex.at(aNeighborId);

                        if (m_fvmMesh.controlVolume(aNeighborId).containsFace(aFaceId))
                        {

                            area += m_fvmMesh.controlVolume(aNeighborId).faceArea(aFaceId);
                            area *= 0.5;
                            dist += m_fvmMesh.controlVolume(aNeighborId).faceDist(aFaceId);

                            dens += m_dens.at(aNeighborId);
                            dens *= 0.5;
                            visc += m_visc.at(aNeighborId);
                            visc *= 0.5;

                            Real xsi = 0.5; // CDS

                            // Convection
                            A.coeffRef(anIndexP, anIndexP) += (1.0 - xsi) * dens * flux / volume;

                            // Diffusion
                            A.coeffRef(anIndexP, anIndexP) += +visc * area / dist / volume;

                            // Convection
                            A.coeffRef(anIndexP, anIndexN) += xsi * dens * flux / volume;

                            // Diffusion
                            A.coeffRef(anIndexP, anIndexN) += -visc * area / dist / volume;
                        }
                    }
                    else
                    {

                        // Diffusion
                        A.coeffRef(anIndexP, anIndexP) += visc * area / dist / volume;

                        bu[anIndexP] += visc * m_uf[aFaceId] * area / dist / volume;
                        bv[anIndexP] += visc * m_vf[aFaceId] * area / dist / volume;
                        bw[anIndexP] += visc * m_wf[aFaceId] * area / dist / volume;

                        // Convection
                        bu[anIndexP] += -dens * m_uf[aFaceId] * area * flux * aNormal.x() / volume;
                        bv[anIndexP] += -dens * m_vf[aFaceId] * area * flux * aNormal.y() / volume;
                        bw[anIndexP] += -dens * m_wf[aFaceId] * area * flux * aNormal.z() / volume;
                    }
                }

                // Source - gravity
                bu[anIndexP] += m_dens.at(aControlVolumeId) * m_gx;
                bv[anIndexP] += m_dens.at(aControlVolumeId) * m_gy;
                bw[anIndexP] += m_dens.at(aControlVolumeId) * m_gz;

                // Source - pressure
                CGeoVector<Real> gradp = this->gradient(m_p, m_pf, aControlVolumeId);

                bu[anIndexP] += gradp.x();
                bv[anIndexP] += gradp.y();
                bw[anIndexP] += gradp.z();

                if (m_dt > 0.0)
                {

                    A.coeffRef(anIndexP, anIndexP) += m_dens.at(aControlVolumeId) / m_dt;

                    bu[anIndexP] += m_dens.at(aControlVolumeId) / m_dt * m_u0.at(aControlVolumeId);
                    bv[anIndexP] += m_dens.at(aControlVolumeId) / m_dt * m_v0.at(aControlVolumeId);
                    bw[anIndexP] += m_dens.at(aControlVolumeId) / m_dt * m_w0.at(aControlVolumeId);
                }
            }

            A.finalize();

            Eigen::BiCGSTAB<Eigen::SparseMatrix<Real>> solver;
            solver.compute(A);

            Eigen::Matrix<Real, Eigen::Dynamic, 1> u = solver.solve(bu);
            Eigen::Matrix<Real, Eigen::Dynamic, 1> v = solver.solve(bv);
            Eigen::Matrix<Real, Eigen::Dynamic, 1> w = solver.solve(bw);

            for (int k = 0; k < A.outerSize(); ++k)
            {

                Integer aControlVolumeId = m_mapIndexToId.at(k);

                CGeoVector<Real> gradp = this->gradient(m_p, m_pf, aControlVolumeId);

                m_Hu[aControlVolumeId] = bu[k] - gradp.x();
                m_Hv[aControlVolumeId] = bv[k] - gradp.y();
                m_Hw[aControlVolumeId] = bw[k] - gradp.z();

                for (typename Eigen::SparseMatrix<Real>::InnerIterator it(A, k); it; ++it)
                {

                    if (it.row() == it.col())
                        m_ap[aControlVolumeId] = it.value();
                    else
                    {
                        m_Hu.at(aControlVolumeId) += -it.value() * u[k];
                        m_Hv.at(aControlVolumeId) += -it.value() * v[k];
                        m_Hw.at(aControlVolumeId) += -it.value() * w[k];
                    }
                }
            }
        }

        template <typename Real>
        void CFvmPisoSolver<Real>::calculatePressureField()
        {

            Eigen::SparseMatrix<Real> A;
            Eigen::Matrix<Real, Eigen::Dynamic, 1> b;

            A.resize(m_fvmMesh.nbControlVolumes(), m_fvmMesh.nbControlVolumes());
            A.reserve(m_fvmMesh.nbControlVolumes());

            b.resize(m_fvmMesh.nbControlVolumes());
            b.setZero();

            // Assemble pressure matrix
            for (Integer i = 0; i < m_fvmMesh.nbControlVolumes(); ++i)
            {

                Integer aControlVolumeId = m_fvmMesh.controlVolumeId(i);

                Integer anIndexP = m_mapIdToIndex.at(aControlVolumeId);

                for (Integer j = 0; j < m_fvmMesh.controlVolume(aControlVolumeId).nbFaces(); ++j)
                {

                    Integer aFaceId = m_fvmMesh.controlVolume(aControlVolumeId).faceId(j);

                    CGeoNormal<Real>& aNormal = m_fvmMesh.controlVolume(aControlVolumeId).faceNormal(aFaceId);

                    Real area = m_fvmMesh.controlVolume(aControlVolumeId).faceArea(aFaceId);
                    Real dist = m_fvmMesh.controlVolume(aControlVolumeId).faceDist(aFaceId);

                    Real apj = m_ap.at(aControlVolumeId);

                    Real Huj = m_Hu.at(aControlVolumeId);
                    Real Hvj = m_Hv.at(aControlVolumeId);
                    Real Hwj = m_Hw.at(aControlVolumeId);

                    if (m_fvmMesh.face(aFaceId).hasPair())
                    {

                        Integer aNeighborId = m_fvmMesh.face(aFaceId).neighborId(aControlVolumeId);

                        Integer anIndexN = m_mapIdToIndex.at(aNeighborId);

                        if (m_fvmMesh.controlVolume(aNeighborId).containsFace(aFaceId))
                        {

                            area += m_fvmMesh.controlVolume(aNeighborId).faceArea(aFaceId);
                            area *= 0.5;
                            dist += m_fvmMesh.controlVolume(aNeighborId).faceDist(aFaceId);

                            apj += m_ap.at(aNeighborId);
                            apj *= 0.5;

                            Huj += m_Hu.at(aNeighborId);
                            Huj *= 0.5;
                            Hvj += m_Hv.at(aNeighborId);
                            Hvj *= 0.5;
                            Hwj += m_Hw.at(aNeighborId);
                            Hwj *= 0.5;

                            Real Hf = Huj * aNormal.x() + Hvj * aNormal.y() + Hwj * aNormal.z();

                            A.coeffRef(anIndexP, anIndexP) += -1.0 / (apj * dist) * area;
                            A.coeffRef(anIndexP, anIndexN) += +1.0 / (apj * dist) * area;

                            m_flux[aFaceId] = -Hf / apj * area;

                            b[anIndexP] -= m_flux.at(aFaceId);
                        }
                    }
                    else
                    {

                        Real Hf = Huj * aNormal.x() + Hvj * aNormal.y() + Hwj * aNormal.z();

                        if (m_fvmMesh.face(aFaceId).boundaryType() == BT_WALL_NO_SLIP)
                        {

                            // specified flux = 0
                            // pressure gradient = 0
                        }
                        else if (m_fvmMesh.face(aFaceId).boundaryType() == BT_INLET_FLOW)
                        {

                            // specified flux
                            // pressure gradient = 0
                            b[anIndexP] -= m_flux.at(aFaceId);

                            m_flux[aFaceId] = -(m_uf[aFaceId] * aNormal.x() + m_vf[aFaceId] * aNormal.y() + m_wf[aFaceId] * aNormal.z()) * area;
                        }
                        else if (m_fvmMesh.face(aFaceId).boundaryType() == BT_INLET_PRESSURE)
                        {

                            // specified pressure
                            // velocity gradient = 0

                            A.coeffRef(anIndexP, anIndexP) += -1.0 / (apj * dist) * area;
                            b[anIndexP] += -1.0 / (apj * dist) * area * m_pf.at(aFaceId);

                            m_flux[aFaceId] = -Hf / apj * area;

                            b[anIndexP] -= m_flux.at(aFaceId);

                            m_uf[aFaceId] = -m_flux[aFaceId] / area * aNormal.x();
                            m_vf[aFaceId] = -m_flux[aFaceId] / area * aNormal.y();
                            m_wf[aFaceId] = -m_flux[aFaceId] / area * aNormal.z();
                        }
                        else if (m_fvmMesh.face(aFaceId).boundaryType() == BT_OUTLET)
                        {

                            // specified pressure = 0
                            // velocity gradient = 0

                            A.coeffRef(anIndexP, anIndexP) += -1.0 / (apj * dist) * area;

                            m_flux[aFaceId] = -Hf / apj * area;

                            b[anIndexP] -= m_flux.at(aFaceId);

                            m_uf[aFaceId] = -m_flux[aFaceId] / area * aNormal.x();
                            m_vf[aFaceId] = -m_flux[aFaceId] / area * aNormal.y();
                            m_wf[aFaceId] = -m_flux[aFaceId] / area * aNormal.z();
                        }
                    }
                }
            }

            A.finalize();

            Eigen::ConjugateGradient<Eigen::SparseMatrix<Real>> solver;
            solver.compute(A);

            Eigen::Matrix<Real, Eigen::Dynamic, 1> p = solver.solve(b);

            for (int i = 0; i < p.rows(); ++i)
            {

                Integer aControlVolumeId = m_mapIndexToId.at(i);

                m_p[aControlVolumeId] = p[i];
            }
        }

        template <typename Real>
        void CFvmPisoSolver<Real>::correctFlux()
        {

            for (Integer i = 0; i < m_fvmMesh.nbControlVolumes(); ++i)
            {

                Integer aControlVolumeId = m_fvmMesh.controlVolumeId(i);

                for (Integer j = 0; j < m_fvmMesh.controlVolume(aControlVolumeId).nbFaces(); ++j)
                {

                    Integer aFaceId = m_fvmMesh.controlVolume(aControlVolumeId).faceId(j);

                    Real area = m_fvmMesh.controlVolume(aControlVolumeId).faceArea(aFaceId);
                    Real dist = m_fvmMesh.controlVolume(aControlVolumeId).faceDist(aFaceId);

                    Real apj = m_ap.at(aControlVolumeId);

                    if (m_fvmMesh.face(aFaceId).hasPair())
                    {

                        Integer aNeighborId = m_fvmMesh.face(aFaceId).neighborId(aControlVolumeId);

                        if (m_fvmMesh.controlVolume(aNeighborId).containsFace(aFaceId))
                        {

                            area += m_fvmMesh.controlVolume(aNeighborId).faceArea(aFaceId);
                            area *= 0.5;
                            dist += m_fvmMesh.controlVolume(aNeighborId).faceDist(aFaceId);

                            apj += m_ap.at(aNeighborId);
                            apj *= 0.5;

                            if (m_fvmMesh.face(aFaceId).controlVolumeId() == aControlVolumeId)
                                m_flux[aFaceId] -= area / (apj * dist) * (m_p[aNeighborId] - m_p.at(aControlVolumeId));
                        }
                    }
                    else
                    {

                        if (m_fvmMesh.face(aFaceId).boundaryType() == BT_WALL_NO_SLIP)
                        {

                            m_flux[aFaceId] = 0.0;
                            m_pf[aFaceId] = m_p.at(aControlVolumeId);
                        }
                        else if (m_fvmMesh.face(aFaceId).boundaryType() == BT_INLET_FLOW)
                        {

                            m_pf[aFaceId] = m_p.at(aControlVolumeId);
                        }
                        else if (m_fvmMesh.face(aFaceId).boundaryType() == BT_INLET_PRESSURE)
                        {

                            m_flux[aFaceId] += area / (apj * dist) * (m_pf[aFaceId] - m_p.at(aControlVolumeId));
                        }
                        else if (m_fvmMesh.face(aFaceId).boundaryType() == BT_OUTLET)
                        {

                            m_flux[aFaceId] += area / (apj * dist) * m_p.at(aControlVolumeId);
                            m_pf[aFaceId] = 0.0;
                        }
                    }
                }
            }
        }

        template <typename Real>
        void CFvmPisoSolver<Real>::correctVelocityField()
        {

            for (Integer i = 0; i < m_fvmMesh.nbControlVolumes(); ++i)
            {

                Integer aControlVolumeId = m_fvmMesh.controlVolumeId(i);

                CGeoVector<Real> gradp = this->gradient(m_p, m_pf, aControlVolumeId);

                m_u[aControlVolumeId] = (m_Hu.at(aControlVolumeId) - gradp.x()) / m_ap.at(aControlVolumeId);
                m_v[aControlVolumeId] = (m_Hv.at(aControlVolumeId) - gradp.y()) / m_ap.at(aControlVolumeId);
                m_w[aControlVolumeId] = (m_Hw.at(aControlVolumeId) - gradp.z()) / m_ap.at(aControlVolumeId);
            }
        }

        template <typename Real>
        void CFvmPisoSolver<Real>::correctPressureField()
        {

            Real p_min = std::numeric_limits<Real>::max();

            for (Integer i = 0; i < m_fvmMesh.nbControlVolumes(); ++i)
            {

                Integer aControlVolumeId = m_fvmMesh.controlVolumeId(i);

                p_min = std::min(p_min, m_p.at(aControlVolumeId));
            }

            for (Integer i = 0; i < m_fvmMesh.nbControlVolumes(); ++i)
            {

                Integer aControlVolumeId = m_fvmMesh.controlVolumeId(i);

                m_p.at(aControlVolumeId) -= p_min;
            }

            for (Integer i = 0; i < m_fvmMesh.nbFaces(); ++i)
            {

                Integer aFaceId = m_fvmMesh.faceId(i);

                m_pf[aFaceId] -= p_min;
            }
        }

        template <typename Real>
        void CFvmPisoSolver<Real>::iterate(const Real dt, const bool bInit)
        {

            this->setTimeInterval(dt);

            this->storePreviousQuantities();

            this->calculateVelocityField();
            this->calculatePressureField();

            this->correctFlux();
            this->correctVelocityField();
            this->correctPressureField();
        }

        template <typename Real>
        Real CFvmPisoSolver<Real>::checkMassConservationCell(const Integer aControlVolumeId)
        {

            Real sumCellFlux = 0.0;

            for (Integer j = 0; j < m_fvmMesh.controlVolume(aControlVolumeId).nbFaces(); ++j)
            {

                Integer aFaceId = m_fvmMesh.controlVolume(aControlVolumeId).faceId(j);

                Real flux;

                if (m_fvmMesh.face(aFaceId).controlVolumeId() == aControlVolumeId)
                    flux = +m_flux.at(aFaceId);
                else
                    flux = -m_flux.at(aFaceId);

                sumCellFlux += flux;
            }

            return sumCellFlux;
        }

        template <typename Real>
        void CFvmPisoSolver<Real>::checkMassConservation(Real& aMassError)
        {

            Real sumFlux = 0.0;

            for (Integer i = 0; i < m_fvmMesh.nbControlVolumes(); ++i)
            {

                Integer aControlVolumeId = m_fvmMesh.controlVolumeId(i);

                Real sumCellFlux = this->checkMassConservationCell(aControlVolumeId);

                m_massError[aControlVolumeId] = sumCellFlux;

                sumFlux += fabs(sumCellFlux);
            }

            if (m_fvmMesh.nbControlVolumes() > 0)
                sumFlux /= m_fvmMesh.nbControlVolumes();

            aMassError = sumFlux;
        }

        template <typename Real>
        void CFvmPisoSolver<Real>::residual(Real& ru, Real& rv, Real& rw, Real& rp)
        {

            ru = 0;
            rv = 0;
            rw = 0;
            rp = 0;

            for (Integer i = 0; i < m_fvmMesh.nbControlVolumes(); ++i)
            {

                Integer aControlVolumeId = m_fvmMesh.controlVolumeId(i);

                ru += (m_u.at(aControlVolumeId) - m_u0.at(aControlVolumeId)) * (m_u.at(aControlVolumeId) - m_u0.at(aControlVolumeId));
                rv += (m_v.at(aControlVolumeId) - m_v0.at(aControlVolumeId)) * (m_v.at(aControlVolumeId) - m_v0.at(aControlVolumeId));
                rw += (m_w.at(aControlVolumeId) - m_w0.at(aControlVolumeId)) * (m_w.at(aControlVolumeId) - m_w0.at(aControlVolumeId));
                rp += (m_p.at(aControlVolumeId) - m_p0.at(aControlVolumeId)) * (m_p.at(aControlVolumeId) - m_p0.at(aControlVolumeId));
            }
        }

        template <typename Real>
        Real CFvmPisoSolver<Real>::u(const Integer aControlVolumeId)
        {

            return m_u.at(aControlVolumeId);
        }

        template <typename Real>
        Real CFvmPisoSolver<Real>::v(const Integer aControlVolumeId)
        {

            return m_v.at(aControlVolumeId);
        }

        template <typename Real>
        Real CFvmPisoSolver<Real>::w(const Integer aControlVolumeId)
        {

            return m_w.at(aControlVolumeId);
        }

        template <typename Real>
        Real CFvmPisoSolver<Real>::p(const Integer aControlVolumeId)
        {

            return m_p.at(aControlVolumeId);
        }

        template <typename Real>
        Real CFvmPisoSolver<Real>::uf(const Integer aFaceId)
        {

            return m_uf.at(aFaceId);
        }

        template <typename Real>
        Real CFvmPisoSolver<Real>::vf(const Integer aFaceId)
        {

            return m_vf.at(aFaceId);
        }

        template <typename Real>
        Real CFvmPisoSolver<Real>::wf(const Integer aFaceId)
        {

            return m_wf.at(aFaceId);
        }

        template <typename Real>
        Real CFvmPisoSolver<Real>::pf(const Integer aFaceId)
        {

            return m_pf.at(aFaceId);
        }

        template <typename Real>
        Real CFvmPisoSolver<Real>::flux(const Integer aFaceId)
        {

            return m_flux.at(aFaceId);
        }

        template <typename Real>
        CGeoVector<Real> CFvmPisoSolver<Real>::gradu(const Integer aControlVolumeId)
        {

            return this->gradient(m_u, m_uf, aControlVolumeId);
        }

        template <typename Real>
        CGeoVector<Real> CFvmPisoSolver<Real>::gradv(const Integer aControlVolumeId)
        {

            return this->gradient(m_v, m_vf, aControlVolumeId);
        }

        template <typename Real>
        CGeoVector<Real> CFvmPisoSolver<Real>::gradw(const Integer aControlVolumeId)
        {

            return this->gradient(m_w, m_wf, aControlVolumeId);
        }

        template <typename Real>
        CGeoVector<Real> CFvmPisoSolver<Real>::gradp(const Integer aControlVolumeId)
        {

            return this->gradient(m_p, m_pf, aControlVolumeId);
        }

        template <typename Real>
        Real CFvmPisoSolver<Real>::massError(const Integer aControlVolumeId)
        {

            return m_massError.at(aControlVolumeId);
        }
    }
}
