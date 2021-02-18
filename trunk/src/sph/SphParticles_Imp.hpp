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
namespace sph {
    template <typename Real>
    CSphParticles<Real>::CSphParticles(CSphKernel<Real>& kernel) : m_bCyclic(false)
    {
        m_kernel = &kernel;
    }

    template <typename Real>
    CSphParticles<Real>::~CSphParticles()
    {
    }

    template <typename Real>
    void CSphParticles<Real>::init(CPdeField<Real>& aField, Real mass, Real rho, Real diff, Real h, Real dt, bool bCyclic)
    {
        m_particles.resize(aField.mesh().nbNodes());

        for (Integer i = 0; i < aField.mesh().nbNodes(); ++i) {
            m_particles[i].mass = mass;
            m_particles[i].density = rho;
            m_particles[i].conductivity = diff;

            m_particles[i].velocity = CGeoVector<Real>(0, 0, 0);
        }

        m_bCyclic = bCyclic;

        m_h = h;
        m_dt = dt;
    }

    template <typename Real>
    void CSphParticles<Real>::setBoundary(const CGeoBoundingBox<Real>& aBoundary)
    {
        m_boundary = aBoundary;
    }

    template <typename Real>
    void CSphParticles<Real>::setInitialVelocity(CPdeField<Real>& aField, const CGeoVector<Real>& aVelocity)
    {
        for (Integer i = 0; i < aField.mesh().nbNodes(); ++i) {
            m_particles[i].velocity = aVelocity;
        }
    }

    template <typename Real>
    void CSphParticles<Real>::buildHashGrid(CPdeField<Real>& aField, CGeoHashGrid<Real>& aHashGrid)
    {
        // Build node hash grid
        for (Integer i = 0; i < aField.mesh().nbNodes(); ++i) {
            Integer aNodeId = aField.mesh().nodeId(i);
            CMshNode<Real>& aNode = aField.mesh().node(aNodeId);

            aHashGrid.addGeometricObject(aNodeId, aNode);
        }

        aHashGrid.build();
    }

    template <typename Real>
    void CSphParticles<Real>::advectParticles(CPdeField<Real>& aField)
    {
        for (Integer i = 0; i < aField.mesh().nbNodes(); ++i) {
            Integer aNodeId1 = aField.mesh().nodeId(i);
            CMshNode<Real>& aNode1 = aField.mesh().node(aNodeId1);

            aNode1.x() += m_particles[i].velocity.x() * m_dt;
            aNode1.y() += m_particles[i].velocity.z() * m_dt;
            aNode1.z() += m_particles[i].velocity.y() * m_dt;

            if (m_bCyclic) {
                if (aNode1.x() >= m_boundary.max().x())
                    aNode1.x() = m_boundary.min().x();
                else if (aNode1.x() <= m_boundary.min().x())
                    aNode1.x() = m_boundary.max().x();

                if (aNode1.y() >= m_boundary.max().y())
                    aNode1.y() = m_boundary.min().y();
                else if (aNode1.y() <= m_boundary.min().y())
                    aNode1.y() = m_boundary.max().y();

                if (aNode1.z() >= m_boundary.max().z())
                    aNode1.z() = m_boundary.min().z();
                else if (aNode1.z() <= m_boundary.min().z())
                    aNode1.z() = m_boundary.max().z();
            }
        }
    }

    template <typename Real>
    void CSphParticles<Real>::addDiffusion(CPdeField<Real>& aField, CGeoHashGrid<Real>& aHashGrid)
    {
        for (Integer i = 0; i < aField.mesh().nbNodes(); ++i) {

            if (aField.uFixed.find(i) == aField.uFixed.end()) {
                Real fd = 0.0;

                Integer aNodeId1 = aField.mesh().nodeId(i);
                CMshNode<Real>& aNode1 = aField.mesh().node(aNodeId1);

                std::vector<Integer> sNodeIds;
                aHashGrid.find(sNodeIds, aNode1, m_h);

                for (Integer n = 0; n < static_cast<Integer>(sNodeIds.size()); ++n) {
                    Integer aNodeId2 = sNodeIds[n];
                    CMshNode<Real>& aNode2 = aField.mesh().node(aNodeId2);

                    Integer j = aField.mesh().nodeIndex(aNodeId2);
                    Real Tj = aField.u(j);

                    CGeoVector<Real> r = aNode2 - aNode1;

                    fd += (Tj - aField.u(i)) * m_kernel->laplacianW(r, m_h);
                }

                aField.u(i) += m_particles[i].conductivity * fd * m_dt;

            } else
                aField.u(i) = aField.uFixed.at(i);
        }
    }

    template <typename Real>
    void CSphParticles<Real>::solve(CPdeField<Real>& aField)
    {
        this->advectParticles(aField);

        CGeoHashGrid<Real> aHashGrid;
        this->buildHashGrid(aField, aHashGrid);

        this->addDiffusion(aField, aHashGrid);
    }
}
}
