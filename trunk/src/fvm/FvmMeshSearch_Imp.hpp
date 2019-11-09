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

namespace fvm {

    template <typename Real>
    CFvmMeshSearch<Real>::CFvmMeshSearch()
    {
    }

    template <typename Real>
    CFvmMeshSearch<Real>::CFvmMeshSearch(CFvmMesh<Real>& aMesh)
    {

        this->set(aMesh);
    }

    template <typename Real>
    CFvmMeshSearch<Real>::~CFvmMeshSearch()
    {
    }

    template <typename Real>
    void CFvmMeshSearch<Real>::set(CFvmMesh<Real>& aMesh)
    {

        m_mesh = &aMesh;
        m_boundaryFaceHashGrid.reset();
    }

    template <typename Real>
    void CFvmMeshSearch<Real>::build()
    {

        for (Integer i = 0; i < m_mesh->nbFaces(); ++i) {

            m_mesh->face(m_mesh->faceId(i)).calculateCentroid();

            if (!m_mesh->face(m_mesh->faceId(i)).hasPair())
                m_boundaryFaceHashGrid.addGeometricObject(m_mesh->faceId(i), m_mesh->face(m_mesh->faceId(i)).centroid());
        }

        m_boundaryFaceHashGrid.build();
    }

    template <typename Real>
    void CFvmMeshSearch<Real>::findClosestBoundaryFace(CGeoCoordinate<Real>& aCoordinate, Integer& aFaceId, const Real aTolerance)
    {

        Integer wFaceId = 0;

        std::vector<Integer> sCoordinates;

        m_boundaryFaceHashGrid.find(sCoordinates, aCoordinate, aTolerance);

        if (sCoordinates.size() > 0) {

            Real minDist = std::numeric_limits<Real>::max();

            for (Integer i = 0; i < static_cast<Integer>(sCoordinates.size()); ++i) {

                aFaceId = sCoordinates[i];

                m_mesh->face(aFaceId).calculateCentroid();

                CGeoVector<Real> d = m_mesh->face(aFaceId).centroid() - aCoordinate;

                Real thisDist = d.norm();

                if (thisDist < minDist) {
                    wFaceId = aFaceId;
                    minDist = thisDist;
                }
            }
        }

        aFaceId = wFaceId;
    }
}
}
