// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

namespace ENigMA
{

    namespace fvm
    {

        template <typename Real>
        CFvmMesh<Real>::CFvmMesh() : m_faceIndex(0), m_controlVolumeIndex(0)
        {

        }

        template <typename Real>
        CFvmMesh<Real>::CFvmMesh(CMshMesh<Real>& aMesh)
        {

            this->set(aMesh);

        }

        template <typename Real>
        CFvmMesh<Real>::~CFvmMesh()
        {

        }

        template <typename Real>
        void CFvmMesh<Real>::set(CMshMesh<Real>& aMesh)
        {

            this->reset();

            m_mesh = aMesh;

            std::map<Integer, Integer> newFaceIds;

            std::map<Integer, bool> faceAdded;

            for (Integer i = 0; i < aMesh.nbFaces(); ++i)
                faceAdded[aMesh.faceId(i)] = false;

            Integer aNewFaceId = 0;

            // MshFaces->FvmFaces
            for (Integer i = 0; i < aMesh.nbFaces(); ++i)
            {

                if (faceAdded[aMesh.faceId(i)])
                    continue;

                CMshFace<Real> aMshFace = aMesh.face(aMesh.faceId(i));

                CFvmFace<Real> aFvmFace;

                for (Integer j = 0; j < aMshFace.nbNodeIds(); ++j)
                {

                    Integer aNodeId = aMshFace.nodeId(j);

                    aFvmFace.addNode(CFvmNode<Real>(aMesh.node(aNodeId)));

                }

                aFvmFace.setHasPair(aMshFace.hasPair());

                aFvmFace.setControlVolumeId(aMshFace.elementId());

                if (aMshFace.hasPair())
                    aFvmFace.setNeighborId(aMesh.face(aMshFace.pairFaceId()).elementId());

                aFvmFace.close();

                this->addFace(aNewFaceId, aFvmFace);

                faceAdded[aMesh.faceId(i)] = true;
                newFaceIds[aMesh.faceId(i)] = aNewFaceId;

                if (aMshFace.hasPair())
                {
                    faceAdded[aMshFace.pairFaceId()] = true;
                    newFaceIds[aMshFace.pairFaceId()] = aNewFaceId;
                }

                aNewFaceId++;

            }

            // MshElements->FvmControlVolumes
            for (Integer i = 0; i < aMesh.nbElements(); ++i)
            {

                Integer anElementId = aMesh.elementId(i);

                CMshElement<Real> aMshElement = aMesh.element(anElementId);

                Integer aControlVolumeId = anElementId;

                CFvmControlVolume<Real> aControlVolume;

                aControlVolume.setControlVolumeId(aMesh.elementId(i));

                for (Integer j = 0; j < aMshElement.nbFaceIds(); ++j)
                {

                    Integer aFaceId = newFaceIds[aMshElement.faceId(j)];

                    CFvmFace<Real> aFace = m_faces[aFaceId];

                    if (aFace.controlVolumeId() != aControlVolumeId)
                    {
                        // Invert
                        aFace.polyline().invert();
                    }

                    aControlVolume.addFace(aFaceId, aFace);

                }

                this->addControlVolume(aControlVolumeId, aControlVolume);

            }

        }

        template <typename Real>
        CMshMesh<Real>& CFvmMesh<Real>::mesh()
        {

            return m_mesh;

        }

        template <typename Real>
        void CFvmMesh<Real>::reset()
        {

            m_faceIds.clear();
            m_faces.clear();

            m_controlVolumeIds.clear();
            m_controlVolumes.clear();

        }

        template <typename Real>
        Integer CFvmMesh<Real>::nbFaces()
        {

            return static_cast<Integer>(m_faceIds.size());

        }

        template <typename Real>
        Integer CFvmMesh<Real>::nbControlVolumes()
        {

            return static_cast<Integer>(m_controlVolumeIds.size());

        }

        template <typename Real>
        Integer CFvmMesh<Real>::faceId(const Integer aFaceIndex)
        {

            return m_faceIds[aFaceIndex];

        }

        template <typename Real>
        Integer CFvmMesh<Real>::faceIndex(const Integer aFaceId)
        {

            return m_faceIndices[aFaceId];

        }

        template <typename Real>
        Integer CFvmMesh<Real>::controlVolumeId(const Integer aControlVolumeIndex)
        {

            return m_controlVolumeIds[aControlVolumeIndex];

        }

        template <typename Real>
        Integer CFvmMesh<Real>::controlVolumeIndex(const Integer aControlVolumeId)
        {

            return m_controlVolumeIndices[aControlVolumeId];

        }

        template <typename Real>
        void CFvmMesh<Real>::addFace(const Integer aFaceId, const CFvmFace<Real>& aFace)
        {

            m_faces[aFaceId] = aFace;
            m_faceIds.push_back(aFaceId);

            m_faceIndices[aFaceId] = m_faceIndex;
            m_faceIndex++;

        }

        template <typename Real>
        void CFvmMesh<Real>::setFace(const Integer aFaceId, const CFvmFace<Real>& aFace)
        {

            m_faces[aFaceId] = aFace;

        }

        template <typename Real>
        void CFvmMesh<Real>::removeFace(const Integer aFaceId)
        {

            m_faces.erase(aFaceId);

            std::vector<Integer>::iterator it = std::find(m_faceIds.begin(), m_faceIds.end(), aFaceId);
            
            if (it != m_faceIds.end())
                m_faceIds.erase(it);

        }

        template <typename Real>
        void CFvmMesh<Real>::addControlVolume(const Integer aControlVolumeId, const CFvmControlVolume<Real>& aControlVolume)
        {

            m_controlVolumes[aControlVolumeId] = aControlVolume;
            m_controlVolumeIds.push_back(aControlVolumeId);

            m_controlVolumeIndices[aControlVolumeId] = m_controlVolumeIndex;
            m_controlVolumeIndex++;

        }

        template <typename Real>
        CFvmFace<Real>& CFvmMesh<Real>::face(const Integer aFaceId)
        {

            return m_faces[aFaceId];

        }

        template <typename Real>
        CFvmControlVolume<Real>& CFvmMesh<Real>::controlVolume(const Integer aControlVolumeId)
        {

            return m_controlVolumes[aControlVolumeId];

        }

        template <typename Real>
        Integer CFvmMesh<Real>::nextFaceId()
        {

            if (m_faces.size() > 0)
                return m_faces.rbegin()->first + 1;
            else
                return 0;

        }

        template <typename Real>
        Integer CFvmMesh<Real>::nextControlVolumeId()
        {

            if (m_controlVolumes.size() > 0)
                return m_controlVolumes.rbegin()->first + 1;
            else
                return 0;

        }

        template <typename Real>
        Real CFvmMesh<Real>::volume()
        {

            Real volume = 0.0;

            for (Integer i = 0; i < static_cast<Integer>(m_controlVolumeIds.size()); ++i)
            {

                Integer aControlVolumeId = m_controlVolumeIds[i];

                m_controlVolumes[aControlVolumeId].calculateVolume();

                volume += m_controlVolumes[aControlVolumeId].volume();

            }

            return volume;

        }

    }

}


