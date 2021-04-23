// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include <map>
#include <vector>

#include "MshMesh.hpp"

#include "FvmControlVolume.hpp"
#include "FvmFace.hpp"
#include "FvmNode.hpp"

using namespace ENigMA::mesh;

namespace ENigMA
{
    namespace fvm
    {
        template <typename Real>
        class CFvmMesh
        {
        private:
            typedef std::map<Integer, CFvmFace<Real>> mapFace;
            typedef std::map<Integer, Integer> mapFaceIndex;

            typedef std::map<Integer, CFvmControlVolume<Real>> mapControlVolume;
            typedef std::map<Integer, Integer> mapControlVolumeIndex;

            std::vector<Integer> m_faceIds;
            std::vector<Integer> m_controlVolumeIds;

            Integer m_faceIndex;
            Integer m_controlVolumeIndex;

            mapFaceIndex m_faceIndices;
            mapControlVolumeIndex m_controlVolumeIndices;

            mapFace m_faces;
            mapControlVolume m_controlVolumes;

            CMshMesh<Real> m_mesh;

        public:
            CFvmMesh();
            explicit CFvmMesh(const CMshMesh<Real>& aMesh);
            virtual ~CFvmMesh();

            void set(const CMshMesh<Real>& aMesh);
            CMshMesh<Real>& mesh();

            void reset();

            Integer nbFaces();
            void addFace(const Integer aFaceId, const CFvmFace<Real>& aFace);
            void setFace(const Integer aFaceId, const CFvmFace<Real>& aFace);
            void removeFace(const Integer aFaceId);
            Integer faceId(const Integer aFaceIndex);
            Integer faceIndex(const Integer aFaceId);
            CFvmFace<Real>& face(const Integer aFaceId);

            Integer nbControlVolumes();
            void addControlVolume(const Integer aControlVolumeId, const CFvmControlVolume<Real>& aControlVolume);
            Integer controlVolumeId(const Integer aControlVolumeIndex);
            Integer controlVolumeIndex(const Integer aControlVolumeId);
            CFvmControlVolume<Real>& controlVolume(const Integer aControlVolumeId);

            Integer nextFaceId();
            Integer nextControlVolumeId();

            Real volume();
        };
    }
}

#include "FvmMesh_Imp.hpp"
