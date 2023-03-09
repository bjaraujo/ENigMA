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
    namespace mesh
    {
        template <typename Real>
        CMshFace<Real>::CMshFace()
            : m_faceType(FT_NONE)
            , m_hasPair(false)
            , m_pairFaceId(0)
            , m_elementId(0)
        {
        }

        template <typename Real>
        CMshFace<Real>::CMshFace(const CMshFace<Real>& aFace)
        {
            m_nodeIds = aFace.m_nodeIds;
            m_faceType = aFace.m_faceType;
            m_hasPair = aFace.m_hasPair;
            m_pairFaceId = aFace.m_pairFaceId;
            m_elementId = aFace.m_elementId;
        }

        template <typename Real>
        CMshFace<Real>::CMshFace(EFaceType aFaceType)
            : m_faceType(aFaceType)
            , m_hasPair(false)
            , m_pairFaceId(0)
            , m_elementId(0)
        {
        }

        template <typename Real>
        CMshFace<Real>::~CMshFace()
        {
        }

        template <typename Real>
        Integer CMshFace<Real>::nbNodeIds() const
        {
            return static_cast<Integer>(m_nodeIds.size());
        }

        template <typename Real>
        void CMshFace<Real>::addNodeId(const Integer aNodeId)
        {
            m_nodeIds.emplace_back(aNodeId);
        }

        template <typename Real>
        Integer CMshFace<Real>::nodeId(const Integer aNodeIndex) const
        {
            return m_nodeIds.at(aNodeIndex);
        }

        template <typename Real>
        void CMshFace<Real>::setNodeId(const Integer aNodeIndex, const Integer aNodeId)
        {
            m_nodeIds[aNodeIndex] = aNodeId;
        }

        template <typename Real>
        void CMshFace<Real>::setPairFaceId(const Integer aPairFaceId)
        {
            m_pairFaceId = aPairFaceId;
            m_hasPair = true;
        }

        template <typename Real>
        Integer CMshFace<Real>::pairFaceId() const
        {
            return m_pairFaceId;
        }

        template <typename Real>
        bool CMshFace<Real>::hasPair() const
        {
            return m_hasPair;
        }

        template <typename Real>
        void CMshFace<Real>::setHasPair(bool hasPair)
        {
            m_hasPair = hasPair;
        }

        template <typename Real>
        void CMshFace<Real>::setElementId(const Integer anElementId)
        {
            m_elementId = anElementId;
        }

        template <typename Real>
        Integer CMshFace<Real>::elementId() const
        {
            return m_elementId;
        }

        template <typename Real>
        void CMshFace<Real>::setFaceType(EFaceType aFaceType)
        {
            m_faceType = aFaceType;
        }

        template <typename Real>
        EFaceType CMshFace<Real>::faceType()
        {
            return m_faceType;
        }

        template <typename Real>
        void CMshFace<Real>::reset()
        {
            m_nodeIds.clear();

            m_faceType = FT_NONE;
            m_hasPair = false;
            m_pairFaceId = 0;
            m_elementId = 0;
        }

        template <typename Real>
        std::ostream& operator<<(std::ostream& output, CMshFace<Real>& aFace)
        {
            for (Integer i = 0; i < aFace.nbNodeIds(); ++i)
                output << aFace.nodeId(i) << " ";

            return output;
        }
    }
}
