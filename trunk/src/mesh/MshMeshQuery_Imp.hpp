// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include <algorithm>
#include <iterator>

namespace ENigMA {
namespace mesh {
    template <typename Real>
    CMshMeshQuery<Real>::CMshMeshQuery(CMshMesh<Real>& aMesh)
    {
        for (Integer i = 0; i < aMesh.nbElements(); ++i) {
            Integer anElementId = aMesh.elementId(i);
            CMshElement<Real>& anElement = aMesh.element(anElementId);

            for (Integer j = 0; j < anElement.nbNodeIds(); ++j) {
                Integer aNodeId = anElement.nodeId(j);

                if (m_nodeToElement.find(aNodeId) != m_nodeToElement.end()) {
                    if (std::find(m_nodeToElement[aNodeId].begin(), m_nodeToElement[aNodeId].end(), anElementId) == m_nodeToElement[aNodeId].end())
                        m_nodeToElement[aNodeId].push_back(anElementId);
                } else
                    m_nodeToElement[aNodeId].push_back(anElementId);
            }
        }
    }

    template <typename Real>
    CMshMeshQuery<Real>::~CMshMeshQuery()
    {
    }

    template <typename Real>
    void CMshMeshQuery<Real>::elementsSharingNode(const Integer aNodeId, std::vector<Integer>& sElementIds)
    {
        sElementIds.resize(m_nodeToElement[aNodeId].size());
        std::copy(m_nodeToElement[aNodeId].begin(), m_nodeToElement[aNodeId].end(), sElementIds.begin());
    }

    template <typename Real>
    void CMshMeshQuery<Real>::elementsSharingNodes(const Integer aNodeId1, const Integer aNodeId2, std::vector<Integer>& sElementIds)
    {
        sElementIds.clear();
        std::set_intersection(m_nodeToElement[aNodeId1].begin(), m_nodeToElement[aNodeId1].end(), m_nodeToElement[aNodeId2].begin(), m_nodeToElement[aNodeId2].end(), std::back_inserter(sElementIds));
    }
}
}
