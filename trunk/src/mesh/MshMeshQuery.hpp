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

namespace ENigMA {

namespace mesh {

    template <typename Real>
    class CMshMeshQuery {
    private:
        typedef std::map<Integer, std::vector<Integer>> mapNodeToElement;

        mapNodeToElement m_nodeToElement;

    public:
        CMshMeshQuery(CMshMesh<Real>& aMesh);
        virtual ~CMshMeshQuery();

        void elementsSharingNode(const Integer aNodeId, std::vector<Integer>& sElementIds);
        void elementsSharingNodes(const Integer aNodeId1, const Integer aNodeId2, std::vector<Integer>& sElementIds);
    };
}
}

#include "MshMeshQuery_Imp.hpp"
