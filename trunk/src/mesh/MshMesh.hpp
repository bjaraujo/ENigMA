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

#include "CmnTypes.hpp"
#include "MshElement.hpp"
#include "MshFace.hpp"
#include "MshNode.hpp"

namespace ENigMA {

namespace mesh {

    template <typename Real>
    class CMshMesh {
    private:
        typedef std::map<Integer, CMshNode<Real>> mapNode;
        typedef std::map<Integer, CMshFace<Real>> mapFace;
        typedef std::map<Integer, CMshElement<Real>> mapElement;

        typedef std::map<Integer, Integer> mapNodeIndex;
        typedef std::map<Integer, Integer> mapFaceIndex;
        typedef std::map<Integer, Integer> mapElementIndex;

        std::vector<Integer> m_nodeIds;
        std::vector<Integer> m_faceIds;
        std::vector<Integer> m_elementIds;

        Integer m_nodeIndex;
        Integer m_faceIndex;
        Integer m_elementIndex;

        mapNodeIndex m_nodeIndices;
        mapFaceIndex m_faceIndices;
        mapElementIndex m_elementIndices;

        mapNode m_nodes;
        mapFace m_faces;
        mapElement m_elements;

        Integer m_nbBoundaryFaces;

        Real m_dx, m_dy, m_dz;
        
        std::vector<ENigMA::geometry::CGeoCoordinate<Real>> m_faceCentroid;
        std::vector<ENigMA::geometry::CGeoCoordinate<Real>> m_elementCentroid;
       
    public:
        CMshMesh();
        CMshMesh(const CMshMesh<Real>& aMesh);
        virtual ~CMshMesh();

        void reset();

        Integer nbNodes() const;
        Integer nbFaces() const;
        Integer nbElements() const;

        Integer nodeId(Integer aNodeIndex) const;
        Integer faceId(Integer aFaceIndex) const;
        Integer elementId(Integer anElementIndex) const;

        Integer nodeIndex(Integer aNodeId) const;
        Integer faceIndex(Integer aFaceId) const;
        Integer elementIndex(Integer anElementId) const;

        void addNode(const Integer aNodeId, const ENigMA::mesh::CMshNode<Real>& aNode);
        void addFace(const Integer aFaceId, const ENigMA::mesh::CMshFace<Real>& aFace);
        void addElement(const Integer anElementId, const ENigMA::mesh::CMshElement<Real>& anElement);

        void removeNode(const Integer aNodeId);
        void removeFace(const Integer aFaceId);
        void removeElement(const Integer anElementId);

        void addMesh(CMshMesh<Real>& aMesh);

        ENigMA::mesh::CMshNode<Real>& node(const Integer aNodeId);
        ENigMA::mesh::CMshFace<Real>& face(const Integer aFaceId);
        ENigMA::mesh::CMshElement<Real>& element(const Integer anElementId);

        void generateFaces(const Real aTolerance = 0.0);

        CMshMesh<Real> extractBoundary(const Real aTolerance = 0.0);
        Integer nbBoundaryFaces() const;

        void setDx(const Real aValue);
        Real dx() const;

        void setDy(const Real aValue);
        Real dy() const;

        void setDz(const Real aValue);
        Real dz() const;

        void calculateFaceCentroid();
        void calculateElementCentroid();

        ENigMA::geometry::CGeoCoordinate<Real>& faceCentroid(const Integer aFaceId);
        ENigMA::geometry::CGeoCoordinate<Real>& elementCentroid(const Integer anElementId);

        void scale(const Real aFactor);

        Integer nextNodeId();
        Integer nextFaceId();
        Integer nextElementId();

        void mergeNodes(const Real aTolerance = 0.0);
        void removeInvalidElements();
        void rebuildIndices();
        void removeDanglingNodes();
        void collapseNakedEdges(const Real aTolerance = 0.0);

        void renumber();

        void invert();

        void meshQuality(Real& aMinQ, Real& aMaxQ, Real& aAveQ);

        ENigMA::geometry::CGeoBoundingBox<Real> boundingBox(const Integer anElementId);
        ENigMA::geometry::CGeoBoundingBox<Real> boundingBox();
    };
}
}

#include "MshMesh_Imp.hpp"
