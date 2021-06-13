// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "GeoHashGrid.hpp"

#include "MshTriangle.hpp"
#include "MshQuadrilateral.hpp"
#include "MshTetrahedron.hpp"

namespace ENigMA
{
    namespace mesh
    {
        template <typename Real>
        CMshMesh<Real>::CMshMesh()
            : m_nodeIndex(0)
            , m_faceIndex(0)
            , m_elementIndex(0)
            , m_nbBoundaryFaces(0)
            , m_dx(0.0)
            , m_dy(0.0)
            , m_dz(0.0)
        {
        }

        template <typename Real>
        CMshMesh<Real>::CMshMesh(const CMshMesh<Real>& aMesh)
        {
            reset();

            for (Integer i = 0; i < static_cast<Integer>(aMesh.m_nodeIds.size()); ++i)
            {
                Integer aNodeId = aMesh.m_nodeIds.at(i);
                CMshNode<Real> aNode = aMesh.m_nodes.at(aNodeId);

                addNode(aNodeId, aNode);
            }

            for (Integer i = 0; i < static_cast<Integer>(aMesh.m_faceIds.size()); ++i)
            {
                Integer aFaceId = aMesh.m_faceIds.at(i);
                CMshFace<Real> aFace = aMesh.m_faces.at(aFaceId);

                addFace(aFaceId, aFace);
            }

            for (Integer i = 0; i < static_cast<Integer>(aMesh.m_elementIds.size()); ++i)
            {
                Integer anElementId = aMesh.m_elementIds.at(i);
                CMshElement<Real> anElement = aMesh.m_elements.at(anElementId);

                addElement(anElementId, anElement);
            }
        }

        template <typename Real>
        CMshMesh<Real>::~CMshMesh() { }

        template <typename Real>
        void CMshMesh<Real>::reset()
        {
            m_nodeIds.clear();
            m_nodes.clear();

            m_faceIds.clear();
            m_faces.clear();

            m_elementIds.clear();
            m_elements.clear();

            m_faceIndices.clear();
            m_nodeIndices.clear();
            m_elementIndices.clear();

            m_nodeIndex = 0;
            m_faceIndex = 0;
            m_elementIndex = 0;

            m_dx = 0.0;
            m_dy = 0.0;
            m_dz = 0.0;

            m_nbBoundaryFaces = 0;
        }

        template <typename Real>
        Integer CMshMesh<Real>::nbNodes() const { return static_cast<Integer>(m_nodeIds.size()); }

        template <typename Real>
        Integer CMshMesh<Real>::nbFaces() const { return static_cast<Integer>(m_faceIds.size()); }

        template <typename Real>
        Integer CMshMesh<Real>::nbElements() const { return static_cast<Integer>(m_elementIds.size()); }

        template <typename Real>
        void CMshMesh<Real>::addNode(const Integer aNodeId, const ENigMA::mesh::CMshNode<Real>& aNode)
        {
            m_nodes[aNodeId] = aNode;
            m_nodeIds.push_back(aNodeId);

            m_nodeIndices[aNodeId] = m_nodeIndex;
            m_nodeIndex++;
        }

        template <typename Real>
        void CMshMesh<Real>::addFace(const Integer aFaceId, const ENigMA::mesh::CMshFace<Real>& aFace)
        {
            m_faces[aFaceId] = aFace;
            m_faceIds.push_back(aFaceId);

            m_faceIndices[aFaceId] = m_faceIndex;
            m_faceIndex++;
        }

        template <typename Real>
        void CMshMesh<Real>::addElement(const Integer anElementId, const ENigMA::mesh::CMshElement<Real>& anElement)
        {
            m_elements[anElementId] = anElement;
            m_elementIds.push_back(anElementId);

            m_elementIndices[anElementId] = m_elementIndex;
            m_elementIndex++;
        }

        template <typename Real>
        void CMshMesh<Real>::removeNode(const Integer aNodeId)
        {
            m_nodes.erase(aNodeId);
            m_nodeIndices.erase(aNodeId);
            m_nodeIds.erase(std::find(m_nodeIds.begin(), m_nodeIds.end(), aNodeId));
        }

        template <typename Real>
        void CMshMesh<Real>::removeFace(const Integer aFaceId)
        {
            m_faces.erase(aFaceId);
            m_faceIndices.erase(aFaceId);
            m_faceIds.erase(std::find(m_faceIds.begin(), m_faceIds.end(), aFaceId));
        }

        template <typename Real>
        void CMshMesh<Real>::removeElement(const Integer anElementId)
        {
            m_elements.erase(anElementId);
            m_elementIndices.erase(anElementId);
            m_elementIds.erase(std::find(m_elementIds.begin(), m_elementIds.end(), anElementId));
        }

        template <typename Real>
        void CMshMesh<Real>::addMesh(const CMshMesh<Real>& aMesh)
        {
            mapNodeIndex sNewNodeIds;

            for (Integer i = 0; i < aMesh.nbNodes(); ++i)
            {
                Integer aNodeId = aMesh.nodeId(i);
                const CMshNode<Real>& aNode = aMesh.node(aNodeId);

                Integer aNewNodeId = nextNodeId();
                addNode(aNewNodeId, aNode);
                sNewNodeIds[aNodeId] = aNewNodeId;
            }

            for (Integer i = 0; i < aMesh.nbFaces(); ++i)
            {
                Integer aFaceId = aMesh.faceId(i);
                const CMshFace<Real>& aFace = aMesh.face(aFaceId);

                Integer aNewFaceId = nextFaceId();
                addFace(aNewFaceId, aFace);
                CMshFace<Real>& aNewFace = face(aNewFaceId);

                for (Integer j = 0; j < aNewFace.nbNodeIds(); ++j)
                {
                    Integer aNodeId = aNewFace.nodeId(j);
                    aNewFace.setNodeId(j, sNewNodeIds.at(aNodeId));
                }
            }

            for (Integer i = 0; i < aMesh.nbElements(); ++i)
            {
                Integer anElementId = aMesh.elementId(i);
                const CMshElement<Real>& anElement = aMesh.element(anElementId);

                Integer aNewElementId = nextElementId();
                addElement(aNewElementId, anElement);
                CMshElement<Real>& aNewElement = element(aNewElementId);

                for (Integer j = 0; j < aNewElement.nbNodeIds(); ++j)
                {
                    Integer aNodeId = aNewElement.nodeId(j);
                    aNewElement.setNodeId(j, sNewNodeIds.at(aNodeId));
                }
            }
        }

        template <typename Real>
        ENigMA::mesh::CMshNode<Real>& CMshMesh<Real>::node(const Integer aNodeId) { return m_nodes.at(aNodeId); }

        template <typename Real>
        const ENigMA::mesh::CMshNode<Real>& CMshMesh<Real>::node(const Integer aNodeId) const { return m_nodes.at(aNodeId); }

        template <typename Real>
        ENigMA::mesh::CMshFace<Real>& CMshMesh<Real>::face(const Integer aFaceId) { return m_faces.at(aFaceId); }

        template <typename Real>
        const ENigMA::mesh::CMshFace<Real>& CMshMesh<Real>::face(const Integer aFaceId) const { return m_faces.at(aFaceId); }

        template <typename Real>
        ENigMA::mesh::CMshElement<Real>& CMshMesh<Real>::element(const Integer anElementId) { return m_elements.at(anElementId); }

        template <typename Real>
        const ENigMA::mesh::CMshElement<Real>& CMshMesh<Real>::element(const Integer anElementId) const { return m_elements.at(anElementId); }

        template <typename Real>
        Integer CMshMesh<Real>::nodeId(const Integer aNodeIndex) const { return m_nodeIds.at(aNodeIndex); }

        template <typename Real>
        Integer CMshMesh<Real>::faceId(const Integer aFaceIndex) const { return m_faceIds.at(aFaceIndex); }

        template <typename Real>
        Integer CMshMesh<Real>::elementId(const Integer anElementIndex) const { return m_elementIds.at(anElementIndex); }

        template <typename Real>
        Integer CMshMesh<Real>::nodeIndex(const Integer aNodeId) const { return m_nodeIndices.at(aNodeId); }

        template <typename Real>
        Integer CMshMesh<Real>::faceIndex(const Integer aFaceId) const { return m_faceIndices.at(aFaceId); }

        template <typename Real>
        Integer CMshMesh<Real>::elementIndex(const Integer anElementId) const { return m_elementIndices.at(anElementId); }

        template <typename Real>
        void CMshMesh<Real>::generateFaces(const Real aTolerance)
        {
            m_faceIds.clear();
            m_faces.clear();

            m_faceIndices.clear();

            m_faceIndex = 0;

            for (Integer i = 0; i < static_cast<Integer>(m_elementIds.size()); ++i)
            {
                Integer anElementId = m_elementIds.at(i);

                std::vector<CMshFace<Real>> sFaces;

                m_elements.at(anElementId).generateFaces(sFaces);

                for (Integer j = 0; j < static_cast<Integer>(sFaces.size()); ++j)
                {
                    sFaces.at(j).setElementId(anElementId);
                    Integer aFaceId = nextFaceId();
                    addFace(aFaceId, sFaces.at(j));
                    m_elements.at(anElementId).addFaceId(aFaceId);
                }
            }

            // Discover double faces
            CGeoHashGrid<Real> aHashGrid;

            std::vector<CGeoCoordinate<Real>> sCenterCoordinates;

            for (Integer i = 0; i < static_cast<Integer>(m_faceIds.size()); ++i)
            {
                Integer aFaceId = m_faceIds.at(i);

                CGeoCoordinate<Real> aCenterCoordinate(0.0, 0.0, 0.0);

                for (Integer j = 0; j < m_faces.at(aFaceId).nbNodeIds(); ++j)
                {
                    Integer aNodeId = m_faces.at(aFaceId).nodeId(j);
                    CMshNode<Real> aNode = m_nodes.at(aNodeId);

                    aCenterCoordinate += aNode;
                }

                if (m_faces.at(m_faceIds.at(i)).nbNodeIds() > 0)
                    aCenterCoordinate /= static_cast<Real>(m_faces.at(aFaceId).nbNodeIds());

                sCenterCoordinates.push_back(aCenterCoordinate);

                aHashGrid.addGeometricObject(aFaceId, aCenterCoordinate);
            }

            aHashGrid.build();

            for (Integer i = 0; i < static_cast<Integer>(m_faceIds.size()); ++i)
            {
                Integer aFaceId = m_faceIds.at(i);

                std::vector<Integer> sCoordinates;

                aHashGrid.find(sCoordinates, sCenterCoordinates.at(i), aTolerance);

                for (Integer j = 0; j < static_cast<Integer>(sCoordinates.size()); ++j)
                {
                    if (sCoordinates.at(j) != aFaceId)
                    {
                        Integer aPairFaceId = sCoordinates.at(j);

                        std::vector<Integer> n1;
                        for (Integer k = 0; k < m_faces.at(aFaceId).nbNodeIds(); ++k)
                        {
                            Integer aNodeId = m_faces.at(aFaceId).nodeId(k);
                            n1.push_back(aNodeId);
                        }
                        std::sort(n1.begin(), n1.end());

                        std::vector<Integer> n2;
                        for (Integer k = 0; k < m_faces.at(aPairFaceId).nbNodeIds(); ++k)
                        {
                            Integer aNodeId = m_faces.at(aPairFaceId).nodeId(k);
                            n2.push_back(aNodeId);
                        }
                        std::sort(n2.begin(), n2.end());

                        if (m_faces.at(aFaceId).faceType() == m_faces.at(aPairFaceId).faceType() && n1 == n2)
                        {
                            m_faces.at(aFaceId).setPairFaceId(aPairFaceId);
                            m_faces.at(aPairFaceId).setPairFaceId(aFaceId);
                            break;
                        }
                    }
                }
            }

            m_nbBoundaryFaces = 0;

            for (typename mapFace::iterator it = m_faces.begin(); it != m_faces.end(); ++it)
            {
                if (!it->second.hasPair())
                    m_nbBoundaryFaces++;
            }
        }

        template <typename Real>
        Integer CMshMesh<Real>::nbBoundaryFaces() const { return m_nbBoundaryFaces; }

        template <typename Real>
        ENigMA::mesh::CMshMesh<Real> CMshMesh<Real>::extractBoundary(const Real aTolerance)
        {
            CMshMesh<Real> aBoundaryMesh;

            generateFaces(aTolerance);

            Integer anElementId = 0;

            for (Integer i = 0; i < static_cast<Integer>(m_faceIds.size()); ++i)
            {
                Integer aFaceId = m_faceIds.at(i);

                CMshFace<Real> aFace = m_faces.at(aFaceId);

                if (!aFace.hasPair())
                {
                    CMshElement<Real> anElement;

                    if (aFace.faceType() == FT_LINE)
                        anElement.setElementType(ET_BEAM);

                    if (aFace.faceType() == FT_TRIANGLE)
                        anElement.setElementType(ET_TRIANGLE);

                    if (aFace.faceType() == FT_QUADRILATERAL)
                        anElement.setElementType(ET_QUADRILATERAL);

                    for (Integer j = 0; j < aFace.nbNodeIds(); ++j)
                        anElement.addNodeId(aFace.nodeId(j));

                    aBoundaryMesh.addElement(anElementId, anElement);
                    anElementId++;
                }
            }

            for (Integer i = 0; i < static_cast<Integer>(m_nodeIds.size()); ++i)
            {
                Integer aNodeId = m_nodeIds.at(i);

                CMshNode<Real> aNode = m_nodes.at(aNodeId);

                aBoundaryMesh.addNode(aNodeId, aNode);
            }

            return aBoundaryMesh;
        }

        template <typename Real>
        void CMshMesh<Real>::setDx(const Real aValue) { m_dx = aValue; }

        template <typename Real>
        Real CMshMesh<Real>::dx() const { return m_dx; }

        template <typename Real>
        void CMshMesh<Real>::setDy(const Real aValue) { m_dy = aValue; }

        template <typename Real>
        Real CMshMesh<Real>::dy() const { return m_dy; }

        template <typename Real>
        void CMshMesh<Real>::setDz(const Real aValue) { m_dz = aValue; }

        template <typename Real>
        Real CMshMesh<Real>::dz() const { return m_dz; }

        template <typename Real>
        void CMshMesh<Real>::calculateFaceCentroid()
        {
            m_faceCentroid.clear();

            for (Integer i = 0; i < nbFaces(); ++i)
            {
                Integer aFaceId = faceId(i);

                CMshFace<Real> aFace = face(aFaceId);

                CGeoCoordinate<Real> aCentroid;

                aCentroid << 0.0, 0.0, 0.0;

                for (Integer j = 0; j < aFace.nbNodeIds(); ++j)
                {
                    aCentroid += m_nodes.at(aFace.nodeId(j));
                }

                if (aFace.nbNodeIds() > 0)
                    aCentroid /= static_cast<Real>(aFace.nbNodeIds());

                m_faceCentroid.push_back(aCentroid);
            }
        }

        template <typename Real>
        void CMshMesh<Real>::calculateElementCentroid()
        {
            m_elementCentroid.clear();

            for (Integer i = 0; i < nbElements(); ++i)
            {
                Integer anElementId = m_elementIds.at(i);

                CMshElement<Real>& anElement = m_elements.at(anElementId);

                CGeoCoordinate<Real> aCentroid;

                aCentroid << 0.0, 0.0, 0.0;

                for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                {
                    aCentroid += m_nodes.at(anElement.nodeId(j));
                }

                if (anElement.nbNodeIds() > 0)
                    aCentroid /= static_cast<Real>(anElement.nbNodeIds());

                m_elementCentroid.push_back(aCentroid);
            }
        }

        template <typename Real>
        const CGeoCoordinate<Real>& CMshMesh<Real>::faceCentroid(const Integer aFaceId) const { return m_faceCentroid.at(faceIndex(aFaceId)); }

        template <typename Real>
        const CGeoCoordinate<Real>& CMshMesh<Real>::elementCentroid(const Integer anElementId) const
        {
            return m_elementCentroid.at(m_elementIndices.at(anElementId));
        }

        template <typename Real>
        void CMshMesh<Real>::scale(const Real aFactor)
        {
            for (typename mapNode::iterator it = m_nodes.begin(); it != m_nodes.end(); ++it)
            {
                it->second *= aFactor;
            }
        }

        template <typename Real>
        Integer CMshMesh<Real>::nextNodeId() const
        {
            if (m_nodes.size() > 0)
                return m_nodes.rbegin()->first + 1;
            else
                return 0;
        }

        template <typename Real>
        Integer CMshMesh<Real>::nextFaceId() const
        {
            if (m_faces.size() > 0)
                return m_faces.rbegin()->first + 1;
            else
                return 0;
        }

        template <typename Real>
        Integer CMshMesh<Real>::nextElementId() const
        {
            if (m_elements.size() > 0)
                return m_elements.rbegin()->first + 1;
            else
                return 0;
        }

        template <typename Real>
        void CMshMesh<Real>::mergeNodes(const Real aTolerance)
        {
            CGeoHashGrid<Real> aHashGrid;

            for (Integer i = 0; i < static_cast<Integer>(m_nodeIds.size()); ++i)
            {
                Integer aNodeId = m_nodeIds.at(i);

                CMshNode<Real>& aNode = m_nodes.at(aNodeId);

                aHashGrid.addGeometricObject(aNodeId, aNode);
            }

            aHashGrid.build();

            std::map<Integer, bool> bDeleteNode;
            mapNodeIndex sNewNodeIds;

            for (Integer i = 0; i < static_cast<Integer>(m_nodeIds.size()); ++i)
            {
                Integer aNodeId = m_nodeIds.at(i);

                bDeleteNode[aNodeId] = false;
                sNewNodeIds[aNodeId] = aNodeId;
            }

            for (Integer i = 0; i < static_cast<Integer>(m_nodeIds.size()); ++i)
            {
                Integer aNodeId = m_nodeIds.at(i);

                if (bDeleteNode.at(aNodeId))
                    continue;

                CMshNode<Real>& aNode = m_nodes.at(aNodeId);

                std::vector<Integer> sNodes;

                aHashGrid.find(sNodes, aNode, aTolerance);

                if (sNodes.size() > 1)
                {
                    for (Integer k = 0; k < static_cast<Integer>(sNodes.size()); ++k)
                    {
                        if (sNodes.at(k) > aNodeId)
                        {
                            // Duplicate node found
                            bDeleteNode[sNodes.at(k)] = true;
                            sNewNodeIds[sNodes.at(k)] = aNodeId;
                        }
                    }
                }
            }

            // Renumber elements
            for (Integer i = 0; i < static_cast<Integer>(m_elementIds.size()); ++i)
            {
                Integer anElementId = m_elementIds.at(i);

                CMshElement<Real>& anElement = m_elements.at(anElementId);

                for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                {
                    Integer aNodeId = anElement.nodeId(j);

                    // Set to new node id
                    anElement.setNodeId(j, sNewNodeIds.at(aNodeId));
                }
            }

            // Renumber faces
            for (Integer i = 0; i < static_cast<Integer>(m_faceIds.size()); ++i)
            {
                Integer aFaceId = m_faceIds.at(i);

                CMshFace<Real>& aFace = m_faces.at(aFaceId);

                for (Integer j = 0; j < aFace.nbNodeIds(); ++j)
                {
                    Integer aNodeId = aFace.nodeId(j);

                    // Set to new node id
                    aFace.setNodeId(j, sNewNodeIds.at(aNodeId));
                }
            }

            // Delete nodes
            for (typename std::map<Integer, bool>::const_iterator itr = bDeleteNode.begin(); itr != bDeleteNode.end(); ++itr)
            {
                Integer aNodeId = itr->first;

                if (itr->second)
                    removeNode(aNodeId);
            }

            removeInvalidElements();

            rebuildIndices();
        }

        template <typename Real>
        void CMshMesh<Real>::removeInvalidElements()
        {
            std::map<Integer, bool> bDeleteElement;

            for (Integer i = 0; i < static_cast<Integer>(m_elementIds.size()); ++i)
            {
                Integer anElementId = m_elementIds.at(i);

                bDeleteElement[anElementId] = false;

                CMshElement<Real>& anElement = m_elements.at(anElementId);

                std::vector<Integer> sNodeIds;

                for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                    sNodeIds.push_back(anElement.nodeId(j));

                std::sort(sNodeIds.begin(), sNodeIds.end());

                if (std::adjacent_find(sNodeIds.begin(), sNodeIds.end()) != sNodeIds.end())
                    bDeleteElement[anElementId] = true;
            }

            // Delete elements
            for (typename std::map<Integer, bool>::const_iterator itr = bDeleteElement.begin(); itr != bDeleteElement.end(); ++itr)
            {
                Integer anElementId = itr->first;

                if (itr->second)
                    removeElement(anElementId);
            }
        }

        template <typename Real>
        void CMshMesh<Real>::rebuildIndices()
        {
            // Rebuild node indices
            m_nodeIndices.clear();

            for (Integer i = 0; i < static_cast<Integer>(m_nodeIds.size()); ++i)
            {
                Integer aNodeId = m_nodeIds.at(i);

                m_nodeIndices[aNodeId] = i;
            }

            // Rebuild face indices
            m_faceIndices.clear();

            for (Integer i = 0; i < static_cast<Integer>(m_faceIds.size()); ++i)
            {
                Integer aFaceId = m_faceIds.at(i);

                m_faceIndices[aFaceId] = i;
            }

            // Rebuild element indices
            m_elementIndices.clear();

            for (Integer i = 0; i < static_cast<Integer>(m_elementIds.size()); ++i)
            {
                Integer anElementId = m_elementIds.at(i);

                m_elementIndices[anElementId] = i;
            }
        }

        template <typename Real>
        void CMshMesh<Real>::removeDanglingNodes()
        {
            std::map<Integer, bool> bDeleteNode;

            for (Integer i = 0; i < static_cast<Integer>(m_nodeIds.size()); ++i)
            {
                Integer aNodeId = m_nodeIds.at(i);

                bDeleteNode[aNodeId] = true;
            }

            for (Integer i = 0; i < static_cast<Integer>(m_faceIds.size()); ++i)
            {
                Integer aFaceId = m_faceIds.at(i);

                CMshFace<Real>& aFace = m_faces.at(aFaceId);

                for (Integer k = 0; k < aFace.nbNodeIds(); ++k)
                {
                    Integer aNodeId = aFace.nodeId(k);
                    bDeleteNode[aNodeId] = false;
                }
            }

            for (Integer i = 0; i < static_cast<Integer>(m_elementIds.size()); ++i)
            {
                Integer anElementId = m_elementIds.at(i);

                CMshElement<Real>& anElement = m_elements.at(anElementId);

                for (Integer k = 0; k < anElement.nbNodeIds(); ++k)
                {
                    Integer aNodeId = anElement.nodeId(k);
                    bDeleteNode[aNodeId] = false;
                }
            }

            // Delete nodes
            for (typename std::map<Integer, bool>::const_iterator itr = bDeleteNode.begin(); itr != bDeleteNode.end(); ++itr)
            {
                Integer aNodeId = itr->first;

                if (itr->second)
                    removeNode(aNodeId);
            }
        }

        template <typename Real>
        void CMshMesh<Real>::collapseNakedEdges(const Real aTolerance)
        {
            Real minDistance = std::numeric_limits<Real>::max();

            // Get naked edges
            std::vector<Integer> sFaces;

            for (Integer i = 0; i < static_cast<Integer>(m_faceIds.size()); ++i)
            {
                Integer aFaceId = m_faceIds.at(i);

                CMshFace<Real>& aFace = m_faces.at(aFaceId);

                if (aFace.faceType() != FT_LINE)
                    continue;

                if (!aFace.hasPair())
                {
                    sFaces.push_back(aFaceId);

                    Integer aNodeId1 = aFace.nodeId(0);
                    Integer aNodeId2 = aFace.nodeId(1);

                    CMshNode<Real>& aNode1 = m_nodes.at(aNodeId1);
                    CMshNode<Real>& aNode2 = m_nodes.at(aNodeId2);

                    Real aDistance = (aNode2 - aNode1).norm();

                    minDistance = std::min(minDistance, aDistance);
                }
            }

            minDistance += aTolerance;

            for (Integer i = 0; i < static_cast<Integer>(sFaces.size()); ++i)
            {
                Integer aFaceId = sFaces.at(i);

                CMshFace<Real>& aFace = m_faces.at(aFaceId);

                Integer aNodeId1 = aFace.nodeId(0);
                Integer aNodeId2 = aFace.nodeId(1);

                CMshNode<Real>& aNode1 = m_nodes.at(aNodeId1);
                CMshNode<Real>& aNode2 = m_nodes.at(aNodeId2);

                Real aDistance = (aNode2 - aNode1).norm();

                if (aDistance < minDistance)
                {
                    aNode1 = aNode2;
                }
            }
        }

        template <typename Real>
        void CMshMesh<Real>::renumber()
        {
            mapNodeIndex sNewNodeIds;
            mapNode sNewNodes;

            for (Integer i = 0; i < static_cast<Integer>(m_nodeIds.size()); ++i)
            {

                Integer anOldNodeId = m_nodeIds.at(i);
                Integer aNewNodeId = i;

                sNewNodeIds[anOldNodeId] = aNewNodeId;

                CMshNode<Real> aNode = m_nodes.at(anOldNodeId);

                sNewNodes[aNewNodeId] = aNode;
            }

            mapNodeIndex sNewFaceIds;
            mapFace sNewFaces;

            for (Integer i = 0; i < static_cast<Integer>(m_faceIds.size()); ++i)
            {

                Integer anOldFaceId = m_faceIds.at(i);
                Integer aNewFaceId = i;

                sNewFaceIds[anOldFaceId] = aNewFaceId;

                CMshFace<Real> aFace = m_faces.at(anOldFaceId);

                for (Integer j = 0; j < aFace.nbNodeIds(); ++j)
                    aFace.setNodeId(j, sNewNodeIds.at(aFace.nodeId(j)));

                sNewFaces[aNewFaceId] = aFace;
            }

            mapElement sNewElements;

            for (Integer i = 0; i < static_cast<Integer>(m_elementIds.size()); ++i)
            {

                Integer anOldElementId = m_elementIds.at(i);
                Integer aNewElementId = i;

                CMshElement<Real>& anElement = m_elements.at(anOldElementId);

                for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                    anElement.setNodeId(j, sNewNodeIds.at(anElement.nodeId(j)));

                for (Integer j = 0; j < anElement.nbFaceIds(); ++j)
                    anElement.setFaceId(j, sNewFaceIds.at(anElement.faceId(j)));

                sNewElements[aNewElementId] = anElement;
            }

            Integer n = static_cast<Integer>(m_nodeIds.size());
            m_nodeIds.clear();
            for (Integer i = 0; i < n; ++i)
                m_nodeIds.push_back(i);

            Integer f = static_cast<Integer>(m_faceIds.size());
            m_faceIds.clear();
            for (Integer i = 0; i < f; ++i)
                m_faceIds.push_back(i);

            Integer e = static_cast<Integer>(m_elementIds.size());
            m_elementIds.clear();
            for (Integer i = 0; i < e; ++i)
                m_elementIds.push_back(i);

            m_nodes = sNewNodes;
            m_faces = sNewFaces;
            m_elements = sNewElements;

            rebuildIndices();
        }

        template <typename Real>
        void CMshMesh<Real>::invert()
        {
            for (Integer i = 0; i < static_cast<Integer>(m_elementIds.size()); ++i)
            {
                Integer anElementId = m_elementIds.at(i);

                m_elements.at(anElementId).invert();
            }
        }

        template <typename Real>
        Real CMshMesh<Real>::quality(Integer anElementId)
        {
            Real q = 0.0;

            CMshElement<Real>& anElement = m_elements.at(anElementId);

            if (anElement.elementType() == ET_TRIANGLE)
            {
                CMshTriangle<Real> aTriangle;

                for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                {
                    Integer aNodeId = anElement.nodeId(j);
                    CMshNode<Real>& aNode = m_nodes.at(aNodeId);

                    aTriangle.addVertex(aNode);
                }

                aTriangle.calculateQuality();

                q = aTriangle.quality();
            }
            else if (anElement.elementType() == ET_QUADRILATERAL)
            {
                CMshQuadrilateral<Real> aQuadrilateral;

                for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                {
                    Integer aNodeId = anElement.nodeId(j);
                    CMshNode<Real>& aNode = m_nodes.at(aNodeId);

                    aQuadrilateral.addVertex(aNode);
                }

                aQuadrilateral.calculateQuality();

                q = aQuadrilateral.quality();
            }
            else if (anElement.elementType() == ET_TETRAHEDRON)
            {
                CMshTetrahedron<Real> aTetrahedron;

                for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                {
                    Integer aNodeId = anElement.nodeId(j);
                    CMshNode<Real>& aNode = m_nodes.at(aNodeId);

                    aTetrahedron.addVertex(aNode);
                }

                aTetrahedron.calculateQuality();

                q = aTetrahedron.quality();
            }

            return q;
        }

        template <typename Real>
        void CMshMesh<Real>::meshQuality(Real& aMinQ, Real& aMaxQ, Real& aAveQ)
        {
            aMinQ = 1.0;
            aMaxQ = 0.0;
            aAveQ = 0.0;

            Real aSumQ = 0.0;

            Integer n = 0;

            for (Integer i = 0; i < nbElements(); ++i)
            {
                Integer anElementId = elementId(i);

                q = this->quality(anElementId);

                aSumQ += q;

                aMinQ = std::min(q, aMinQ);
                aMaxQ = std::max(q, aMaxQ);

                n++;
            }

            if (n > 0)
                aAveQ = aSumQ / n;
        }

        template <typename Real>
        CGeoBoundingBox<Real> CMshMesh<Real>::boundingBox(const Integer anElementId)
        {
            CGeoBoundingBox<Real> aBoundingBox;

            for (Integer j = 0; j < m_elements.at(anElementId).nbNodeIds(); ++j)
            {
                Integer aNodeId = m_elements.at(anElementId).nodeId(j);
                CMshNode<Real>& aNode = m_nodes.at(aNodeId);

                aBoundingBox.addCoordinate(aNode);
            }

            return aBoundingBox;
        }

        template <typename Real>
        CGeoBoundingBox<Real> CMshMesh<Real>::boundingBox()
        {
            CGeoBoundingBox<Real> aBoundingBox;

            for (Integer i = 0; i < static_cast<Integer>(m_nodes.size()); ++i)
            {
                Integer aNodeId = m_nodeIds.at(i);
                CMshNode<Real>& aNode = m_nodes.at(aNodeId);

                aBoundingBox.addCoordinate(aNode);
            }

            return aBoundingBox;
        }
    }
}
