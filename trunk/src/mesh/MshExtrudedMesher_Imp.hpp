// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#include "MshMeshQuery.hpp"

namespace ENigMA
{
    namespace mesh
    {
        template <typename Real>
        CMshExtrudedMesher<Real>::CMshExtrudedMesher()
        {
        }

        template <typename Real>
        CMshExtrudedMesher<Real>::~CMshExtrudedMesher()
        {
        }

        template <typename Real>
        bool CMshExtrudedMesher<Real>::generate(CMshMesh<Real>& aSurfaceMesh, const Real aDelta, const Real aTolerance)
        {
            CMshMeshQuery<Real> aMeshQuery(aSurfaceMesh);
            std::vector<Integer> sElementIds;

            CGeoNormal<Real> aNormal;
            std::map<Integer, CGeoNormal<Real>> sNormals;
            CGeoTriangle<Real> aTriangle;

            for (Integer i = 0; i < aSurfaceMesh.nbNodes(); i++)
            {
                Integer aNodeId = aSurfaceMesh.nodeId(i);
                aMeshQuery.elementsSharingNode(aNodeId, sElementIds);

                aNormal << 0, 0, 0;
                for (Integer k = 0; k < sElementIds.size(); k++)
                {
                    CMshElement<Real> anElement = aSurfaceMesh.element(sElementIds[k]);

                    aTriangle.reset();
                    for (Integer j = 0; j < anElement.nbNodeIds(); j++)
                    {
                        aTriangle.addVertex(aSurfaceMesh.node(anElement.nodeId(j)));
                    }
                    aTriangle.calculateNormal(true);

                    aNormal += aTriangle.normal();
                }
                aNormal.normalize();

                sNormals[aNodeId] = -aNormal;
            }

            // Extrude elements
            for (Integer i = 0; i < aSurfaceMesh.nbElements(); i++)
            {
                Integer anElementId = aSurfaceMesh.elementId(i);
                CMshElement<Real> anElement = aSurfaceMesh.element(anElementId);

                if (anElement.elementType() == ET_TRIANGLE)
                {
                    // Bottom
                    CMshNode<Real> aNode1 = aSurfaceMesh.node(anElement.nodeId(0));
                    Integer aNodeId1 = m_mesh.nextNodeId();
                    m_mesh.addNode(aNodeId1, aNode1);

                    CMshNode<Real> aNode2 = aSurfaceMesh.node(anElement.nodeId(1));
                    Integer aNodeId2 = m_mesh.nextNodeId();
                    m_mesh.addNode(aNodeId2, aNode2);

                    CMshNode<Real> aNode3 = aSurfaceMesh.node(anElement.nodeId(2));
                    Integer aNodeId3 = m_mesh.nextNodeId();
                    m_mesh.addNode(aNodeId3, aNode3);

                    // Extrude
                    aNode1 += sNormals.at(anElement.nodeId(0)) * aDelta;
                    aNode2 += sNormals.at(anElement.nodeId(1)) * aDelta;
                    aNode3 += sNormals.at(anElement.nodeId(2)) * aDelta;

                    // Top
                    Integer aNodeId4 = m_mesh.nextNodeId();
                    m_mesh.addNode(aNodeId4, aNode1);

                    Integer aNodeId5 = m_mesh.nextNodeId();
                    m_mesh.addNode(aNodeId5, aNode2);

                    Integer aNodeId6 = m_mesh.nextNodeId();
                    m_mesh.addNode(aNodeId6, aNode3);

                    CMshElement<Real> aSolidElement(ET_TRIANGULAR_PRISM);

                    aSolidElement.addNodeId(aNodeId1);
                    aSolidElement.addNodeId(aNodeId2);
                    aSolidElement.addNodeId(aNodeId3);
                    aSolidElement.addNodeId(aNodeId4);
                    aSolidElement.addNodeId(aNodeId5);
                    aSolidElement.addNodeId(aNodeId6);

                    m_mesh.addElement(m_mesh.nextElementId(), aSolidElement);
                }
                else if (anElement.elementType() == ET_QUADRILATERAL)
                {
                    // Bottom
                    CMshNode<Real> aNode1 = aSurfaceMesh.node(anElement.nodeId(0));
                    Integer aNodeId1 = m_mesh.nextNodeId();
                    m_mesh.addNode(aNodeId1, aNode1);

                    CMshNode<Real> aNode2 = aSurfaceMesh.node(anElement.nodeId(1));
                    Integer aNodeId2 = m_mesh.nextNodeId();
                    m_mesh.addNode(aNodeId2, aNode2);

                    CMshNode<Real> aNode3 = aSurfaceMesh.node(anElement.nodeId(2));
                    Integer aNodeId3 = m_mesh.nextNodeId();
                    m_mesh.addNode(aNodeId3, aNode3);

                    CMshNode<Real> aNode4 = aSurfaceMesh.node(anElement.nodeId(3));
                    Integer aNodeId4 = m_mesh.nextNodeId();
                    m_mesh.addNode(aNodeId4, aNode4);

                    // Extrude
                    aNode1 += sNormals.at(anElement.nodeId(0)) * aDelta;
                    aNode2 += sNormals.at(anElement.nodeId(1)) * aDelta;
                    aNode3 += sNormals.at(anElement.nodeId(2)) * aDelta;
                    aNode4 += sNormals.at(anElement.nodeId(3)) * aDelta;

                    // Top
                    Integer aNodeId5 = m_mesh.nextNodeId();
                    m_mesh.addNode(aNodeId5, aNode1);

                    Integer aNodeId6 = m_mesh.nextNodeId();
                    m_mesh.addNode(aNodeId6, aNode2);

                    Integer aNodeId7 = m_mesh.nextNodeId();
                    m_mesh.addNode(aNodeId7, aNode3);

                    Integer aNodeId8 = m_mesh.nextNodeId();
                    m_mesh.addNode(aNodeId8, aNode4);

                    CMshElement<Real> aSolidElement(ET_HEXAHEDRON);

                    aSolidElement.addNodeId(aNodeId1);
                    aSolidElement.addNodeId(aNodeId2);
                    aSolidElement.addNodeId(aNodeId3);
                    aSolidElement.addNodeId(aNodeId4);
                    aSolidElement.addNodeId(aNodeId5);
                    aSolidElement.addNodeId(aNodeId8);
                    aSolidElement.addNodeId(aNodeId7);
                    aSolidElement.addNodeId(aNodeId6);

                    m_mesh.addElement(m_mesh.nextElementId(), aSolidElement);
                }
            }

            // Update mesh
            for (Integer i = 0; i < aSurfaceMesh.nbNodes(); i++)
            {
                Integer aNodeId = aSurfaceMesh.nodeId(i);
                CMshNode<Real>& aNode = aSurfaceMesh.node(aNodeId);

                aNode += sNormals.at(aNodeId) * aDelta;
            }

            return true;
        }

        template <typename Real>
        CMshMesh<Real>& CMshExtrudedMesher<Real>::mesh()
        {
            return m_mesh;
        }
    }
}
