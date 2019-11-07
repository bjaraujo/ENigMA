// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

namespace ENigMA {

namespace mesh {

    template <typename Real>
    CMshExtrudedMesher<Real>::CMshExtrudedMesher()
    {
    }

    template <typename Real>
    CMshExtrudedMesher<Real>::~CMshExtrudedMesher()
    {
    }

    template <typename Real>
    bool CMshExtrudedMesher<Real>::generate(CMshMesh<Real>& aPlanarMesh, const Integer nw, Real dw, const Real aTolerance)
    {

        for (Integer k = 0; k < nw; k++) {

            // Extrude elements
            for (Integer i = 0; i < aPlanarMesh.nbElements(); i++) {
                Integer anElementId = aPlanarMesh.elementId(i);
                CMshElement<Real>& aPlanarElement = aPlanarMesh.element(anElementId);

                if (aPlanarElement.elementType() == ET_TRIANGLE) {

                    // Bottom
                    CMshNode<Real> aNode1 = aPlanarMesh.node(aPlanarElement.nodeId(0));
                    Integer aNodeId1 = m_mesh.nextNodeId();
                    m_mesh.addNode(aNodeId1, aNode1);
                    m_mesh.node(aNodeId1).z() += (k + 0) * dw;

                    CMshNode<Real> aNode2 = aPlanarMesh.node(aPlanarElement.nodeId(1));
                    Integer aNodeId2 = m_mesh.nextNodeId();
                    m_mesh.addNode(aNodeId2, aNode2);
                    m_mesh.node(aNodeId2).z() += (k + 0) * dw;

                    CMshNode<Real> aNode3 = aPlanarMesh.node(aPlanarElement.nodeId(2));
                    Integer aNodeId3 = m_mesh.nextNodeId();
                    m_mesh.addNode(aNodeId3, aNode3);
                    m_mesh.node(aNodeId3).z() += (k + 0) * dw;

                    // Top
                    CMshNode<Real> aNode4 = aNode1;
                    Integer aNodeId4 = m_mesh.nextNodeId();
                    m_mesh.addNode(aNodeId4, aNode4);
                    m_mesh.node(aNodeId4).z() += (k + 1) * dw;

                    CMshNode<Real> aNode5 = aNode2;
                    Integer aNodeId5 = m_mesh.nextNodeId();
                    m_mesh.addNode(aNodeId5, aNode5);
                    m_mesh.node(aNodeId5).z() += (k + 1) * dw;

                    CMshNode<Real> aNode6 = aNode3;
                    Integer aNodeId6 = m_mesh.nextNodeId();
                    m_mesh.addNode(aNodeId6, aNode6);
                    m_mesh.node(aNodeId6).z() += (k + 1) * dw;

                    CMshElement<Real> aSolidElement(ET_TRIANGULAR_PRISM);

                    aSolidElement.addNodeId(aNodeId1);
                    aSolidElement.addNodeId(aNodeId2);
                    aSolidElement.addNodeId(aNodeId3);
                    aSolidElement.addNodeId(aNodeId4);
                    aSolidElement.addNodeId(aNodeId5);
                    aSolidElement.addNodeId(aNodeId6);

                    m_mesh.addElement(m_mesh.nextElementId(), aSolidElement);

                } else if (aPlanarElement.elementType() == ET_QUADRILATERAL) {

                    // Bottom
                    CMshNode<Real> aNode1 = aPlanarMesh.node(aPlanarElement.nodeId(0));
                    Integer aNodeId1 = m_mesh.nextNodeId();
                    m_mesh.addNode(aNodeId1, aNode1);
                    m_mesh.node(aNodeId1).z() += (k + 0) * dw;

                    CMshNode<Real> aNode2 = aPlanarMesh.node(aPlanarElement.nodeId(1));
                    Integer aNodeId2 = m_mesh.nextNodeId();
                    m_mesh.addNode(aNodeId2, aNode2);
                    m_mesh.node(aNodeId2).z() += (k + 0) * dw;

                    CMshNode<Real> aNode3 = aPlanarMesh.node(aPlanarElement.nodeId(2));
                    Integer aNodeId3 = m_mesh.nextNodeId();
                    m_mesh.addNode(aNodeId3, aNode3);
                    m_mesh.node(aNodeId3).z() += (k + 0) * dw;

                    CMshNode<Real> aNode4 = aPlanarMesh.node(aPlanarElement.nodeId(3));
                    Integer aNodeId4 = m_mesh.nextNodeId();
                    m_mesh.addNode(aNodeId4, aNode4);
                    m_mesh.node(aNodeId4).z() += (k + 0) * dw;

                    // Top
                    CMshNode<Real> aNode5 = aNode1;
                    Integer aNodeId5 = m_mesh.nextNodeId();
                    m_mesh.addNode(aNodeId5, aNode5);
                    m_mesh.node(aNodeId5).z() += (k + 1) * dw;

                    CMshNode<Real> aNode6 = aNode2;
                    Integer aNodeId6 = m_mesh.nextNodeId();
                    m_mesh.addNode(aNodeId6, aNode6);
                    m_mesh.node(aNodeId6).z() += (k + 1) * dw;

                    CMshNode<Real> aNode7 = aNode3;
                    Integer aNodeId7 = m_mesh.nextNodeId();
                    m_mesh.addNode(aNodeId7, aNode7);
                    m_mesh.node(aNodeId7).z() += (k + 1) * dw;

                    CMshNode<Real> aNode8 = aNode4;
                    Integer aNodeId8 = m_mesh.nextNodeId();
                    m_mesh.addNode(aNodeId8, aNode8);
                    m_mesh.node(aNodeId8).z() += (k + 1) * dw;

                    CMshElement<Real> aSolidElement(ET_HEXAHEDRON);

                    aSolidElement.addNodeId(aNodeId1);
                    aSolidElement.addNodeId(aNodeId2);
                    aSolidElement.addNodeId(aNodeId3);
                    aSolidElement.addNodeId(aNodeId4);
                    aSolidElement.addNodeId(aNodeId5);
                    aSolidElement.addNodeId(aNodeId6);
                    aSolidElement.addNodeId(aNodeId7);
                    aSolidElement.addNodeId(aNodeId8);

                    m_mesh.addElement(m_mesh.nextElementId(), aSolidElement);
                }
            }

            // Connect elements
            m_mesh.mergeNodes(aTolerance);
            m_mesh.renumber();
            m_mesh.generateFaces(aTolerance);
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
