// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "GeoNormal.hpp"

using namespace ENigMA::geometry;

namespace ENigMA
{
    namespace mesh
    {
        template <typename Real>
        CMshBasicMesher<Real>::CMshBasicMesher()
        {
        }

        template <typename Real>
        CMshBasicMesher<Real>::~CMshBasicMesher()
        {
        }

        template <typename Real>
        bool CMshBasicMesher<Real>::generate(const CGeoLine<Real>& aLine, const Integer nu)
        {
            Real du;

            CGeoVector<Real> aVectorU = aLine.vector();
            CGeoNormal<Real> aNormalU = aVectorU;
            aNormalU.normalize();

            if (nu > 0)
                du = aVectorU.norm() / nu;
            else
                du = 0.0;

            m_mesh.reset();

            m_mesh.setDx(du);

            Integer nn = 0;

            for (Integer i = 0; i < nu + 1; ++i)
            {
                CMshNode<Real> aNode;

                aNode = aLine.startPoint() + aNormalU * du * static_cast<Real>(i);

                m_mesh.addNode(nn, aNode);
                nn++;
            }

            Integer ne = 0;

            for (Integer i = 0; i < nu; ++i)
            {
                CMshElement<Real> anElement;

                anElement.addNodeId(i);
                anElement.addNodeId(i + 1);

                anElement.setElementType(ET_BEAM);

                m_mesh.addElement(ne, anElement);
                ne++;
            }

            return true;
        }

        template <typename Real>
        bool CMshBasicMesher<Real>::generate(const CGeoQuadrilateral<Real>& aQuadrilateral, const Integer nu, const Integer nv, bool decimate)
        {
            Real du, dv;

            CGeoVector<Real> aVectorU = aQuadrilateral.vertex(1) - aQuadrilateral.vertex(0);
            CGeoNormal<Real> aNormalU = aVectorU;
            aNormalU.normalize();

            CGeoVector<Real> aVectorV = aQuadrilateral.vertex(3) - aQuadrilateral.vertex(0);
            CGeoNormal<Real> aNormalV = aVectorV;
            aNormalV.normalize();

            if (nu > 0)
                du = aVectorU.norm() / nu;
            else
                du = 0.0;

            if (nv > 0)
                dv = aVectorV.norm() / nv;
            else
                dv = 0.0;

            m_mesh.reset();

            m_mesh.setDx(du);
            m_mesh.setDy(dv);

            Integer nn = 0;

            for (Integer j = 0; j < nv + 1; ++j)
            {
                for (Integer i = 0; i < nu + 1; ++i)
                {
                    CMshNode<Real> aNode;

                    aNode = aQuadrilateral.vertex(0) + aNormalU * du * static_cast<Real>(i) + aNormalV * dv * static_cast<Real>(j);

                    m_mesh.addNode(nn, aNode);
                    nn++;
                }
            }

            Integer ne = 0;

            for (Integer j = 0; j < nv; ++j)
            {
                for (Integer i = 0; i < nu; ++i)
                {
                    if (decimate)
                    {
                        CMshElement<Real> anElement;

                        // Triangle 1
                        anElement.reset();
                        anElement.addNodeId((nu + 1) * j + i);
                        anElement.addNodeId((nu + 1) * j + i + 1);
                        anElement.addNodeId((nu + 1) * (j + 1) + i + 1);

                        anElement.setElementType(ET_TRIANGLE);

                        m_mesh.addElement(ne, anElement);
                        ne++;

                        // Triangle 2
                        anElement.reset();
                        anElement.addNodeId((nu + 1) * j + i);
                        anElement.addNodeId((nu + 1) * (j + 1) + i + 1);
                        anElement.addNodeId((nu + 1) * (j + 1) + i);

                        anElement.setElementType(ET_TRIANGLE);

                        m_mesh.addElement(ne, anElement);
                        ne++;
                    }
                    else
                    {
                        // Quadrilateral
                        CMshElement<Real> anElement;

                        anElement.addNodeId((nu + 1) * j + i);
                        anElement.addNodeId((nu + 1) * j + i + 1);
                        anElement.addNodeId((nu + 1) * (j + 1) + i + 1);
                        anElement.addNodeId((nu + 1) * (j + 1) + i);

                        anElement.setElementType(ET_QUADRILATERAL);

                        m_mesh.addElement(ne, anElement);
                        ne++;
                    }
                }
            }

            return true;
        }

        template <typename Real>
        bool CMshBasicMesher<Real>::generate(const CGeoHexahedron<Real>& aHexahedron, const Integer nu, const Integer nv, const Integer nw, bool decimate)
        {
            CGeoVector<Real> aVectorU = (aHexahedron.vertex(1) - aHexahedron.vertex(0));
            CGeoNormal<Real> aNormalU = aVectorU;
            aNormalU.normalize();

            CGeoVector<Real> aVectorV = aHexahedron.vertex(3) - aHexahedron.vertex(0);
            CGeoNormal<Real> aNormalV = aVectorV;
            aNormalV.normalize();

            CGeoVector<Real> aVectorW = aHexahedron.vertex(4) - aHexahedron.vertex(0);
            CGeoNormal<Real> aNormalW = aVectorW;
            aNormalW.normalize();

            Real du, dv, dw;

            du = dv = dw = 0.0;

            if (nu > 0)
                du = aVectorU.norm() / nu;

            if (nv > 0)
                dv = aVectorV.norm() / nv;

            if (nw > 0)
                dw = aVectorW.norm() / nw;

            m_mesh.reset();

            m_mesh.setDx(du);
            m_mesh.setDy(dv);
            m_mesh.setDz(dw);

            Integer nn = 0;

            for (Integer k = 0; k < nw + 1; ++k)
            {
                for (Integer j = 0; j < nv + 1; ++j)
                {
                    for (Integer i = 0; i < nu + 1; ++i)
                    {
                        CMshNode<Real> aNode;

                        aNode = aHexahedron.vertex(0) + aNormalU * du * static_cast<Real>(i) + aNormalV * dv * static_cast<Real>(j) + aNormalW * dw * static_cast<Real>(k);

                        m_mesh.addNode(nn, aNode);
                        nn++;
                    }
                }
            }

            Integer ne = 0;

            for (Integer k = 0; k < nw; ++k)
            {
                for (Integer j = 0; j < nv; ++j)
                {
                    for (Integer i = 0; i < nu; ++i)
                    {
                        Integer nodeId[8];

                        nodeId[0] = (nu + 1) * (nv + 1) * k + (nu + 1) * j + i;
                        nodeId[1] = (nu + 1) * (nv + 1) * k + (nu + 1) * j + i + 1;
                        nodeId[2] = (nu + 1) * (nv + 1) * k + (nu + 1) * (j + 1) + i + 1;
                        nodeId[3] = (nu + 1) * (nv + 1) * k + (nu + 1) * (j + 1) + i;
                        nodeId[4] = (nu + 1) * (nv + 1) * (k + 1) + (nu + 1) * j + i;
                        nodeId[5] = (nu + 1) * (nv + 1) * (k + 1) + (nu + 1) * j + i + 1;
                        nodeId[6] = (nu + 1) * (nv + 1) * (k + 1) + (nu + 1) * (j + 1) + i + 1;
                        nodeId[7] = (nu + 1) * (nv + 1) * (k + 1) + (nu + 1) * (j + 1) + i;

                        if (decimate)
                        {
                            CMshElement<Real> anElement;

                            // Tetrahedron 1
                            anElement.reset();
                            anElement.addNodeId(nodeId[1]);
                            anElement.addNodeId(nodeId[0]);
                            anElement.addNodeId(nodeId[3]);
                            anElement.addNodeId(nodeId[7]);

                            anElement.setElementType(ET_TETRAHEDRON);

                            m_mesh.addElement(ne, anElement);
                            ne++;

                            // Tetrahedron 2
                            anElement.reset();
                            anElement.addNodeId(nodeId[4]);
                            anElement.addNodeId(nodeId[0]);
                            anElement.addNodeId(nodeId[1]);
                            anElement.addNodeId(nodeId[7]);

                            anElement.setElementType(ET_TETRAHEDRON);

                            m_mesh.addElement(ne, anElement);
                            ne++;

                            // Tetrahedron 3
                            anElement.reset();
                            anElement.addNodeId(nodeId[4]);
                            anElement.addNodeId(nodeId[1]);
                            anElement.addNodeId(nodeId[5]);
                            anElement.addNodeId(nodeId[7]);

                            anElement.setElementType(ET_TETRAHEDRON);

                            m_mesh.addElement(ne, anElement);
                            ne++;

                            // Tetrahedron 4
                            anElement.reset();
                            anElement.addNodeId(nodeId[2]);
                            anElement.addNodeId(nodeId[1]);
                            anElement.addNodeId(nodeId[3]);
                            anElement.addNodeId(nodeId[7]);

                            anElement.setElementType(ET_TETRAHEDRON);

                            m_mesh.addElement(ne, anElement);
                            ne++;

                            // Tetrahedron 5
                            anElement.reset();
                            anElement.addNodeId(nodeId[6]);
                            anElement.addNodeId(nodeId[1]);
                            anElement.addNodeId(nodeId[2]);
                            anElement.addNodeId(nodeId[7]);

                            anElement.setElementType(ET_TETRAHEDRON);

                            m_mesh.addElement(ne, anElement);
                            ne++;

                            // Tetrahedron 6
                            anElement.reset();
                            anElement.addNodeId(nodeId[5]);
                            anElement.addNodeId(nodeId[1]);
                            anElement.addNodeId(nodeId[6]);
                            anElement.addNodeId(nodeId[7]);

                            anElement.setElementType(ET_TETRAHEDRON);

                            m_mesh.addElement(ne, anElement);
                            ne++;
                        }
                        else
                        {
                            // Hexahedron
                            CMshElement<Real> anElement;

                            anElement.addNodeId(nodeId[0]);
                            anElement.addNodeId(nodeId[1]);
                            anElement.addNodeId(nodeId[2]);
                            anElement.addNodeId(nodeId[3]);
                            anElement.addNodeId(nodeId[4]);
                            anElement.addNodeId(nodeId[5]);
                            anElement.addNodeId(nodeId[6]);
                            anElement.addNodeId(nodeId[7]);

                            anElement.setElementType(ET_HEXAHEDRON);

                            m_mesh.addElement(ne, anElement);
                            ne++;
                        }
                    }
                }
            }

            return true;
        }

        template <typename Real>
        bool CMshBasicMesher<Real>::generate(const CGeoBoundingBox<Real>& aBoundingBox, const Real meshSize, bool decimate)
        {
            CGeoVector<Real> aVector = aBoundingBox.max() - aBoundingBox.min();

            CGeoVector<Real> aVectorU = CGeoVector<Real>(aVector.x(), 0.0, 0.0);
            CGeoNormal<Real> aNormalU = aVectorU;
            aNormalU.normalize();

            CGeoVector<Real> aVectorV = CGeoVector<Real>(0.0, aVector.y(), 0.0);
            CGeoNormal<Real> aNormalV = aVectorV;
            aNormalV.normalize();

            CGeoVector<Real> aVectorW = CGeoVector<Real>(0.0, 0.0, aVector.z());
            CGeoNormal<Real> aNormalW = aVectorW;
            aNormalW.normalize();

            Integer nu, nv, nw;

            nu = static_cast<Integer>(aVectorU.norm() / meshSize);
            nv = static_cast<Integer>(aVectorV.norm() / meshSize);
            nw = static_cast<Integer>(aVectorW.norm() / meshSize);

            m_mesh.reset();

            m_mesh.setDx(meshSize);
            m_mesh.setDy(meshSize);
            m_mesh.setDz(meshSize);

            Integer nn = 0;

            for (Integer k = 0; k < nw + 1; ++k)
            {
                for (Integer j = 0; j < nv + 1; ++j)
                {
                    for (Integer i = 0; i < nu + 1; ++i)
                    {
                        CMshNode<Real> aNode;

                        aNode = aBoundingBox.min() + aNormalU * meshSize * static_cast<Real>(i) + aNormalV * meshSize * static_cast<Real>(j) + aNormalW * meshSize * static_cast<Real>(k);

                        m_mesh.addNode(nn, aNode);
                        nn++;
                    }
                }
            }

            Integer ne = 0;

            for (Integer k = 0; k < nw; ++k)
            {
                for (Integer j = 0; j < nv; ++j)
                {
                    for (Integer i = 0; i < nu; ++i)
                    {
                        Integer nodeId[8];

                        nodeId[0] = (nu + 1) * (nv + 1) * k + (nu + 1) * j + i;
                        nodeId[1] = (nu + 1) * (nv + 1) * k + (nu + 1) * j + i + 1;
                        nodeId[2] = (nu + 1) * (nv + 1) * k + (nu + 1) * (j + 1) + i + 1;
                        nodeId[3] = (nu + 1) * (nv + 1) * k + (nu + 1) * (j + 1) + i;
                        nodeId[4] = (nu + 1) * (nv + 1) * (k + 1) + (nu + 1) * j + i;
                        nodeId[5] = (nu + 1) * (nv + 1) * (k + 1) + (nu + 1) * j + i + 1;
                        nodeId[6] = (nu + 1) * (nv + 1) * (k + 1) + (nu + 1) * (j + 1) + i + 1;
                        nodeId[7] = (nu + 1) * (nv + 1) * (k + 1) + (nu + 1) * (j + 1) + i;

                        if (decimate)
                        {
                            CMshElement<Real> anElement;

                            // Tetrahedron 1
                            anElement.reset();
                            anElement.addNodeId(nodeId[0]);
                            anElement.addNodeId(nodeId[3]);
                            anElement.addNodeId(nodeId[1]);
                            anElement.addNodeId(nodeId[7]);

                            anElement.setElementType(ET_TETRAHEDRON);

                            m_mesh.addElement(ne, anElement);
                            ne++;

                            // Tetrahedron 2
                            anElement.reset();
                            anElement.addNodeId(nodeId[0]);
                            anElement.addNodeId(nodeId[1]);
                            anElement.addNodeId(nodeId[4]);
                            anElement.addNodeId(nodeId[7]);

                            anElement.setElementType(ET_TETRAHEDRON);

                            m_mesh.addElement(ne, anElement);
                            ne++;

                            // Tetrahedron 3
                            anElement.reset();
                            anElement.addNodeId(nodeId[1]);
                            anElement.addNodeId(nodeId[5]);
                            anElement.addNodeId(nodeId[4]);
                            anElement.addNodeId(nodeId[7]);

                            anElement.setElementType(ET_TETRAHEDRON);

                            m_mesh.addElement(ne, anElement);
                            ne++;

                            // Tetrahedron 4
                            anElement.reset();
                            anElement.addNodeId(nodeId[1]);
                            anElement.addNodeId(nodeId[3]);
                            anElement.addNodeId(nodeId[2]);
                            anElement.addNodeId(nodeId[7]);

                            anElement.setElementType(ET_TETRAHEDRON);

                            m_mesh.addElement(ne, anElement);
                            ne++;

                            // Tetrahedron 5
                            anElement.reset();
                            anElement.addNodeId(nodeId[1]);
                            anElement.addNodeId(nodeId[2]);
                            anElement.addNodeId(nodeId[6]);
                            anElement.addNodeId(nodeId[7]);

                            anElement.setElementType(ET_TETRAHEDRON);

                            m_mesh.addElement(ne, anElement);
                            ne++;

                            // Tetrahedron 6
                            anElement.reset();
                            anElement.addNodeId(nodeId[1]);
                            anElement.addNodeId(nodeId[6]);
                            anElement.addNodeId(nodeId[5]);
                            anElement.addNodeId(nodeId[7]);

                            anElement.setElementType(ET_TETRAHEDRON);

                            m_mesh.addElement(ne, anElement);
                            ne++;
                        }
                        else
                        {
                            // Hexahedron
                            CMshElement<Real> anElement;

                            anElement.addNodeId(nodeId[0]);
                            anElement.addNodeId(nodeId[1]);
                            anElement.addNodeId(nodeId[2]);
                            anElement.addNodeId(nodeId[3]);
                            anElement.addNodeId(nodeId[4]);
                            anElement.addNodeId(nodeId[5]);
                            anElement.addNodeId(nodeId[6]);
                            anElement.addNodeId(nodeId[7]);

                            anElement.setElementType(ET_HEXAHEDRON);

                            m_mesh.addElement(ne, anElement);
                            ne++;
                        }
                    }
                }
            }

            return true;
        }

        template <typename Real>
        CMshMesh<Real>& CMshBasicMesher<Real>::mesh()
        {
            return m_mesh;
        }
    }
}
