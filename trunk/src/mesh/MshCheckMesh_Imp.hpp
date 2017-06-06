// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// <Copyright> Copyright (c) 2012, All Rights Reserved </Copyright>
// *****************************************************************************

#pragma once

namespace ENigMA
{

    namespace mesh
    {

        template <typename Real>
        CMshCheckMesh<Real>::CMshCheckMesh()
        {

        }

        template <typename Real>
        CMshCheckMesh<Real>::~CMshCheckMesh()
        {

        }

        template <typename Real>
        bool CMshCheckMesh<Real>::checkOpen(CMshMesh<Real>& aMesh)
        {

            for (Integer i = 0; i < aMesh.nbElements(); ++i)
            {

                Integer anElementId = aMesh.elementId(i);

                CMshElement<Real>& anElement = aMesh.element(anElementId);

                for (Integer j = 0; j < anElement.nbFaceIds(); ++j)
                {

                    Integer aFaceId = anElement.faceId(j);
                    CMshFace<Real>& aFace = aMesh.face(aFaceId);

                    if (!aFace.hasPair())
                        return true;

                }

            }

            return false;

        }

        template <typename Real>
        bool CMshCheckMesh<Real>::checkIntersections(CMshMesh<Real>& aMesh, const Real aTolerance)
        {

            CGeoRtree<Real> aRtree;

            for (Integer i = 0; i < aMesh.nbElements(); ++i)
            {

                Integer anElementId = aMesh.elementId(i);

                aRtree.addGeometricObject(anElementId, aMesh.boundingBox(anElementId));

            }

            std::vector<Integer> sElements;

            for (Integer i = 0; i < aMesh.nbElements(); ++i)
            {

                Integer anElementId = aMesh.elementId(i);

                CMshElement<Real> anElement = aMesh.element(anElementId);

                Integer nLines = 0;
                Integer nTriangles = 0;

                CGeoLine<Real> aLine1;
                CGeoTriangle<Real> aTriangle1;

                if (anElement.elementType() == ET_BEAM)
                {

                    aLine1.setStartPoint(aMesh.node(anElement.nodeId(0)));
                    aLine1.setEndPoint(aMesh.node(anElement.nodeId(1)));
                    nLines++;

                } 
                else if (anElement.elementType() == ET_TRIANGLE)
                {

                    aTriangle1.addVertex(aMesh.node(anElement.nodeId(0)));
                    aTriangle1.addVertex(aMesh.node(anElement.nodeId(1)));
                    aTriangle1.addVertex(aMesh.node(anElement.nodeId(2)));
                    nTriangles++;

                }

                aRtree.find(sElements, aMesh.boundingBox(anElementId));

                for (Integer j = 0; j < sElements.size(); ++j)
                {

                    Integer anotherElementId = sElements[j];

                    if (anElementId == anotherElementId)
                        continue;

                    CMshElement<Real> anotherElement = aMesh.element(anotherElementId);

                    CGeoLine<Real> aLine2;
                    CGeoTriangle<Real> aTriangle2;

                    if (anotherElement.elementType() == ET_BEAM)
                    {

                        aLine2.setStartPoint(aMesh.node(anotherElement.nodeId(0)));
                        aLine2.setEndPoint(aMesh.node(anotherElement.nodeId(1)));
                        nLines++;

                    } 
                    else if (anotherElement.elementType() == ET_TRIANGLE)
                    {

                        aTriangle2.addVertex(aMesh.node(anotherElement.nodeId(0)));
                        aTriangle2.addVertex(aMesh.node(anotherElement.nodeId(1)));
                        aTriangle2.addVertex(aMesh.node(anotherElement.nodeId(2)));
                        nTriangles++;

                    }

                    CGeoIntersectionType anIntersectionType;

                    if (nLines == 2)
                    {
                        CGeoCoordinate<Real> aPoint;
                        if (aLine1.intersects(aLine2, aPoint, anIntersectionType, aTolerance))
                            if (anIntersectionType == IT_EDGE || anIntersectionType == IT_COINCIDENT || 
                                anIntersectionType == IT_INTERNAL || anIntersectionType == IT_SWAP)
                            return true;
                    }

                    if (nTriangles == 2)
                    {
                        if (aTriangle1.intersects(aTriangle2, anIntersectionType, aTolerance))
                            if (anIntersectionType == IT_EDGE || anIntersectionType == IT_COINCIDENT || 
                                anIntersectionType == IT_INTERNAL || anIntersectionType == IT_SWAP)
                            return true;
                    }

                }

            }

            return false;

        }

    }

}


