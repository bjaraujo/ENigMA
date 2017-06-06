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

    namespace geometry
    {

        template <typename Real>
        CGeoOctree<Real>::CGeoOctree()
        {

        }

        template <typename Real>
        CGeoOctree<Real>::~CGeoOctree()
        {

        }

        template <typename Real>
        void CGeoOctree<Real>::reset()
        {

            CGeoContainer<CGeoCoordinate<Real>, Real>::reset();

        }

        template <typename Real>
        void CGeoOctree<Real>::createNode(CGeoOctreeNode<Real>& anOctreeNode, const CGeoOctreeNode<Real>& anOctreeParent, Integer i)
        {

            try
            {

                double xmin = anOctreeParent.xmin;
                double ymin = anOctreeParent.ymin;
                double zmin = anOctreeParent.zmin;
                double xmax = anOctreeParent.xmax;
                double ymax = anOctreeParent.ymax;
                double zmax = anOctreeParent.zmax;
                double xmid = anOctreeParent.xmid;
                double ymid = anOctreeParent.ymid;
                double zmid = anOctreeParent.zmid;

                switch (i)
                {
                case OL_BACK_BOTTOM_LEFT:
                    anOctreeNode.xmin = xmin;
                    anOctreeNode.xmax = xmid;
                    anOctreeNode.ymin = ymin;
                    anOctreeNode.ymax = ymid;
                    anOctreeNode.zmin = zmin;
                    anOctreeNode.zmax = zmid;
                    break;
                case OL_BACK_TOP_LEFT:
                    anOctreeNode.xmin = xmin;
                    anOctreeNode.xmax = xmid;
                    anOctreeNode.ymin = ymin;
                    anOctreeNode.ymax = ymid;
                    anOctreeNode.zmin = zmid;
                    anOctreeNode.zmax = zmax;
                    break;
                case OL_FRONT_BOTTOM_LEFT:
                    anOctreeNode.xmin = xmin;
                    anOctreeNode.xmax = xmid;
                    anOctreeNode.ymin = ymid;
                    anOctreeNode.ymax = ymax;
                    anOctreeNode.zmin = zmin;
                    anOctreeNode.zmax = zmid;
                    break;
                case OL_FRONT_TOP_LEFT:
                    anOctreeNode.xmin = xmin;
                    anOctreeNode.xmax = xmid;
                    anOctreeNode.ymin = ymid;
                    anOctreeNode.ymax = ymax;
                    anOctreeNode.zmin = zmid;
                    anOctreeNode.zmax = zmax;
                    break;
                case OL_BACK_BOTTOM_RIGHT:
                    anOctreeNode.xmin = xmid;
                    anOctreeNode.xmax = xmax;
                    anOctreeNode.ymin = ymin;
                    anOctreeNode.ymax = ymid;
                    anOctreeNode.zmin = zmin;
                    anOctreeNode.zmax = zmid;
                    break;
                case OL_BACK_TOP_RIGHT:
                    anOctreeNode.xmin = xmid;
                    anOctreeNode.xmax = xmax;
                    anOctreeNode.ymin = ymin;
                    anOctreeNode.ymax = ymid;
                    anOctreeNode.zmin = zmid;
                    anOctreeNode.zmax = zmax;
                    break;
                case OL_FRONT_BOTTOM_RIGHT:
                    anOctreeNode.xmin = xmid;
                    anOctreeNode.xmax = xmax;
                    anOctreeNode.ymin = ymid;
                    anOctreeNode.ymax = ymax;
                    anOctreeNode.zmin = zmin;
                    anOctreeNode.zmax = zmid;
                    break;
                case OL_FRONT_TOP_RIGHT:
                    anOctreeNode.xmin = xmid;
                    anOctreeNode.xmax = xmax;
                    anOctreeNode.ymin = ymid;
                    anOctreeNode.ymax = ymax;
                    anOctreeNode.zmin = zmid;
                    anOctreeNode.zmax = zmax;
                    break;
                }

                anOctreeNode.xmid = (anOctreeNode.xmin + anOctreeNode.xmax) * 0.5;
                anOctreeNode.ymid = (anOctreeNode.ymin + anOctreeNode.ymax) * 0.5;
                anOctreeNode.zmid = (anOctreeNode.zmin + anOctreeNode.zmax) * 0.5;

            }
            catch (const std::exception& e)
            {
                std::cout << "Error: std exception: " << e.what() << " in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
            }
            catch (...)
            {
                std::cout << "Error: unknown exception in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
            }

        }

        template <typename Real>
        EOctreeLocation CGeoOctree<Real>::dispatch(CGeoOctreeNode<Real>& anOctreeNode, Real x, Real y, Real z)
        {

            Real xmid = anOctreeNode.xmid;
            Real ymid = anOctreeNode.ymid;
            Real zmid = anOctreeNode.zmid;

            if(x < xmid)
            {
                if(y < ymid)
                {
                    if(z < zmid) 
                        return OL_BACK_BOTTOM_LEFT;
                    else       
                        return OL_BACK_TOP_LEFT;
                }
                else
                {
                    if (z < zmid) 
                        return OL_FRONT_BOTTOM_LEFT;
                    else        
                        return OL_BACK_TOP_LEFT;
                }
            }
            else
            {
                if(y < ymid)
                {
                    if(z < zmid) 
                        return OL_BACK_BOTTOM_RIGHT;
                    else       
                        return OL_BACK_TOP_RIGHT;
                }
                else
                {
                    if(z < zmid) 
                        return OL_FRONT_BOTTOM_RIGHT;
                    else       
                        return OL_FRONT_TOP_RIGHT;
                }
            }

        }

        template <typename Real>
        void CGeoOctree<Real>::createRecursive(CGeoOctreeNode<Real>& anOctreeNode)
        {

            try
            {

                int ok[8];

                Integer number_of_the_node = 0;

                anOctreeNode.nodes.resize(8);
                for (Integer i = 0; i < 8; ++i)
                {
                    createNode(anOctreeNode.nodes[i], anOctreeNode, i);
                }

                for (Integer k = 0; k < anOctreeNode.nbCoordinates(); ++k)
                {

                    for (Integer l = 0; l < 8; l++)
                        ok[l] = 0;

                    Integer i = anOctreeNode.coordinate(k);

                    const Real C_EPS = std::numeric_limits<Real>::epsilon();

                    double xm = CGeoContainer<CGeoCoordinate<Real>, Real>::m_geometricObjects[i].x() - C_EPS;
                    double xp = CGeoContainer<CGeoCoordinate<Real>, Real>::m_geometricObjects[i].x() + C_EPS;

                    double ym = CGeoContainer<CGeoCoordinate<Real>, Real>::m_geometricObjects[i].y() - C_EPS;
                    double yp = CGeoContainer<CGeoCoordinate<Real>, Real>::m_geometricObjects[i].y() + C_EPS;

                    double zm = CGeoContainer<CGeoCoordinate<Real>, Real>::m_geometricObjects[i].z() - C_EPS;
                    double zp = CGeoContainer<CGeoCoordinate<Real>, Real>::m_geometricObjects[i].z() + C_EPS;

                    number_of_the_node = static_cast<Integer> (dispatch(anOctreeNode, xm, ym, zm));
                    if (ok[number_of_the_node] == 0)
                        ok[number_of_the_node] = 1;

                    number_of_the_node = static_cast<Integer> (dispatch(anOctreeNode, xm, ym, zp));
                    if (ok[number_of_the_node] == 0)
                        ok[number_of_the_node] = 1;

                    number_of_the_node = static_cast<Integer> (dispatch(anOctreeNode, xm, yp, zm));
                    if (ok[number_of_the_node] == 0)
                        ok[number_of_the_node] = 1;

                    number_of_the_node = static_cast<Integer> (dispatch(anOctreeNode, xm, yp, zp));
                    if (ok[number_of_the_node] == 0)
                        ok[number_of_the_node] = 1;

                    number_of_the_node = static_cast<Integer> (dispatch(anOctreeNode, xp, ym, zm));
                    if (ok[number_of_the_node] == 0)
                        ok[number_of_the_node] = 1;

                    number_of_the_node = static_cast<Integer> (dispatch(anOctreeNode, xp, ym, zp));
                    if (ok[number_of_the_node] == 0)
                        ok[number_of_the_node] = 1;

                    number_of_the_node = static_cast<Integer> (dispatch(anOctreeNode, xp, yp, zm));
                    if (ok[number_of_the_node] == 0)
                        ok[number_of_the_node] = 1;

                    number_of_the_node = static_cast<Integer> (dispatch(anOctreeNode, xp, yp, zp));
                    if (ok[number_of_the_node] == 0)
                        ok[number_of_the_node] = 1;

                    for (Integer m = 0; m < 8; m++)
                    {
                        if (ok[m] == 1)
                        {
                            anOctreeNode.nodes[m].addCoordinate(i);
                        }
                    }
                }

                for (Integer i = 0; i < 8; ++i)
                {

                    if (anOctreeNode.nodes[i].nbCoordinates() != anOctreeNode.nbCoordinates())
                    {
                        if (anOctreeNode.nodes[i].nbCoordinates() > MAX_NB_ELEMENTS_PER_NODE)
                        {
                            // Then we need to construct another node starting from this one
                            createRecursive(anOctreeNode.nodes[i]);
                            anOctreeNode.setIsLeaf(false);
                            anOctreeNode.reset();
                        }
                    }

                }
            }
            catch (const std::exception& e)
            {
                std::cout << "Error: std exception: " << e.what() << " in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
            }
            catch (...)
            {
                std::cout << "Error: unknown exception in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
            }

        }

        template <typename Real>
        void CGeoOctree<Real>::build()
        {

            try
            {

                // Get bounding box
                for (Integer i = 0; i < static_cast<Integer> (CGeoContainer<CGeoCoordinate<Real>, Real>::m_geometricObjects.size()); ++i)
                    m_boundingBox.addCoordinate(CGeoContainer<CGeoCoordinate<Real>, Real>::m_geometricObjects[i]);

                const Real C_EPS = std::numeric_limits<Real>::epsilon();

                m_root.xmin = m_boundingBox.min().x() - C_EPS;
                m_root.xmax = m_boundingBox.min().x() + C_EPS;

                m_root.ymin = m_boundingBox.min().y() - C_EPS;
                m_root.ymax = m_boundingBox.min().y() + C_EPS;

                m_root.zmin = m_boundingBox.min().z() - C_EPS;
                m_root.zmax = m_boundingBox.min().z() + C_EPS;

                m_root.xmid = (m_boundingBox.min().x() + m_boundingBox.max().x()) * 0.5;
                m_root.ymid = (m_boundingBox.min().y() + m_boundingBox.max().y()) * 0.5;
                m_root.zmid = (m_boundingBox.min().z() + m_boundingBox.max().z()) * 0.5;

                for (Integer i = 0; i < static_cast<Integer> (CGeoContainer<CGeoCoordinate<Real>, Real>::m_geometricObjects.size()); ++i)
                {
                    m_root.addCoordinate(i);
                }

                createRecursive(m_root);

            }
            catch (const std::exception& e)
            {
                std::cout << "Error: std exception: " << e.what() << " in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
            }
            catch (...)
            {
                std::cout << "Error: unknown exception in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
            }

        }

        template <typename Real>
        CGeoOctreeNode<Real>& CGeoOctree<Real>::findRecursive(CGeoOctreeNode<Real>& anOctreeNode, CGeoCoordinate<Real>& aCoordinate)
        {

            Integer octreeLocation = static_cast<Integer> (dispatch(anOctreeNode, aCoordinate.x(), aCoordinate.y(), aCoordinate.z()));
    
            if (!anOctreeNode.isLeaf())
                findRecursive(anOctreeNode.nodes[octreeLocation], aCoordinate);
                
            return anOctreeNode.nodes[octreeLocation];

        }

        template <typename Real>
        void CGeoOctree<Real>::find(std::vector<Integer>& coordinateIds, CGeoCoordinate<Real>& aCoordinate, const Real aTolerance)
        {

            try
            {

                CGeoOctreeNode<Real> anOctreeNode = findRecursive(m_root, aCoordinate);

                for (Integer i = 0; i < anOctreeNode.nbCoordinates(); ++i)
                {

                    Integer aCoordinateIndex = anOctreeNode.coordinate(i);

                    Real dist = (CGeoContainer<CGeoCoordinate<Real>, Real>::m_geometricObjects[aCoordinateIndex] - aCoordinate).norm();

                    if (dist <= aTolerance)
                        coordinateIds.push_back(CGeoContainer<CGeoCoordinate<Real>, Real>::m_geometricObjectIds[aCoordinateIndex]);

                }

            }
            catch (const std::exception& e)
            {
                std::cout << "Error: std exception: " << e.what() << " in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
            }
            catch (...)
            {
                std::cout << "Error: unknown exception in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
            }

        }

    }

}

