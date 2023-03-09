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
    namespace geometry
    {
        template <typename Real>
        CGeoOctreeNode<Real>::CGeoOctreeNode()
            : m_isLeaf(true)
        {
        }

        template <typename Real>
        CGeoOctreeNode<Real>::~CGeoOctreeNode()
        {
        }

        template <typename Real>
        void CGeoOctreeNode<Real>::addCoordinate(Integer aCoordinateIndex)
        {
            m_coordinateList.emplace_back(aCoordinateIndex);
        }

        template <typename Real>
        Integer CGeoOctreeNode<Real>::nbCoordinates() const
        {
            return static_cast<Integer>(m_coordinateList.size());
        }

        template <typename Real>
        Integer CGeoOctreeNode<Real>::coordinate(Integer aCoordinateIndex)
        {
            return m_coordinateList[aCoordinateIndex];
        }

        template <typename Real>
        void CGeoOctreeNode<Real>::setIsLeaf(bool isLeaf)
        {
            m_isLeaf = isLeaf;
        }

        template <typename Real>
        bool CGeoOctreeNode<Real>::isLeaf()
        {
            return m_isLeaf;
        }

        template <typename Real>
        void CGeoOctreeNode<Real>::reset()
        {
            m_coordinateList.clear();
        }
    }
}
