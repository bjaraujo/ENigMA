// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

namespace ENigMA {
namespace geometry {
    template <typename Real>
    CGeoRtree<Real>::CGeoRtree()
    {
    }

    template <typename Real>
    CGeoRtree<Real>::~CGeoRtree()
    {
    }

    template <typename Real>
    void CGeoRtree<Real>::reset()
    {
        m_tree.RemoveAll();
        CGeoContainer<CGeoBoundingBox<Real>, Real>::reset();
    }

    template <typename Real>
    void CGeoRtree<Real>::addGeometricObject(const Integer aBoundingBoxId, CGeoBoundingBox<Real>& aBoundingBox)
    {
        Real b_min[3];
        Real b_max[3];

        b_min[0] = aBoundingBox.min().x();
        b_min[1] = aBoundingBox.min().y();
        b_min[2] = aBoundingBox.min().z();

        b_max[0] = aBoundingBox.max().x();
        b_max[1] = aBoundingBox.max().y();
        b_max[2] = aBoundingBox.max().z();

        m_tree.Insert(b_min, b_max, aBoundingBoxId);
    }

    template <typename Real>
    void CGeoRtree<Real>::removeGeometricObject(const Integer aBoundingBoxId, CGeoBoundingBox<Real>& aBoundingBox)
    {
        Real b_min[3];
        Real b_max[3];

        b_min[0] = aBoundingBox.min().x();
        b_min[1] = aBoundingBox.min().y();
        b_min[2] = aBoundingBox.min().z();

        b_max[0] = aBoundingBox.max().x();
        b_max[1] = aBoundingBox.max().y();
        b_max[2] = aBoundingBox.max().z();

        m_tree.Remove(b_min, b_max, aBoundingBoxId);
    }

    template <typename Real>
    void CGeoRtree<Real>::build()
    {
    }

    template <typename Real>
    void CGeoRtree<Real>::find(std::vector<Integer>& boundingBoxIds, CGeoBoundingBox<Real>& aBoundingBox, const Real aTolerance)
    {
        Real b_min[3];
        Real b_max[3];

        b_min[0] = aBoundingBox.min().x();
        b_min[1] = aBoundingBox.min().y();
        b_min[2] = aBoundingBox.min().z();

        b_max[0] = aBoundingBox.max().x();
        b_max[1] = aBoundingBox.max().y();
        b_max[2] = aBoundingBox.max().z();

        boundingBoxIds.clear();

        m_tree.Search(b_min, b_max, boundingBoxIds, NULL);
    }
}
}
