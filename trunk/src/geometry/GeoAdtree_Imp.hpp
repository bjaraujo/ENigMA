// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#include <array>

#pragma once

namespace ENigMA {
namespace geometry {
    template <typename Real>
    CGeoAdtreeNode6<Real>::CGeoAdtreeNode6()
    {
        id = -1;

        leftNode = NULL;
        rightNode = NULL;
        fatherNode = NULL;

        nbChildren = 0;
    }

    template <typename Real>
    void CGeoAdtreeNode6<Real>::deleteChildren()
    {
        if (leftNode) {
            leftNode->deleteChildren();
            delete leftNode;
            leftNode = NULL;
        }

        if (rightNode) {
            rightNode->deleteChildren();
            delete rightNode;
            rightNode = NULL;
        }
    }

    template <typename Real>
    CGeoAdtree6<Real>::CGeoAdtree6()
    {
    }

    template <typename Real>
    void CGeoAdtree6<Real>::set(const Real acmin[6], const Real acmax[6])
    {
        memcpy(m_cmin, acmin, 6 * sizeof(Real));
        memcpy(m_cmax, acmax, 6 * sizeof(Real));

        m_root = new CGeoAdtreeNode6<Real>();
        m_root->sep = (m_cmin[0] + m_cmax[0]) / 2;
    }

    template <typename Real>
    void CGeoAdtree6<Real>::reset()
    {
        m_root->deleteChildren();
        delete m_root;
        m_root = new CGeoAdtreeNode6<Real>();
        m_root->sep = (m_cmin[0] + m_cmax[0]) / 2;
        m_nodes.clear();
    }

    template <typename Real>
    CGeoAdtree6<Real>::~CGeoAdtree6()
    {
        m_root->deleteChildren();
        delete m_root;
    }

    template <typename Real>
    void CGeoAdtree6<Real>::insert(const Real p[6], const Integer anId)
    {
        CGeoAdtreeNode6<Real>* node;
        CGeoAdtreeNode6<Real>* next;

        Integer dir;
        Integer lr(0);

        Real bmin[6];
        Real bmax[6];

        memcpy(bmin, m_cmin, 6 * sizeof(Real));
        memcpy(bmax, m_cmax, 6 * sizeof(Real));

        next = m_root;
        dir = 0;

        while (next) {
            node = next;

            if (node->id == -1) {
                memcpy(node->data, p, 6 * sizeof(Real));
                node->id = static_cast<Integer>(anId);

                m_nodes.push_back(node);

                return;
            }

            if (node->sep > p[dir]) {
                next = node->leftNode;
                bmax[dir] = node->sep;
                lr = 0;
            } else {
                next = node->rightNode;
                bmin[dir] = node->sep;
                lr = 1;
            }

            dir++;

            if (dir == 6)
                dir = 0;
        }

        next = new CGeoAdtreeNode6<Real>();

        memcpy(next->data, p, 6 * sizeof(Real));
        next->id = static_cast<Integer>(anId);
        next->sep = static_cast<Real>((bmin[dir] + bmax[dir]) / 2.0);

        m_nodes.push_back(next);

        if (lr)
            node->rightNode = next;
        else
            node->leftNode = next;

        next->fatherNode = node;

        while (node) {
            node->nbChildren++;
            node = node->fatherNode;
        }
    }

    template <typename Real>
    void CGeoAdtree6<Real>::remove(Integer anId)
    {
        CGeoAdtreeNode6<Real>* node = m_nodes[anId];

        node->id = -1;

        node = node->fatherNode;

        while (node) {
            node->nbChildren--;
            node = node->fatherNode;
        }
    }

    template <typename Real>
    void CGeoAdtree6<Real>::getIntersecting(const Real bmin[6], const Real bmax[6], std::vector<Integer>& sIds)
    {
        const int maxDepth = 1024;

        struct node_dir {
            CGeoAdtreeNode6<Real>* node;
            unsigned char dir;
        };

        std::array<node_dir, maxDepth> nodes;

        nodes[0].node = m_root;
        nodes[0].dir = 0;

        Integer stacks = 0;

        while (stacks >= 0) {
            node_dir a_node_dir = nodes.at(stacks);

            CGeoAdtreeNode6<Real>* node = a_node_dir.node;
            Integer dir = a_node_dir.dir;

            if (node->id != -1) {
                if (node->data[0] <= bmax[0] && node->data[1] <= bmax[1] && node->data[2] <= bmax[2] && node->data[3] >= bmin[3] && node->data[4] >= bmin[4] && node->data[5] >= bmin[5]) {
                    sIds.push_back(node->id);
                }
            }

            unsigned char ndir = (dir + 1) % 6;

            stacks--;

            if (node->leftNode && bmin[dir] <= node->sep) {
                stacks++;

                if (stacks >= maxDepth) {
                    std::cout << "Error: max depth of " << maxDepth << " reached!" << std::endl;
                    break;
                }

                nodes.at(stacks).node = node->leftNode;
                nodes.at(stacks).dir = ndir;
            }

            if (node->rightNode && bmax[dir] >= node->sep) {
                stacks++;

                if (stacks >= maxDepth) {
                    std::cout << "Error: max depth of " << maxDepth << " reached!" << std::endl;
                    break;
                }

                nodes.at(stacks).node = node->rightNode;
                nodes.at(stacks).dir = ndir;
            }
        }
    }

    template <typename Real>
    CGeoAdtree<Real>::CGeoAdtree()
    {
    }

    template <typename Real>
    CGeoAdtree<Real>::~CGeoAdtree()
    {
    }

    template <typename Real>
    void CGeoAdtree<Real>::set(CGeoBoundingBox<Real>& aBoundingBox)
    {
        m_boundingBox = aBoundingBox;

        Real acmin[6];
        Real acmax[6];

        acmin[0] = acmin[3] = m_boundingBox.min().x();
        acmin[1] = acmin[4] = m_boundingBox.min().y();
        acmin[2] = acmin[5] = m_boundingBox.min().z();

        acmax[0] = acmax[3] = m_boundingBox.max().x();
        acmax[1] = acmax[4] = m_boundingBox.max().y();
        acmax[2] = acmax[5] = m_boundingBox.max().z();

        m_tree.set(acmin, acmax);
    }
    template <typename Real>
    void CGeoAdtree<Real>::reset()
    {
        CGeoContainer<CGeoBoundingBox<Real>, Real>::reset();
        m_tree.reset();
    }

    template <typename Real>
    void CGeoAdtree<Real>::addGeometricObject(const Integer aBoundingBoxId, CGeoBoundingBox<Real>& aBoundingBox)
    {
        Real b[6];

        b[0] = aBoundingBox.min().x();
        b[1] = aBoundingBox.min().y();
        b[2] = aBoundingBox.min().z();

        b[3] = aBoundingBox.max().x();
        b[4] = aBoundingBox.max().y();
        b[5] = aBoundingBox.max().z();

        m_tree.insert(b, aBoundingBoxId);
    }

    template <typename Real>
    void CGeoAdtree<Real>::removeGeometricObject(const Integer aBoundingBoxId)
    {
        m_tree.remove(aBoundingBoxId);
    }

    template <typename Real>
    void CGeoAdtree<Real>::build()
    {
    }

    template <typename Real>
    void CGeoAdtree<Real>::find(std::vector<Integer>& boundingBoxIds, CGeoBoundingBox<Real>& aBoundingBox, const Real aTolerance)
    {
        Real b_min[6];
        Real b_max[6];

        b_min[0] = m_boundingBox.min().x();
        b_min[1] = m_boundingBox.min().y();
        b_min[2] = m_boundingBox.min().z();

        b_min[3] = aBoundingBox.min().x();
        b_min[4] = aBoundingBox.min().y();
        b_min[5] = aBoundingBox.min().z();

        b_max[0] = aBoundingBox.max().x();
        b_max[1] = aBoundingBox.max().y();
        b_max[2] = aBoundingBox.max().z();

        b_max[3] = m_boundingBox.max().x();
        b_max[4] = m_boundingBox.max().y();
        b_max[5] = m_boundingBox.max().z();

        boundingBoxIds.clear();

        m_tree.getIntersecting(b_min, b_max, boundingBoxIds);
    }
}
}
