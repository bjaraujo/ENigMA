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

#include "GeoTriangularPrism.hpp"
#include "IntTriangularPrism.hpp"
#include "FemElement.hpp"
#include "FemSolid.hpp"

using namespace ENigMA::geometry;
using namespace ENigMA::integration;

namespace ENigMA
{

    namespace fem
    {

        /*
                       5*
                 eta   /|\
                  |   / |  \
                  |  /  |    \
                  | /   |      \
                  |/    |        \
              +1.0-     |          \
                 /|     |            \
                / |    3*-------------*4
               /  |    /             /
              /   |   /             /
             /    |  /             /
           2*     | /             /
            |\    |/             /
            |  \  O-------------|----> xi
            |    \             /+1.0
            |   /  \          /
            |  /     \       /
            | /        \    /
            |/           \ /
           0*-------------*1
           /+1.0
          /
         zeta
           */

        template <typename Real, Integer NbNodes, Integer Dof, Integer Order>
        class CFemTriangularPrism : public CFemElement<Real>, public CFemSolid<Real>, public CIntTriangularPrism<Real>, public CGeoTriangularPrism<Real>
        {



        };

        template <typename Real>
        class CFemTriangularPrism<Real, 6, 1, 1> : public CFemElement<Real>, public CFemSolid<Real>, public CIntTriangularPrism<Real>, public CGeoTriangularPrism<Real>
        {
        protected:

            Eigen::Matrix<Real, 6, 6> m_jacobian;

            Real m_x1, m_y1, m_z1;
            Real m_x2, m_y2, m_z2;
            Real m_x3, m_y3, m_z3;
            Real m_x4, m_y4, m_z4;
            Real m_x5, m_y5, m_z5;
            Real m_x6, m_y6, m_z6;
    
            void setTransientTerm();
            void setDiffusionTerm();
            void setConvectiveTerm();

        public:

            CFemTriangularPrism();
            ~CFemTriangularPrism();

            void rebuild();
            void update();

            void setSourceOnNode(const Integer aNodeIndex, const Real aValue);
            void setSourceOnEdge(const Integer anEdgeIndex, const Real aValue);
            void setSourceOnFace(const Integer anFaceIndex, const Real aValue);

        };

    }

}

#include "FemTriangularPrism_Imp.hpp"


