// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "FemElement.hpp"
#include "FemSolid.hpp"
#include "GeoHexahedron.hpp"
#include "IntHexahedron.hpp"

using namespace ENigMA::geometry;
using namespace ENigMA::integration;

namespace ENigMA {

namespace fem {

    /*
        Hexahedral element.

                        zeta
                        |
               4 *------|-----------* 5
                /|      |          /|
              /  |      |        /  |
          7 /    |      |     6/    |
           *------------|-----*     |
           |     |      |     |     |
           |     |      |----------------> xi
           |     |     /      |     |
           |   0 *----/-------|-----* 1
           |    /    /        |     /
           |  /    eta        |   /
           |/                 | /
           *------------------*
          3                   2

       */

    template <typename Real, Integer NbNodes, Integer Dof, Integer Order>
    class CFemHexahedron : public CFemElement<Real>, public CFemSolid<Real>, public CIntHexahedron<Real>, public CGeoHexahedron<Real> {
    };

    template <typename Real>
    class CFemHexahedron<Real, 8, 1, 1> : public CFemElement<Real>, public CFemSolid<Real>, public CIntHexahedron<Real>, public CGeoHexahedron<Real> {
    protected:
        Eigen::Matrix<Real, 8, 8> m_jacobian;

        Real m_x1, m_y1, m_z1;
        Real m_x2, m_y2, m_z2;
        Real m_x3, m_y3, m_z3;
        Real m_x4, m_y4, m_z4;
        Real m_x5, m_y5, m_z5;
        Real m_x6, m_y6, m_z6;
        Real m_x7, m_y7, m_z7;
        Real m_x8, m_y8, m_z8;

        void setTransientTerm();
        void setDiffusionTerm();
        void setConvectiveTerm();

    public:
        CFemHexahedron();
        ~CFemHexahedron();

        void rebuild();
        void update();

        void setSourceOnNode(const Integer aNodeIndex, const Real aValue);
        void setSourceOnEdge(const Integer anEdgeIndex, const Real aValue);
        void setSourceOnFace(const Integer anFaceIndex, const Real aValue);
    };
}
}

#include "FemHexahedron_Imp.hpp"
