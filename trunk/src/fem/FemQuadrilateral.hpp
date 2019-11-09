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
#include "FemFace.hpp"
#include "GeoQuadrilateral.hpp"
#include "IntQuadrilateral.hpp"

using namespace ENigMA::geometry;
using namespace ENigMA::integration;

namespace ENigMA {

namespace fem {

    /*
        Quadrilateral element.

                 eta
                 |
        *3-------|--------*2
        |        |        |
        |        |        |
        |        |        |
        |        |-------------> xi
        |                 |
        |                 |
        |                 |
        *0----------------*1
        */

    template <typename Real, Integer NbNodes, Integer Dof, Integer Order>
    class CFemQuadrilateral : public CFemElement<Real>, public CFemFace<Real>, public CIntQuadrilateral<Real>, public CGeoQuadrilateral<Real> {
    };

    template <typename Real>
    class CFemQuadrilateral<Real, 4, 1, 1> : public CFemElement<Real>, public CFemFace<Real>, public CIntQuadrilateral<Real>, public CGeoQuadrilateral<Real> {
    protected:
        Eigen::Matrix<Real, 4, 4> m_jacobian;

        Real m_x1, m_x2, m_x3, m_x4;
        Real m_y1, m_y2, m_y3, m_y4;

        Real m_x21, m_x32, m_x43, m_x14;
        Real m_y21, m_y32, m_y43, m_y14;

        void setTransientTerm();
        void setDiffusionTerm();
        void setConvectiveTerm();

    public:
        CFemQuadrilateral();
        virtual ~CFemQuadrilateral();

        void rebuild();
        void update();

        void setSourceOnNode(const Integer aNodeIndex, const Real aValue);
        void setSourceOnEdge(const Integer anEdgeIndex, const Real aValue);
    };
}
}

#include "FemQuadrilateral_Imp.hpp"
