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
#include "GeoTetrahedron.hpp"

using namespace ENigMA::geometry;

namespace ENigMA {
namespace fem {
    /*
        Tetrahedral element.

              eta

              |
              *2
             /| \  
              |   \
            / |     \
              |       \ 1
           /  *0-------*-----> xi
             /      - 
          / /    - 
           /  -
          *3
         zeta    

       */

    template <typename Real, Integer NbNodes, Integer Dof, Integer Order>
    class CFemTetrahedron : public CFemElement<Real>, public CFemSolid<Real>, public CGeoTetrahedron<Real> {
    };

    template <typename Real>
    class CFemTetrahedron<Real, 4, 1, 1> : public CFemElement<Real>, public CFemSolid<Real>, public CGeoTetrahedron<Real> {
    protected:
        Real m_x1, m_y1, m_z1;
        Real m_x2, m_y2, m_z2;
        Real m_x3, m_y3, m_z3;
        Real m_x4, m_y4, m_z4;

        Real m_x14, m_x24, m_x34;
        Real m_y14, m_y24, m_y34;
        Real m_z14, m_z24, m_z34;

        void setTransientTerm();
        void setDiffusionTerm();
        void setConvectiveTerm();

        void calculateB(Eigen::Matrix<Real, 3, 4>& B);

    public:
        CFemTetrahedron();
        virtual ~CFemTetrahedron();

        void rebuild();
        void update();

        void setSourceOnNode(const Integer aNodeIndex, const Real aValue);
        void setSourceOnEdge(const Integer anEdgeIndex, const Real aValue);
        void setSourceOnFace(const Integer anFaceIndex, const Real aValue);
    };

    template <typename Real>
    class CFemTetrahedron<Real, 4, 3, 1> : public CFemTetrahedron<Real, 4, 1, 1> {
    protected:
        void calculateB(Eigen::Matrix<Real, 6, 12>& B);
    };
}
}

#include "FemTetrahedron_Imp.hpp"
