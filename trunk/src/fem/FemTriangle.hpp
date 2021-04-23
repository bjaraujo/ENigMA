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
#include "GeoTriangle.hpp"
#include "IntTriangle.hpp"

using namespace ENigMA::geometry;
using namespace ENigMA::integration;

namespace ENigMA
{
    namespace fem
    {
        /*
        Triangle element.

        eta

        |
        *2
        | \
        |   \
        |     \
        |       \ 1
        *0-------*-----> xi
        */

        template <typename Real, Integer NbNodes, Integer Dof, Integer Order>
        class CFemTriangle : public CFemElement<Real>, public CFemFace<Real>, public CGeoTriangle<Real>
        {
        };

        template <typename Real>
        class CFemTriangle<Real, 3, 1, 1> : public CFemElement<Real>, public CFemFace<Real>, public CGeoTriangle<Real>
        {
        protected:
            Real m_x1, m_x2, m_x3;
            Real m_y1, m_y2, m_y3;

            Real m_x21, m_y21;
            Real m_x32, m_y32;
            Real m_x13, m_y13;

            void calculateB(Eigen::Matrix<Real, 2, 3>& B);

            void setTransientTerm();
            void setDiffusionTerm();
            void setConvectiveTerm();

        public:
            CFemTriangle();
            virtual ~CFemTriangle();

            void rebuild();
            void update();

            void setSourceOnNode(const Integer aNodeIndex, const Real aValue);
            void setSourceOnEdge(const Integer anEdgeIndex, const Real aValue);
        };

        template <typename Real>
        class CFemTriangle<Real, 3, 2, 1> : public CFemTriangle<Real, 3, 1, 1>
        {
        protected:
            void calculateB(Eigen::Matrix<Real, 3, 6>& B);
        };
    }
}

#include "FemTriangle_Imp.hpp"
