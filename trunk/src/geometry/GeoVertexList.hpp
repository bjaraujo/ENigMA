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

#include <vector>

#include "GeoCoordinate.hpp"

namespace ENigMA
{

    namespace geometry
    {

        template <typename Real>
        class CGeoVertexList
        {
        protected:

            std::vector<CGeoCoordinate<Real> > m_vertices;

        public:
            CGeoVertexList();
            ~CGeoVertexList();

            Integer nbVertices();

            virtual void reset();

            void addVertex(CGeoCoordinate<Real>& aVertex);
            void insertVertex(const Integer aVertexIndex, CGeoCoordinate<Real>& aVertex);
            void removeVertex(const Integer aVertexIndex);

            CGeoCoordinate<Real>& vertex(const Integer aVertexIndex);

            void removeDuplicates();
            void removeCollinear(const Real aTolerance = 0.0);

            void invert();

        };

    }

}

#include "GeoVertexList_Imp.hpp"
