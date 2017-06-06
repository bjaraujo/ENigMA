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
#include "GeoBoundingBox.hpp"
#include "GeoContainer.hpp"

namespace ENigMA
{

    namespace geometry
    {

        template <typename Real>
        class CGeoHashGrid : public CGeoContainer<CGeoCoordinate<Real>, Real>
        {
        private:

            Integer m_nbCellsX, m_nbCellsY, m_nbCellsZ;
            Integer m_nbCellsXY, m_nbCells;

            std::vector<Integer> m_coordinateList;
            std::vector<Integer> m_coordinateListPtr;

            CGeoBoundingBox<Real> m_boundingBox;

            CGeoVector<Real> m_adOrig;
            CGeoVector<Real> m_adDelta;

        public:
            CGeoHashGrid();
            ~CGeoHashGrid();

            void reset();

            void build();

            void find(std::vector<Integer>& coordinateIds, CGeoCoordinate<Real>& aCoordinate, const Real aTolerance = 0.0);

        };

    }

}

#include "GeoHashGrid_Imp.hpp"
