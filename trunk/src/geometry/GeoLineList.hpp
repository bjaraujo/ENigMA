// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include <vector>

#include "GeoLine.hpp"

namespace ENigMA {

namespace geometry {

    template <typename Real>
    class CGeoLineList : public CGeoLength<Real> {
    protected:
        std::vector<CGeoLine<Real>> m_lines;

    public:
        CGeoLineList();
        virtual ~CGeoLineList();

        Integer nbLines();

        void reset();

        void addLine(CGeoLine<Real>& aLine);
        CGeoLine<Real>& line(const Integer aLineIndex);

        inline void calculateLength(bool bReCalculate = false);
        inline void calculateBoundingBox(bool bReCalculate = false);

        void sort(const Real aTolerance = 0.0);

        void invert();
    };
}
}

#include "GeoLineList_Imp.hpp"
