// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include <Eigen/Dense>

#include "CmnTypes.hpp"
#include "GeoCoordinateSystem.hpp"

namespace ENigMA {
namespace geometry {
    template <typename Real>
    struct CGeoCoordinate : public Eigen::Matrix<Real, 3, 1> {
    public:
        typedef Eigen::Matrix<Real, 3, 1> Base;

        inline CGeoCoordinate(Real value = (Real)0) { this->setConstant(value); }

        inline CGeoCoordinate(const Real x, const Real y, const Real z)
            : Base(x, y, z)
        {
        }

        template <typename Derived>
        inline CGeoCoordinate(const Eigen::MatrixBase<Derived>& p)
            : Base(p)
        {
        }

        template <typename Derived>
        CGeoCoordinate& operator=(const Eigen::MatrixBase<Derived>& p)
        {
            this->Base::operator=(p);
            return *this;
        }

        inline void transform(const CGeoCoordinateSystem<Real>& cs);
    };

    template <typename Real>
    std::ostream& operator<<(std::ostream& output, CGeoCoordinate<Real>& aCoordinate);
}
}

#include "GeoCoordinate_Imp.hpp"
