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

namespace ENigMA {

namespace geometry {

    template <typename Real>
    struct CGeoVector : public Eigen::Matrix<Real, 3, 1> {
    public:
        typedef Eigen::Matrix<Real, 3, 1> Base;

        inline CGeoVector(Real value = (Real)0) { this->setConstant(value); }

        inline CGeoVector(const Real x, const Real y, const Real z)
            : Base(x, y, z)
        {
        }

        template <typename Derived>
        inline CGeoVector(const Eigen::MatrixBase<Derived>& p)
            : Base(p)
        {
        }

        template <typename Derived>
        CGeoVector& operator=(const Eigen::MatrixBase<Derived>& p)
        {
            this->Base::operator=(p);
            return *this;
        }

        inline Real angle(const CGeoVector<Real>& vec);
        inline void rotate(const Real angle);
    };

    template <typename Real>
    std::ostream& operator<<(std::ostream& output, CGeoVector<Real>& aVector);
}
}

#include "GeoVector_Imp.hpp"
