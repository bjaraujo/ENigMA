// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "GeoCoordinate.hpp"

namespace ENigMA {

namespace mesh {

    template <typename Real>
    struct CMshNode : public ENigMA::geometry::CGeoCoordinate<Real> {
    public:
        typedef ENigMA::geometry::CGeoCoordinate<Real> Base;

        inline CMshNode(Real value = (Real)0) { this->setConstant(value); }

        inline CMshNode(const Real x, const Real y, const Real z)
            : Base(x, y, z)
        {
        }

        template <typename Derived>
        inline CMshNode(const Eigen::MatrixBase<Derived>& p)
            : Base(p)
        {
        }

        template <typename Derived>
        CMshNode& operator=(const Eigen::MatrixBase<Derived>& p)
        {
            this->Base::operator=(p);
            return *this;
        }
    };
}
}
