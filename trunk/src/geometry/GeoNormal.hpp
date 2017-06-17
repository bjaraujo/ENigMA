// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "GeoVector.hpp"

namespace ENigMA
{

    namespace geometry
    {

        template <typename Real>
        struct CGeoNormal : public CGeoVector<Real>
        {
        public:

            typedef CGeoVector<Real> Base;

            inline CGeoNormal(Real value = (Real) 0) { this->setConstant(value); }

            inline CGeoNormal(const Real x, const Real y, const Real z) : Base(x, y, z) { Base::normalize(); }

            template <typename Derived>
            inline CGeoNormal(const Eigen::MatrixBase<Derived>& p) : Base(p) { }

            template <typename Derived>
            CGeoNormal &operator=(const Eigen::MatrixBase<Derived>& p)
            {
                this->Base::operator=(p);
                return *this;
            }

        };

    }

}
