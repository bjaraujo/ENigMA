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

#include <Eigen/Dense>

namespace ENigMA
{

    namespace geometry
    {

        template <typename Real>
        struct CGeoCoordinateSystem : public Eigen::Matrix<Real, 3, 3>
        {
        public:

            typedef Eigen::Matrix<Real, 3, 3> Base;

            inline CGeoCoordinateSystem(Real value = (Real) 0) { this->setConstant(value); }

            template <typename Derived>
            inline CGeoCoordinateSystem(const Eigen::MatrixBase<Derived>& p) : Base(p) { }

            template <typename Derived>
            CGeoCoordinateSystem &operator=(const Eigen::MatrixBase<Derived>& p)
            {
                this->Base::operator=(p);
                return *this;
            }

        };

    }

}



