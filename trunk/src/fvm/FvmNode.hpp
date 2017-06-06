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

#include "GeoCoordinate.hpp"

using namespace ENigMA::geometry;

namespace ENigMA
{

    namespace fvm
    {

        template <typename Real>
        class CFvmNode : public CGeoCoordinate<Real>
        {
        public:

            typedef CGeoCoordinate<Real> Base;

            inline CFvmNode(Real value = (Real) 0) { this->setConstant(value); }

            inline CFvmNode(const Real x, const Real y, const Real z) : Base(x, y, z) { }

            template <typename Derived>
            inline CFvmNode(const Eigen::MatrixBase<Derived>& p) : Base(p) { }

            template <typename Derived>
            CFvmNode &operator=(const Eigen::MatrixBase<Derived>& p)
            {
                this->Base::operator=(p);
                return *this;
            }

        };

    }

}

