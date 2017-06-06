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

namespace ENigMA
{

    namespace geometry
    {

        template <typename Real>
        inline void CGeoCoordinate<Real>::transform(const CGeoCoordinateSystem<Real>& cs)
        {

            *this = (cs * *this);

        }

        template <typename Real>
        std::ostream &operator<< (std::ostream &output, CGeoCoordinate<Real> &aCoordinate)
        {

            output << aCoordinate.x() << ", " << aCoordinate.y() << ", " << aCoordinate.z();
            
            return output;
        }

    }

}
