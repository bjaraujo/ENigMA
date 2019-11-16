// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

namespace ENigMA {

namespace geometry {

    template <typename Real>
    inline Real CGeoVector<Real>::angle(const CGeoVector<Real>& vec)
    {
        Real rad = 0.0;
        
        if (this->squaredNorm() > 0 && vec.squaredNorm() > 0)
        {
            rad = this->dot(vec) / std::sqrt(this->squaredNorm() * vec.squaredNorm());
        }

        if (rad < -1.0)
            rad = -1.0;
        if (rad > +1.0)
            rad = +1.0;

        return (acos(rad));
    }

    template <typename Real>
    inline void CGeoVector<Real>::rotate(const Real angle)
    {
        Real r;
        Real theta;

        r = std::sqrt((this->x() * this->x()) + (this->y() * this->y()));
        theta = atan2(this->y(), this->x());
        this->x() = r * std::cos(theta + angle);
        this->y() = r * std::sin(theta + angle);
    }

    template <typename Real>
    std::ostream& operator<<(std::ostream& output, CGeoVector<Real>& aVector)
    {
        output << aVector.x() << ", " << aVector.y() << ", " << aVector.z();
        return output;
    }
}
}
