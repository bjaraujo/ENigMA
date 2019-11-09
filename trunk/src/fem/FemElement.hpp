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

namespace fem {

    template <typename Real>
    class CFemElement {
    protected:
        Real m_dt; // delta t
        Real m_diffusionCoefficient; // ex: viscosity
        Real m_convectionCoefficient; // ex: velocity

        bool m_transient;

        virtual void setTransientTerm() = 0;
        virtual void setDiffusionTerm() = 0;
        virtual void setConvectiveTerm() = 0;

    public:
        CFemElement();
        virtual ~CFemElement();

        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> ddt;
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> laplacian;
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> divergence;
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> gradient;
        Eigen::Matrix<Real, Eigen::Dynamic, 1> source;

        void setDt(const Real aValue);
        Real dt() const;

        void setDiffusionCoefficient(const Real aValue);
        Real diffusionCoefficient() const;

        void setConvectionCoefficient(const Real aValue);
        Real convectionCoefficient() const;

        void setTransient(const bool aValue);
        bool transient() const;

        virtual void rebuild() = 0;
        virtual void update() = 0;
    };
}
}

#include "FemElement_Imp.hpp"
