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

    namespace fem
    {

        template <typename Real>
        class CFemElement
        {
        protected:

            Real m_dt;                            // delta t
            Real m_diffusionCoefficient;          // ex: viscosity
            Real m_convectionCoefficient;         // ex: velocity

            bool m_transient;

            virtual void setTransientTerm() = 0;
            virtual void setDiffusionTerm() = 0;
            virtual void setConvectiveTerm() = 0;

        public:

            CFemElement();
            ~CFemElement();

            Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> ddt;
            Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> laplacian;
            Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> divergence;
            Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> gradient;
            Eigen::Matrix<Real, Eigen::Dynamic, 1> source;

            Real& dt();
            Real& diffusionCoefficient();
            Real& convectionCoefficient();

            bool& transient();

            virtual void rebuild() = 0;
            virtual void update() = 0;

        };

    }

}

#include "FemElement_Imp.hpp"
