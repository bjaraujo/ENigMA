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
namespace analytical {
    template <typename Real>
    CAnaTemperature<Real>::CAnaTemperature()
    {
    }

    template <typename Real>
    CAnaTemperature<Real>::~CAnaTemperature()
    {
    }

    template <typename Real>
    void CAnaTemperature<Real>::steadyStateHeatConduction1D(Real x, Real& T)
    {
        // Source: Integration of heat equation

        // Case: steady-state heat conduction in a beam with assigned temperature T=0 on x=0 and T=1 on x=1.

        // x: [0, 1]

        // T: [0, 1]

        T = x;
    }

    template <typename Real>
    void CAnaTemperature<Real>::steadyStateHeatConduction2D(Real x, Real y, Real& T)
    {
        // Source: Myers, G. (1971). Analytical methods in conduction heat transfer. Mc-Graw Hill.

        // Case: steady-state heat conduction in a 1x1 square with assigned temperature T=1 on top side y=1. Other sides T=0.

        // x: [0, 1]
        // y: [0, 1]

        // T: [0, 1]

        const Real pi = std::acos(-1.0);

        T = 0.0;

        for (Integer p = 0; p < 25; p++) {
            T += 4.0 / pi * sin((2 * p + 1) * pi * x) * sinh((2 * p + 1) * pi * y) / ((2 * p + 1) * sinh((2 * p + 1) * pi));
        }
    }

    template <typename Real>
    void CAnaTemperature<Real>::steadyStateHeatConduction3D(Real x, Real y, Real z, Real& T)
    {
        // Source: Carslaw, H. and Jaeger, J. (1959). Conduction of heat in solids. Oxford Clarendon Press.

        // Case: steady-state heat conduction 1x1x1 cube with assigned temperature T=1 on side x=1 and T=0.5 on side x=0. Other sides T=0.

        // x: [0, 1]
        // y: [0, 1]
        // z: [0, 1]

        // T: [0, 1]

        const Real pi = std::acos(-1.0);

        T = 0.0;

        for (Integer p = 0; p < 10; p++) {
            for (Integer q = 0; q < 10; q++) {
                Real lambda = sqrt((2 * p + 1) * (2 * p + 1) * pi * pi + (2 * q + 1) * (2 * q + 1) * pi * pi);

                T += 16.0 / (pi * pi) * ((0.5 * sinh(lambda * (1.0 - x)) + sinh(lambda * x)) * sin((2 * p + 1) * pi * y) * sin((2 * q + 1) * pi * z)) / ((2 * p + 1) * (2 * q + 1) * sinh(lambda));
            }
        }
    }

    template <typename Real>
    void CAnaTemperature<Real>::steadyStateHeatConduction1D(Real x, Real Ta, Real Tb, Real h, Real Tinf, Real k, Real perimeter, Real sectionArea, Real length, Real& T)
    {
        // Source: http://www.thermoanalytics.com/support/validation/example003.html

        // Case: steady-state heat conduction in beam with end temperatures Ta, Tb and convection

        Real m = sqrt(h * perimeter / (k * sectionArea));

        Real c1 = -(exp(m * length) * Tb - exp(2 * m * length) * Ta + Tinf * exp(2 * m * length) - Tinf * exp(m * length)) / expm1(2 * m * length);
        Real c2 = (Ta - Tinf) - c1;

        T = Tinf + c1 * exp(-m * x) + c2 * exp(m * x);
    }

    template <typename Real>
    void CAnaTemperature<Real>::steadyStateHeatConvectionRadiation1D(Real x, Real Tb, Real h, Real e, Real k, Real perimeter, Real sectionArea, Real& T)
    {
        // Case: steady-state heat conduction in beam with temperature Tb at x=0 and radiation on sides

        // Stefan-boltzmann constant
        const Real sigma = 5.6704E-8;

        // Tinf = 0

        Real G, M, x1, Tb1;

        G = 2.0 / 3.0 * perimeter * e * sigma / (k * sectionArea);
        M = h * perimeter / (k * sectionArea);
        x1 = x;
        Tb1 = Tb;

        CAnaFunction<Real>::defineVariable("G", G);
        CAnaFunction<Real>::defineVariable("M", M);
        CAnaFunction<Real>::defineVariable("x1", x1);
        CAnaFunction<Real>::defineVariable("Tb1", Tb1);

        CAnaFunction<Real>::set("1/3*M^(-1/2)*(ln(((G*Tb^3+M)^(1/2)-M^(1/2))/((G*Tb^3+M)^(1/2)+M^(1/2)))-ln(((G*T^3+M)^(1/2)-M^(1/2))/((G*T^3+M)^(1/2)+M^(1/2))))-x");

        Integer nIter;

        T = CAnaFunction<Real>::root("T", 1E-12, Tb, nIter, 100, 1E-12);
    }

    template <typename Real>
    void CAnaTemperature<Real>::transientHeatConduction1D(Real x, Real t, Real alpha, Real& T)
    {
        // Case: transient heat conduction in beam length = 1.
        // Initial temperature Ti=1. Side x=0 is adiabatic, x=1 is set to T=0.

        const Real pi = std::acos(-1.0);

        T = 0.0;

        for (int p = 1; p < 100; p++) {
            Real lambda = (2 * p - 1) * pi * 0.5;

            T += 4.0 / pi * pow(-1.0, (p + 1)) / (2 * p - 1) * exp(-alpha * lambda * lambda * t) * cos(lambda * x);
        }
    }

    template <typename Real>
    void CAnaTemperature<Real>::transientHeatConduction2D(Real x, Real y, Real t, Real alpha, Real& T)
    {
        // Case: transient heat conduction in 1x1 square with center at 0,0.
        // All sides are at T=0. Initial temperature Ti=1.

        const Real pi = std::acos(-1.0);

        T = 0.0;

        for (int p = 0; p < 50; p++) {
            for (int q = 0; q < 50; q++) {
                Real lambda = alpha * pi * pi * ((2 * p + 1) * (2 * p + 1) + (2 * q + 1) * (2 * q + 1));

                Real phi = cos((2 * p + 1) * pi * x) * cos((2 * q + 1) * pi * y);

                T += 16.0 / (pi * pi) * pow(-1.0, (p + q)) / ((2 * p + 1) * (2 * q + 1)) * phi * exp(-lambda * t);
            }
        }
    }

    template <typename Real>
    void CAnaTemperature<Real>::transientHeatConduction3D(Real x, Real y, Real z, Real t, Real alpha, Real& T)
    {
        // Case: transient heat conduction in 1x1x1 cube with center at 0,0,0.
        // All sides at T=0. Initial temperature Ti=1.

        const Real pi = std::acos(-1.0);

        T = 0.0;

        for (int p = 0; p < 50; p++) {
            for (int q = 0; q < 50; q++) {
                for (int r = 0; r < 50; r++) {
                    Real lambda = alpha * pi * pi * ((2 * p + 1) * (2 * p + 1) + (2 * q + 1) * (2 * q + 1) + (2 * r + 1) * (2 * r + 1));

                    T += 64.0 / (pi * pi * pi) * pow(-1.0, p + q) / ((2 * p + 1) * (2 * q + 1) * (2 * r + 1)) * cos((2 * p + 1) * pi * x) * cos((2 * q + 1) * pi * y) * cos((2 * r + 1) * pi * z) * exp(-lambda * t);
                }
            }
        }
    }
}
}
