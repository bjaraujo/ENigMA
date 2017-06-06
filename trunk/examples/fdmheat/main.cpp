// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// <Copyright> Copyright (c) 2012, All Rights Reserved </Copyright>
// *****************************************************************************

#include <iostream>
#include <iomanip>

#include "MshBasicMesher.hpp"
#include "SleSystem.hpp"
#include "PdeEquation.hpp"
#include "PosGnuplot.hpp"
#include "PosGmsh.hpp"
#include "AnaTemperature.hpp"

using namespace ENigMA::mesh;
using namespace ENigMA::post;
using namespace ENigMA::analytical;

void steadyHeatConduction1D()
{

    /***************************************************************
     * One-dimensional steady-state heat conduction                *
     ***************************************************************/

    std::cout << "**** One-dimensional steady-state heat conduction ****" << std::endl;

    CGeoCoordinate<double> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<double> aPoint2(1.0, 0.0, 0.0);

    CGeoLine<double> aLine(aPoint1, aPoint2);

    CMshBasicMesher<double> aBasicMesher;

    const Integer nu = 100;

    aBasicMesher.generate(aLine, nu);

    CMshMesh<double> aMesh;
    
    aMesh = aBasicMesher.mesh();

    CPdeField<double> T;

    T.setMesh(aMesh);
    T.setDiscretMethod(DM_FDM);
    T.setDiscretOrder(DO_LINEAR);
    T.setDiscretLocation(DL_NODE);
    T.setSimulationType(ST_THERMAL);
    T.setNbDofs(1);

    // Set BC
    for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
    {

        Integer aNodeId = T.mesh().nodeId(i);
        CMshNode<double> aNode = T.mesh().node(aNodeId);

        if (fabs(aNode.x() - 0.0) < 1E-6)
            T.setFixedValue(i,  0.0);

        if (fabs(aNode.x() - 1.0) < 1E-6)
            T.setFixedValue(i,  1.0);

    }

    // Steady-state conduction in a line
    CPdeEquation<double> aPdeEquation(laplacian<double>(T) = 0);

    aPdeEquation.solve(T);

    CPosGnuplot<double> aPosGnuplot;
    aPosGnuplot.save(T, "fdm_heat_cond_st_1D.dat");

    CAnaTemperature<double> aAnaTemperature;

    for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
    {

        Integer aNodeId = T.mesh().nodeId(i);
        CMshNode<double> aNode = T.mesh().node(aNodeId);

        double x = aNode.x();

        double Ti;

        aAnaTemperature.steadyStateHeatConduction1D(x, Ti);

        T.u(i) = Ti;

    }

    aPosGnuplot.save(T, "ana_heat_cond_st_1D.dat");

    std::ofstream plotFile;

    plotFile.open("heat_cond_st_1D.plt", std::ios_base::out | std::ios_base::trunc);

    if (plotFile.is_open())
    {

        plotFile << "set yrange [0:1]" << std::endl;
        plotFile << "plot \"fdm_heat_cond_st_1D.dat\" using 2:5 title \"FDM\" with points, \\" << std::endl;
        plotFile << "      \"ana_heat_cond_st_1D.dat\" using 2:5 title \"Analytical\" with lines" << std::endl;
        plotFile << "pause -1" << std::endl;

        plotFile.close();

    }

}

void unsteadyHeatConduction1D()
{

    /***************************************************************
     * One-dimensional unsteady heat conduction                    *
     ***************************************************************/

    std::cout << "**** One-dimensional unsteady heat conduction ****" << std::endl;

    CGeoCoordinate<double> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<double> aPoint2(1.0, 0.0, 0.0);

    CGeoLine<double> aLine(aPoint1, aPoint2);

    CMshBasicMesher<double> aBasicMesher;

    const Integer nu = 100;

    aBasicMesher.generate(aLine, nu);

    CMshMesh<double> aMesh = aBasicMesher.mesh();

    CPdeField<double> T;

    T.setMesh(aMesh);
    T.setDiscretMethod(DM_FDM);
    T.setDiscretOrder(DO_LINEAR);
    T.setDiscretLocation(DL_NODE);
    T.setSimulationType(ST_THERMAL);
    T.setNbDofs(1);

    // Set BC and initial conditions
    T.u.resize(T.mesh().nbNodes());

    Integer ii;

    for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
    {

        Integer aNodeId = T.mesh().nodeId(i);
        CMshNode<double> aNode = T.mesh().node(aNodeId);

        if (fabs(aNode.x() - 1.0) < 1E-6)
        {
            T.setFixedValue(i, 0.0);
            ii = i;
        }

        T.u(i) = 1.0;

    }

    double rho = 1.0;    // density
    double Cp = 1.0;    // specific heat
    double k = 1.0;        // conductivity

    double dt = 0.01;
    Integer nIter = 10;
    double time = dt * nIter;

    // Unsteady conduction in a line
    for (Integer i = 0; i < nIter; ++i)
    {

        std::cout << "Time = " << dt * (i + 1) << std::endl;

        CPdeEquation<double> aPdeEquation(rho * Cp / dt * ddt<double>(T) - k * laplacian<double>(T) = 0);

        aPdeEquation.solve(T);

    }

    CPosGnuplot<double> aPosGnuplot;
    aPosGnuplot.save(T, "fdm_heat_cond_un_1D.dat");

    CAnaTemperature<double> aAnaTemperature;

    double Ti;

    for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
    {

        Integer aNodeId = T.mesh().nodeId(i);
        CMshNode<double> aNode = T.mesh().node(aNodeId);

        double x = aNode.x();

        aAnaTemperature.transientHeatConduction1D(x, time, 1.0, Ti);

        T.u(i) = Ti;

    }

    aPosGnuplot.save(T, "ana_heat_cond_un_1D.dat");

    std::ofstream plotFile;

    plotFile.open("heat_cond_un_1D.plt", std::ios_base::out | std::ios_base::trunc);

    if (plotFile.is_open())
    {

        plotFile << "set yrange [0:1]" << std::endl;
        plotFile << "plot \"fdm_heat_cond_un_1D.dat\" using 2:5 title \"FDM\" with linespoints, \\" << std::endl;
        plotFile << "      \"ana_heat_cond_un_1D.dat\" using 2:5 title \"Analytical\" with lines" << std::endl;
        plotFile << "pause -1" << std::endl;

        plotFile.close();

    }

}

void unsteadyHeatConvection1D()
{

    /***************************************************************
     * One-dimensional unsteady heat convection                    *
     ***************************************************************/

    std::cout << "**** One-dimensional unsteady heat convection ****" << std::endl;

    CGeoCoordinate<double> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<double> aPoint2(1.0, 0.0, 0.0);

    CGeoLine<double> aLine(aPoint1, aPoint2);

    CMshBasicMesher<double> aBasicMesher;

    const Integer nu = 100;

    aBasicMesher.generate(aLine, nu);

    CMshMesh<double> aMesh = aBasicMesher.mesh();

    CPdeField<double> T;

    T.setMesh(aMesh);
    T.setDiscretMethod(DM_FDM);
    T.setDiscretOrder(DO_LINEAR);
    T.setDiscretLocation(DL_NODE);
    T.setSimulationType(ST_THERMAL);
    T.setNbDofs(1);

    // Set BC and initial conditions
    T.u.resize(T.mesh().nbNodes());

    for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
    {

        Integer aNodeId = T.mesh().nodeId(i);
        CMshNode<double> aNode = T.mesh().node(aNodeId);

        double x = aNode.x();
        double p = (x-0.5)/0.08;

        T.u(i) = exp(-0.5*p*p);

    }

    double rho = 1.0;    // density
    double Cp = 1.0;    // specific heat

    double dt = 0.01;
    Integer nIter = 10;
    double time = dt * nIter;

    double u = 1.0;

    for (Integer i = 0; i < nIter; ++i)
    {

        std::cout << "Time = " << dt * (i + 1) << std::endl;

        CPdeEquation<double> aPdeEquation(rho * Cp / dt * ddt<double>(T) + u * divergence<double>(T) = 0);

        aPdeEquation.solve(T);

    }

    CPosGnuplot<double> aPosGnuplot;
    aPosGnuplot.save(T, "fdm_heat_conv_un_1D.dat");

    for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
    {

        Integer aNodeId = T.mesh().nodeId(i);
        CMshNode<double> aNode = T.mesh().node(aNodeId);

        double x = aNode.x() - u * time;
        double p = (x-0.5)/0.08;

        T.u(i) = exp(-0.5*p*p);

    }

    aPosGnuplot.save(T, "ini_heat_conv_un_1D.dat");

    std::ofstream plotFile;

    plotFile.open("heat_conv_un_1D.plt", std::ios_base::out | std::ios_base::trunc);

    if (plotFile.is_open())
    {

        plotFile << "plot \"fdm_heat_conv_un_1D.dat\" using 2:5 title \"FDM\" with linespoints, \\" << std::endl;
        plotFile << "      \"ini_heat_conv_un_1D.dat\" using 2:5 title \"Analytical\" with lines" << std::endl;
        plotFile << "pause -1" << std::endl;

        plotFile.close();

    }

}

int main(int argc, char** argv)
{

    steadyHeatConduction1D();    
    unsteadyHeatConduction1D();
    unsteadyHeatConvection1D();

}
