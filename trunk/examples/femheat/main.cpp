// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#include <iostream>
#include <iomanip>

#include "FemTriangle.hpp"
#include "MshBasicMesher.hpp"
#include "SleSystem.hpp"
#include "PdeEquation.hpp"
#include "PosGnuplot.hpp"
#include "PosGmsh.hpp"
#include "AnaTemperature.hpp"

using namespace ENigMA::fem;
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

    const Integer nu = 30;

    aBasicMesher.generate(aLine, nu);

    CMshMesh<double> aMesh = aBasicMesher.mesh();

    CPdeField<double> T;

    T.setMesh(aMesh);
    T.setDiscretMethod(DM_FEM);
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

    aPdeEquation.setElimination(T);

    aPdeEquation.solve(T);

    CPosGnuplot<double> aPosGnuplot;
    aPosGnuplot.save(T, "fem_heat_cond_st_1D.dat");

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
        plotFile << "plot \"fem_heat_cond_st_1D.dat\" using 2:5 title \"FEM\" with points, \\" << std::endl;
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

    const Integer nu = 30;

    aBasicMesher.generate(aLine, nu);

    CMshMesh<double> aMesh = aBasicMesher.mesh();

    CPdeField<double> T;

    T.setMesh(aMesh);
    T.setDiscretMethod(DM_FEM);
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

        if (fabs(aNode.x() - 1.0) < 1E-6)
            T.setFixedValue(i,  0.0);

        T.u(i) = 1.0;

    }

    double rho = 1.0;    // density
    double Cp = 1.0;    // specific heat
    double k = 1.0;        // conductivity

    double time = 0.1;
    Integer nIter = 10;
    double dt = time / nIter;

    // Unsteady conduction in a line
    for (Integer i = 0; i < nIter; ++i)
    {

        std::cout << "Time = " << dt * (i + 1) << std::endl;

        CPdeEquation<double> aPdeEquation(rho * Cp / dt * ddt<double>(T) - k * laplacian<double>(T) = 0);

        aPdeEquation.setElimination(T);
 
        aPdeEquation.solve(T);

    }

    CPosGnuplot<double> aPosGnuplot;
    aPosGnuplot.save(T, "fem_heat_cond_un_1D.dat");

    CAnaTemperature<double> aAnaTemperature;

    for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
    {

        Integer aNodeId = T.mesh().nodeId(i);
        CMshNode<double> aNode = T.mesh().node(aNodeId);

        double x = aNode.x();

        double Ti;

        aAnaTemperature.transientHeatConduction1D(x, time, 1.0, Ti);

        T.u(i) = Ti;

    }

    aPosGnuplot.save(T, "ana_heat_cond_un_1D.dat");

    std::ofstream plotFile;

    plotFile.open("heat_cond_un_1D.plt", std::ios_base::out | std::ios_base::trunc);

    if (plotFile.is_open())
    {

        plotFile << "set yrange [0:1]" << std::endl;
        plotFile << "plot \"fem_heat_cond_un_1D.dat\" using 2:5 title \"FEM\" with points, \\" << std::endl;
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

    const Integer nu = 50;

    aBasicMesher.generate(aLine, nu);

    CMshMesh<double> aMesh = aBasicMesher.mesh();

    CPdeField<double> T;

    T.setMesh(aMesh);
    T.setDiscretMethod(DM_FEM);
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
    aPosGnuplot.save(T, "fem_heat_conv_un_1D.dat");

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

        plotFile << "plot \"fem_heat_conv_un_1D.dat\" using 2:5 title \"FEM\" with linespoints, \\" << std::endl;
        plotFile << "      \"ini_heat_conv_un_1D.dat\" using 2:5 title \"Analytical\" with lines" << std::endl;
        plotFile << "pause -1" << std::endl;

        plotFile.close();

    }

}

void steadyHeatConduction2D()
{

    /***************************************************************
     * Two-dimensional steady-state heat conduction                *
     ***************************************************************/

    std::cout << "**** Two-dimensional steady-state heat conduction ****" << std::endl;
    
    CGeoCoordinate<double> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<double> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<double> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<double> aVertex4(0.0, 1.0, 0.0);

    CGeoQuadrilateral<double> aQuadrilateral;

    aQuadrilateral.addVertex(aVertex1);
    aQuadrilateral.addVertex(aVertex2);
    aQuadrilateral.addVertex(aVertex3);
    aQuadrilateral.addVertex(aVertex4);

    CMshBasicMesher<double> aBasicMesher;

    const Integer nu = 15;
    const Integer nv = 15;

    aBasicMesher.generate(aQuadrilateral, nu, nv, true);

    CMshMesh<double> aMesh = aBasicMesher.mesh();

    CPdeField<double> T;

    T.setMesh(aMesh);
    T.setDiscretMethod(DM_FEM);
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
            T.setFixedValue(i, 0.0);

        if (fabs(aNode.x() - 1.0) < 1E-6)
            T.setFixedValue(i, 0.0);

        if (fabs(aNode.y() - 0.0) < 1E-6)
            T.setFixedValue(i, 0.0);

        if (fabs(aNode.y() - 1.0) < 1E-6)
            T.setFixedValue(i, 1.0);

    }

    // Steady-state conduction in a plate
    CPdeEquation<double> aPdeEquation(laplacian<double>(T) = 0);

    aPdeEquation.setElimination(T);

    aPdeEquation.solve(T);

    CPosGmsh<double> aPosGmsh;
    aPosGmsh.save(T, "fem_heat_cond_st_2D.msh", "Temperature");

    CAnaTemperature<double> aAnaTemperature;

    double x, y, Ti;

    for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
    {

        Integer aNodeId = T.mesh().nodeId(i);
        CMshNode<double> aNode = T.mesh().node(aNodeId);

        x = aNode.x();
        y = aNode.y();

        aAnaTemperature.steadyStateHeatConduction2D(x, y, Ti);

        T.u(i) = Ti;

    }

    aPosGmsh.save(T, "ana_heat_cond_st_2D.msh", "Temperature");

}

void unsteadyHeatConduction2D()
{

    /***************************************************************
     * Two-dimensional unsteady heat conduction                    *
     ***************************************************************/

    std::cout << "**** Two-dimensional unsteady heat conduction ****" << std::endl;

    CGeoCoordinate<double> aVertex1(-0.5, -0.5, 0.0);
    CGeoCoordinate<double> aVertex2(+0.5, -0.5, 0.0);
    CGeoCoordinate<double> aVertex3(+0.5, +0.5, 0.0);
    CGeoCoordinate<double> aVertex4(-0.5, +0.5, 0.0);

    CGeoQuadrilateral<double> aQuadrilateral;

    aQuadrilateral.addVertex(aVertex1);
    aQuadrilateral.addVertex(aVertex2);
    aQuadrilateral.addVertex(aVertex3);
    aQuadrilateral.addVertex(aVertex4);

    CMshBasicMesher<double> aBasicMesher;

    const Integer nu = 15;
    const Integer nv = 15;

    aBasicMesher.generate(aQuadrilateral, nu, nv, false);

    CMshMesh<double> aMesh = aBasicMesher.mesh();

    CPdeField<double> T;

    T.setMesh(aMesh);
    T.setDiscretMethod(DM_FEM);
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

        if (fabs(aNode.x() + 0.5) < 1E-6)
            T.setFixedValue(i, 0.0);

        if (fabs(aNode.x() - 0.5) < 1E-6)
            T.setFixedValue(i, 0.0);

        if (fabs(aNode.y() + 0.5) < 1E-6)
            T.setFixedValue(i, 0.0);

        if (fabs(aNode.y() - 0.5) < 1E-6)
            T.setFixedValue(i, 0.0);

        T.u(i) = 1.0;

    }

    double rho = 1.0;    // density
    double Cp = 1.0;    // specific heat
    double k = 1.0;        // conductivity

    double time = 0.1;
    Integer nIter = 10;
    double dt = time / nIter;

    // Unsteady conduction in a plate
    for (Integer i = 0; i < nIter; ++i)
    {

        std::cout << "Time = " << dt * (i + 1) << std::endl;

        CPdeEquation<double> aPdeEquation(rho * Cp / dt * ddt<double>(T) - k * laplacian<double>(T) = 0);

        //aPdeEquation.setPenaltyFactor(T, 1000);
        aPdeEquation.setElimination(T);

        aPdeEquation.solve(T);

    }

    CPosGmsh<double> aPosGmsh;
    aPosGmsh.save(T, "fem_heat_un_2D.msh", "Temperature");

    CAnaTemperature<double> aAnaTemperature;

    double x, y, Ti;

    for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
    {

        Integer aNodeId = T.mesh().nodeId(i);
        CMshNode<double> aNode = T.mesh().node(aNodeId);

        x = aNode.x();
        y = aNode.y();

        aAnaTemperature.transientHeatConduction2D(x, y, time, 1.0, Ti);

        T.u(i) = Ti;

    }

    aPosGmsh.save(T, "ana_heat_un_2D.msh", "Temperature");

}

void steadyHeatConduction3D()
{

    /***************************************************************
     * Three-dimensional steady-state heat conduction              *
     ***************************************************************/

    std::cout << "**** Three-dimensional steady-state heat conduction ****" << std::endl;

    CGeoCoordinate<double> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<double> aVertex2(1.0, 0.0, 0.0);
    CGeoCoordinate<double> aVertex3(1.0, 1.0, 0.0);
    CGeoCoordinate<double> aVertex4(0.0, 1.0, 0.0);
    CGeoCoordinate<double> aVertex5(0.0, 0.0, 1.0);
    CGeoCoordinate<double> aVertex6(1.0, 0.0, 1.0);
    CGeoCoordinate<double> aVertex7(1.0, 1.0, 1.0);
    CGeoCoordinate<double> aVertex8(0.0, 1.0, 1.0);

    CGeoHexahedron<double> aHexahedron;

    aHexahedron.addVertex(aVertex1);
    aHexahedron.addVertex(aVertex2);
    aHexahedron.addVertex(aVertex3);
    aHexahedron.addVertex(aVertex4);
    aHexahedron.addVertex(aVertex5);
    aHexahedron.addVertex(aVertex6);
    aHexahedron.addVertex(aVertex7);
    aHexahedron.addVertex(aVertex8);

    CMshBasicMesher<double> aBasicMesher;

    const Integer nu = 5;
    const Integer nv = 5;
    const Integer nw = 5;

    aBasicMesher.generate(aHexahedron, nu, nv, nw, false);

    CMshMesh<double> aMesh = aBasicMesher.mesh();

    CPdeField<double> T;

    T.setMesh(aMesh);
    T.setDiscretMethod(DM_FEM);
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
            T.setFixedValue(i, 0.5);

        if (fabs(aNode.x() - 1.0) < 1E-6)
            T.setFixedValue(i, 1.0);

        if (fabs(aNode.y() - 0.0) < 1E-6)
            T.setFixedValue(i, 0.0);

        if (fabs(aNode.y() - 1.0) < 1E-6)
            T.setFixedValue(i, 0.0);

        if (fabs(aNode.z() - 0.0) < 1E-6)
            T.setFixedValue(i, 0.0);

        if (fabs(aNode.z() - 1.0) < 1E-6)
            T.setFixedValue(i, 0.0);

    }

    // Steady-state conduction in a cube
    CPdeEquation<double> aPdeEquation(laplacian<double>(T) = 0);

    //aPdeEquation.setPenaltyFactor(T, 1000);
    aPdeEquation.setElimination(T);

    aPdeEquation.solve(T);

    CPosGmsh<double> aPosGmsh;
    aPosGmsh.save(T, "fem_heat_cond_st_3D.msh", "Temperature");

    CAnaTemperature<double> aAnaTemperature;

    double x, y, z, Ti;

    for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
    {

        Integer aNodeId = T.mesh().nodeId(i);
        CMshNode<double> aNode = T.mesh().node(aNodeId);

        x = aNode.x();
        y = aNode.y();
        z = aNode.z();

        aAnaTemperature.steadyStateHeatConduction3D(x, y, z, Ti);

        T.u(i) = Ti;

    }

    aPosGmsh.save(T, "ana_heat_cond_st_3D.msh", "Temperature");

}

void unsteadyHeatConduction3D()
{

    /***************************************************************
     * Three-dimensional unsteady heat conduction                  *
     ***************************************************************/

    std::cout << "**** Three-dimensional unsteady heat conduction ****" << std::endl;

    CGeoCoordinate<double> aVertex1(-0.5, -0.5, -0.5);
    CGeoCoordinate<double> aVertex2(+0.5, -0.5, -0.5);
    CGeoCoordinate<double> aVertex3(+0.5, +0.5, -0.5);
    CGeoCoordinate<double> aVertex4(-0.5, +0.5, -0.5);
    CGeoCoordinate<double> aVertex5(-0.5, -0.5, +0.5);
    CGeoCoordinate<double> aVertex6(+0.5, -0.5, +0.5);
    CGeoCoordinate<double> aVertex7(+0.5, +0.5, +0.5);
    CGeoCoordinate<double> aVertex8(-0.5, +0.5, +0.5);

    CGeoHexahedron<double> aHexahedron;

    aHexahedron.addVertex(aVertex1);
    aHexahedron.addVertex(aVertex2);
    aHexahedron.addVertex(aVertex3);
    aHexahedron.addVertex(aVertex4);
    aHexahedron.addVertex(aVertex5);
    aHexahedron.addVertex(aVertex6);
    aHexahedron.addVertex(aVertex7);
    aHexahedron.addVertex(aVertex8);

    CMshBasicMesher<double> aBasicMesher;

    const Integer nu = 5;
    const Integer nv = 5;
    const Integer nw = 5;

    aBasicMesher.generate(aHexahedron, nu, nv, nw, false);

    CMshMesh<double> aMesh = aBasicMesher.mesh();

    CPdeField<double> T;

    T.setMesh(aMesh);
    T.setDiscretMethod(DM_FEM);
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

        if (fabs(aNode.x() + 0.5) < 1E-6)
            T.setFixedValue(i, 0.0);

        if (fabs(aNode.x() - 0.5) < 1E-6)
            T.setFixedValue(i, 0.0);

        if (fabs(aNode.y() + 0.5) < 1E-6)
            T.setFixedValue(i, 0.0);

        if (fabs(aNode.y() - 0.5) < 1E-6)
            T.setFixedValue(i,  0.0);

        if (fabs(aNode.z() + 0.5) < 1E-6)
            T.setFixedValue(i,  0.0);

        if (fabs(aNode.z() - 0.5) < 1E-6)
            T.setFixedValue(i,  0.0);

        T.u(i) = 1.0;

    }

    double rho = 1.0;    // density
    double Cp = 1.0;    // specific heat
    double k = 1.0;        // conductivity

    double time = 0.1;
    Integer nIter = 10;
    double dt = time / nIter;

    // Unsteady conduction in a cube
    for (Integer i = 0; i < nIter; ++i)
    {

        std::cout << "Time = " << dt * (i + 1) << std::endl;

        CPdeEquation<double> aPdeEquation(rho * Cp / dt * ddt<double>(T) - k * laplacian<double>(T) = 0);

        //aPdeEquation.setPenaltyFactor(T, 1000);
        aPdeEquation.setElimination(T);

        aPdeEquation.solve(T);

    }

    CPosGmsh<double> aPosGmsh;
    aPosGmsh.save(T, "fem_heat_un_3D.msh", "Temperature");

    CAnaTemperature<double> aAnaTemperature;

    double x, y, z, Ti;

    for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
    {

        Integer aNodeId = T.mesh().nodeId(i);
        CMshNode<double> aNode = T.mesh().node(aNodeId);

        x = aNode.x();
        y = aNode.y();
        z = aNode.z();

        aAnaTemperature.transientHeatConduction3D(x, y, z, time, 1.0, Ti);

        T.u(i) = Ti;

    }

    aPosGmsh.save(T, "ana_heat_un_3D.msh", "Temperature");

}

int main(int argc, char** argv)
{

    steadyHeatConduction1D();    
    unsteadyHeatConduction1D();

    unsteadyHeatConvection1D();

    steadyHeatConduction2D();
    unsteadyHeatConduction2D();

    steadyHeatConduction3D();
    unsteadyHeatConduction3D();

}
