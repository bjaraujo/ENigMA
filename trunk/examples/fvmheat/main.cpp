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

    CGeoCoordinate<double> aVertex1(+0.00, -0.05, -0.05);
    CGeoCoordinate<double> aVertex2(+1.00, -0.05, -0.05);
    CGeoCoordinate<double> aVertex3(+1.00, +0.05, -0.05);
    CGeoCoordinate<double> aVertex4(+0.00, +0.05, -0.05);
    CGeoCoordinate<double> aVertex5(+0.00, -0.05, +0.05);
    CGeoCoordinate<double> aVertex6(+1.00, -0.05, +0.05);
    CGeoCoordinate<double> aVertex7(+1.00, +0.05, +0.05);
    CGeoCoordinate<double> aVertex8(+0.00, +0.05, +0.05);

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

    const Integer nu = 100;
    const Integer nv = 1;
    const Integer nw = 1;

    aBasicMesher.generate(aHexahedron, nu, nv, nw);

    CMshMesh<double> aMesh = aBasicMesher.mesh();

    aMesh.generateFaces(1E-12);

    aMesh.calculateFaceCentroid();
    aMesh.calculateElementCentroid();

    CPdeField<double> T;

    T.setMesh(aMesh);
    T.setDiscretMethod(DM_FVM);
    T.setDiscretOrder(DO_LINEAR);
    T.setDiscretLocation(DL_ELEMENT_CENTER);
    T.setSimulationType(ST_THERMAL);
    T.setNbDofs(1);

    // Set BC and initial conditions
    T.u.resize(T.mesh().nbElements());

    for (Integer i = 0; i < T.mesh().nbElements(); ++i)
    {
        T.u(i) = 0.0;
    }

    for (Integer i = 0; i < T.mesh().nbFaces(); ++i)
    {

        Integer aFaceId = T.mesh().faceId(i);

        if (fabs(T.mesh().faceCentroid(aFaceId).x() - 0.0) < 1E-6)
        {
            CPdeBoundaryCondition<double> aFixedTemperature(BT_GENERIC_FIXED_VALUE);
            aFixedTemperature.addCondition(CT_GENERIC_FIXED_VALUE, 0.0);
            T.addBCFace(aFaceId, aFixedTemperature);
        }

        if (fabs(T.mesh().faceCentroid(aFaceId).x() - 1.0) < 1E-6)
        {
            CPdeBoundaryCondition<double> aFixedTemperature(BT_GENERIC_FIXED_VALUE);
            aFixedTemperature.addCondition(CT_GENERIC_FIXED_VALUE, 1.0);
            T.addBCFace(aFaceId, aFixedTemperature);
        }

    }

    double dt = 0.01;
    Integer nIter = 100;

    // Steady conduction in a line
    for (Integer i = 0; i < nIter; ++i)
    {

        // Steady-state conduction in a line
        CPdeEquation<double> aPdeEquation(1.0 / dt * ddt<double>(T) -laplacian<double>(T) = 0);

        aPdeEquation.solve(T);

    }

    CPosGnuplot<double> aPosGnuplot;
    aPosGnuplot.save(T, "fvm_heat_cond_st_1D.dat");

    CAnaTemperature<double> aAnaTemperature;

    for (Integer i = 0; i < T.mesh().nbElements(); ++i)
    {

        Integer anElementId = T.mesh().elementId(i);

        double x = T.mesh().elementCentroid(anElementId).x();

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
        plotFile << "plot \"fvm_heat_cond_st_1D.dat\" using 2:5 title \"FVM\" with points, \\" << std::endl;
        plotFile << "     \"ana_heat_cond_st_1D.dat\" using 2:5 title \"Analytical\" with lines" << std::endl;
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

    CGeoCoordinate<double> aVertex1(+0.00, -0.05, -0.05);
    CGeoCoordinate<double> aVertex2(+1.00, -0.05, -0.05);
    CGeoCoordinate<double> aVertex3(+1.00, +0.05, -0.05);
    CGeoCoordinate<double> aVertex4(+0.00, +0.05, -0.05);
    CGeoCoordinate<double> aVertex5(+0.00, -0.05, +0.05);
    CGeoCoordinate<double> aVertex6(+1.00, -0.05, +0.05);
    CGeoCoordinate<double> aVertex7(+1.00, +0.05, +0.05);
    CGeoCoordinate<double> aVertex8(+0.00, +0.05, +0.05);

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

    const Integer nu = 100;
    const Integer nv = 1;
    const Integer nw = 1;

    aBasicMesher.generate(aHexahedron, nu, nv, nw);

    CMshMesh<double> aMesh = aBasicMesher.mesh();

    aMesh.generateFaces(1E-12);

    aMesh.calculateFaceCentroid();
    aMesh.calculateElementCentroid();

    CPdeField<double> T;

    T.setMesh(aMesh);
    T.setDiscretMethod(DM_FVM);
    T.setDiscretOrder(DO_LINEAR);
    T.setDiscretLocation(DL_ELEMENT_CENTER);
    T.setSimulationType(ST_THERMAL);
    T.setNbDofs(1);

    // Set BC and initial conditions
    T.u.resize(T.mesh().nbElements());

    for (Integer i = 0; i < T.mesh().nbElements(); ++i)
    {
        T.u(i) = 1.0;
    }

    for (Integer i = 0; i < T.mesh().nbFaces(); ++i)
    {

        Integer aFaceId = T.mesh().faceId(i);

        if (fabs(T.mesh().faceCentroid(aFaceId).x() - 1.0) < 1E-6)
        {
            CPdeBoundaryCondition<double> aFixedTemperature(BT_GENERIC_FIXED_VALUE);
            aFixedTemperature.addCondition(CT_GENERIC_FIXED_VALUE, 0.0);
            T.addBCFace(aFaceId, aFixedTemperature);
        }

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

        CPdeEquation<double> aPdeEquation(rho * Cp / dt * ddt<double>(T) -k * laplacian<double>(T) = 0);

        aPdeEquation.solve(T);

    }

    CPosGnuplot<double> aPosGnuplot;
    aPosGnuplot.save(T, "fvm_heat_cond_un_1D.dat");

    CAnaTemperature<double> aAnaTemperature;

    for (Integer i = 0; i < T.mesh().nbElements(); ++i)
    {

        Integer anElementId = T.mesh().elementId(i);
        double x = T.mesh().elementCentroid(anElementId).x();

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
        plotFile << "plot \"fvm_heat_cond_un_1D.dat\" using 2:5 title \"FVM\" with linespoints, \\" << std::endl;
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

    CGeoCoordinate<double> aVertex1(+0.00, -0.05, -0.05);
    CGeoCoordinate<double> aVertex2(+1.00, -0.05, -0.05);
    CGeoCoordinate<double> aVertex3(+1.00, +0.05, -0.05);
    CGeoCoordinate<double> aVertex4(+0.00, +0.05, -0.05);
    CGeoCoordinate<double> aVertex5(+0.00, -0.05, +0.05);
    CGeoCoordinate<double> aVertex6(+1.00, -0.05, +0.05);
    CGeoCoordinate<double> aVertex7(+1.00, +0.05, +0.05);
    CGeoCoordinate<double> aVertex8(+0.00, +0.05, +0.05);

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

    const Integer nu = 100;
    const Integer nv = 1;
    const Integer nw = 1;

    aBasicMesher.generate(aHexahedron, nu, nv, nw);

    CMshMesh<double> aMesh = aBasicMesher.mesh();

    aMesh.generateFaces(1E-12);

    aMesh.calculateFaceCentroid();
    aMesh.calculateElementCentroid();

    CPdeField<double> T;

    T.setMesh(aMesh);
    T.setDiscretMethod(DM_FVM);
    T.setDiscretOrder(DO_LINEAR);
    T.setDiscretLocation(DL_ELEMENT_CENTER);
    T.setSimulationType(ST_THERMAL);
    T.setNbDofs(1);

    // Set BC and initial conditions
    T.u.resize(T.mesh().nbElements());

    for (Integer i = 0; i < T.mesh().nbElements(); ++i)
    {

        Integer anElementId = T.mesh().elementId(i);
        double x = T.mesh().elementCentroid(anElementId).x();

        double p = (x - 0.5) / 0.08;

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

        CPdeEquation<double> aPdeEquation(rho * Cp / dt * ddt<double>(T) +u * divergence<double>(T) = 0);

        aPdeEquation.solve(T);

    }

    CPosGnuplot<double> aPosGnuplot;
    aPosGnuplot.save(T, "fvm_heat_conv_un_1D.dat");

    for (Integer i = 0; i < T.mesh().nbElements(); ++i)
    {

        Integer anElementId = T.mesh().elementId(i);
        double x = T.mesh().elementCentroid(anElementId).x() - u * time;

        double p = (x - 0.5) / 0.08;

        T.u(i) = exp(-0.5*p*p);

    }

    aPosGnuplot.save(T, "ini_heat_conv_un_1D.dat");

    std::ofstream plotFile;

    plotFile.open("heat_conv_un_1D.plt", std::ios_base::out | std::ios_base::trunc);

    if (plotFile.is_open())
    {

        plotFile << "plot \"fvm_heat_conv_un_1D.dat\" using 2:5 title \"FVM\" with linespoints, \\" << std::endl;
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
