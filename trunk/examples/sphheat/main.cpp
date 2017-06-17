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
#include "SphParticles.hpp"

using namespace ENigMA::mesh;
using namespace ENigMA::post;
using namespace ENigMA::analytical;
using namespace ENigMA::sph;

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

    CMshMesh<double> aMesh = aBasicMesher.mesh();

    CPdeField<double> T;

    T.setMesh(aMesh);
    T.setDiscretMethod(DM_SPH);
    T.setDiscretOrder(DO_LINEAR);
    T.setDiscretLocation(DL_NODE);
    T.setSimulationType(ST_THERMAL);
    T.setNbDofs(1);

    // Set BC
    T.u.resize(T.mesh().nbNodes());

    for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
    {

        Integer aNodeId = T.mesh().nodeId(i);
        CMshNode<double> aNode = T.mesh().node(aNodeId);

        if (fabs(aNode.x() - 0.0) < 1E-6)
            T.setFixedValue(i, 0.0);

        if (fabs(aNode.x() - 1.0) < 1E-6)
            T.setFixedValue(i, 1.0);

        T.u(i) = 0.0;

    }

    CSphSpiky<double> aKernel(1);

    CSphParticles<double> sParticles(aKernel);

    double dt = 0.01;
    Integer nIter = 100;

    double h = 0.1;
    double diff = 1.0;
    double mass = 1.0;
    double rho = 1.0;

    sParticles.init(T, mass, rho, diff, h, dt);

    for (Integer i = 0; i < nIter; ++i)
    {

        // Solve
        sParticles.solve(T);

    }

    double sumError = 0.0;

    CAnaTemperature<double> aAnaTemperature;

    for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
    {

        Integer aNodeId = T.mesh().nodeId(i);
        CMshNode<double> aNode = T.mesh().node(aNodeId);
        double x = aNode.x();

        double Ti;

        aAnaTemperature.steadyStateHeatConduction1D(x, Ti);

        sumError += (T.u(i) - Ti) * (T.u(i) - Ti);

    }

    sumError /= T.mesh().nbNodes();

    std::cout << "Error 1: " << sumError << std::endl;

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

    const Integer nu = 25;
    const Integer nv = 25;

    aBasicMesher.generate(aQuadrilateral, nu, nv, false);

    CMshMesh<double> aMesh = aBasicMesher.mesh();

    CPdeField<double> T;

    T.setMesh(aMesh);
    T.setDiscretMethod(DM_SPH);
    T.setDiscretOrder(DO_LINEAR);
    T.setDiscretLocation(DL_NODE);
    T.setSimulationType(ST_THERMAL);
    T.setNbDofs(1);

    // Set BC
    T.u.resize(T.mesh().nbNodes());

    for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
    {

        Integer aNodeId = T.mesh().nodeId(i);
        CMshNode<double> aNode = T.mesh().node(aNodeId);

        if (fabs(aNode.x() - 0.0) < 1E-4)
            T.setFixedValue(i, 0.0);

        if (fabs(aNode.x() - 1.0) < 1E-4)
            T.setFixedValue(i, 0.0);

        if (fabs(aNode.y() - 0.0) < 1E-4)
            T.setFixedValue(i, 0.0);

        if (fabs(aNode.y() - 1.0) < 1E-4)
            T.setFixedValue(i, 1.0);

        T.u(i) = 0.0;

    }

    CSphSpiky<double> aKernel(2);

    CSphParticles<double> sParticles(aKernel);

    double dt = 0.01;
    Integer nIter = 10;

    double h = 1.0 / nu * 5.0;
    double diff = 1.0;
    double mass = 1.0;
    double rho = 1.0;

    sParticles.init(T, mass, rho, diff, h, dt);

    for (Integer i = 0; i < nIter; ++i)
    {

        // Solve
        sParticles.solve(T);

    }

    double sumError = 0.0;

    CAnaTemperature<double> aAnaTemperature;

    for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
    {

        Integer aNodeId = T.mesh().nodeId(i);
        CMshNode<double> aNode = T.mesh().node(aNodeId);
        double x = aNode.x();
        double y = aNode.y();

        double Ti;

        aAnaTemperature.steadyStateHeatConduction2D(x, y, Ti);

        sumError += (T.u(i) - Ti) * (T.u(i) - Ti);

    }

    sumError /= T.mesh().nbNodes();

    std::cout << "Error 2: " << sumError << std::endl;

}

void steadyHeatConduction3D()
{

    /***************************************************************
    * Three-dimensional steady-state heat conduction                *
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

    const Integer nu = 25;
    const Integer nv = 25;
    const Integer nw = 25;

    aBasicMesher.generate(aHexahedron, nu, nv, nw, false);

    CMshMesh<double> aMesh = aBasicMesher.mesh();

    CPdeField<double> T;

    T.setMesh(aMesh);
    T.setDiscretMethod(DM_SPH);
    T.setDiscretOrder(DO_LINEAR);
    T.setDiscretLocation(DL_NODE);
    T.setSimulationType(ST_THERMAL);
    T.setNbDofs(1);

    // Set BC
    T.u.resize(T.mesh().nbNodes());

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

        T.u(i) = 0.0;

    }

    CSphSpiky<double> aKernel(3);

    CSphParticles<double> sParticles(aKernel);

    double dt = 0.01;
    Integer nIter = 100;

    double h = 0.1;
    double diff = 1.0;
    double mass = 1.0;
    double rho = 1.0;

    sParticles.init(T, mass, rho, diff, h, dt);

    for (Integer i = 0; i < nIter; ++i)
    {

        // Solve
        sParticles.solve(T);

    }

    double sumError = 0.0;

    CAnaTemperature<double> aAnaTemperature;

    for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
    {

        Integer aNodeId = T.mesh().nodeId(i);
        CMshNode<double> aNode = T.mesh().node(aNodeId);
        double x = aNode.x();
        double y = aNode.y();
        double z = aNode.y();

        double Ti;

        aAnaTemperature.steadyStateHeatConduction3D(x, y, z, Ti);

        sumError += (T.u(i) - Ti) * (T.u(i) - Ti);

    }

    sumError /= T.mesh().nbNodes();

    std::cout << "Error 3: " << sumError << std::endl;

}

int main(int argc, char** argv)
{

    steadyHeatConduction1D();
    steadyHeatConduction2D();
    steadyHeatConduction3D();

}
