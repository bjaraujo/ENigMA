// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#include <QtGui>
#include "MainWindow.h"

using namespace ENigMA::mesh;
using namespace ENigMA::post;
using namespace ENigMA::analytical;
using namespace ENigMA::sph;

MainWindow::MainWindow(QWidget *parent) : m_step(0), m_plotId(0)
{

    Ui::MainWindow::setupUi(this);

    customPlot->xAxis->setLabel("x");
    customPlot->yAxis->setLabel("T");

    customPlot->xAxis->setRange(0.0, 1.0);
    customPlot->yAxis->setRange(0.0, 1.0);

    m_plotTitle = new QCPTextElement(customPlot, "", QFont("sans", 12, QFont::Bold));

    customPlot->plotLayout()->insertRow(0);
    customPlot->plotLayout()->addElement(0, 0, m_plotTitle);

    m_pointStyles.append(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::black, Qt::white, 6));
    m_pointStyles.append(QCPScatterStyle(QCPScatterStyle::ssSquare, Qt::red, Qt::white, 6));
    m_pointStyles.append(QCPScatterStyle(QCPScatterStyle::ssDiamond, Qt::green, Qt::white, 6));
    m_pointStyles.append(QCPScatterStyle(QCPScatterStyle::ssPlus, Qt::blue, Qt::white, 6));
    m_pointStyles.append(QCPScatterStyle(QCPScatterStyle::ssStar, Qt::darkYellow, Qt::white, 6));

    m_lineStyles.append(QPen(Qt::black));

}

void MainWindow::resetPlots()
{

    customPlot->clearGraphs();

    customPlot->legend->setVisible(true);

    customPlot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignLeft | Qt::AlignTop);

    customPlot->replot();
    customPlot->repaint();

    m_plotId = 0;

}

void MainWindow::addPlot(CPdeField<double> T, QString plotName, const bool bLine)
{

    QCPGraph* aGraph = customPlot->addGraph();

    aGraph->setName(plotName);

    if (bLine)
    {
        aGraph->setLineStyle(QCPGraph::lsLine);
        aGraph->setPen(m_lineStyles[m_plotId % 1]);
    }
    else
    {
        aGraph->setLineStyle(QCPGraph::lsNone);
        aGraph->setScatterStyle(m_pointStyles[m_plotId % 5]);
    }

    QVector<double> x, y;

    if (T.discretLocation() == DL_NODE)
    {

        for (Integer i = 0; i < T.mesh().nbNodes(); i++)
        {

            Integer id = T.mesh().nodeId(i);

            x.push_back(T.mesh().node(id).x());
            y.push_back(T.u(i));

        }

    }
    else if (T.discretLocation() == DL_ELEMENT_CENTER)
    {

        for (Integer i = 0; i < T.mesh().nbElements(); i++)
        {

            Integer id = T.mesh().elementId(i);

            CGeoCoordinate<double> aCentroid = T.mesh().elementCentroid(id);

            x.push_back(aCentroid.x());
            y.push_back(T.u(i));

        }

    }

    aGraph->setData(x, y);

    customPlot->replot();
    customPlot->repaint();

    m_plotId++;

}

void MainWindow::solveHeatConduction1D()
{

    std::cout << "**** One-dimensional steady-state heat conduction ****" << std::endl;

    m_plotTitle->setText("One-dimensional steady-state heat conduction");

    CGeoCoordinate<double> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<double> aPoint2(1.0, 0.0, 0.0);

    CGeoLine<double> aLine(aPoint1, aPoint2);

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

    const Integer nu = 50;
    const Integer nv = 1;
    const Integer nw = 1;

    CMshBasicMesher<double> aBasicMesher;
    
    aBasicMesher.generate(aLine, nu);
    CMshMesh<double> aLineMesh = aBasicMesher.mesh();

    aBasicMesher.generate(aHexahedron, nu, nv, nw);
    CMshMesh<double> aVolumeMesh = aBasicMesher.mesh();

    aVolumeMesh.generateFaces(1E-12);
    aVolumeMesh.calculateFaceCentroid();
    aVolumeMesh.calculateElementCentroid();

    double dt = 0.01;
    Integer nIter = 100;

    // Analytical
    {

        CPdeField<double> Ta;

        Ta.setMesh(aLineMesh);
        Ta.setDiscretLocation(DL_NODE);
        Ta.setDiscretOrder(DO_LINEAR);
        Ta.setSimulationType(ST_THERMAL);
        Ta.setNbDofs(1);

        Ta.u.resize(Ta.mesh().nbNodes());

        CAnaTemperature<double> aAnaTemperature;

        for (Integer i = 0; i < Ta.mesh().nbNodes(); ++i)
        {

            Integer aNodeId = Ta.mesh().nodeId(i);
            CMshNode<double> aNode = Ta.mesh().node(aNodeId);
            double x = aNode.x();

            double Ti;
            aAnaTemperature.steadyStateHeatConduction1D(x, Ti);

            Ta.u(i) = Ti;

        }

        this->addPlot(Ta, "Analytical", true);

    }

    // FDM
    {

        CPdeField<double> T;

        T.setMesh(aLineMesh);
        T.setDiscretMethod(DM_FDM);
        T.setDiscretLocation(DL_NODE);
        T.setDiscretOrder(DO_LINEAR);
        T.setSimulationType(ST_THERMAL);
        T.setNbDofs(1);

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

        for (Integer i = 0; i < nIter; ++i)
        {

            CPdeEquation<double> aPdeEquation(1.0 / dt * ddt<double>(T) - laplacian<double>(T) = 0);

            aPdeEquation.solve(T);

        }

        this->addPlot(T, "FDM");

    }

    // FEM
    {

        CPdeField<double> T;

        T.setMesh(aLineMesh);
        T.setDiscretMethod(DM_FEM);
        T.setDiscretLocation(DL_NODE);
        T.setDiscretOrder(DO_LINEAR);
        T.setSimulationType(ST_THERMAL);
        T.setNbDofs(1);

        T.u.resize(T.mesh().nbNodes());
        for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
        {

            Integer aNodeId = T.mesh().nodeId(i);
            CMshNode<double> aNode = T.mesh().node(aNodeId);

            if (fabs(aNode.x() - 0.0) < 1E-6)
                T.setFixedValue(i, 0.0);

            if (fabs(aNode.x() - 1.0) < 1E-6)
                T.setFixedValue(i, 1.0);

        }

        CPdeEquation<double> aPdeEquation(laplacian<double>(T) = 0);

        aPdeEquation.setElimination(T);

        aPdeEquation.solve(T);

        this->addPlot(T, "FEM");

    }

    // FVM
    {

        CPdeField<double> T;

        T.setMesh(aVolumeMesh);
        T.setDiscretMethod(DM_FVM);
        T.setDiscretLocation(DL_ELEMENT_CENTER);
        T.setDiscretOrder(DO_LINEAR);
        T.setSimulationType(ST_THERMAL);
        T.setNbDofs(1);

        T.u.resize(T.mesh().nbElements());
        for (Integer i = 0; i < T.mesh().nbElements(); ++i)
            T.u(i) = 0.0;

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

        for (Integer i = 0; i < nIter; ++i)
        {

            CPdeEquation<double> aPdeEquation(1.0 / dt * ddt<double>(T) - laplacian<double>(T) = 0);

            aPdeEquation.solve(T);

        }

        this->addPlot(T, "FVM");

    }

    // SPH
    {

        CPdeField<double> T;

        T.setMesh(aLineMesh);
        T.setDiscretMethod(DM_SPH);
        T.setDiscretLocation(DL_NODE);
        T.setDiscretOrder(DO_LINEAR);
        T.setSimulationType(ST_THERMAL);
        T.setNbDofs(1);

        // Ghost points
        CMshNode<double> aNode;

        // Side 1
        for (Integer i = 0; i < 5; ++i)
        {
            aNode << 0.0 - (double)(i + 1) / nu, 0.0, 0.0;
            T.mesh().addNode(T.mesh().nextNodeId(), aNode);
        }

        // Side 2
        for (Integer i = 0; i < 5; ++i)
        {
            aNode << 1.0 + (double)(i + 1) / nu, 0.0, 0.0;
            T.mesh().addNode(T.mesh().nextNodeId(), aNode);
        }

        T.u.resize(T.mesh().nbNodes());
        for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
        {

            Integer aNodeId = T.mesh().nodeId(i);
            CMshNode<double> aNode = T.mesh().node(aNodeId);

            if (aNode.x() < 0.0 + 1E-6)
                T.setFixedValue(i, 0.0);

            if (aNode.x() > 1.0 - 1E-6)
                T.setFixedValue(i, 1.0);

            T.u(i) = 0.0;

        }

        CSphSpiky<double> aKernel(1);
        CSphParticles<double> sParticles(aKernel);

        double h = 0.25;
        double diff = 1.0;
        double mass = 1.0;
        double rho = 1.0;

        sParticles.init(T, mass, rho, diff, h, dt);

        for (Integer i = 0; i < nIter; ++i)
        {
            sParticles.solve(T);
        }

        this->addPlot(T, "SPH");

    }

    std::cout << "Done." << std::endl;

}

void MainWindow::unsteadyHeatConduction1D()
{

    // One-dimensional unsteady heat conduction

    std::cout << "**** One-dimensional unsteady heat conduction ****" << std::endl;

    m_plotTitle->setText("One-dimensional unsteady heat conduction");

    CGeoCoordinate<double> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<double> aPoint2(1.0, 0.0, 0.0);

    CGeoLine<double> aLine(aPoint1, aPoint2);

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

    const Integer nu = 50;
    const Integer nv = 1;
    const Integer nw = 1;

    CMshBasicMesher<double> aBasicMesher;

    aBasicMesher.generate(aLine, nu);
    CMshMesh<double> aLineMesh = aBasicMesher.mesh();

    aBasicMesher.generate(aHexahedron, nu, nv, nw);
    CMshMesh<double> aVolumeMesh = aBasicMesher.mesh();

    aVolumeMesh.generateFaces(1E-12);
    aVolumeMesh.calculateFaceCentroid();
    aVolumeMesh.calculateElementCentroid();

    double rho = 1.0;   // density
    double Cp = 1.0;    // specific heat
    double k = 1.0;     // conductivity

    double dt = 0.001;
    Integer nIter = 100;

    // Analytical
    {

        CPdeField<double> Ta;

        Ta.setMesh(aLineMesh);
        Ta.setDiscretLocation(DL_NODE);
        Ta.setDiscretOrder(DO_LINEAR);
        Ta.setSimulationType(ST_THERMAL);
        Ta.setNbDofs(1);

        Ta.u.resize(Ta.mesh().nbNodes());

        CAnaTemperature<double> aAnaTemperature;

        for (Integer i = 0; i < Ta.mesh().nbNodes(); ++i)
        {

            Integer aNodeId = Ta.mesh().nodeId(i);
            CMshNode<double> aNode = Ta.mesh().node(aNodeId);
            double x = aNode.x();

            double Ti;

            aAnaTemperature.transientHeatConduction1D(x, dt * nIter, 1.0, Ti);

            Ta.u(i) = Ti;

        }

        this->addPlot(Ta, "Analytical", true);

    }

    // FDM
    {

        CPdeField<double> T;

        T.setMesh(aLineMesh);
        T.setDiscretMethod(DM_FDM);
        T.setDiscretLocation(DL_NODE);
        T.setDiscretOrder(DO_LINEAR);
        T.setSimulationType(ST_THERMAL);
        T.setNbDofs(1);

        T.u.resize(T.mesh().nbNodes());
        for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
            T.u(i) = 0.0;

        for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
        {

            Integer aNodeId = T.mesh().nodeId(i);
            CMshNode<double> aNode = T.mesh().node(aNodeId);

            if (fabs(aNode.x() - 1.0) < 1E-6)
                T.setFixedValue(i, 0.0);

            T.u(i) = 1.0;

        }

        for (Integer i = 0; i < nIter; ++i)
        {

            CPdeEquation<double> aPdeEquation(rho * Cp / dt * ddt<double>(T) - k * laplacian<double>(T) = 0);

            aPdeEquation.solve(T);

        }

        this->addPlot(T, "FDM");

    }

    // FEM
    {

        CPdeField<double> T;

        T.setMesh(aLineMesh);
        T.setDiscretMethod(DM_FEM);
        T.setDiscretLocation(DL_NODE);
        T.setDiscretOrder(DO_LINEAR);
        T.setSimulationType(ST_THERMAL);
        T.setNbDofs(1);

        T.u.resize(T.mesh().nbNodes());
        for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
        {

            Integer aNodeId = T.mesh().nodeId(i);
            CMshNode<double> aNode = T.mesh().node(aNodeId);

            if (fabs(aNode.x() - 1.0) < 1E-6)
                T.setFixedValue(i, 0.0);

            T.u(i) = 1.0;

        }

        for (Integer i = 0; i < nIter; ++i)
        {

            CPdeEquation<double> aPdeEquation(rho * Cp / dt * ddt<double>(T) - k * laplacian<double>(T) = 0);

            aPdeEquation.setElimination(T);

            aPdeEquation.solve(T);

        }

        this->addPlot(T, "FEM");

    }

    // FVM
    {

        CPdeField<double> T;

        T.setMesh(aVolumeMesh);
        T.setDiscretMethod(DM_FVM);
        T.setDiscretLocation(DL_ELEMENT_CENTER);
        T.setDiscretOrder(DO_LINEAR);
        T.setSimulationType(ST_THERMAL);
        T.setNbDofs(1);

        T.u.resize(T.mesh().nbElements());
        for (Integer i = 0; i < T.mesh().nbElements(); ++i)
            T.u(i) = 1.0;

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

        for (Integer i = 0; i < nIter; ++i)
        {

            CPdeEquation<double> aPdeEquation(rho * Cp / dt * ddt<double>(T) - k * laplacian<double>(T) = 0);

            aPdeEquation.solve(T);

        }

        this->addPlot(T, "FVM");

    }

    // SPH
    {

        CPdeField<double> T;

        T.setMesh(aLineMesh);
        T.setDiscretMethod(DM_SPH);
        T.setDiscretLocation(DL_NODE);
        T.setDiscretOrder(DO_LINEAR);
        T.setSimulationType(ST_THERMAL);
        T.setNbDofs(1);

        CMshNode<double> aNode;

        // Side 1
        for (Integer i = 0; i < 5; ++i)
        {
            aNode << 1.0 + (double)(i + 1) / nu, 0.0, 0.0;
            T.mesh().addNode(T.mesh().nextNodeId(), aNode);
        }

        T.u.resize(T.mesh().nbNodes());
        for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
        {

            Integer aNodeId = T.mesh().nodeId(i);
            CMshNode<double> aNode = T.mesh().node(aNodeId);

            if (aNode.x() > 1.0 - 1E-6)
                T.setFixedValue(i, 0.0);

            T.u(i) = 1.0;

        }

        CSphSpiky<double> aKernel(1);
        CSphParticles<double> sParticles(aKernel);

        double h = 0.25;
        double diff = 1.0;
        double mass = 1.0;
        double rho = 1.0;

        sParticles.init(T, mass, rho, diff, h, dt);

        for (Integer i = 0; i < nIter; ++i)
        {
            sParticles.solve(T);
        }

        this->addPlot(T, "SPH");

    }

    std::cout << "Done." << std::endl;

}

void MainWindow::unsteadyHeatConvection1D()
{

    // One-dimensional unsteady heat convection

    std::cout << "**** One-dimensional unsteady heat convection ****" << std::endl;

    m_plotTitle->setText("One-dimensional unsteady heat convection");

    CGeoCoordinate<double> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<double> aPoint2(1.0, 0.0, 0.0);

    CGeoLine<double> aLine(aPoint1, aPoint2);

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

    const Integer nu = 100;
    const Integer nv = 1;
    const Integer nw = 1;

    CMshBasicMesher<double> aBasicMesher;

    aBasicMesher.generate(aLine, nu);
    CMshMesh<double> aLineMesh = aBasicMesher.mesh();

    aBasicMesher.generate(aHexahedron, nu, nv, nw);
    CMshMesh<double> aVolumeMesh = aBasicMesher.mesh();

    aVolumeMesh.generateFaces(1E-12);
    aVolumeMesh.calculateFaceCentroid();
    aVolumeMesh.calculateElementCentroid();

    // Cyclic
    Integer aFaceId1, aFaceId2;
    for (Integer i = 0; i < aVolumeMesh.nbFaces(); ++i)
    {

        Integer aFaceId = aVolumeMesh.faceId(i);

        if (fabs(aVolumeMesh.faceCentroid(aFaceId).x() - 0.0) < 1E-6)
            aFaceId1 = aFaceId;

        if (fabs(aVolumeMesh.faceCentroid(aFaceId).x() - 1.0) < 1E-6)
            aFaceId2 = aFaceId;

    }
    aVolumeMesh.face(aFaceId1).setPairFaceId(aFaceId2);
    aVolumeMesh.face(aFaceId2).setPairFaceId(aFaceId1);

    double dt = 0.005;
    Integer nIter = 400;
    double time = dt * nIter;

    double u = 1.0;

    // Analytical
    {

        CPdeField<double> Ta;

        Ta.setMesh(aLineMesh);
        Ta.setDiscretLocation(DL_NODE);
        Ta.setDiscretOrder(DO_LINEAR);
        Ta.setSimulationType(ST_THERMAL);
        Ta.setNbDofs(1);

        Ta.u.resize(Ta.mesh().nbNodes());
        for (Integer i = 0; i < Ta.mesh().nbNodes(); ++i)
        {

            Integer aNodeId = Ta.mesh().nodeId(i);
            CMshNode<double> aNode = Ta.mesh().node(aNodeId);
            double x = aNode.x();

            double p = (x - 0.5) / 0.08;

            Ta.u(i) = exp(-0.5*p*p);

        }

        this->addPlot(Ta, "Analytical", true);

    }

    // FDM
    {

        CPdeField<double> T;

        T.setMesh(aLineMesh);
        T.setDiscretMethod(DM_FDM);
        T.setDiscretLocation(DL_NODE);
        T.setDiscretOrder(DO_LINEAR);
        T.setSimulationType(ST_THERMAL);
        T.setNbDofs(1);

        T.u.resize(T.mesh().nbNodes());
        for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
        {

            Integer aNodeId = T.mesh().nodeId(i);
            CMshNode<double> aNode = T.mesh().node(aNodeId);
            double x = aNode.x();

            double p = (x - 0.5) / 0.08;

            T.u(i) = exp(-0.5*p*p);

        }

        for (Integer i = 0; i < nIter; ++i)
        {

            CPdeEquation<double> aPdeEquation(1.0 / dt * ddt<double>(T) + u * divergence<double>(T) = 0);

            aPdeEquation.solve(T);

            // Zero gradient (Upwind)
            T.u(T.mesh().nbNodes() - 1) = T.u(T.mesh().nbNodes() - 2);

            // Cyclic BC
            T.u(0) = T.u(T.mesh().nbNodes() - 1);

        }

        this->addPlot(T, "FDM");

    }

    // FEM
    {

        CPdeField<double> T;

        T.setMesh(aLineMesh);
        T.setDiscretMethod(DM_FEM);
        T.setDiscretLocation(DL_NODE);
        T.setDiscretOrder(DO_LINEAR);
        T.setSimulationType(ST_THERMAL);
        T.setNbDofs(1);

        T.u.resize(T.mesh().nbNodes());
        for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
        {

            Integer aNodeId = T.mesh().nodeId(i);
            CMshNode<double> aNode = T.mesh().node(aNodeId);
            double x = aNode.x();

            double p = (x - 0.5) / 0.08;

            T.u(i) = exp(-0.5*p*p);

        }

        for (Integer i = 0; i < nIter; ++i)
        {

            CPdeEquation<double> aPdeEquation(1.0 / dt * ddt<double>(T) + u * divergence<double>(T) = 0);

            aPdeEquation.solve(T);

            // Cyclic BC
            T.u(0) = T.u(T.mesh().nbNodes() - 1);

        }

        this->addPlot(T, "FEM");

    }

    // FVM
    {

        CPdeField<double> T;

        T.setMesh(aVolumeMesh);
        T.setDiscretMethod(DM_FVM);
        T.setDiscretLocation(DL_ELEMENT_CENTER);
        T.setDiscretOrder(DO_LINEAR);
        T.setSimulationType(ST_THERMAL);
        T.setNbDofs(1);

        T.u.resize(T.mesh().nbElements());
        for (Integer i = 0; i < T.mesh().nbElements(); ++i)
        {

            Integer anElementId = T.mesh().elementId(i);
            double x = T.mesh().elementCentroid(anElementId).x();

            double p = (x - 0.5) / 0.08;

            T.u(i) = exp(-0.5*p*p);

        }

        for (Integer i = 0; i < nIter; ++i)
        {

            CPdeEquation<double> aPdeEquation(1.0 / dt * ddt<double>(T) + u * divergence<double>(T) = 0);

            aPdeEquation.solve(T);

            // Cyclic BC
            T.u(0) = T.u(T.mesh().nbElements() - 1);

        }

        this->addPlot(T, "FVM");

    }

    // SPH
    {

        CPdeField<double> T;

        T.setMesh(aLineMesh);
        T.setDiscretMethod(DM_SPH);
        T.setDiscretLocation(DL_NODE);
        T.setDiscretOrder(DO_LINEAR);
        T.setSimulationType(ST_THERMAL);
        T.setNbDofs(1);

        T.u.resize(T.mesh().nbNodes());
        for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
        {

            Integer aNodeId = T.mesh().nodeId(i);
            CMshNode<double> aNode = T.mesh().node(aNodeId);
            double x = aNode.x();

            double p = (x - 0.5) / 0.08;

            T.u(i) = exp(-0.5*p*p);

        }

        CSphSpiky<double> aKernel(1);
        CSphParticles<double> sParticles(aKernel);

        double h = 0.25;
        double mass = 1.0;
        double rho = 1.0;

        sParticles.init(T, mass, rho, 0.0, h, dt, true);

        const CGeoBoundingBox<double> aBoundingBox(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
        sParticles.setBoundary(aBoundingBox);

        const CGeoVector<double> aVelocity(u, 0, 0);
        sParticles.setInitialVelocity(T, aVelocity);

        for (Integer i = 0; i < nIter; ++i)
        {
            sParticles.solve(T);
        }

        this->addPlot(T, "SPH");

    }

    std::cout << "Done." << std::endl;

}

void MainWindow::on_btnCalculate_clicked()
{

    this->resetPlots();

    if (m_step % 3 == 0)
        this->solveHeatConduction1D();
    else if (m_step % 3 == 1)
        this->unsteadyHeatConduction1D();
    else if (m_step % 3 == 2)
        this->unsteadyHeatConvection1D();

    customPlot->replot();
    customPlot->repaint();

    m_step++;

}
