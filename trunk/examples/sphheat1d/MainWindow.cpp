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

void MainWindow::addPlot(CPdeField<double>& T, QString plotName, const bool bLine)
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

void MainWindow::solveHeatConduction1D1()
{

    std::cout << "**** One-dimensional steady-state heat conduction ****" << std::endl;

    m_plotTitle->setText("One-dimensional steady-state heat conduction");

    CGeoCoordinate<double> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<double> aPoint2(1.0, 0.0, 0.0);

    CGeoLine<double> aLine(aPoint1, aPoint2);

    const Integer nu = 100;

    CMshBasicMesher<double> aBasicMesher;
    
    aBasicMesher.generate(aLine, nu);
    CMshMesh<double> aLineMesh = aBasicMesher.mesh();

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

        Ta.init();

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

    // SPH 1
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

        double dx = 1.0 / nu;

        // Side 1
        for (Integer i = 0; i < 20; ++i)
        {
            aNode << 0.0 - (double)(i + 1) * dx, 0.0, 0.0;
            T.mesh().addNode(T.mesh().nextNodeId(), aNode);
        }

        // Side 2
        for (Integer i = 0; i < 20; ++i)
        {
            aNode << 1.0 + (double)(i + 1) * dx, 0.0, 0.0;
            T.mesh().addNode(T.mesh().nextNodeId(), aNode);
        }

        T.u.resize(T.mesh().nbNodes());
        T.u.fill(0.0);
        for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
        {

            Integer aNodeId = T.mesh().nodeId(i);
            CMshNode<double> aNode = T.mesh().node(aNodeId);

            if (aNode.x() <= 0.0)
                T.setFixedValue(i, 0.0);

            if (aNode.x() >= 1.0)
                T.setFixedValue(i, 1.0);

            T.u(i) = 0.0;

        }

        double h = 0.25;
        double diff = 1.0;
        double mass = 1.0;
        double rho = 1.0;

        {
            CSphSpiky<double> aKernel(1);
            CSphParticles<double> sParticles(aKernel);

            T.u.fill(0.0);
            sParticles.init(T, mass, rho, diff, h, dt);

            for (Integer i = 0; i < nIter; ++i)
            {
                sParticles.solve(T);
            }

            this->addPlot(T, "SPH - Spiky");
        }

        {
            CSphGaussian<double> aKernel(1);
            CSphParticles<double> sParticles(aKernel);

            T.u.fill(0.0);
            sParticles.init(T, mass, rho, diff, h, dt);

            for (Integer i = 0; i < nIter; ++i)
            {
                sParticles.solve(T);
            }

            this->addPlot(T, "SPH - Gaussian");
        }

        {
            CSphConvex<double> aKernel(1);
            CSphParticles<double> sParticles(aKernel);

            T.u.fill(0.0);
            sParticles.init(T, mass, rho, diff, h, dt);

            for (Integer i = 0; i < nIter; ++i)
            {
                sParticles.solve(T);
            }

            this->addPlot(T, "SPH - Convex");
        }

        {
            CSphQuintic<double> aKernel(1);
            CSphParticles<double> sParticles(aKernel);

            T.u.fill(0.0);
            sParticles.init(T, mass, rho, diff, h, dt);

            for (Integer i = 0; i < nIter; ++i)
            {
                sParticles.solve(T);
            }

            this->addPlot(T, "SPH - Quintic");
        }

        {
            CSphCubicSpline<double> aKernel(1);
            CSphParticles<double> sParticles(aKernel);

            T.u.fill(0.0);
            sParticles.init(T, mass, rho, diff, h, dt);

            for (Integer i = 0; i < nIter; ++i)
            {
                sParticles.solve(T);
            }

            this->addPlot(T, "SPH - CubicSpline");
        }
    }

    std::cout << "Done." << std::endl;

}

void MainWindow::solveHeatConduction1D2()
{

    std::cout << "**** One-dimensional steady-state heat conduction ****" << std::endl;

    m_plotTitle->setText("One-dimensional steady-state heat conduction");

    CGeoCoordinate<double> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<double> aPoint2(1.0, 0.0, 0.0);

    CGeoLine<double> aLine(aPoint1, aPoint2);

    const Integer nu = 100;

    CMshBasicMesher<double> aBasicMesher;

    aBasicMesher.generate(aLine, nu);
    CMshMesh<double> aLineMesh = aBasicMesher.mesh();

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

        Ta.init();

        CAnaTemperature<double> aAnaTemperature;

        for (Integer i = 0; i < Ta.mesh().nbNodes(); ++i)
        {

            Integer aNodeId = Ta.mesh().nodeId(i);
            CMshNode<double> aNode = Ta.mesh().node(aNodeId);
            double x = aNode.x();

            double Ti;
            aAnaTemperature.steadyStateHeatConduction1D(1.0 - x, Ti);

            Ta.u(i) = Ti;

        }

        this->addPlot(Ta, "Analytical", true);

    }

    // SPH 1
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

        double dx = 1.0 / nu;

        // Side 1
        for (Integer i = 0; i < 20; ++i)
        {
            aNode << 0.0 - (double)(i + 1) * dx, 0.0, 0.0;
            T.mesh().addNode(T.mesh().nextNodeId(), aNode);
        }

        // Side 2
        for (Integer i = 0; i < 20; ++i)
        {
            aNode << 1.0 + (double)(i + 1) * dx, 0.0, 0.0;
            T.mesh().addNode(T.mesh().nextNodeId(), aNode);
        }

        T.u.resize(T.mesh().nbNodes());
        T.u.fill(0.0);
        for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
        {

            Integer aNodeId = T.mesh().nodeId(i);
            CMshNode<double> aNode = T.mesh().node(aNodeId);

            if (aNode.x() <= 0.0)
                T.setFixedValue(i, 1.0);

            if (aNode.x() >= 1.0)
                T.setFixedValue(i, 0.0);

            T.u(i) = 0.0;

        }

        double h = 0.25;
        double diff = 1.0;
        double mass = 1.0;
        double rho = 1.0;

        {
            CSphSpiky<double> aKernel(1);
            CSphParticles<double> sParticles(aKernel);

            T.u.fill(0.0);
            sParticles.init(T, mass, rho, diff, h, dt);

            for (Integer i = 0; i < nIter; ++i)
            {
                sParticles.solve(T);
            }

            this->addPlot(T, "SPH - Spiky");
        }

        {
            CSphGaussian<double> aKernel(1);
            CSphParticles<double> sParticles(aKernel);

            T.u.fill(0.0);
            sParticles.init(T, mass, rho, diff, h, dt);

            for (Integer i = 0; i < nIter; ++i)
            {
                sParticles.solve(T);
            }

            this->addPlot(T, "SPH - Gaussian");
        }

        {
            CSphConvex<double> aKernel(1);
            CSphParticles<double> sParticles(aKernel);

            T.u.fill(0.0);
            sParticles.init(T, mass, rho, diff, h, dt);

            for (Integer i = 0; i < nIter; ++i)
            {
                sParticles.solve(T);
            }

            this->addPlot(T, "SPH - Convex");
        }

        {
            CSphQuintic<double> aKernel(1);
            CSphParticles<double> sParticles(aKernel);

            T.u.fill(0.0);
            sParticles.init(T, mass, rho, diff, h, dt);

            for (Integer i = 0; i < nIter; ++i)
            {
                sParticles.solve(T);
            }

            this->addPlot(T, "SPH - Quintic");
        }

        {
            CSphCubicSpline<double> aKernel(1);
            CSphParticles<double> sParticles(aKernel);

            T.u.fill(0.0);
            sParticles.init(T, mass, rho, diff, h, dt);

            for (Integer i = 0; i < nIter; ++i)
            {
                sParticles.solve(T);
            }

            this->addPlot(T, "SPH - CubicSpline");
        }
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

    const Integer nu = 100;

    CMshBasicMesher<double> aBasicMesher;

    aBasicMesher.generate(aLine, nu);
    CMshMesh<double> aLineMesh = aBasicMesher.mesh();

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

        Ta.init();

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

        double dx = 1.0 / nu;

        // Side 1
        for (Integer i = 0; i < 20; ++i)
        {
            aNode << 1.0 + (double)(i + 1) * dx, 0.0, 0.0;
            T.mesh().addNode(T.mesh().nextNodeId(), aNode);
        }

        T.u.resize(T.mesh().nbNodes());
        for (Integer i = 0; i < T.mesh().nbNodes(); ++i)
        {

            Integer aNodeId = T.mesh().nodeId(i);
            CMshNode<double> aNode = T.mesh().node(aNodeId);

            if (aNode.x() >= 1.0)
                T.setFixedValue(i, 0.0);

            T.u(i) = 1.0;

        }

        double h = 0.25;
        double diff = 1.0;
        double mass = 1.0;
        double rho = 1.0;

        {
            CSphSpiky<double> aKernel(1);
            CSphParticles<double> sParticles(aKernel);

            T.u.fill(1.0);
            sParticles.init(T, mass, rho, diff, h, dt);

            for (Integer i = 0; i < nIter; ++i)
            {
                sParticles.solve(T);
            }

            this->addPlot(T, "SPH - Spiky");
        }

        {
            CSphGaussian<double> aKernel(1);
            CSphParticles<double> sParticles(aKernel);

            T.u.fill(1.0);
            sParticles.init(T, mass, rho, diff, h, dt);

            for (Integer i = 0; i < nIter; ++i)
            {
                sParticles.solve(T);
            }

            this->addPlot(T, "SPH - Gaussian");
        }

        {
            CSphConvex<double> aKernel(1);
            CSphParticles<double> sParticles(aKernel);

            T.u.fill(1.0);
            sParticles.init(T, mass, rho, diff, h, dt);

            for (Integer i = 0; i < nIter; ++i)
            {
                sParticles.solve(T);
            }

            this->addPlot(T, "SPH - Convex");
        }

        {
            CSphQuintic<double> aKernel(1);
            CSphParticles<double> sParticles(aKernel);

            T.u.fill(1.0);
            sParticles.init(T, mass, rho, diff, h, dt);

            for (Integer i = 0; i < nIter; ++i)
            {
                sParticles.solve(T);
            }

            this->addPlot(T, "SPH - Quintic");
        }

        {
            CSphCubicSpline<double> aKernel(1);
            CSphParticles<double> sParticles(aKernel);

            T.u.fill(1.0);
            sParticles.init(T, mass, rho, diff, h, dt);

            for (Integer i = 0; i < nIter; ++i)
            {
                sParticles.solve(T);
            }

            this->addPlot(T, "SPH - Cubic Spline");
        }

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

    const Integer nu = 100;

    CMshBasicMesher<double> aBasicMesher;

    aBasicMesher.generate(aLine, nu);
    CMshMesh<double> aLineMesh = aBasicMesher.mesh();

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

        sParticles.init(T, mass, 0.0, h, dt, true);

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

    if (m_step % 4 == 0)
        this->solveHeatConduction1D1();
    else if (m_step % 4 == 1)
        this->solveHeatConduction1D2();
    else if (m_step % 4 == 2)
        this->unsteadyHeatConduction1D();
    else if (m_step % 4 == 3)
        this->unsteadyHeatConvection1D();

    customPlot->replot();
    customPlot->repaint();

    m_step++;

}
