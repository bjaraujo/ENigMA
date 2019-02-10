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
    customPlot->yAxis->setLabel("u");

    m_plotTitle = new QCPTextElement(customPlot, "", QFont("sans", 12, QFont::Bold));

    customPlot->plotLayout()->insertRow(0);
    customPlot->plotLayout()->addElement(0, 0, m_plotTitle);

    m_pointStyles.append(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::black, Qt::white, 6));
    m_pointStyles.append(QCPScatterStyle(QCPScatterStyle::ssSquare, Qt::red, Qt::white, 6));
    m_pointStyles.append(QCPScatterStyle(QCPScatterStyle::ssDiamond, Qt::green, Qt::white, 6));
    m_pointStyles.append(QCPScatterStyle(QCPScatterStyle::ssPlus, Qt::blue, Qt::white, 6));
    m_pointStyles.append(QCPScatterStyle(QCPScatterStyle::ssStar, Qt::darkYellow, Qt::white, 6));

    m_lineStyles.append(QPen(Qt::black));

    for (int i = 0; i < 10; i++)
    {       
        m_lineStyles.append(QPen(QColor::fromHslF(1.0 - i * 0.1, 0.5, 0.5)));
    }

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

void MainWindow::addPlot(CPdeField<double>& u, QString plotName, const bool bLine)
{

    QCPGraph* aGraph = customPlot->addGraph();

    aGraph->setName(plotName);

    if (bLine)
    {
        aGraph->setLineStyle(QCPGraph::lsLine);
        aGraph->setPen(m_lineStyles[m_plotId % m_lineStyles.count()]);
    }
    else
    {
        aGraph->setLineStyle(QCPGraph::lsNone);
        aGraph->setScatterStyle(m_pointStyles[m_plotId % m_pointStyles.count()]);
    }

    QVector<double> x, y;

    if (u.discretLocation() == DL_NODE)
    {

        for (Integer i = 0; i < u.mesh().nbNodes(); i++)
        {

            Integer id = u.mesh().nodeId(i);

            x.push_back(u.mesh().node(id).x());
            y.push_back(u.u(i));

        }

    }
    else if (u.discretLocation() == DL_ELEMENT_CENTER)
    {

        for (Integer i = 0; i < u.mesh().nbElements(); i++)
        {

            Integer id = u.mesh().elementId(i);

            CGeoCoordinate<double> aCentroid = u.mesh().elementCentroid(id);

            x.push_back(aCentroid.x());
            y.push_back(u.u(i));

        }

    }

    aGraph->setData(x, y);

    customPlot->replot();
    customPlot->repaint();

    m_plotId++;

}

void MainWindow::solveViscousBurgersEquation1D()
{

    // One-dimensional viscous Burgers equation

    std::cout << "**** One-dimensional viscous Burgers equation ****" << std::endl;

    m_plotTitle->setText("One-dimensional viscous Burgers equation");

    double L = 2 * M_PI;

    customPlot->xAxis->setRange(0.0, L);
    customPlot->yAxis->setRange(0.0, 10);

    CGeoCoordinate<double> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<double> aPoint2(L, 0.0, 0.0);

    CGeoLine<double> aLine(aPoint1, aPoint2);

    CGeoCoordinate<double> aVertex1(+0.00, -0.05, -0.05);
    CGeoCoordinate<double> aVertex2(L, -0.05, -0.05);
    CGeoCoordinate<double> aVertex3(L, +0.05, -0.05);
    CGeoCoordinate<double> aVertex4(+0.00, +0.05, -0.05);
    CGeoCoordinate<double> aVertex5(+0.00, -0.05, +0.05);
    CGeoCoordinate<double> aVertex6(L, -0.05, +0.05);
    CGeoCoordinate<double> aVertex7(L, +0.05, +0.05);
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

    const Integer nu = 1000;
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

        if (fabs(aVolumeMesh.faceCentroid(aFaceId).x() - L) < 1E-6)
            aFaceId2 = aFaceId;

    }
    aVolumeMesh.face(aFaceId1).setPairFaceId(aFaceId2);
    aVolumeMesh.face(aFaceId2).setPairFaceId(aFaceId1);

    double dt = 2.5e-4;
    double dx = L / nu;

    Integer nIter = 5*400;
    double time = dt * nIter;

    const double mu = 0.05;

    CAnaFlow<double> aAnaFlow;

    // Analytical
    {

        CPdeField<double> ua;

        ua.setMesh(aLineMesh);
        ua.setDiscretLocation(DL_NODE);
        ua.setDiscretOrder(DO_LINEAR);
        ua.setSimulationType(ST_THERMAL);
        ua.setNbDofs(1);

        ua.init();

        for (Integer i = 0; i < ua.mesh().nbNodes(); ++i)
        {

            Integer aNodeId = ua.mesh().nodeId(i);
            CMshNode<double> aNode = ua.mesh().node(aNodeId);
            double x = aNode.x();

            double ui;
            aAnaFlow.viscousBurgersEquation(x, dx, 0.5, L, mu, ui);

            ua.u(i) = ui;
        }

        this->addPlot(ua, "Analytical", true);

    }

    // FDM
    {

        CPdeField<double> u;

        u.setMesh(aLineMesh);
        u.setDiscretMethod(DM_FDM);
        u.setDiscretLocation(DL_NODE);
        u.setDiscretOrder(DO_LINEAR);
        u.setSimulationType(ST_GENERIC);
        u.setNbDofs(1);

        u.init();

        for (Integer i = 0; i < u.mesh().nbNodes(); ++i)
        {

            Integer aNodeId = u.mesh().nodeId(i);
            CMshNode<double> aNode = u.mesh().node(aNodeId);
            double x = aNode.x();

            double ui;
            aAnaFlow.viscousBurgersEquation(x, dx, 0, L, mu, ui);

            u.u(i) = ui;
        }

        for (Integer j = 0; j < nIter; ++j)
        {
            std::cout << j << std::endl;

            CPdeEquation<double> aPdeEquation(1.0 / dt * ddt<double>(u) + u.u * divergence<double>(u) - mu * laplacian<double>(u) = 0);

            aPdeEquation.solve(u);

            // Cyclic BC
            u.u(0) = u.u(u.mesh().nbNodes() - 1);
        }

        this->addPlot(u, "FDM");
    }

    // FEM
    {

        CPdeField<double> u;

        u.setMesh(aLineMesh);
        u.setDiscretMethod(DM_FEM);
        u.setDiscretLocation(DL_NODE);
        u.setDiscretOrder(DO_LINEAR);
        u.setSimulationType(ST_GENERIC);
        u.setNbDofs(1);

        u.init();

        for (Integer i = 0; i < u.mesh().nbNodes(); ++i)
        {

            Integer aNodeId = u.mesh().nodeId(i);
            CMshNode<double> aNode = u.mesh().node(aNodeId);
            double x = aNode.x();

            double ui;
            aAnaFlow.viscousBurgersEquation(x, dx, 0, L, mu, ui);

            u.u(i) = ui;
        }

        for (Integer j = 0; j < nIter; ++j)
        {
            std::cout << j << std::endl;

            CPdeEquation<double> aPdeEquation(1.0 / dt * ddt<double>(u) + u.u * divergence<double>(u) - mu * laplacian<double>(u) = 0);

            aPdeEquation.solve(u);

            // Cyclic BC
            u.u(0) = u.u(u.mesh().nbNodes() - 1);

        }

        this->addPlot(u, "FEM");
    }

    // FVM
    {

        CPdeField<double> u;

        u.setMesh(aVolumeMesh);
        u.setDiscretMethod(DM_FVM);
        u.setDiscretLocation(DL_ELEMENT_CENTER);
        u.setDiscretOrder(DO_LINEAR);
        u.setSimulationType(ST_GENERIC);

        u.init();

        for (Integer i = 0; i < u.mesh().nbElements(); ++i)
        {

            Integer anElementId = u.mesh().elementId(i);
            double x = u.mesh().elementCentroid(anElementId).x();

            double ui;
            aAnaFlow.viscousBurgersEquation(x, dx, 0, L, mu, ui);

            u.u(i) = ui;
        }

        for (Integer j = 0; j < nIter; ++j)
        {
            std::cout << j << std::endl;

            CPdeEquation<double> aPdeEquation(1.0 / dt * ddt<double>(u) + u.u * divergence<double>(u) - mu * laplacian<double>(u) = 0);

            aPdeEquation.solve(u);

            // Cyclic BC
            u.u(0) = u.u(u.mesh().nbElements() - 1);
        }

        this->addPlot(u, "FVM");
    }

    std::cout << "Done." << std::endl;

}

void MainWindow::solveInviscidBurgersEquation1D()
{

    // One-dimensional inviscid Burgers equation

    std::cout << "**** One-dimensional inviscid Burgers equation ****" << std::endl;

    m_plotTitle->setText("One-dimensional inviscid Burgers equation");

    double L = 1.0;

    customPlot->xAxis->setRange(0.0, L);
    customPlot->yAxis->setRange(0.0, 1.5);

    CGeoCoordinate<double> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<double> aPoint2(L, 0.0, 0.0);

    CGeoLine<double> aLine(aPoint1, aPoint2);

    const Integer nu = 500;
    const Integer nv = 1;
    const Integer nw = 1;

    CMshBasicMesher<double> aBasicMesher;

    aBasicMesher.generate(aLine, nu);
    CMshMesh<double> aLineMesh = aBasicMesher.mesh();

    double dt = 5e-4;
    double dx = L / nu;

    Integer nIter = 1000;
    double time = dt * nIter;

    m_plotId++;

    const double b = 0.5;
    const double c = 0.1;

    CAnaFlow<double> aAnaFlow;

    // FDM
    {

        CPdeField<double> u;

        u.setMesh(aLineMesh);
        u.setDiscretMethod(DM_FDM);
        u.setDiscretLocation(DL_NODE);
        u.setDiscretOrder(DO_LINEAR);
        u.setSimulationType(ST_GENERIC);
        u.setNbDofs(1);

        u.init();

        for (Integer i = 0; i < nu; ++i)
        {
            double x = i * dx;
            u.u[i] = exp(-(x - b) * (x - b) / (2 * c * c));
        }

        double t = 0.0;

        for (Integer j = 0; j < nIter; ++j)
        {
            std::cout << j << std::endl;

            for (int k = 1; k < u.u.size() - 1; k++)
            {
                u.u[k] = u.u[k] - u.u[k] * (dt / dx)*(u.u[k] - u.u[k - 1]);
            }
            u.u[0] = u.u[nu - 1] = 0;

            if (j % 100 == 0)
            {
                this->addPlot(u, "t = " + QString::number(t, 'f', 3), true);
            }

            t += dt;
        }

        this->addPlot(u, "t = " + QString::number(t, 'f', 3), true);
    }

    std::cout << "Done." << std::endl;

}

void MainWindow::on_btnCalculate_clicked()
{

    this->resetPlots();

    //this->solveViscousBurgersEquation1D();
    this->solveInviscidBurgersEquation1D();

    customPlot->replot();
    customPlot->repaint();

    m_step++;

}
