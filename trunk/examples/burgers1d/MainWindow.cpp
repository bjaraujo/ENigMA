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

    m_length = 2 * M_PI;

    customPlot->xAxis->setLabel("x");
    customPlot->yAxis->setLabel("T");

    customPlot->xAxis->setRange(0.0, m_length);
    customPlot->yAxis->setRange(0.0, 10);

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

double MainWindow::phi(const double x, const double t, const double mu)
{

    return exp(-(x - 4 * t) * (x - 4 * t) / (4 * mu * (t + 1))) + exp(-(x - 4 * t - m_length) * (x - 4 * t - m_length) / (4 * mu * (t + 1)));

}

void MainWindow::solveBurgersEquation1D()
{

    // One-dimensional unsteady heat convection

    std::cout << "**** One-dimensional burgers equation ****" << std::endl;

    m_plotTitle->setText("One-dimensional burgers equation");

    CGeoCoordinate<double> aPoint1(0.0, 0.0, 0.0);
    CGeoCoordinate<double> aPoint2(m_length, 0.0, 0.0);

    CGeoLine<double> aLine(aPoint1, aPoint2);

    const Integer nu = 1000;

    CMshBasicMesher<double> aBasicMesher;

    aBasicMesher.generate(aLine, nu);
    CMshMesh<double> aLineMesh = aBasicMesher.mesh();

    double dt = 2.5e-4;
    double dx = m_length / nu;

    Integer nIter = 400;
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

        ua.u.resize(ua.mesh().nbNodes());
        for (Integer i = 0; i < ua.mesh().nbNodes(); ++i)
        {

            Integer aNodeId = ua.mesh().nodeId(i);
            CMshNode<double> aNode = ua.mesh().node(aNodeId);
            double x = aNode.x();

            double ui;
            aAnaFlow.burgersEquation(x, dx, 0.1, m_length, mu, ui);

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
        u.setSimulationType(ST_THERMAL);
        u.setNbDofs(1);

        u.u.resize(u.mesh().nbNodes());
        for (Integer i = 0; i < u.mesh().nbNodes(); ++i)
        {

            Integer aNodeId = u.mesh().nodeId(i);
            CMshNode<double> aNode = u.mesh().node(aNodeId);
            double x = aNode.x();

            double ui;
            aAnaFlow.burgersEquation(x, dx, 0, m_length, mu, ui);

            u.u(i) = ui;

        }

        for (Integer i = 0; i < nIter; ++i)
        {

            std::cout << i << std::endl;

            CPdeEquation<double> aPdeEquation(1.0 / dt * ddt<double>(u) + u.u * divergence<double>(u) - mu * laplacian<double>(u) = 0);

            aPdeEquation.solve(u);

            // Cyclic BC
            u.u(0) = u.u(u.mesh().nbNodes() - 1);

        }

        this->addPlot(u, "FDM");

    }

    /*
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

            CPdeEquation<double> aPdeEquation(1.0 / dt * ddt<double>(T) +u * divergence<double>(T) = 0);

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

            CPdeEquation<double> aPdeEquation(1.0 / dt * ddt<double>(T) +u * divergence<double>(T) = 0);

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
    */

    std::cout << "Done." << std::endl;

}

void MainWindow::on_btnCalculate_clicked()
{

    this->resetPlots();

    this->solveBurgersEquation1D();

    customPlot->replot();
    customPlot->repaint();

    m_step++;

}
