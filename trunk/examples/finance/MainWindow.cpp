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

MainWindow::MainWindow(QWidget *parent)
{

    Ui::MainWindow::setupUi(this);

    QObject::connect(&m_blackScholes, SIGNAL(plotChanged(Integer, QVector<double>, QVector<double>)), this, SLOT(updatePlot(Integer, QVector<double>, QVector<double>)));

    customPlot->xAxis->setLabel("Stock price");
    customPlot->yAxis->setLabel("Option value");

    customPlot->xAxis->setRange(txtMinSharePrice->text().toDouble(), 1.1*txtMaxSharePrice->text().toDouble());
    customPlot->yAxis->setRange(0, 1.1*std::max(txtStrikePrice->text().toDouble(), txtSharePrice->text().toDouble()));
    customPlot->replot();
    customPlot->repaint();

}

void MainWindow::updatePlot(Integer plotId, QVector<double> x, QVector<double> y)
{

    customPlot->graph(plotId)->setData(x, y);
    customPlot->replot();
    customPlot->repaint();

}

void MainWindow::on_btnCalculate_clicked()
{

    m_blackScholes.setInterestRate(txtInterestRate->text().toDouble());
    m_blackScholes.setStrikePrice(txtStrikePrice->text().toDouble());
    m_blackScholes.setVolatility(txtVolatility->text().toDouble());
    m_blackScholes.setMaturity(txtMaturity->text().toDouble());
    m_blackScholes.setInitialPrice(txtSharePrice->text().toDouble());
    m_blackScholes.setInterval(txtMinSharePrice->text().toDouble(), txtMaxSharePrice->text().toDouble());

    m_blackScholes.enableDiffTerm(chkDiffTerm->isChecked());
    m_blackScholes.enableAdvectTerm(chkAdvectTerm->isChecked());
    m_blackScholes.enableReactTerm(chkReactTerm->isChecked());

    customPlot->xAxis->setRange(txtMinSharePrice->text().toDouble(), 1.1*txtMaxSharePrice->text().toDouble());
    customPlot->yAxis->setRange(0, 1.1*std::max(txtStrikePrice->text().toDouble(), txtSharePrice->text().toDouble()));
    customPlot->replot();
    customPlot->repaint();

    std::string strFVMResult = "fvm_black_scholes.dat";
    std::string strFEMResult = "fem_black_scholes.dat";
    std::string strFDMResult = "fdm_black_scholes.dat";

    customPlot->legend->setVisible(true);

    std::cout << "**** Black scholes equation - FVM ****" << std::endl;    
    if (customPlot->graphCount() == 0)
        customPlot->addGraph();

    customPlot->graph(0)->setName("FVM");
    customPlot->graph(0)->setPen(QPen(Qt::blue, 2.0));
    m_blackScholes.setPlotId(0);
    m_blackScholes.calculate(DM_FVM, strFVMResult);

    std::cout << "**** Black scholes equation - FEM ****" << std::endl;
    if (customPlot->graphCount() == 1)
        customPlot->addGraph();

    customPlot->graph(1)->setName("FEM");
    customPlot->graph(1)->setPen(QPen(Qt::red, 2.0));
    m_blackScholes.setPlotId(1);
    m_blackScholes.calculate(DM_FEM, strFEMResult);

    std::cout << "**** Black scholes equation - FDM ****" << std::endl;
    if (customPlot->graphCount() == 2)
        customPlot->addGraph();

    customPlot->graph(2)->setName("FDM");
    customPlot->graph(2)->setPen(QPen(Qt::green, 2.0));
    m_blackScholes.setPlotId(2);
    m_blackScholes.calculate(DM_FDM, strFDMResult);

}
