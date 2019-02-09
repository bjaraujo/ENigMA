
#pragma once

#include "ui_MainWindow.h"

#include "MshBasicMesher.hpp"
#include "PdeEquation.hpp"
#include "PosGnuplot.hpp"
#include "PosGmsh.hpp"
#include "AnaFlow.hpp"
#include "SphParticles.hpp"

class MainWindow : public QMainWindow, private Ui::MainWindow
{
    Q_OBJECT

private:
    double m_length;
    Integer m_step;
    Integer m_plotId;
    QList<QCPScatterStyle> m_pointStyles;
    QList<QPen> m_lineStyles;

    QCPTextElement* m_plotTitle;

    void solveBurgersEquation1D();

    void resetPlots();
    void addPlot(CPdeField<double> T, QString plotName, const bool bLine = false);
    double phi(const double x, const double t, const double mu);

public:
    MainWindow(QWidget *parent = 0);

private slots:
    void on_btnCalculate_clicked();

};

