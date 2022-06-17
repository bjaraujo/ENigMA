
#pragma once

#include "ui_MainWindow.h"

#include "MshBasicMesher.hpp"
#include "PdeEquation.hpp"
#include "PosGnuplot.hpp"
#include "PosGmsh.hpp"
#include "AnaTemperature.hpp"
#include "SphParticles.hpp"

class MainWindow : public QMainWindow, private Ui::MainWindow
{
    Q_OBJECT

private:
    Integer m_step;
    Integer m_plotId;
    QList<QCPScatterStyle> m_pointStyles;
    QList<QPen> m_lineStyles;

    QCPTextElement* m_plotTitle;

    void solveHeatConduction1D1();
    void solveHeatConduction1D2();
    void unsteadyHeatConduction1D();
    void unsteadyHeatConvection1D();

    void resetPlots();
    void addPlot(CPdeField<double>& T, QString plotName, const bool bLine = false);

public:
    MainWindow(QWidget *parent = 0);

private slots:
    void on_btnCalculate_clicked();

};

