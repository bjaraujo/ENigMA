
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
    Integer m_step;
    Integer m_plotId;
    QList<QCPScatterStyle> m_pointStyles;
    QList<QPen> m_lineStyles;

    QCPTextElement* m_plotTitle;

    void solveViscousBurgersEquation1D();
    void solveInviscidBurgersEquation1D();

    void resetPlots();
    void addPlot(CPdeField<double>& u, QString plotName, const bool bLine = false);

public:
    MainWindow(QWidget *parent = 0);

private slots:
    void on_btnCalculate_clicked();

};

