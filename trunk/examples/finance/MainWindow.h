
#pragma once

#include "ui_MainWindow.h"

#include "BlackScholes.h"

class MainWindow : public QMainWindow, private Ui::MainWindow
{
    Q_OBJECT

private:

    CBlackScholes m_blackScholes;

public:
    MainWindow(QWidget *parent = 0);

private slots:
    void on_btnCalculate_clicked();

public slots:
    void updatePlot(Integer plotId, QVector<double> x, QVector<double> y);

};

