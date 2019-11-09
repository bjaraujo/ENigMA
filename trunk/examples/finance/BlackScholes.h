
#pragma once

#include <iostream>

#include <QtCore>

#include "MshBasicMesher.hpp"
#include "PdeEquation.hpp"
#include "PdeBoundaryCondition.hpp"
#include "PosGnuplot.hpp"

using namespace ENigMA::mesh;
using namespace ENigMA::pde;
using namespace ENigMA::post;

class CBlackScholes : public QObject
{
    Q_OBJECT
private:

    double m_K;               // Strike price
    double m_r;               // Interest rate
    double m_sigma;           // Volatility
    double m_T;               // Maturity (years)
    double m_S0;              // Initial price
    double m_Smin;            // Minimum share price
    double m_Smax;            // Maximum share price

    Integer m_N;       // Number of timesteps
    Integer m_M;       // Number of Share prices to look at

    bool m_bDiffTerm;
    bool m_bAdvectTerm;
    bool m_bReactTerm;

    Integer m_plotId;

public:
    CBlackScholes();
    ~CBlackScholes();

    void setInterestRate(double r);
    void setVolatility(double sigma);
    void setMaturity(double T);
    void setStrikePrice(double K);
    void setInitialPrice(double S0);
    void setInterval(double Smin, double Smax);

    void enableDiffTerm(bool bDiffTerm);
    void enableAdvectTerm(bool bAvectTerm);
    void enableReactTerm(bool bReactTerm);

    void setPlotId(Integer plotId);

    void calculate(EDiscretMethod aDiscretMethod, std::string strOutputFileName);

signals:
    void plotChanged(Integer plotId, QVector<double> x, QVector<double> y);

};




