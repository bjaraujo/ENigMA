// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#include "BlackScholes.h"

CBlackScholes::CBlackScholes()
{

    m_N = 100;
    m_M = 50;

}

CBlackScholes::~CBlackScholes()
{


}

void CBlackScholes::setInterestRate(double r)
{

    m_r = r;

}

void CBlackScholes::setVolatility(double sigma)
{

    m_sigma = sigma;
}


void CBlackScholes::setMaturity(double T)
{

    m_T = T;
}

void CBlackScholes::setStrikePrice(double K)
{

    m_K = K;

}

void CBlackScholes::setInitialPrice(double S0)
{

    m_S0 = S0;

}

void CBlackScholes::setInterval(double Smin, double Smax)
{

    m_Smin = Smin;
    m_Smax = Smax;

}

void CBlackScholes::setPlotId(Integer plotId)
{

    m_plotId = plotId;

}

void CBlackScholes::enableDiffTerm(bool bDiffTerm)
{

    m_bDiffTerm = bDiffTerm;

}

void CBlackScholes::enableAdvectTerm(bool bAdvectTerm)
{

    m_bAdvectTerm = bAdvectTerm;

}

void CBlackScholes::enableReactTerm(bool bReactTerm)
{

    m_bReactTerm = bReactTerm;

}

void CBlackScholes::calculate(EDiscretMethod aDiscretMethod, std::string strOutputFileName)
{

    CMshMesh<double> aMesh;
    CMshBasicMesher<double> aBasicMesher;

    const Integer nu = m_M;

    if (aDiscretMethod == DM_FDM || aDiscretMethod == DM_FEM)
    {

       CGeoCoordinate<double> aPoint1(m_Smin, 0.0, 0.0);
       CGeoCoordinate<double> aPoint2(m_Smax, 0.0, 0.0);

       CGeoLine<double> aLine(aPoint1, aPoint2);

       aBasicMesher.generate(aLine, nu);
       aMesh = aBasicMesher.mesh();

    }
    else if (aDiscretMethod == DM_FVM)
    {

       CGeoCoordinate<double> aVertex1(+m_Smin, -0.05, -0.05);
       CGeoCoordinate<double> aVertex2(+m_Smax, -0.05, -0.05);
       CGeoCoordinate<double> aVertex3(+m_Smax, +0.05, -0.05);
       CGeoCoordinate<double> aVertex4(+m_Smin, +0.05, -0.05);
       CGeoCoordinate<double> aVertex5(+m_Smin, -0.05, +0.05);
       CGeoCoordinate<double> aVertex6(+m_Smax, -0.05, +0.05);
       CGeoCoordinate<double> aVertex7(+m_Smax, +0.05, +0.05);
       CGeoCoordinate<double> aVertex8(+m_Smin, +0.05, +0.05);

       CGeoHexahedron<double> aBox;

       aBox.addVertex(aVertex1);
       aBox.addVertex(aVertex2);
       aBox.addVertex(aVertex3);
       aBox.addVertex(aVertex4);
       aBox.addVertex(aVertex5);
       aBox.addVertex(aVertex6);
       aBox.addVertex(aVertex7);
       aBox.addVertex(aVertex8);

       aBasicMesher.generate(aBox, nu, 1, 1);

       aMesh = aBasicMesher.mesh();

       aMesh.generateFaces(1E-12);

       aMesh.calculateFaceCentroid();
       aMesh.calculateElementCentroid();

    }

    CPdeField<double> V;
    CPdeField<double> S, Ssqr;

    V.setMesh(aMesh);
    V.setDiscretMethod(aDiscretMethod);
    V.setDiscretOrder(DO_LINEAR);

    if (aDiscretMethod == DM_FDM || aDiscretMethod == DM_FEM)
       V.setDiscretLocation(DL_NODE);
    else if (aDiscretMethod == DM_FVM)
       V.setDiscretLocation(DL_ELEMENT_CENTER);

    V.setSimulationType(ST_GENERIC);
    V.setNbDofs(1);

    int wi = -1;

    // Set BC and initial conditions
    if (aDiscretMethod == DM_FDM || aDiscretMethod == DM_FEM)
    {
       V.u.resize(V.mesh().nbNodes());

       S.u.resize(V.mesh().nbNodes());
       S.setDiscretMethod(DM_NONE);

       Ssqr.u.resize(V.mesh().nbNodes());
       Ssqr.setDiscretMethod(DM_NONE);

       for (Integer i = 0; i < V.mesh().nbNodes(); ++i)
       {

          Integer aNodeId = V.mesh().nodeId(i);
          double x = V.mesh().node(aNodeId).x();

          S.u(i) = x;
          Ssqr.u(i) = x * x;

          V.u(i) = std::max(m_K - x, 0.0); //  Value of option at maturity

          if (x > m_S0 && wi == -1)
             wi = i;

       }

    }
    else if (aDiscretMethod == DM_FVM)
    {

       V.u.resize(V.mesh().nbElements());

       S.u.resize(V.mesh().nbElements());
       S.setDiscretMethod(DM_NONE);

       Ssqr.u.resize(V.mesh().nbElements());
       Ssqr.setDiscretMethod(DM_NONE);

       for (Integer i = 0; i < V.mesh().nbElements(); ++i)
       {

          Integer anElementId = V.mesh().elementId(i);
          double x = V.mesh().elementCentroid(anElementId).x();

          S.u(i) = x;
          Ssqr.u(i) = x * x;

          V.u(i) = std::max(m_K - x, 0.0); //  Value of option at maturity

          if (x > m_S0 && wi == -1)
             wi = i;

       }

    }

    double dt = m_T / m_N;

    // Black-scholes equation
    //CPdeEquation<double> aPdeEquation(-1.0 / dt * ddt<double>(V) + 0.5 * sigma * sigma * (Ssqr.u * laplacian<double>(V)) + r * (S.u * divergence<double>(V)) - r * V = 0);

    // Faster to store systems in memory and then reuse them for each time step

    double fd;
    if (m_bDiffTerm) fd = 1.0; else fd = 0.0;

    double fa;
    if (m_bAdvectTerm) fa = 1.0; else fa = 0.0;

    double fr;
    if (m_bReactTerm) fr = 1.0; else fr = 0.0;

    CSleSystem<double> DiffusionTerm = +fd * 0.5 * m_sigma * m_sigma * (Ssqr.u * laplacian<double>(V));
    CSleSystem<double> AdvectionTerm = +fa * m_r * (S.u * divergence<double>(V));
    CSleSystem<double> ReactionTerm  = -fr * m_r * V;

    for (Integer i = 0; i < m_N; ++i)
    {

       //std::cout << "Time = " << dt * (i + 1) << std::endl;

       // Transient term is inverted: it is backwards in time
       CPdeEquation<double> aPdeEquation(-1.0 / dt * ddt<double>(V) + DiffusionTerm + AdvectionTerm + ReactionTerm = 0);

       aPdeEquation.solve(V);

       QVector<double> x, y;

       x.resize(S.u.size());
       y.resize(V.u.size());

       for (int i = 0; i < S.u.size(); ++i)
           x[i] = S.u[i];

       for (int i = 0; i < V.u.size(); ++i)
           y[i] = V.u[i];

       emit plotChanged(m_plotId, x, y);

    }

    if (wi >= 0)
    {
       // Linear interpolation
       double P = V.u(wi-1) + (V.u(wi) - V.u(wi-1)) * (m_S0 - S.u(wi-1)) / (S.u(wi) - S.u(wi-1));
       std::cout << "Option price = " << P << std::endl;
    }

}

