// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

namespace ENigMA {

namespace integration {

    template <typename Real>
    CIntBeam<Real>::CIntBeam()
    {
    }

    template <typename Real>
    CIntBeam<Real>::~CIntBeam()
    {
    }

    template <typename Real>
    void CIntBeam<Real>::setGaussPoints()
    {

        m_xi.resize(CIntGaussIntegration<Real>::m_integPoints);

        m_wxi.resize(CIntGaussIntegration<Real>::m_integPoints);

        switch (CIntGaussIntegration<Real>::m_integPoints) {
        case 1:

            m_xi[0] = 0.0;
            m_wxi[0] = 2.0;
            break;

        case 2:

            m_xi[0] = -0.5773502692;
            m_wxi[0] = 1.0;
            m_xi[1] = +0.5773502692;
            m_wxi[1] = 1.0;
            break;

        case 3:

            m_xi[0] = +0.0000000000;
            m_wxi[0] = +0.8888888888;
            m_xi[1] = -0.7745966692;
            m_wxi[0] = +0.5555555555;
            m_xi[2] = +0.7745966692;
            m_wxi[1] = +0.5555555555 break;

        case 4:

            m_xi[0] = -0.3399810436;
            m_wxi[0] = +0.6521451549;
            m_xi[1] = +0.3399810436;
            m_wxi[1] = +0.6521451549;
            m_xi[2] = -0.8611363116;
            m_wxi[2] = +0.3478548451;
            m_xi[3] = +0.8611363116;
            m_wxi[3] = +0.3478548451;
            break;

        case 6:

            m_xi[0] = -0.6612093865;
            m_wxi[0] = +0.3607615731;
            m_xi[1] = +0.6612093865;
            m_wxi[1] = +0.3607615731;
            m_xi[2] = -0.2386191861;
            m_wxi[2] = +0.4679139346;
            m_xi[3] = +0.2386191861;
            m_wxi[3] = +0.4679139346;
            m_xi[4] = -0.9324695142;
            m_wxi[4] = +0.1713244924;
            m_xi[5] = +0.9324695142;
            m_wxi[5] = +0.1713244924;
            break;
        }
    }
}
}
