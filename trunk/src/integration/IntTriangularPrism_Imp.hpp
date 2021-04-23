// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

namespace ENigMA
{
    namespace integration
    {
        template <typename Real>
        CIntTriangularPrism<Real>::CIntTriangularPrism()
        {
        }

        template <typename Real>
        CIntTriangularPrism<Real>::~CIntTriangularPrism()
        {
        }

        template <typename Real>
        void CIntTriangularPrism<Real>::setGaussPoints()
        {
            m_xi.resize(CIntGaussIntegration<Real>::m_integPoints);
            m_eta.resize(CIntGaussIntegration<Real>::m_integPoints);
            m_zeta.resize(CIntGaussIntegration<Real>::m_integPoints);

            m_wxi.resize(CIntGaussIntegration<Real>::m_integPoints);
            m_weta.resize(CIntGaussIntegration<Real>::m_integPoints);
            m_wzeta.resize(CIntGaussIntegration<Real>::m_integPoints);

            switch (CIntGaussIntegration<Real>::m_integPoints)
            {
            case 1:

                m_xi[0] = 0.0;
                m_eta[0] = 0.0;
                m_zeta[0] = 0.0;
                m_wxi[0] = 2.0;
                m_weta[0] = 2.0;
                m_wzeta[0] = 2.0;
                break;

            case 6:

                m_xi[0] = 1.0 / 6.0;
                m_eta[0] = 1.0 / 6.0;
                m_zeta[0] = -1.0 / sqrt(3.0);
                m_wxi[0] = 0.5;
                m_weta[0] = 0.5;
                m_wzeta[0] = 0.5;
                m_xi[1] = 1.0 / 3.0;
                m_eta[1] = 1.0 / 6.0;
                m_zeta[1] = -1.0 / sqrt(3.0);
                m_wxi[1] = 0.5;
                m_weta[1] = 0.5;
                m_wzeta[1] = 0.5;
                m_xi[2] = 1.0 / 6.0;
                m_eta[2] = 1.0 / 3.0;
                m_zeta[2] = -1.0 / sqrt(3.0);
                m_wxi[2] = 0.5;
                m_weta[2] = 0.5;
                m_wzeta[2] = 0.5;
                m_xi[3] = 1.0 / 6.0;
                m_eta[3] = 1.0 / 6.0;
                m_zeta[3] = 1.0 / sqrt(3.0);
                m_wxi[3] = 0.5;
                m_weta[3] = 0.5;
                m_wzeta[3] = 0.5;
                m_xi[4] = 1.0 / 3.0;
                m_eta[4] = 1.0 / 6.0;
                m_zeta[4] = 1.0 / sqrt(3.0);
                m_wxi[4] = 0.5;
                m_weta[4] = 0.5;
                m_wzeta[4] = 0.5;
                m_xi[5] = 1.0 / 6.0;
                m_eta[5] = 1.0 / 3.0;
                m_zeta[5] = 1.0 / sqrt(3.0);
                m_wxi[5] = 0.5;
                m_weta[5] = 0.5;
                m_wzeta[5] = 0.5;
                break;

            case 9:

                m_xi[0] = 1.0 / 6.0;
                m_eta[0] = 1.0 / 6.0;
                m_zeta[0] = -sqrt(3.0) / sqrt(5.0);
                m_wxi[0] = 5.0 / 18.0;
                m_weta[0] = 5.0 / 18.0;
                m_wzeta[0] = 5.0 / 18.0;
                m_xi[1] = 1.0 / 3.0;
                m_eta[1] = 1.0 / 6.0;
                m_zeta[1] = -sqrt(3.0) / sqrt(5.0);
                m_wxi[1] = 5.0 / 18.0;
                m_weta[1] = 5.0 / 18.0;
                m_wzeta[1] = 5.0 / 18.0;
                m_xi[2] = 1.0 / 6.0;
                m_eta[2] = 1.0 / 3.0;
                m_zeta[2] = -sqrt(3.0) / sqrt(5.0);
                m_wxi[2] = 5.0 / 18.0;
                m_weta[2] = 5.0 / 18.0;
                m_wzeta[2] = 5.0 / 18.0;
                m_xi[3] = 1.0 / 6.0;
                m_eta[3] = 1.0 / 6.0;
                m_zeta[3] = 0.0;
                m_wxi[3] = 8.0 / 18.0;
                m_weta[3] = 8.0 / 18.0;
                m_wzeta[3] = 8.0 / 18.0;
                m_xi[4] = 1.0 / 3.0;
                m_eta[4] = 1.0 / 6.0;
                m_zeta[4] = 0.0;
                m_wxi[4] = 8.0 / 18.0;
                m_weta[4] = 8.0 / 18.0;
                m_wzeta[4] = 8.0 / 18.0;
                m_xi[5] = 1.0 / 6.0;
                m_eta[5] = 1.0 / 3.0;
                m_zeta[5] = 0.0;
                m_wxi[5] = 8.0 / 18.0;
                m_weta[5] = 8.0 / 18.0;
                m_wzeta[5] = 8.0 / 18.0;
                m_xi[6] = 1.0 / 6.0;
                m_eta[6] = 1.0 / 6.0;
                m_zeta[6] = +sqrt(3.0) / sqrt(5.0);
                m_wxi[6] = 5.0 / 18.0;
                m_weta[6] = 5.0 / 18.0;
                m_wzeta[6] = 5.0 / 18.0;
                m_xi[7] = 1.0 / 3.0;
                m_eta[7] = 1.0 / 6.0;
                m_zeta[7] = +sqrt(3.0) / sqrt(5.0);
                m_wxi[7] = 5.0 / 18.0;
                m_weta[7] = 5.0 / 18.0;
                m_wzeta[7] = 5.0 / 18.0;
                m_xi[8] = 1.0 / 6.0;
                m_eta[8] = 1.0 / 3.0;
                m_zeta[8] = +sqrt(3.0) / sqrt(5.0);
                m_wxi[8] = 5.0 / 18.0;
                m_weta[8] = 5.0 / 18.0;
                m_wzeta[8] = 5.0 / 18.0;
                break;
            }
        }
    }
}
