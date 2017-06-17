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
        CIntQuadrilateral<Real>::CIntQuadrilateral()
        {
        
        }

        template <typename Real>
        CIntQuadrilateral<Real>::~CIntQuadrilateral()
        {
        
        }

        template <typename Real>
        void CIntQuadrilateral<Real>::setGaussPoints()
        {

            m_xi.resize(CIntGaussIntegration<Real>::m_integPoints);
            m_eta.resize(CIntGaussIntegration<Real>::m_integPoints);

            m_wxi.resize(CIntGaussIntegration<Real>::m_integPoints);
            m_weta.resize(CIntGaussIntegration<Real>::m_integPoints);

            switch(CIntGaussIntegration<Real>::m_integPoints)
            {
            case 1:

                m_xi[0] = 0.0; m_eta[0] = 0.0; m_wxi[0] = 2.0; m_weta[0] = 2.0;
                break;

            case 4:

                m_xi[0] = -0.5773502692; m_eta[0] = -0.5773502692; m_wxi[0] = 1.0; m_weta[0] = 1.0;
                m_xi[1] = +0.5773502692; m_eta[1] = -0.5773502692; m_wxi[1] = 1.0; m_weta[1] = 1.0;
                m_xi[2] = -0.5773502692; m_eta[2] = +0.5773502692; m_wxi[2] = 1.0; m_weta[2] = 1.0;
                m_xi[3] = +0.5773502692; m_eta[3] = +0.5773502692; m_wxi[3] = 1.0; m_weta[3] = 1.0;
                break;

            case 8:

                // TODO: CHECK THIS!
                m_xi[0] = -0.7745966692; m_eta[0] = -0.7745966692; m_wxi[0] = 0.5555555556; m_weta[0] = 0.5555555556;
                m_xi[1] = +0.7745966692; m_eta[1] = -0.7745966692; m_wxi[1] = 0.5555555556; m_weta[1] = 0.5555555556;
                m_xi[2] = -0.7745966692; m_eta[2] = +0.7745966692; m_wxi[2] = 0.5555555556; m_weta[2] = 0.5555555556;
                m_xi[3] = +0.7745966692; m_eta[3] = +0.7745966692; m_wxi[3] = 0.5555555556; m_weta[3] = 0.5555555556;
                m_xi[4] = -0.7745966692; m_eta[4] = +0.0000000000; m_wxi[4] = 0.5555555556; m_weta[4] = 0.8888888889;
                m_xi[5] = +0.7745966692; m_eta[5] = +0.0000000000; m_wxi[5] = 0.5555555556; m_weta[5] = 0.8888888889;
                m_xi[6] = +0.0000000000; m_eta[6] = -0.7745966692; m_wxi[6] = 0.8888888889; m_weta[6] = 0.5555555556;
                m_xi[7] = +0.0000000000; m_eta[7] = -0.7745966692; m_wxi[7] = 0.8888888889; m_weta[7] = 0.5555555556;
                break;

            }

        }

    }
        
}
