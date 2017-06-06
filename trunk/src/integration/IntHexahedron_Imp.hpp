// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// <Copyright> Copyright (c) 2012, All Rights Reserved </Copyright>
// *****************************************************************************

namespace ENigMA
{

    namespace integration
    {

        template <typename Real>
        CIntHexahedron<Real>::CIntHexahedron()
        {
        
        }

        template <typename Real>
        CIntHexahedron<Real>::~CIntHexahedron()
        {
        
        }


        template <typename Real>
        void CIntHexahedron<Real>::setGaussPoints()
        {

            m_xi.resize(CIntGaussIntegration<Real>::m_integPoints);
            m_eta.resize(CIntGaussIntegration<Real>::m_integPoints);
            m_zeta.resize(CIntGaussIntegration<Real>::m_integPoints);

            m_wxi.resize(CIntGaussIntegration<Real>::m_integPoints);
            m_weta.resize(CIntGaussIntegration<Real>::m_integPoints);
            m_wzeta.resize(CIntGaussIntegration<Real>::m_integPoints);

            switch(CIntGaussIntegration<Real>::m_integPoints)
            {
            case 1:
                    
                m_xi[0] = 0.0; m_eta[0] = 0.0; m_zeta[0] = 0.0; m_wxi[0] = 2.0; m_weta[0] = 2.0; m_wzeta[0] = 2.0;
                break;

            case 8:

                m_xi[0] = -0.5773502692; m_eta[0] = -0.5773502692; m_zeta[0] = -0.5773502692; m_wxi[0] = 1.0; m_weta[0] = 1.0; m_wzeta[0] = 1.0;
                m_xi[1] = +0.5773502692; m_eta[1] = -0.5773502692; m_zeta[1] = -0.5773502692; m_wxi[1] = 1.0; m_weta[1] = 1.0; m_wzeta[1] = 1.0;
                m_xi[2] = -0.5773502692; m_eta[2] = +0.5773502692; m_zeta[2] = -0.5773502692; m_wxi[2] = 1.0; m_weta[2] = 1.0; m_wzeta[2] = 1.0;
                m_xi[3] = +0.5773502692; m_eta[3] = +0.5773502692; m_zeta[3] = -0.5773502692; m_wxi[3] = 1.0; m_weta[3] = 1.0; m_wzeta[3] = 1.0;
                m_xi[4] = -0.5773502692; m_eta[4] = -0.5773502692; m_zeta[4] = +0.5773502692; m_wxi[4] = 1.0; m_weta[4] = 1.0; m_wzeta[4] = 1.0;
                m_xi[5] = +0.5773502692; m_eta[5] = -0.5773502692; m_zeta[5] = +0.5773502692; m_wxi[5] = 1.0; m_weta[5] = 1.0; m_wzeta[5] = 1.0;
                m_xi[6] = -0.5773502692; m_eta[6] = +0.5773502692; m_zeta[6] = +0.5773502692; m_wxi[6] = 1.0; m_weta[6] = 1.0; m_wzeta[6] = 1.0;
                m_xi[7] = +0.5773502692; m_eta[7] = +0.5773502692; m_zeta[7] = +0.5773502692; m_wxi[7] = 1.0; m_weta[7] = 1.0; m_wzeta[7] = 1.0;
                break;

            case 27:

                m_xi[0]  = -0.7745966692; m_eta[0]  = -0.7745966692; m_zeta[0]  = -0.7745966692; m_wxi[0]  = 0.5555555556; m_weta[0]  = 0.5555555556; m_wzeta[0]  = 0.5555555556;
                m_xi[1]  = +0.7745966692; m_eta[1]  = -0.7745966692; m_zeta[1]  = -0.7745966692; m_wxi[1]  = 0.5555555556; m_weta[1]  = 0.5555555556; m_wzeta[1]  = 0.5555555556;
                m_xi[2]  = -0.7745966692; m_eta[2]  = +0.7745966692; m_zeta[2]  = -0.7745966692; m_wxi[2]  = 0.5555555556; m_weta[2]  = 0.5555555556; m_wzeta[2]  = 0.5555555556;
                m_xi[3]  = +0.7745966692; m_eta[3]  = +0.7745966692; m_zeta[3]  = -0.7745966692; m_wxi[3]  = 0.5555555556; m_weta[3]  = 0.5555555556; m_wzeta[3]  = 0.5555555556;
                m_xi[4]  = -0.7745966692; m_eta[4]  = +0.0000000000; m_zeta[4]  = -0.7745966692; m_wxi[4]  = 0.5555555556; m_weta[4]  = 0.8888888889; m_wzeta[4]  = 0.5555555556;
                m_xi[5]  = +0.7745966692; m_eta[5]  = +0.0000000000; m_zeta[5]  = -0.7745966692; m_wxi[5]  = 0.5555555556; m_weta[5]  = 0.8888888889; m_wzeta[5]  = 0.5555555556;
                m_xi[6]  = +0.0000000000; m_eta[6]  = -0.7745966692; m_zeta[6]  = -0.7745966692; m_wxi[6]  = 0.8888888889; m_weta[6]  = 0.5555555556; m_wzeta[6]  = 0.5555555556;
                m_xi[7]  = +0.0000000000; m_eta[7]  = -0.7745966692; m_zeta[7]  = -0.7745966692; m_wxi[7]  = 0.8888888889; m_weta[7]  = 0.5555555556; m_wzeta[7]  = 0.5555555556;
                m_xi[8]  = +0.0000000000; m_eta[8]  = +0.0000000000; m_zeta[8]  = -0.7745966692; m_wxi[8]  = 0.8888888889; m_weta[8]  = 0.8888888889; m_wzeta[8]  = 0.5555555556;
                m_xi[9]  = -0.7745966692; m_eta[9]  = -0.7745966692; m_zeta[9]  = +0.7745966692; m_wxi[9]  = 0.5555555556; m_weta[9]  = 0.5555555556; m_wzeta[9]  = 0.5555555556;
                m_xi[10] = +0.7745966692; m_eta[10] = -0.7745966692; m_zeta[10] = +0.7745966692; m_wxi[10] = 0.5555555556; m_weta[10] = 0.5555555556; m_wzeta[10] = 0.5555555556;
                m_xi[11] = -0.7745966692; m_eta[11] = +0.7745966692; m_zeta[11] = +0.7745966692; m_wxi[11] = 0.5555555556; m_weta[11] = 0.5555555556; m_wzeta[11] = 0.5555555556;
                m_xi[12] = +0.7745966692; m_eta[12] = +0.7745966692; m_zeta[12] = +0.7745966692; m_wxi[12] = 0.5555555556; m_weta[12] = 0.5555555556; m_wzeta[12] = 0.5555555556;
                m_xi[13] = -0.7745966692; m_eta[13] = +0.0000000000; m_zeta[13] = +0.7745966692; m_wxi[13] = 0.5555555556; m_weta[13] = 0.8888888889; m_wzeta[13] = 0.5555555556;
                m_xi[14] = +0.7745966692; m_eta[14] = +0.0000000000; m_zeta[14] = +0.7745966692; m_wxi[14] = 0.5555555556; m_weta[14] = 0.8888888889; m_wzeta[14] = 0.5555555556;
                m_xi[15] = +0.0000000000; m_eta[15] = -0.7745966692; m_zeta[15] = +0.7745966692; m_wxi[15] = 0.8888888889; m_weta[15] = 0.5555555556; m_wzeta[15] = 0.5555555556;
                m_xi[16] = +0.0000000000; m_eta[16] = -0.7745966692; m_zeta[16] = +0.7745966692; m_wxi[16] = 0.8888888889; m_weta[16] = 0.5555555556; m_wzeta[16] = 0.5555555556;
                m_xi[17] = +0.0000000000; m_eta[17] = +0.0000000000; m_zeta[17] = +0.7745966692; m_wxi[17] = 0.8888888889; m_weta[17] = 0.8888888889; m_wzeta[17] = 0.5555555556;
                m_xi[18] = -0.7745966692; m_eta[18] = -0.7745966692; m_zeta[18] = +0.0000000000; m_wxi[18] = 0.5555555556; m_weta[18] = 0.5555555556; m_wzeta[18] = 0.8888888889;
                m_xi[19] = +0.7745966692; m_eta[19] = -0.7745966692; m_zeta[19] = +0.0000000000; m_wxi[19] = 0.5555555556; m_weta[19] = 0.5555555556; m_wzeta[19] = 0.8888888889;
                m_xi[20] = -0.7745966692; m_eta[20] = +0.7745966692; m_zeta[20] = +0.0000000000; m_wxi[20] = 0.5555555556; m_weta[20] = 0.5555555556; m_wzeta[20] = 0.8888888889;
                m_xi[21] = +0.7745966692; m_eta[21] = +0.7745966692; m_zeta[21] = +0.0000000000; m_wxi[21] = 0.5555555556; m_weta[21] = 0.5555555556; m_wzeta[21] = 0.8888888889;
                m_xi[22] = -0.7745966692; m_eta[22] = +0.0000000000; m_zeta[22] = +0.0000000000; m_wxi[22] = 0.5555555556; m_weta[22] = 0.8888888889; m_wzeta[22] = 0.8888888889;
                m_xi[23] = +0.7745966692; m_eta[23] = +0.0000000000; m_zeta[23] = +0.0000000000; m_wxi[23] = 0.5555555556; m_weta[23] = 0.8888888889; m_wzeta[23] = 0.8888888889;
                m_xi[24] = +0.0000000000; m_eta[24] = -0.7745966692; m_zeta[24] = +0.0000000000; m_wxi[24] = 0.8888888889; m_weta[24] = 0.5555555556; m_wzeta[24] = 0.8888888889;
                m_xi[25] = +0.0000000000; m_eta[25] = -0.7745966692; m_zeta[25] = +0.0000000000; m_wxi[25] = 0.8888888889; m_weta[25] = 0.5555555556; m_wzeta[25] = 0.8888888889;
                m_xi[26] = +0.0000000000; m_eta[26] = +0.0000000000; m_zeta[26] = +0.0000000000; m_wxi[26] = 0.8888888889; m_weta[26] = 0.8888888889; m_wzeta[26] = 0.8888888889;
                break;

            }

        }

    }
        
}
