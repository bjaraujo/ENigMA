// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

using namespace ENigMA::geometry;

namespace ENigMA
{

    namespace bem
    {

        template <typename Real>
        CBemTriangle<Real>::CBemTriangle()
        {

            CIntGaussIntegration<Real>::m_integPoints = 20;

            CBemTriangle<Real>::setBarycentricGaussPoints();

        }

        template <typename Real>
        CBemTriangle<Real>::~CBemTriangle()
        {

        }

        template <typename Real>
        bool CBemTriangle<Real>::getIntegrationPoints(std::vector<CGeoCoordinate<Real> >& sPoints, std::vector<Real>& sWeights)
        {

            CGeoCoordinate<Real> n1 = this->m_vertices[0];
            CGeoCoordinate<Real> n2 = this->m_vertices[1];
            CGeoCoordinate<Real> n3 = this->m_vertices[2];

            for (Integer k = 0; k < CIntGaussIntegration<Real>::m_integPoints; ++k)
            {

                CGeoCoordinate<Real> aPoint;

                aPoint.x() = CIntTriangle<Real>::m_beta1[k] * n1.x() + CIntTriangle<Real>::m_beta2[k] * n2.x() + CIntTriangle<Real>::m_beta3[k] * n3.x();
                aPoint.y() = CIntTriangle<Real>::m_beta1[k] * n1.y() + CIntTriangle<Real>::m_beta2[k] * n2.y() + CIntTriangle<Real>::m_beta3[k] * n3.y();
                aPoint.z() = CIntTriangle<Real>::m_beta1[k] * n1.z() + CIntTriangle<Real>::m_beta2[k] * n2.z() + CIntTriangle<Real>::m_beta3[k] * n3.z();

                sPoints.push_back(aPoint);
                sWeights.push_back(CIntTriangle<Real>::m_wbeta[k]);

            }

            return true;

        }

        template <typename Real>
        void CBemTriangle<Real>::laplacianCoeff(const Integer i, const Integer j, CBemTriangle<Real>& aBemTriangle, Real& h, Real& g)
        {

            CGeoTriangle<Real>::calculateCentroid();
            CGeoCoordinate<Real> ci = CGeoTriangle<Real>::centroid();

            Real xi = ci.x();
            Real yi = ci.y();
            Real zi = ci.z();

            aBemTriangle.calculateCentroid();
            CGeoCoordinate<Real> cj = aBemTriangle.centroid();

            Real h_sum = 0.0;
            Real g_sum = 0.0;

            aBemTriangle.calculateArea();
            Real area = aBemTriangle.area();

            if (i == j)
            {

                CGeoCoordinate<Real> n1 = aBemTriangle.vertex(0);
                CGeoCoordinate<Real> n2 = aBemTriangle.vertex(1);
                CGeoCoordinate<Real> n3 = aBemTriangle.vertex(2);

                CGeoVector<Real> r1, r2, r3;
                CGeoVector<Real> r21, r32, r13;

                r21 = n2 - n1;
                r32 = n3 - n2;
                r13 = n1 - n3;

                r1 = cj - n1;
                r2 = cj - n2;
                r3 = cj - n3;

                Real teta1 = r2.angle(r3);
                Real teta2 = r3.angle(r1);
                Real teta3 = r1.angle(r2);

                Real alpha1 = r1.angle(r21);
                Real alpha2 = r2.angle(r32);
                Real alpha3 = r3.angle(r13);

                h_sum = 0.0;
                g_sum = 2.0 / 3.0 * area * (1.0 / r32.norm() * log(tan((teta1 + alpha2) / 2.0)/tan(alpha2 / 2.0)) + 
                                            1.0 / r13.norm() * log(tan((teta2 + alpha3) / 2.0)/tan(alpha3 / 2.0)) +
                                            1.0 / r21.norm() * log(tan((teta3 + alpha1) / 2.0)/tan(alpha1 / 2.0)));

            }
            else
            {

                std::vector<CGeoCoordinate<Real> > sPoints;
                std::vector<Real> sWeights;

                aBemTriangle.getIntegrationPoints(sPoints, sWeights);

                aBemTriangle.calculateNormal();

                Real nx = aBemTriangle.normal().x();
                Real ny = aBemTriangle.normal().y();
                Real nz = aBemTriangle.normal().z();

                h_sum = 0.0f;
                g_sum = 0.0f;

                for (Integer k = 0; k < CIntGaussIntegration<Real>::m_integPoints; ++k)
                {

                    Real xk = sPoints[k].x();
                    Real yk = sPoints[k].y();
                    Real zk = sPoints[k].z();

                    Real wk = sWeights[k];

                    h_sum += 2.0 * area * wk * ((xk-xi)*nx+(yk-yi)*ny+(zk-zi)*nz) / pow((xk-xi)*(xk-xi)+(yk-yi)*(yk-yi)+(zk-zi)*(zk-zi), 1.5);
                    g_sum += 2.0 * area * wk / sqrt((xk-xi)*(xk-xi)+(yk-yi)*(yk-yi)+(zk-zi)*(zk-zi));

                }

            }

            const Real pi = std::acos(-1.0);

            h = -1.0 / (4.0 * pi) * h_sum;
            g = +1.0 / (4.0 * pi) * g_sum;

            return;

        }

    }

}
