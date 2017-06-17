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

    namespace fem
    {

        namespace thermal
        {

            template <typename Real>
            CFemLinearTemperatureTriangle<Real, 3, 1, 1>::CFemLinearTemperatureTriangle()
            {
                
            }

            template <typename Real>
            CFemLinearTemperatureTriangle<Real, 3, 1, 1>::~CFemLinearTemperatureTriangle()
            {

            }

            template <typename Real>
            void CFemLinearTemperatureTriangle<Real, 3, 1, 1>::setDiffusionTerm()
            {

                CFemTriangle<Real, 3, 1, 1>::setDiffusionTerm();

                CFemElement<Real>::laplacian *= CFemThermalElement<Real>::thermalConductivity();

            }

            template <typename Real>
            void CFemLinearTemperatureTriangle<Real, 3, 1, 1>::setConvectionOnEdge(const Integer anEdgeIndex, const Real h, const Real Tinf)
            {

                Integer aNodeIndex1 = (anEdgeIndex + 0) % 3; 
                Integer aNodeIndex2 = (anEdgeIndex + 1) % 3; 

                CGeoCoordinate<Real> n1 = CGeoVertexList<Real>::vertex(aNodeIndex1);
                CGeoCoordinate<Real> n2 = CGeoVertexList<Real>::vertex(aNodeIndex2);

                CGeoVector<Real> v = n2 - n1;

                Real edgeLenth = v.norm();

                CFemElement<Real>::source(aNodeIndex1) += h * Tinf * edgeLenth * this->m_thickness * 0.5;
                CFemElement<Real>::source(aNodeIndex2) += h * Tinf * edgeLenth * this->m_thickness * 0.5;

                CFemElement<Real>::laplacian(aNodeIndex1, aNodeIndex1) += h * edgeLenth * this->m_thickness / 3.0;
                CFemElement<Real>::laplacian(aNodeIndex1, aNodeIndex2) += h * edgeLenth * this->m_thickness / 6.0;

                CFemElement<Real>::laplacian(aNodeIndex2, aNodeIndex1) += h * edgeLenth * this->m_thickness / 6.0;
                CFemElement<Real>::laplacian(aNodeIndex2, aNodeIndex2) += h * edgeLenth * this->m_thickness / 3.0;

            }

            template <typename Real>
            void CFemLinearTemperatureTriangle<Real, 3, 1, 1>::setConvectionOnEdge(const Integer anEdgeIndex, const Real e, const Real teta, const Real Tinf)
            {

                // Stefan-boltzmann constant
                Real sigma = 5.6704E-8;

                CFemThermalElement<Real>::setSourceOnEdge(anEdgeIndex, sigma * e * teta, Tinf);

            }

        }

    }

}
