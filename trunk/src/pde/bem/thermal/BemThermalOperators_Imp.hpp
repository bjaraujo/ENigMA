// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// <Copyright> Copyright (c) 2012, All Rights Reserved </Copyright>
// *****************************************************************************

#pragma once

#include "BemTriangle.hpp"
#include "MshElement.hpp"

using namespace ENigMA::bem;
using namespace ENigMA::mesh;

namespace ENigMA
{

    namespace pde
    {

        namespace bem
        {

            namespace thermal
            {

                template <typename Real>
                void ddt(CSleSystem<Real>& aSystem, CPdeField<Real>& aField)
                {

                }

                template <typename Real>
                void laplacian(CSleSystem<Real>& aSystem, CPdeField<Real>& aField)
                {

                    aSystem.matrixType = MT_DENSE;

                    aSystem.matrixA.resize(aField.mesh().nbElements(), aField.mesh().nbElements());
                    aSystem.matrixA.setZero();

                    aSystem.vectorB.resize(aField.mesh().nbElements());
                    aSystem.vectorB.setZero();

                    std::vector<bool> isFixed;

                    isFixed.resize(aField.mesh().nbElements(), false);

                    for(typename std::map<Integer, Real>::const_iterator itr = aField.uFixed.begin(); itr != aField.uFixed.end(); ++itr)
                    {

                        Integer anIndex = itr->first;
                        isFixed[anIndex] = true;

                    }

                    for (Integer eli = 0; eli < aField.mesh().nbElements(); ++eli)
                    {

                        CMshElement<Real>& anElementi = aField.mesh().element(eli);

                        if (anElementi.elementType() != ET_TRIANGLE)
                            continue;

                        CBemTriangle<Real> aTrianglei;

                        for (Integer i = 0; i < anElementi.nbNodeIds(); ++i)
                        {

                            CGeoCoordinate<Real> aVertex(aField.mesh().node(anElementi.nodeId(i)));
                            aTrianglei.addVertex(aVertex);

                        }

                        for (Integer elj = 0; elj < aField.mesh().nbElements(); ++elj)
                        {

                            CMshElement<Real>& anElementj = aField.mesh().element(elj);

                            if (anElementj.elementType() != ET_TRIANGLE)
                                continue;

                            CBemTriangle<Real> aTrianglej;

                            for (Integer j = 0; j < anElementj.nbNodeIds(); ++j)
                            {

                                CGeoCoordinate<Real> aVertex(aField.mesh().node(anElementj.nodeId(j)));
                                aTrianglej.addVertex(aVertex);

                            }

                            Real hij, gij;

                            aTrianglei.laplacianCoeff(eli, elj, aTrianglej, hij, gij);

                            Real Delta = aTrianglej.area();

                            // H
                            if (eli == elj)
                                aSystem.matrixA.coeffRef(eli, elj) += hij - 0.5;
                            else
                                aSystem.matrixA.coeffRef(eli, elj) += hij;

                            if (isFixed[elj])
                            {
                                // H^
                                aSystem.matrixA.coeffRef(eli, elj) += -gij / Delta;

                                aSystem.vectorB[eli] += -gij / Delta * aField.uFixed[elj];
                            }

                        }

                    }

                }

                template <typename Real>
                void divergence(CSleSystem<Real>& aSystem, CPdeField<Real>& aField)
                {

                    // TODO:

                }

                template <typename Real>
                void source(Eigen::Matrix<Real, Eigen::Dynamic, 1>& aVectorB, CPdeField<Real>& aField, Real aSource)
                {

                }

            }

        }

    }

}

