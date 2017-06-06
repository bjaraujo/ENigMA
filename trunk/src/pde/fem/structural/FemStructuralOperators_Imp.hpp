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

#include "../../src/fem/structural/FemConstantStrainTriangle.hpp"
#include "../../src/fem/structural/FemConstantStrainTetrahedron.hpp"

using namespace ENigMA::fem;
using namespace ENigMA::fem::structural;

namespace ENigMA
{

    namespace pde
    {

        namespace fem
        {

            namespace structural
            {

                template <typename Real>
                void ddt(CSleSystem<Real>& aSystem, CPdeField<Real>& aField)
                {

                    // TODO:

                }

                template <typename Real>
                void laplacian(CSleSystem<Real>& aSystem, CPdeField<Real>& aField)
                {

                    aSystem.matrixType = MT_SPARSE_SYMMETRIC;

                    aSystem.matrixA.resize(aField.mesh().nbNodes() * aField.nbDofs(), aField.mesh().nbNodes() * aField.nbDofs());

                    aSystem.matrixA.reserve(aField.mesh().nbNodes() * aField.nbDofs());

                    aSystem.vectorB.resize(aField.mesh().nbNodes() * aField.nbDofs());
                    aSystem.vectorB.setZero();

                    for (Integer el = 0; el < aField.mesh().nbElements(); el++)
                    {

                        Integer anElementId = aField.mesh().elementId(el);

                        CMshElement<Real> anElement = aField.mesh().element(anElementId);

                        if (anElement.elementType() == ET_BEAM && anElement.nbNodeIds() == 2)
                        {

                            // TODO:

                        } else if (anElement.elementType() == ET_TRIANGLE && anElement.nbNodeIds() == 3)
                        {

                            CFemConstantStrainTriangle<Real, 3, 2, 1> aTriangle;

                            ENigMA::material::CMatMaterial<Real> aMaterial = aField.material();

                            aTriangle.thickness() = anElement.thickness();

                            aTriangle.elasticModulus() = aMaterial.propertyValue(ENigMA::material::PT_ELASTIC_MODULUS);
                            aTriangle.coeffPoisson() = aMaterial.propertyValue(ENigMA::material::PT_POISSON_COEFFICIENT);

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {

                                CGeoCoordinate<Real> aVertex(aField.mesh().node(anElement.nodeId(i)));
                                aTriangle.addVertex(aVertex);

                            }

                            aTriangle.calculateArea();

                            aTriangle.transient() = false;

                            aTriangle.update();

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {
                                for (Integer k = 0; k < aField.nbDofs(); ++k)
                                {
                                    for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                                    {
                                        for (Integer l = 0; l < aField.nbDofs(); l++)
                                        {

                                            aSystem.matrixA.coeffRef(
                                                aField.mesh().nodeIndex(anElement.nodeId(i)) * aField.nbDofs() + k, 
                                                aField.mesh().nodeIndex(anElement.nodeId(j)) * aField.nbDofs() + l) += aTriangle.laplacian(i * aField.nbDofs() + k, j * aField.nbDofs() + l);
                                        }

                                    }

                                }

                                aSystem.vectorB(aField.mesh().nodeIndex(anElement.nodeId(i))) += aTriangle.source(i);

                            }

                        } else if (anElement.elementType() == ET_QUADRILATERAL && anElement.nbNodeIds() == 4)
                        {

                            // TODO:

                        } else if (anElement.elementType() == ET_TETRAHEDRON && anElement.nbNodeIds() == 4)
                        {

                            CFemConstantStrainTetrahedron<Real, 4, 3, 1> aTetrahedron;

                            ENigMA::material::CMatMaterial<Real> aMaterial = aField.material();

                            aTetrahedron.elasticModulus() = aMaterial.propertyValue(ENigMA::material::PT_ELASTIC_MODULUS);
                            aTetrahedron.coeffPoisson() = aMaterial.propertyValue(ENigMA::material::PT_POISSON_COEFFICIENT);

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {

                                CGeoCoordinate<Real> aVertex(aField.mesh().node(anElement.nodeId(i)));
                                aTetrahedron.addVertex(aVertex);

                            }

                            aTetrahedron.calculateVolume();

                            aTetrahedron.transient() = false;

                            aTetrahedron.update();

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {
                                for (Integer k = 0; k < aField.nbDofs(); ++k)
                                {
                                    for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                                    {
                                        for (Integer l = 0; l < aField.nbDofs(); l++)
                                        {

                                            aSystem.matrixA.coeffRef(
                                                aField.mesh().nodeIndex(anElement.nodeId(i)) * aField.nbDofs() + k, 
                                                aField.mesh().nodeIndex(anElement.nodeId(j)) * aField.nbDofs() + l) += aTetrahedron.laplacian(i * aField.nbDofs() + k, j * aField.nbDofs() + l);
                                        }

                                    }

                                }

                                aSystem.vectorB(aField.mesh().nodeIndex(anElement.nodeId(i))) += aTetrahedron.source(i);

                            }

                        } else if (anElement.elementType() == ET_HEXAHEDRON && anElement.nbNodeIds() == 8)
                        {

                            // TODO:

                        }

                    }

                    aSystem.matrixA.finalize();

                }

                template <typename Real>
                void divergence(CSleSystem<Real>& aSystem, CPdeField<Real>& aField)
                {

                    // TODO:

                }

                template <typename Real>
                void source(Eigen::Matrix<Real, Eigen::Dynamic, 1>& aVectorB, CPdeField<Real>& aField, const Real aSource)
                {

                    // TODO:

                }

            }

        }

    }

}

