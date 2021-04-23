// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "flow/FemFlowTetrahedron.hpp"
#include "flow/FemFlowTriangle.hpp"

using namespace ENigMA::fem;
using namespace ENigMA::fem::flow;

namespace ENigMA
{
    namespace pde
    {
        namespace fem
        {
            namespace flow
            {
                template <typename Real>
                void ddt(CSleSystem<Real>& aSystem, CPdeField<Real>& aField)
                {
                    aSystem.matrixType = MT_SPARSE_SYMMETRIC;

                    aSystem.matrixA.resize(aField.mesh().nbNodes(), aField.mesh().nbNodes());
                    aSystem.matrixA.reserve(aField.mesh().nbNodes());

                    aSystem.vectorB.resize(aField.mesh().nbNodes());

                    aSystem.vectorB.setZero();

                    for (Integer el = 0; el < aField.mesh().nbElements(); el++)
                    {
                        CMshElement<Real> anElement = aField.mesh().element(el);

                        if (anElement.elementType() == ET_BEAM && anElement.nbNodeIds() == 2)
                        {
                            // TODO:
                        }
                        else if (anElement.elementType() == ET_TRIANGLE && anElement.nbNodeIds() == 3)
                        {
                            CFemFlowTriangle<Real, 3, 1, 1> aTriangle;

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {
                                CGeoCoordinate<Real> aVertex(aField.mesh().node(anElement.nodeId(i)));
                                aTriangle.addVertex(aVertex);
                            }

                            aTriangle.calculateArea();

                            aTriangle.calculateTransientTerm();

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {
                                for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                                {
                                    aSystem.matrixA.coeffRef(
                                        aField.mesh().nodeIndex(anElement.nodeId(i)),
                                        aField.mesh().nodeIndex(anElement.nodeId(j)))
                                        += aTriangle.ddt(i, j);
                                }

                                aSystem.vectorB(anElement.nodeId(i)) += aTriangle.source(i);
                            }
                        }
                        else if (anElement.elementType() == ET_QUADRILATERAL && anElement.nbNodeIds() == 4)
                        {
                            // TODO:
                        }
                        else if (anElement.elementType() == ET_TETRAHEDRON && anElement.nbNodeIds() == 4)
                        {
                            CFemFlowTetrahedron<Real, 4, 1, 1> aTetrahedron;

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {
                                CGeoCoordinate<Real> aVertex(aField.mesh().node(anElement.nodeId(i)));
                                aTetrahedron.addVertex(aVertex);
                            }

                            aTetrahedron.calculateVolume();

                            aTetrahedron.calculateTransientTerm();

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {
                                for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                                {
                                    aSystem.matrixA.coeffRef(
                                        aField.mesh().nodeIndex(anElement.nodeId(i)),
                                        aField.mesh().nodeIndex(anElement.nodeId(j)))
                                        += aTetrahedron.ddt(i, j);
                                }

                                aSystem.vectorB(anElement.nodeId(i)) += aTetrahedron.source(i);
                            }
                        }
                        else if (anElement.elementType() == ET_HEXAHEDRON && anElement.nbNodeIds() == 8)
                        {
                            // TODO:
                        }
                    }

                    aSystem.matrixA.finalize();

                    aSystem.vectorB += aSystem.matrixA * aField.u;
                }

                template <typename Real>
                void laplacian(CSleSystem<Real>& aSystem, CPdeField<Real>& aField)
                {
                    aSystem.matrixType = MT_SPARSE_SYMMETRIC;

                    aSystem.matrixA.resize(aField.mesh().nbNodes(), aField.mesh().nbNodes());

                    aSystem.matrixA.reserve(aField.mesh().nbNodes());

                    aSystem.vectorB.resize(aField.mesh().nbNodes());
                    aSystem.vectorB.setZero();

                    for (Integer el = 0; el < aField.mesh().nbElements(); el++)
                    {
                        Integer anElementId = aField.mesh().elementId(el);

                        CMshElement<Real> anElement = aField.mesh().element(anElementId);

                        if (anElement.elementType() == ET_BEAM && anElement.nbNodeIds() == 2)
                        {
                            // TODO:
                        }
                        else if (anElement.elementType() == ET_TRIANGLE && anElement.nbNodeIds() == 3)
                        {
                            CFemFlowTriangle<Real, 3, 1, 1> aTriangle;

                            ENigMA::material::CMatMaterial<Real> aMaterial = aField.material();

                            aTriangle.setThickness(anElement.thickness());

                            aTriangle.setDensity(aMaterial.propertyValue(ENigMA::material::PT_DENSITY));
                            aTriangle.setViscosity(aMaterial.propertyValue(ENigMA::material::PT_DENSITY));

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {
                                CGeoCoordinate<Real> aVertex(aField.mesh().node(anElement.nodeId(i)));
                                aTriangle.addVertex(aVertex);
                            }

                            aTriangle.calculateArea();

                            aTriangle.calculateDiffusiveTerm();

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
                                                aField.mesh().nodeIndex(anElement.nodeId(j)) * aField.nbDofs() + l)
                                                += aTriangle.laplacian(i * aField.nbDofs() + k, j * aField.nbDofs() + l);
                                        }
                                    }
                                }

                                aSystem.vectorB(aField.mesh().nodeIndex(anElement.nodeId(i))) += aTriangle.source(i);
                            }
                        }
                        else if (anElement.elementType() == ET_QUADRILATERAL && anElement.nbNodeIds() == 4)
                        {
                            // TODO:
                        }
                        else if (anElement.elementType() == ET_TETRAHEDRON && anElement.nbNodeIds() == 4)
                        {
                            CFemFlowTetrahedron<Real, 4, 1, 1> aTetrahedron;

                            ENigMA::material::CMatMaterial<Real> aMaterial = aField.material();

                            aTetrahedron.setDensity(aMaterial.propertyValue(ENigMA::material::PT_DENSITY));
                            aTetrahedron.setViscosity(aMaterial.propertyValue(ENigMA::material::PT_DENSITY));

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {
                                CGeoCoordinate<Real> aVertex(aField.mesh().node(anElement.nodeId(i)));
                                aTetrahedron.addVertex(aVertex);
                            }

                            aTetrahedron.calculateVolume();

                            aTetrahedron.calculateDiffusiveTerm();

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
                                                aField.mesh().nodeIndex(anElement.nodeId(j)) * aField.nbDofs() + l)
                                                += aTetrahedron.laplacian(i * aField.nbDofs() + k, j * aField.nbDofs() + l);
                                        }
                                    }
                                }

                                aSystem.vectorB(aField.mesh().nodeIndex(anElement.nodeId(i))) += aTetrahedron.source(i);
                            }
                        }
                        else if (anElement.elementType() == ET_TRIANGULAR_PRISM && anElement.nbNodeIds() == 6)
                        {
                            // TODO:
                        }
                        else if (anElement.elementType() == ET_HEXAHEDRON && anElement.nbNodeIds() == 8)
                        {
                            // TODO:
                        }
                    }

                    aSystem.matrixA.finalize();
                }

                template <typename Real>
                void divergence(CSleSystem<Real>& aSystem, CPdeField<Real>& aField)
                {
                    aSystem.matrixType = MT_SPARSE;

                    aSystem.matrixA.resize(aField.mesh().nbNodes(), aField.mesh().nbNodes());

                    aSystem.matrixA.reserve(aField.mesh().nbNodes());

                    aSystem.vectorB.resize(aField.mesh().nbNodes());
                    aSystem.vectorB.setZero();

                    for (Integer el = 0; el < aField.mesh().nbElements(); el++)
                    {
                        Integer anElementId = aField.mesh().elementId(el);

                        CMshElement<Real> anElement = aField.mesh().element(anElementId);

                        if (anElement.elementType() == ET_BEAM && anElement.nbNodeIds() == 2)
                        {
                            // TODO:
                        }
                        else if (anElement.elementType() == ET_TRIANGLE && anElement.nbNodeIds() == 3)
                        {
                            // TODO:
                        }
                        else if (anElement.elementType() == ET_QUADRILATERAL && anElement.nbNodeIds() == 4)
                        {
                            // TODO:
                        }
                        else if (anElement.elementType() == ET_TETRAHEDRON && anElement.nbNodeIds() == 4)
                        {
                            // TODO:
                        }
                        else if (anElement.elementType() == ET_TRIANGULAR_PRISM && anElement.nbNodeIds() == 6)
                        {
                            // TODO:
                        }
                        else if (anElement.elementType() == ET_HEXAHEDRON && anElement.nbNodeIds() == 8)
                        {
                            // TODO:
                        }
                    }

                    aSystem.matrixA.finalize();
                }

                template <typename Real>
                void divergence(CSleSystem<Real>& aSystem, CPdeField<Real>& aField1, CPdeField<Real>& aField2, Real dt)
                {
                    aSystem.matrixType = MT_SPARSE;

                    aSystem.matrixA.resize(aField1.mesh().nbNodes(), aField1.mesh().nbNodes());

                    aSystem.matrixA.reserve(aField1.mesh().nbNodes());

                    aSystem.vectorB.resize(aField1.mesh().nbNodes());
                    aSystem.vectorB.setZero();

                    for (Integer el = 0; el < aField1.mesh().nbElements(); el++)
                    {
                        Integer anElementId = aField1.mesh().elementId(el);

                        CMshElement<Real> anElement = aField1.mesh().element(anElementId);

                        if (anElement.elementType() == ET_BEAM && anElement.nbNodeIds() == 2)
                        {
                            // TODO:
                        }
                        else if (anElement.elementType() == ET_TRIANGLE && anElement.nbNodeIds() == 3)
                        {
                            CFemFlowTriangle<Real, 3, 1, 1> aTriangle;

                            ENigMA::material::CMatMaterial<Real> aMaterial = aField1.material();

                            aTriangle.setThickness(anElement.thickness());

                            aTriangle.setDt(dt);

                            aTriangle.setDensity(aMaterial.propertyValue(ENigMA::material::PT_DENSITY));
                            aTriangle.setViscosity(aMaterial.propertyValue(ENigMA::material::PT_VISCOSITY));

                            double ue[3], ve[3];

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {
                                Integer aNodeId = anElement.nodeId(i);

                                CGeoCoordinate<Real> aVertex(aField1.mesh().node(aNodeId));
                                aTriangle.addVertex(aVertex);

                                Integer aNodeIndex = aField1.mesh().nodeIndex(aNodeId);

                                ue[i] = aField1.u(aNodeIndex);
                                ve[i] = aField2.u(aNodeIndex);
                            }

                            aTriangle.calculateArea();

                            aTriangle.calculateConvectiveTerm(ue, ve);

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {
                                for (Integer k = 0; k < aField1.nbDofs(); ++k)
                                {
                                    for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                                    {
                                        for (Integer l = 0; l < aField1.nbDofs(); l++)
                                        {
                                            aSystem.matrixA.coeffRef(
                                                aField1.mesh().nodeIndex(anElement.nodeId(i)) * aField1.nbDofs() + k,
                                                aField1.mesh().nodeIndex(anElement.nodeId(j)) * aField1.nbDofs() + l)
                                                += aTriangle.divergence(i * aField1.nbDofs() + k, j * aField1.nbDofs() + l);
                                        }
                                    }
                                }

                                aSystem.vectorB(aField1.mesh().nodeIndex(anElement.nodeId(i))) += aTriangle.source(i);
                            }
                        }
                        else if (anElement.elementType() == ET_QUADRILATERAL && anElement.nbNodeIds() == 4)
                        {
                            // TODO:
                        }
                        else if (anElement.elementType() == ET_TETRAHEDRON && anElement.nbNodeIds() == 4)
                        {
                            // TODO:
                        }
                        else if (anElement.elementType() == ET_TRIANGULAR_PRISM && anElement.nbNodeIds() == 6)
                        {
                            // TODO:
                        }
                        else if (anElement.elementType() == ET_HEXAHEDRON && anElement.nbNodeIds() == 8)
                        {
                            // TODO:
                        }
                    }

                    aSystem.matrixA.finalize();
                }

                template <typename Real>
                void divergence(CSleSystem<Real>& aSystem, CPdeField<Real>& aField1, CPdeField<Real>& aField2, CPdeField<Real>& aField3, Real dt)
                {
                    aSystem.matrixType = MT_SPARSE;

                    aSystem.matrixA.resize(aField1.mesh().nbNodes(), aField1.mesh().nbNodes());

                    aSystem.matrixA.reserve(aField1.mesh().nbNodes());

                    aSystem.vectorB.resize(aField1.mesh().nbNodes());
                    aSystem.vectorB.setZero();

                    for (Integer el = 0; el < aField1.mesh().nbElements(); el++)
                    {
                        Integer anElementId = aField1.mesh().elementId(el);

                        CMshElement<Real> anElement = aField1.mesh().element(anElementId);

                        if (anElement.elementType() == ET_BEAM && anElement.nbNodeIds() == 2)
                        {
                            // TODO:
                        }
                        else if (anElement.elementType() == ET_TRIANGLE && anElement.nbNodeIds() == 3)
                        {
                            // TODO:
                        }
                        else if (anElement.elementType() == ET_QUADRILATERAL && anElement.nbNodeIds() == 4)
                        {
                            // TODO:
                        }
                        else if (anElement.elementType() == ET_TETRAHEDRON && anElement.nbNodeIds() == 4)
                        {
                            CFemFlowTetrahedron<Real, 4, 1, 1> aTetrahedron;

                            ENigMA::material::CMatMaterial<Real> aMaterial = aField1.material();

                            aTetrahedron.setDt(dt);

                            aTetrahedron.setDensity(aMaterial.propertyValue(ENigMA::material::PT_DENSITY));
                            aTetrahedron.setViscosity(aMaterial.propertyValue(ENigMA::material::PT_VISCOSITY));

                            double ue[4], ve[4], we[4];

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {
                                Integer aNodeId = anElement.nodeId(i);

                                CGeoCoordinate<Real> aVertex(aField1.mesh().node(aNodeId));
                                aTetrahedron.addVertex(aVertex);

                                Integer aNodeIndex = aField1.mesh().nodeIndex(aNodeId);

                                ue[i] = aField1.u(aNodeIndex);
                                ve[i] = aField2.u(aNodeIndex);
                                we[i] = aField3.u(aNodeIndex);
                            }

                            aTetrahedron.calculateVolume();

                            aTetrahedron.calculateConvectiveTerm(ue, ve, we);

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {
                                for (Integer k = 0; k < aField1.nbDofs(); ++k)
                                {
                                    for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                                    {
                                        for (Integer l = 0; l < aField1.nbDofs(); l++)
                                        {
                                            aSystem.matrixA.coeffRef(
                                                aField1.mesh().nodeIndex(anElement.nodeId(i)) * aField1.nbDofs() + k,
                                                aField1.mesh().nodeIndex(anElement.nodeId(j)) * aField1.nbDofs() + l)
                                                += aTetrahedron.divergence(i * aField1.nbDofs() + k, j * aField1.nbDofs() + l);
                                        }
                                    }
                                }

                                aSystem.vectorB(aField1.mesh().nodeIndex(anElement.nodeId(i))) += aTetrahedron.source(i);
                            }
                        }
                        else if (anElement.elementType() == ET_TRIANGULAR_PRISM && anElement.nbNodeIds() == 6)
                        {
                            // TODO:
                        }
                        else if (anElement.elementType() == ET_HEXAHEDRON && anElement.nbNodeIds() == 8)
                        {
                            // TODO:
                        }
                    }

                    aSystem.matrixA.finalize();
                }

                template <typename Real>
                void gradient(CSleSystem<Real>& aSystem, CPdeField<Real>& aField, const EComponent aComponent)
                {
                    aSystem.matrixType = MT_SPARSE;

                    aSystem.matrixA.resize(aField.mesh().nbNodes(), aField.mesh().nbNodes());

                    aSystem.matrixA.reserve(aField.mesh().nbNodes());

                    aSystem.vectorB.resize(aField.mesh().nbNodes());
                    aSystem.vectorB.setZero();

                    for (Integer el = 0; el < aField.mesh().nbElements(); el++)
                    {
                        Integer anElementId = aField.mesh().elementId(el);

                        CMshElement<Real> anElement = aField.mesh().element(anElementId);

                        if (anElement.elementType() == ET_BEAM && anElement.nbNodeIds() == 2)
                        {
                            // TODO:
                        }
                        else if (anElement.elementType() == ET_TRIANGLE && anElement.nbNodeIds() == 3)
                        {
                            CFemFlowTriangle<Real, 3, 1, 1> aTriangle;

                            ENigMA::material::CMatMaterial<Real> aMaterial = aField.material();

                            aTriangle.setThickness(anElement.thickness());

                            aTriangle.setDensity(aMaterial.propertyValue(ENigMA::material::PT_DENSITY));
                            aTriangle.setViscosity(aMaterial.propertyValue(ENigMA::material::PT_DENSITY));

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {
                                CGeoCoordinate<Real> aVertex(aField.mesh().node(anElement.nodeId(i)));
                                aTriangle.addVertex(aVertex);
                            }

                            aTriangle.calculateArea();

                            aTriangle.calculateGradientTerm(aComponent);

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
                                                aField.mesh().nodeIndex(anElement.nodeId(j)) * aField.nbDofs() + l)
                                                += aTriangle.gradient(i * aField.nbDofs() + k, j * aField.nbDofs() + l);
                                        }
                                    }
                                }

                                aSystem.vectorB(aField.mesh().nodeIndex(anElement.nodeId(i))) += aTriangle.source(i);
                            }
                        }
                        else if (anElement.elementType() == ET_QUADRILATERAL && anElement.nbNodeIds() == 4)
                        {
                            // TODO:
                        }
                        else if (anElement.elementType() == ET_TETRAHEDRON && anElement.nbNodeIds() == 4)
                        {
                            CFemFlowTetrahedron<Real, 4, 1, 1> aTetrahedron;

                            ENigMA::material::CMatMaterial<Real> aMaterial = aField.material();

                            aTetrahedron.setDensity(aMaterial.propertyValue(ENigMA::material::PT_DENSITY));
                            aTetrahedron.setViscosity(aMaterial.propertyValue(ENigMA::material::PT_DENSITY));

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {
                                CGeoCoordinate<Real> aVertex(aField.mesh().node(anElement.nodeId(i)));
                                aTetrahedron.addVertex(aVertex);
                            }

                            aTetrahedron.calculateVolume();

                            aTetrahedron.calculateGradientTerm(aComponent);

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
                                                aField.mesh().nodeIndex(anElement.nodeId(j)) * aField.nbDofs() + l)
                                                += aTetrahedron.gradient(i * aField.nbDofs() + k, j * aField.nbDofs() + l);
                                        }
                                    }
                                }

                                aSystem.vectorB(aField.mesh().nodeIndex(anElement.nodeId(i))) += aTetrahedron.source(i);
                            }
                        }
                        else if (anElement.elementType() == ET_TRIANGULAR_PRISM && anElement.nbNodeIds() == 6)
                        {
                            // TODO:
                        }
                        else if (anElement.elementType() == ET_HEXAHEDRON && anElement.nbNodeIds() == 8)
                        {
                            // TODO:
                        }
                    }

                    aSystem.matrixA.finalize();
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
