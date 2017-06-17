// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "../../src/fem/FemBeam.hpp"
#include "../../src/fem/FemTriangle.hpp"
#include "../../src/fem/FemQuadrilateral.hpp"
#include "../../src/fem/FemTetrahedron.hpp"
#include "../../src/fem/FemTriangularPrism.hpp"
#include "../../src/fem/FemHexahedron.hpp"

using namespace ENigMA::fem;

namespace ENigMA
{

    namespace pde
    {

        namespace fem
        {

            namespace generic
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

                            CFemBeam<Real, 2, 1, 1> aBeam;

                            CGeoCoordinate<Real> p1(aField.mesh().node(anElement.nodeId(0)));
                            CGeoCoordinate<Real> p2(aField.mesh().node(anElement.nodeId(1)));

                            aBeam.setStartPoint(p1);
                            aBeam.setEndPoint(p2);

                            aBeam.calculateLength();

                            aBeam.transient() = true;

                            aBeam.update();

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {
                                for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                                {

                                    aSystem.matrixA.coeffRef(
                                        aField.mesh().nodeIndex(anElement.nodeId(i)), 
                                        aField.mesh().nodeIndex(anElement.nodeId(j))) += aBeam.ddt(i, j);

                                }

                                aSystem.vectorB(anElement.nodeId(i)) += aBeam.source(i);

                            }

                        } else if (anElement.elementType() == ET_TRIANGLE && anElement.nbNodeIds() == 3)
                        {

                            CFemTriangle<Real, 3, 1, 1> aTriangle;

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {

                                CGeoCoordinate<Real> aVertex(aField.mesh().node(anElement.nodeId(i)));
                                aTriangle.addVertex(aVertex);

                            }

                            aTriangle.calculateArea();

                            aTriangle.transient() = true;

                            aTriangle.update();

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {
                                for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                                {

                                    aSystem.matrixA.coeffRef(
                                        aField.mesh().nodeIndex(anElement.nodeId(i)), 
                                        aField.mesh().nodeIndex(anElement.nodeId(j))) += aTriangle.ddt(i, j);

                                }

                                aSystem.vectorB(anElement.nodeId(i)) += aTriangle.source(i);

                            }

                        } else if (anElement.elementType() == ET_QUADRILATERAL && anElement.nbNodeIds() == 4)
                        {

                            CFemQuadrilateral<Real, 4, 1, 1> aQuadrilateral;

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {

                                CGeoCoordinate<Real> aVertex(aField.mesh().node(anElement.nodeId(i)));
                                aQuadrilateral.addVertex(aVertex);

                            }

                            aQuadrilateral.calculateArea();

                            aQuadrilateral.transient() = true;

                            aQuadrilateral.update();

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {
                                for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                                {

                                    aSystem.matrixA.coeffRef(
                                        aField.mesh().nodeIndex(anElement.nodeId(i)), 
                                        aField.mesh().nodeIndex(anElement.nodeId(j))) += aQuadrilateral.ddt(i, j);

                                }

                                aSystem.vectorB(anElement.nodeId(i)) += aQuadrilateral.source(i);

                            }

                        } else if (anElement.elementType() == ET_TETRAHEDRON && anElement.nbNodeIds() == 4)
                        {

                            CFemTetrahedron<Real, 4, 1, 1> aTetrahedron;

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {

                                CGeoCoordinate<Real> aVertex(aField.mesh().node(anElement.nodeId(i)));
                                aTetrahedron.addVertex(aVertex);

                            }

                            aTetrahedron.calculateVolume();

                            aTetrahedron.transient() = true;

                            aTetrahedron.update();

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {
                                for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                                {

                                    aSystem.matrixA.coeffRef(
                                        aField.mesh().nodeIndex(anElement.nodeId(i)), 
                                        aField.mesh().nodeIndex(anElement.nodeId(j))) += aTetrahedron.ddt(i, j);

                                }

                                aSystem.vectorB(anElement.nodeId(i)) += aTetrahedron.source(i);

                            }

                        } else if (anElement.elementType() == ET_HEXAHEDRON && anElement.nbNodeIds() == 8)
                        {

                            CFemHexahedron<Real, 8, 1, 1> aHexahedron;

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {

                                CGeoCoordinate<Real> aVertex(aField.mesh().node(anElement.nodeId(i)));
                                aHexahedron.addVertex(aVertex);

                            }

                            aHexahedron.calculateVolume();

                            aHexahedron.transient() = true;

                            aHexahedron.update();

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {
                                for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                                {

                                    aSystem.matrixA.coeffRef(
                                        aField.mesh().nodeIndex(anElement.nodeId(i)), 
                                        aField.mesh().nodeIndex(anElement.nodeId(j))) += aHexahedron.ddt(i, j);

                                }

                                aSystem.vectorB(anElement.nodeId(i)) += aHexahedron.source(i);

                            }

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

                            CFemBeam<Real, 2, 1, 1> aBeam;

                            CGeoCoordinate<Real> p1(aField.mesh().node(anElement.nodeId(0)));
                            CGeoCoordinate<Real> p2(aField.mesh().node(anElement.nodeId(1)));

                            aBeam.setStartPoint(p1);
                            aBeam.setEndPoint(p2);

                            aBeam.calculateLength();

                            aBeam.transient() = false;

                            aBeam.update();

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
                                                aField.mesh().nodeIndex(anElement.nodeId(j)) * aField.nbDofs() + l) += aBeam.laplacian(i * aField.nbDofs() + k, j * aField.nbDofs() + l);
                                        }

                                    }

                                }

                                aSystem.vectorB(aField.mesh().nodeIndex(anElement.nodeId(i))) += aBeam.source(i);

                            }

                        } else if (anElement.elementType() == ET_TRIANGLE && anElement.nbNodeIds() == 3)
                        {

                            CFemTriangle<Real, 3, 1, 1> aTriangle;

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

                            CFemQuadrilateral<Real, 4, 1, 1> aQuadrilateral;

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {

                                CGeoCoordinate<Real> aVertex(aField.mesh().node(anElement.nodeId(i)));
                                aQuadrilateral.addVertex(aVertex);

                            }

                            aQuadrilateral.calculateArea();

                            aQuadrilateral.transient() = false;

                            aQuadrilateral.update();

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
                                                aField.mesh().nodeIndex(anElement.nodeId(j)) * aField.nbDofs() + l) += aQuadrilateral.laplacian(i * aField.nbDofs() + k, j * aField.nbDofs() + l);
                                        }

                                    }

                                }

                                aSystem.vectorB(aField.mesh().nodeIndex(anElement.nodeId(i))) += aQuadrilateral.source(i);

                            }

                        } else if (anElement.elementType() == ET_TETRAHEDRON && anElement.nbNodeIds() == 4)
                        {

                            CFemTetrahedron<Real, 4, 1, 1> aTetrahedron;

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

                        } 
                        else if (anElement.elementType() == ET_TRIANGULAR_PRISM && anElement.nbNodeIds() == 6)
                        {

                            CFemTriangularPrism<Real, 6, 1, 1> aTriangularPrism;

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {

                                CGeoCoordinate<Real> aVertex(aField.mesh().node(anElement.nodeId(i)));
                                aTriangularPrism.addVertex(aVertex);

                            }

                            aTriangularPrism.calculateVolume();

                            aTriangularPrism.transient() = false;

                            aTriangularPrism.update();

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
                                                aField.mesh().nodeIndex(anElement.nodeId(j)) * aField.nbDofs() + l) += aTriangularPrism.laplacian(i * aField.nbDofs() + k, j * aField.nbDofs() + l);
                                        }

                                    }

                                }

                                aSystem.vectorB(aField.mesh().nodeIndex(anElement.nodeId(i))) += aTriangularPrism.source(i);

                            }

                        }
                        else if (anElement.elementType() == ET_HEXAHEDRON && anElement.nbNodeIds() == 8)
                        {

                            CFemHexahedron<Real, 8, 1, 1> aHexahedron;

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {

                                CGeoCoordinate<Real> aVertex(aField.mesh().node(anElement.nodeId(i)));
                                aHexahedron.addVertex(aVertex);

                            }

                            aHexahedron.calculateVolume();

                            aHexahedron.transient() = false;

                            aHexahedron.update();

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
                                                aField.mesh().nodeIndex(anElement.nodeId(j)) * aField.nbDofs() + l) += aHexahedron.laplacian(i * aField.nbDofs() + k, j * aField.nbDofs() + l);
                                        }

                                    }

                                }

                                aSystem.vectorB(aField.mesh().nodeIndex(anElement.nodeId(i))) += aHexahedron.source(i);

                            }

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

                            CFemBeam<Real, 2, 1, 1> aBeam;

                            CGeoCoordinate<Real> p1(aField.mesh().node(anElement.nodeId(0)));
                            CGeoCoordinate<Real> p2(aField.mesh().node(anElement.nodeId(1)));

                            aBeam.setStartPoint(p1);
                            aBeam.setEndPoint(p2);

                            aBeam.calculateLength();

                            aBeam.transient() = false;

                            aBeam.update();

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
                                                aField.mesh().nodeIndex(anElement.nodeId(j)) * aField.nbDofs() + l) += aBeam.divergence(i * aField.nbDofs() + k, j * aField.nbDofs() + l);
                                        }

                                    }

                                }

                                aSystem.vectorB(aField.mesh().nodeIndex(anElement.nodeId(i))) += aBeam.source(i);

                            }

                        } else if (anElement.elementType() == ET_TRIANGLE && anElement.nbNodeIds() == 3)
                        {

                            // TODO:

                        } else if (anElement.elementType() == ET_QUADRILATERAL && anElement.nbNodeIds() == 4)
                        {

                            // TODO:

                        } else if (anElement.elementType() == ET_TETRAHEDRON && anElement.nbNodeIds() == 4)
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
                void source(Eigen::Matrix<Real, Eigen::Dynamic, 1>& aVectorB, CPdeField<Real>& aField, const Real aSource)
                {

                    aVectorB.resize(aField.mesh().nbNodes());
                    aVectorB.setZero();

                    for (Integer el = 0; el < aField.mesh().nbElements(); el++)
                    {

                        Integer anElementId = aField.mesh().elementId(el);

                        CMshElement<Real> anElement = aField.mesh().element(anElementId);

                        if (anElement.elementType() == ET_BEAM && anElement.nbNodeIds() == 2)
                        {

                            CFemBeam<Real, 2, 1, 1> aBeam;

                            CGeoCoordinate<Real> p1(aField.mesh().node(anElement.nodeId(0)));
                            CGeoCoordinate<Real> p2(aField.mesh().node(anElement.nodeId(1)));

                            aBeam.setStartPoint(p1);
                            aBeam.setEndPoint(p2);

                            aBeam.calculateLength();

                            aBeam.transient() = false;

                            aBeam.update();

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                                aBeam.setSourceOnNode(i, aSource);
        
                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {
                                aVectorB(aField.mesh().nodeIndex(anElement.nodeId(i))) += aBeam.source(i);
                            }

                        } else if (anElement.elementType() == ET_TRIANGLE && anElement.nbNodeIds() == 3)
                        {

                            CFemTriangle<Real, 3, 1, 1> aTriangle;

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {

                                CGeoCoordinate<Real> aVertex(aField.mesh().node(anElement.nodeId(i)));
                                aTriangle.addVertex(aVertex);

                            }

                            aTriangle.calculateArea();

                            aTriangle.transient() = false;

                            aTriangle.update();

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                                aTriangle.setSourceOnNode(i, aSource);
                                
                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {
                                aVectorB(aField.mesh().nodeIndex(anElement.nodeId(i))) += aTriangle.source(i);
                            }

                        } else if (anElement.elementType() == ET_QUADRILATERAL && anElement.nbNodeIds() == 4)
                        {

                            CFemQuadrilateral<Real, 4, 1, 1> aQuadrilateral;

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {

                                CGeoCoordinate<Real> aVertex(aField.mesh().node(anElement.nodeId(i)));
                                aQuadrilateral.addVertex(aVertex);

                            }

                            aQuadrilateral.calculateArea();

                            aQuadrilateral.transient() = false;

                            aQuadrilateral.update();

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                                aQuadrilateral.setSourceOnNode(i, aSource);

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {
                                aVectorB(aField.mesh().nodeIndex(anElement.nodeId(i))) += aQuadrilateral.source(i);
                            }
                            
                        } else if (anElement.elementType() == ET_TETRAHEDRON && anElement.nbNodeIds() == 4)
                        {

                            CFemTetrahedron<Real, 4, 1, 1> aTetrahedron;

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {

                                CGeoCoordinate<Real> aVertex(aField.mesh().node(anElement.nodeId(i)));
                                aTetrahedron.addVertex(aVertex);

                            }

                            aTetrahedron.calculateVolume();

                            aTetrahedron.transient() = false;

                            aTetrahedron.update();

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                                aTetrahedron.setSourceOnNode(i, aSource);

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {
                                aVectorB(aField.mesh().nodeIndex(anElement.nodeId(i))) += aTetrahedron.source(i);
                            }

                        } else if (anElement.elementType() == ET_HEXAHEDRON && anElement.nbNodeIds() == 8)
                        {

                            CFemHexahedron<Real, 8, 1, 1> aHexahedron;

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {

                                CGeoCoordinate<Real> aVertex(aField.mesh().node(anElement.nodeId(i)));
                                aHexahedron.addVertex(aVertex);

                            }

                            aHexahedron.calculateVolume();

                            aHexahedron.transient() = false;

                            aHexahedron.update();

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                                aHexahedron.setSourceOnNode(i, aSource);

                            for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
                            {
                                aVectorB(aField.mesh().nodeIndex(anElement.nodeId(i))) += aHexahedron.source(i);
                            }

                        }

                    }

                }
                    
            }

        }

    }

}

