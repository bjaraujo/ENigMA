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

#include "../../src/fvm/FvmControlVolume.hpp"

using namespace ENigMA::fvm;

namespace ENigMA
{

    namespace pde
    {

        namespace fvm
        {

            namespace generic
            {
    
                template <typename Real>
                void ddt(CSleSystem<Real>& aSystem, CPdeField<Real>& aField)
                {

                    aSystem.matrixType = MT_SPARSE;

                    aSystem.matrixA.resize(aField.mesh().nbElements(), aField.mesh().nbElements());
                    aSystem.matrixA.reserve(aField.mesh().nbElements());

                    aSystem.vectorB.resize(aField.mesh().nbElements());

                    aSystem.vectorB.setZero();

                    std::vector<bool> isFixed;

                    isFixed.resize(aField.mesh().nbElements(), false);

                    for(typename std::map<Integer, Real>::const_iterator itr = aField.uFixed.begin(); itr != aField.uFixed.end(); ++itr)
                    {

                        Integer anIndex = itr->first;
                        isFixed[anIndex] = true;

                    }

                    for (Integer el = 0; el < aField.mesh().nbElements(); el++)
                    {

                        aSystem.matrixA.coeffRef(el, el) += 1.0;

                        if (isFixed[el])
                            aSystem.vectorB[el] += aField.uFixed[el];
                        else
                            aSystem.vectorB[el] += aField.u[el];

                    }

                    aSystem.matrixA.finalize();

                }

                template <typename Real>
                void laplacian(CSleSystem<Real>& aSystem, CPdeField<Real>& aField)
                {

                    aSystem.matrixType = MT_SPARSE;

                    aSystem.matrixA.resize(aField.mesh().nbElements(), aField.mesh().nbElements());

                    aSystem.matrixA.reserve(aField.mesh().nbElements());

                    aSystem.vectorB.resize(aField.mesh().nbElements());
                    aSystem.vectorB.setZero();

                    for (Integer el = 0; el < aField.mesh().nbElements(); el++)
                    {

                        Integer anElementId = aField.mesh().elementId(el);

                        CMshElement<Real> anElement = aField.mesh().element(anElementId);

                        CFvmControlVolume<Real> aControlVolume;

                        for (Integer i = 0; i < anElement.nbFaceIds(); ++i)
                        {

                            Integer aFaceId = anElement.faceId(i);

                            CFvmFace<Real> aFace;

                            for (Integer j = 0; j < aField.mesh().face(aFaceId).nbNodeIds(); ++j)
                            {

                                Integer aNodeId = aField.mesh().face(aFaceId).nodeId(j);

                                aFace.addNode(CFvmNode<Real>(aField.mesh().node(aNodeId)));

                            }

                            aControlVolume.addFace(aFaceId, aFace);

                        }

                        aControlVolume.calculateCentroid();
                        aControlVolume.calculateVolume();

                        Real Vp = aControlVolume.volume();

                        for (Integer i = 0; i < aControlVolume.nbFaces(); ++i)
                        {

                            Integer aFaceId = aControlVolume.faceId(i);

                            aControlVolume.calculateFaceArea(aFaceId);

                            Real Af = aControlVolume.faceArea(aFaceId);
                            CGeoNormal<Real> nf = aControlVolume.faceNormal(aFaceId);

                            if (aField.mesh().face(aFaceId).hasPair())
                            {

                                Integer aPairFaceId = aField.mesh().face(aFaceId).pairFaceId();
                                Integer aNeighborId = aField.mesh().face(aPairFaceId).elementId();

                                Integer nl = aField.mesh().elementIndex(aNeighborId);

                                Real df = (aField.mesh().elementCentroid(aNeighborId) - aField.mesh().elementCentroid(anElementId)).norm();

                                aSystem.matrixA.coeffRef(el, el) += -Af / df / Vp;
                                aSystem.matrixA.coeffRef(el, nl) += +Af / df / Vp;

                            }
                            else
                            {

                                // Check if it has boundary condition
                                if (aField.faceHasBC(aFaceId))
                                {
                                    
                                    CPdeBoundaryCondition<Real> aCondition = aField.faceBC(aFaceId);

                                    if (aCondition.boundaryConditionType() == BT_GENERIC_FIXED_VALUE)
                                    {

                                        Real df = (aField.mesh().faceCentroid(aFaceId) - aField.mesh().elementCentroid(anElementId)).norm();

                                        aSystem.matrixA.coeffRef(el, el) += -Af / df / Vp;
                                        aSystem.vectorB[el] += -Af / df / Vp * aCondition.conditionValue(CT_GENERIC_FIXED_VALUE);

                                    }

                                }

                            }
        
                        }

                    }

                    aSystem.matrixA.finalize();

                }

                template <typename Real>
                void divergence(CSleSystem<Real>& aSystem, CPdeField<Real>& aField)
                {

                    aSystem.matrixType = MT_SPARSE;

                    aSystem.matrixA.resize(aField.mesh().nbElements(), aField.mesh().nbElements());

                    aSystem.matrixA.reserve(aField.mesh().nbElements());

                    aSystem.vectorB.resize(aField.mesh().nbElements());
                    aSystem.vectorB.setZero();

                    for (Integer el = 0; el < aField.mesh().nbElements(); el++)
                    {

                        Integer anElementId = aField.mesh().elementId(el);

                        CMshElement<Real> anElement = aField.mesh().element(anElementId);

                        CFvmControlVolume<Real> aControlVolume;

                        for (Integer i = 0; i < anElement.nbFaceIds(); ++i)
                        {

                            Integer aFaceId = anElement.faceId(i);

                            CFvmFace<Real> aFace;

                            for (Integer j = 0; j < aField.mesh().face(aFaceId).nbNodeIds(); ++j)
                            {

                                Integer aNodeId = aField.mesh().face(aFaceId).nodeId(j);

                                aFace.addNode(CFvmNode<Real>(aField.mesh().node(aNodeId)));

                            }

                            aControlVolume.addFace(aFaceId, aFace);

                        }

                        aControlVolume.calculateCentroid();
                        aControlVolume.calculateVolume();

                        Real Vp = aControlVolume.volume();

                        for (Integer i = 0; i < aControlVolume.nbFaces(); ++i)
                        {

                            Integer aFaceId = aControlVolume.faceId(i);

                            aControlVolume.calculateFaceArea(aFaceId);

                            Real Af = aControlVolume.faceArea(aFaceId);
                            CGeoNormal<Real> nf = aControlVolume.faceNormal(aFaceId);

                            if (aField.mesh().face(aFaceId).hasPair())
                            {

                                Integer aPairFaceId = aField.mesh().face(aFaceId).pairFaceId();
                                Integer aNeighborId = aField.mesh().face(aPairFaceId).elementId();

                                Integer nl = aField.mesh().elementIndex(aNeighborId);

                                aSystem.matrixA.coeffRef(el, el) += +0.5 * Af / Vp * (nf.x() + nf.y() + nf.z());
                                aSystem.matrixA.coeffRef(el, nl) += +0.5 * Af / Vp * (nf.x() + nf.y() + nf.z());

                            }
                            else
                            {

                                // Check if it has boundary condition
                                if (aField.faceHasBC(aFaceId))
                                {
                                    
                                    CPdeBoundaryCondition<Real> aCondition = aField.faceBC(aFaceId);

                                    if (aCondition.boundaryConditionType() == BT_GENERIC_FIXED_VALUE)
                                    {

                                        aSystem.vectorB[el] += Af / Vp * (nf.x() + nf.y() + nf.z()) * aCondition.conditionValue(CT_GENERIC_FIXED_VALUE);

                                    }

                                }

                            }
        
                        }

                    }

                    aSystem.matrixA.finalize();

                }

                template <typename Real>
                void source(Eigen::Matrix<Real, Eigen::Dynamic, 1>& aVectorB, CPdeField<Real>& aField, Real aSource)
                {

                }

            }

        }

    }
    
}

