// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

namespace ENigMA {
namespace pde {
    namespace fdm {
        namespace generic {
            template <typename Real>
            void ddt(CSleSystem<Real>& aSystem, CPdeField<Real>& aField)
            {
                aSystem.matrixType = MT_SPARSE;

                aSystem.matrixA.resize(aField.mesh().nbNodes(), aField.mesh().nbNodes());
                aSystem.matrixA.reserve(aField.mesh().nbNodes());
                aSystem.matrixA.setZero();

                aSystem.vectorB.resize(aField.mesh().nbNodes());
                aSystem.vectorB.setZero();

                std::vector<bool> isFixed;

                isFixed.resize(aField.mesh().nbNodes(), false);

                for (typename std::map<Integer, Real>::const_iterator itr = aField.uFixed.begin(); itr != aField.uFixed.end(); ++itr) {
                    Integer anIndex = itr->first;
                    isFixed[anIndex] = true;
                }

                for (Integer ii = 0; ii < aField.mesh().nbNodes(); ++ii) {
                    aSystem.matrixA.coeffRef(ii, ii) += 1.0;

                    if (isFixed[ii])
                        aSystem.vectorB[ii] += aField.uFixed[ii];
                    else
                        aSystem.vectorB[ii] += aField.u[ii];
                }

                aSystem.matrixA.finalize();
            }

            template <typename Real>
            void laplacian(CSleSystem<Real>& aSystem, CPdeField<Real>& aField)
            {
                aSystem.matrixType = MT_SPARSE;

                aSystem.matrixA.resize(aField.mesh().nbNodes(), aField.mesh().nbNodes());
                aSystem.matrixA.reserve(aField.mesh().nbNodes());
                aSystem.matrixA.setZero();

                aSystem.vectorB.resize(aField.mesh().nbNodes());
                aSystem.vectorB.setZero();

                std::vector<bool> isFixed;

                isFixed.resize(aField.mesh().nbNodes(), false);

                for (typename std::map<Integer, Real>::const_iterator itr = aField.uFixed.begin(); itr != aField.uFixed.end(); ++itr) {
                    Integer anIndex = itr->first;
                    isFixed[anIndex] = true;
                }

                // Assume constant spacing
                Real dx = aField.mesh().dx();

                for (Integer ii = 0; ii < aField.mesh().nbNodes(); ++ii) {
                    if (ii == 0) {
                        Integer ip;

                        ip = ii + 1;

                        if (isFixed[ii]) {
                            aSystem.matrixA.coeffRef(ii, ii) += -2.0 / (dx * dx);
                            aSystem.matrixA.coeffRef(ii, ip) += +1.0 / (dx * dx);

                            aSystem.vectorB[ii] += -1.0 / (dx * dx) * aField.uFixed[ii];
                        } else {
                            // Fixed gradient
                            aSystem.matrixA.coeffRef(ii, ii) += -1.0 / (dx * dx);
                            aSystem.matrixA.coeffRef(ii, ip) += +1.0 / (dx * dx);
                        }

                    } else if (ii < aField.mesh().nbNodes() - 1) {
                        Integer im, ip;

                        im = ii - 1;
                        ip = ii + 1;

                        aSystem.matrixA.coeffRef(ii, im) += +1.0 / (dx * dx);
                        aSystem.matrixA.coeffRef(ii, ii) += -2.0 / (dx * dx);
                        aSystem.matrixA.coeffRef(ii, ip) += +1.0 / (dx * dx);

                    } else if (ii == aField.mesh().nbNodes() - 1) {
                        Integer im;

                        im = ii - 1;

                        if (isFixed[ii]) {
                            aSystem.matrixA.coeffRef(ii, im) += +1.0 / (dx * dx);
                            aSystem.matrixA.coeffRef(ii, ii) += -2.0 / (dx * dx);

                            aSystem.vectorB[ii] += -1.0 / (dx * dx) * aField.uFixed[ii];
                        } else {
                            // Fixed gradient
                            aSystem.matrixA.coeffRef(ii, im) += +1.0 / (dx * dx);
                            aSystem.matrixA.coeffRef(ii, ii) += -1.0 / (dx * dx);
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
                aSystem.matrixA.setZero();

                aSystem.vectorB.resize(aField.mesh().nbNodes());
                aSystem.vectorB.setZero();

                std::vector<bool> isFixed;

                isFixed.resize(aField.mesh().nbNodes(), false);

                for (typename std::map<Integer, Real>::const_iterator itr = aField.uFixed.begin(); itr != aField.uFixed.end(); ++itr) {
                    Integer anIndex = itr->first;
                    isFixed[anIndex] = true;
                }

                // Assume constant spacing
                Real dx = aField.mesh().dx();

                for (Integer ii = 0; ii < aField.mesh().nbNodes(); ++ii) {
                    if (ii == 0) {
                        Integer ip;

                        ip = ii + 1;

                        aSystem.matrixA.coeffRef(ii, ii) = 1.0;

                        if (isFixed[ii])
                            aSystem.vectorB[ii] += -(aField.u(ip) - aField.uFixed[ii]) / dx;

                    } else if (ii < aField.mesh().nbNodes() - 1) {
                        Integer im, ip;

                        im = ii - 1;
                        ip = ii + 1;

                        aSystem.matrixA.coeffRef(ii, im) += -1.0 / (2.0 * dx);
                        aSystem.matrixA.coeffRef(ii, ip) += +1.0 / (2.0 * dx);

                    } else if (ii == aField.mesh().nbNodes() - 1) {
                        Integer im;

                        im = ii - 1;

                        aSystem.matrixA.coeffRef(ii, ii) = 1.0;

                        if (isFixed[ii])
                            aSystem.vectorB[ii] += -(aField.uFixed[ii] - aField.u(im)) / dx;
                    }
                }

                aSystem.matrixA.finalize();
            }

            template <typename Real>
            void source(Eigen::Matrix<Real, Eigen::Dynamic, 1>& aVectorB, CPdeField<Real>& aField, Real aSource)
            {
                // TODO:
            }
        }
    }
}
}
