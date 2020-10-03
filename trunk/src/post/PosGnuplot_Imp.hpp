// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include <iomanip>
#include <fstream>

namespace ENigMA {
namespace post {
    template <typename Real>
    CPosGnuplot<Real>::CPosGnuplot()
    {
    }

    template <typename Real>
    CPosGnuplot<Real>::~CPosGnuplot()
    {
    }

    template <typename Real>
    bool CPosGnuplot<Real>::save(CPdeField<Real>& aField, const std::string& strFileName)
    {
        std::ofstream fileGnuplot;

        fileGnuplot.open(strFileName.c_str(), std::ios_base::out | std::ios_base::trunc);

        if (fileGnuplot.is_open()) {
            fileGnuplot << "# id x y z u" << std::endl;

            if (aField.discretLocation() == DL_NODE) {
                for (Integer i = 0; i < aField.mesh().nbNodes(); ++i) {
                    Integer id = aField.mesh().nodeIndex(i);
                    CMshNode<Real> aNode = aField.mesh().node(id);

                    fileGnuplot << id + 1 << std::setprecision(16)
                                << " " << aNode.x()
                                << " " << aNode.y()
                                << " " << aNode.z()
                                << " " << aField.u(i) << std::endl;
                }
            }

            if (aField.discretLocation() == DL_ELEMENT_CENTER) {
                for (Integer i = 0; i < aField.mesh().nbElements(); ++i) {
                    Integer id = aField.mesh().elementIndex(i);
                    CGeoCoordinate<Real> aCentroid = aField.mesh().elementCentroid(id);

                    fileGnuplot << id + 1 << std::setprecision(16)
                                << " " << aCentroid.x()
                                << " " << aCentroid.y()
                                << " " << aCentroid.z()
                                << " " << aField.u(i)
                                << std::endl;
                }
            }

            fileGnuplot.close();

        } else
            return false;

        return true;
    }
}
}
