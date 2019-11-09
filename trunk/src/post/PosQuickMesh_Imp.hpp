// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include <fstream>

namespace ENigMA {

namespace post {

    template <typename Real>
    CPosQuickMesh<Real>::CPosQuickMesh()
    {
    }

    template <typename Real>
    CPosQuickMesh<Real>::~CPosQuickMesh()
    {
    }

    template <typename Real>
    bool CPosQuickMesh<Real>::load(CPdeField<Real>& aField, std::string strFileName)
    {

        std::ifstream fileQuickMesh;

        fileQuickMesh.open(strFileName.c_str(), std::ios_base::in);

        if (fileQuickMesh.is_open()) {

            while (!fileQuickMesh.eof()) {

                std::string line;

                std::getline(fileQuickMesh, line);

                if (line.find("C") != std::string::npos) {

                    std::string desc;
                    Integer id;
                    Real x, y, z;

                    std::sscanf(line, "%s %d %f %f %f", &desc, &id, &x, &y, &z);

                    CMshNode<Real> aNode(x, y, z);

                    aField.mesh().addNode(id - 1, aNode);

                } else if (line.find("t") != std::string::npos) {

                    std::string desc;
                    Integer id;
                    Integer n1, n2, n3, n4;
                    Integer p1, p2, p3, p4;
                    Real v;

                    std::sscanf(line, "%s %d %d %d %d %d %d %d %d %d %f", &desc, &id, &n1, &n2, &n3, &n4, &p1, &p2, &p3, &p4, &v);

                    CMshElement<Real> anElement;

                    anElement.setElementType(ET_TETRAHEDRON);

                    anElement.addNodeId(n1 - 1);
                    anElement.addNodeId(n2 - 1);
                    anElement.addNodeId(n3 - 1);
                    anElement.addNodeId(n4 - 1);

                    aField.mesh().addElement(id - 1, anElement);
                }
            }

            fileQuickMesh.close();

        } else
            return false;

        return true;
    }

    template <typename Real>
    bool CPosQuickMesh<Real>::save(CPdeField<Real>& aField, std::string strFileName)
    {

        std::ofstream fileQuickMesh;

        fileQuickMesh.open(strFileName.c_str(), std::ios_base::out | std::ios_base::trunc);

        if (fileQuickMesh.is_open()) {

            fileQuickMesh << "P " << aField.mesh().nbNodes() << std::endl;
            fileQuickMesh << "T " << aField.mesh().nbElements() << std::endl;

            for (Integer i = 0; i < aField.mesh().nbNodes(); ++i) {

                Integer id = aField.mesh().nodeId(i);
                Integer index = aField.mesh().nodeIndex(id);

                fileQuickMesh << "C " << index + 1 << std::setprecision(16) << " " << aField.mesh().node(id).x() << " " << aField.mesh().node(id).y() << " " << aField.mesh().node(id).z() << std::endl;
            }

            for (Integer i = 0; i < aField.mesh().nbElements(); ++i) {

                Integer id = aField.mesh().elementId(i);
                Integer index = aField.mesh().elementIndex(id);

                fileQuickMesh << "t " << index + 1;

                for (Integer j = 0; j < aField.mesh().element(id).nbNodeIds(); ++j) {
                    Integer aNodeId = aField.mesh().element(id).nodeId(j);
                    fileQuickMesh << " " << aField.mesh().nodeIndex(aNodeId) + 1;
                }

                fileQuickMesh << " 1" << std::endl;
            }

            fileQuickMesh.close();

        } else
            return false;

        return true;
    }
}
}
