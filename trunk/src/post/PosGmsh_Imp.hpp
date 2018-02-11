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

namespace ENigMA
{

    namespace post
    {

        template <typename Real>
        CPosGmsh<Real>::CPosGmsh()
        {

        }

        template <typename Real>
        CPosGmsh<Real>::~CPosGmsh()
        {

        }

        template <typename Real>
        bool CPosGmsh<Real>::findSection(std::ifstream& fileStream, std::string sectionName)
        {

            while (!fileStream.eof())
            {

                std::string line;

                std::getline(fileStream, line);

                if (line.find(sectionName) != std::string::npos)
                {
                    return true;
                }

            }

            return false;

        }

        template <typename Real>
        bool CPosGmsh<Real>::load(ENigMA::pde::CPdeField<Real>& aField, const std::string strFileName)
        {

            std::ifstream fileGmsh;

            fileGmsh.open(strFileName.c_str(), std::ios_base::in);

            if (fileGmsh.is_open())
            {

                if (findSection(fileGmsh, "$MeshFormat"))
                {
                    findSection(fileGmsh, "$EndMeshFormat");
                }

                aField.mesh().reset();

                if (findSection(fileGmsh, "$Nodes"))
                {

                    Integer nbNodes;

                    fileGmsh >> nbNodes;

                    for (Integer i = 0; i < nbNodes; ++i)
                    {

                        Integer id;
                        Real x, y, z;

                        fileGmsh >> id >> x >> y >> z;

                        ENigMA::mesh::CMshNode<Real> aNode(x, y, z);

                        aField.mesh().addNode(id - 1, aNode);

                    }

                    //findSection(fileGmsh, "$EndNodes");

                }

                if (findSection(fileGmsh, "$Elements"))
                {

                    Integer nbElements;

                    fileGmsh >> nbElements;

                    for (Integer i = 0; i < nbElements; ++i)
                    {

                        Integer id, nType, nTags, nAux;

                        Integer nbNodeIds = 0;

                        fileGmsh >> id >> nType >> nTags;

                        if (nType < 1 || nType > 6)
                        {
                            std::string line;
                            std::getline(fileGmsh, line);
                            continue;
                        }

                        // TODO: Use tags
                        for (Integer j = 0; j < nTags; ++j)
                        {
                            fileGmsh >> nAux;
                        }

                        ENigMA::mesh::CMshElement<Real> anElement;

                        switch (nType)
                        {
                        case 1:
                            // line
                            anElement.setElementType(ENigMA::mesh::ET_BEAM);
                            nbNodeIds = 2;
                            break;
                        case 2:
                            // triangle
                            anElement.setElementType(ENigMA::mesh::ET_TRIANGLE);
                            nbNodeIds = 3;
                            break;
                        case 3:
                            // quad
                            anElement.setElementType(ENigMA::mesh::ET_QUADRILATERAL);
                            nbNodeIds = 4;
                            break;
                        case 4:
                            // tetra
                            anElement.setElementType(ENigMA::mesh::ET_TETRAHEDRON);
                            nbNodeIds = 4;
                            break;
                        case 5:
                            // hexa
                            anElement.setElementType(ENigMA::mesh::ET_HEXAHEDRON);
                            nbNodeIds = 8;
                            break;
                        case 6:
                            // triangular prism
                            anElement.setElementType(ENigMA::mesh::ET_TRIANGULAR_PRISM);
                            nbNodeIds = 6;
                            break;
                        }

                        Integer aNodeId;

                        for (Integer j = 0; j < nbNodeIds; ++j)
                        {

                            fileGmsh >> aNodeId;

                            anElement.addNodeId(aNodeId - 1);

                        }

                        if (anElement.elementType() == ENigMA::mesh::ET_TETRAHEDRON)
                            anElement.invert();

                        aField.mesh().addElement(id - 1, anElement);

                    }

                    //findSection(fileGmsh, "$EndElements");

                }

                fileGmsh.close();

            }
            else
                return false;

            return true;

        }

        template <typename Real>
        bool CPosGmsh<Real>::save(ENigMA::pde::CPdeField<Real>& aField, const std::string strFileName, const std::string strViewName)
        {

            std::ofstream fileGmsh;

            fileGmsh.open(strFileName.c_str(), std::ios_base::out | std::ios_base::trunc);

            if (fileGmsh.is_open())
            {

                fileGmsh << "$MeshFormat" << std::endl;
                fileGmsh << "2.2 0 8" << std::endl;
                fileGmsh << "$EndMeshFormat" << std::endl;

                fileGmsh << "$Nodes" << std::endl;
                fileGmsh << aField.mesh().nbNodes() << std::endl;

                for (Integer i = 0; i < aField.mesh().nbNodes(); ++i)
                {

                    Integer id = aField.mesh().nodeId(i);
                    Integer index = aField.mesh().nodeIndex(id);

                    fileGmsh << index + 1 << std::setprecision(16) << " " << aField.mesh().node(id).x() << " " << aField.mesh().node(id).y() << " " << aField.mesh().node(id).z() << std::endl;

                }

                fileGmsh << "$EndNodes" << std::endl;

                fileGmsh << "$Elements" << std::endl;
                fileGmsh << aField.mesh().nbElements() << std::endl;

                for (Integer i = 0; i < aField.mesh().nbElements(); ++i)
                {

                    Integer id = aField.mesh().elementId(i);
                    Integer index = aField.mesh().elementIndex(id);

                    ENigMA::mesh::CMshElement<Real> anElement = aField.mesh().element(id);

                    if (anElement.elementType() == ENigMA::mesh::ET_TETRAHEDRON)
                        anElement.invert();

                    Integer nElemType;

                    switch (anElement.elementType())
                    {
                    case ENigMA::mesh::ET_BEAM:
                        nElemType = 1;
                        break;
                    case ENigMA::mesh::ET_TRIANGLE:
                        nElemType = 2;
                        break;
                    case ENigMA::mesh::ET_QUADRILATERAL:
                        nElemType = 3;
                        break;
                    case ENigMA::mesh::ET_TETRAHEDRON:
                        nElemType = 4;
                        break;
                    case ENigMA::mesh::ET_HEXAHEDRON:
                        nElemType = 5;
                        break;
                    case ENigMA::mesh::ET_TRIANGULAR_PRISM:
                        nElemType = 6;
                        break;
                    default:
                        nElemType = 0;
                        break;
                    }

                    fileGmsh << index + 1 << " " << nElemType << " 0";

                    for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                    {
                        Integer aNodeId = anElement.nodeId(j);
                        fileGmsh << " " << aField.mesh().nodeIndex(aNodeId) + 1;
                    }

                    fileGmsh << std::endl;

                }

                fileGmsh << "$EndElements" << std::endl;

                if (aField.u.size() > 0)
                {

                    if (aField.discretLocation() == ENigMA::pde::DL_NODE)
                    {

                        fileGmsh << "$NodeData" << std::endl;
                        fileGmsh << "1" << std::endl;
                        fileGmsh << """" << strViewName << """" << std::endl;

                        fileGmsh << "1" << std::endl;
                        fileGmsh << "0.0" << std::endl;

                        fileGmsh << "3" << std::endl;
                        fileGmsh << "0" << std::endl;

                        fileGmsh << aField.nbDofs() << std::endl;

                        fileGmsh << aField.mesh().nbNodes() << std::endl;

                        for (Integer i = 0; i < aField.mesh().nbNodes(); ++i)
                        {

                            Integer index = aField.mesh().nodeIndex(i);

                            if (aField.nbDofs() == 1)
                                fileGmsh << index + 1 << std::setprecision(16) << " " << aField.u(i) << std::endl;
                            else if (aField.nbDofs() == 2)
                                fileGmsh << index + 1 << std::setprecision(16) << " " << aField.u(i * 2 + 0) << " " << aField.u(i * 2 + 1) << std::endl;
                            else if (aField.nbDofs() == 3)
                                fileGmsh << index + 1 << std::setprecision(16) << " " << aField.u(i * 3 + 0) << " " << aField.u(i * 3 + 1) << " " << aField.u(i * 3 + 2) << std::endl;

                        }

                        fileGmsh << "$EndNodeData" << std::endl;

                    }
                    else if (aField.discretLocation() == ENigMA::pde::DL_ELEMENT_CENTER)
                    {

                        fileGmsh << "$ElementData" << std::endl;
                        fileGmsh << "1" << std::endl;
                        fileGmsh << """" << strViewName << """" << std::endl;

                        fileGmsh << "1" << std::endl;
                        fileGmsh << "0.0" << std::endl;

                        fileGmsh << "3" << std::endl;
                        fileGmsh << "0" << std::endl;

                        fileGmsh << aField.nbDofs() << std::endl;

                        fileGmsh << aField.mesh().nbElements() << std::endl;

                        for (Integer i = 0; i < aField.mesh().nbElements(); ++i)
                        {

                            Integer index = aField.mesh().elementIndex(i);

                            if (aField.nbDofs() == 1)
                                fileGmsh << index + 1 << std::setprecision(16) << " " << aField.u(i) << std::endl;
                            else if (aField.nbDofs() == 2)
                                fileGmsh << index + 1 << std::setprecision(16) << " " << aField.u(i * 2 + 0) << " " << aField.u(i * 2 + 1) << std::endl;
                            else if (aField.nbDofs() == 3)
                                fileGmsh << index + 1 << std::setprecision(16) << " " << aField.u(i * 3 + 0) << " " << aField.u(i * 3 + 1) << " " << aField.u(i * 3 + 2) << std::endl;

                        }

                        fileGmsh << "$EndElementData" << std::endl;

                    }

                }

                fileGmsh.close();

            }
            else
                return false;

            return true;

        }

    }

}

