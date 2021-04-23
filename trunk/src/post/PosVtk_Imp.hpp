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
#include <iomanip>

namespace ENigMA
{
    namespace post
    {
        template <typename Real>
        CPosVtk<Real>::CPosVtk()
        {
        }

        template <typename Real>
        CPosVtk<Real>::~CPosVtk()
        {
        }

        template <typename Real>
        bool CPosVtk<Real>::load(CPdeField<Real>& aField, const std::string& strFileName)
        {
            std::ifstream fileVtk;

            fileVtk.open(strFileName.c_str(), std::ios_base::in);

            if (fileVtk.is_open())
            {
                // TODO:

                fileVtk.close();
            }
            else
                return false;

            return true;
        }

        template <typename Real>
        bool CPosVtk<Real>::save(CPdeField<Real>& aField, const std::string& strFileName)
        {
            std::ofstream fileVtk;

            fileVtk.open(strFileName.c_str(), std::ios_base::out | std::ios_base::trunc);

            if (fileVtk.is_open())
            {
                fileVtk << "# vtk DataFile Version 2.0" << std::endl;
                fileVtk << " . " << std::endl;
                fileVtk << "ASCII" << std::endl;
                fileVtk << "DATASET UNSTRUCTURED_GRID" << std::endl;

                fileVtk << "POINTS " << aField.mesh().nbNodes() << " float" << std::endl;

                for (Integer i = 0; i < aField.mesh().nbNodes(); ++i)
                {
                    Integer id = aField.mesh().nodeId(i);

                    fileVtk << aField.mesh().node(id).x() << std::setprecision(16) << " " << aField.mesh().node(id).y() << " " << aField.mesh().node(id).z() << std::endl;
                }

                Integer nSize = 0;

                for (Integer i = 0; i < aField.mesh().nbElements(); ++i)
                {
                    Integer nElemSize = 0;

                    Integer id = aField.mesh().elementId(i);

                    switch (aField.mesh().element(id).elementType())
                    {
                    case ET_BEAM:
                        nElemSize = 2;
                        break;
                    case ET_TRIANGLE:
                        nElemSize = 3;
                        break;
                    case ET_QUADRILATERAL:
                        nElemSize = 4;
                        break;
                    case ET_TETRAHEDRON:
                        nElemSize = 4;
                        break;
                    case ET_HEXAHEDRON:
                        nElemSize = 8;
                        break;
                    case ET_TRIANGULAR_PRISM:
                        nElemSize = 6;
                        break;
                    default:
                        nElemSize = 0;
                        break;
                    }

                    nSize += (nElemSize + 1);
                }

                fileVtk << "CELLS " << aField.mesh().nbElements() << " " << nSize << std::endl;

                for (Integer i = 0; i < aField.mesh().nbElements(); ++i)
                {
                    Integer nElemSize = 0;

                    Integer id = aField.mesh().elementId(i);

                    switch (aField.mesh().element(id).elementType())
                    {
                    case ET_BEAM:
                        nElemSize = 2;
                        fileVtk << nElemSize;
                        fileVtk << " " << aField.mesh().nodeIndex(aField.mesh().element(id).nodeId(0));
                        fileVtk << " " << aField.mesh().nodeIndex(aField.mesh().element(id).nodeId(1));
                        break;
                    case ET_TRIANGLE:
                        nElemSize = 3;
                        fileVtk << nElemSize;
                        fileVtk << " " << aField.mesh().nodeIndex(aField.mesh().element(id).nodeId(0));
                        fileVtk << " " << aField.mesh().nodeIndex(aField.mesh().element(id).nodeId(1));
                        fileVtk << " " << aField.mesh().nodeIndex(aField.mesh().element(id).nodeId(2));
                        break;
                    case ET_QUADRILATERAL:
                        nElemSize = 4;
                        break;
                    case ET_TETRAHEDRON:
                        nElemSize = 4;
                        fileVtk << nElemSize;
                        fileVtk << " " << aField.mesh().nodeIndex(aField.mesh().element(id).nodeId(0));
                        fileVtk << " " << aField.mesh().nodeIndex(aField.mesh().element(id).nodeId(2));
                        fileVtk << " " << aField.mesh().nodeIndex(aField.mesh().element(id).nodeId(1));
                        fileVtk << " " << aField.mesh().nodeIndex(aField.mesh().element(id).nodeId(3));
                        break;
                    case ET_HEXAHEDRON:
                        nElemSize = 8;
                        fileVtk << nElemSize;
                        fileVtk << " " << aField.mesh().nodeIndex(aField.mesh().element(id).nodeId(0));
                        fileVtk << " " << aField.mesh().nodeIndex(aField.mesh().element(id).nodeId(1));
                        fileVtk << " " << aField.mesh().nodeIndex(aField.mesh().element(id).nodeId(2));
                        fileVtk << " " << aField.mesh().nodeIndex(aField.mesh().element(id).nodeId(3));
                        fileVtk << " " << aField.mesh().nodeIndex(aField.mesh().element(id).nodeId(4));
                        fileVtk << " " << aField.mesh().nodeIndex(aField.mesh().element(id).nodeId(5));
                        fileVtk << " " << aField.mesh().nodeIndex(aField.mesh().element(id).nodeId(6));
                        fileVtk << " " << aField.mesh().nodeIndex(aField.mesh().element(id).nodeId(7));
                        break;
                    case ET_TRIANGULAR_PRISM:
                        nElemSize = 6;
                        fileVtk << nElemSize;
                        fileVtk << " " << aField.mesh().nodeIndex(aField.mesh().element(id).nodeId(0));
                        fileVtk << " " << aField.mesh().nodeIndex(aField.mesh().element(id).nodeId(1));
                        fileVtk << " " << aField.mesh().nodeIndex(aField.mesh().element(id).nodeId(2));
                        fileVtk << " " << aField.mesh().nodeIndex(aField.mesh().element(id).nodeId(3));
                        fileVtk << " " << aField.mesh().nodeIndex(aField.mesh().element(id).nodeId(4));
                        fileVtk << " " << aField.mesh().nodeIndex(aField.mesh().element(id).nodeId(5));
                        break;
                    default:
                        nElemSize = 0;
                        break;
                    }

                    fileVtk << std::endl;
                }

                fileVtk << "CELL_TYPES " << aField.mesh().nbElements() << std::endl;

                for (Integer i = 0; i < aField.mesh().nbElements(); ++i)
                {
                    Integer id = aField.mesh().elementId(i);

                    Integer nElemType;

                    switch (aField.mesh().element(id).elementType())
                    {
                    case ET_BEAM:
                        nElemType = 3;
                        break;
                    case ET_TRIANGLE:
                        nElemType = 5;
                        break;
                    case ET_QUADRILATERAL:
                        nElemType = 9;
                        break;
                    case ET_TETRAHEDRON:
                        nElemType = 10;
                        break;
                    case ET_HEXAHEDRON:
                        nElemType = 12;
                        break;
                    case ET_TRIANGULAR_PRISM:
                        nElemType = 13;
                        break;
                    default:
                        nElemType = 0;
                        break;
                    }

                    fileVtk << nElemType << std::endl;
                }

                if (aField.u.size() > 0)
                {
                    if (aField.discretLocation() == DL_NODE)
                    {
                        fileVtk << "POINT_DATA " << aField.mesh().nbNodes() << std::endl;
                        fileVtk << "SCALARS fixed float" << std::endl;
                        fileVtk << "LOOKUP_TABLE default" << std::endl;

                        for (Integer i = 0; i < aField.mesh().nbNodes(); ++i)
                        {
                            Integer index = aField.mesh().nodeIndex(i);

                            if (aField.nbDofs() == 1)
                                fileVtk << aField.u(i) << std::endl;
                        }
                    }
                }

                fileVtk.close();
            }
            else
                return false;

            return true;
        }
    }
}
