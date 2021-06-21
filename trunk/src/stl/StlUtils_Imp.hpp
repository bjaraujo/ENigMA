// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#pragma warning(disable : 4996)

#include <fstream>

#include "GeoHashGrid.hpp"

using namespace ENigMA::geometry;

namespace ENigMA
{
    namespace stl
    {
        template <typename Real>
        CStlUtils<Real>::CStlUtils()
        {
            m_stlFile.stats.type = FT_ASCII;
        }

        template <typename Real>
        CStlUtils<Real>::CStlUtils(ENigMA::mesh::CMshMesh<Real>& aTriMesh)
        {
            m_stlFile.stats.type = FT_ASCII;

            this->set(aTriMesh);
        }

        template <typename Real>
        CStlUtils<Real>::~CStlUtils()
        {
        }

        template <typename Real>
        ENigMA::stl::CStlFile<Real>& CStlUtils<Real>::stlFile()
        {
            return m_stlFile;
        }

        template <typename Real>
        void CStlUtils<Real>::set(ENigMA::mesh::CMshMesh<Real>& aTriMesh)
        {
            m_stlFile.reset();

            for (Integer i = 0; i < aTriMesh.nbElements(); ++i)
            {
                Integer anElementId = aTriMesh.elementId(i);

                ENigMA::mesh::CMshElement<Real> anElement = aTriMesh.element(anElementId);

                if (anElement.elementType() != ENigMA::mesh::ET_TRIANGLE)
                    continue;

                CStlFacet<Real> aFacet;

                for (Integer j = 0; j < anElement.nbNodeIds(); ++j)
                {
                    Integer aNodeId = anElement.nodeId(j);

                    CGeoCoordinate<Real> aVertex = aTriMesh.node(aNodeId);

                    aFacet.addVertex(aVertex);
                }

                Integer aFacetId = anElementId;

                m_stlFile.addFacet(aFacetId, aFacet);
            }
        }

        template <typename Real>
        bool CStlUtils<Real>::load(const std::string& strFileName)
        {
            // Open the file
            m_stlFile.fileName = strFileName;

            std::ifstream fileSTL;

            fileSTL.open(strFileName.c_str(), std::ios_base::in | std::ios::binary);

            std::string line;

            if (fileSTL.is_open())
            {
                // Find size of file
                fileSTL.seekg(0, std::ios::end);

                m_stlFile.stats.fileSize = fileSTL.tellg();

                // Check for binary or ASCII file
                fileSTL.seekg(HEADER_SIZE);

                char chtest[256];
                fileSTL.read(chtest, sizeof(chtest));

                if (fileSTL)
                {
                    m_stlFile.stats.type = FT_ASCII;

                    char* pch;
                    pch = strstr(chtest, "vertex");

                    if (pch == NULL)
                        m_stlFile.stats.type = FT_BINARY;

                    fileSTL.seekg(0);

                    Integer numFacets = 0;

                    // Get the header and the number of facets in the .STL file
                    // If the .STL file is binary, then do the following
                    if (m_stlFile.stats.type == FT_BINARY)
                    {
                        // Test if the STL file has the right size
                        if (((m_stlFile.stats.fileSize - HEADER_SIZE) % SIZEOF_STL_FACET != 0)
                            || (m_stlFile.stats.fileSize < STL_MIN_FILE_SIZE))
                            return false;

                        numFacets = static_cast<Integer>((m_stlFile.stats.fileSize - HEADER_SIZE) / SIZEOF_STL_FACET);

                        // Read the header
                        fileSTL.read(m_stlFile.stats.header, LABEL_SIZE);

                        if (fileSTL)
                        {
                            m_stlFile.stats.header[80] = '\0';

                            // Read the int following the header.  This should contain # of facets
                            //Integer headerNumFacets = this->getInt(fileSTL);
                            //std::cout << headerNumFacets << std::endl;
                        }
                    }
                    else
                    {
                        // Otherwise, if the .STL file is ASCII, then do the following

                        // Get the header
                        std::getline(fileSTL, line);

                        char fhtest[5];

                        // Discard spaces
                        Integer c = 0;

                        while (c < 5)
                        {
                            char letter = fileSTL.get();

                            if (letter != 32)
                            {
                                fhtest[c] = letter;
                                c++;
                            }
                        }

#ifdef WIN32
                        if (strnicmp(fhtest, "facet", sizeof(fhtest)) != 0)
                            return false;
#else
                        if (strncasecmp(fhtest, "facet", sizeof(fhtest)) != 0)
                            return false;
#endif

                        fileSTL.seekg(0);

                        // Find the number of facets
                        Integer numLines = 1;

                        for (Integer i = 0, j = 0; i < m_stlFile.stats.fileSize; i++, ++j)
                        {
                            if (fileSTL.get() == '\n')
                            {
                                if (j > 4) // don't count short lines
                                    numLines++;

                                j = 0;
                            }
                        }

                        fileSTL.seekg(0);

                        // Read the header
                        for (Integer i = 0; (i < 80) && (m_stlFile.stats.header[i] = static_cast<char>(fileSTL.get())) != '\n'; ++i)
                            ;

                        m_stlFile.stats.header[80] = '\0';

                        numFacets = numLines / ASCII_LINES_PER_FACET;
                    }

                    if (m_stlFile.stats.type == FT_BINARY)
                    {
                        fileSTL.seekg(HEADER_SIZE);
                    }
                    else
                    {
                        fileSTL.seekg(0);
                        // Skip the first line of the file
                        std::getline(fileSTL, line);
                    }

                    m_stlFile.reset();

                    for (Integer i = 0; i < numFacets; ++i)
                    {
                        CStlFacet<Real> aFacet;

                        if (m_stlFile.stats.type == FT_BINARY)
                        {
                            // Read a single facet from a binary .STL file
                            aFacet.normal().x() = this->getFloat(fileSTL);
                            aFacet.normal().y() = this->getFloat(fileSTL);
                            aFacet.normal().z() = this->getFloat(fileSTL);

                            for (Integer j = 0; j < 3; ++j)
                            {
                                CGeoCoordinate<Real> aVertex;

                                aVertex.x() = this->getFloat(fileSTL);
                                aVertex.y() = this->getFloat(fileSTL);
                                aVertex.z() = this->getFloat(fileSTL);

                                aFacet.addVertex(aVertex);
                            }

                            aFacet.extra[0] = static_cast<char>(fileSTL.get());
                            aFacet.extra[1] = static_cast<char>(fileSTL.get());
                        }
                        else
                        {
                            std::string strFacet, strNormal;
                            std::string strOuter, strLoop;
                            std::string strVertex;
                            std::string strEndLoop;
                            std::string strEndFacet;

                            // Read a single facet from an ASCII .STL file
                            fileSTL >> strFacet >> strNormal >> aFacet.normal().x() >> aFacet.normal().y() >> aFacet.normal().z();

                            fileSTL >> strOuter >> strLoop;

                            for (Integer j = 0; j < 3; ++j)
                            {
                                CGeoCoordinate<Real> aVertex;

                                fileSTL >> strVertex >> aVertex.x() >> aVertex.y() >> aVertex.z();

                                aFacet.addVertex(aVertex);
                            }

                            fileSTL >> strEndLoop;
                            fileSTL >> strEndFacet;
                        }

                        // Add the facet.
                        m_stlFile.addFacet(i, aFacet);
                    }
                }

                fileSTL.close();

                this->calculateStatistics();

                return true;
            }
            else
                return false;
        }

        template <typename Real>
        bool CStlUtils<Real>::saveAscii(const std::string& strFileName)
        {
            std::ofstream fileSTL;

            fileSTL.open(strFileName.c_str(), std::ios_base::out);

            if (fileSTL.is_open())
            {
                fileSTL << "solid" << std::endl;

                fileSTL.precision(8);
                fileSTL.setf(std::ios::scientific);

                for (Integer i = 0; i < m_stlFile.nbFacets(); ++i)
                {
                    Integer aFacetId = m_stlFile.facetId(i);

                    fileSTL << "facet normal "
                            << m_stlFile.facet(aFacetId).normal().x() << " "
                            << m_stlFile.facet(aFacetId).normal().y() << " "
                            << m_stlFile.facet(aFacetId).normal().z() << std::endl;

                    fileSTL << "outer loop" << std::endl;

                    for (Integer j = 0; j < 3; ++j)
                    {
                        fileSTL << "vertex "
                                << m_stlFile.facet(aFacetId).vertex(j).x() << " "
                                << m_stlFile.facet(aFacetId).vertex(j).y() << " "
                                << m_stlFile.facet(aFacetId).vertex(j).z() << std::endl;
                    }

                    fileSTL << "endloop" << std::endl;
                    fileSTL << "endfacet" << std::endl;
                }

                fileSTL << "endsolid" << std::endl;

                fileSTL.close();
            }

            return true;
        }

        template <typename Real>
        bool CStlUtils<Real>::writeFloat(std::ofstream& stream, const Real aValue)
        {
            float fvalue = static_cast<float>(aValue);
            stream.write(reinterpret_cast<const char*>(&fvalue), sizeof(float));
            return true;
        }

        template <typename Real>
        bool CStlUtils<Real>::writeInt(std::ofstream& stream, const Integer value)
        {
            stream.write(reinterpret_cast<const char*>(&value), sizeof(int));
            return true;
        }

        template <typename Real>
        Real CStlUtils<Real>::getFloat(std::ifstream& stream)
        {
            float fvalue;
            stream.read(reinterpret_cast<char*>(&fvalue), sizeof(float));
            return static_cast<Real>(fvalue);
        }

        template <typename Real>
        Integer CStlUtils<Real>::getInt(std::ifstream& stream)
        {
            Integer ivalue;
            stream.read(reinterpret_cast<char*>(&ivalue), sizeof(int));
            return ivalue;
        }

        template <typename Real>
        bool CStlUtils<Real>::saveBinary(const std::string& strFileName)
        {
            std::ofstream fileSTL;

            fileSTL.open(strFileName.c_str(), std::ios_base::out | std::ios::binary);

            if (fileSTL.is_open())
            {
                for (Integer i = 0; i < LABEL_SIZE; ++i)
                    fileSTL.put(0);

                this->writeInt(fileSTL, m_stlFile.nbFacets());

                for (Integer i = 0; i < m_stlFile.nbFacets(); ++i)
                {
                    Integer aFacetId = m_stlFile.facetId(i);

                    this->writeFloat(fileSTL, m_stlFile.facet(aFacetId).normal().x());
                    this->writeFloat(fileSTL, m_stlFile.facet(aFacetId).normal().y());
                    this->writeFloat(fileSTL, m_stlFile.facet(aFacetId).normal().z());

                    for (Integer j = 0; j < 3; ++j)
                    {
                        this->writeFloat(fileSTL, m_stlFile.facet(aFacetId).vertex(j).x());
                        this->writeFloat(fileSTL, m_stlFile.facet(aFacetId).vertex(j).y());
                        this->writeFloat(fileSTL, m_stlFile.facet(aFacetId).vertex(j).z());
                    }

                    fileSTL.put(m_stlFile.facet(aFacetId).extra[0]);
                    fileSTL.put(m_stlFile.facet(aFacetId).extra[1]);
                }

                fileSTL.close();
            }

            return true;
        }

        template <typename Real>
        bool CStlUtils<Real>::save(const std::string& strFileName)
        {
            if (m_stlFile.stats.type == FT_ASCII)
                return this->saveAscii(strFileName);
            else
                return this->saveBinary(strFileName);
        }

        template <typename Real>
        bool CStlUtils<Real>::add(CStlUtils<Real>& aStl)
        {
            for (Integer i = 0; i < aStl.stlFile().nbFacets(); ++i)
            {
                Integer aFacetId = aStl.stlFile().facetId(i);
                Integer aNewFacetId = m_stlFile.nextFacetId();

                this->addFacet(aNewFacetId, aStl.stlFile().facet(aFacetId));
            }

            this->calculateStatistics();

            return true;
        }

        template <typename Real>
        bool CStlUtils<Real>::calculateStatistics()
        {
            // Initialize the max and min values the first time through
            if (m_stlFile.nbFacets() > 0)
            {
                Integer aFacetId = m_stlFile.facetId(0);

                CStlFacet<Real> aFacet = m_stlFile.facet(aFacetId);

                m_stlFile.stats.max = aFacet.vertex(0);
                m_stlFile.stats.min = aFacet.vertex(0);
            }

            for (Integer i = 0; i < m_stlFile.nbFacets(); ++i)
            {
                Integer aFacetId = m_stlFile.facetId(i);

                CStlFacet<Real> aFacet = m_stlFile.facet(aFacetId);

                // while we are going through all of the facets, let's find the
                // maximum and minimum values for x, y, and z

                // Now find the max and min values
                for (Integer j = 0; j < 3; ++j)
                {
                    m_stlFile.stats.max.x() = std::max(m_stlFile.stats.max.x(), aFacet.vertex(j).x());
                    m_stlFile.stats.min.x() = std::min(m_stlFile.stats.min.x(), aFacet.vertex(j).x());
                    m_stlFile.stats.max.y() = std::max(m_stlFile.stats.max.y(), aFacet.vertex(j).y());
                    m_stlFile.stats.min.y() = std::min(m_stlFile.stats.min.y(), aFacet.vertex(j).y());
                    m_stlFile.stats.max.z() = std::max(m_stlFile.stats.max.z(), aFacet.vertex(j).z());
                    m_stlFile.stats.min.z() = std::min(m_stlFile.stats.min.z(), aFacet.vertex(j).z());
                }
            }

            m_stlFile.stats.size = m_stlFile.stats.max - m_stlFile.stats.min;

            return true;
        }

        template <typename Real>
        bool CStlUtils<Real>::addFacet(const Integer aFacetId, ENigMA::stl::CStlFacet<Real> aFacet)
        {
            m_stlFile.addFacet(aFacetId, aFacet);
            return true;
        }

        template <typename Real>
        bool CStlUtils<Real>::removeFacet(const Integer aFacetId)
        {
            m_stlFile.removeFacet(aFacetId);
            return true;
        }

        template <typename Real>
        bool CStlUtils<Real>::removeDuplicateFacets(const Real aTolerance)
        {
            CGeoHashGrid<Real> aHashGrid;

            for (Integer i = 0; i < m_stlFile.nbFacets(); ++i)
            {
                Integer aFacetId = m_stlFile.facetId(i);

                m_stlFile.facet(aFacetId).calculateCentroid();
                CGeoCoordinate<Real> aCenterCoordinate = m_stlFile.facet(aFacetId).centroid();

                aHashGrid.addGeometricObject(aFacetId, aCenterCoordinate);
            }

            aHashGrid.build();

            for (Integer i = 0; i < m_stlFile.nbFacets(); ++i)
            {
                Integer aFacetId = m_stlFile.facetId(i);

                m_stlFile.facet(aFacetId).calculateArea();

                std::vector<Integer> sCoordinates;

                CGeoCoordinate<Real> aCenterCoordinate = m_stlFile.facet(aFacetId).centroid();

                aHashGrid.find(sCoordinates, aCenterCoordinate, aTolerance);

                for (Integer k = 1; k < static_cast<Integer>(sCoordinates.size()); ++k)
                {
                    Integer wFacetId = sCoordinates.at(k);

                    m_stlFile.facet(wFacetId).calculateArea();

                    Real area_diff = std::fabs(m_stlFile.facet(aFacetId).area() - m_stlFile.facet(wFacetId).area());

                    if (area_diff < aTolerance * aTolerance)
                    {
                        // Remove this facet
                        this->removeFacet(wFacetId);
                        i--;
                    }
                }
            }

            return true;
        }

        template <typename Real>
        bool CStlUtils<Real>::removeInvalidFacets(const Real aTolerance)
        {
            for (Integer i = 0; i < m_stlFile.nbFacets(); ++i)
            {
                Integer aFacetId = m_stlFile.facetId(i);

                Real d1 = (m_stlFile.facet(aFacetId).vertex(1) - m_stlFile.facet(aFacetId).vertex(0)).norm();
                Real d2 = (m_stlFile.facet(aFacetId).vertex(2) - m_stlFile.facet(aFacetId).vertex(1)).norm();
                Real d3 = (m_stlFile.facet(aFacetId).vertex(0) - m_stlFile.facet(aFacetId).vertex(2)).norm();

                bool naked1 = m_stlFile.facet(aFacetId).edge(0).naked();
                bool naked2 = m_stlFile.facet(aFacetId).edge(1).naked();
                bool naked3 = m_stlFile.facet(aFacetId).edge(2).naked();

                // Remove invalid facets
                if (d1 < aTolerance || d2 < aTolerance || d3 < aTolerance || (naked1 && naked2 && naked3))
                {
                    this->removeFacet(aFacetId);
                    i--;
                }
            }

            return true;
        }

        template <typename Real>
        bool CStlUtils<Real>::generateConnectivity(const Real aTolerance)
        {
            // Discover double edges
            CGeoHashGrid<Real> aHashGrid;

            std::vector<CGeoCoordinate<Real>> sCenterCoordinates;

            Integer anEdgeId;

            for (Integer i = 0; i < m_stlFile.nbFacets(); ++i)
            {
                Integer aFacetId = m_stlFile.facetId(i);

                for (Integer j = 0; j < 3; ++j)
                {
                    CGeoCoordinate<Real> aCenterCoordinate = (m_stlFile.facet(aFacetId).vertex((j + 0) % 3) + m_stlFile.facet(aFacetId).vertex((j + 1) % 3)) * 0.5;

                    sCenterCoordinates.push_back(aCenterCoordinate);

                    anEdgeId = aFacetId * 3 + j;

                    aHashGrid.addGeometricObject(anEdgeId, aCenterCoordinate);
                }
            }

            aHashGrid.build();

            for (Integer i = 0; i < m_stlFile.nbFacets(); ++i)
            {
                Integer aFacetId = m_stlFile.facetId(i);

                // initialize neighbors list to -1 to mark unconnected edges
                for (Integer j = 0; j < 3; ++j)
                {
                    m_stlFile.facet(aFacetId).edge(j).setNeighbor(-1);
                    m_stlFile.facet(aFacetId).edge(j).setWhichVertexNot(-1);
                }
            }

            for (Integer i = 0; i < m_stlFile.nbFacets(); ++i)
            {
                Integer aFacetId = m_stlFile.facetId(i);

                for (Integer j = 0; j < 3; ++j)
                {
                    CGeoCoordinate<Real> aCenterCoordinate = (m_stlFile.facet(aFacetId).vertex((j + 0) % 3) + m_stlFile.facet(aFacetId).vertex((j + 1) % 3)) * 0.5;

                    std::vector<Integer> sCoordinates;

                    aHashGrid.find(sCoordinates, aCenterCoordinate, aTolerance);

                    for (Integer k = 0; k < static_cast<Integer>(sCoordinates.size()); ++k)
                    {
                        if (sCoordinates.at(k) != aFacetId * 3 + j)
                        {
                            Integer ni = sCoordinates.at(k) / 3;
                            Integer nj = sCoordinates.at(k) % 3;

                            m_stlFile.facet(aFacetId).edge(j).setNeighbor(ni);
                            m_stlFile.facet(aFacetId).edge(j).setWhichVertexNot((nj + 2) % 3);

                            break;
                        }
                    }
                }
            }

            return true;
        }

        template <typename Real>
        ENigMA::mesh::CMshMesh<Real> CStlUtils<Real>::mesh()
        {
            ENigMA::mesh::CMshMesh<Real> aMesh;

            Integer aMaxNodeId = std::numeric_limits<Integer>::max();

            for (Integer i = 0; i < m_stlFile.nbFacets(); ++i)
            {
                Integer aFacetId = m_stlFile.facetId(i);

                ENigMA::mesh::CMshElement<Real> anElement;

                anElement.setElementType(ENigMA::mesh::ET_TRIANGLE);

                anElement.addNodeId(aMaxNodeId);
                anElement.addNodeId(aMaxNodeId);
                anElement.addNodeId(aMaxNodeId);

                Integer anElementId = aFacetId;

                aMesh.addElement(anElementId, anElement);
            }

            Integer aNodeId = 0;

            for (Integer i = 0; i < m_stlFile.nbFacets(); ++i)
            {
                Integer aFacetId = m_stlFile.facetId(i);

                Integer aFirstFacet = aFacetId;

                Integer anElementId = aFacetId;

                for (Integer j = 0; j < 3; ++j)
                {
                    if (aMesh.element(anElementId).nodeId(j) != aMaxNodeId)
                        continue;

                    ENigMA::mesh::CMshNode<Real> aNode;

                    aNode = m_stlFile.facet(aFacetId).vertex(j);

                    aMesh.addNode(aNodeId, aNode);

                    Integer aDirection = 0;
                    Integer bReversed = 0;
                    Integer aFacetNum = aFacetId;
                    Integer aVertexNot = (j + 2) % 3;

                    Integer aPivotVertex, aNextEdge;

                    for (Integer k = 0; k < 100; k++)
                    {
                        if (aVertexNot > 2)
                        {
                            if (aDirection == 0)
                            {
                                aPivotVertex = (aVertexNot + 2) % 3;
                                aNextEdge = aPivotVertex;
                                aDirection = 1;
                            }
                            else
                            {
                                aPivotVertex = (aVertexNot + 1) % 3;
                                aNextEdge = aVertexNot % 3;
                                aDirection = 0;
                            }
                        }
                        else
                        {
                            if (aDirection == 0)
                            {
                                aPivotVertex = (aVertexNot + 1) % 3;
                                aNextEdge = aVertexNot;
                            }
                            else
                            {
                                aPivotVertex = (aVertexNot + 2) % 3;
                                aNextEdge = aPivotVertex;
                            }
                        }

                        aMesh.element(aFacetNum).setNodeId(aPivotVertex, aNodeId);

                        Integer aNextFacet = m_stlFile.facet(aFacetNum).edge(aNextEdge).neighbor();

                        if (aNextFacet == -1)
                        {
                            if (bReversed)
                            {
                                break;
                            }
                            else
                            {
                                aDirection = 1;
                                aVertexNot = (j + 1) % 3;
                                bReversed = 1;
                                aFacetNum = aFirstFacet;
                            }
                        }
                        else if (aNextFacet != aFirstFacet)
                        {
                            aVertexNot = m_stlFile.facet(aFacetNum).edge(aNextEdge).whichVertexNot();
                            aFacetNum = aNextFacet;
                        }
                        else
                        {
                            break;
                        }
                    }

                    aNodeId++;
                }
            }

            return aMesh;
        }

        template <typename Real>
        bool CStlUtils<Real>::getFacetQuality(ENigMA::stl::CStlFacet<Real>& aFacet, Integer& smallEdge, Integer& bigEdge, Real& dmin, Real& dmax, Real& q)
        {
            CGeoVector<Real> v0 = aFacet.vertex(1) - aFacet.vertex(0);
            CGeoVector<Real> v1 = aFacet.vertex(2) - aFacet.vertex(1);
            CGeoVector<Real> v2 = aFacet.vertex(0) - aFacet.vertex(2);

            Real d0 = v0.norm();
            Real d1 = v1.norm();
            Real d2 = v2.norm();

            dmin = d0;
            dmax = d0;

            bigEdge = 0;

            if (d1 > dmax)
            {
                dmax = d1;
                bigEdge = 1;
            }

            if (d2 > dmax)
            {
                dmax = d2;
                bigEdge = 2;
            }

            smallEdge = 0;

            if (d1 < dmin)
            {
                dmin = d1;
                smallEdge = 1;
            }

            if (d2 < dmin)
            {
                dmin = d2;
                smallEdge = 2;
            }

            Real a = v0.angle(-v2);
            Real b = v1.angle(-v0);
            Real c = v2.angle(-v1);

            Real d = sin(a) + sin(b) + sin(c);

            if (d > std::numeric_limits<Real>::epsilon())
                q = 4.0 * sin(a) * sin(b) * sin(c) / d;
            else
                q = 0.0;

            return true;
        }

        template <typename Real>
        bool CStlUtils<Real>::splitFacets(const Real splitSize, const Real fq)
        {
            CGeoCoordinate<Real> pn;
            CGeoCoordinate<Real> p3;
            CGeoCoordinate<Real> p[3];

            Real fr;
            Real dx, dy, dz;

            Integer ai;
            Integer neighbor, vertexNot;

            Integer neighbor_c, vertex_not_c;

            for (Integer i = 0; i < m_stlFile.nbFacets(); ++i)
            {
                Integer aFacetId = m_stlFile.facetId(i);

                Integer sedg, bedg;
                Real dmin, dmax;
                Real q;

                this->getFacetQuality(m_stlFile.facet(aFacetId), sedg, bedg, dmin, dmax, q);

                for (Integer j = 0; j < 3; ++j)
                    p[j] = m_stlFile.facet(aFacetId).vertex(j);

                if (q >= fq && dmax > splitSize)
                {
                    ai = (Integer)(dmax / splitSize);

                    if (ai > 1 && ai < 10)
                    {
                        if (ai % 2 != 0)
                            fr = 0.5f - 0.5f / ai;
                        else
                            fr = 0.5f;
                    }
                    else
                        fr = 0.5;

                    if (fr == 1.0)
                        continue;

                    dx = p[(bedg + 1) % 3].x() - p[(bedg + 0) % 3].x();
                    dy = p[(bedg + 1) % 3].y() - p[(bedg + 0) % 3].y();
                    dz = p[(bedg + 1) % 3].z() - p[(bedg + 0) % 3].z();

                    pn.x() = p[(bedg + 0) % 3].x() + dx * fr;
                    pn.y() = p[(bedg + 0) % 3].y() + dy * fr;
                    pn.z() = p[(bedg + 0) % 3].z() + dz * fr;

                    neighbor = m_stlFile.facet(aFacetId).edge((bedg + 0) % 3).neighbor();
                    vertexNot = m_stlFile.facet(aFacetId).edge((bedg + 0) % 3).whichVertexNot();

                    if (vertexNot > 2 || vertexNot < 0)
                        continue;

                    // New facet ids
                    Integer aNewFacetId1 = m_stlFile.nextFacetId() + 0;
                    Integer aNewFacetId2 = m_stlFile.nextFacetId() + 1;

                    // New facet
                    CStlFacet<Real> aNewFacet1;

                    aNewFacet1.addVertex(pn);
                    aNewFacet1.addVertex(p[(bedg + 1) % 3]);
                    aNewFacet1.addVertex(p[(bedg + 2) % 3]);

                    aNewFacet1.calculateNormal();

                    if (neighbor != -1)
                        aNewFacet1.edge(0).setNeighbor(aNewFacetId2);
                    else
                        aNewFacet1.edge(0).setNeighbor(-1);

                    aNewFacet1.edge(1).setNeighbor(m_stlFile.facet(aFacetId).edge((bedg + 1) % 3).neighbor());
                    aNewFacet1.edge(2).setNeighbor(i);

                    aNewFacet1.edge(0).setWhichVertexNot(1);
                    aNewFacet1.edge(1).setWhichVertexNot(m_stlFile.facet(aFacetId).edge((bedg + 1) % 3).whichVertexNot());
                    aNewFacet1.edge(2).setWhichVertexNot((bedg + 0) % 3);

                    // Correct neighbor
                    neighbor_c = m_stlFile.facet(aFacetId).edge((bedg + 1) % 3).neighbor();
                    vertex_not_c = m_stlFile.facet(aFacetId).edge((bedg + 1) % 3).whichVertexNot();

                    if (neighbor_c != -1)
                    {
                        m_stlFile.facet(neighbor_c).edge((vertex_not_c + 1) % 3).setNeighbor(m_stlFile.nbFacets());
                        m_stlFile.facet(neighbor_c).edge((vertex_not_c + 1) % 3).setWhichVertexNot(0);
                    }

                    // Facet: i
                    m_stlFile.facet(aFacetId).setVertex((bedg + 1) % 3, pn);
                    m_stlFile.facet(aFacetId).edge((bedg + 1) % 3).setNeighbor(m_stlFile.nbFacets());
                    m_stlFile.facet(aFacetId).edge((bedg + 1) % 3).setWhichVertexNot(1);

                    // Add new facet 1
                    this->addFacet(aNewFacetId1, aNewFacet1);

                    if (neighbor != -1)
                    {
                        p3 = m_stlFile.facet(neighbor).vertex(vertexNot);

                        // New facet
                        CStlFacet<Real> aNewFacet2;

                        aNewFacet2.addVertex(pn);
                        aNewFacet2.addVertex(p3);
                        aNewFacet2.addVertex(p[(bedg + 1) % 3]);

                        aNewFacet2.calculateNormal();

                        aNewFacet2.edge(0).setNeighbor(neighbor);
                        aNewFacet2.edge(1).setNeighbor(m_stlFile.facet(neighbor).edge((vertexNot + 0) % 3).neighbor());
                        aNewFacet2.edge(2).setNeighbor(aNewFacetId1);

                        aNewFacet2.edge(0).setWhichVertexNot((vertexNot + 2) % 3);
                        aNewFacet2.edge(1).setWhichVertexNot(m_stlFile.facet(neighbor).edge((vertexNot + 0) % 3).whichVertexNot());
                        aNewFacet2.edge(2).setWhichVertexNot(2);

                        // Correct neighbor
                        neighbor_c = m_stlFile.facet(neighbor).edge((vertexNot + 0) % 3).neighbor();
                        vertex_not_c = m_stlFile.facet(neighbor).edge((vertexNot + 0) % 3).whichVertexNot();

                        if (neighbor_c != -1)
                        {
                            m_stlFile.facet(neighbor_c).edge((vertex_not_c + 1) % 3).setNeighbor(m_stlFile.nbFacets());
                            m_stlFile.facet(neighbor_c).edge((vertex_not_c + 1) % 3).setWhichVertexNot(0);
                        }

                        // Facet: neighbor
                        m_stlFile.facet(neighbor).setVertex((vertexNot + 1) % 3, pn);
                        m_stlFile.facet(neighbor).edge((vertexNot + 0) % 3).setNeighbor(m_stlFile.nbFacets());
                        m_stlFile.facet(neighbor).edge((vertexNot + 0) % 3).setWhichVertexNot(2);

                        // Add new facet 2
                        this->addFacet(aNewFacetId2, aNewFacet2);
                    }
                }
            }

            return true;
        }

        template <typename Real>
        bool CStlUtils<Real>::relaxVertices()
        {
            std::vector<CGeoCoordinate<Real>> sVertices;
            std::vector<CStlFacet<Real>> aFacet_i;
            std::vector<CStlFacet<Real>> aFacet_f;
            std::vector<Integer> whichVertex;
            std::vector<Integer> whichFacet;

            const Integer nMaxSize = 100;

            for (Integer i = 0; i < m_stlFile.nbFacets(); ++i)
            {
                Integer aFacetId = m_stlFile.facetId(i);

                m_stlFile.facet(aFacetId).extra[0] = '?';
            }

            for (Integer i = 0; i < m_stlFile.nbFacets(); ++i)
            {
                Integer aFacetId = m_stlFile.facetId(i);

                for (Integer j = 0; j < 3; ++j)
                {
                    if (m_stlFile.facet(aFacetId).edge(j).outline() == true)
                        continue;

                    Integer aNextFacet = aFacetId;
                    Integer aNextEdge = j;

                    sVertices.clear();
                    aFacet_i.clear();
                    aFacet_f.clear();
                    whichFacet.clear();
                    whichVertex.clear();

                    Integer n = 0;

                    bool fixedNode = false;

                    do
                    {
                        if (aNextFacet == -1)
                        {
                            fixedNode = true;
                            break;
                        }

                        if (m_stlFile.facet(aNextFacet).edge(aNextEdge).outline() == true)
                        {
                            fixedNode = true;
                            break;
                        }

                        Integer aNeighbor = m_stlFile.facet(aNextFacet).edge(aNextEdge).neighbor();
                        Integer aVertexNot = m_stlFile.facet(aNextFacet).edge(aNextEdge).whichVertexNot();

                        if (aNeighbor == -1)
                        {
                            fixedNode = true;
                            break;
                        }

                        sVertices.push_back(m_stlFile.facet(aNeighbor).vertex(aVertexNot));
                        aFacet_i.push_back(m_stlFile.facet(aNeighbor));
                        aFacet_f.push_back(m_stlFile.facet(aNeighbor));
                        whichFacet.push_back(aNeighbor);
                        whichVertex.push_back((aVertexNot + 1) % 3);
                        n++;

                        if (n >= nMaxSize)
                        {
                            fixedNode = true;
                            break;
                        }

                        aNextFacet = aNeighbor;
                        aNextEdge = (aVertexNot + 2) % 3;

                    } while (aNextFacet != i);

                    if (fixedNode == true)
                        continue;

                    CGeoCoordinate<Real> aNewVertex(0.0, 0.0, 0.0);

                    for (Integer k = 0; k < static_cast<Integer>(sVertices.size()); ++k)
                    {
                        aNewVertex += sVertices[k];
                    }

                    if (sVertices.size() > 0)
                        aNewVertex /= static_cast<Integer>(sVertices.size());

                    fixedNode = false;

                    for (Integer k = 0; k < static_cast<Integer>(whichVertex.size()); ++k)
                    {
                        CGeoNormal<Real> normal_i, normal_f;

                        aFacet_i[k].calculateNormal();
                        normal_i = aFacet_i[k].normal();

                        aFacet_f[k].setVertex((whichVertex[k] + 1) % 3, aNewVertex);

                        aFacet_f[k].calculateNormal();
                        normal_f = aFacet_f[k].normal();

                        if (std::fabs(normal_f.angle(normal_i)) > 0.5)
                            fixedNode = true;
                    }

                    if (fixedNode == false)
                    {
                        for (Integer k = 0; k < static_cast<Integer>(whichFacet.size()); ++k)
                            m_stlFile.facet(whichFacet[k]).setVertex((whichVertex[k] + 1) % 3, aNewVertex);
                    }
                }
            }

            return true;
        }

        template <typename Real>
        bool CStlUtils<Real>::checkEdges(const Real angle)
        {
            for (Integer i = 0; i < m_stlFile.nbFacets(); ++i)
            {
                Integer aFacetId = m_stlFile.facetId(i);

                m_stlFile.facet(aFacetId).calculateNormal();
            }

            for (Integer i = 0; i < m_stlFile.nbFacets(); ++i)
            {
                Integer aFacetId = m_stlFile.facetId(i);

                for (Integer j = 0; j < 3; ++j)
                {
                    CGeoCoordinate<Real> p1 = m_stlFile.facet(aFacetId).vertex((j + 0) % 3);
                    CGeoCoordinate<Real> p2 = m_stlFile.facet(aFacetId).vertex((j + 1) % 3);

                    m_stlFile.facet(aFacetId).edge(j).setStartPoint(p1);
                    m_stlFile.facet(aFacetId).edge(j).setEndPoint(p2);

                    m_stlFile.facet(aFacetId).edge(j).setOutline(false);

                    if (m_stlFile.facet(aFacetId).edge(j).neighbor() != -1)
                    {
                        m_stlFile.facet(aFacetId).edge(j).setNaked(false);

                        if (std::fabs(m_stlFile.facet(aFacetId).normal().angle(m_stlFile.facet(m_stlFile.facet(aFacetId).edge(j).neighbor()).normal())) > angle)
                            m_stlFile.facet(aFacetId).edge(j).setOutline(true);
                    }
                    else
                        m_stlFile.facet(aFacetId).edge(j).setNaked(true);
                }
            }

            return true;
        }

        template <typename Real>
        bool CStlUtils<Real>::flipEdges(const Real swapAngle)
        {
            for (Integer i = 0; i < m_stlFile.nbFacets(); ++i)
            {
                Integer aFacetId = m_stlFile.facetId(i);

                m_stlFile.facet(aFacetId).extra[0] = '?';
            }

            for (Integer i = 0; i < m_stlFile.nbFacets(); ++i)
            {
                Integer aFacetId = m_stlFile.facetId(i);

                if (m_stlFile.facet(aFacetId).extra[0] == 'c')
                    continue;

                CGeoCoordinate<Real> p0 = m_stlFile.facet(aFacetId).vertex(0);
                CGeoCoordinate<Real> p1 = m_stlFile.facet(aFacetId).vertex(1);
                CGeoCoordinate<Real> p2 = m_stlFile.facet(aFacetId).vertex(2);

                m_stlFile.facet(aFacetId).calculateNormal();
                CGeoNormal<Real> normal_i = m_stlFile.facet(aFacetId).normal();

                for (Integer j = 0; j < 3; ++j)
                {
                    Integer neighbor = m_stlFile.facet(aFacetId).edge(j).neighbor();
                    Integer vertexNot = m_stlFile.facet(aFacetId).edge(j).whichVertexNot();

                    if (neighbor == -1)
                        continue;

                    if (m_stlFile.facet(neighbor).extra[0] == 'c')
                        continue;

                    if (vertexNot < 0 || vertexNot > 2)
                        continue;

                    CGeoCoordinate<Real> p3 = m_stlFile.facet(neighbor).vertex(vertexNot);

                    m_stlFile.facet(neighbor).calculateNormal();
                    CGeoNormal<Real> normal_j = m_stlFile.facet(neighbor).normal();

                    Real ang = normal_i.angle(normal_j);

                    m_stlFile.facet(aFacetId).calculateArea();
                    Real area = m_stlFile.facet(aFacetId).area();

                    if (std::fabs(ang) < swapAngle || area < std::numeric_limits<Real>::min())
                    {
                        Real dmin, dmax;

                        Integer sedg, bedg;
                        Real qi, qj;

                        this->getFacetQuality(m_stlFile.facet(aFacetId), sedg, bedg, dmin, dmax, qi);

                        m_stlFile.facet(neighbor).calculateNormal();
                        CGeoNormal<Real> iNormal_i = m_stlFile.facet(aFacetId).normal();

                        this->getFacetQuality(m_stlFile.facet(neighbor), sedg, bedg, dmin, dmax, qj);

                        m_stlFile.facet(neighbor).calculateNormal();
                        CGeoNormal<Real> iNormal_j = m_stlFile.facet(neighbor).normal();

                        Real qmini = std::min(qi, qj);

                        CStlFacet<Real> newFacet_i, newFacet_j;

                        switch (j)
                        {
                        case 0:

                            newFacet_i = m_stlFile.facet(aFacetId);
                            newFacet_i.addVertex(p0);
                            newFacet_i.addVertex(p3);
                            newFacet_i.addVertex(p2);

                            newFacet_j = m_stlFile.facet(neighbor);
                            newFacet_j.addVertex(p1);
                            newFacet_j.addVertex(p2);
                            newFacet_j.addVertex(p3);

                            break;

                        case 1:

                            newFacet_i = m_stlFile.facet(aFacetId);
                            newFacet_i.addVertex(p0);
                            newFacet_i.addVertex(p1);
                            newFacet_i.addVertex(p3);

                            newFacet_j = m_stlFile.facet(neighbor);
                            newFacet_j.addVertex(p0);
                            newFacet_j.addVertex(p3);
                            newFacet_j.addVertex(p2);

                            break;

                        case 2:

                            newFacet_i = m_stlFile.facet(aFacetId);
                            newFacet_i.addVertex(p0);
                            newFacet_i.addVertex(p1);
                            newFacet_i.addVertex(p3);

                            newFacet_j = m_stlFile.facet(neighbor);
                            newFacet_j.addVertex(p1);
                            newFacet_j.addVertex(p2);
                            newFacet_j.addVertex(p3);

                            break;
                        }

                        this->getFacetQuality(newFacet_i, sedg, bedg, dmin, dmax, qi);

                        newFacet_i.calculateNormal();
                        CGeoNormal<Real> fNormal_i = newFacet_i.normal();

                        this->getFacetQuality(newFacet_j, sedg, bedg, dmin, dmax, qj);

                        newFacet_j.calculateNormal();
                        CGeoNormal<Real> fNormal_j = newFacet_j.normal();

                        Real qminf = std::min(qi, qj);

                        m_stlFile.facet(aFacetId).calculateArea();
                        m_stlFile.facet(neighbor).calculateArea();

                        newFacet_i.calculateArea();
                        newFacet_j.calculateArea();

                        Real areai = m_stlFile.facet(aFacetId).area() + m_stlFile.facet(neighbor).area();
                        Real areaf = newFacet_i.area() + newFacet_j.area();

                        if (qminf > qmini)
                        {
                            Real cangi;

                            if (qmini > 1E-15)
                                cangi = iNormal_i.angle(fNormal_i);
                            else
                                cangi = 0.0f;

                            Real cangj = iNormal_j.angle(fNormal_j);

                            if (std::fabs(cangi) < 2.5 && std::fabs(cangj) < 2.5 && std::fabs(areai - areaf) < areai * 0.01)
                            {
                                m_stlFile.facet(aFacetId) = newFacet_i;
                                m_stlFile.facet(neighbor) = newFacet_j;

                                m_stlFile.facet(aFacetId).extra[0] = 'c';
                                m_stlFile.facet(neighbor).extra[0] = 'c';

                                break;
                            }
                        }
                    }
                }
            }

            return true;
        }

        template <typename Real>
        bool CStlUtils<Real>::collapseEdges(const Real minQuality, const Real maxAngle)
        {
            std::vector<CStlFacet<Real>> aFacet_i;
            std::vector<CStlFacet<Real>> aFacet_f;
            std::vector<Integer> whichVertex;
            std::vector<Integer> whichFacet;

            const Integer nMaxSize = 100;

            for (Integer i = 0; i < m_stlFile.nbFacets(); ++i)
            {
                Integer aFacetId = m_stlFile.facetId(i);

                m_stlFile.facet(aFacetId).extra[0] = '?';
            }

            for (Integer i = 0; i < m_stlFile.nbFacets(); ++i)
            {
                Integer aFacetId = m_stlFile.facetId(i);

                if (m_stlFile.facet(aFacetId).extra[0] == 'c')
                    continue;

                Real dmin, dmax;

                Integer sedgi, bedgi;
                Real qi;

                this->getFacetQuality(m_stlFile.facet(aFacetId), sedgi, bedgi, dmin, dmax, qi);

                if (qi > minQuality)
                    continue;

                Integer pair = m_stlFile.facet(aFacetId).edge(sedgi).neighbor();

                if (pair == -1)
                    continue;

                Integer sedgj, bedgj;
                Real qj;

                this->getFacetQuality(m_stlFile.facet(pair), sedgj, bedgj, dmin, dmax, qj);

                // collapse edge
                Integer aNextFacet = aFacetId;
                Integer aNextEdge = (sedgi + 1) % 3;

                aFacet_i.clear();
                aFacet_f.clear();
                whichFacet.clear();
                whichVertex.clear();

                Integer n = 0;

                bool fixedNode = false;

                do
                {
                    Integer aNeighborId = m_stlFile.facet(aNextFacet).edge(aNextEdge).neighbor();
                    Integer aVertexNot = m_stlFile.facet(aNextFacet).edge(aNextEdge).whichVertexNot();

                    if (aNeighborId != -1)
                    {
                        aFacet_i.push_back(m_stlFile.facet(aNeighborId));
                        aFacet_f.push_back(m_stlFile.facet(aNeighborId));
                        whichFacet.push_back(aNeighborId);
                        whichVertex.push_back((aVertexNot + 1) % 3);
                        n++;

                        if (n >= nMaxSize)
                        {
                            fixedNode = true;
                            break;
                        }
                    }
                    else
                    {
                        fixedNode = true;
                        break;
                    }

                    aNextFacet = aNeighborId;
                    aNextEdge = (aVertexNot + 2) % 3;

                    if (aNextFacet == -1)
                    {
                        fixedNode = true;
                        break;
                    }

                } while (aNextFacet != i);

                if (fixedNode)
                    continue;

                CGeoCoordinate<Real> aNewVertex = m_stlFile.facet(aFacetId).vertex(sedgi);

                fixedNode = false;

                for (Integer k = 0; k < static_cast<Integer>(aFacet_f.size()); ++k)
                {
                    aFacet_f[k].setVertex((whichVertex[k] + 1) % 3, aNewVertex);

                    aFacet_f[k].calculateArea(true);

                    if (aFacet_f[k].area() > 0.0)
                    {
                        aFacet_f[k].calculateNormal(true);
                        aFacet_i[k].calculateNormal(true);

                        Real cang = aFacet_f[k].normal().angle(aFacet_i[k].normal());

                        if (std::fabs(cang) > maxAngle)
                            fixedNode = true;
                    }
                }

                if (!fixedNode)
                {
                    for (Integer k = 0; k < static_cast<Integer>(whichFacet.size()); ++k)
                        m_stlFile.facet(whichFacet[k]).setVertex((whichVertex[k] + 1) % 3, aNewVertex);

                    m_stlFile.facet(aFacetId).extra[0] = 'c';
                    m_stlFile.facet(pair).extra[0] = 'c';
                }
            }

            return true;
        }

        template <typename Real>
        void CStlUtils<Real>::invertFacet(Integer aFacetId, Integer aVertexNot)
        {
            // Fix facet
            CGeoCoordinate<Real> anAux = m_stlFile.facet(aFacetId).vertex((aVertexNot + 2) % 3);
            m_stlFile.facet(aFacetId).setVertex((aVertexNot + 2) % 3, m_stlFile.facet(aFacetId).vertex((aVertexNot + 1) % 3));
            m_stlFile.facet(aFacetId).setVertex((aVertexNot + 1) % 3, anAux);

            Integer aNeighborId2 = m_stlFile.facet(aFacetId).edge((aVertexNot + 2) % 3).neighbor();
            Integer aNeighborId3 = m_stlFile.facet(aFacetId).edge((aVertexNot + 3) % 3).neighbor();

            Integer aVertexNot2 = m_stlFile.facet(aFacetId).edge((aVertexNot + 2) % 3).whichVertexNot();
            Integer aVertexNot3 = m_stlFile.facet(aFacetId).edge((aVertexNot + 3) % 3).whichVertexNot();

            // Fix the vnots of the neighboring facets
            if (aVertexNot2 != -1)
            {
                if (aNeighborId2 != -1)
                    m_stlFile.facet(aNeighborId2).edge((aVertexNot2 + 1) % 3).setWhichVertexNot((m_stlFile.facet(aNeighborId2).edge((aVertexNot2 + 1) % 3).whichVertexNot() + 4) % 3);
                else
                    m_stlFile.facet(aNeighborId2).edge((aVertexNot2 + 1) % 3).setWhichVertexNot(-1);
            }

            if (aVertexNot3 != -1)
            {
                if (aNeighborId3 != -1)
                    m_stlFile.facet(aNeighborId3).edge((aVertexNot3 + 1) % 3).setWhichVertexNot((m_stlFile.facet(aNeighborId3).edge((aVertexNot3 + 1) % 3).whichVertexNot() + 2) % 3);
                else
                    m_stlFile.facet(aNeighborId3).edge((aVertexNot3 + 1) % 3).setWhichVertexNot(-1);
            }

            // Swap the neighbors of the facet that is being reversed
            m_stlFile.facet(aFacetId).edge((aVertexNot + 2) % 3).setNeighbor(aNeighborId3);
            m_stlFile.facet(aFacetId).edge((aVertexNot + 3) % 3).setNeighbor(aNeighborId2);

            // Swap the vnots of the facet that is being reversed
            m_stlFile.facet(aFacetId).edge((aVertexNot + 2) % 3).setWhichVertexNot(aVertexNot3);
            m_stlFile.facet(aFacetId).edge((aVertexNot + 3) % 3).setWhichVertexNot(aVertexNot2);
        }

        template <typename Real>
        bool CStlUtils<Real>::setOrientation(Integer aPivotFacetId, const Real aTolerance)
        {
            std::vector<Integer> sFacetsToFix;

            sFacetsToFix.push_back(aPivotFacetId);

            for (Integer i = 0; i < static_cast<Integer>(sFacetsToFix.size()); ++i)
            {
                Integer aFacetId = sFacetsToFix[i];

                for (Integer j = 0; j < 3; ++j)
                {
                    Integer aNeighborId = m_stlFile.facet(aFacetId).edge(j).neighbor();
                    Integer aVertexNot = m_stlFile.facet(aFacetId).edge(j).whichVertexNot();

                    if (aNeighborId != -1)
                    {
                        if (std::find(sFacetsToFix.begin(), sFacetsToFix.end(), aNeighborId) == sFacetsToFix.end())
                        {
                            sFacetsToFix.push_back(aNeighborId);

                            Real d1 = (m_stlFile.facet(aNeighborId).vertex((aVertexNot + 1) % 3) - m_stlFile.facet(aFacetId).vertex((j + 0) % 3)).norm();
                            Real d2 = (m_stlFile.facet(aNeighborId).vertex((aVertexNot + 2) % 3) - m_stlFile.facet(aFacetId).vertex((j + 1) % 3)).norm();

                            if (d1 <= aTolerance && d2 <= aTolerance)
                                this->invertFacet(aNeighborId, aVertexNot);
                        }
                    }
                }
            }

            return true;
        }

        template <typename Real>
        bool CStlUtils<Real>::setOrientation(bool bPointOut, const Real aTolerance)
        {
            Integer aFacetId = m_stlFile.facetId(0);
            CStlFacet<Real>& aFacet = m_stlFile.facet(aFacetId);

            aFacet.calculateNormal();
            aFacet.calculateCentroid();

            CGeoCoordinate<Real> aStartPoint = aFacet.centroid();
            CGeoCoordinate<Real> anEndPoint = aStartPoint + aFacet.normal() * m_stlFile.stats.size.norm();

            CGeoIntersectionType anIntersectionType;
            CGeoLine<Real> aLine(aStartPoint, anEndPoint);
            CGeoCoordinate<Real> aNewPoint;

            Integer nbIntersections = 0;

            for (Integer i = 1; i < m_stlFile.nbFacets(); ++i)
            {
                aFacetId = m_stlFile.facetId(i);
                aFacet = m_stlFile.facet(aFacetId);

                aFacet.calculateNormal();

                if (aFacet.intersects(aLine, aNewPoint, anIntersectionType, aTolerance))
                    nbIntersections++;
            }

            //std::cout << nbIntersections << std::endl;

            if (nbIntersections % 2 != 0)
                this->invertFacet(aFacetId, 0);

            return this->setOrientation(aFacetId, aTolerance);
        }
    }
}
