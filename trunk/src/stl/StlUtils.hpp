// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "MshMesh.hpp"
#include "StlFile.hpp"

namespace ENigMA {

namespace stl {

#define LABEL_SIZE 80
#define NUM_FACET_SIZE 4
#define HEADER_SIZE 84
#define STL_MIN_FILE_SIZE 134
#define ASCII_LINES_PER_FACET 7
#define SIZEOF_EDGE_SORT 24
#define SIZEOF_STL_FACET 50

    template <typename Real>
    struct CDisjointPoint {

        ENigMA::geometry::CGeoCoordinate<Real> aStartPoint;
        ENigMA::geometry::CGeoCoordinate<Real> aPoint;

        bool operator<(const CDisjointPoint& a) const
        {
            return (aPoint - aStartPoint).norm() < (a.aPoint - aStartPoint).norm();
        }
    };

    template <typename Real>
    struct CDisjointLine {
        ENigMA::geometry::CGeoLine<Real> aLine;
        std::vector<CDisjointPoint<Real>> sPoints;

        Integer wFacet;
        Integer wEdge;
    };

    template <typename Real>
    class CStlUtils {
    private:
        CStlFile<Real> m_stlFile;

        bool saveAscii(std::string strFileName);
        bool saveBinary(std::string strFileName);

        bool writeFloat(std::ofstream& stream, const Real aValue);
        bool writeInt(std::ofstream& stream, const Integer aValue);

        Real getFloat(std::ifstream& stream);
        Integer getInt(std::ifstream& stream);

        void invertFacet(Integer aFacetId, Integer aVertexNot);

    public:
        CStlUtils();
        CStlUtils(ENigMA::mesh::CMshMesh<Real>& aTriMesh);

        ~CStlUtils();

        ENigMA::stl::CStlFile<Real>& stlFile();

        void set(ENigMA::mesh::CMshMesh<Real>& aTriMesh);

        bool load(std::string strFileName);
        bool save(std::string strFileName);

        bool add(CStlUtils<Real>& aStl);

        bool addFacet(const Integer aFacetId, ENigMA::stl::CStlFacet<Real> aFacet);
        bool removeFacet(const Integer aFacetId);

        bool calculateStatistics();

        bool removeDuplicateFacets(const Real aTolerance = 0.0);
        bool removeInvalidFacets(const Real aTolerance = 0.0);

        bool generateConnectivity(const Real aTolerance = 0.0);
        bool setOrientation(Integer aPivotFacetId, const Real aTolerance = 0.0);
        bool setOrientation(bool bPointOut = true, const Real aTolerance = 0.0);

        ENigMA::mesh::CMshMesh<Real> mesh();

        bool getFacetQuality(ENigMA::stl::CStlFacet<Real>& aFacet, Integer& smallEdge, Integer& bigEdge, Real& dmin, Real& dmax, Real& q);

        bool splitFacets(const Real splitSize, const Real fq);
        bool relaxVertices();
        bool checkEdges(const Real angle);
        bool flipEdges(const Real swapAngle);
        bool collapseEdges(const Real minQuality, const Real maxAngle);
    };
}
}

#include "StlUtils_Imp.hpp"
