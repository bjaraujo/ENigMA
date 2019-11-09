// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#include "gtest/gtest.h"

#include "TypeDef.hpp"

#include "GeoConvexHull.hpp"

using namespace ENigMA::geometry;

class CTestGeoConvexHull : public ::testing::Test {
protected:

    virtual void SetUp() {

    }

    virtual void TearDown() {

    }

};

TEST_F(CTestGeoConvexHull, convexHull) {


    std::vector<CGeoCoordinate<decimal> > sVertices;

    CGeoCoordinate<decimal> aVertex1(0.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex2(0.25, 0.2, 0.0);
    CGeoCoordinate<decimal> aVertex3(1.0, 0.0, 0.0);
    CGeoCoordinate<decimal> aVertex4(1.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex5(0.0, 1.0, 0.0);
    CGeoCoordinate<decimal> aVertex6(0.5, 0.5, 0.0);
    CGeoCoordinate<decimal> aVertex7(0.1, 0.1, 0.0);
    
    sVertices.push_back(aVertex1);
    sVertices.push_back(aVertex2);
    sVertices.push_back(aVertex3);
    sVertices.push_back(aVertex4);
    sVertices.push_back(aVertex5);
    sVertices.push_back(aVertex6);
    sVertices.push_back(aVertex7);

    CGeoConvexHull<decimal> aConvexHull(sVertices);
    
    EXPECT_EQ(4, aConvexHull.nbVertices());

}

