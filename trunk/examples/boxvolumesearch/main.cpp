// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include <time.h>

#include "GeoRtree.hpp"
#include "GeoAdtree.hpp"

using namespace ENigMA::geometry;

int main(int argc, char *argv[])
{

    CGeoRtree<double> aRtree;
    CGeoAdtree<double> aAdtree;

    Integer maxBound = 100000;

    CGeoBoundingBox<double> aBoundingBox(0, 0, 0, maxBound, maxBound, maxBound);

    aAdtree.set(aBoundingBox);

    for (Integer n = 0; n < 20; n++)
    {

        std::cout << "************************" << std::endl;
        std::cout << n + 1 << std::endl;

        aRtree.reset();
        aAdtree.reset();

        srand((Integer)time(0));

        std::vector<CGeoBoundingBox<double> > sBoundingBoxes;

        for (Integer i = 0; i < maxBound; ++i)
        {

            CGeoCoordinate<double> aCoordinate1;
            CGeoCoordinate<double> aCoordinate2;

            aCoordinate1.x() = (rand() % maxBound) / (double) maxBound;
            aCoordinate1.y() = (rand() % maxBound) / (double) maxBound;
            aCoordinate1.z() = (rand() % maxBound) / (double) maxBound;

            aCoordinate2.x() = (rand() % maxBound) / (double) maxBound;
            aCoordinate2.y() = (rand() % maxBound) / (double) maxBound;
            aCoordinate2.z() = (rand() % maxBound) / (double) maxBound;

            CGeoBoundingBox<double> aBoundingBox1;

            aBoundingBox1.addCoordinate(aCoordinate1);
            aBoundingBox1.addCoordinate(aCoordinate2);

            sBoundingBoxes.push_back(aBoundingBox1);

            aRtree.addGeometricObject(i, aBoundingBox1);
            aAdtree.addGeometricObject(i, aBoundingBox1);

        }

        CGeoCoordinate<double> aCoordinate3(0.199, 0.199, 0.199);
        CGeoCoordinate<double> aCoordinate4(0.200, 0.200, 0.200);

        CGeoBoundingBox<double> aBoundingBox2;

        aBoundingBox2.addCoordinate(aCoordinate3);
        aBoundingBox2.addCoordinate(aCoordinate4);

        std::vector<Integer> sBoundingBoxes1, sBoundingBoxes2;

        std::cout << "Searching..." << std::endl;

        clock_t start1 = clock();
        for (Integer i = 0; i < 100; ++i)
            aRtree.find(sBoundingBoxes1, aBoundingBox2);
        std::cout << "Time R-Tree: " << static_cast<double>(clock() - start1) / CLOCKS_PER_SEC << " seconds." << std::endl;

        clock_t start2 = clock();
        for (Integer i = 0; i < 100; ++i)
            aAdtree.find(sBoundingBoxes2, aBoundingBox2);
        std::cout << "Time ADTree: " << static_cast<double>(clock() - start2) / CLOCKS_PER_SEC << " seconds." << std::endl;

        std::cout << "Number of elements found by R-Tee: " << sBoundingBoxes1.size() << std::endl;
        std::cout << "Number of elements found by ADTree: " << sBoundingBoxes2.size() << std::endl;

        std::sort(sBoundingBoxes1.begin(), sBoundingBoxes1.end());
        std::sort(sBoundingBoxes2.begin(), sBoundingBoxes2.end());

        if (sBoundingBoxes1 != sBoundingBoxes2)
        {

            int size1 = static_cast<int>(sBoundingBoxes1.size());
            int size2 = static_cast<int>(sBoundingBoxes2.size());

            for (int i = 0; i < std::min(5, size1); ++i)
                std::cout << sBoundingBoxes1[i] << std::endl;

            std::cout << "--" << std::endl;

            for (int i = 0; i < std::min(5, size2); ++i)
                std::cout << sBoundingBoxes2[i] << std::endl;

            std::cout << "-> 0" << std::endl;

            std::cout << sBoundingBoxes[0].min() << std::endl;
            std::cout << sBoundingBoxes[0].max() << std::endl;

            std::cout << "0 Intersects?: ";

            if (aBoundingBox2.intersects(sBoundingBoxes[0]))
                std::cout << "yes" << std::endl;
            else
                std::cout << "no" << std::endl;

            break;

        }

    }

    std::cout << "Done." << std::endl;

}

