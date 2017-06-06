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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef CPU86
#include <float.h>
#endif /* CPU86 */
#ifdef LINUX
#include <fpu_control.h>
#endif /* LINUX */

#include "GeoCoordinate.hpp"

namespace ENigMA
{

    namespace geometry
    {

        template <typename Real>
        class CGeoPredicates
        {
        public:

            static Real exactInit();
            static Real inCircle(CGeoCoordinate<Real>& aPoint1, CGeoCoordinate<Real>& aPoint2, CGeoCoordinate<Real>& aPoint3, CGeoCoordinate<Real>& aPoint);
            static Real inSphere(CGeoCoordinate<Real>& aPoint1, CGeoCoordinate<Real>& aPoint2, CGeoCoordinate<Real>& aPoint3, CGeoCoordinate<Real>& aPoint4, CGeoCoordinate<Real>& aPoint);

        private:

            static int growExpansion(int elen, Real *e, Real b, Real *h);
            static int growExpansionZeroelim(int elen, Real *e, Real b, Real *h);
            static int expansionSum(int elen, Real *e, int flen, Real *f, Real *h);
            static int expansionSumZeroelim1(int elen, Real *e, int flen, Real *f, Real *h);
            static int expansionSumZeroelim2(int elen, Real *e, int flen, Real *f, Real *h);
            static int fastExpansionSum(int elen, Real *e, int flen, Real *f, Real *h);
            static int fastExpansionSumZeroelim(int elen, Real *e, int flen, Real *f, Real *h);
            static int linearExpansionSum(int elen, Real *e, int flen, Real *f, Real *h);
            static int linearExpansionSumZeroelim(int elen, Real *e, int flen, Real *f, Real *h);
            static int scaleExpansion(int elen, Real *e, Real b, Real *h);
            static int scaleExpansionZeroelim(int elen, Real *e, Real b, Real *h);
            static int compress(int elen, Real *e, Real *h);
            static Real estimate(int elen, Real *e);
            static Real orient2dFast(Real *pa, Real *pb, Real *pc);
            static Real orient2dExact(Real *pa, Real *pb, Real *pc);
            static Real orient2dSlow(Real *pa, Real *pb, Real *pc);
            static Real orient2dAdapt(Real *pa, Real *pb, Real *pc, Real detsum);
            static Real orient2d(Real *pa, Real *pb, Real *pc);
            static Real orient3dFast(Real *pa, Real *pb, Real *pc, Real *pd);
            static Real orient3dExact(Real *pa, Real *pb, Real *pc, Real *pd);
            static Real orient3dSlow(Real *pa, Real *pb, Real *pc, Real *pd);
            static Real orient3dAdapt(Real *pa, Real *pb, Real *pc, Real *pd, Real permanent);
            static Real orient3d(Real *pa, Real *pb, Real *pc, Real *pd);
            static Real inCircleFast(Real *pa, Real *pb, Real *pc, Real *pd);
            static Real inCircleExact(Real *pa, Real *pb, Real *pc, Real *pd);
            static Real inCircleSlow(Real *pa, Real *pb, Real *pc, Real *pd);
            static Real inCircleAdapt(Real *pa, Real *pb, Real *pc, Real *pd, Real permanent);
            static Real inCircle(Real *pa, Real *pb, Real *pc, Real *pd);
            static Real inSphereFast(Real *pa, Real *pb, Real *pc, Real *pd, Real *pe);
            static Real inSphereExact(Real *pa, Real *pb, Real *pc, Real *pd, Real *pe);
            static Real inSphereSlow(Real *pa, Real *pb, Real *pc, Real *pd, Real *pe);
            static Real inSphereAdapt(Real *pa, Real *pb, Real *pc, Real *pd, Real *pe, Real permanent);
            static Real inSphere(Real *pa, Real *pb, Real *pc, Real *pd, Real *pe);

        };

    }

}

#include "GeoPredicates_Imp.hpp"
