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

#include "GeoPolyhedron.hpp"
#include "FvmCell.hpp"
#include "FvmFace.hpp"

using namespace ENigMA::geometry;

namespace ENigMA
{

    namespace fvm
    {

        template <typename Real>
        class CFvmControlVolume : public CFvmCell<Real>, public CGeoVolume<Real>
        {
        private:

            bool m_clipped;

            Integer m_clippedFaceId;
            CFvmFace<Real> m_clippedFace;

            CGeoPolyhedron<Real> m_polyhedron;
            CGeoPolyhedron<Real> m_clippedPolyhedron;

            Integer m_controlVolumeId;

            Real m_originalVolume;

        public:
            CFvmControlVolume();
            CFvmControlVolume(CGeoPolyhedron<Real>& aPolyhedron);
            ~CFvmControlVolume();

            void setControlVolumeId(const Integer aControlVolumeId);
            Integer controlVolumeId();

            Integer nbFaces();
            void addFace(const Integer aFaceId, CFvmFace<Real>& aFace);
            Integer faceId(const Integer aFaceIndex);

            bool containsFace(const Integer aFaceId);

            void calculateFaceArea(const Integer aFaceId, bool bRecalculate = false);
            void calculateFaceCentroid(const Integer aFaceId, bool bRecalculate = false);

            Real faceArea(const Integer aFaceId);
            CGeoNormal<Real>& faceNormal(const Integer aFaceId);
            Real faceDist(const Integer aFaceId);

            void reset();

            // Clip control volume
            void setClippedFaceId(Integer aFaceId);
            Integer clippedFaceId();
            CFvmFace<Real>& clippedFace();
            inline void clip(CGeoNormal<Real>& aNormal, const Real volumeFractionReq, Real& volumeFractionAct, Integer& nIterations, const Integer nMaxIterations, const Real aTolerance);

            void calculateOriginalVolume(bool bReCalculate = false);
            Real originalVolume();
            Real originalAvFaceDist(const Integer aFaceId);

            void calculateCentroid(bool bReCalculate = false);
            void calculateSurfaceArea(bool bReCalculate = false);
            void calculateVolume(bool bReCalculate = false);
            void calculateBoundingBox(bool bReCalculate = false);

            bool isClipped();
            void setClipped(bool clipped);

        };

    }

}

#include "FvmControlVolume_Imp.hpp"
