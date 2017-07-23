// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "GeoPolyhedron.hpp"
#include "FvmCell.hpp"
#include "FvmFace.hpp"

namespace ENigMA
{

    namespace fvm
    {

        template <typename Real>
        class CFvmControlVolume : public CFvmCell<Real>, public ENigMA::geometry::CGeoVolume<Real>
        {
        private:

            bool m_clipped;

            Integer m_clippedFaceId;
            CFvmFace<Real> m_clippedFace;

            ENigMA::geometry::CGeoPolyhedron<Real> m_polyhedron;
            ENigMA::geometry::CGeoPolyhedron<Real> m_clippedPolyhedron;

            Integer m_controlVolumeId;

            Real m_originalVolume;

        public:
            CFvmControlVolume();
            CFvmControlVolume(ENigMA::geometry::CGeoPolyhedron<Real>& aPolyhedron);
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
            ENigMA::geometry::CGeoNormal<Real>& faceNormal(const Integer aFaceId);
            Real faceDist(const Integer aFaceId);

            void reset();

            // Clip control volume
            void setClippedFaceId(Integer aFaceId);
            Integer clippedFaceId();
            CFvmFace<Real>& clippedFace();
            inline void clip(ENigMA::geometry::CGeoNormal<Real>& aNormal, const Real volumeFractionReq, Real& volumeFractionAct, Integer& nIterations, const Integer nMaxIterations, const Real aTolerance);

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
