// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

using namespace ENigMA::geometry;

namespace ENigMA
{

    namespace fvm
    {

        template <typename Real>
        CFvmControlVolume<Real>::CFvmControlVolume()
            : m_clipped(false)
            , m_clippedFaceId(-1)
            , m_controlVolumeId(0)
            , m_originalVolume(0.0)
        {

        }

        template <typename Real>
        CFvmControlVolume<Real>::CFvmControlVolume(CGeoPolyhedron<Real>& aPolyhedron) : CFvmControlVolume() 
        {

        }

        template <typename Real>
        CFvmControlVolume<Real>::~CFvmControlVolume()
        {

        }

        template <typename Real>
        void CFvmControlVolume<Real>::setControlVolumeId(const Integer aControlVolumeId)
        {

            m_controlVolumeId = aControlVolumeId;

        }

        template <typename Real>
        Integer CFvmControlVolume<Real>::controlVolumeId()
        {

            return m_controlVolumeId;

        }

        template <typename Real>
        Integer CFvmControlVolume<Real>::nbFaces()
        {

            if (m_clipped)
                return m_clippedPolyhedron.nbPolygons();
            else
                return m_polyhedron.nbPolygons();

        }

        template <typename Real>
        void CFvmControlVolume<Real>::addFace(const Integer aFaceId, CFvmFace<Real>& aFace)
        {

            CGeoPolygon<Real> aPolygon;

            aPolygon.setPolyline(aFace.polyline());

            m_polyhedron.addPolygon(aFaceId, aPolygon);

        }

        template <typename Real>
        Integer CFvmControlVolume<Real>::faceId(const Integer aFaceIndex)
        {

            const Integer aPolygonIndex = aFaceIndex;

            if (m_clipped)
                return m_clippedPolyhedron.polygonId(aPolygonIndex);
            else
                return m_polyhedron.polygonId(aPolygonIndex);

        }

        template <typename Real>
        bool CFvmControlVolume<Real>::containsFace(const Integer aFaceId)
        {

            Integer aPolygonId = aFaceId;

            if (m_clipped)
                return m_clippedPolyhedron.containsPolygon(aPolygonId);
            else
                return m_polyhedron.containsPolygon(aPolygonId);

        }

        template <typename Real>
        void CFvmControlVolume<Real>::calculateFaceArea(const Integer aFaceId, bool bRecalculate)
        {

            const Integer aPolygonId = aFaceId;

            if (m_clipped)
                m_clippedPolyhedron.polygon(aPolygonId).calculateArea(bRecalculate);
            else
                m_polyhedron.polygon(aPolygonId).calculateArea(bRecalculate);

        }

        template <typename Real>
        Real CFvmControlVolume<Real>::faceArea(const Integer aFaceId)
        {

            const Integer aPolygonId = aFaceId;

            if (m_clipped)
                return m_clippedPolyhedron.polygon(aPolygonId).area();
            else
                return m_polyhedron.polygon(aPolygonId).area();

        }

        template <typename Real>
        void CFvmControlVolume<Real>::calculateFaceNormal(const Integer aFaceId, bool bRecalculate)
        {

            const Integer aPolygonId = aFaceId;

            if (m_clipped)
                m_clippedPolyhedron.polygon(aPolygonId).calculateNormal(bRecalculate);
            else
                m_polyhedron.polygon(aPolygonId).calculateNormal(bRecalculate);

        }

        template <typename Real>
        CGeoNormal<Real>& CFvmControlVolume<Real>::faceNormal(const Integer aFaceId)
        {

            const Integer aPolygonId = aFaceId;

            if (m_clipped)
                return m_clippedPolyhedron.polygon(aPolygonId).normal();
            else
                return m_polyhedron.polygon(aPolygonId).normal();

        }

        template <typename Real>
        void CFvmControlVolume<Real>::calculateFaceCentroid(const Integer aFaceId, bool bRecalculate)
        {

            const Integer aPolygonId = aFaceId;

            if (m_clipped)
                m_clippedPolyhedron.polygon(aPolygonId).calculateCentroid(bRecalculate);
            else
                m_polyhedron.polygon(aPolygonId).calculateCentroid(bRecalculate);

        }

        template <typename Real>
        CGeoVector<Real> CFvmControlVolume<Real>::faceDist(const Integer aFaceId)
        {
            const Integer aPolygonId = aFaceId;

            if (m_clipped)
                return CGeoVector<Real>(m_clippedPolyhedron.polygon(aPolygonId).centroid() - this->m_centroid);
            else
                return CGeoVector<Real>(m_polyhedron.polygon(aPolygonId).centroid() - this->m_centroid);
        }

        template <typename Real>
        void CFvmControlVolume<Real>::reset()
        {

        }

        template <typename Real>
        void CFvmControlVolume<Real>::setClippedFaceId(Integer aFaceId)
        {
            m_clippedFaceId = aFaceId;
        }

        template <typename Real>
        Integer CFvmControlVolume<Real>::clippedFaceId()
        {
            return m_clippedFaceId;
        }

        template <typename Real>
        CFvmFace<Real>& CFvmControlVolume<Real>::clippedFace()
        {
            return m_clippedFace;
        }

        template <typename Real>
        void CFvmControlVolume<Real>::clip(const CGeoNormal<Real>& aNormal, const Real volumeFractionReq, Real& volumeFractionAct, Integer& nIterations, const Integer nMaxIterations,  const Real aVolumeFractionTolerance, const Real aTolerance)
        {

            CGeoPolygon<Real> aPolygon;
            CGeoPlane<Real> aPlane(aNormal, 0.0);

            m_clippedPolyhedron = m_polyhedron.clip(aPolygon, m_clippedFaceId, aPlane, volumeFractionReq, volumeFractionAct, nIterations, nMaxIterations, aVolumeFractionTolerance, aTolerance);

            m_clippedFace.set(aPolygon);
            m_clipped = true;

        }

        template <typename Real>
        void CFvmControlVolume<Real>::calculateSurfaceArea(bool bReCalculate)
        {

            if (!this->m_bSurfaceArea || bReCalculate)
            {

                if (m_clipped)
                {
                    m_clippedPolyhedron.calculateSurfaceArea();
                    this->m_surfaceArea = m_clippedPolyhedron.surfaceArea();
                }
                else
                {
                    m_polyhedron.calculateSurfaceArea();
                    this->m_surfaceArea = m_polyhedron.surfaceArea();
                }

                this->m_bSurfaceArea = true;

            }

        }

        template <typename Real>
        void CFvmControlVolume<Real>::calculateCentroid(bool bReCalculate)
        {

            if (!this->m_bCentroid || bReCalculate)
            {

                if (m_clipped)
                {
                    m_clippedPolyhedron.calculateCentroid();
                    this->m_centroid = m_clippedPolyhedron.centroid();
                }
                else
                {
                    m_polyhedron.calculateCentroid();
                    this->m_centroid = m_polyhedron.centroid();
                }

                this->m_bCentroid = true;

            }

        }

        template <typename Real>
        void CFvmControlVolume<Real>::calculateVolume(bool bReCalculate)
        {

            if (!this->m_bVolume || bReCalculate)
            {

                if (m_clipped)
                {
                    m_clippedPolyhedron.calculateVolume(bReCalculate);
                    this->m_volume = m_clippedPolyhedron.volume();
                }
                else
                {
                    m_polyhedron.calculateVolume(bReCalculate);
                    this->m_volume = m_polyhedron.volume();
                }

                this->m_bVolume = true;

            }

        }

        template <typename Real>
        void CFvmControlVolume<Real>::calculateBoundingBox(bool bReCalculate)
        {


            // TODO:

        }

        template <typename Real>
        void CFvmControlVolume<Real>::calculateOriginalVolume(bool bReCalculate)
        {

            m_polyhedron.calculateVolume(bReCalculate);

            m_originalVolume = m_polyhedron.volume();

        }

        template <typename Real>
        Real CFvmControlVolume<Real>::originalVolume()
        {

            return m_originalVolume;

        }

        template <typename Real>
        Real CFvmControlVolume<Real>::originalAvFaceDist(const Integer aFaceId)
        {

            return pow(m_originalVolume, 1.0 / 3.0) * 0.5;

        }

        template <typename Real>
        bool CFvmControlVolume<Real>::isClipped()
        {

            return m_clipped;

        }

        template <typename Real>
        void CFvmControlVolume<Real>::setClipped(bool clipped)
        {

            m_clipped = clipped;

        }

    }

}


