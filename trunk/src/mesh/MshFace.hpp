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

#include <vector>

#include "GeoArea.hpp"

namespace ENigMA
{

    namespace mesh
    {

        enum EFaceType 
        {
            FT_NONE = 0,
            FT_POINT,
            FT_LINE,
            FT_TRIANGLE,
            FT_QUADRILATERAL
        };

        template <typename Real>
        class CMshFace
        {
        private:

            std::vector<Integer> m_nodeIds;

            EFaceType m_faceType;

            bool m_hasPair;

            Integer m_pairFaceId;
            Integer m_elementId;

        public:
            CMshFace();
            CMshFace(EFaceType aFaceType);
            ~CMshFace();

            Integer nbNodeIds();

            void addNodeId(Integer aNodeId);
            Integer nodeId(Integer aNodeIndex);
            void setNodeId(const Integer aNodeIndex, const Integer aNodeId);

            void setPairFaceId(Integer aPairFace);
            Integer pairFaceId();

            void setHasPair(bool hasPair);
            bool hasPair();

            void setElementId(Integer anElementId);
            Integer elementId();

            void setFaceType(EFaceType aFaceType);
            EFaceType faceType();

            void reset();

        };

        template <typename Real>
        std::ostream &operator<< (std::ostream &output, CMshFace<Real> &aFace);

    }

}

#include "MshFace_Imp.hpp"
