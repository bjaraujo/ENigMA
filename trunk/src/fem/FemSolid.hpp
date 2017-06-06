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

using namespace ENigMA::geometry;

namespace ENigMA
{

    namespace fem
    {

        template <typename Real>
        class CFemSolid
        {

            virtual void setSourceOnNode(const Integer aNodeIndex, const Real aValue) = 0;
            virtual void setSourceOnEdge(const Integer anEdgeIndex, const Real aValue) = 0;
            virtual void setSourceOnFace(const Integer anFaceIndex, const Real aValue) = 0;

        };

    }

}
