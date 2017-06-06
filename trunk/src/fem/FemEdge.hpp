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

namespace ENigMA
{

    namespace fem
    {

        template <typename Real>
        class CFemEdge
        {
        protected:

            Real m_sectionArea;
            Real m_perimeter;
                
        public:

            CFemEdge();
            ~CFemEdge();

            Real& sectionArea();
            Real& perimeter();

            virtual void setSourceOnNode(const Integer aNodeIndex, const Real aValue) = 0;

        };

    }

}

#include "FemEdge_Imp.hpp"
