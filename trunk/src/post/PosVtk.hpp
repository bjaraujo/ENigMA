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

#include "PdeField.hpp"

using namespace ENigMA::pde;

namespace ENigMA
{

    namespace post
    {

        template <typename Real>
        class CPosVtk
        {
        public:

            CPosVtk();
            ~CPosVtk();

            bool load(CPdeField<Real>& aField, std::string strFileName);
            bool save(CPdeField<Real>& aField, std::string strFileName);

        };

    }

}

#include "PosVtk_Imp.hpp"

