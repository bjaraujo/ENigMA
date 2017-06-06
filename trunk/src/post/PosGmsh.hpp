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

namespace ENigMA
{

    namespace post
    {

        template <typename Real>
        class CPosGmsh
        {
        private:

            bool findSection(std::ifstream& fileStream, std::string sectionName);

        public:

            CPosGmsh();
            ~CPosGmsh();

            bool load(ENigMA::pde::CPdeField<Real>& aField, const std::string strFileName);
            bool save(ENigMA::pde::CPdeField<Real>& aField, const std::string strFileName, const std::string strViewName);

        };

    }

}

#include "PosGmsh_Imp.hpp"

