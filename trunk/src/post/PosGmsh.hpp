// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
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
            bool findSection(std::ifstream& fileStream, const std::string& sectionName);

        public:
            CPosGmsh();
            virtual ~CPosGmsh();

            bool load(ENigMA::pde::CPdeField<Real>& aField, const std::string& strFileName);
            bool save(ENigMA::pde::CPdeField<Real>& aField, const std::string& strFileName, const std::string& strViewName);
        };
    }
}

#include "PosGmsh_Imp.hpp"
