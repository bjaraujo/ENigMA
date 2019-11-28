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

using namespace ENigMA::pde;

namespace ENigMA {
namespace post {
    template <typename Real>
    class CPosQuickMesh {
    public:
        CPosQuickMesh();
        virtual ~CPosQuickMesh();

        bool load(CPdeField<Real>& aField, const std::string& strFileName);
        bool save(CPdeField<Real>& aField, const std::string& strFileName);
    };
}
}

#include "PosQuickMesh_Imp.hpp"
