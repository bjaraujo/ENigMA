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
    class CPosGnuplot {
    public:
        CPosGnuplot();
        virtual ~CPosGnuplot();

        bool save(CPdeField<Real>& aField, std::string strFileName);
    };
}
}

#include "PosGnuplot_Imp.hpp"
