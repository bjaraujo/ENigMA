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

#include "AnaFunction.hpp"

using namespace ENigMA::analytical;

namespace ENigMA
{

    namespace pde
    {

        enum EBoundaryConditionType 
        {
            BT_GENERIC_FIXED_VALUE = 0,
            BT_GENERIC_ZERO_GRADIENT,
            BT_GENERIC_CYCLIC,
            BT_HEAT_CONVECTIVE
        };

        enum EBoundaryLocation 
        {
            BL_NODE = 0,
            BL_EDGE,
            BL_FACE
        };

        enum EConditionType 
        {
            CT_GENERIC_FIXED_VALUE = 0,
            CT_GENERIC_ZERO_GRADIENT,
            CT_HEAT_TRANSFER_COEFFICIENT,
            CT_HEAT_INFINITESIMAL_TEMPERATURE
        };

        template <typename Real>
        class CPdeBoundaryCondition
        {
        private:

            EBoundaryConditionType m_boundaryConditionType;
            EBoundaryLocation m_boundaryLocation;
            
            std::map<EConditionType, Real> m_condition;
            std::map<EConditionType, CAnaFunction<Real> > m_conditionFunc;

        public:

            CPdeBoundaryCondition(EBoundaryConditionType aBoundaryType);
            CPdeBoundaryCondition();
            ~CPdeBoundaryCondition();

            void setBoundaryConditionType(EBoundaryConditionType aBoundaryType);
            EBoundaryConditionType boundaryConditionType();

            void setLocation(EBoundaryLocation aBoundaryLocation);
            EBoundaryLocation location();

            void addCondition(EConditionType aConditionType, const Real aConditionValue);
            void addCondition(EConditionType aConditionType, CAnaFunction<Real>& aPropertyFunction);

            Real conditionValue(EConditionType aConditionType);
            Real conditionValue(EConditionType aConditionType, std::string aVariable1, Real& aValue1);

        };

    }

}

#include "PdeBoundaryCondition_Imp.hpp"

