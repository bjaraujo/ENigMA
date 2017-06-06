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

#include <string>

#include "exprtk.hpp"

#include "CmnTypes.hpp"

namespace ENigMA
{

    namespace analytical
    {

        template <typename Real>
        class CAnaFunction
        {
        private:

            std::string m_expressionString;
            exprtk::symbol_table<Real> m_symbolTable;
            Real m_value;
            bool m_constant;

        public:

            CAnaFunction();
            CAnaFunction(Real aConstant);
            CAnaFunction(std::string anExpression);
            ~CAnaFunction();

            void set(Real aConstant);
            void set(std::string anExpression);

            void defineVariable(const std::string aVariable, Real& aValue);
            void undefineVariable(const std::string aVariable);
            void removeAllVariables();

            Real evaluate();
            Real evaluate(const std::string aVariable, Real& aValue);
            Real evaluate(const std::string aVariable1, Real& aValue1, const std::string aVariable2, Real& aValue2);

            Real bisection(const std::string strVar, Real lowerBnd, Real upperBnd, Integer& nIterations, const Integer nMaxIterations, const Real aTolerance);
            Real brent(const std::string strVar, Real lowerBnd, Real upperBnd, Integer& nIterations, const Integer nMaxIterations, const Real aTolerance);
            Real root(const std::string strVar, Real lowerBnd, Real upperBnd, Integer& nIterations, const Integer nMaxIterations, const Real aTolerance);

        };

    }

}

#include "AnaFunction_Imp.hpp"

