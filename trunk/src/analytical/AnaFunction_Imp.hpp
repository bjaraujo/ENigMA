// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

namespace ENigMA
{
    namespace analytical
    {
        template <typename Real>
        CAnaFunction<Real>::CAnaFunction()
            : m_constant(false)
        {
        }

        template <typename Real>
        CAnaFunction<Real>::CAnaFunction(Real aConstant) { this->set(aConstant); }

        template <typename Real>
        CAnaFunction<Real>::CAnaFunction(const std::string& anExpression) { this->set(anExpression); }

        template <typename Real>
        CAnaFunction<Real>::~CAnaFunction() { }

        template <typename Real>
        void CAnaFunction<Real>::set(Real aConstant)
        {
            m_value = aConstant;
            m_constant = true;
        }

        template <typename Real>
        void CAnaFunction<Real>::set(const std::string& anExpression)
        {
            m_expressionString = anExpression;
            m_constant = false;
        }

        template <typename Real>
        void CAnaFunction<Real>::defineVariable(const std::string& aVariable, Real& aValue) { m_symbolTable.add_variable(aVariable, aValue); }

        template <typename Real>
        void CAnaFunction<Real>::undefineVariable(const std::string& aVariable) { m_symbolTable.remove_variable(aVariable); }

        template <typename Real>
        void CAnaFunction<Real>::removeAllVariables() { m_symbolTable.clear(); }

        template <typename Real>
        Real CAnaFunction<Real>::evaluate()
        {
            if (m_constant)
                return m_value;

            exprtk::expression<Real> expression;
            exprtk::parser<Real> parser;

            m_symbolTable.add_constants();

            expression.register_symbol_table(m_symbolTable);
            parser.compile(m_expressionString, expression);

            return expression.value();
        }

        template <typename Real>
        Real CAnaFunction<Real>::evaluate(const std::string& aVariable, Real& aValue)
        {
            exprtk::expression<Real> expression;
            exprtk::parser<Real> parser;

            m_symbolTable.clear();

            m_symbolTable.add_variable(aVariable, aValue);
            m_symbolTable.add_constants();

            expression.register_symbol_table(m_symbolTable);
            parser.compile(m_expressionString, expression);

            return expression.value();
        }

        template <typename Real>
        Real CAnaFunction<Real>::evaluate(const std::string& aVariable1, Real& aValue1, const std::string& aVariable2, Real& aValue2)
        {
            exprtk::expression<Real> expression;
            exprtk::parser<Real> parser;

            m_symbolTable.clear();

            m_symbolTable.add_variable(aVariable1, aValue1);
            m_symbolTable.add_variable(aVariable2, aValue2);
            m_symbolTable.add_constants();

            expression.register_symbol_table(m_symbolTable);
            parser.compile(m_expressionString, expression);

            return expression.value();
        }

        template <typename Real>
        Real CAnaFunction<Real>::bisection(const std::string& strVar, Real lowerBnd, Real upperBnd, Integer& nIterations, const Integer nMaxIterations, const Real aTolerance)
        {
            Real a = lowerBnd;
            Real b = upperBnd;

            Integer i = 0;

            // Start loop
            while (fabs(b - a) > aTolerance)
            {
                // calculate midpoint
                Real s = (a + b) * 0.5;

                // find f(midpoint)
                if (evaluate(strVar, a) * evaluate(strVar, s) > 0)
                {
                    // throw away left half
                    a = s;
                }
                else
                {
                    // throw away right half
                    b = s;
                }

                i++;

                if (i > nMaxIterations)
                    break;
            };

            nIterations = i;

            return (a + b) * 0.5;
        }

        template <typename Real>
        Real CAnaFunction<Real>::brent(const std::string& strVar, Real lowerBnd, Real upperBnd, Integer& nIterations, const Integer nMaxIterations, const Real aTolerance)
        {
            // http://en.wikipedia.org/wiki/Brent's_method

            Real a = lowerBnd;
            Real b = upperBnd;
            Real c = 0.0;
            Real d = std::numeric_limits<Real>::max();

            Real fa = evaluate(strVar, a);
            Real fb = evaluate(strVar, b);

            Real fc = 0.0;

            // if f(a) f(b) >= 0 then error-exit
            if (fa * fb >= 0)
            {
                if (fa < fb)
                    return a;
                else
                    return b;
            }

            // if |f(a)| < |f(b)| then swap (a,b) end if
            if (fabs(fa) < fabs(fb))
            {
                Real tmp = a;
                a = b;
                b = tmp;
                tmp = fa;
                fa = fb;
                fb = tmp;
            }

            c = a;
            fc = fa;

            bool mflag = true;

            Integer i = 0;

            while ((fb != 0) && (fabs(a - b) > aTolerance))
            {
                Real s;

                if ((fa != fc) && (fb != fc))
                    // Inverse quadratic interpolation
                    s = a * fb * fc / (fa - fb) / (fa - fc) + b * fa * fc / (fb - fa) / (fb - fc) + c * fa * fb / (fc - fa) / (fc - fb);
                else
                    // Secant Rule
                    s = b - fb * (b - a) / (fb - fa);

                Real tmp2 = (3 * a + b) / 4;

                if ((!(((s > tmp2) && (s < b)) || ((s < tmp2) && (s > b)))) || (mflag && (fabs(s - b) >= (fabs(b - c) / 2))) || (!mflag && (fabs(s - b) >= (fabs(c - d) / 2))))
                {
                    s = (a + b) * 0.5;
                    mflag = true;
                }
                else
                {
                    if ((mflag && (fabs(b - c) < aTolerance)) || (!mflag && (fabs(c - d) < aTolerance)))
                    {
                        s = (a + b) / 2;
                        mflag = true;
                    }
                    else
                        mflag = false;
                }

                Real fs = evaluate(strVar, s);
                d = c;
                c = b;
                fc = fb;

                if (fa * fs < 0)
                {
                    b = s;
                    fb = fs;
                }
                else
                {
                    a = s;
                    fa = fs;
                }

                // if |f(a)| < |f(b)| then swap (a,b) end if
                if (fabs(fa) < fabs(fb))
                {
                    Real tmp = a;
                    a = b;
                    b = tmp;
                    tmp = fa;
                    fa = fb;
                    fb = tmp;
                }

                i++;

                if (i > nMaxIterations)
                    break;
            }

            nIterations = i;

            return b;
        }

        template <typename Real>
        Real CAnaFunction<Real>::root(const std::string& strVar, Real lowerBnd, Real upperBnd, Integer& nIterations, const Integer nMaxIterations, const Real aTolerance) { return bisection(strVar, lowerBnd, upperBnd, nIterations, nMaxIterations, aTolerance); }
    }
}
