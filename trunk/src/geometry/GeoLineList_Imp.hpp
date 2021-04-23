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
    namespace geometry
    {
        template <typename Real>
        CGeoLineList<Real>::CGeoLineList()
        {
        }

        template <typename Real>
        CGeoLineList<Real>::~CGeoLineList()
        {
            reset();
        }

        template <typename Real>
        void CGeoLineList<Real>::reset()
        {
            m_lines.clear();
        }

        template <typename Real>
        Integer CGeoLineList<Real>::nbLines() const
        {
            return static_cast<Integer>(m_lines.size());
        }

        template <typename Real>
        void CGeoLineList<Real>::addLine(const CGeoLine<Real>& aLine)
        {
            m_lines.push_back(aLine);
        }

        template <typename Real>
        CGeoLine<Real>& CGeoLineList<Real>::line(const Integer aLineIndex)
        {
            return m_lines[aLineIndex];
        }

        template <typename Real>
        void CGeoLineList<Real>::calculateLength(bool bReCalculate)
        {
            if (!this->m_bLength || bReCalculate)
            {
                CGeoLength<Real>::m_length = 0;

                for (Integer i = 0; i < static_cast<Integer>(m_lines.size()); ++i)
                {
                    m_lines[i].calculateLength();
                    CGeoLength<Real>::m_length += m_lines[i].length();
                }

                this->m_bLength = true;
            }
        }

        template <typename Real>
        void CGeoLineList<Real>::calculateBoundingBox(bool bReCalculate)
        {
            // TODO:
        }

        template <typename Real>
        void CGeoLineList<Real>::sort(const Real aTolerance)
        {
            for (Integer i = 0; i < static_cast<Integer>(m_lines.size()); ++i)
            {
                for (Integer j = i + 1; j < static_cast<Integer>(m_lines.size()); ++j)
                {
                    // Start point
                    if ((m_lines[j].startPoint() - m_lines[i].endPoint()).norm() <= aTolerance)
                    {
                        if (i + 1 != j)
                            std::swap(m_lines[i + 1], m_lines[j]);
                    }

                    // End point
                    if ((m_lines[j].endPoint() - m_lines[i].endPoint()).norm() <= aTolerance)
                    {
                        m_lines[j].invert();

                        if (i + 1 != j)
                            std::swap(m_lines[i + 1], m_lines[j]);
                    }
                }
            }
        }

        template <typename Real>
        void CGeoLineList<Real>::invert()
        {
            std::reverse(m_lines.begin(), m_lines.end());

            for (Integer i = 0; i < static_cast<Integer>(m_lines.size()); ++i)
                m_lines[i].invert();
        }
    }
}
