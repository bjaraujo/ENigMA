// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

// <ProjectName> $ProjectName$ </ProjectName>

#pragma once

namespace ENigMA {

namespace geometry {

    template <typename Real>
    CGeoPolyline<Real>::CGeoPolyline()
    {
    }

    template <typename Real>
    CGeoPolyline<Real>::CGeoPolyline(CGeoLineList<Real>& aLineList, bool sort, const Real aTolerance)
    {

        set(aLineList, sort, aTolerance);
    }

    template <typename Real>
    CGeoPolyline<Real>::~CGeoPolyline()
    {
    }

    template <typename Real>
    void CGeoPolyline<Real>::set(CGeoLineList<Real>& aLineList, bool sort, const Real aTolerance)
    {

        CGeoVertexList<Real>::reset();

        if (sort)
            aLineList.sort(aTolerance);

        for (Integer i = 0; i < aLineList.nbLines(); ++i) {

            aLineList.line(i).calculateLength();

            if (aLineList.line(i).length() > 0) {

                CGeoCoordinate<Real> aVertex;

                aVertex = aLineList.line(i).startPoint();
                CGeoVertexList<Real>::addVertex(aVertex);

                aVertex = aLineList.line(i).endPoint();
                CGeoVertexList<Real>::addVertex(aVertex);
            }
        }

        for (Integer i = 1; i < CGeoVertexList<Real>::nbVertices(); ++i) {

            Real d = (CGeoVertexList<Real>::vertex(i) - CGeoVertexList<Real>::vertex(i - 1)).norm();

            if (d <= aTolerance) {
                // fuse vertices
                CGeoVertexList<Real>::vertex(i) = CGeoVertexList<Real>::vertex(i - 1);
            }
        }

        CGeoVertexList<Real>::removeDuplicates();
    }

    template <typename Real>
    Integer CGeoPolyline<Real>::nbLines()
    {

        if (CGeoVertexList<Real>::nbVertices() == 0)
            return 0;

        return CGeoVertexList<Real>::nbVertices() - 1;
    }

    template <typename Real>
    CGeoLine<Real> CGeoPolyline<Real>::line(Integer aLineIndex)
    {

        CGeoLine<Real> aLine;

        CGeoCoordinate<Real> p1 = CGeoVertexList<Real>::vertex((aLineIndex + 0) % CGeoVertexList<Real>::nbVertices());
        CGeoCoordinate<Real> p2 = CGeoVertexList<Real>::vertex((aLineIndex + 1) % CGeoVertexList<Real>::nbVertices());

        aLine.setStartPoint(p1);
        aLine.setEndPoint(p2);

        return aLine;
    }

    template <typename Real>
    bool CGeoPolyline<Real>::isClosed(const Real aTolerance)
    {

        if (CGeoVertexList<Real>::nbVertices() < 2)
            return false;

        return (CGeoVertexList<Real>::vertex(CGeoVertexList<Real>::nbVertices() - 1) - this->m_vertices[0]).norm() <= aTolerance;
    }

    template <typename Real>
    void CGeoPolyline<Real>::close(const Real aTolerance)
    {

        if (CGeoVertexList<Real>::nbVertices() > 2) {

            Real d = (CGeoVertexList<Real>::vertex(CGeoVertexList<Real>::nbVertices() - 1) - this->m_vertices[0]).norm();

            if (d > aTolerance)
                CGeoVertexList<Real>::addVertex(this->m_vertices[0]);
        }
    }

    template <typename Real>
    void CGeoPolyline<Real>::calculateLength(bool bReCalculate)
    {

        if (!this->m_bLength || bReCalculate) {

            CGeoLength<Real>::m_length = 0.0;

            for (Integer i = 0; i < CGeoVertexList<Real>::nbVertices() - 1; ++i) {
                CGeoLength<Real>::m_length += (CGeoVertexList<Real>::vertex(i + 1) - CGeoVertexList<Real>::vertex(i)).norm();
            }

            this->m_bLength = true;
        }
    }

    template <typename Real>
    void CGeoPolyline<Real>::calculateBoundingBox(bool bReCalculate)
    {

        // TODO:
    }

    template <typename Real>
    CGeoLineList<Real> CGeoPolyline<Real>::clip(CGeoPlane<Real>& aPlane)
    {

        CGeoLineList<Real> aLineList;

        for (Integer i = 0; i < CGeoVertexList<Real>::nbVertices() - 1; ++i) {

            CGeoLine<Real> aLine;

            CGeoCoordinate<Real> p1, p2;

            p1 = CGeoVertexList<Real>::vertex(i);
            p2 = CGeoVertexList<Real>::vertex(i + 1);

            aLine.setStartPoint(p1);
            aLine.setEndPoint(p2);

            CGeoLine<Real> aClippedLine = aLine.clip(aPlane);

            aClippedLine.calculateLength();

            if (aClippedLine.length() > 0) {
                aLineList.addLine(aClippedLine);
            }
        }

        return aLineList;
    }

    template <typename Real>
    bool CGeoPolyline<Real>::intersects(CGeoPlane<Real>& aPlane)
    {

        try {

            CGeoCoordinate<Real> aPoint;
            CGeoIntersectionType anIntersectionType;

            for (Integer i = 0; i < CGeoVertexList<Real>::nbVertices() - 1; ++i) {

                CGeoLine<Real> aLine;

                CGeoCoordinate<Real> p1, p2;

                p1 = CGeoVertexList<Real>::vertex(i);
                p2 = CGeoVertexList<Real>::vertex(i + 1);

                aLine.setStartPoint(p1);
                aLine.setEndPoint(p2);

                if (aLine.intersects(aPlane, aPoint, anIntersectionType))
                    return true;
            }

        } catch (const std::exception& e) {
            std::cout << "Error: std exception: " << e.what() << " in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
        } catch (...) {
            std::cout << "Error: unknown exception in function: " << ENIGMA_CURRENT_FUNCTION << std::endl;
        }

        return false;
    }

    template <typename Real>
    void CGeoPolyline<Real>::invert()
    {

        CGeoVertexList<Real>::invert();
    }
}
}
