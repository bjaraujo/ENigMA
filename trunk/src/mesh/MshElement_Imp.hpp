// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

namespace ENigMA {
namespace mesh {
    template <typename Real>
    CMshElement<Real>::CMshElement()
        : m_elementType(ET_NONE)
        , m_thickness(1.0)
    {
    }

    template <typename Real>
    CMshElement<Real>::CMshElement(const CMshElement<Real>& anElement)
    {
        m_elementType = anElement.m_elementType;
        m_thickness = anElement.m_thickness;
        m_nodeIds = anElement.m_nodeIds;
        m_faceIds = anElement.m_faceIds;
    }

    template <typename Real>
    CMshElement<Real>::CMshElement(EElementType anElementType)
        : m_elementType(anElementType)
    {
    }

    template <typename Real>
    CMshElement<Real>::~CMshElement()
    {
    }

    template <typename Real>
    Integer CMshElement<Real>::nbNodeIds() const
    {
        return static_cast<Integer>(m_nodeIds.size());
    }

    template <typename Real>
    void CMshElement<Real>::addNodeId(const Integer aNodeId)
    {
        m_nodeIds.push_back(aNodeId);
    }

    template <typename Real>
    void CMshElement<Real>::setNodeId(const Integer aNodeIndex, const Integer aNodeId)
    {
        m_nodeIds[aNodeIndex] = aNodeId;
    }

    template <typename Real>
    Integer CMshElement<Real>::nodeId(const Integer aNodeIndex)
    {
        return m_nodeIds.at(aNodeIndex);
    }

    template <typename Real>
    Integer CMshElement<Real>::nbFaceIds() const
    {
        return static_cast<Integer>(m_faceIds.size());
    }

    template <typename Real>
    void CMshElement<Real>::addFaceId(const Integer aFaceId)
    {
        m_faceIds.push_back(aFaceId);
    }

    template <typename Real>
    void CMshElement<Real>::setFaceId(const Integer aFaceIndex, const Integer aFaceId)
    {
        m_faceIds[aFaceIndex] = aFaceId;
    }

    template <typename Real>
    Integer CMshElement<Real>::faceId(const Integer aFaceIndex)
    {
        return m_faceIds.at(aFaceIndex);
    }

    template <typename Real>
    void CMshElement<Real>::setElementType(EElementType anElementType)
    {
        m_elementType = anElementType;
    }

    template <typename Real>
    EElementType CMshElement<Real>::elementType()
    {
        return m_elementType;
    }

    template <typename Real>
    void CMshElement<Real>::reset()
    {
        m_elementType = ET_NONE;
        m_nodeIds.clear();
    }

    template <typename Real>
    void CMshElement<Real>::generateFaces(std::vector<ENigMA::mesh::CMshFace<Real>>& sFaces)
    {
        CMshFace<Real> aFace;

        switch (m_elementType) {
        case ET_NONE:
            break;
        case ET_BEAM:

            // face 1 (node)
            aFace.addNodeId(m_nodeIds[0]);
            aFace.setFaceType(FT_POINT);

            sFaces.push_back(aFace);
            aFace.reset();

            // face 2 (node)
            aFace.addNodeId(m_nodeIds[1]);
            aFace.setFaceType(FT_POINT);

            sFaces.push_back(aFace);
            aFace.reset();

            break;
        case ET_TRIANGLE:

            // face 1 (edge)
            aFace.addNodeId(m_nodeIds[0]);
            aFace.addNodeId(m_nodeIds[1]);
            aFace.setFaceType(FT_LINE);

            sFaces.push_back(aFace);
            aFace.reset();

            // face 2 (edge)
            aFace.addNodeId(m_nodeIds[1]);
            aFace.addNodeId(m_nodeIds[2]);
            aFace.setFaceType(FT_LINE);

            sFaces.push_back(aFace);
            aFace.reset();

            // face 3 (edge)
            aFace.addNodeId(m_nodeIds[2]);
            aFace.addNodeId(m_nodeIds[0]);
            aFace.setFaceType(FT_LINE);

            sFaces.push_back(aFace);
            aFace.reset();

            break;
        case ET_QUADRILATERAL:

            // face 1 (edge)
            aFace.addNodeId(m_nodeIds[0]);
            aFace.addNodeId(m_nodeIds[1]);
            aFace.setFaceType(FT_LINE);

            sFaces.push_back(aFace);
            aFace.reset();

            // face 2 (edge)
            aFace.addNodeId(m_nodeIds[1]);
            aFace.addNodeId(m_nodeIds[2]);
            aFace.setFaceType(FT_LINE);

            sFaces.push_back(aFace);
            aFace.reset();

            // face 3 (edge)
            aFace.addNodeId(m_nodeIds[2]);
            aFace.addNodeId(m_nodeIds[3]);
            aFace.setFaceType(FT_LINE);

            sFaces.push_back(aFace);
            aFace.reset();

            // face 4 (edge)
            aFace.addNodeId(m_nodeIds[3]);
            aFace.addNodeId(m_nodeIds[0]);
            aFace.setFaceType(FT_LINE);

            sFaces.push_back(aFace);
            aFace.reset();

            break;
        case ET_TETRAHEDRON:

            // face 1
            aFace.addNodeId(m_nodeIds[0]);
            aFace.addNodeId(m_nodeIds[1]);
            aFace.addNodeId(m_nodeIds[2]);
            aFace.setFaceType(FT_TRIANGLE);

            sFaces.push_back(aFace);
            aFace.reset();

            // face 2
            aFace.addNodeId(m_nodeIds[1]);
            aFace.addNodeId(m_nodeIds[3]);
            aFace.addNodeId(m_nodeIds[2]);
            aFace.setFaceType(FT_TRIANGLE);

            sFaces.push_back(aFace);
            aFace.reset();

            // face 3
            aFace.addNodeId(m_nodeIds[0]);
            aFace.addNodeId(m_nodeIds[3]);
            aFace.addNodeId(m_nodeIds[1]);
            aFace.setFaceType(FT_TRIANGLE);

            sFaces.push_back(aFace);
            aFace.reset();

            // face 4
            aFace.addNodeId(m_nodeIds[0]);
            aFace.addNodeId(m_nodeIds[2]);
            aFace.addNodeId(m_nodeIds[3]);
            aFace.setFaceType(FT_TRIANGLE);

            sFaces.push_back(aFace);
            aFace.reset();

            break;

        case ET_TRIANGULAR_PRISM:

            // face 1
            aFace.addNodeId(m_nodeIds[0]);
            aFace.addNodeId(m_nodeIds[2]);
            aFace.addNodeId(m_nodeIds[1]);
            aFace.setFaceType(FT_TRIANGLE);

            sFaces.push_back(aFace);
            aFace.reset();

            // face 2
            aFace.addNodeId(m_nodeIds[3]);
            aFace.addNodeId(m_nodeIds[4]);
            aFace.addNodeId(m_nodeIds[5]);
            aFace.setFaceType(FT_TRIANGLE);

            sFaces.push_back(aFace);
            aFace.reset();

            // face 3
            aFace.addNodeId(m_nodeIds[0]);
            aFace.addNodeId(m_nodeIds[3]);
            aFace.addNodeId(m_nodeIds[4]);
            aFace.addNodeId(m_nodeIds[1]);
            aFace.setFaceType(FT_QUADRILATERAL);

            sFaces.push_back(aFace);
            aFace.reset();

            // face 4
            aFace.addNodeId(m_nodeIds[2]);
            aFace.addNodeId(m_nodeIds[1]);
            aFace.addNodeId(m_nodeIds[4]);
            aFace.addNodeId(m_nodeIds[5]);
            aFace.setFaceType(FT_QUADRILATERAL);

            sFaces.push_back(aFace);
            aFace.reset();

            // face 5
            aFace.addNodeId(m_nodeIds[0]);
            aFace.addNodeId(m_nodeIds[2]);
            aFace.addNodeId(m_nodeIds[5]);
            aFace.addNodeId(m_nodeIds[3]);
            aFace.setFaceType(FT_QUADRILATERAL);

            sFaces.push_back(aFace);
            aFace.reset();

            break;

        case ET_HEXAHEDRON:

            // face 1
            aFace.addNodeId(m_nodeIds[0]);
            aFace.addNodeId(m_nodeIds[3]);
            aFace.addNodeId(m_nodeIds[2]);
            aFace.addNodeId(m_nodeIds[1]);
            aFace.setFaceType(FT_QUADRILATERAL);

            sFaces.push_back(aFace);
            aFace.reset();

            // face 2
            aFace.addNodeId(m_nodeIds[4]);
            aFace.addNodeId(m_nodeIds[5]);
            aFace.addNodeId(m_nodeIds[6]);
            aFace.addNodeId(m_nodeIds[7]);
            aFace.setFaceType(FT_QUADRILATERAL);

            sFaces.push_back(aFace);
            aFace.reset();

            // face 3
            aFace.addNodeId(m_nodeIds[4]);
            aFace.addNodeId(m_nodeIds[0]);
            aFace.addNodeId(m_nodeIds[1]);
            aFace.addNodeId(m_nodeIds[5]);
            aFace.setFaceType(FT_QUADRILATERAL);

            sFaces.push_back(aFace);
            aFace.reset();

            // face 4
            aFace.addNodeId(m_nodeIds[7]);
            aFace.addNodeId(m_nodeIds[6]);
            aFace.addNodeId(m_nodeIds[2]);
            aFace.addNodeId(m_nodeIds[3]);
            aFace.setFaceType(FT_QUADRILATERAL);

            sFaces.push_back(aFace);
            aFace.reset();

            // face 5
            aFace.addNodeId(m_nodeIds[5]);
            aFace.addNodeId(m_nodeIds[1]);
            aFace.addNodeId(m_nodeIds[2]);
            aFace.addNodeId(m_nodeIds[6]);
            aFace.setFaceType(FT_QUADRILATERAL);

            sFaces.push_back(aFace);
            aFace.reset();

            // face 6
            aFace.addNodeId(m_nodeIds[0]);
            aFace.addNodeId(m_nodeIds[4]);
            aFace.addNodeId(m_nodeIds[7]);
            aFace.addNodeId(m_nodeIds[3]);
            aFace.setFaceType(FT_QUADRILATERAL);

            sFaces.push_back(aFace);
            aFace.reset();

            break;

        case ET_POLYHEDRON:

            break;
        }

        m_faceIds.clear();
    }

    template <typename Real>
    void CMshElement<Real>::invert()
    {
        Integer idAux;

        switch (m_elementType) {
        case ET_NONE:
            // TODO:
            break;
        case ET_BEAM:
            idAux = m_nodeIds[0];
            m_nodeIds[0] = m_nodeIds[1];
            m_nodeIds[1] = idAux;
            break;
        case ET_TRIANGLE:
            idAux = m_nodeIds[1];
            m_nodeIds[1] = m_nodeIds[2];
            m_nodeIds[2] = idAux;
            break;
        case ET_QUADRILATERAL:
            idAux = m_nodeIds[1];
            m_nodeIds[1] = m_nodeIds[3];
            m_nodeIds[3] = idAux;
            break;
        case ET_TETRAHEDRON:
            idAux = m_nodeIds[1];
            m_nodeIds[1] = m_nodeIds[2];
            m_nodeIds[2] = idAux;
            break;
        case ET_TRIANGULAR_PRISM:
            // TODO:
            break;
        case ET_HEXAHEDRON:
            // TODO:
            break;
        case ET_POLYHEDRON:
            // TODO:
            break;
        }
    }

    template <typename Real>
    void CMshElement<Real>::setThickness(Real aThickness)
    {
        m_thickness = aThickness;
    }

    template <typename Real>
    Real CMshElement<Real>::thickness()
    {
        return m_thickness;
    }

    template <typename Real>
    std::ostream& operator<<(std::ostream& output, CMshElement<Real>& anElement)
    {
        for (Integer i = 0; i < anElement.nbNodeIds(); ++i)
            output << anElement.nodeId(i) << " ";

        return output;
    }
}
}
