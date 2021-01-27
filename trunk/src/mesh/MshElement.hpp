// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "MshFace.hpp"

namespace ENigMA {
namespace mesh {
    enum EElementType {
        ET_NONE = 0,
        ET_BEAM,
        ET_TRIANGLE,
        ET_QUADRILATERAL,
        ET_TETRAHEDRON,
        ET_TRIANGULAR_PRISM,
        ET_HEXAHEDRON,
        ET_POLYHEDRON
    };

    template <typename Real>
    class CMshElement {
    private:
        EElementType m_elementType;
        Real m_thickness;

        std::vector<Integer> m_nodeIds;
        std::vector<Integer> m_faceIds;

    public:
        CMshElement();
        CMshElement(const CMshElement<Real>& anElement);
        explicit CMshElement(EElementType anElementType);
        virtual ~CMshElement();

        Integer nbNodeIds() const;

        void addNodeId(const Integer aNodeId);
        Integer nodeId(const Integer aNodeIndex);
        void setNodeId(const Integer aNodeIndex, const Integer aNodeId);

        Integer nbFaceIds() const;

        void addFaceId(const Integer aFaceId);
        Integer faceId(const Integer aFaceIndex);
        void setFaceId(const Integer aFaceIndex, const Integer aFaceId);

        void setElementType(EElementType anElementType);
        EElementType elementType();

        void setThickness(Real aThickness);
        Real thickness();

        void reset();

        void generateFaces(std::vector<ENigMA::mesh::CMshFace<Real>>& sFaces);

        void invert();
    };

    template <typename Real>
    std::ostream& operator<<(std::ostream& output, ENigMA::mesh::CMshElement<Real>& anElement);
}
}

#include "MshElement_Imp.hpp"
