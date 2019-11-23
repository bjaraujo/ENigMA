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

namespace pde {

    template <typename Real>
    CPdeField<Real>::CPdeField()
        : m_nbDofs(0),
        m_discretMethod(DM_NONE),
        m_discretOrder(DO_LINEAR),
        m_discretLocation(DL_NODE),
        m_simulationType(ST_GENERIC)
    {
    }

    template <typename Real>
    CPdeField<Real>::~CPdeField()
    {
    }

    template <typename Real>
    void CPdeField<Real>::init()
    {
        if (m_discretLocation == DL_NODE) {
            u.resize(m_mesh.nbNodes() * m_nbDofs);
        } else if (m_discretLocation == DL_ELEMENT_CENTER) {
            u.resize(m_mesh.nbElements() * m_nbDofs);
        }
    }

    template <typename Real>
    void CPdeField<Real>::setMesh(const ENigMA::mesh::CMshMesh<Real>& aMesh)
    {

        m_mesh = aMesh;
    }

    template <typename Real>
    ENigMA::mesh::CMshMesh<Real>& CPdeField<Real>::mesh()
    {

        return m_mesh;
    }

    template <typename Real>
    void CPdeField<Real>::setMaterial(ENigMA::material::CMatMaterial<Real>& aMaterial)
    {

        m_material = aMaterial;
    }

    template <typename Real>
    ENigMA::material::CMatMaterial<Real>& CPdeField<Real>::material()
    {

        return m_material;
    }

    template <typename Real>
    void CPdeField<Real>::setNbDofs(const Integer nbDofs)
    {

        m_nbDofs = nbDofs;
    }

    template <typename Real>
    Integer CPdeField<Real>::nbDofs() const
    {

        return m_nbDofs;
    }

    template <typename Real>
    void CPdeField<Real>::setDiscretMethod(EDiscretMethod aDiscretMethod)
    {

        m_discretMethod = aDiscretMethod;
    }

    template <typename Real>
    EDiscretMethod CPdeField<Real>::discretMethod()
    {

        return m_discretMethod;
    }

    template <typename Real>
    void CPdeField<Real>::setDiscretOrder(EDiscretOrder aDiscretOrder)
    {

        m_discretOrder = aDiscretOrder;
    }

    template <typename Real>
    EDiscretOrder CPdeField<Real>::discretOrder()
    {

        return m_discretOrder;
    }

    template <typename Real>
    void CPdeField<Real>::setDiscretLocation(EDiscretLocation aDiscretLocation)
    {

        m_discretLocation = aDiscretLocation;
    }

    template <typename Real>
    EDiscretLocation CPdeField<Real>::discretLocation()
    {

        return m_discretLocation;
    }

    template <typename Real>
    void CPdeField<Real>::setSimulationType(ESimulationType aSimulationType)
    {

        m_simulationType = aSimulationType;
    }

    template <typename Real>
    ESimulationType CPdeField<Real>::simulationType()
    {

        return m_simulationType;
    }

    template <typename Real>
    void CPdeField<Real>::setSize(const Integer aSize)
    {

        u.resize(aSize);
    }

    template <typename Real>
    void CPdeField<Real>::setValue(const Integer anIndex, Real aValue)
    {

        u(anIndex) = aValue;
    }

    template <typename Real>
    Real CPdeField<Real>::value(const Integer anIndex)
    {

        return u(anIndex);
    }

    template <typename Real>
    void CPdeField<Real>::setFixedValue(const Integer anIndex, Real aValue)
    {

        uFixed[anIndex] = aValue;
    }

    template <typename Real>
    void CPdeField<Real>::setFixedValue(const Integer anIndex, const Integer aDof, Real aValue)
    {

        uFixed[anIndex * m_nbDofs + aDof] = aValue;
    }

    template <typename Real>
    void CPdeField<Real>::setSource(const Integer anIndex, Real aValue)
    {

        uSource[anIndex] = aValue;
    }

    template <typename Real>
    void CPdeField<Real>::setSource(const Integer anIndex, const Integer aDof, Real aValue)
    {

        uSource[anIndex * m_nbDofs + aDof] = aValue;
    }

    template <typename Real>
    void CPdeField<Real>::addMaterial(const Integer anElementIndex, ENigMA::material::CMatMaterial<Real> aMatMaterial)
    {

        m_materials[anElementIndex] = aMatMaterial;
    }

    template <typename Real>
    bool CPdeField<Real>::hasMaterial(const Integer anElementIndex)
    {

        if (m_materials.find(anElementIndex) != m_materials.end())
            return true;
        else
            return false;
    }

    template <typename Real>
    ENigMA::material::CMatMaterial<Real>& CPdeField<Real>::material(const Integer anElementIndex)
    {

        return m_materials.at(anElementIndex);
    }

    // Boundary conditions

    template <typename Real>
    void CPdeField<Real>::addBCNode(const Integer aNodeId, const ENigMA::pde::CPdeBoundaryCondition<Real>& aBoundaryCondition)
    {

        m_bcNode[aNodeId] = aBoundaryCondition;
    }

    template <typename Real>
    bool CPdeField<Real>::nodeHasBC(const Integer aNodeId)
    {

        if (m_bcNode.find(aNodeId) != m_bcNode.end())
            return true;
        else
            return false;
    }

    template <typename Real>
    ENigMA::pde::CPdeBoundaryCondition<Real>& CPdeField<Real>::nodeBC(const Integer aNodeId)
    {

        return m_bcNode.at(aNodeId);
    }

    template <typename Real>
    void CPdeField<Real>::addBCFace(const Integer aFaceId, const ENigMA::pde::CPdeBoundaryCondition<Real>& aBoundaryCondition)
    {

        m_bcFace[aFaceId] = aBoundaryCondition;
    }

    template <typename Real>
    bool CPdeField<Real>::faceHasBC(const Integer aFaceId)
    {

        if (m_bcFace.find(aFaceId) != m_bcFace.end())
            return true;
        else
            return false;
    }

    template <typename Real>
    ENigMA::pde::CPdeBoundaryCondition<Real>& CPdeField<Real>::faceBC(const Integer aFaceId)
    {

        return m_bcFace.at(aFaceId);
    }

    template <typename Real>
    void CPdeField<Real>::addBCElement(const Integer anElementId, const Integer anIndex, const ENigMA::pde::CPdeBoundaryCondition<Real>& aBoundaryCondition)
    {

        m_bcElement[std::make_pair(anElementId, anIndex)] = aBoundaryCondition;
    }

    template <typename Real>
    bool CPdeField<Real>::elementHasBC(const Integer anElementId, const Integer anIndex)
    {

        if (m_bcElement.find(std::make_pair(anElementId, anIndex)) != m_bcElement.end())
            return true;
        else
            return false;
    }

    template <typename Real>
    ENigMA::pde::CPdeBoundaryCondition<Real>& CPdeField<Real>::elementBC(const Integer anElementId, const Integer anIndex)
    {

        return m_bcElement[std::make_pair(anElementId, anIndex)];
    }
}
}
