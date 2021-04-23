// *****************************************************************************
// <ProjectName> ENigMA </ProjectName>
// <Description> Extended Numerical Multiphysics Analysis </Description>
// <HeadURL> $HeadURL$ </HeadURL>
// <LastChangedDate> $LastChangedDate$ </LastChangedDate>
// <LastChangedRevision> $LastChangedRevision$ </LastChangedRevision>
// <Author> Billy Araujo </Author>
// *****************************************************************************

#pragma once

#include "MatMaterial.hpp"
#include "MshMesh.hpp"
#include "PdeBoundaryCondition.hpp"

namespace ENigMA
{
    namespace pde
    {
        enum EDiscretMethod
        {
            DM_NONE = 0,
            DM_FDM,
            DM_BEM,
            DM_FEM,
            DM_FVM,
            DM_LBM,
            DM_SPH
        };

        enum EDiscretOrder
        {
            DO_LINEAR = 1,
            DO_QUADRATIC = 2,
            DO_CUBIC = 3
        };

        enum EDiscretLocation
        {
            DL_NODE = 1,
            DL_ELEMENT_CENTER
        };

        enum ESimulationType
        {
            ST_GENERIC = 0,
            ST_THERMAL,
            ST_FLOW,
            ST_STRUCTURAL
        };

        template <typename Real>
        class CPdeField
        {
        private:
            ENigMA::mesh::CMshMesh<Real> m_mesh;

            ENigMA::material::CMatMaterial<Real> m_material;
            std::map<Integer, ENigMA::material::CMatMaterial<Real>> m_materials;

            std::map<Integer, CPdeBoundaryCondition<Real>> m_bcNode;
            std::map<Integer, CPdeBoundaryCondition<Real>> m_bcFace;
            std::map<std::pair<Integer, Integer>, CPdeBoundaryCondition<Real>> m_bcElement;

            Integer m_nbDofs;

            EDiscretMethod m_discretMethod;
            EDiscretOrder m_discretOrder;
            EDiscretLocation m_discretLocation;
            ESimulationType m_simulationType;

        public:
            CPdeField();
            virtual ~CPdeField();

            Eigen::Matrix<Real, Eigen::Dynamic, 1> u;

            std::map<Integer, Real> uFixed;
            std::map<Integer, Real> uSource;

            void init();

            void setMesh(const ENigMA::mesh::CMshMesh<Real>& aMesh);
            ENigMA::mesh::CMshMesh<Real>& mesh();

            void setNbDofs(const Integer nbDofs);
            Integer nbDofs() const;

            void setDiscretMethod(EDiscretMethod aDiscretMethod);
            EDiscretMethod discretMethod();

            void setDiscretOrder(EDiscretOrder aDiscretOrder);
            EDiscretOrder discretOrder();

            void setDiscretLocation(EDiscretLocation aDiscretLocation);
            EDiscretLocation discretLocation();

            void setSimulationType(ESimulationType aSimulationType);
            ESimulationType simulationType();

            void setSize(const Integer aSize);

            void setValue(const Integer anIndex, Real aValue);
            Real value(const Integer anIndex);

            void setFixedValue(const Integer anIndex, Real aValue);
            void setFixedValue(const Integer anIndex, const Integer aDof, Real aValue);

            void setSource(const Integer anIndex, Real aValue);
            void setSource(const Integer anIndex, const Integer aDof, Real aValue);

            void setMaterial(ENigMA::material::CMatMaterial<Real>& aMaterial);
            ENigMA::material::CMatMaterial<Real>& material();
            void addMaterial(const Integer anElementIndex, ENigMA::material::CMatMaterial<Real> aMatMaterial);
            bool hasMaterial(const Integer anElementIndex);
            ENigMA::material::CMatMaterial<Real>& material(const Integer anElementIndex);

            void addBCNode(const Integer aNodeId, const ENigMA::pde::CPdeBoundaryCondition<Real>& aBoundaryCondition);
            bool nodeHasBC(const Integer aNodeId);
            ENigMA::pde::CPdeBoundaryCondition<Real>& nodeBC(const Integer aNodeId);

            void addBCFace(const Integer aFaceId, const ENigMA::pde::CPdeBoundaryCondition<Real>& aBoundaryCondition);
            bool faceHasBC(const Integer aFaceId);
            ENigMA::pde::CPdeBoundaryCondition<Real>& faceBC(const Integer aFaceId);

            void addBCElement(const Integer anElementId, const Integer anIndex, const ENigMA::pde::CPdeBoundaryCondition<Real>& aBoundaryCondition);
            bool elementHasBC(const Integer anElementId, const Integer anIndex);
            ENigMA::pde::CPdeBoundaryCondition<Real>& elementBC(const Integer anElementId, const Integer anIndex);
        };
    }
}

#include "PdeField_Imp.hpp"
