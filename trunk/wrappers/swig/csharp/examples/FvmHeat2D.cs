using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Kitware.VTK;

namespace Demo
{
    class FvmHeat3D
    {

        static public void GenerateMesh(CMshBasicMesher aMesher, vtkRenderer aRenderer)
        {

            var aVertex1 = new CGeoCoordinate(+0.00, +0.00, -0.1);
            var aVertex2 = new CGeoCoordinate(+1.00, +0.00, -0.1);
            var aVertex3 = new CGeoCoordinate(+1.00, +1.00, -0.1);
            var aVertex4 = new CGeoCoordinate(+0.00, +1.00, -0.1);
            var aVertex5 = new CGeoCoordinate(+0.00, +0.00, +0.1);
            var aVertex6 = new CGeoCoordinate(+1.00, +0.00, +0.1);
            var aVertex7 = new CGeoCoordinate(+1.00, +1.00, +0.1);
            var aVertex8 = new CGeoCoordinate(+0.00, +1.00, +0.1);

            var aHexahedron = new CGeoHexahedron();

            aHexahedron.addVertex(aVertex1);
            aHexahedron.addVertex(aVertex2);
            aHexahedron.addVertex(aVertex3);
            aHexahedron.addVertex(aVertex4);
            aHexahedron.addVertex(aVertex5);
            aHexahedron.addVertex(aVertex6);
            aHexahedron.addVertex(aVertex7);
            aHexahedron.addVertex(aVertex8);

            aMesher.generate(aHexahedron, 40, 40, 1, false);

            aMesher.mesh().generateFaces(1E-3);

        }

        static public void DrawMesh(CMshBasicMesher aMesher, vtkUnstructuredGrid anUnstructuredGrid)
        {

            vtkPoints sPoints = vtkPoints.New();

            for (int i = 0; i < aMesher.mesh().nbNodes(); ++i)
            {

                int aNodeId = aMesher.mesh().nodeId(i);

                var aNode = aMesher.mesh().node(aNodeId);

                sPoints.InsertNextPoint(aNode.x(), aNode.y(), aNode.z());

            }

            anUnstructuredGrid.SetPoints(sPoints);

            for (int i = 0; i < aMesher.mesh().nbElements(); ++i)
            {

                int anElementId = aMesher.mesh().elementId(i);

                var anElement = aMesher.mesh().element(anElementId);

                if (anElement.elementType() == EElementType.ET_HEXAHEDRON)
                {

                    var aHexahedron = vtkHexahedron.New();
                    aHexahedron.GetPointIds().SetId(0, aMesher.mesh().nodeIndex(anElement.nodeId(0)));
                    aHexahedron.GetPointIds().SetId(1, aMesher.mesh().nodeIndex(anElement.nodeId(1)));
                    aHexahedron.GetPointIds().SetId(2, aMesher.mesh().nodeIndex(anElement.nodeId(2)));
                    aHexahedron.GetPointIds().SetId(3, aMesher.mesh().nodeIndex(anElement.nodeId(3)));
                    aHexahedron.GetPointIds().SetId(4, aMesher.mesh().nodeIndex(anElement.nodeId(4)));
                    aHexahedron.GetPointIds().SetId(5, aMesher.mesh().nodeIndex(anElement.nodeId(5)));
                    aHexahedron.GetPointIds().SetId(6, aMesher.mesh().nodeIndex(anElement.nodeId(6)));
                    aHexahedron.GetPointIds().SetId(7, aMesher.mesh().nodeIndex(anElement.nodeId(7)));

                    anUnstructuredGrid.InsertNextCell(aHexahedron.GetCellType(), aHexahedron.GetPointIds());

                }

            }

        }

        static public void Simulate(CMshBasicMesher aMesher, vtkUnstructuredGrid anUnstructuredGrid, vtkRenderer aRenderer)
        {

            var T = new CPdeField();

            aMesher.mesh().generateFaces(1E-12);

            aMesher.mesh().calculateFaceCentroid();
            aMesher.mesh().calculateElementCentroid();

            T.setMesh(aMesher.mesh());
            T.setDiscretMethod(EDiscretMethod.DM_FVM);
            T.setDiscretOrder(EDiscretOrder.DO_LINEAR);
            T.setDiscretLocation(EDiscretLocation.DL_ELEMENT_CENTER);
            T.setSimulationType(ESimulationType.ST_THERMAL);
            T.setNbDofs(1);

            // Set BC and initial conditions
            T.setSize(T.mesh().nbElements());

            for (int i = 0; i < T.mesh().nbFaces(); ++i)
            {

                int aFaceId = T.mesh().faceId(i);

                if (Math.Abs(T.mesh().faceCentroid(aFaceId).x() - 0.0) < 1E-6)
                {
                    var aFixedTemperature = new CPdeBoundaryCondition(EBoundaryConditionType.BT_GENERIC_FIXED_VALUE);
                    aFixedTemperature.addCondition(EConditionType.CT_GENERIC_FIXED_VALUE, 0.0);
                    T.addBCFace(aFaceId, aFixedTemperature);
                }

                if (Math.Abs(T.mesh().faceCentroid(aFaceId).x() - 1.0) < 1E-6)
                {
                    var aFixedTemperature = new CPdeBoundaryCondition(EBoundaryConditionType.BT_GENERIC_FIXED_VALUE);
                    aFixedTemperature.addCondition(EConditionType.CT_GENERIC_FIXED_VALUE, 0.0);
                    T.addBCFace(aFaceId, aFixedTemperature);
                }

                if (Math.Abs(T.mesh().faceCentroid(aFaceId).y() - 0.0) < 1E-6)
                {
                    var aFixedTemperature = new CPdeBoundaryCondition(EBoundaryConditionType.BT_GENERIC_FIXED_VALUE);
                    aFixedTemperature.addCondition(EConditionType.CT_GENERIC_FIXED_VALUE, 0.0);
                    T.addBCFace(aFaceId, aFixedTemperature);
                }

                if (Math.Abs(T.mesh().faceCentroid(aFaceId).y() - 1.0) < 1E-6)
                {
                    var aFixedTemperature = new CPdeBoundaryCondition(EBoundaryConditionType.BT_GENERIC_FIXED_VALUE);
                    aFixedTemperature.addCondition(EConditionType.CT_GENERIC_FIXED_VALUE, 1.0);
                    T.addBCFace(aFaceId, aFixedTemperature);
                }

            }

            // Steady conduction in a rectangle
            var aPdeEquation = new CPdeEquation(ENigMA.laplacian(T));

            aPdeEquation.solve(T);

            var sScalars = vtkFloatArray.New();

            sScalars.SetNumberOfTuples(T.mesh().nbElements());

            for (int i = 0; i < T.mesh().nbElements(); ++i)
            {

                int anElementId = T.mesh().elementId(i);

                sScalars.SetTuple1(anElementId, T.value(i));

            }

            anUnstructuredGrid.GetCellData().SetScalars(sScalars);
            
            var aGeometryFilter = vtkGeometryFilter.New();
            aGeometryFilter.SetInput(anUnstructuredGrid);
            aGeometryFilter.Update();

            var aCellDataToPointDataFilter1 = vtkCellDataToPointData.New();
            aCellDataToPointDataFilter1.SetInputConnection(aGeometryFilter.GetOutputPort());
            aCellDataToPointDataFilter1.Update();

            var aBoundedFilter = vtkBandedPolyDataContourFilter.New();
            aBoundedFilter.SetInput(aCellDataToPointDataFilter1.GetOutput());
            aBoundedFilter.GenerateValues(24, sScalars.GetRange()[0], sScalars.GetRange()[1]);

            var aLookupTable = vtkLookupTable.New();
            aLookupTable.SetNumberOfColors(256);
            aLookupTable.SetHueRange(0.667, 0.0);
            aLookupTable.Build();

            var aBandedMapper = vtkPolyDataMapper.New();
            aBandedMapper.SetInputConnection(aBoundedFilter.GetOutputPort());
            aBandedMapper.SetScalarModeToUsePointData();
            aBandedMapper.SetScalarRange(sScalars.GetRange()[0], sScalars.GetRange()[1]);
            aBandedMapper.SetLookupTable(aLookupTable);

            // add actor to the renderer
            var aMeshActor = vtkActor.New();
            aRenderer.AddActor(aMeshActor);

            aMeshActor.SetMapper(aBandedMapper);

            var aScalarBarActor = vtkScalarBarActor.New();
            aRenderer.AddActor(aScalarBarActor);

            aScalarBarActor.SetLookupTable(aLookupTable);
            aScalarBarActor.SetTitle("Temperature");
            aScalarBarActor.SetNumberOfLabels(6);
            aScalarBarActor.Modified();

        }

        static public void Show(RenderWindowControl aRenderControl)
        {

            // get a reference to the renderwindow of our renderWindowControl1
            var aRenderWindow = aRenderControl.RenderWindow;

            // get a reference to the renderer
            var aRenderer = aRenderWindow.GetRenderers().GetFirstRenderer();

            aRenderer.Clear();
            aRenderer.RemoveAllViewProps(); 

            // set background color
            aRenderer.SetBackground(0.2, 0.3, 0.4);

            var aMesher = new CMshBasicMesher();

            FvmHeat3D.GenerateMesh(aMesher, aRenderer);

            var anUnstructuredGrid = vtkUnstructuredGrid.New();

            FvmHeat3D.DrawMesh(aMesher, anUnstructuredGrid);

            FvmHeat3D.Simulate(aMesher, anUnstructuredGrid, aRenderer);

            aRenderer.ResetCamera();

            aRenderControl.Refresh();

        }

    }
}
