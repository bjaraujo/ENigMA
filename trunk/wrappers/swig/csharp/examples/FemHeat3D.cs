using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Kitware.VTK;

namespace Demo
{
    class FemHeat3D
    {

        static public void GenerateMesh(CMshBasicMesherDouble aMesher, vtkRenderer aRenderer)
        {

            var aVertex1 = new CGeoCoordinateDouble(+0.00, +0.00, -0.05);
            var aVertex2 = new CGeoCoordinateDouble(+1.00, +0.00, -0.05);
            var aVertex3 = new CGeoCoordinateDouble(+1.00, +1.00, -0.05);
            var aVertex4 = new CGeoCoordinateDouble(+0.00, +1.00, -0.05);
            var aVertex5 = new CGeoCoordinateDouble(+0.00, +0.00, +0.05);
            var aVertex6 = new CGeoCoordinateDouble(+1.00, +0.00, +0.05);
            var aVertex7 = new CGeoCoordinateDouble(+1.00, +1.00, +0.05);
            var aVertex8 = new CGeoCoordinateDouble(+0.00, +1.00, +0.05);

            var aHexahedron = new CGeoHexahedronDouble();

            aHexahedron.addVertex(aVertex1);
            aHexahedron.addVertex(aVertex2);
            aHexahedron.addVertex(aVertex3);
            aHexahedron.addVertex(aVertex4);
            aHexahedron.addVertex(aVertex5);
            aHexahedron.addVertex(aVertex6);
            aHexahedron.addVertex(aVertex7);
            aHexahedron.addVertex(aVertex8);

            aMesher.generate(aHexahedron, 40, 40, 1, true);

            aMesher.mesh().generateFaces(1E-3);

        }

        static public void DrawMesh(CMshBasicMesherDouble aMesher, vtkUnstructuredGrid anUnstructuredGrid)
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

                if (anElement.elementType() == EElementType.ET_TETRAHEDRON)
                {

                    var aTetrahedron = vtkTetra.New();
                    aTetrahedron.GetPointIds().SetId(0, aMesher.mesh().nodeIndex(anElement.nodeId(0)));
                    aTetrahedron.GetPointIds().SetId(1, aMesher.mesh().nodeIndex(anElement.nodeId(1)));
                    aTetrahedron.GetPointIds().SetId(2, aMesher.mesh().nodeIndex(anElement.nodeId(2)));
                    aTetrahedron.GetPointIds().SetId(3, aMesher.mesh().nodeIndex(anElement.nodeId(3)));

                    anUnstructuredGrid.InsertNextCell(aTetrahedron.GetCellType(), aTetrahedron.GetPointIds());

                }

            }

        }

        static public void Simulate(CMshBasicMesherDouble aMesher, vtkUnstructuredGrid anUnstructuredGrid, vtkRenderer aRenderer)
        {

            var T = new CPdeFieldDouble();

            T.setMesh(aMesher.mesh());
            T.setDiscretMethod(EDiscretMethod.DM_FEM);
            T.setDiscretOrder(EDiscretOrder.DO_LINEAR);
            T.setDiscretLocation(EDiscretLocation.DL_NODE);
            T.setSimulationType(ESimulationType.ST_THERMAL);
            T.setNbDofs(1);

            // Set BC and initial conditions
            T.setSize(T.mesh().nbNodes());

            for (int i = 0; i < T.mesh().nbNodes(); ++i)
            {

                int aNodeId = T.mesh().nodeId(i);
                var aNode = T.mesh().node(aNodeId);

                if (Math.Abs(aNode.x() + 0.0) < 1E-6)
                    T.setFixedValue(i, 0.0);

                if (Math.Abs(aNode.x() - 1.0) < 1E-6)
                    T.setFixedValue(i, 0.0);

                if (Math.Abs(aNode.y() + 0.0) < 1E-6)
                    T.setFixedValue(i, 0.0);

                if (Math.Abs(aNode.y() - 1.0) < 1E-6)
                    T.setFixedValue(i, 1.0);

                T.setValue(i, 0.0);

            }

            // Steady conduction in a rectangle
            var aPdeEquation = new CPdeEquationDouble(ENigMA.laplacian(T));

            aPdeEquation.setElimination(T);

            aPdeEquation.solve(T);

            var sScalars = vtkFloatArray.New();

            sScalars.SetNumberOfTuples(T.mesh().nbNodes());

            for (int i = 0; i < T.mesh().nbNodes(); ++i)
            {

                int aNodeId = T.mesh().nodeId(i);

                sScalars.SetTuple1(aNodeId, T.value(i));

            }

            anUnstructuredGrid.GetPointData().SetScalars(sScalars);

            var aGeometryFilter = vtkGeometryFilter.New();
            aGeometryFilter.SetInput(anUnstructuredGrid);
            aGeometryFilter.Update();

            var aBoundedFilter = vtkBandedPolyDataContourFilter.New();
            aBoundedFilter.SetInput(aGeometryFilter.GetOutput());
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

            var aMesher = new CMshBasicMesherDouble();

            FemHeat3D.GenerateMesh(aMesher, aRenderer);

            var anUnstructuredGrid = vtkUnstructuredGrid.New();

            FemHeat3D.DrawMesh(aMesher, anUnstructuredGrid);

            FemHeat3D.Simulate(aMesher, anUnstructuredGrid, aRenderer);

            aRenderer.ResetCamera();

            aRenderControl.Refresh();

        }

    }
}
