using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Kitware.VTK;

namespace Demo
{
    class FvmLidDriven2D
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

        static public void Simulate(CMshBasicMesherDouble aMesher, vtkUnstructuredGrid anUnstructuredGrid, vtkRenderer aRenderer, RenderWindowControl aRenderControl)
        {

            aMesher.mesh().generateFaces(1E-12);

            aMesher.mesh().calculateFaceCentroid();
            aMesher.mesh().calculateElementCentroid();

            var aFvmMesh = new CFvmMesh(aMesher.mesh());

            double U = 1.0; // Lid velocity
            double mu = 0.001; // dynamic viscosity
            double rho = 1.0; // density

            double nu = mu / rho; // kinematic viscosity 

            var aPisoSolver = new CFvmPisoSolver(aFvmMesh);

            aPisoSolver.setGravity(0.0, 0.0, 0.0);
            aPisoSolver.setMaterialProperties(rho, mu);

            var sFaceIds = new StdVectorInt();

            sFaceIds.Clear();
            for (int i = 0; i < aFvmMesh.nbFaces(); ++i)
            {   

                int aFaceId = aFvmMesh.faceId(i);

                var aFace = aFvmMesh.face(aFaceId);

                aFace.calculateCentroid();

                if (aFace.centroid().x() == 0.0 ||
                    aFace.centroid().x() == 1.0 ||
                    aFace.centroid().y() == 0.0 ||
                    aFace.centroid().y() == 1.0)
                {
                    sFaceIds.Add(aFaceId);
                }

            }

            aPisoSolver.setBoundaryVelocity(sFaceIds, EBoundaryType.BT_WALL_NO_SLIP, 0.0, 0.0, 0.0);
            aPisoSolver.setBoundaryPressure(sFaceIds, EBoundaryType.BT_WALL_NO_SLIP, 0.0);

            sFaceIds.Clear();
            for (int i = 0; i < aFvmMesh.nbFaces(); ++i)
            {

                int aFaceId = aFvmMesh.faceId(i);

                var aFace = aFvmMesh.face(aFaceId);

                aFace.calculateCentroid();

                if (aFace.centroid().y() == 1.0)
                {
                    sFaceIds.Add(aFaceId);
                }

            }

            aPisoSolver.setBoundaryVelocity(sFaceIds, EBoundaryType.BT_WALL_NO_SLIP, U, 0.0, 0.0);

            double dt = U / 40;			// Courant < 1

            var sScalars = vtkFloatArray.New();
            sScalars.SetNumberOfTuples(aFvmMesh.nbControlVolumes());
            var aGeometryFilter = vtkGeometryFilter.New();
            var aCellDataToPointDataFilter1 = vtkCellDataToPointData.New();
            var aBoundedFilter = vtkBandedPolyDataContourFilter.New();
            var aLookupTable = vtkLookupTable.New();
            aLookupTable.SetNumberOfColors(256);
            aLookupTable.SetHueRange(0.667, 0.0);
            aLookupTable.Build();
            var aBandedMapper = vtkPolyDataMapper.New();

            // add actor to the renderer
            var aMeshActor = vtkActor.New();
            aRenderer.AddActor(aMeshActor);

            var aScalarBarActor = vtkScalarBarActor.New();
            aRenderer.AddActor(aScalarBarActor);

            int nIter = 20;

            // Flow in a rectangle
            for (int ii = 0; ii < nIter; ++ii)
            {

                aPisoSolver.iterate(dt);

                for (int i = 0; i < aFvmMesh.nbControlVolumes(); ++i)
                {

                    var aControlVolumeId = aFvmMesh.controlVolumeId(i);

                    double u = aPisoSolver.u(aControlVolumeId);
                    double v = aPisoSolver.v(aControlVolumeId);

                    double vel = Math.Sqrt(u * u + v * v);

                    sScalars.SetTuple1(aControlVolumeId, vel);

                }

                sScalars.Modified();

                anUnstructuredGrid.GetCellData().SetScalars(sScalars);
                aGeometryFilter.SetInput(anUnstructuredGrid);
                aGeometryFilter.Update();

                aCellDataToPointDataFilter1.SetInputConnection(aGeometryFilter.GetOutputPort());
                aCellDataToPointDataFilter1.Update();

                aBoundedFilter.SetInput(aCellDataToPointDataFilter1.GetOutput());
                aBoundedFilter.GenerateValues(24, sScalars.GetRange()[0], sScalars.GetRange()[1]);

                aBandedMapper.SetInputConnection(aBoundedFilter.GetOutputPort());
                aBandedMapper.SetScalarModeToUsePointData();
                aBandedMapper.SetScalarRange(sScalars.GetRange()[0], sScalars.GetRange()[1]);
                aBandedMapper.SetLookupTable(aLookupTable);

                aMeshActor.SetMapper(aBandedMapper);

                aScalarBarActor.SetLookupTable(aLookupTable);
                aScalarBarActor.SetTitle("Velocity");
                aScalarBarActor.SetNumberOfLabels(6);

                aRenderer.Render();

                if (ii == 0)
                {
                    aRenderer.ResetCamera();
                }

                aRenderControl.Refresh();
            }

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

            FvmLidDriven2D.GenerateMesh(aMesher, aRenderer);

            var anUnstructuredGrid = vtkUnstructuredGrid.New();

            FvmLidDriven2D.DrawMesh(aMesher, anUnstructuredGrid);

            FvmLidDriven2D.Simulate(aMesher, anUnstructuredGrid, aRenderer, aRenderControl);

            aRenderer.ResetCamera();

            aRenderControl.Refresh();
        }
    }
}
