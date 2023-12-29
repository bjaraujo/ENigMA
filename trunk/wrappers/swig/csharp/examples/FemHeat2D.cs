﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Kitware.VTK;

namespace Demo
{
    class FemHeat2D
    {

        static public void GenerateMesh(CMshTriangleMesher aMesher, vtkRenderer aRenderer)
        {

            var aNode1 = new CMshNode(0.0, 0.0, 0.0);
            var aNode2 = new CMshNode(1.0, 0.0, 0.0);
            var aNode3 = new CMshNode(1.0, 1.0, 0.0);
            var aNode4 = new CMshNode(0.0, 1.0, 0.0);

            var anEdgeMesh = new CMshMesh();

            anEdgeMesh.addNode(1, aNode1);
            anEdgeMesh.addNode(2, aNode2);
            anEdgeMesh.addNode(3, aNode3);
            anEdgeMesh.addNode(4, aNode4);

            var anElement1 = new CMshElement(EElementType.ET_BEAM);
            anElement1.addNodeId(1);
            anElement1.addNodeId(2);

            var anElement2 = new CMshElement(EElementType.ET_BEAM);
            anElement2.addNodeId(2);
            anElement2.addNodeId(3);

            var anElement3 = new CMshElement(EElementType.ET_BEAM);
            anElement3.addNodeId(3);
            anElement3.addNodeId(4);

            var anElement4 = new CMshElement(EElementType.ET_BEAM);
            anElement4.addNodeId(4);
            anElement4.addNodeId(1);

            anEdgeMesh.addElement(1, anElement1);
            anEdgeMesh.addElement(2, anElement2);
            anEdgeMesh.addElement(3, anElement3);
            anEdgeMesh.addElement(4, anElement4);

            var aMeshSize = 0.02;

            anEdgeMesh.generateFaces(1E-3);

            aMesher.remesh(anEdgeMesh, aMeshSize);
            aMesher.generate(anEdgeMesh, 99999, aMeshSize, 0.1);

            for (int i = 0; i < 3; i++)
            {
                aMesher.relaxNodes(1E-5);
                aMesher.flipEdges(1E-5);
            }

        }

        static public void DrawMesh(CMshTriangleMesher aMesher, vtkUnstructuredGrid anUnstructuredGrid)
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

                if (anElement.elementType() == EElementType.ET_TRIANGLE)
                {

                    var aTriangle = vtkTriangle.New();
                    aTriangle.GetPointIds().SetId(0, aMesher.mesh().nodeIndex(anElement.nodeId(0)));
                    aTriangle.GetPointIds().SetId(1, aMesher.mesh().nodeIndex(anElement.nodeId(1)));
                    aTriangle.GetPointIds().SetId(2, aMesher.mesh().nodeIndex(anElement.nodeId(2)));

                    anUnstructuredGrid.InsertNextCell(aTriangle.GetCellType(), aTriangle.GetPointIds());

                }

            }

        }

        static public void Simulate(CMshTriangleMesher aMesher, vtkUnstructuredGrid anUnstructuredGrid, vtkRenderer aRenderer)
        {

            var T = new CPdeField();

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
            var aPdeEquation = new CPdeEquation(ENigMA.laplacian(T));

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

            var aMesher = new CMshTriangleMesher();

            FemHeat2D.GenerateMesh(aMesher, aRenderer);

            var anUnstructuredGrid = vtkUnstructuredGrid.New();

            FemHeat2D.DrawMesh(aMesher, anUnstructuredGrid);

            FemHeat2D.Simulate(aMesher, anUnstructuredGrid, aRenderer);

            aRenderer.ResetCamera();

            aRenderControl.Refresh();

        }

    }
}
