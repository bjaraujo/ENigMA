using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Kitware.VTK;

namespace Demo
{
    class Mesh2D
    {

        static public void AddCircle(CMshMeshDouble anEdgeMesh, double cx, double cy, double radius)
        {

            int aFirstNodeId = anEdgeMesh.nextNodeId();

            for (int i = 0; i < 16; i++)
            {

                var teta = i / 16.0 * 2 * Math.PI;

                var aNode5 = new CMshNodeDouble(cx - radius * Math.Cos(teta), cy + radius * Math.Sin(teta), 0.0);
                anEdgeMesh.addNode(anEdgeMesh.nextNodeId(), aNode5);

                if (i > 0)
                {
                    var anElement = new CMshElementDouble(EElementType.ET_BEAM);
                    anElement.addNodeId(anEdgeMesh.nextNodeId() - 2);
                    anElement.addNodeId(anEdgeMesh.nextNodeId() - 1);
                    anEdgeMesh.addElement(anEdgeMesh.nextElementId(), anElement);
                }

            }

            {
                var anElement = new CMshElementDouble(EElementType.ET_BEAM);
                anElement.addNodeId(anEdgeMesh.nextNodeId() - 1);
                anElement.addNodeId(aFirstNodeId);
                anEdgeMesh.addElement(anEdgeMesh.nextElementId(), anElement);
            }

        }

        static public void GenerateMesh(CMshTriangleMesherDouble aMesher, vtkRenderer aRenderer)
        {

            var aNode1 = new CMshNodeDouble(0.0, 0.0, 0.0);
            var aNode2 = new CMshNodeDouble(1.0, 0.0, 0.0);
            var aNode3 = new CMshNodeDouble(1.0, 1.0, 0.0);
            var aNode4 = new CMshNodeDouble(0.0, 1.0, 0.0);

            var anEdgeMesh = new CMshMeshDouble();

            anEdgeMesh.addNode(1, aNode1);
            anEdgeMesh.addNode(2, aNode2);
            anEdgeMesh.addNode(3, aNode3);
            anEdgeMesh.addNode(4, aNode4);

            var anElement1 = new CMshElementDouble(EElementType.ET_BEAM);
            anElement1.addNodeId(1);
            anElement1.addNodeId(2);

            var anElement2 = new CMshElementDouble(EElementType.ET_BEAM);
            anElement2.addNodeId(2);
            anElement2.addNodeId(3);

            var anElement3 = new CMshElementDouble(EElementType.ET_BEAM);
            anElement3.addNodeId(3);
            anElement3.addNodeId(4);

            var anElement4 = new CMshElementDouble(EElementType.ET_BEAM);
            anElement4.addNodeId(4);
            anElement4.addNodeId(1);

            anEdgeMesh.addElement(1, anElement1);
            anEdgeMesh.addElement(2, anElement2);
            anEdgeMesh.addElement(3, anElement3);
            anEdgeMesh.addElement(4, anElement4);

            Mesh2D.AddCircle(anEdgeMesh, 0.3, 0.3, 0.1);
            Mesh2D.AddCircle(anEdgeMesh, 0.7, 0.3, 0.1);
            Mesh2D.AddCircle(anEdgeMesh, 0.7, 0.7, 0.1);
            Mesh2D.AddCircle(anEdgeMesh, 0.3, 0.7, 0.1);

            var sLinesPolyData = vtkPolyData.New();

            vtkPoints sPoints = vtkPoints.New();

            for (int i = 0; i < anEdgeMesh.nbNodes(); ++i)
            {

                int aNodeId = anEdgeMesh.nodeId(i);

                var aNode = anEdgeMesh.node(aNodeId);

                sPoints.InsertNextPoint(aNode.x(), aNode.y(), aNode.z());

            }

            sLinesPolyData.SetPoints(sPoints);

            var sLines = vtkCellArray.New();

            for (int i = 0; i < anEdgeMesh.nbElements(); ++i)
            {
                int anElementId = anEdgeMesh.elementId(i);

                var anElement = anEdgeMesh.element(anElementId);

                var aLine = vtkLine.New();
                aLine.GetPointIds().SetId(0, anEdgeMesh.nodeIndex(anElement.nodeId(0)));
                aLine.GetPointIds().SetId(1, anEdgeMesh.nodeIndex(anElement.nodeId(1)));

                sLines.InsertNextCell(aLine);
            }

            sLinesPolyData.SetLines(sLines);

            // add actor to the renderer
            var anEdgeActor = vtkActor.New();
            aRenderer.AddActor(anEdgeActor);

            var anEdgeMapper = vtkDataSetMapper.New();
            anEdgeMapper.SetInput(sLinesPolyData);

            anEdgeActor.SetMapper(anEdgeMapper);

            var aMeshSize = 0.026;

            anEdgeMesh.generateFaces(1E-3);

            aMesher.remesh(anEdgeMesh, aMeshSize);
            aMesher.generate(anEdgeMesh, 99999, aMeshSize, 0.1, 1E-6);

            for (int i = 0; i < 3; i++)
            {
                aMesher.relaxNodes(1E-5);
                aMesher.flipEdges(1E-5);
            }

            aMesher.mesh().renumber();

        }

        static public void DrawMesh(CMshTriangleMesherDouble aMesher, vtkUnstructuredGrid anUnstructuredGrid, vtkRenderer aRenderer)
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

            // add actor to the renderer
            var aMeshActor = vtkActor.New();
            aRenderer.AddActor(aMeshActor);

            var aMeshMapper = vtkDataSetMapper.New();
            aMeshMapper.SetInput(anUnstructuredGrid);

            aMeshActor.SetMapper(aMeshMapper);

            aMeshActor.GetProperty().SetRepresentationToWireframe();

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

            var aMesher = new CMshTriangleMesherDouble();

            Mesh2D.GenerateMesh(aMesher, aRenderer);

            var anUnstructuredGrid = vtkUnstructuredGrid.New();

            Mesh2D.DrawMesh(aMesher, anUnstructuredGrid, aRenderer);

            aRenderer.ResetCamera();

            aRenderControl.Refresh();

        }

    }
}
