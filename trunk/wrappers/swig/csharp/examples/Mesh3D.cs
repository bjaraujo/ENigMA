using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Kitware.VTK;

namespace Demo
{
    class Mesh3D
    {

        static public void GenerateMesh(CMshTetrahedronMesherDouble aMesher, vtkRenderer aRenderer)
        {

            var aVertex1 = new CGeoCoordinateDouble(0.0, 0.0, 0.0);
            var aVertex2 = new CGeoCoordinateDouble(1.0, 0.0, 0.0);
            var aVertex3 = new CGeoCoordinateDouble(1.0, 1.0, 0.0);
            var aVertex4 = new CGeoCoordinateDouble(0.0, 1.0, 0.0);
            var aVertex5 = new CGeoCoordinateDouble(0.0, 0.0, 1.0);
            var aVertex6 = new CGeoCoordinateDouble(1.0, 0.0, 1.0);
            var aVertex7 = new CGeoCoordinateDouble(1.0, 1.0, 1.0);
            var aVertex8 = new CGeoCoordinateDouble(0.0, 1.0, 1.0);

            var aHexahedron = new CGeoHexahedronDouble();

            aHexahedron.addVertex(aVertex1);
            aHexahedron.addVertex(aVertex2);
            aHexahedron.addVertex(aVertex3);
            aHexahedron.addVertex(aVertex4);
            aHexahedron.addVertex(aVertex5);
            aHexahedron.addVertex(aVertex6);
            aHexahedron.addVertex(aVertex7);
            aHexahedron.addVertex(aVertex8);

            var aBasicMesher = new CMshBasicMesherDouble();

            aBasicMesher.generate(aHexahedron, 10, 10, 10, true);

            CMshMeshDouble aBoundaryMesh = aBasicMesher.mesh().extractBoundary(1E-6);

            aBoundaryMesh.generateFaces(1E-4);

            var aMeshSize = 0.1;

            aMesher.generate(aBoundaryMesh, 99999, aMeshSize, 0.1, 1E-4);

            aMesher.mesh().renumber();

        }

        static public void DrawMesh(CMshTetrahedronMesherDouble aMesher, vtkUnstructuredGrid anUnstructuredGrid, vtkRenderer aRenderer)
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

            var aMesher = new CMshTetrahedronMesherDouble();

            Mesh3D.GenerateMesh(aMesher, aRenderer);

            var anUnstructuredGrid = vtkUnstructuredGrid.New();

            Mesh3D.DrawMesh(aMesher, anUnstructuredGrid, aRenderer);

            aRenderer.ResetCamera();

            aRenderControl.Refresh();

        }

    }
}
