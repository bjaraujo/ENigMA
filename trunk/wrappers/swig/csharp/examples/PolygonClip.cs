using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Kitware.VTK;

namespace Demo
{
    class PolygonClip
    {

        static public void DrawPolygon(vtkRenderer aRenderer)
        {

            var aVertex1 = new CGeoCoordinate(+0.00, +0.00, -0.5);
            var aVertex2 = new CGeoCoordinate(+1.00, +0.00, -0.5);
            var aVertex3 = new CGeoCoordinate(+1.00, +1.00, -0.5);
            var aVertex4 = new CGeoCoordinate(+0.00, +1.00, -0.5);
            var aVertex5 = new CGeoCoordinate(+0.00, +0.00, +0.5);
            var aVertex6 = new CGeoCoordinate(+1.00, +0.00, +0.5);
            var aVertex7 = new CGeoCoordinate(+1.00, +1.00, +0.5);
            var aVertex8 = new CGeoCoordinate(+0.00, +1.00, +0.5);

            var aHexahedron = new CGeoHexahedron();

            aHexahedron.addVertex(aVertex1);
            aHexahedron.addVertex(aVertex2);
            aHexahedron.addVertex(aVertex3);
            aHexahedron.addVertex(aVertex4);
            aHexahedron.addVertex(aVertex5);
            aHexahedron.addVertex(aVertex6);
            aHexahedron.addVertex(aVertex7);
            aHexahedron.addVertex(aVertex8);

            var aPolyhedron = new CGeoPolyhedron(aHexahedron);

            double aFracAct;

            int aNewPolygonId = 999;
            var aNewPolygon = new CGeoPolygon();

            var aNewNormal = new CGeoNormal(0.45, 0.45, 0.1);
            aNewNormal.normalize();

            var aNewPlane = new CGeoPlane();
            double d;

            int nIter;

            aPolyhedron = aPolyhedron.clip(aNewPolygon, aNewPolygonId, aNewNormal, out d, 0.6, out aFracAct, out nIter, 100, 1E-5);

            var sPoints = vtkPoints.New();
            var sPolygons = vtkCellArray.New();

            var sColors = vtkUnsignedCharArray.New();

            sColors.SetNumberOfComponents(3);
            sColors.SetName("Colors");

            int np = 0;

            for (int i = 0; i < aPolyhedron.nbPolygons(); ++i)
            {

                var aPolygon = vtkPolygon.New();

                int aPolygonId = aPolyhedron.polygonId(i);

                aPolygon.GetPointIds().SetNumberOfIds(aPolyhedron.polygon(aPolygonId).polyline().nbVertices() - 1);

                for (int j = 0; j < aPolyhedron.polygon(aPolygonId).polyline().nbVertices() - 1; ++j)
                {

                    var aVertex = aPolyhedron.polygon(aPolygonId).polyline().vertex(j);

                    aPolygon.GetPointIds().SetId(j, np);

                    sPoints.InsertNextPoint(aVertex.x(), aVertex.y(), aVertex.z());

                    if (aPolygonId == aNewPolygonId)
                    {
                        sColors.InsertNextTuple3(255, 255, 0);
                    }
                    else
                    {
                        sColors.InsertNextTuple3(0, 255, 255);
                    }

                    np++;

                }

                // Add the polygon to a list of polygons
                sPolygons.InsertNextCell(aPolygon);

            }

            var aPolygonPolyData = vtkPolyData.New();

            aPolygonPolyData.SetPoints(sPoints);
            aPolygonPolyData.SetPolys(sPolygons);

            aPolygonPolyData.GetPointData().SetScalars(sColors);

            // add actor to the renderer
            var anActor = vtkActor.New();
            aRenderer.AddActor(anActor);

            var aMapper = vtkPolyDataMapper.New();
            aMapper.SetInput(aPolygonPolyData);

            anActor.SetMapper(aMapper);

            anActor.SetMapper(aMapper);

        }

        static public void Show(RenderWindowControl aRenderControl)
        {

            // get a reference to the renderwindow of our renderWindowControl1
            var aRenderWindow = aRenderControl.RenderWindow;

            // get a reference to the renderer
            var aRenderer = aRenderWindow.GetRenderers().GetFirstRenderer();

            aRenderer.Clear();
            aRenderer.RemoveAllViewProps();

            PolygonClip.DrawPolygon(aRenderer);

            aRenderer.ResetCamera();

            aRenderControl.Refresh();
        }
    }
}
