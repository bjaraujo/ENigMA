using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Kitware.VTK;

namespace Demo
{
    class CSG
    {

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

            var aCenter = new CGeoCoordinateDouble(0, 0, 0);

            var aRadius = new CGeoVectorDouble(1.0);

            var a = CCsgCubeDouble.create(aCenter, aRadius);
            var b = CCsgSphereDouble.create(aCenter, 1.35, 16, 12);

            var aVertex_x1 = new CGeoCoordinateDouble(-2, 0, 0);
            var aVertex_x2 = new CGeoCoordinateDouble(+2, 0, 0);

            var c = CCsgCylinderDouble.create(aVertex_x1, aVertex_x2, 0.7, 16);

            var aVertex_y1 = new CGeoCoordinateDouble(0, -2, 0);
            var aVertex_y2 = new CGeoCoordinateDouble(0, +2, 0);

            var d = CCsgCylinderDouble.create(aVertex_y1, aVertex_y2, 0.7, 16);

            var aVertex_z1 = new CGeoCoordinateDouble(0, 0, -2);
            var aVertex_z2 = new CGeoCoordinateDouble(0, 0, +2);

            var e = CCsgCylinderDouble.create(aVertex_z1, aVertex_z2, 0.7, 16);

            a.setTolerance(1E-6);
            b.setTolerance(1E-6);
            c.setTolerance(1E-6);
            d.setTolerance(1E-6);
            e.setTolerance(1E-6);

            var f = a.intersect(b).subtract(c.add(d).add(e));

            var aFileName = "demo.stl";

            f.toStl().save(aFileName);

            var aReader = vtkSTLReader.New();

            aReader.SetFileName(aFileName);
            aReader.Update();

            var aMapper = vtkPolyDataMapper.New();
            aMapper.SetInputConnection(aReader.GetOutputPort());

            var anActor = vtkActor.New();
            anActor.SetMapper(aMapper);
            aRenderer.AddActor(anActor);

            aRenderer.ResetCamera();

            aRenderControl.Refresh();

        }

    }
}
