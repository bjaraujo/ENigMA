using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace CSDemo
{
    public partial class Form1 : Form
    {

        public Form1()
        {
            InitializeComponent();
        }

        private void Form1_Load(object sender, EventArgs e)
        {

            comboBox1.Items.Add("Line Intersection");
            comboBox1.Items.Add("Line Distance");
            comboBox1.Items.Add("Polygon Clipping");
            comboBox1.Items.Add("Mesh Generation");

            comboBox1.SelectedIndex = 0;

        }

        private void DrawPoint(PaintEventArgs e, GeoCoordinate aPoint, Color aColor, int aRadius)
        {

            Pen aPen = new Pen(aColor);
            e.Graphics.DrawEllipse(aPen, (float)aPoint.x() - aRadius * 0.5f, (float)aPoint.y() - aRadius * 0.5f, aRadius, aRadius);

        }

        private void DrawLine(PaintEventArgs e, GeoLine aLine, Color aColor)
        {

            Pen aPen = new Pen(aColor);
            e.Graphics.DrawLine(aPen, (float)aLine.startPoint().x(), (float)aLine.startPoint().y(), (float)aLine.endPoint().x(), (float)aLine.endPoint().y());

        }

        private void DrawTriangle(PaintEventArgs e, GeoTriangle aTriangle, Color aColor)
        {

            Pen aPen = new Pen(aColor);
            PointF[] sPoints = new PointF[4];

            sPoints[0] = new PointF((float)aTriangle.vertex(0).x(), (float)aTriangle.vertex(0).y());
            sPoints[1] = new PointF((float)aTriangle.vertex(1).x(), (float)aTriangle.vertex(1).y());
            sPoints[2] = new PointF((float)aTriangle.vertex(2).x(), (float)aTriangle.vertex(2).y());
            sPoints[3] = new PointF((float)aTriangle.vertex(0).x(), (float)aTriangle.vertex(0).y());

            e.Graphics.DrawLines(aPen, sPoints);

        }

        private void DrawQuadrilateral(PaintEventArgs e, GeoQuadrilateral aQuadrilateral, Color aColor)
        {

            Pen aPen = new Pen(aColor);
            PointF[] sPoints = new PointF[5];

            sPoints[0] = new PointF((float)aQuadrilateral.vertex(0).x(), (float)aQuadrilateral.vertex(0).y());
            sPoints[1] = new PointF((float)aQuadrilateral.vertex(1).x(), (float)aQuadrilateral.vertex(1).y());
            sPoints[2] = new PointF((float)aQuadrilateral.vertex(2).x(), (float)aQuadrilateral.vertex(2).y());
            sPoints[3] = new PointF((float)aQuadrilateral.vertex(3).x(), (float)aQuadrilateral.vertex(3).y());
            sPoints[4] = new PointF((float)aQuadrilateral.vertex(0).x(), (float)aQuadrilateral.vertex(0).y());

            e.Graphics.DrawLines(aPen, sPoints);

        }

        private void DrawPolygon(PaintEventArgs e, GeoPolygon aPolygon, Color aColor)
        {

            var aPolyline = aPolygon.polyline();

            for (int i = 0; i < aPolyline.nbLines(); i++)
            {

                var aLine = aPolyline.line(i);

                DrawLine(e, aLine, aColor);

            }

        }

        private void DrawLineIntersection(PaintEventArgs e)
        {

            var aPoint1 = new GeoCoordinate(10.0, 2.0 * trackBar1.Value + 10.0, 0.0);
            var aPoint2 = new GeoCoordinate(400.0, 400.0, 0.0);

            var aPoint3 = new GeoCoordinate(10.0, 305.0, 0.0);
            var aPoint4 = new GeoCoordinate(400.0, 10.0, 0.0);

            var aLine1 = new GeoLine(aPoint1, aPoint2);
            var aLine2 = new GeoLine(aPoint3, aPoint4);

            DrawLine(e, aLine1, Color.FromArgb(255, 255, 0, 0));
            DrawLine(e, aLine2, Color.FromArgb(255, 0, 255, 0));

            var aPoint = new GeoCoordinate();

            if (aLine1.intersects(aLine2, aPoint, 1E-6))
            {
                DrawPoint(e, aPoint, Color.FromArgb(255, 0, 0, 255), 10);
            }
        }

        private void DrawLineDistance(PaintEventArgs e)
        {

            var aPoint1 = new GeoCoordinate(10.0, 2.0 * trackBar1.Value + 10.0, 0.0);
            var aPoint2 = new GeoCoordinate(400.0, 400.0, 0.0);

            var aPoint3 = new GeoCoordinate(10.0, 10.0, 0.0);
            var aPoint4 = new GeoCoordinate(400.0, 3.5 * trackBar1.Value + 10.0, 0.0);

            var aLine1 = new GeoLine(aPoint1, aPoint2);
            var aLine2 = new GeoLine(aPoint3, aPoint4);

            DrawLine(e, aLine1, Color.FromArgb(255, 255, 0, 0));
            DrawLine(e, aLine2, Color.FromArgb(255, 0, 255, 0));

            var aPoint5 = new GeoCoordinate();
            var aPoint6 = new GeoCoordinate();

            double aDistance;

            if (aLine1.distance(aLine2, aPoint5, aPoint6, out aDistance))
            {

                DrawPoint(e, aPoint5, Color.FromArgb(255, 0, 0, 255), 10);
                DrawPoint(e, aPoint6, Color.FromArgb(255, 0, 0, 255), 10);

                label1.Text = aDistance.ToString();
            }

        }

        private void DrawMesh(PaintEventArgs e)
        {

            var aNode1 = new MshNode(10, 10, 0);
            var aNode2 = new MshNode(400, 10, 0);
            var aNode3 = new MshNode(400, 400, 0);
            var aNode4 = new MshNode(10, 400, 0);

            var anEdgeMesh = new MshMesh();

            anEdgeMesh.addNode(1, aNode1);
            anEdgeMesh.addNode(2, aNode2);
            anEdgeMesh.addNode(3, aNode3);
            anEdgeMesh.addNode(4, aNode4);

            var anElement1 = new MshElement(eElementType.ET_BEAM);
            anElement1.addNodeId(1);
            anElement1.addNodeId(2);

            var anElement2 = new MshElement(eElementType.ET_BEAM);
            anElement2.addNodeId(2);
            anElement2.addNodeId(3);

            var anElement3 = new MshElement(eElementType.ET_BEAM);
            anElement3.addNodeId(3);
            anElement3.addNodeId(4);

            var anElement4 = new MshElement(eElementType.ET_BEAM);
            anElement4.addNodeId(4);
            anElement4.addNodeId(1);

            anEdgeMesh.addElement(1, anElement1);
            anEdgeMesh.addElement(2, anElement2);
            anEdgeMesh.addElement(3, anElement3);
            anEdgeMesh.addElement(4, anElement4);

            var aMesher = new MshTriangleMesher();
            //var aMesher = new MshQuadrilateralMesher();

            var aMeshSize = (100 - trackBar1.Value) + 15.0;

            anEdgeMesh.generateFaces(1E-3);

            aMesher.remesh(anEdgeMesh, aMeshSize);
            aMesher.generate(anEdgeMesh, 99999, aMeshSize, 0.1);

            for (int i = 0; i < 3; i++)
            {
                aMesher.relaxNodes(1E-5);
                aMesher.flipEdges(1E-5);
            }

            for (int i = 0; i < aMesher.mesh().nbElements(); i++)
            {

                var anElementId = aMesher.mesh().elementId(i);
                var anElement = aMesher.mesh().element(anElementId);

                if (anElement.elementType() != eElementType.ET_TRIANGLE)
                {
                    continue;
                }

                var aTriangle = new GeoTriangle();

                aTriangle.addVertex(aMesher.mesh().node(anElement.nodeId(0)));
                aTriangle.addVertex(aMesher.mesh().node(anElement.nodeId(1)));
                aTriangle.addVertex(aMesher.mesh().node(anElement.nodeId(2)));

                DrawTriangle(e, aTriangle, Color.FromArgb(255, 0, 0, 255));

                /*
                if (anElement.elementType() != eElementType.ET_QUADRILATERAL)
                    continue;

                var aQuadrilateral = new GeoQuadrilateral();

                aQuadrilateral.addVertex(aMesher.mesh().node(anElement.nodeId(0)));
                aQuadrilateral.addVertex(aMesher.mesh().node(anElement.nodeId(1)));
                aQuadrilateral.addVertex(aMesher.mesh().node(anElement.nodeId(2)));
                aQuadrilateral.addVertex(aMesher.mesh().node(anElement.nodeId(3)));

                DrawQuadrilateral(e, aQuadrilateral, Color.FromArgb(255, 0, 0, 255));
                */
            }

        }

        private void DrawPolygonClip(PaintEventArgs e)
        {

            var aPoint1 = new GeoCoordinate(20.0, 20.0, 0.0);
            var aPoint2 = new GeoCoordinate(30.0, 200.0, 0.0);
            var aPoint3 = new GeoCoordinate(100.0, 300.0, 0.0);
            var aPoint4 = new GeoCoordinate(400.0, 350.0, 0.0);
            var aPoint5 = new GeoCoordinate(350.0, 150.0, 0.0);

            var aPolyline = new GeoPolyline();

            aPolyline.addVertex(aPoint1);
            aPolyline.addVertex(aPoint2);
            aPolyline.addVertex(aPoint3);
            aPolyline.addVertex(aPoint4);
            aPolyline.addVertex(aPoint5);
            aPolyline.close();

            var aPolygon = new GeoPolygon();

            aPolygon.setPolyline(aPolyline);

            DrawPolygon(e, aPolygon, Color.FromArgb(255, 0, 0, 255));

            var aNormal = new GeoNormal(1, 1, 0);

            var aPlane = new GeoPlane(aNormal, 3.0 * trackBar1.Value + 100.0);

            var aClippedPolygon = aPolygon.clip(aPlane);

            DrawPolygon(e, aClippedPolygon, Color.FromArgb(255, 255, 0, 0));


        }

        private void comboBox1_SelectedIndexChanged(object sender, EventArgs e)
        {

            Refresh();

        }

        private void trackBar1_Scroll(object sender, EventArgs e)
        {

            Refresh();

        }

        private void pictureBox1_Paint(object sender, PaintEventArgs e)
        {

            switch (comboBox1.SelectedIndex)
            {
                case 0:
                    DrawLineIntersection(e);
                    break;
                case 1:
                    DrawLineDistance(e);
                    break;
                case 2:
                    DrawPolygonClip(e);
                    break;
                case 3:
                    DrawMesh(e);
                    break;
                default:
                    break;
            }

        }


    }
}
