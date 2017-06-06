using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Diagnostics;

using Kitware.VTK;

namespace CSDemo
{
    public partial class Form2 : Form
    {

        public Form2()
        {
            InitializeComponent();
        }

        private void GenerateMesh(CMshTriangleMesherDouble aMesher)
        {

            var aNode1 = new CMshNodeDouble(0, 0, 0);
            var aNode2 = new CMshNodeDouble(1, 0, 0);
            var aNode3 = new CMshNodeDouble(1, 1, 0);
            var aNode4 = new CMshNodeDouble(0, 1, 0);

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

        private void DrawMesh(CMshTriangleMesherDouble aMesher, vtkUnstructuredGrid anUnstructuredGrid)
        {

            vtkPoints points = vtkPoints.New();

            for (int i = 0; i < aMesher.mesh().nbNodes(); ++i)
            {

                int aNodeId = aMesher.mesh().nodeId(i);

                var aNode = aMesher.mesh().node(aNodeId);

                points.InsertNextPoint(aNode.x(), aNode.y(), aNode.z());

            }

            anUnstructuredGrid.SetPoints(points);

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

        private void Simulate(CMshTriangleMesherDouble aMesher, vtkUnstructuredGrid anUnstructuredGrid, vtkRenderer aRenderer)
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
                    T.setFixedValue(i, 0.0);

                T.setValue(i, 1.0);

            }

            double rho = 1.0;   // density
            double Cp = 1.0;    // specific heat
            double k = 1.0;     // conductivity

            double time = 0.1;
            int nIter = 10;
            double dt = time / nIter;

            // Unsteady conduction in a plate
            for (int i = 0; i < nIter; ++i)
            {

                Debug.WriteLine("Time = " + dt * (i + 1));

                var aTransientTerm = new CSleSystemDouble();
                aTransientTerm = ENigMA.ddt(T);
                aTransientTerm.Multiply(rho * Cp / dt);

                var aDiffusionTerm = new CSleSystemDouble();
                aDiffusionTerm = ENigMA.laplacian(T);
                aDiffusionTerm.Multiply(k);

                var aSystem = new CSleSystemDouble();
                aSystem = aTransientTerm;
                aSystem.Plus(aDiffusionTerm);

                var aPdeEquation = new CPdeEquationDouble(aSystem);

                aPdeEquation.setElimination(T);

                aPdeEquation.solve(T);

            }

            var sScalars = vtkFloatArray.New();

            sScalars.SetNumberOfTuples(T.mesh().nbNodes());
    
            for (int i = 0; i < T.mesh().nbNodes(); ++i)
            {
            
                int aControlVolumeId = T.mesh().nodeId(i);

                sScalars.SetTuple1(aControlVolumeId, T.value(i));

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
            var MeshActor = vtkActor.New();
            aRenderer.AddActor(MeshActor);

            var MeshMapper = vtkDataSetMapper.New();
            MeshMapper.SetInput(anUnstructuredGrid);

            MeshActor.SetMapper(aBandedMapper);

        }

        private void Form1_Activated(object sender, EventArgs e)
        {

            var aMesher = new CMshTriangleMesherDouble();

            GenerateMesh(aMesher);

            // get a reference to the renderwindow of our renderWindowControl1
            var aRenderWindow = renderWindowControl1.RenderWindow;

            // get a reference to the renderer
            var aRenderer = aRenderWindow.GetRenderers().GetFirstRenderer();

            // set background color
            aRenderer.SetBackground(0.2, 0.3, 0.4);

            var anUnstructuredGrid = vtkUnstructuredGrid.New();

            DrawMesh(aMesher, anUnstructuredGrid);

            Simulate(aMesher, anUnstructuredGrid, aRenderer);

            aRenderer.ResetCamera();

        }

    }
}
