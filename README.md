|Build Status|License|
|:--:|:--:|
|[![Build Status](https://github.com/bjaraujo/ENigMA/workflows/CMake/badge.svg)](https://github.com/bjaraujo/ENigMA/actions)|[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/gpl-2.0)|

# ENigMA - Extended Numerical Multiphysics Analysis #

ENigMA is an object-oriented C++ template library which aim is to provide multi-physics simulation in a multi-domain environment. The code implements several numerical methods such as Finite Volume Methods (FVM), Finite Difference Methods (FDM), Finite Element Methods (FEM), Boundary Element Methods (BEM), Smoothed Particle Hydrodynamics (SPH), etc. for numerical approximation of Partial Differential Equations (PDE) in each domain. It also provides classes for robust mesh generation (triangular, block, constrained tetrahedral, etc), intersection and clipping operations and implements R-tree, octree and hashgrid methods for spatial searching. It can be used for three-dimensional flow, thermal and structural analysis.

It was developed to be cross-platform using STL, Eigen (for vectors and matrices) and exprtk (for math expressions). The SWIG tool is used to expose ENigMA's classes and methods to other languages such as Python, C#, etc. It also uses Gtest for unit testing (> 160 unit tests), CMake is used for cross-platform building and Git is used for source code management.

![mesh](https://github.com/bjaraujo/ENigMA/blob/master/images/mesh.png)

ENigMA implements several operators such as divergence, laplacian, gradient, etc. so, for example, a heat equation such as: 

![equation](https://github.com/bjaraujo/ENigMA/blob/master/images/equation.gif)

can be represented in ENigMA as:

```cpp
CPdeEquation<double> aPdeEquation{1.0 / dt * ddt<double>(u) - alpha * laplacian<double>(u) = 0};
```

The argument passed to the CPdeEquation accepts an object of type CSleSystem which is a linear system of equations. The solution of the PDE is achieved by simply calling solve method:

```cpp
aPdeEquation.solve(u);
```

### Repository ###

https://github.com/bjaraujo/ENigMA

### Structure ###

ENigMA has the following namespaces:

- **analytical**: analytical solutions
- **geometry**: geometry, intersection, clipping, spatial search (R-tree, octree, hashgrid), etc.
- **integration**: numerical integration
- **bem**: Boundary Element Method (BEM)
- **fdm**: Finite Difference Method (FDM)
- **fem**: Finite Element Method (FEM)
- **fvm**: Finite Volume Method (FVM)
- **sph**: Smoothed Particle Hydrodynamics (SPH)
- **pde**: Partial Differential Equations
- **sle**: System of Linear Equations
- **mesh**: mesh generation (block, triangular, tetrahedral*, hexahedral, etc.)
- **post**: post-processing
- **stl**: STL file processing

### Capabilities ###

- Mesh generation (2D/3D), 
- Structural mechanics (linear elasticity),
- Electro-magnetics (Laplace equation),
- Fluid flow and heat transfer (Navier-Stokes & Laplace equation),
- Finance (Black-Scholes equation).

### Dependencies ###

#### Core ####

- Eigen: http://eigen.tuxfamily.org
- Exprtk: https://github.com/ArashPartow/exprtk or dependencies folder
- RTree: https://github.com/nushoin/RTree (slightly modified) or dependencies folder
- ViennaCL (optional): http://viennacl.sourceforge.net

#### Tests ####

- Gtest (optional): https://github.com/google/googletest

#### Examples ####

- Qt (optional): https://www.qt.io
- VTK (optional): http://www.vtk.org
- OpenGL/GLUT (optional): http://freeglut.sourceforge.net
- OpenCASCADE (optional): https://www.opencascade.com
- Gmsh (optional): http://gmsh.info

### Install ###

Configure correctly the environment variables (EIGEN_DIR, EXPRTK_DIR, RTREE_DIR, etc.) and run [CMake](https://cmake.org).

- git clone https://github.com/bjaraujo/ENigMA.git
- cd ENigMA
- cmake-gui&
- setup source/build directory
- set paths to Exprtk  (dependencies folder), RTree (dependencies folder)
- configure/generate
- make

### Usage ###

- Check the tests and examples folders,
- Visit the [Wiki](https://github.com/bjaraujo/ENigMA/wiki),
- Run demo (releases).

### Quick Start ###

- Install python 3 32/64bit (if you don't have already), 
- Download the latest release,
- Install using pip:
```bash
pip3 install ENigMApy-VERSION-py3-none-DISTRO.whl
```
- Create file named "Example1.py"
- Go to wiki and copy example "Creating a mesh using python"
- Run command "python Example1.py"

### Examples ###

These examples show the mesh generation capability of the ENigMA library. 
The [pythonocc](https://github.com/tpaviot/pythonocc) library is used to build the CAD geometries and the ENigMA python wrapper to perform surface mesh generation. These examples use [Gmsh](http://gmsh.info/) for post-processing.

To reproduce these examples:

1. Download miniconda (python3.6 64bit): https://conda.io/miniconda.html

2. Install pythonocc: 
```bash
conda install -c conda-forge -c dlr-sc -c pythonocc -c oce pythonocc-core==0.18.1 python=3
```

3. Download ENigMA 0.1.2: [ENigMA_python3_64bit_0.1.2.0.zip](https://github.com/bjaraujo/ENigMA/releases/download/v0.1.2.0/ENigMA_python3_64bit_0.1.2.0.zip)

4. Download the examples: [ENigMA_python_examples.zip](https://github.com/bjaraujo/ENigMA/releases/download/v0.1.2.0/ENigMA_python_examples.zip)

#### Mesh Generation ####

##### A Cylinder #####

![cylinder](https://github.com/bjaraujo/ENigMA/blob/master/images/occ_01.png)

<details><summary>Code</summary>
<p>
        
```python
from OCC.BRepPrimAPI import BRepPrimAPI_MakeCylinder

import ENigMAocc
        
cylinder = BRepPrimAPI_MakeCylinder(3.0, 10.0).Shape() 

mesh = ENigMAocc.meshShape(cylinder, 0.5, 1E-3)
ENigMAocc.saveMeshFile(mesh, "occ_01.msh")
```

</p>
</details>

##### A Torus #####

![torus](https://github.com/bjaraujo/ENigMA/blob/master/images/occ_02.png)

<details><summary>Code</summary>
<p>

```python
from OCC.BRepPrimAPI import BRepPrimAPI_MakeTorus

import ENigMAocc
    
sphere = BRepPrimAPI_MakeTorus(4.0, 1.5).Shape() 

mesh = ENigMAocc.meshShape(sphere, 0.5, 1E-3)
ENigMAocc.saveMeshFile(mesh, "occ_02.msh")
```

</p>
</details>

##### A Cut-out Sphere #####

![sphere](https://github.com/bjaraujo/ENigMA/blob/master/images/occ_03.png)

<details><summary>Code</summary>
<p>

```python
from OCC.gp import gp_Pnt2d, gp_Pnt, gp_Vec, gp_Trsf
from OCC.BRepPrimAPI import BRepPrimAPI_MakeSphere
from OCC.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCC.BRepAlgoAPI import BRepAlgoAPI_Cut

import ENigMAocc
        
sphere1 = BRepPrimAPI_MakeSphere(3.0).Shape() 
sphere2 = BRepPrimAPI_MakeSphere(6.0).Shape() 

move = gp_Trsf()
move.SetTranslation(gp_Vec(4.0, 0.0, 0.0))
sphere2 = BRepBuilderAPI_Transform(sphere2, move, False).Shape()

shape = BRepAlgoAPI_Cut(sphere1, sphere2).Shape()

mesh = ENigMAocc.meshShape(shape, 0.5, 1E-3)
ENigMAocc.saveMeshFile(mesh, "occ_03.msh")
```

</p>
</details>

##### A Cut-out Shape #####

![cutout](https://github.com/bjaraujo/ENigMA/blob/master/images/occ_04a.png)

<details><summary>Code</summary>
<p>

```python
from OCC.gp import gp_Vec, gp_Trsf
from OCC.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
from OCC.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCC.BRepAlgoAPI import BRepAlgoAPI_Cut

import ENigMAocc
        
box = BRepPrimAPI_MakeBox(8.0, 8.0, 8.0).Shape()
cylinder = BRepPrimAPI_MakeCylinder(3.0, 6.0).Shape() 

move = gp_Trsf()
move.SetTranslation(gp_Vec(-4.0, -4.0, -4.0))
box = BRepBuilderAPI_Transform(box, move, False).Shape()

shape = BRepAlgoAPI_Cut(box, cylinder).Shape()

mesh = ENigMAocc.meshShape(shape, 0.5, 1E-3)
ENigMAocc.saveMeshFile(mesh, "occ_04a.msh")
```

</p>
</details>

##### A Fused Shape #####

![fused](https://github.com/bjaraujo/ENigMA/blob/master/images/occ_04b.png)

<details><summary>Code</summary>
<p>

```python
from OCC.gp import gp_Vec, gp_Trsf
from OCC.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
from OCC.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCC.BRepAlgoAPI import BRepAlgoAPI_Fuse

import ENigMAocc
        
box = BRepPrimAPI_MakeBox(8.0, 8.0, 8.0).Shape()
cylinder = BRepPrimAPI_MakeCylinder(3.0, 6.0).Shape() 

move = gp_Trsf()
move.SetTranslation(gp_Vec(-4.0, -4.0, -4.0))
box = BRepBuilderAPI_Transform(box, move, False).Shape()

shape = BRepAlgoAPI_Fuse(box, cylinder).Shape()

mesh = ENigMAocc.meshShape(shape, 0.5, 1E-3)
ENigMAocc.saveMeshFile(mesh, "occ_04b.msh")
```

</p>
</details>

##### More Cut-outs #####

![mcutout](https://github.com/bjaraujo/ENigMA/blob/master/images/occ_05.png)

<details><summary>Code</summary>
<p>

```python
from OCC.gp import gp_Vec, gp_Trsf
from OCC.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.BRepPrimAPI import BRepPrimAPI_MakeCylinder
from OCC.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCC.BRepAlgoAPI import BRepAlgoAPI_Cut

import ENigMAocc
        
box = BRepPrimAPI_MakeBox(100.0, 50.0, 4.0).Shape() 
cylinder1 = BRepPrimAPI_MakeCylinder(3.5, 10.0).Shape() 

move = gp_Trsf()

cut = box

for i in range(0, 9):
    for j in range(0, 4):
        move.SetTranslation(gp_Vec(i*10 + 10, j*10 + 10, 0.0))
        cylinder_copy = BRepBuilderAPI_Transform(cylinder1, move, True).Shape()
        cut = BRepAlgoAPI_Cut(cut, cylinder_copy).Shape()
    
cylinder2 = BRepPrimAPI_MakeCylinder(6.0, 10.0).Shape()  

move.SetTranslation(gp_Vec(0, 0, 0.0))
cylinder_copy = BRepBuilderAPI_Transform(cylinder2, move, True).Shape()
cut = BRepAlgoAPI_Cut(cut, cylinder_copy).Shape()

move.SetTranslation(gp_Vec(100, 0, 0.0))
cylinder_copy = BRepBuilderAPI_Transform(cylinder2, move, True).Shape()
cut = BRepAlgoAPI_Cut(cut, cylinder_copy).Shape()

move.SetTranslation(gp_Vec(100, 50, 0.0))
cylinder_copy = BRepBuilderAPI_Transform(cylinder2, move, True).Shape()
cut = BRepAlgoAPI_Cut(cut, cylinder_copy).Shape()

move.SetTranslation(gp_Vec(0, 50, 0.0))
cylinder_copy = BRepBuilderAPI_Transform(cylinder2, move, True).Shape()
shape = BRepAlgoAPI_Cut(cut, cylinder_copy).Shape()
        
mesh = ENigMAocc.meshShape(shape, 2.0, 1E-3)
ENigMAocc.saveMeshFile(mesh, "occ_05.msh")
```

</p>
</details>

##### A STEP File #####

![step](https://github.com/bjaraujo/ENigMA/blob/master/images/occ_06.png)

<details><summary>Code</summary>
<p>

```python
from OCC.STEPControl import STEPControl_Reader

import ENigMAocc
        
step_reader = STEPControl_Reader()
step_reader.ReadFile("lego.step")

step_reader.TransferRoot(1)
shape = step_reader.Shape(1)
    
mesh = ENigMAocc.meshShape(shape, 1.0, 1E-3)
ENigMAocc.saveMeshFile(mesh, "occ_06.msh")
```

</p>
</details>

##### The OCC Bottle #####

![occbottle](https://github.com/bjaraujo/ENigMA/blob/master/images/occ_07.png)

<details><summary>Code</summary>
<p>

```python

import math

from OCC.gp import gp_Pnt, gp_OX, gp_Vec, gp_Trsf, gp_DZ, gp_Ax2, gp_Ax3, gp_Pnt2d, gp_Dir2d, gp_Ax2d
from OCC.GC import GC_MakeArcOfCircle, GC_MakeSegment
from OCC.GCE2d import GCE2d_MakeSegment
from OCC.Geom import Geom_Plane, Geom_Curve, Geom_CylindricalSurface, Handle_Geom_Plane, Handle_Geom_Surface
from OCC.Geom2d import Geom2d_Ellipse, Geom2d_TrimmedCurve, Handle_Geom2d_Ellipse, Handle_Geom2d_Curve
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeFace, BRepBuilderAPI_Transform
from OCC.BRepPrimAPI import BRepPrimAPI_MakePrism, BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeBox
from OCC.BRepFilletAPI import BRepFilletAPI_MakeFillet
from OCC.BRepAlgoAPI import BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut
from OCC.BRepOffsetAPI import BRepOffsetAPI_MakeThickSolid, BRepOffsetAPI_ThruSections
from OCC.BRepLib import breplib
from OCC.BRep import BRep_Builder, BRep_Tool, BRep_Tool_Curve, BRep_Tool_Surface, BRep_Tool_CurveOnSurface
from OCC.BRepGProp import brepgprop_LinearProperties
from OCC.BRepAdaptor import BRepAdaptor_Curve
from OCC.TopExp import TopExp_Explorer
from OCC.TopTools import TopTools_ListOfShape
from OCC.TopExp import TopExp_Explorer
from OCC.TopoDS import topods, TopoDS_Compound, topods_Face, topods_Wire, topods_Edge
from OCC.TopAbs import TopAbs_FACE, TopAbs_WIRE, TopAbs_EDGE, TopAbs_REVERSED

import ENigMAocc

def faceIsPlane(face):
    hs = BRep_Tool_Surface(face)
    downcast_result = Handle_Geom_Plane.DownCast(hs)
    # The handle is null if downcast failed or is not possible, that is to say the face is not a plane
    if downcast_result.IsNull():
        return False
    else:
        return True

def geomPlaneFromFace(aFace):
    return Handle_Geom_Plane.DownCast(BRep_Tool_Surface(aFace)).GetObject()

height = 70
width = 50
thickness = 30

# The points we'll use to create the profile of the bottle's body
aPnt1 = gp_Pnt(-width / 2.0, 0, 0)
aPnt2 = gp_Pnt(-width / 2.0, -thickness / 4.0, 0)
aPnt3 = gp_Pnt(0, -thickness / 2.0, 0)
aPnt4 = gp_Pnt(width / 2.0, -thickness / 4.0, 0)
aPnt5 = gp_Pnt(width / 2.0, 0, 0)

aArcOfCircle = GC_MakeArcOfCircle(aPnt2, aPnt3, aPnt4)
aSegment1 = GC_MakeSegment(aPnt1, aPnt2)
aSegment2 = GC_MakeSegment(aPnt4, aPnt5)

# Could also construct the line edges directly using the points instead of the resulting line
aEdge1 = BRepBuilderAPI_MakeEdge(aSegment1.Value())
aEdge2 = BRepBuilderAPI_MakeEdge(aArcOfCircle.Value())
aEdge3 = BRepBuilderAPI_MakeEdge(aSegment2.Value())

# Create a wire out of the edges
aWire = BRepBuilderAPI_MakeWire(aEdge1.Edge(), aEdge2.Edge(), aEdge3.Edge())

# Quick way to specify the X axis
xAxis = gp_OX()

# Set up the mirror
aTrsf = gp_Trsf()
aTrsf.SetMirror(xAxis)

# Apply the mirror transformation
aBRespTrsf = BRepBuilderAPI_Transform(aWire.Wire(), aTrsf)

# Get the mirrored shape back out of the transformation and convert back to a wire
aMirroredShape = aBRespTrsf.Shape()

# A wire instead of a generic shape now
aMirroredWire = topods.Wire(aMirroredShape)

# Combine the two constituent wires
mkWire = BRepBuilderAPI_MakeWire()
mkWire.Add(aWire.Wire())
mkWire.Add(aMirroredWire)
myWireProfile = mkWire.Wire()

# The face that we'll sweep to make the prism
myFaceProfile = BRepBuilderAPI_MakeFace(myWireProfile)

# We want to sweep the face along the Z axis to the height
aPrismVec = gp_Vec(0, 0, height)
myBody = BRepPrimAPI_MakePrism(myFaceProfile.Face(), aPrismVec)

# Add fillets to all edges through the explorer
mkFillet = BRepFilletAPI_MakeFillet(myBody.Shape())
anEdgeExplorer = TopExp_Explorer(myBody.Shape(), TopAbs_EDGE)

while anEdgeExplorer.More():
    anEdge = topods.Edge(anEdgeExplorer.Current())
    mkFillet.Add(thickness / 12.0, anEdge)

    anEdgeExplorer.Next()

myBody = mkFillet

# Create the neck of the bottle
neckLocation = gp_Pnt(0, 0, height)
neckAxis = gp_DZ()
neckAx2 = gp_Ax2(neckLocation, neckAxis)

myNeckRadius = thickness / 4.0
myNeckHeight = height / 10.0

mkCylinder = BRepPrimAPI_MakeCylinder(neckAx2, myNeckRadius, myNeckHeight)

myBody = BRepAlgoAPI_Fuse(myBody.Shape(), mkCylinder.Shape())

# Our goal is to find the highest Z face and remove it
faceToRemove = None
zMax = -1

# We have to work our way through all the faces to find the highest Z face so we can remove it for the shell
aFaceExplorer = TopExp_Explorer(myBody.Shape(), TopAbs_FACE)
while aFaceExplorer.More():
    aFace = topods.Face(aFaceExplorer.Current())

    if faceIsPlane(aFace):
        aPlane = geomPlaneFromFace(aFace)

        # We want the highest Z face, so compare this to the previous faces
        aPnt = aPlane.Location()
        aZ = aPnt.Z()
        if aZ > zMax:
            zMax = aZ
            faceToRemove = aFace

    aFaceExplorer.Next()

facesToRemove = TopTools_ListOfShape()
facesToRemove.Append(faceToRemove)

myBody = BRepOffsetAPI_MakeThickSolid(myBody.Shape(), facesToRemove, -thickness / 50.0, 0.001)

# Set up our surfaces for the threading on the neck
neckAx2_Ax3 = gp_Ax3(neckLocation, gp_DZ())
aCyl1 = Geom_CylindricalSurface(neckAx2_Ax3, myNeckRadius * 0.99)
aCyl2 = Geom_CylindricalSurface(neckAx2_Ax3, myNeckRadius * 1.05)

# Set up the curves for the threads on the bottle's neck
aPnt = gp_Pnt2d(2.0 * math.pi, myNeckHeight / 2.0)
aDir = gp_Dir2d(2.0 * math.pi, myNeckHeight / 4.0)
anAx2d = gp_Ax2d(aPnt, aDir)

aMajor = 2.0 * math.pi
aMinor = myNeckHeight / 10.0

anEllipse1 = Geom2d_Ellipse(anAx2d, aMajor, aMinor)
anEllipse2 = Geom2d_Ellipse(anAx2d, aMajor, aMinor / 4.0)

anArc1 = Geom2d_TrimmedCurve(Handle_Geom2d_Ellipse(anEllipse1), 0, math.pi)
anArc2 = Geom2d_TrimmedCurve(Handle_Geom2d_Ellipse(anEllipse2), 0, math.pi)

anEllipsePnt1 = anEllipse1.Value(0)
anEllipsePnt2 = anEllipse1.Value(math.pi)

aSegment = GCE2d_MakeSegment(anEllipsePnt1, anEllipsePnt2)

# Build edges and wires for threading
anEdge1OnSurf1 = BRepBuilderAPI_MakeEdge(Handle_Geom2d_Curve(anArc1), Handle_Geom_Surface(aCyl1))
anEdge2OnSurf1 = BRepBuilderAPI_MakeEdge(aSegment.Value(), Handle_Geom_Surface(aCyl1))
anEdge1OnSurf2 = BRepBuilderAPI_MakeEdge(Handle_Geom2d_Curve(anArc2), Handle_Geom_Surface(aCyl2))
anEdge2OnSurf2 = BRepBuilderAPI_MakeEdge(aSegment.Value(), Handle_Geom_Surface(aCyl2))

threadingWire1 = BRepBuilderAPI_MakeWire(anEdge1OnSurf1.Edge(), anEdge2OnSurf1.Edge())
threadingWire2 = BRepBuilderAPI_MakeWire(anEdge1OnSurf2.Edge(), anEdge2OnSurf2.Edge())

# Compute the 3D representations of the edges/wires
breplib.BuildCurves3d(threadingWire1.Shape())
breplib.BuildCurves3d(threadingWire2.Shape())

# Create the surfaces of the threading
aTool = BRepOffsetAPI_ThruSections(True)
aTool.AddWire(threadingWire1.Wire())
aTool.AddWire(threadingWire2.Wire())
aTool.CheckCompatibility(False)
myThreading = aTool.Shape()

# Build the resulting compound
aRes = TopoDS_Compound()
aBuilder = BRep_Builder()
aBuilder.MakeCompound(aRes)
aBuilder.Add(aRes, myBody.Shape())
aBuilder.Add(aRes, myThreading)

shape = aRes

mesh = ENigMAocc.meshShape(shape, 1.5, 1E-3)
ENigMAocc.saveMeshFile(mesh, "occ_07.msh")
```

</p>
</details>

#### Structural Analysis ####

##### A Cantilever #####

A cantilever beam fixed at x=0 with an applied force of 1000 at x=1.

```
Max deflection (calculated)  = -0.021174580524850536
Max deflection (theoretical) = -0.021333333333333326
```

![cantilever](https://github.com/bjaraujo/ENigMA/blob/master/images/fem_01.png)

<details><summary>Code</summary>
<p>

```python
import math
import ENigMA

b = 0.05
h = 0.05

F = -1000.0
L = 1.0
E = 30E+9
I = b * h * h * h / 12

vertex1 = ENigMA.CMshNodeDouble(0.0, 0.0, 0.0)
vertex2 = ENigMA.CMshNodeDouble(L, 0.0, 0.0)
vertex3 = ENigMA.CMshNodeDouble(L, h, 0.0)
vertex4 = ENigMA.CMshNodeDouble(0.0, h, 0.0)

quadrilateral = ENigMA.CGeoQuadrilateralDouble()

quadrilateral.addVertex(vertex1)
quadrilateral.addVertex(vertex2)
quadrilateral.addVertex(vertex3)
quadrilateral.addVertex(vertex4)

basicMesher = ENigMA.CMshBasicMesherDouble()

basicMesher.generate(quadrilateral, 400, 20, True);

surfaceMesh = basicMesher.mesh();

for i in range(0, surfaceMesh.nbElements()):
    elementId = surfaceMesh.elementId(i)
    surfaceMesh.element(elementId).setThickness(b);

# Material
material = ENigMA.CMatMaterialDouble()

material.addProperty(ENigMA.PT_ELASTIC_MODULUS, E)
material.addProperty(ENigMA.PT_POISSON_COEFFICIENT, 0.3)

# Stress field
u = ENigMA.CPdeFieldDouble()

u.setMesh(surfaceMesh)
u.setMaterial(material)
u.setSimulationType(ENigMA.ST_STRUCTURAL)
u.setDiscretMethod(ENigMA.DM_FEM)
u.setDiscretOrder(ENigMA.DO_LINEAR)
u.setDiscretLocation(ENigMA.DL_NODE)
u.setNbDofs(2)

index = 0;

for i in range(0, surfaceMesh.nbNodes()):
    nodeId = surfaceMesh.nodeId(i)
    node = surfaceMesh.node(nodeId)
    
    if (math.fabs(node.x() - 0.0) < 1E-6):
        u.setFixedValue(surfaceMesh.nodeIndex(nodeId), 0, 0.0)
        u.setFixedValue(surfaceMesh.nodeIndex(nodeId), 1, 0.0)
        
    if (math.fabs(node.x() - L) < 1E-6 and math.fabs(node.y() - h) < 1E-6):
        index = i;
        u.setSource(surfaceMesh.nodeIndex(nodeId), 1, F)

sleSystem = ENigMA.laplacian(u)
sleSystem.setRhs(0)
        
pdeEquation = ENigMA.CPdeEquationDouble(sleSystem)

pdeEquation.setSources(u);

pdeEquation.setElimination(u)

pdeEquation.solve(u)

posGmsh = ENigMA.CPosGmshDouble()
posGmsh.save(u, "fem_01.msh", "Structural")

# Theoretical displacement
deflection = (F * L * L * L) / (3 * E * I)

print('Max deflection (calculated)  = ' + str(u.value(index * 2 + 1)))
print('Max deflection (theoretical) = ' + str(deflection))
```

</p>
</details>

#### Thermal Analysis ####

##### A Plate #####

Two-dimensional steady-state thermal analysis in a rectangular plate heated at the top (T = 1), other sides kept at 0.

```
Temperature at point (0.5, 0.5) (calculated)  = 0.25
Temperature at point (0.5, 0.5) (theoretical) = 0.25
```

![plate](https://github.com/bjaraujo/ENigMA/blob/master/images/fem_02.png)

<details><summary>Code</summary>
<p>

```python
import math
import ENigMA

vertex1 = ENigMA.CMshNodeDouble(0.0, 0.0, 0.0)
vertex2 = ENigMA.CMshNodeDouble(1.0, 0.0, 0.0)
vertex3 = ENigMA.CMshNodeDouble(1.0, 1.0, 0.0)
vertex4 = ENigMA.CMshNodeDouble(0.0, 1.0, 0.0)

quadrilateral = ENigMA.CGeoQuadrilateralDouble()

quadrilateral.addVertex(vertex1)
quadrilateral.addVertex(vertex2)
quadrilateral.addVertex(vertex3)
quadrilateral.addVertex(vertex4)

basicMesher = ENigMA.CMshBasicMesherDouble()

basicMesher.generate(quadrilateral, 16, 16, True);

surfaceMesh = basicMesher.mesh();

# Temperature field
u = ENigMA.CPdeFieldDouble()

u.setMesh(surfaceMesh)
u.setSimulationType(ENigMA.ST_THERMAL)
u.setDiscretMethod(ENigMA.DM_FEM)
u.setDiscretOrder(ENigMA.DO_LINEAR)
u.setDiscretLocation(ENigMA.DL_NODE)
u.setNbDofs(1)

index = 0;

for i in range(0, surfaceMesh.nbNodes()):
    nodeId = surfaceMesh.nodeId(i)
    node = surfaceMesh.node(nodeId)
    
    if (math.fabs(node.x() - 0.0) < 1E-6):
        u.setFixedValue(surfaceMesh.nodeIndex(nodeId), 0.0)

    if (math.fabs(node.x() - 1.0) < 1E-6):
        u.setFixedValue(surfaceMesh.nodeIndex(nodeId), 0.0)

    if (math.fabs(node.y() - 0.0) < 1E-6):
        u.setFixedValue(surfaceMesh.nodeIndex(nodeId), 0.0)
        
    if (math.fabs(node.y() - 1.0) < 1E-6):
        u.setFixedValue(surfaceMesh.nodeIndex(nodeId), 1.0)

sleSystem = ENigMA.laplacian(u)
sleSystem.setRhs(0)
        
pdeEquation = ENigMA.CPdeEquationDouble(sleSystem)

pdeEquation.setSources(u);

pdeEquation.setElimination(u)

pdeEquation.solve(u)

posGmsh = ENigMA.CPosGmshDouble()
posGmsh.save(u, "fem_02.msh", "Thermal")

for i in range(0, surfaceMesh.nbNodes()):
    nodeId = surfaceMesh.nodeId(i)
    node = surfaceMesh.node(nodeId)
    
    if (math.fabs(node.x() - 0.5) < 1E-6 and math.fabs(node.y() - 0.5) < 1E-6):
        index = surfaceMesh.nodeIndex(nodeId)

anaTemperature = ENigMA.CAnaTemperatureDouble()

ut = anaTemperature.steadyStateHeatConduction2D(0.5, 0.5)

print('Temperature at point (0.5, 0.5) (theoretical) = ' + str(ut))
print('Temperature at point (0.5, 0.5) (calculated)  = ' + str(u.value(index)))
```

</p>
</details>

#### Flow Analysis ####

Lid-driven cavity (Re = 1000).

![liddriven](https://github.com/bjaraujo/ENigMA/blob/master/images/fvm_01.png)

<details><summary>Code</summary>
<p>

```python
import math
import ENigMA

vertex1 = ENigMA.CGeoCoordinateDouble(+0.00, +0.00, -0.1)
vertex2 = ENigMA.CGeoCoordinateDouble(+1.00, +0.00, -0.1)
vertex3 = ENigMA.CGeoCoordinateDouble(+1.00, +1.00, -0.1)
vertex4 = ENigMA.CGeoCoordinateDouble(+0.00, +1.00, -0.1)
vertex5 = ENigMA.CGeoCoordinateDouble(+0.00, +0.00, +0.1)
vertex6 = ENigMA.CGeoCoordinateDouble(+1.00, +0.00, +0.1)
vertex7 = ENigMA.CGeoCoordinateDouble(+1.00, +1.00, +0.1)
vertex8 = ENigMA.CGeoCoordinateDouble(+0.00, +1.00, +0.1)

hexahedron = ENigMA.CGeoHexahedronDouble()

hexahedron.addVertex(vertex1)
hexahedron.addVertex(vertex2)
hexahedron.addVertex(vertex3)
hexahedron.addVertex(vertex4)
hexahedron.addVertex(vertex5)
hexahedron.addVertex(vertex6)
hexahedron.addVertex(vertex7)
hexahedron.addVertex(vertex8)

mesher = ENigMA.CMshBasicMesherDouble()
mesher.generate(hexahedron, 40, 40, 1, False)
mesher.mesh().generateFaces(1E-3)
mesher.mesh().calculateFaceCentroid()
mesher.mesh().calculateElementCentroid()

volumeMesh = mesher.mesh()
fvmMesh = ENigMA.CFvmMeshDouble(volumeMesh)

U = 1.0    # lid velocity
mu = 0.001 # dynamic viscosity
rho = 1.0  # density

nu = mu / rho # kinematic viscosity 

pisoSolver = ENigMA.CFvmPisoSolverDouble(fvmMesh)

pisoSolver.setGravity(0.0, 0.0, 0.0)
pisoSolver.setMaterialProperties(rho, mu)

faceIds = ENigMA.StdVectorInt()

faceIds.clear()
for i in range(0, fvmMesh.nbFaces()):

    faceId = fvmMesh.faceId(i)
    face = fvmMesh.face(faceId)
    face.calculateCentroid()

    if (face.centroid().x() == 0.0 or
        face.centroid().x() == 1.0 or
        face.centroid().y() == 0.0 or
        face.centroid().y() == 1.0):
        faceIds.push_back(faceId)

pisoSolver.setBoundaryVelocity(faceIds, ENigMA.BT_WALL_NO_SLIP, 0.0, 0.0, 0.0)
pisoSolver.setBoundaryPressure(faceIds, ENigMA.BT_WALL_NO_SLIP, 0.0)

faceIds.clear()
for i in range(0, fvmMesh.nbFaces()):

    faceId = fvmMesh.faceId(i)
    face = fvmMesh.face(faceId)
    face.calculateCentroid()

    if (face.centroid().y() == 1.0):
        faceIds.push_back(faceId)

pisoSolver.setBoundaryVelocity(faceIds, ENigMA.BT_WALL_NO_SLIP, U, 0.0, 0.0)

# Courant < 1
dt = U / 40			

iter = 100

# Flow in a rectangle
for ii in range(0, iter):
    print('Iter = ' + str(ii + 1) + ' of ' + str(iter))
    pisoSolver.iterate(dt)

u = ENigMA.CPdeFieldDouble()
u.setMesh(volumeMesh)
u.setSimulationType(ENigMA.ST_FLOW)
u.setDiscretMethod(ENigMA.DM_FVM)
u.setDiscretOrder(ENigMA.DO_LINEAR)
u.setDiscretLocation(ENigMA.DL_ELEMENT_CENTER)
u.setNbDofs(2)
u.setSize(fvmMesh.nbControlVolumes() * 2)
   
for i in range(0, fvmMesh.nbControlVolumes()):

    controlVolumeId = fvmMesh.controlVolumeId(i)

    u.setValue(i * 2 + 0, pisoSolver.u(controlVolumeId))
    u.setValue(i * 2 + 1, pisoSolver.v(controlVolumeId))

posGmsh = ENigMA.CPosGmshDouble()
posGmsh.save(u, "fvm_01.msh", "Flow")
```

</p>
</details>




