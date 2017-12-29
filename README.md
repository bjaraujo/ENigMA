[![License: GPL v3](https://img.shields.io/badge/license-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://travis-ci.org/bjaraujo/ENigMA.svg?branch=master)](https://travis-ci.org/bjaraujo/ENigMA/build)


# ENigMA - Extended Numerical Multiphysics Analysis #

ENigMA is an object-oriented C++ template library which aim is to provide multi-physics simulation in a multi-domain environment. The code implements several numerical methods such as Finite Volume Methods (FVM), Finite Difference Methods (FDM), Finite Element Methods (FEM), Boundary Element Methods (BEM), Smoothed Particle Hydrodynamics (SPH), etc. for numerical approximation of Partial Differential Equations (PDE) in each domain. 

ENigMA provides classes for mesh generation (triangular, block, constrained tetrahedral, etc), intersection and clipping operations and implements R-tree, octree and hashgrid methods for spatial searching. It can be used for three-dimensional flow, thermal and structural analysis.

It was developed to be cross-platform using STL, Eigen (for vectors and matrices) and exprtk (for math expressions). The SWIG tool is used to expose ENigMA's classes and methods to other languages such as Python, C#, etc. It also uses Gtest for unit testing (> 160 unit tests), CMake is used for cross-platform building and Git is used for source code management.

![mesh](https://github.com/bjaraujo/ENigMA/blob/master/images/mesh.png)

ENigMA implements several operators such as divergence, laplacian, gradient, etc. so, for example, a heat equation such as: 

![equation](https://github.com/bjaraujo/ENigMA/blob/master/images/equation.gif)

can be represented in ENigMA as:

```cpp
CPdeEquation<double> aPdeEquation(1.0 / dt * ddt<double>(u) - alpha * laplacian<double>(u) = 0);
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

Note: the constrained advancing-front tetrahedral mesher is part of ENigMA+. Will be added to ENigMA source code base when this project reaches 100 stars. 

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
- OpenCASCADE (optional): https://www.opencascade.com/
- Gmsh (optional): http://gmsh.info/

### Install ###

Configure correctly the environment variables (EIGEN_DIR, EXPRTK_DIR, RTREE_DIR, etc.) and run [CMake](https://cmake.org).

- git clone https://github.com/bjaraujo/ENigMA.git
- cd ENigMA
- cmake-gui&
- setup source/build directory
- download Eigen and GTest
- set paths to Eigen, GTest, Exprtk  (dependencies folder), RTree (dependencies folder)
- configure/generate
- make

### Usage ###

- Check the tests and examples folders,
- Visit the [Wiki](https://github.com/bjaraujo/ENigMA/wiki),
- Run demo (releases).

### Quick Start ###

- Install python 3 64bit (if you don't have already), 
- Download the latest release (python wrapper),
- Create file named "Example1.py"
- Go to wiki and copy example "Creating a mesh using python"
- Run command "python Example1.py"

### Examples ###

#### Mesh Generation ####

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

##### A Cylinder #####

![cylinder](https://github.com/bjaraujo/ENigMA/blob/master/images/occ_01.png)
```python
from OCC.BRepPrimAPI import BRepPrimAPI_MakeCylinder

import ENigMAocc
        
cylinder = BRepPrimAPI_MakeCylinder(3.0, 10.0).Shape() 

mesh = ENigMAocc.meshShape(cylinder, 0.5, 1E-3)
ENigMAocc.saveMeshFile(mesh, "occ_01.msh")
```

##### A Torus #####

![torus](https://github.com/bjaraujo/ENigMA/blob/master/images/occ_02.png)
```python
from OCC.BRepPrimAPI import BRepPrimAPI_MakeTorus

import ENigMAocc
    
sphere = BRepPrimAPI_MakeTorus(4.0, 1.5).Shape() 

mesh = ENigMAocc.meshShape(sphere, 0.5, 1E-3)
ENigMAocc.saveMeshFile(mesh, "occ_02.msh")
```

##### A Cut-out Sphere #####

![sphere](https://github.com/bjaraujo/ENigMA/blob/master/images/occ_03.png)
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

##### A Cut-out Shape #####

![cutout](https://github.com/bjaraujo/ENigMA/blob/master/images/occ_04a.png)
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

##### A Fused Shape #####

![fused](https://github.com/bjaraujo/ENigMA/blob/master/images/occ_04b.png)
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

##### More Cut-outs #####

![mcutout](https://github.com/bjaraujo/ENigMA/blob/master/images/occ_05.png)
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

##### A STEP File #####

![step](https://github.com/bjaraujo/ENigMA/blob/master/images/occ_06.png)
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
