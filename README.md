# ENigMA - Extended Numerical Multiphysics Analysis #

ENigMA is an object-oriented C++ template library which aim is to provide multi-physics simulation in a multi-domain environment. The code implements several numerical methods such as Finite Volume Methods (FVM), Finite Difference Methods (FDM), Finite Element Methods (FEM), Boundary Element Methods (BEM), etc. for numerical approximation of Partial Differential Equations (PDE) in each domain. It can be used for three-dimensional flow, thermal and structural analysis and also provides a series of automatic mesh generators (triangular, constrained tetrahedral, etc). ENigMA provides 3D geometry classes for intersection and clipping operations and implements RTree, Octree and HashGrid methods for spatial searching. It was developed to be cross-platform using STL, Eigen (for vectors and matrices) and exprtk (for math expressions). The SWIG tool is used to expose ENigMA's classes and methods to other languages such as Python, C#, etc. It also uses Gtest for unit testing (> 160 unit tests), CMake is used for cross-platform building and Git is used for source code management.

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

### Dependencies ###

Core:
- Eigen: http://eigen.tuxfamily.org/index.php?title=Main_Page
- Exprtk: https://github.com/ArashPartow/exprtk or dependencies folder
- ViennaCL (optional): http://viennacl.sourceforge.net

Tests:
- Gtest: https://github.com/google/googletest

Examples:
- Qt: https://www.qt.io
- VTK: http://www.vtk.org
- OpenGL/GLUT: http://freeglut.sourceforge.net

### Install ###

Configure correct environment variables (GTEST_DIR, EXPRTK_DIR, EIGEN_DIR, VIENNACL_DIR, etc.) and run CMake.

### Usage ###

Check the tests and examples folders.
