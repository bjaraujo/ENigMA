import subprocess

import ENigMA
 
vertex1 = ENigMA.CGeoCoordinateDouble(0, 0, 0)
vertex2 = ENigMA.CGeoCoordinateDouble(1, 0, 0)
vertex3 = ENigMA.CGeoCoordinateDouble(1, 1, 0)
vertex4 = ENigMA.CGeoCoordinateDouble(0, 1, 0)

quadrilateral = ENigMA.CGeoQuadrilateralDouble()

quadrilateral.addVertex(vertex1)
quadrilateral.addVertex(vertex2)
quadrilateral.addVertex(vertex3)
quadrilateral.addVertex(vertex4)

basicMesher = ENigMA.CMshBasicMesherDouble()

basicMesher.generate(quadrilateral, 15, 15, True)

surfaceMesh = basicMesher.mesh()

pdeField = ENigMA.CPdeFieldDouble()

pdeField.setMesh(surfaceMesh)
pdeField.setDiscretMethod(ENigMA.DM_FEM)
pdeField.setDiscretOrder(ENigMA.DO_LINEAR)
pdeField.setDiscretLocation(ENigMA.DL_NODE)
pdeField.setSimulationType(ENigMA.ST_THERMAL)
pdeField.setNbDofs(1)

for i in range(0, surfaceMesh.nbNodes()):
    
    nodeId = surfaceMesh.nodeId(i)
    node = surfaceMesh.node(nodeId)
    
    if (abs(node.x() - 0) < 1E-6):
        pdeField.setFixedValue(i, 0.0)

    if (abs(node.x() - 1) < 1E-6):
        pdeField.setFixedValue(i, 0.0)

    if (abs(node.y() - 0) < 1E-6):
        pdeField.setFixedValue(i, 0.0)

    if (abs(node.y() - 1) < 1E-6):
        pdeField.setFixedValue(i, 1.0)
        
pdeEquation = ENigMA.CPdeEquationDouble(ENigMA.laplacian(pdeField) == 0)

#pdeEquation.setElimination(pdeField);

pdeEquation.solve(pdeField);

posGmsh = ENigMA.CPosGmshDouble()

posGmsh.save(pdeField, "pde1.msh", "tris");

#subprocess.call(['gmsh.exe', 'pde1.msh'])
