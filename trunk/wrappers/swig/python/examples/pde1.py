import subprocess

from ENigMA import ENigMA
 
vertex1 = ENigMA.CGeoCoordinate(0, 0, 0)
vertex2 = ENigMA.CGeoCoordinate(1, 0, 0)
vertex3 = ENigMA.CGeoCoordinate(1, 1, 0)
vertex4 = ENigMA.CGeoCoordinate(0, 1, 0)

quadrilateral = ENigMA.CGeoQuadrilateral()

quadrilateral.addVertex(vertex1)
quadrilateral.addVertex(vertex2)
quadrilateral.addVertex(vertex3)
quadrilateral.addVertex(vertex4)

basicMesher = ENigMA.CMshBasicMesher()

basicMesher.generate(quadrilateral, 15, 15, True)

surfaceMesh = basicMesher.mesh()

T = ENigMA.CPdeField()

T.setMesh(surfaceMesh)
T.setDiscretMethod(ENigMA.DM_FEM)
T.setDiscretOrder(ENigMA.DO_LINEAR)
T.setDiscretLocation(ENigMA.DL_NODE)
T.setSimulationType(ENigMA.ST_THERMAL)
T.setNbDofs(1)

T.setSize(surfaceMesh.nbNodes())

for i in range(0, surfaceMesh.nbNodes()):
    
    nodeId = surfaceMesh.nodeId(i)
    node = surfaceMesh.node(nodeId)
    
    if (abs(node.x() - 0) < 1E-3):
        T.setFixedValue(i, 0.0)

    if (abs(node.x() - 1) < 1E-3):
        T.setFixedValue(i, 0.0)

    if (abs(node.y() - 0) < 1E-3):
        T.setFixedValue(i, 0.0)

    if (abs(node.y() - 1) < 1E-3):
        T.setFixedValue(i, 1.0)
        
    T.setValue(i, 0.0)
    
pdeEquation = ENigMA.CPdeEquation(ENigMA.laplacian(T))

pdeEquation.setElimination(T)

pdeEquation.solve(T)

posGmsh = ENigMA.CPosGmsh()

posGmsh.save(T, "pde1.msh", "tris")

#subprocess.call(['gmsh.exe', 'pde1.msh'])
