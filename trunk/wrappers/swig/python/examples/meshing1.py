
from timeit import default_timer as timer
#import subprocess

from ENigMA import ENigMA

edgeMesh = ENigMA.CMshMesh()

node1 = ENigMA.CMshNode(10, 10, 0)
node2 = ENigMA.CMshNode(400, 10, 0)
node3 = ENigMA.CMshNode(400, 400, 0)
node4 = ENigMA.CMshNode(10, 400, 0)

edgeMesh.addNode(1, node1)
edgeMesh.addNode(2, node2)
edgeMesh.addNode(3, node3)
edgeMesh.addNode(4, node4)

element1 = ENigMA.CMshElement(ENigMA.ET_BEAM)
element1.addNodeId(1)
element1.addNodeId(2)

element2 = ENigMA.CMshElement(ENigMA.ET_BEAM)
element2.addNodeId(2)
element2.addNodeId(3)

element3 = ENigMA.CMshElement(ENigMA.ET_BEAM)
element3.addNodeId(3)
element3.addNodeId(4)

element4 = ENigMA.CMshElement(ENigMA.ET_BEAM)
element4.addNodeId(4)
element4.addNodeId(1)

edgeMesh.addElement(1, element1)
edgeMesh.addElement(2, element2)
edgeMesh.addElement(3, element3)
edgeMesh.addElement(4, element4)

triangleMesher = ENigMA.CMshTriangleMesher()
interiorPoints = ENigMA.StdVectorCGeoCoordinate()

meshSize = 20

edgeMesh.generateFaces(1E-3)
triangleMesher.remesh(edgeMesh, meshSize)

start = timer()
triangleMesher.generate(edgeMesh, 9999, interiorPoints, meshSize, meshSize * 0.1, meshSize * 10.0, 1E-3)
end = timer()

print("Elapsed time: " + str(end - start))

surfaceMesh = triangleMesher.mesh()

triangleMesher.flipEdges(surfaceMesh)
triangleMesher.relaxNodes(surfaceMesh)
triangleMesher.flipEdges(surfaceMesh)

print(surfaceMesh.nbNodes())
print(surfaceMesh.nbElements())

pdeField = ENigMA.CPdeField()

pdeField.setMesh(surfaceMesh)

posGmsh = ENigMA.CPosGmsh()

posGmsh.save(pdeField, "mesh1.msh", "tris")

#subprocess.call(['gmsh.exe', 'mesh1.msh'])
