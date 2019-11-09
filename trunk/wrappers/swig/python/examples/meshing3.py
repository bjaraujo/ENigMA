
from timeit import default_timer as timer
#import subprocess

import ENigMA

d = 0.125

nu = 12
nv = 12
nw = 12

vertex1 = ENigMA.CGeoCoordinateDouble(0.0, 0.0, 0.0)
vertex2 = ENigMA.CGeoCoordinateDouble(nu * d, 0.0, 0.0)
vertex3 = ENigMA.CGeoCoordinateDouble(nu * d, nv * d, 0.0)
vertex4 = ENigMA.CGeoCoordinateDouble(0.0, nv * d, 0.0)
vertex5 = ENigMA.CGeoCoordinateDouble(0.0, 0.0, nw * d)
vertex6 = ENigMA.CGeoCoordinateDouble(nu * d, 0.0, nw * d)
vertex7 = ENigMA.CGeoCoordinateDouble(nu * d, nv * d, nw * d)
vertex8 = ENigMA.CGeoCoordinateDouble(0.0, nv * d, nw * d)

hexahedron = ENigMA.CGeoHexahedronDouble()

hexahedron.addVertex(vertex1)
hexahedron.addVertex(vertex2)
hexahedron.addVertex(vertex3)
hexahedron.addVertex(vertex4)
hexahedron.addVertex(vertex5)
hexahedron.addVertex(vertex6)
hexahedron.addVertex(vertex7)
hexahedron.addVertex(vertex8)

basicMesher = ENigMA.CMshBasicMesherDouble()

basicMesher.generate(hexahedron, nu, nv, nw, True)

surfaceMesh = basicMesher.mesh().extractBoundary(1E-3)

surfaceMesh.generateFaces(1E-5)

tetrahedronMesher = ENigMA.CMshTetrahedronMesherDouble()

start = timer()
tetrahedronMesher.generate(surfaceMesh, 99999, d, 0.1, 1E-3);
end = timer()

print("Elapsed time: " + str(end - start))

volumeMesh = tetrahedronMesher.mesh()

print(volumeMesh.nbNodes())
print(volumeMesh.nbElements())

pdeField = ENigMA.CPdeFieldDouble()

pdeField.setMesh(volumeMesh)

posGmsh = ENigMA.CPosGmshDouble()

posGmsh.save(pdeField, "mesh3.msh", "tetras");

#subprocess.call(['gmsh.exe', 'mesh2.msh'])
