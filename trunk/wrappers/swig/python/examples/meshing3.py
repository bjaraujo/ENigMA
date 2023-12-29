
from timeit import default_timer as timer
#import subprocess

from ENigMA import ENigMA

d = 0.125

nu = 10
nv = 10
nw = 10

vertex1 = ENigMA.CGeoCoordinate(0.0, 0.0, 0.0)
vertex2 = ENigMA.CGeoCoordinate(nu * d, 0.0, 0.0)
vertex3 = ENigMA.CGeoCoordinate(nu * d, nv * d, 0.0)
vertex4 = ENigMA.CGeoCoordinate(0.0, nv * d, 0.0)
vertex5 = ENigMA.CGeoCoordinate(0.0, 0.0, nw * d)
vertex6 = ENigMA.CGeoCoordinate(nu * d, 0.0, nw * d)
vertex7 = ENigMA.CGeoCoordinate(nu * d, nv * d, nw * d)
vertex8 = ENigMA.CGeoCoordinate(0.0, nv * d, nw * d)

hexahedron = ENigMA.CGeoHexahedron()

hexahedron.addVertex(vertex1)
hexahedron.addVertex(vertex2)
hexahedron.addVertex(vertex3)
hexahedron.addVertex(vertex4)
hexahedron.addVertex(vertex5)
hexahedron.addVertex(vertex6)
hexahedron.addVertex(vertex7)
hexahedron.addVertex(vertex8)

basicMesher = ENigMA.CMshBasicMesher()

basicMesher.generate(hexahedron, nu, nv, nw, True)

surfaceMesh = basicMesher.mesh().extractBoundary(1E-3)

surfaceMesh.generateFaces(1E-5)

tetrahedronMesher = ENigMA.CMshTetrahedronMesher()
interiorPoints = ENigMA.StdVectorCGeoCoordinate()

start = timer()
tetrahedronMesher.generate(surfaceMesh, 99999, interiorPoints, d, d * 0.5, d * 2.0, 1E-6)
end = timer()

print("Elapsed time: " + str(end - start))

volumeMesh = tetrahedronMesher.mesh()

print(volumeMesh.nbNodes())
print(volumeMesh.nbElements())

pdeField = ENigMA.CPdeField()

pdeField.setMesh(volumeMesh)

posGmsh = ENigMA.CPosGmsh()

posGmsh.save(pdeField, "mesh3.msh", "tetras")

#subprocess.call(['gmsh.exe', 'mesh2.msh'])
