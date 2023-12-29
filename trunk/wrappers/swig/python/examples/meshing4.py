
from ENigMA import ENigMA

pdeField = ENigMA.CPdeField()
posGmsh = ENigMA.CPosGmsh()

posGmsh.load(pdeField, "surf.msh")

surfaceMesh = pdeField.mesh()
surfaceMesh.generateFaces(0.02)

tetrahedronMesher = ENigMA.CMshTetrahedronMesher()
interiorPoints = ENigMA.StdVectorCGeoCoordinate()

try:
    tetrahedronMesher.generate(surfaceMesh, 99999, interiorPoints, 2.6, 2.6, 2.6, 0.02)
except:
    pass
    
print(tetrahedronMesher.mesh().nbNodes())
print(tetrahedronMesher.mesh().nbElements())

pdeField.setMesh(tetrahedronMesher.mesh())
posGmsh.save(pdeField, "vol.msh", "tetras")


