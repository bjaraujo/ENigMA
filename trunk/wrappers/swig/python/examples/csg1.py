
import subprocess

from ENigMA import ENigMA
 
center = ENigMA.CGeoCoordinate(0,0,0)
radius = ENigMA.CGeoVector(1,1,1)

a = ENigMA.CCsgCubeDouble.create(center, radius)
b = ENigMA.CCsgSphereDouble.create(center, 1.35, 16, 12)

vertex_x1 = ENigMA.CGeoCoordinate(-2,0,0)
vertex_x2 = ENigMA.CGeoCoordinate(+2,0,0)

c = ENigMA.CCsgCylinderDouble.create(vertex_x1, vertex_x2, 0.7, 16)

vertex_y1 = ENigMA.CGeoCoordinate(0,-2,0)
vertex_y2 = ENigMA.CGeoCoordinate(0,+2,0)

d = ENigMA.CCsgCylinderDouble.create(vertex_y1, vertex_y2, 0.7, 16)

vertex_z1 = ENigMA.CGeoCoordinate(0,0,-2)
vertex_z2 = ENigMA.CGeoCoordinate(0,0,+2)

e = ENigMA.CCsgCylinderDouble.create(vertex_z1, vertex_z2, 0.7, 16)

tol = 1E-6

a.setTolerance(tol)
b.setTolerance(tol)
c.setTolerance(tol)
d.setTolerance(tol)
e.setTolerance(tol)

f = a.intersect(b).subtract(c.add(d).add(e))

print(f.nbPolygons())

f.toStl().save("f.stl")

subprocess.call(['gmsh.exe', 'f.stl'])
