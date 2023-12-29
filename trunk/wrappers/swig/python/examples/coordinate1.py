
from ENigMA import ENigMA

a = ENigMA.CGeoCoordinate(1, 2, 3)

b = ENigMA.CGeoCoordinateSystem()

# set value
b(0, 0, 1)
b(1, 1, 1)
b(2, 2, 1)

# get value
print(b(0, 0))

a.transform(b)

print(a.x())
print(a.y())
print(a.z())
