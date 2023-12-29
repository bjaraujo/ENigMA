
from ENigMA import ENigMA

a = ENigMA.CGeoCoordinate(0, 0, 0)
b = ENigMA.CGeoCoordinate(1, 0, 0)

l1 = ENigMA.CGeoLine(a, b)

l1.calculateLength()

print('> l1.length')
print(l1.length())

n = ENigMA.CGeoNormal(1, 0, 0)
p = ENigMA.CGeoPlane(n, 0.2)

l2 = l1.clip(p)

l2.calculateLength()

print('> l2.length')
print(l2.length())

