
import ENigMA

a = ENigMA.CGeoCoordinateDouble(0, 0, 0)
b = ENigMA.CGeoCoordinateDouble(1, 0, 0)

l1 = ENigMA.CGeoLineDouble(a, b)

l1.calculateLength()

print '> l1.length'
print l1.length()

n = ENigMA.CGeoNormalDouble(1, 0, 0)
p = ENigMA.CGeoPlaneDouble(n, 0.2)

l2 = l1.clip(p)

l2.calculateLength()

print '> l2.length'
print l2.length()

