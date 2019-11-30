
from ENigMA import ENigMA

a = ENigMA.CGeoVectorDouble(1, 0, 0)
b = ENigMA.CGeoVectorDouble(0, 1, 0)

print '> a'
print a.x()
print a.y()
print a.z()

print '> b'
print b.x()
print b.y()
print b.z()

print '> c = a + b'
c = a + b
print c.x()
print c.y()
print c.z()

print '> d = a - b'
d = a - b
print d.x()
print d.y()
print d.z()

print '> e = a * 3'
e = a * 3
print e.x()
print e.y()
print e.z()

print '> f = a / 3'
f = a / 3
print f.x()
print f.y()
print f.z()

print '> a += b'
a += b
print a.x()
print a.y()
print a.z()

print a.angle(b)
