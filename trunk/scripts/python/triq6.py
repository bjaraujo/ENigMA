from sympy.solvers import solve
from sympy import Symbol

u1 = Symbol('u1')
u2 = Symbol('u2')
u3 = Symbol('u3')
u4 = Symbol('u4')
u5 = Symbol('u5')
u6 = Symbol('u6')

alpha1 = Symbol('alpha1')
alpha2 = Symbol('alpha2')
alpha3 = Symbol('alpha3')
alpha4 = Symbol('alpha4')
alpha5 = Symbol('alpha5')
alpha6 = Symbol('alpha6')

x1 = Symbol('x1')
x2 = Symbol('x2')
x3 = Symbol('x3')
x4 = Symbol('x4')
x5 = Symbol('x5')
x6 = Symbol('x6')

y1 = Symbol('y1')
y2 = Symbol('y2')
y3 = Symbol('y3')
y4 = Symbol('y4')
y5 = Symbol('y5')
y6 = Symbol('y6')

f1 = alpha1 + alpha2*x1 + alpha3*y1 + alpha4*x1*y1 + alpha5*x1**2 + alpha6*y1**2 - u1
f2 = alpha1 + alpha2*x2 + alpha3*y2 + alpha4*x2*y2 + alpha5*x2**2 + alpha6*y2**2 - u2
f3 = alpha1 + alpha2*x3 + alpha3*y3 + alpha4*x3*y3 + alpha5*x3**2 + alpha6*y3**2 - u3
f4 = alpha1 + alpha2*x4 + alpha3*y4 + alpha4*x4*y4 + alpha5*x4**2 + alpha6*y4**2 - u4
f5 = alpha1 + alpha2*x5 + alpha3*y5 + alpha4*x5*y5 + alpha5*x5**2 + alpha6*y5**2 - u5
f6 = alpha1 + alpha2*x6 + alpha3*y6 + alpha4*x6*y6 + alpha5*x6**2 + alpha6*y6**2 - u6

solution = solve([f1, f2, f3, f4, f5, f6], [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6])
print(solution)
