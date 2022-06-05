from sympy.solvers import solve
from sympy import Symbol

u1 = Symbol('u1')
u2 = Symbol('u2')
u3 = Symbol('u3')
u4 = Symbol('u4')

alpha1 = Symbol('alpha1')
alpha2 = Symbol('alpha2')
alpha3 = Symbol('alpha3')
alpha4 = Symbol('alpha4')

x1 = Symbol('x1')
x2 = Symbol('x2')
x3 = Symbol('x3')
x4 = Symbol('x4')

y1 = Symbol('y1')
y2 = Symbol('y2')
y3 = Symbol('y3')
y4 = Symbol('y4')

f1 = alpha1 + alpha2*x1 + alpha3*y1 + alpha4*x1*y1 - u1
f2 = alpha1 + alpha2*x2 + alpha3*y2 + alpha4*x2*y2 - u2
f3 = alpha1 + alpha2*x3 + alpha3*y3 + alpha4*x3*y3 - u3
f4 = alpha1 + alpha2*x4 + alpha3*y4 + alpha4*x4*y4 - u4

solution = solve([f1, f2, f3, f4], [alpha1, alpha2, alpha3, alpha4])
print(solution)
