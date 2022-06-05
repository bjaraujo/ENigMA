from sympy.solvers import solve
from sympy import Symbol

u1 = Symbol('u1')
u2 = Symbol('u2')
u3 = Symbol('u3')
u4 = Symbol('u4')
u5 = Symbol('u5')
u6 = Symbol('u6')
u7 = Symbol('u7')
u8 = Symbol('u8')

alpha1 = Symbol('alpha1')
alpha2 = Symbol('alpha2')
alpha3 = Symbol('alpha3')
alpha4 = Symbol('alpha4')
alpha5 = Symbol('alpha5')
alpha6 = Symbol('alpha6')
alpha7 = Symbol('alpha7')
alpha8 = Symbol('alpha8')

x1 = Symbol('x1')
x2 = Symbol('x2')
x3 = Symbol('x3')
x4 = Symbol('x4')
x5 = Symbol('x5')
x6 = Symbol('x6')
x7 = Symbol('x7')
x8 = Symbol('x8')

y1 = Symbol('y1')
y2 = Symbol('y2')
y3 = Symbol('y3')
y4 = Symbol('y4')
y5 = Symbol('y5')
y6 = Symbol('y6')
y7 = Symbol('y7')
y8 = Symbol('y8')

z1 = Symbol('z1')
z2 = Symbol('z2')
z3 = Symbol('z3')
z4 = Symbol('z4')
z5 = Symbol('z5')
z6 = Symbol('z6')
z7 = Symbol('z7')
z8 = Symbol('z8')

f1 = alpha1 + alpha2*x1 + alpha3*y1 + alpha4*z1 + alpha5*x1*y1 + alpha6*y1*z1 + alpha7*z1*x1 + alpha8*x1*y1*z1 - u1
f2 = alpha1 + alpha2*x2 + alpha3*y2 + alpha4*z2 + alpha5*x2*y2 + alpha6*y2*z2 + alpha7*z2*x2 + alpha8*x2*y2*z2 - u2
f3 = alpha1 + alpha2*x3 + alpha3*y3 + alpha4*z3 + alpha5*x3*y3 + alpha6*y3*z3 + alpha7*z3*x3 + alpha8*x3*y3*z3 - u3
f4 = alpha1 + alpha2*x4 + alpha3*y4 + alpha4*z4 + alpha5*x4*y4 + alpha6*y4*z4 + alpha7*z4*x4 + alpha8*x4*y4*z4 - u4
f5 = alpha1 + alpha2*x5 + alpha3*y5 + alpha4*z5 + alpha5*x5*y5 + alpha6*y5*z5 + alpha7*z5*x5 + alpha8*x5*y5*z5 - u5
f6 = alpha1 + alpha2*x6 + alpha3*y6 + alpha4*z6 + alpha5*x6*y6 + alpha6*y6*z6 + alpha7*z6*x6 + alpha8*x6*y6*z6 - u6
f7 = alpha1 + alpha2*x7 + alpha3*y7 + alpha4*z7 + alpha5*x7*y7 + alpha6*y7*z7 + alpha7*z7*x7 + alpha8*x7*y7*z7 - u7
f8 = alpha1 + alpha2*x8 + alpha3*y8 + alpha4*z8 + alpha5*x8*y8 + alpha6*y8*z8 + alpha7*z8*x8 + alpha8*x8*y8*z8 - u8

solution = solve([f1, f2, f3, f4, f5, f6, f7, f8], [alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8])
print(solution)
