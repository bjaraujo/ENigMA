from sympy.solvers import solve
from sympy import Symbol

u1 = Symbol('u1')
u2 = Symbol('u2')
u3 = Symbol('u3')

alpha1 = Symbol('alpha1')
alpha2 = Symbol('alpha2')
alpha3 = Symbol('alpha3')

x1 = Symbol('x1')
x2 = Symbol('x2')
x3 = Symbol('x3')

f1 = alpha1 + alpha2*x1 + alpha3*x1**2 - u1
f2 = alpha1 + alpha2*x2 + alpha3*x2**2 - u2
f3 = alpha1 + alpha2*x3 + alpha3*x3**2 - u3

solution = solve([f1, f2, f3], [alpha1, alpha2, alpha3])
print(solution)
