from sympy.solvers import solve
from sympy.plotting import plot
from sympy import Symbol
 
u1 = Symbol('u1')
u2 = Symbol('u2')

alpha1 = Symbol('alpha1')
alpha2 = Symbol('alpha2')

x1 = Symbol('x1')
x2 = Symbol('x2')

f1 = alpha1 + alpha2*x1 - u1
f2 = alpha1 + alpha2*x2 - u2
 
solution = solve([f1, f2], [alpha1, alpha2])
#print(solution)

alpha1 = solution[alpha1]
alpha2 = solution[alpha2]

x = Symbol('x')

u = alpha1 + alpha2 * x
#print(u)

n1 = u.subs({u1:1, u2:0})
n2 = u.subs({u1:0, u2:1})

xi = Symbol('xi')

n1 = n1.subs({x:xi, x1:0, x2:1})
n2 = n2.subs({x:xi, x1:0, x2:1})

print(n1)
print(n2)

#plot(n1, n2, (xi, 0, 1))
#plot(n1)
#plot(n2)
