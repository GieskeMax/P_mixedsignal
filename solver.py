import math
from sympy.solvers import solve
from sympy import symbols
from sympy.codegen.cfunctions import Sqrt
from sympy.simplify import sqrtdenest



x = symbols('x')
# solve(2*Sqrt(2*c.i_q/(c.k_0n*x))+c.gamma*(Sqrt(c.phi+Sqrt(2*c.i_q/(c.k_0n*x)))-Sqrt(c.phi))-1)
# print(x)


# print(solve(2 * Sqrt(2 * c.i_q / (c.k_0n * x)) + c.gamma * (sqrtdenest(Sqrt(c.phi+Sqrt(2*c.i_q/(c.k_0n*x)))) - Sqrt(c.phi)) - 1))

print(solve(2*x+5))
