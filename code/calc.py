from sympy import * 

g, l, s, h = symbols('g l s h')

system = [ s - 1/(g+l), 
            h - g/(g+l)
            ]

nonlinsolve(system, [g,l])