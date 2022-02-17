# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 10:17:41 2022

@author: TH0LXL1
"""

from sympy import symbols, integrate,simplify
h,y,L = symbols('h y L')
expr = integrate(y**2/(y**2+h**2)**2,y)
expr2 = expr.subs(y,L/2)-expr.subs(y,-L/2)
expr2 = simplify(expr2)
print(expr2)