#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 10:33:26 2019

@author: arthursantiago
"""

import numpy as np #bbt para uso de computação numérica
import sympy as sp #bbt para uso de computação simbólica
sp.init_printing()


#criando as funções que retornam cálculo de gradiante e laplaciano de uma função

def grad(f,x):
    return sp.Matrix([f.diff(i) for i in x])

def laplacian(f,x):
    return sum([f.diff(i,2) for i in x]).simplify()


#crianção das variáveis presentes no exercício

M_t, a, psi, G, theta, z, x, y = sp.var('M_t, a, psi, G, theta, Z, X, Y')
C=sp.var(['C_'+str(i) for i in range(10)])

#criando as bordas da seção transversal (triângulo)

f=sp.Matrix([
        y + sp.sqrt(3)*x - 2*a,
        y - sp.sqrt(3)*x - 2*a,
        y -a,
        ])
    

p2=[]

for i in f:
    p1=sp.plot_implicit(
            sp.Eq(i.subs({a:1}),0),
            (x, -4, 4), (y, -4,4)
            )
    p2.append(p1)

[p2[0].extend(i) for i in p2]
p2[0].show()


#criação de psi dada no exercício

psi, alpha, mu, g = sp.var('psi, alpha, mu, G')
theta = sp.Function('theta')(x,y,z)

psi = -g*theta.diff(z)*((x**2+y**2)/2-(x**3-3*x*y**2)/(2*a) -2*a**2/27)


T = sp.zeros(3)
T[0,1]=T[1,0]=+psi.diff(y,).simplify()
T[0,2]=T[2,0]=-psi.diff(x,).simplify()
display(T)