#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 10:03:07 2019

@author: arthursantiago
"""

import numpy as np #bbt computacao numerica
import sympy as sp #bbt computacao simbólica
import math
sp.init_printing()

#%%

# Cria as variáveis do problema
theta, phi, x_1, x_2, x_3 = sp.var('theta, phi, X_1, X_2, X_3')
lamb, mu = sp.var('lambda, mu')

a = sp.var('a', real=True, positive=True)
b = sp.var('b', real=True, positive=True)

# Theta em função de x1
theta = sp.Function('theta')(x_1)

# Vetores canônicos
e_1 = sp.Matrix([1,0,0])
e_2 = sp.Matrix([0,1,0])
e_3 = sp.Matrix([0,0,1])
r = x_1*e_1 + x_2*e_2 + x_3*e_3

T = sp.Matrix([[  0 ,-2*x_3  , x_2],
                 [-2*x_3 ,0 , 0],
                 [x_2 , 0 , 0]
                 ])
    
#%%

#estudo de caso: Elipse com eixos a E b

def gradf(f,x):
    return sp.Matrix([f.diff(i) for i in x])

f = x_2**2/1**2 + x_3**2/0.5 -1
grad_f = gradf(f,r)
n = grad_f / sp.sqrt(grad_f.dot(grad_f))



var = T*n
display(var)


t = T*(-e_1)

#eq_equilibrio = t[0] + t[1] + t[2]

fr = []
for i in range(3):
    fr.append(  sp.integrate (t[i],  (x_2,
                     -sp.sqrt(1**2 - 2*(x_3**2)),
                     sp.sqrt(1**2 - 2*(x_3**2))
                     ),
                    (x_3,-1,1)
    )
    )
    
    
print("A FORCA RESULTANTE EM CADA UM DOS SENTIDOS EH: ","e1: ", fr[0] , "e2: ", fr[1] ,"e3: ", fr[2])
print("A FORCA RESULTANTE EH NULA POIS A ORIGEM ESTA SOBRE A LINHA NEUTRA: ")


da = sp.var('da')
M = []
r= r.subs(x_1,0)
vec = r.cross(t)


# REPRESENTACAO DO MOMENTO DE INERCIA
vec[0] = -sp.Integral(vec[0],a)
display(vec[0])

# SABENDO QUE O MOMENTO DE INERCIA EM UMA ELIPSE EH DADO POR 1/4 * PI*A*Bˆ3 NO EIXO X
# E O MOMENTO NO EIXO Y EH DADO POR 1/4 * PI*B*Aˆ3
# TEMOS QUE :

# CALCULANDO O MOMENTO DE INERCIA EM REALACAO A X E Y
parte_1 = 1/4 * (math.pi*1*math.sqrt(0.5))
parte_2 = 1/4 * (math.pi*1*math.sqrt(0.5)**3)

display (-parte_1 - 2*parte_2)


