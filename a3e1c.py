#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 11:38:48 2019

@author: arthursantiago
"""

import numpy as np #bbt computacao numerica
import sympy as sp #bbt computacao simbólica
sp.init_printing()

#%%

# Cria as variáveis do problema
theta,alpha,v,x_1, x_2, x_3 = sp.var('theta,alpha,v,X_1, X_2, X_3')

#theta como funcao de x1 e x2
theta = sp.Function('theta')(x_2,x_3)


# Vetores canônicos
e_1 = sp.Matrix([1,0,0])
e_2 = sp.Matrix([0,1,0])
e_3 = sp.Matrix([0,0,1])
r = x_1*e_1 + x_2*e_2 + x_3*e_3

def gradf(f,x):
    return sp.Matrix([f.diff(i) for i in x])

# Montando o Tensor de tensoes dado no problema
E_1 = sp.Matrix([[  (1/alpha)*theta ,0  , 0],
                 [0 ,(-v/alpha)*theta  , 0],
                 [0 , 0 , (-v/alpha)*theta]
                 ])
    

# Montando as condicoes de compatibilidade
c1_parte_1 = E_1[0,0].diff(x_2).diff(x_2)
c1_parte_2 = E_1[1,1].diff(x_1).diff(x_1)
c1_parte_3 = 2 *  E_1[0,1].diff(x_2).diff(x_1)
display(c1_parte_3)
print(" ")
print("como a Derivada segunda em relaxao a X_2 e X_3 eh igual a zero, a funcao eh linear ")




c2_parte_1 = E_1[0,0].diff(x_3).diff(x_3)
c2_parte_2 = E_1[2,2].diff(x_1).diff(x_1)
c2_parte_3 = 2 *  E_1[0,2].diff(x_3).diff(x_1)
display(c2_parte_3)
print(" ")
print("como a Derivada segunda em relaxao a X_2 e X_3 eh igual a zero, a funcao eh linear ")


c3_parte_1 = E_1[1,1].diff(x_3).diff(x_3)
c3_parte_2 = E_1[2,2].diff(x_2).diff(x_2)
c3_parte_3 = 2 *  E_1[1,2].diff(x_3).diff(x_2)
display(c3_parte_3)
print(" ")
print("como a Derivada segunda em relaxao a X_2 e X_3 eh igual a zero, a funcao eh linear ")



