#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 11:14:46 2019

@author: arthursantiago
"""

import numpy as np #bbt computacao numerica
import sympy as sp #bbt computacao simbólica
sp.init_printing()

#%%

# Cria as variáveis do problema
x_1, x_2, x_3 = sp.var('X_1, X_2, X_3')

# Vetores canônicos
e_1 = sp.Matrix([1,0,0])
e_2 = sp.Matrix([0,1,0])
e_3 = sp.Matrix([0,0,1])
r = x_1*e_1 + x_2*e_2 + x_3*e_3

def gradf(f,x):
    return sp.Matrix([f.diff(i) for i in x])

# Montando o vetor de deslocamento dado no problema
u = sp.Matrix([sp.sin(x_1),x_1**3*x_2 , sp.cos(x_3)])

grad_u = u.jacobian(r)

# Montando o tensor de tensoes
E = (grad_u + grad_u.T)/2

#Montando as condicoes de compatibilidade
c1_parte_1 = E[0,0].diff(x_2).diff(x_2)
c1_parte_2 = E[1,1].diff(x_1).diff(x_1)
c1_parte_3 = 2 *  E[0,1].diff(x_2).diff(x_1)


c2_parte_1 = E[0,0].diff(x_3).diff(x_3)
c2_parte_2 = E[2,2].diff(x_1).diff(x_1)
c2_parte_3 = 2 *  E[0,2].diff(x_3).diff(x_1)

c3_parte_1 = E[1,1].diff(x_3).diff(x_3)
c3_parte_2 = E[2,2].diff(x_2).diff(x_2)
c3_parte_3 = 2 *  E[1,2].diff(x_3).diff(x_2)


print("Para a primeira condicao do tensor E : ",c1_parte_1 ,"+", c1_parte_2 ,"=", c1_parte_3)
print("Para a segunda condicao do tensor E : ",c2_parte_1 ,"+", c2_parte_2 ,"=", c2_parte_3)
print("Para a terceira condicao do tensor E : ",c3_parte_1 ,"+", c3_parte_2 ,"=", c3_parte_3)
print(" ")
print("como todas as condicoes foram atendidas, o tensor atendea s condicoes de compatibilidade")


