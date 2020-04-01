#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 08:40:35 2019

@author: arthursantiago
"""

import numpy as np #bbt computacao numerica
import sympy as sp #bbt computacao simbólica
sp.init_printing()

#%%

# Cria as variáveis do problema
theta,x_1, x_2, x_3, c = sp.var('theta,X_1, X_2, X_3,C')

#theta como funcao de x1 e x2
theta = sp.Function('theta')(x_1,x_2)


# Vetores canônicos
e_1 = sp.Matrix([1,0,0])
e_2 = sp.Matrix([0,1,0])
e_3 = sp.Matrix([0,0,1])
r = x_1*e_1 + x_2*e_2 + x_3*e_3
t = (1+x_2)*e_1 + (5-x_2)*e_2

def gradf(f,x):
    return sp.Matrix([f.diff(i) for i in x])

T = sp.Matrix([[  x_1+x_2 ,theta  , 0],
                 [theta , x_1- 2*x_2  , 0],
                 [0 , 0 , x_2]
                 ])
eq_1 = T[0,0].diff(x_1) + T[0,1].diff(x_2) + T[0,2].diff(x_3)   
eq_2 = T[1,0].diff(x_1) + T[1,1].diff(x_2) + T[1,2].diff(x_3)   
eq_3 = T[2,0].diff(x_1) + T[2,1].diff(x_2) + T[2,2].diff(x_3)   

print(" ")
print("Primeira condicao de equilibrio")
display(eq_1)

print(" ")
print("Segunda condicao de equilibrio")
display(eq_2)

print(" ")
print("Terceira condicao de equilibrio")
display(eq_2)

#%%

# Criando funcoes auxiliares
lamb,psi = sp.var('lambda,psi')
lamb = sp.Function('lambda')(x_1)
psi = sp.Function('psi') (x_2)

#T12 TEM UMA FUNCAO LAMBDA DE X1 E UMA FUNCAO DE X2 QUE EH -X2
T12A = -x_2 + lamb

# T12 TEM UMA FUNCAO PSI  DE X2 E UMA FUNCAO DE X1 QUE EH 2X1
T12B = 2*x_1 + psi

# SE JUNTARMOS AS DUAS PARTES TEMOS QUE A FUNCAO DE X1 EH 2X1 E A FUNCAO DE X2 EH -X2, INTEGRANDO TEMOS UMA CONSTANTE C
T12 = 2*x_1 - x_2 + c

#SUBSTITUINDO A FUNCAO TTHETA(X1,X2) POR T12 QUE CALCULAMOS
T = T.subs(theta,T12)
T =T.subs(x_1,1)


eq= t-T*e_1

#eq = -x_2 + 5 - c +x_2 -2
#for i in range(3):
#   print (sp.solve(eq[i],c))
T12 = T12.subs(c,sp.solve(eq[1],c)[0])

display(T12)

