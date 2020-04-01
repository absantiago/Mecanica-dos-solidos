#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 10:16:00 2019

@author: arthursantiago
"""

import numpy as np #bbt computacao numerica
import sympy as sp #bbt computacao simbólica
sp.init_printing()

#%%

# Cria as variáveis do problema
k,x_1, x_2, x_3 = sp.var('K,X_1, X_2, X_3')

# Vetores canônicos
e_1 = sp.Matrix([1,0,0])
e_2 = sp.Matrix([0,1,0])
e_3 = sp.Matrix([0,0,1])
r = x_1*e_1 + x_2*e_2 + x_3*e_3

def gradf(f,x):
    return sp.Matrix([f.diff(i) for i in x])

# Tensores de deformação (E) e de tensão(T)
E_1 = sp.Matrix([[x_1**2,x_3**2 + x_2**2 ,x_1*x_3],
                 [x_3**2 + x_2**2 , 0 ,x_1],
                 [x_1*x_3 , x_1 , x_2**2]
                 ])
    
E_1 = E_1*k
    

# montando a primeira condicao de compatibilidade ---------
c1_parte_1 = E_1[0,0].diff(x_2).diff(x_2)
c1_parte_2 = E_1[1,1].diff(x_1).diff(x_1)
c1_parte_3 = 2 *  E_1[0,1].diff(x_2).diff(x_1)
#display(c1_parte_1,c1_parte_2,c1_parte_3)


# montando a segunda condicao de compatibilidade ---------
c2_parte_1 = E_1[0,0].diff(x_3).diff(x_3)
c2_parte_2 = E_1[2,2].diff(x_1).diff(x_1)
c2_parte_3 = 2 *  E_1[0,2].diff(x_3).diff(x_1)


# montando a terceira condicao de compatibilidade ---------
c3_parte_1 = E_1[1,1].diff(x_3).diff(x_3)
c3_parte_2 = E_1[2,2].diff(x_2).diff(x_2)
c3_parte_3 = 2 *  E_1[1,2].diff(x_3).diff(x_2)


print("Para a primeira condicao do tensor E_1 : ",c1_parte_1 ,"+", c1_parte_2 ,"=", c1_parte_3)
print("Para a segunda condicao do tensor E_1 : ",c2_parte_1 ,"+", c2_parte_2 ,"!=", c2_parte_3)
print("Para a terceira condicao do tensor E_1 : ",c3_parte_1 ,"+", c3_parte_2 ,"!=", c3_parte_3)
print(" ")
print("como a terceira condicao de compatibilidade nao foi satisfeita, o tensor E_1 não atende as condicoes ")




# Montando o segundo Tensor
E_2 = sp.Matrix([[x_1+x_2 , x_1 , x_2],
                 [x_1 ,x_3 + x_2 , x_3],
                 [x_2 , x_3 , x_1 + x_3]
                 ])
    
E_2 = E_2*k


# Montando as condicoes de compatibilidade do segundo tensor
d1_parte_1 = E_2[0,0].diff(x_2).diff(x_2)
d1_parte_2 = E_2[1,1].diff(x_1).diff(x_1)
d1_parte_3 = 2 *  E_2[0,1].diff(x_2).diff(x_1)

d2_parte_1 = E_2[0,0].diff(x_3).diff(x_3)
d2_parte_2 = E_2[2,2].diff(x_1).diff(x_1)
d2_parte_3 = 2 *  E_2[0,2].diff(x_3).diff(x_1)

d3_parte_1 = E_2[1,1].diff(x_3).diff(x_3)
d3_parte_2 = E_2[2,2].diff(x_2).diff(x_2)
d3_parte_3 = 2 *  E_2[1,2].diff(x_3).diff(x_2)

print(" ")
print("Para a primeira condicao do tensor E_2 : ",d1_parte_1 ,"+", d1_parte_2 ,"=", d1_parte_3)
print("Para a segunda condicao do tensor E_2 : ",d2_parte_1 ,"+", d2_parte_2 ,"=", d2_parte_3)
print("Para a terceira condicao do tensor E_2 : ",d3_parte_1 ,"+", d3_parte_2 ,"=", d3_parte_3)
print(" ")
print("como todas as condicoes foram atendidas, o tensor atendea s condicoes de compatibilidade")



