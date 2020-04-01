# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 08:45:34 2019

@author: Tales
"""

import numpy as np #bbt computacao numerica
import sympy as sp #bbt computacao simbólica
sp.init_printing()

# Cria as variáveis do problema
x_1, x_2, x_3 = sp.var('X_1, X_2, X_3')
u, k, r = sp.var('u, k, r')
lamb, mu = sp.var('lambda, mu')

# Vetores canônicos
e_1 = sp.Matrix([1,0,0])
e_2 = sp.Matrix([0,1,0])
e_3 = sp.Matrix([0,0,1])
r = x_1*e_1 + x_2*e_2 + x_3*e_3

u = sp.Matrix([X_2*X_3, X_3*X_1, X_1**2 - X_2**2])
u = k*u
 
# Tensores de deformação (E) e de tensão(T)
grad_u = u.jacobian(r)
E = (grad_u + grad_u.T)/2
T = lamb * E.trace() * sp.eye(3) + 2 * mu * E

eql_1 = T[0,0].diff(x_1) +T[0,1].diff(x_2) + T[0,2].diff(x_3)
eql_2 = T[1,0].diff(x_1) +T[1,1].diff(x_2) + T[1,2].diff(x_3)
eql_3 = T[2,0].diff(x_1) +T[2,1].diff(x_2) + T[2,2].diff(x_3)

eql = sp.Matrix([eql_1, eql_2, eql_3])


print('Na ausência de forças de corpo, B vale zero, portanto, o gradiente de T deve ser nulo.')
print('')
print('Primeiro vetor de deslocamento: ')
display(u)
print('O gradiente do tensor T é: ')
display(eql)
print('Portanto, o campo de Tensão está em equilíbrio.')
print('')
#%%
#---------------------------------------------------------------------------------------

# Cria as variáveis do problema
y_1, y_2, y_3 = sp.var('Y_1, Y_2, Y_3')
u_2, k, r_2 = sp.var('u_2, k, r_2')
lamb, mu = sp.var('lambda, mu')

# Vetores canônicos
e_1 = sp.Matrix([1,0,0])
e_2 = sp.Matrix([0,1,0])
e_3 = sp.Matrix([0,0,1])
r_2 = y_1*e_1 + y_2*e_2 + y_3*e_3

u_2 = sp.Matrix([y_2*y_3, y_3*y_1, y_1*y_2])
u_2 = k*u_2
 
# Tensores de deformação (E) e de tensão(T)
grad_u_2 = u_2.jacobian(r_2)
E_2 = (grad_u_2 + grad_u_2.T)/2
T_2 = lamb * E_2.trace() * sp.eye(3) + 2 * mu * E_2

eql_12 = T_2[0,0].diff(y_1) +T_2[0,1].diff(y_2) + T_2[0,2].diff(y_3)
eql_22 = T_2[1,0].diff(y_1) +T_2[1,1].diff(y_2) + T_2[1,2].diff(y_3)
eql_32 = T_2[2,0].diff(y_1) +T_2[2,1].diff(y_2) + T_2[2,2].diff(y_3)

eql_2 = sp.Matrix([eql_12, eql_22, eql_32])


print('Segundo vetor de deslocamento: ')
display(u_2)
print('O Divergente do tensor T é: ')
display(eql_2)
print('Portanto, o campo de Tensão está em equilíbrio.')
print('')
#%%
#--------------------------------------------------------------------------------

# Cria as variáveis do problema
z_1, z_2, z_3 = sp.var('Z_1, Z_2, Z_3')
u_3, k, r_3 = sp.var('u_3, k, r_3')
lamb, mu = sp.var('lambda, mu')

# Vetores canônicos
e_1 = sp.Matrix([1,0,0])
e_2 = sp.Matrix([0,1,0])
e_3 = sp.Matrix([0,0,1])
r_3 = z_1*e_1 + z_2*e_2 + z_3*e_3

u_3 = sp.Matrix([z_2*z_3, z_1*z_3, z_1*z_2 + z_3**2])
u_3 = k*u_3
 
# Tensores de deformação (E) e de tensão(T)
grad_u_3 = u_3.jacobian(r_3)
E_3 = (grad_u_3 + grad_u_3.T)/2
T_3 = lamb * E_3.trace() * sp.eye(3) + 2 * mu * E_3

eql_13 = T_3[0,0].diff(z_1) +T_3[0,1].diff(z_2) + T_3[0,2].diff(z_3)
eql_23 = T_3[1,0].diff(z_1) +T_3[1,1].diff(z_2) + T_3[1,2].diff(z_3)
eql_33 = T_3[2,0].diff(z_1) +T_3[2,1].diff(z_2) + T_3[2,2].diff(z_3)

eql_3 = sp.Matrix([eql_13, eql_23, eql_33])

print('Terceiro vetor de deslocamento: ')
display(u_3)
print('O gradiente do tensor T é: ')
display(eql_3)
print('Portanto, o campo de Tensão não está em equilíbrio, pois o grad(T) != 0.')
print('')
