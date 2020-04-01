# -*- coding: utf-8 -*-
"""
Spyder Editor

Este é um arquivo de script temporário.
"""

import numpy as np #bbt computacao numerica
import sympy as sp #bbt computacao simbólica
sp.init_printing()


# Cria as variáveis do problema
phi, x_1, x_2, x_3 = sp.var('phi, X_1, X_2, X_3')
alpha_1, alpha_2, alpha_3 = sp.var('alpha_1, alpha_2, alpha_3')

phi = alpha_1*(x_1**2) + alpha_2*x_1*x_2 + alpha_3*(x_2**2)

# Diferenciando 4 vezes phi com relação a X1
d_11 = phi.diff(x_1)
d_12 = d_11.diff(x_1)
d_13 = d_12.diff(x_1)
d_14 = d_13.diff(x_1)

# Diferenciando 4 vezes phi com relação a X2
d_21 = phi.diff(x_2)
d_22 = d_21.diff(x_2)
d_23 = d_22.diff(x_2)
d_24 = d_23.diff(x_2)

# Diferenciando 2 vezes phi com relação a X2 e 2 vezes com relaçao a X1
d_31 = phi.diff(x_2)
d_32 = d_31.diff(x_2)
d_33 = d_32.diff(x_1)
d_34 = d_33.diff(x_1)

print('Primeiro caso:')
print('Equação 1:')
display(phi)
print('Termos do gradiente à quarta:')
display(d_14, d_24, d_34)
print('A função de Airy satisfaz a função bi-harmônica, pois o seu gradiente à quarta é nulo')

#Determinando T11,T12 e T22

T_11, T12, T22 = sp.var('T_11, T_12, T_22')

T_11 = d_22
T_22 = d_12
T_12 = -(phi.diff(x_2).diff(x_1))

print('As Tensões T_11, T_12 e T_22 são respectivamente:')
display(T_11, T_12, T_22)

#%%
#-----------------------------------------------------------------------------

# Segunda fundação phi dada

# Cria as variáveis do problema
theta, y_1, y_2, y_3 = sp.var('theta, Y_1, Y_2, Y_3')
alpha_1 = sp.var('alpha_1')


theta = alpha_1*(y_1**2)*y_2

# Diferenciando 4 vezes phi com relação a X1
di_11 = theta.diff(y_1)
di_12 = di_11.diff(y_1)
di_13 = di_12.diff(y_1)
di_14 = di_13.diff(y_1)

# Diferenciando 4 vezes phi com relação a X2
di_21 = theta.diff(y_2)
di_22 = di_21.diff(y_2)
di_23 = di_22.diff(y_2)
di_24 = di_23.diff(y_2)

# Diferenciando 2 vezes phi com relação a X2 e 2 vezes com relaçao a X1
di_31 = theta.diff(y_2)
di_32 = di_31.diff(y_2)
di_33 = di_32.diff(y_1)
di_34 = di_33.diff(y_1)

print('Segundo caso:')
print('Equação 2: ')
display(theta)
print('Termos do gradiente à quarta:')
display(di_14, di_24, di_34)
print('')
print('A função de Airy satisfaz a função bi-harmônica, pois o seu gradiente à quarta é nulo')

#Determinando T11,T12 e T22

Te_11, Te_12, Te_22 = sp.var('Te_11, Te_12, Te_22')


Te_11 = di_22
Te_22 = di_12
Te_12 = -(theta.diff(y_2).diff(y_1))

print('As Tensões Te_11, Te_12 e Te_22 são respectivamente:')
display(Te_11, Te_12, Te_22)
#%%
#--------------------------------------------------------------------------

# Cria as variáveis do problema
eta, z_1, z_2, z_3 = sp.var('eta, Z_1, Z_2, Z_3')
alpha_1 = sp.var('alpha_1')


eta = alpha_1*(z_1**4 - z_2**4)

# Diferenciando 4 vezes phi com relação a X1
dif_11 = eta.diff(z_1)
dif_12 = dif_11.diff(z_1)
dif_13 = dif_12.diff(z_1)
dif_14 = dif_13.diff(z_1)

# Diferenciando 4 vezes phi com relação a X2
dif_21 = eta.diff(z_2)
dif_22 = dif_21.diff(z_2)
dif_23 = dif_22.diff(z_2)
dif_24 = dif_23.diff(z_2)

# Diferenciando 2 vezes phi com relação a X2 e 2 vezes com relaçao a X1
dif_31 = eta.diff(z_2)
dif_32 = dif_31.diff(z_2)
dif_33 = dif_32.diff(z_1)
dif_34 = dif_33.diff(z_1)

print('Terceiro caso')
print('Equação 3:')
display(eta)
print('Termos do gradiente à quarta:')
display(dif_14, dif_24, dif_34)
print('')
print('A função de Airy satisfaz a função bi-harmônica, pois o seu gradiente à quarta é nulo')

#Determinando T11,T12 e T22

Ten_11, Ten_12, Ten_22 = sp.var('Ten_11, Ten_12, Ten_22')

Ten_11 = dif_22
Ten_22 = dif_12
Ten_12 = -(eta.diff(z_2).diff(z_1))

print('As Tensões Ten_11, Ten_12 e Ten_22 são respectivamente:')
display(Ten_11, Ten_12, Ten_22)


#%%
#--------------------------------------------------------------------------
#Quarta função dada

# Cria as variáveis do problema
psi, w_1, w_2, w_3 = sp.var('psi, W_1, W_2, W_3')
alpha_1 = sp.var('alpha_1')

psi = alpha_1*(w_1*w_2**2) + w_1*w_2**3

# Diferenciando 4 vezes phi com relação a X1
dife_11 = psi.diff(w_1)
dife_12 = dife_11.diff(w_1)
dife_13 = dife_12.diff(w_1)
dife_14 = dife_13.diff(w_1)

# Diferenciando 4 vezes phi com relação a X2
dife_21 = psi.diff(w_2)
dife_22 = dife_21.diff(w_2)
dife_23 = dife_22.diff(w_2)
dife_24 = dife_23.diff(w_2)

# Diferenciando 2 vezes phi com relação a X2 e 2 vezes com relaçao a X1
dife_31 = psi.diff(w_2)
dife_32 = dife_31.diff(w_2)
dife_33 = dife_32.diff(w_1)
dife_34 = dife_33.diff(w_1)

print('Quarto caso:')
print('Equação 4:')
display(psi)
print('Termos do gradiente à quarta:')
display(dife_14, dife_24, dife_34)
print('')
print('A função de Airy satisfaz a função bi-harmônica, pois o seu gradiente à quarta é nulo')

#Determinando T11,T12 e T22

Tens_11, Tens_12, Tens_22 = sp.var('Tens_11, Tens_12, Tens_22')


Tens_11 = dife_22
Tens_22 = dife_12
Tens_12 = -(psi.diff(w_2).diff(w_1))

print('As Tensões Tens_11, Tens_12 e Tens_22 são respectivamente:')
display(Tens_11, Tens_12, Tens_22)