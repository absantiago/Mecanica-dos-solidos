import numpy as np #bbt para uso de computação numérica
import sympy as sp #bbt para uso de computação simbólica
import matplotlib.pyplot as plt
import math
sp.init_printing()

#%%

psi, f_1, f_2, a, b, k, l, x_1, x_2, x_3, pi = sp.var('psi, f_1, f_2, a, b, k, l, x_1, x_2, x_3, pi')
psi = sp.Function('psi')(x_2,x_3)



#%%
   
f_1=sp.Matrix([
     x_2 -a/2,
     x_2 +a/2,
    ])

f_2=sp.Matrix([
     x_3 - b/2,
     x_3 + b/2,
    ])
 
#aux=list( sp.utilities.iterables.subsets(f,3) ) # combina 2 a 2
#display(psi)    

p2=[]
for i in f_1:
   # display(psi.subs({'x_2':sp.solve(i,x_2)[0]}).simplify())
    p1=sp.plot_implicit(sp.Eq(i.subs({a:2}), 0),(x_2, -4, 4), (x_3, -4, 4))
    p2.append(p1)

for i in f_2:
  #  display(psi.subs({'x_3':sp.solve(i,x_3)[0]}).simplify())
    p1=sp.plot_implicit(sp.Eq(i.subs({b:2}), 0),(x_2, -4, 4), (x_3, -4, 4))
    p2.append(p1)

[p2[0].extend(i) for i in p2]
p2[0].show()


#%%
psi_t = ((2**5)*(a**2)*(b**2))/(sp.pi**4) * ((-1)**((k+l)/2 -1))/(k*l*((k**2*b**2)+(l**2*a**2))) * sp.cos((k*x_2*sp.pi)/a) * sp.cos((l*x_3*sp.pi)/b)

psi = ((-1)**((k+l)/2 -1))/(k*l*((k**2*b**2)+(l**2*a**2))) * sp.cos((k*x_2*sp.pi)/a) * sp.cos((l*x_3*sp.pi)/b)

display(psi) 

#%%
   
psi = psi.subs(a,2)
psi = psi.subs(b,2)

#%%

psi1, psi2, psi3, psi4, psi5 = sp.var('psi1, psi2, psi3, psi4, psi5')

psi1 = psi.subs(k,1)
psi1 = psi1.subs(l,1)

psi2 = psi.subs(k,3)
psi2 = psi1 + psi2.subs(l,3)

psi3 = psi.subs(k,5)
psi3 = psi2 + psi3.subs(l,5)

psi4 = psi.subs(k,7)
psi4 = psi3 + psi4.subs(l,7)

psi5 = psi.subs(k,9)
psi5 = psi4 + psi5.subs(l,9)

psi_c = ((2**5)*(a**2)*(b**2))/(sp.pi**4) * psi5

psi_c = psi_c.subs(a,2)
psi_c = psi_c.subs(b,2)

display(psi_c)
psi_diff = psi_c.diff(x_3)

fi = sp.integrate(psi_diff,x_2)
display(fi)

fi1 = []
for i in range (100):
    aux = fi.subs(x_2,i*0.001)
    aux = aux.subs(x_3,i*0.001)
    fi1.append(aux)

plt.plot(fi1)
plt.axis([-100, 100, -0.02, 0.02])
plt.show()
