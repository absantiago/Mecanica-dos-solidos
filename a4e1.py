import numpy as np #bbt para uso de computação numérica
import sympy as sp #bbt para uso de computação simbólica
sp.init_printing()

#%%
#criando as funções que retornam cálculo de gradiante e laplaciano de uma função

def grad(f,x):
    return sp.Matrix([f.diff(i) for i in x])

def laplacian(f,x):
    return sum([f.diff(i,2) for i in x]).simplify()

#%%
#crianção das variáveis presentes no exercício

M_t, a, theta, z, x, y, M = sp.var('M_t, a, theta, z, x, y, M')
C=sp.var(['C_'+str(i) for i in range(10)])

#%%
#criando as bordas da seção transversal (triângulo)

f=sp.Matrix([
        x + sp.sqrt(3)*y - 2*a/3,
        x - sp.sqrt(3)*y - 2*a/3,
        x + a/3,
        ])
   

p2=[]

for i in f:
    p1=sp.plot_implicit(
            sp.Eq(i.subs({a:3}),0),
            (x, -4, 4), (y, -4,4)
            )
    p2.append(p1)

[p2[0].extend(i) for i in p2]
p2[0].show()

#%%
#criação de psi dada no exercício

psi, alpha, mu, G, phi = sp.var('psi, alpha, mu, G, phi')
theta = sp.Function('theta')(z)

psi = -G*theta.diff(z)*((x**2+y**2)/2-(x**3-3*x*y**2)/(2*a) -2*a**2/27)


T = sp.zeros(3)
T[0,1]=T[1,0]=+psi.diff(y,).simplify()
T[0,2]=T[2,0]=-psi.diff(x,).simplify()

print('Os valores de T12=T21 e T13=T31 procurados são mostrados no tensor T abaixo:')
display(T)


M = sp.integrate(-2*psi,
                 (
                     y,
                     sp.solve(f[1],y)[0],
                     sp.solve(f[0],y)[0],
             ),
             (
                     x,
                     -a/3,
                     2*a/3
            ),
        )


eq = M-M_t
resp = sp.solve(eq, theta.diff(z))

print('A derivada do ângulo de torção em relação a z é igual a:')
display(resp)

phi = sp.Function('phi')(x,y,z)

#%%
#determinando a função de empenamento

parte1 = psi.diff(z)
parte2 = -G*z*theta.diff(z)
parte3 = G*phi.diff(y)


empenamento = sp.integrate(z*resp[0],
                           
                           (
                            y,
                            sp.solve(f[1],y)[0],
                            sp.solve(f[0],y)[0],
                           )      
       
        )

print('Função de Empenamento é igual a:')
display(empenamento.simplify())


#%%Parte referente à letra a) -> não chegamos a uma conclusão ainda
respa1 = psi.subs(y, sp.solve(f[1],
                              y)[0])
respa2 = psi.subs(y, sp.solve(f[0],y)[0])
respa3 = psi.subs(x, -a/3)


display()