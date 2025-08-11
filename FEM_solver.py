import numpy as np
import pandas as pd
import math
from sympy import diff, integrate, simplify
import sympy as sp
import matplotlib.pyplot as plt
from scipy.integrate import simpson

#Number of elements
N = 10
M_vector = [3,4]
def u_analytical(x):
    term1 = np.exp(np.pi)*np.sin(x/2)
    term2 = (2*(3-np.exp(2*np.pi))/np.exp(np.pi))*np.cos(x/2)
    return 2*np.exp(-x)*(term1+term2)+np.sin(x)

def L_e(U,u,x,N,nodes): 
    aux = 0
    for i in range(N): 
        y_values =U[2*i:2*i+nodes]-u[2*i:2*i+nodes]
        x_values =x[2*i:2*i+nodes]
        aux+= (simpson(y_values, x_values))**(2)
        
    return (aux)**(0.5)

plt.plot(1)
for M in M_vector: 
    #Number of nodes
    n=((M+1)*N)+1

    xi = sp.Symbol('xi')

    # Define the mapping from x in [0,1] to xi in [-1,1]
    x = (xi + 1)/2

    # Nodes in the physical domain
    x_var = [float(val) for val in np.linspace(0, 1, M+2)]
    x_var = [sp.sympify(val) for val in x_var]
    phi = [sp.S.Zero for _ in range(M+2)]
    for j in range(2+M): 
        aux=1
        for m in range(M+2): 
            if m !=j: 
                aux= aux*(x-x_var[m])/(x_var[j]-x_var[m])
        phi[j]=sp.simplify(aux)

    #length or hight of each element 
    h=(2*math.pi)/(N) 

    #Creating local stiffness matrix 
    K_local = np.zeros((M+2,M+2),dtype=float)

    #Creating local length 
    h_local = (2*math.pi)/(n-1)

    #Creating the global force vector 
    f  = np.zeros(n)

    #Creating the global stiffness matrix 
    K = np.zeros((n,n),dtype=float)

    for e in range(N):
        #Defining local stiffness matrix at each element 
        x_min = -math.pi+h*e
        x_max = -math.pi+h*(e+1)
        je = (x_max-x_min)/2
        x = je*xi+(x_min+x_max)/2
        for i in range(0,(M+2)): 
            for j in range(0,(M+2)): 
                aux = (-4*(1/je**2)*diff(phi[i],xi)*diff(phi[j],xi)+8*(1/je)*phi[i]*diff(phi[j],xi)+5*phi[i]*phi[j])*je
                K_local[i,j] = float(integrate(simplify(aux), (xi, -1,1)))

        #Assembling local force vector
        f_local = np.zeros(((M+2)))
        for i in range(len(f_local)): 
            f_local[i]=float(integrate((sp.sin(x)+8*sp.cos(x))*phi[i]*je,(xi,-1,1)).evalf())

        #Assembling the global stiffness matrix
        global_dofs = np.zeros(((M+2)))
        for i in range(M+2): 
            global_dofs[i] =int((M+1)*e+i)

        for i in range(M+2):
            for j in range(len(K_local[0,:])):
                K[int(global_dofs[i]), int(global_dofs[j])] += K_local[i, j]
        f[int(global_dofs[0]):int(global_dofs[-1]+1)]+=f_local

    
    ##For boundary conditions: 
    f[0]+=-20
    u_pi=2
    for j in range(0,n-1):
        f[j]-= u_pi*K[j,-1]

    u=np.linalg.solve(K[0:-1,0:-1],f[0:-1])
    x_analytical =np.linspace(-np.pi,np.pi,len(u))
    u_an=u_analytical(x_analytical)
    nodes = M+2
    plt.plot(x_analytical,u,label=f'M={M}')

plt.title(f'Quadratic Basis Function ({N} Elements)')
plt.plot(x_analytical,u_an,label='Analytical')
plt.xlabel('Position, x')
plt.ylabel('u')
plt.xlim([-np.pi,np.pi])
plt.legend()
plt.show()
