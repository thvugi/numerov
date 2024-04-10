import numpy as np
import matplotlib.pyplot as plt
import csv 
import math
import pandas as pd

omega = 1
h2m = 41.47
m=1
E=[]
#number of energy levels
n = 0
for i in range(n):
   E.append((i+.5)*h2m*omega)


def V(s):
    R=1.65
    return -60.9*math.exp(-(s/R)**2)
    
LHS = 0
RHS = 20
L = RHS-LHS

N = 1000
M = 300

delx = L/N

U = []
VV = []
for j in range(N+1):
    x = LHS+delx*j
    VV.append(V(x))
for j in range(N+1):
    x = LHS+delx*j
    nmb = (1/(h2m))*(-2-V(x))
    U.append(nmb)

BC = 0

def BC1(x):
   return x

psipos = np.zeros(N+1)
psineg=np.zeros(N+1)
psipos[0]=0
psipos[1]=delx
psineg[N]=math.exp(-0.05*RHS)
psineg[N-1]=math.exp(-0.05*(RHS+delx))

for i in range(1,M+1):
        psipos[i+1] = (-psipos[i-1]*(12+(delx**2)*U[i-1])+2*psipos[i]*((-5*delx**2)*U[i]+12))/((delx**2)*U[i+1]+12)


for j in range(N-1,M-1,-1):
         psineg[j-1] = (-psineg[j+1]*((delx**2)*U[j+1]+12)+2*psineg[j]*(12-(5*delx**2)*U[j]))/(12+(delx**2)*U[j-1])
        
        

#plt.plot(VV)
plt.plot(psipos,linewidth=1.5,color='red')
#plt.plot(psineg,linewidth=1.5,color='blue')
plt.show()

dict = {'pos value': psipos,'neg value': psineg}
df = pd.DataFrame(dict)
df.to_csv('soln.csv')