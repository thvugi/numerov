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
    #return abs(s)#((omega**2)*s**2)/2
    R=1.65
    return -60.9*math.exp(-(s/R)**2)
    #if(s>1):
   #     return 0
   # if(s<=1):
   #     return -100

LHS = 0
RHS = 20
L = RHS-LHS

N = 1000

delx = L/N

U = []
VV = []
for j in range(N+1):
    x = LHS+delx*j
    VV.append(V(x))
for j in range(N+1):
    x = LHS+delx*j
    nmb = (1/(h2m))*(-10+V(x))
    U.append(nmb)

BC = 0

def BC1(x):
   return x

psipos = []
psineg=[]
psipos.append(BC)
psipos.append(BC1(LHS+delx))
psineg.append(BC)
psineg.append(BC1((math.exp(0.01*delx))))

for i in range(1,int(N/2)-1):
        temp = (psipos[i-1]*(12-(delx**2)*U[i-1])-2*psipos[i]*((5*delx**2)*U[i]+12))/((delx**2)*U[i+1]-12)
        psipos.append(temp)
    
num = N


for j in range(1,int(N/2)-1):
        temp = (psineg[j-1]*((delx**2)*U[num]-12)+2*psineg[j]*((5*delx**2)*U[num-1]+12))/(12-(delx**2)*U[num-2])
        num=num-1
        psineg.append(temp)

#plt.plot(VV)
plt.plot(psipos,linewidth=1.5,color='red')
#plt.plot(psineg,linewidth=1.5,color='blue')
plt.show()

dict = {'pos value': psipos,'neg value': psineg}
df = pd.DataFrame(dict)
df.to_csv('soln.csv')