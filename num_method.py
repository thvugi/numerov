import numpy as np
import matplotlib.pyplot as plt
import csv 
import math
import pandas as pd


#parameters 
omega = 1 
h2m = 41.47 #planck constant with mass

E=[] #energy
for i in range(2000):
    E.append((-4 +i*.002)) #initial guess
delE = []
Delta = []


def soln(LHS, delx, E, l, M, N, psipos, psineg, Delta, delE, U):
    for m in range(N+1):
        x = LHS+delx*m
        nmb = (1/(h2m))*(E[l]-V(x))
        U[m]=nmb

    for i in range(1,M+1):
        psipos[i+1] = (-psipos[i-1]*(12+(delx**2)*U[i-1])+2*psipos[i]*((-5*delx**2)*U[i]+12))/((delx**2)*U[i+1]+12)


    for j in range(N-1,M-1,-1):
        psineg[j-1] = (-psineg[j+1]*((delx**2)*U[j+1]+12)+2*psineg[j]*(12-(5*delx**2)*U[j]))/(12+(delx**2)*U[j-1])

    negprime = (psineg[M]-psineg[M-1])/delx
    posprime = (psipos[M+1]-psipos[M])/delx

    Delta.append(np.abs(-(negprime/psineg[M])+(posprime/psipos[M]))) #this computes the delta
    delE.append(abs(E[l]-E[l-1]))
    plot = np.zeros(N+1) #this gives the plot connecting the right and left graphs 
    for i in range(N+1):
        if i >M+1:
            plot[i] = np.flip(psineg[abs(N-i)])
        else:
            plot[i]= psipos[i]

    #plt.plot(VV)
    #plt.plot(psipos,linewidth=1.5,color='red')
    #plt.plot(np.flip(psineg),linewidth=1.5,color='blue')
    plt.plot(plot)
    #plt.plot(Delta) #this plots the Delta 
    plt.show()

def V(s):
    R=1.65
    return -60.9*math.exp(-(s/R)**2)
    
LHS = 0
RHS = 20
L = RHS-LHS

N = 1000
M = 500

delx = L/N

U = np.zeros(N+1) #This is the A vector
VV = [] #This is the graph of the potential 

for d in range(N+1):
    x = LHS+delx*d
    VV.append(V(x))

BC = 0
def BC1(x):
    return math.exp(x)

psipos = np.zeros(N+1)
psineg=np.zeros(N+1)
psipos[0]=BC
psipos[1]=delx
psineg[N]=BC1(-0.05*RHS)
psineg[N-1]=BC1(-0.05*(RHS+delx))


for l in range(2000):
    for m in range(N+1):
        x = LHS+delx*m
        nmb = (1/(h2m))*(E[l]-V(x))
        U[m]=nmb

    for i in range(1,M+1):
            psipos[i+1] = (-psipos[i-1]*(12+(delx**2)*U[i-1])+2*psipos[i]*((-5*delx**2)*U[i]+12))/((delx**2)*U[i+1]+12)


    for j in range(N-1,M-1,-1):
            psineg[j-1] = (-psineg[j+1]*((delx**2)*U[j+1]+12)+2*psineg[j]*(12-(5*delx**2)*U[j]))/(12+(delx**2)*U[j-1])
    
    negprime = (psineg[M]-psineg[M-1])/delx
    posprime = (psipos[M+1]-psipos[M])/delx

    Delta.append(np.abs(-(negprime/psineg[M])+(posprime/psipos[M]))) #this computes the delta
    delE.append(abs(E[l]-E[l-1]))# this computes Delta En and En+1

l=(Delta.index(min(Delta)))#This finds the energy where the Delta is the minimum
soln(LHS, delx, E, l, M, N, psipos, psineg, Delta, delE, U) #

dict = {"Delta": Delta, "Energy": E, "Delta E": delE}
df = pd.DataFrame(dict)
df.to_csv('solution.csv')


'''

   
    if l >0:

        slope = (Delta[l-1] - Delta[l])/(E[l-1]-E[l])
        x1 = E[l-1]
        y1 = Delta[l-1]
        x2=E[l] 
        y2 = Delta[l]
    
        slope = (y2 - y1) / (x2 - x1)
    
    
        E.append(-y1 / slope + x1)

        x_values = [E[l], E[l-1],E[l+1]]
        y_values = [Delta[l], Delta[l-1],0]
    
        #plt.plot(x_values, y_values, marker='o', linestyle='-')
        #plt.xlabel('x')
        #plt.ylabel('y')
        #plt.title('Linear Function and Root')
        #plt.grid(True)

        tol = abs(E[l]-E[l-1])
        print(tol)
        print(E[l])
        print(Delta[l])'''