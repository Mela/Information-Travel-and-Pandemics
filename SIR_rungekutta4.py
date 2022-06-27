#!/usr/bin/env python
# coding: utf-8

# In[39]:


import numpy as np


# In[57]:


#Parameters for the modulation
alpha = 0.1 #information transmission rate
beta = 0.1
gamma = 0.1 #Recovery rate
N =  100 #population size
lam = 0.4 #information fading


# In[70]:


#Starting Parameters
t0 = 0  #start time
tmax = 100 #maximum time
S0 = 100 #number of susceptible
I0 = 1  #number of infected at start
R0 = 0  #number of recovered at start
N0 = [0]*2
y0 = [S0,I0,R0,N0,N0] #starting conditions
n_steps = 100 #number of steps


# In[89]:


def SIR_N(t,arr):
    S = arr[0]
    I = arr[1]
    R = arr[2]
    N_before = arr[3]
    Niarr = arr[4]
    
    #calculate N_{<i}=sum_{j=0}^{i-1}N_j
    Nsmalleri = [0]
    for i in Niarr:
        x = i+sum(Nsmalleri)
        Nsmalleri.append(x)
    Nsmalleri.pop()
    #differential equation
    for i in range(0,len(Niarr)):
        Niarrnew = -alpha*Niarr[i]/N*Nsmalleri[i]+alpha*N_before[i]/N*(N-Nsmalleri[i])-lam*(Niarr[i]-N_before[i])                                                               
    Snew = -beta*S/N*I
    Inew = beta*S/N*I-gamma*I
    Rnew = gamma*I
    print([Snew, Inew,Rnew,Niarr,Niarrnew])
    return([Snew, Inew,Rnew,Niarr,Niarrnew])


# In[90]:


def rungekutta4(t0, tmax,y0,n):
    """
    input starting t0, maximum t tmax, start condition array y0
    and number of steps n
    """
    h = (tmax-t0)/n  #step size
    t = t0
    y = y0
    yarr = y0
    tarr = t0
    for i in range(n):
        print(y)
        k1a = (SIR_N(t,y))
        #print(k1a)
        k1 = [print(z) for z in k1a]
        k2 = h * (SIR_N((t+h/2), [z+x/2 for z in y and x in k1]))
        k3 = h * (SIR_N((t+h/2), (y+k2/2)))
        k4 = h * (SIR_N((t+h), (y+k3)))
        k = (k1+2*k2+2*k3+k4)/6
        ynew = y + k
        y = ynew
        t = t+h
        yarr.append(y)
        tarr.append(t)
    
    return(tarr,yarr)


# In[91]:


t, SIRN_arr = rungekutta4(t0,tmax,y0,n_steps)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




