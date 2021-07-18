"""
                EE2703 - Assignment 6
           Author : EE19B105, Pranav Phatak  

The code can be called in the commandline in the following format:
python Assgn6.py n=value M=value Msig=Value nk=Value u0=Value p=Value seed=Value 
where all the arguments are optional and value is some number                

Inputs: M ( number of electrons injected ) , n ( length of tubelight ) , nk ( time steps to be taken ) , u0 ( Velocity after which collision is possible )
        p ( probability that an electron at velocity higher than u0 collides ) , Msig ( value of standard deviation for normal distribution of number of electrons injected )
        
Outputs: No. of electrons at a position vs position graph
         Intensity of light at a position vs position graph
         Velocity of electron at a position vs position graph
         Tabular data of Intensity vs position                       
"""
import os
import numpy as np    
import matplotlib.pyplot as plt
import sys
from pylab import *
import pandas as pd    # pandas for showing in tabular form


M = 5                                          # number of electrons injected per turn.
n=100                                           # spatial grid size.
nk=500                                          # number of turns to simulate.
u0= 5                                           # threshold velocity.
p=0.25                                          # probability that ionization will occur                
Msig=1                                          # taking sigma of normal distribution to be 1

# Taking inputs from sys.argv if given, if not give use the values mentioned above

args = sys.argv[1:]  
argList = {}         # Dictionary declared
for i in args:
    i = i.split('=')    # 'varname=value' string is split into varname and value
    argList[i[0]] = float(i[1])   # The dictionary is updated with values from commandline


if 'n' in argList:
    n = argList['n']
if 'M' in argList:
    M = argList['M']
if 'nk' in argList:
    nk = argList['nk']
if 'u0' in argList:
    u0 = argList['u0']
if 'p' in argList:
    p = argList['p']
if 'Msig' in argList:
    Msig = argList['Msig']


xx = np.zeros(n*M)      #Electron position
u = np.zeros(n*M)       #Electron velocity  
dx = np.zeros(n*M)      #Electron displacement per time step

I = []          #Stores value of intensity of emitted light at every time-step
V = []          #Stores value of electron velocity at every time-step
X = []          #Stores value of electron position at every time-step    


for i in range(1,nk):
    ii = np.where(xx>0)[0]
    
    dx[ii] = u[ii] + 0.5
    xx[ii] = xx[ii] + dx[ii]
    u[ii] = u[ii] + 1.0
    
    npos = np.where(xx>n)[0]
    xx[npos] = 0.0
    dx[npos] = 0.0
    u[npos] = 0.0
    
    kk=np.where(u>=u0)[0] 
    # array which gives indices of electrons which can cause collision
    ll=np.where(rand(len(kk))<=p)[0]
    kl=kk[ll]                   
    # array to store indices of electrons which collided to give a photon  
    
    P = np.random.rand(len(kl)) 
    # value to reset position by when collision occurs
    xx[kl] = xx[kl]-np.multiply(dx[kl],P)                 
    # position reset by small factor since collision can take place at any time  
    u[kl] = 0                                             
    # inelastic collision implies velocity reset to 0
    
    I.extend(xx[kl].tolist())
    
    m=np.random.randn()*Msig+M
    empty_xpos = np.where(xx==0)[0]
    electrons_generated = min(len(empty_xpos),(int)(m))
    xx[empty_xpos[0:electrons_generated]] = 1.0
    u[empty_xpos[0:electrons_generated]] = 0.0
    
    X.extend(xx[ii].tolist())
    V.extend(u[ii].tolist()) 
       
     
# Plotting graphs now

# Population plot for electron density
plt.figure(0)
plt.hist(X,histtype='bar', bins=np.arange(1,n,n/100),ec='black',alpha=0.5) 
plt.title('Population Plot')
plt.xlabel('x')
plt.ylabel('No. of electrons')
plt.savefig("Electrons_vs_x.jpg")
plt.close()

# Population plot for intensity of emitted light
plt.figure(1)
plt.hist(I,histtype='bar', bins=np.arange(1,n,n/100),ec='black',alpha=0.5)
plt.title('Intensity Plot')
plt.xlabel('x')
plt.ylabel('No. of photons emitted ($\propto$ Intensity)')
plt.savefig("Intensity_vs_x.jpg")    
plt.close()

# Plot for Electron phase space ( pos vs velocity )
plt.figure(2)
plt.plot(X,V,'x')
plt.title('Electron phase space')
plt.xlabel('x')
plt.ylabel('velocity')
plt.savefig("Pos_vs_Velocity.jpg")
plt.close()

# Tabulating data for intensity vs position
bins = plt.hist(I,bins=np.arange(1,n,n/100))[1]    # Bin positions are obtained
count = plt.hist(I,bins=np.arange(1,n,n/100))[0]   # Population counts obtained
xpos = 0.5*(bins[0:-1] + bins[1:])     # As no. of end-points of bins would be 1 more than actual no. of bins, the mean of bin end-points are used to get population of count a particular bin
df = pd.DataFrame()   # A pandas dataframe is initialized to do the tabular plotting of values.
df['Xpos'] = xpos
df['count'] = count

base_filename = 'values.txt'
with open(base_filename,'w') as outfile:
    df.to_string(outfile)

