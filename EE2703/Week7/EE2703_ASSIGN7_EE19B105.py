"""
                       EE2703 - Assignment 7 Sympy

AUTHOR : Pranav Phatak (EE19B105)
PURPOSE : Using the sympy module to perform analysis of a LPF circuit and a HPF circuit
INPUT : NULL
OUTPUT :9 Graphs
"""

import numpy as np 
import pylab as plt
import scipy.signal as sp
from sympy import *
import warnings


warnings.filterwarnings("ignore")                   #To ignore the warnings so that the terminal doesn't become crowded when program is compiled



#Since we have many plots will define a class General Plotter to save code lines and increase efficiency
class General_Plotter():
    ''' Class used for plotting different plots. Shortens the code by quite a bit'''
    
    fig_num=0   #Defined static variable for the figure number
    def __init__(self,xlabel,ylabel,title):
        ''' xlabel,ylabel,title are used in every graph''' 

        self.xlabel=xlabel
        self.ylabel=ylabel
        self.title=title
        self.fig=plt.figure(self.__class__.fig_num)
        self.__class__.fig_num+=1
    
    def general_funcs(self,ax):
        ''' General functions for every graph'''
        
        ax.set_ylabel(self.ylabel)
        ax.set_xlabel(self.xlabel)
        ax.set_title(self.title)
        self.fig.savefig(self.title+".png")

    def general_plot(self,X,Y,legend_txt=None):
        ''' General Plot'''
        axes=self.fig.add_subplot(111)
        axes.plot(X,Y)
        if legend_txt is not None:
            axes.legend(labels=legend_txt)
        self.general_funcs(axes)
        
    
    def semilogx(self,X,Y):
        ''' Semilog scale''' 
        axes=self.fig.add_subplot(111)
        axes.semilogx(X,Y)
        axes.grid()
        self.general_funcs(axes)

    def loglog(self,X,Y):
        ''' Loglog scale'''
        axes=self.fig.add_subplot(111)
        axes.loglog(X,Y)
        axes.grid()
        self.general_funcs(axes)
        
#1st (Low Pass Filter Circuit)        
def lowpass(R1,R2,C1,C2,G,Vi):
    s=symbols('s')
    A=Matrix([[0,0,1,-1/G],[-1/(1+s*R2*C2),1,0,0],[0,-G,G,1],[-1/R1-1/R2-s*C1,1/R2,0,s*C1]])
    b=Matrix([0,0,0,Vi/R1])        
    V=A.inv()*b
    return (A,b,V)
    
#2nd (High Pass Filter Circuit)    
def highpass(R1,R2,C1,C2,G,Vi):
    s=symbols('s')
    A = Matrix([[0,0,1,-1/G],[-1/(1+1/(s*R2*C2)),1,0,0], [0,-G,G,1],[0-s*C1-s*C2-1/R1,s*C2,0,1/R1]])
    b = Matrix([0,0,0,Vi*s*C1])
    V = A.inv()*b
    return (A,b,V)
    
    
#Question1
H = sp.lti([-0.0001586],[2e-14,4.414e-9,0.0002])
w,HS,phi = H.bode()

plt.subplot(2,1,1)
plt.semilogx(w,HS)
plt.xlabel(r"$\omega$")
plt.ylabel(r"$|H(s)|$")
plt.subplot(2,1,2)
plt.semilogx(w,phi)
plt.xlabel(r"$\omega$")
plt.ylabel(r"$\angle(H(s))$")
plt.suptitle("Bode plot of transfer function of Low Pass Filter")
plt.savefig("Q1: Bode plot of transfer function of Low Pass Filter")
plt.close()

A,b,V = lowpass(10000,10000,1e-9,1e-9,1.586,1)
Vo = V[3]
t1 = np.linspace(0,0.001,10000)
v = sp.lti(np.poly1d([0.0001586]),np.poly1d([2e-14,4.414e-9,0.0002,0]))
stepl = sp.impulse(v,None,t1)
p3 = General_Plotter("t", "$V_o(step)$", "Q1: Step Response of Low Pass Filter")
p3.general_plot(stepl[0],stepl[1])

#Question2
t = np.linspace(0,0.001,200000)
u1 = np.sin(2000*np.pi*t) 
u2 = np.cos(200000*np.pi*t) 
Vo = simplify(Vo)
n,d = fraction(Vo)
p4 = General_Plotter("t","$V_i(t)$","Q2: Components of the input signal $V_i$")
p4.general_plot(t,np.array([u1,u2]).T,legend_txt=["sin($2\cdot10^3\pi$t)","cos($2\cdot10^6\pi$t)"])    
    
H = sp.lti([0.0001586],[2e-14,4.414e-9,0.0002])
t,y,_ = sp.lsim(H,u1+u2,t)
tn = np.linspace(0,0.00001,200000)
tn,z,_ = sp.lsim(H,u2,tn)

p5 = General_Plotter("t","Voltage","Q2: Vo(t) - Output of Low Pass Filter with both frequency inputs")
p5.general_plot(t,np.array([u1+u2,y]).T,legend_txt=["Vi(t) = sin($2\cdot10^3\pi$t) + cos($2\cdot10^6\pi$t)", "Vo(t)"])

p6 = General_Plotter("t","Voltage","Q2: Vo(t) - Output of Low Pass Filter with single frequency input")
p6.general_plot(tn,np.array([u2,z]).T,legend_txt=["Vi(t) = cos($2\cdot10^6\pi$t)", "Vo(t)"])
    
#Question3
Hs = sp.lti([1.586e-9,0,0],[2e-9,4.414e-4,20.0])
W,Hhp,phi_Hs = Hs.bode()

plt.subplot(2,1,1)
plt.semilogx(W,Hhp)
plt.xlabel(r"$\omega$")
plt.ylabel(r"$|H(s)|$")
plt.subplot(2,1,2)
plt.semilogx(W,phi_Hs)
plt.xlabel(r"$\omega$")
plt.ylabel(r"$\angle(H(s))$")
plt.suptitle("Bode plot of transfer function of High Pass Filter")
plt.savefig("Q3: Bode plot of transfer function of High Pass Filter")
plt.close()

#Question5
A,b,V = highpass(10000,10000,1e-9,1e-9,1.586,1)
Vo = V[3]
Vo = simplify(Vo)
v = sp.lti(np.poly1d([1.586e-9,0]),np.poly1d([2e-9,4.414e-4,20.0]))
t2 = np.linspace(0,0.001,10000)
steph = sp.impulse(v,None,t2)
p9 = General_Plotter("t","$V_o(step)$","Q5: Step reponse of High Pass Filter")
p9.general_plot(steph[0],steph[1])

#Question4
ta = np.linspace(0,5,200000)
u1 = np.sin(20*np.pi*ta)*np.exp(-5*ta)
ta,y,_ = sp.lsim(Hs,u1,ta)
p10 = General_Plotter("t","Voltage","Q4: Response of HPF to Low-frequency damped sinusoid")
p10.general_plot(ta,np.array([u1,y]).T,legend_txt=["$V_i(t)$ = sin($20\pi$t)$e^{-5t}$","$V_o(t)$"])

tb = np.linspace(0,0.001,200000)
u2 = np.sin(2000000*np.pi*t)*np.exp(-5000*t)
tb,y,_ = sp.lsim(Hs,u2,tb)
p11 = General_Plotter("t","Voltage","Q4: Response of HPF to High-frequency damped sinusoid")
p11.general_plot(tb,np.array([u2,y]).T,legend_txt=["$V_i(t)$ = sin($2\cdot10^6\pi$t)$e^{-5000t}$","$V_o(t)$"])    
