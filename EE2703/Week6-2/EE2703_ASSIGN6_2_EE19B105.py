"""
                       EE2703 - Assignment 6.2 The Laplace Transform 

AUTHOR : Pranav Phatak (EE19B105)
PURPOSE : Using the scipy.signal module to perform analysis of systems with rational polynomial transfer functions
INPUT : NULL
OUTPUT : Plotting the time response of a coupled system of differential equations and a linear electrical circuit which behaves like a low-pass filter
"""


import numpy as np 
import matplotlib.pyplot as plt 
import scipy.signal as sp


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



#Question 1
def input_signal(freq,decay):
    """Transfer function of the given system"""
    
    n = np.poly1d([1,decay])
    d = n*n+freq**2
    return n,d
    
def general_transfer_fcn(wn=1.5,zeta=0,gain=1/2.25):
    """General transfer function"""
    
    n = np.poly1d([wn**2*gain])
    d = np.poly1d([1,2*wn*zeta,wn**2])
    return n,d
    
def LTI(decay,freq=1.5):
    """Find the response to the given system to a decaying cosine."""
    
    input_numerator, input_denominator = input_signal(freq,decay=decay)
    transfer_numerator, transfer_denominator = general_transfer_fcn()

    output_numerator,output_denominator = input_numerator*transfer_numerator, input_denominator*transfer_denominator
    out_s = sp.lti(output_numerator.coeffs, output_denominator.coeffs)

    t = np.linspace(0,50,1000)
    
    return sp.impulse(out_s,None,t)

X1,Y1 = LTI(0.5)    #Decay constant is 0.5 
p1=General_Plotter(r"$t$",r"$x$",r"Q1: System response with decay of $0.5$")
p1.general_plot(X1,Y1)

#Question 2 
X2,Y2 = LTI(0.05)   #Decay constant is 0.05
p2=General_Plotter(r"$t$",r"$x$",r"Q2: System response with decay of $0.05$")
p2.general_plot(X2,Y2)


#Question 3
def input_time(t,decay=0.5,freq=1.5):
    """Exponentially decaying cosine function."""

    u_t = 1*(t>0)
    return np.cos(freq*t)*np.exp(-decay*t) * u_t

time_stamps = np.linspace(0,100,1000)
transfer_function = general_transfer_fcn(wn=1.5,zeta=0,gain=1/2.25)
outputs = []

#Looping through different values of frequency and plotting the output time domain simulation
for freq in np.arange(1.4,1.6,0.05):
    #Simulating
    t,response,_ = sp.lsim(transfer_function,input_time(time_stamps,0.05,freq),time_stamps)
    outputs.append(response)	

p3 = General_Plotter(r"$t$",r"$x$",r"Q3: System responses with variation of input frequency")
p3.general_plot(t,np.array(outputs).T,["Freq = ${:.2f}$".format(f) for f in np.arange(1.4,1.6,0.05)])

w,S,phi=sp.lti(*transfer_function).bode()

p4 = General_Plotter("Frequency in rad/s (log)","Magnitude in dB","Q3: Magnitude plot")
p4.semilogx(w,S)
    
p5 = General_Plotter("Frequency in rad/s (log)","Phase in degrees","Q3: Phase plot")
p5.semilogx(w,phi)
    
#Question 4
def coupled_spring():
    #Laplace function for x and y obtained by decopling the equations
    laplace_function_x = sp.lti(np.poly1d([1,0,2]),np.poly1d([1,0,3,0]))
    laplace_function_y = sp.lti(np.poly1d([2]),np.poly1d([1,0,3,0]))
    time_stamps = np.linspace(0,20,1000)

    #Obtaining the output in time domain for each laplace function
    response_x = sp.impulse(laplace_function_x,None,time_stamps)
    response_y = sp.impulse(laplace_function_y,None,time_stamps)

    p6 = General_Plotter("Time","X","Q4: X vs time")
    p6.general_plot(response_x[0],response_x[1])

    p7 = General_Plotter("Time","Y","Q4: Y vs time")
    p7.general_plot(response_y[0],response_y[1])
	
    p8 = General_Plotter("Time","Displacement","Q4: Responses of coupled system")
    p8.general_plot(t,np.array([response_x[1],response_y[1]]).T,legend_txt=[r"$x(t)$", r"$y(t)$"])
    
coupled_spring()

#Question 5
def two_port():
	# Find the transfer function of the given circuit
    R = 100
    L = 1e-6
    C = 1e-6

    wn = 1/np.sqrt(L*C) # natural frequency
    Q = 1/R * np.sqrt(L/C) # quality factor
    zeta = 1/(2*Q) # damping constant

	# transfer function
    n,d = general_transfer_fcn(gain=1,wn=wn,zeta=zeta)

    # make system
    transfer_function = sp.lti(n,d)

    # get bode plots
    w,S,phi=transfer_function.bode()
    
    p9 = General_Plotter("Frequency in rad/s (log)","Magnitude in dB","Q5: Magnitude plot ")
    p9.semilogx(w,S)
    
    p10 = General_Plotter("Frequency in rad/s (log)","Phase in degrees","Q5: Phase plot ")
    p10.semilogx(w,phi)

    return transfer_function

transfer_function = two_port()

#Question 6
def multiple_frequencies(transfer_function):
	
	#Simulating and plotting Q5 for given input for time uptill 30us 
	time = np.linspace(0,30*0.000001,1000)
	vi = np.multiply(np.cos(1000*time)-np.cos(1000000*time),np.heaviside(time,0.5))
	_,output_time,_ = sp.lsim(transfer_function,vi,time)
	p11 = General_Plotter("Time$\longrightarrow$","Voltage$\longrightarrow$","Q6: Voltage uptill 30mus")
	p11.general_plot(time,output_time)
	
	#Simulating and plotting Q5 for given input for time uptill 1ms
	time = np.linspace(0,10*0.001,100000)
	vi = np.multiply(np.cos(1000*time)-np.cos(1000000*time),np.heaviside(time,0.5))
	_,output_time,_ = sp.lsim(transfer_function,vi,time)
	p12 = General_Plotter("Time$\longrightarrow$","Voltage$\longrightarrow$","Q6: Voltage uptill 10ms")
	p12.general_plot(time,output_time)

multiple_frequencies(transfer_function)


