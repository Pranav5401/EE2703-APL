"""
                EE2703 - Assignment 4
                Pranav Phatak
                EE19B105
Inputs : None from user, will create the functions using lambda function to find fourier series coefficients for exp(x) and cos(cos(x))
Outputs : Error in integral and matrix solving methods for finding fourier series coefficients, and multiple graphs as askedin assignment
"""

import math
import scipy.integrate as integ
from pylab import *
import scipy
import scipy.special as sp
import numpy as np
import warnings


warnings.filterwarnings("ignore")                   #To ignore the warnings so that the terminal doesn't become crowded when program is compiled


f1 = lambda x: np.exp(x)                      #Defining the 2 given functions
f2 = lambda x: np.cos(np.cos(x))

f1_periodic = lambda x: np.exp(x%(2*math.pi))   #Defining the periodically extended exponential function which we use for fourier coefficients

#Plotting aperiodic e^x and periodically extended e^x in log scale from x = -2pi to +4pi by splitting the x range into 300 points from -2pi to 4pi 
def plotting_exponentials():
    x = np.linspace(-2*math.pi,4*math.pi,300)     
    figure(1)
    semilogy(x,f1(x),'m',label='Aperiodic exponential function in log scale')      
    semilogy(x,f1_periodic(x),'g',label='Periodically extended exponential function in log scale')

    grid(True)
    title("Q1 : Plotting original functions")
    ylabel(r'$log(e^{x})$',fontsize=12)
    xlabel(r'$x$',fontsize=12)
    legend()
    
    savefig("Ques1_1.jpg")
    close()

#Plotting cos(cos(x)) in log scale from x = -2pi to +4pi by splitting the x range into 300 points from -2pi to 4pi
def plotting_coscosx():
    x = np.linspace(-2*np.pi,4*np.pi,300)
    figure(2)
    plot(x,f2(x), 'y')
    
    grid(True)
    title(r'y = cos(cos(x))')
    ylabel(r'$cos(cos(x))$', fontsize = 12)
    xlabel(r'$x$', fontsize = 12)
    
    savefig("Ques1_2.jpg")
    close()
                            

#Returns the first k ak and (k-1) bk fourier series coefficients for a function f
def FSC(f,k):                       
    coeff=[]                        #List which stores all the coefficients
    a=[]                            #List of only an's                            
    b=[]                            #List of only bn's
                            
    u = lambda x, n: f(x)*math.cos(n*x)
    v = lambda x, n: f(x)*math.sin(n*x)
    
    b.append(0)                                                                     #Since there is no b(0) will make it 0    
    a.append((1/(2*math.pi))*integ.quad(u, 0, 2*math.pi, args=0)[0])
    coeff.append((1/(2*math.pi))*integ.quad(u, 0, 2*math.pi, args=0)[0])            #Solve for n=0
    for n in range(1,k):
        a.append((1/math.pi)*integ.quad(u, 0, 2*math.pi, args=n)[0])    
        coeff.append((1/math.pi)*integ.quad(u, 0, 2*math.pi, args=n)[0])
        
        b.append((1/math.pi)*integ.quad(v, 0, 2*math.pi, args=n)[0])   
        coeff.append((1/math.pi)*integ.quad(v, 0, 2*math.pi, args=n)[0]) 
    
    return coeff,a,b 

coeff_f1,a_f1,b_f1 = FSC(f1,26)
coeff_f2,a_f2,b_f2 = FSC(f2,26)

#Plots the coefficients in semilogy and loglog scale
def plotting_coefficients():                                                #In all the graphs b(0) will be 0, explained why this is done in report
    figure(3)
    semilogy(np.abs(a_f1), 'ro')
    semilogy(np.abs(b_f1), 'bo')
        
    grid(True)
    title(r'Magnitudes of coefficients in log scale for e^x with $a_n$ in red and $b_n$ in blue', fontsize = 10)
    ylabel(r'log(coeff)')
    xlabel(r'$n$')
    savefig("Ques3_1.jpg")
    close()
    
    
    figure(4)
    loglog(np.abs(a_f1), 'ro')
    loglog(np.abs(b_f1), 'bo')
    
    grid(True)
    title(r'Magnitudes of coefficients in loglog scale for e^x with $a_n$ in red and $b_n$ in blue', fontsize = 10)
    ylabel(r'log(log(coeff))')
    xlabel(r'$n$')
    savefig("Ques3_2.jpg")
    close()
    
    
    figure(5)
    semilogy(np.abs(a_f2), 'ro')
    semilogy(np.abs(b_f2), 'bo')
    
    grid(True)
    title(r'Magnitudes of coefficients in log scale for cos(cos(x)) with $a_n$ in red and $b_n$ in blue', fontsize = 10)
    ylabel(r'log(coeff)')
    xlabel(r'$n$')
    savefig("Ques3_3.jpg")
    close()
    
    
    figure(6)
    loglog(np.abs(a_f2), 'ro')
    loglog(np.abs(b_f2), 'bo')
    
    grid(True)
    title(r'Magnitudes of coefficients in loglog scale for cos(cos(x)) with $a_n$ in red and $b_n$ in blue', fontsize = 10)
    ylabel(r'log(log(coeff))')
    xlabel(r'$n$')
    savefig("Ques3_4.jpg")
    close()



def matrix_method(A,f):                              #Function which return all coefficients and distinct lists for a and b as well.
    coeff = scipy.linalg.lstsq(A,f(x))[0]            #using lstsq for finding value of matrix c that best fits the equation 
    coeff_a = []
    coeff_b = []

    coeff_a.append(coeff[0])
    coeff_b.append(0)
    for i in range(1,51,2):
	    coeff_a.append(coeff[i])

    for i in range(2,51,2):
	    coeff_b.append(coeff[i])

    return coeff, coeff_a, coeff_b
    

#Least Squares Approach
x = np.linspace(0,2*pi,401)
x=x[:-1]  
A = np.zeros((400,51)) 
A[:,0]=1 
for k in range(1,26):
    A[:,2*k-1] = np.cos(k*x) 
    A[:,2*k] = np.sin(k*x) 
#Matrix A (51 x 400) has been defined 

coeff1, coeff1_a , coeff1_b = matrix_method(A,f1)       #Calling above function
coeff2, coeff2_a , coeff2_b = matrix_method(A,f2)
    
    
def plotting_comparing_coeff():
    figure(7)
    fig, axs = plt.subplots(2)                     #To divide axis into 2 to plot an on one and bn on another, similarly done in all 4 plots
    
    axs[0].semilogy(np.abs(coeff1_a), 'bo', label = 'Least Squares Approach')
    axs[0].semilogy(np.abs(a_f1), 'go', label = 'Integration Approach')
    axs[1].semilogy(np.abs(coeff1_b), 'bo', label = 'Least Squares Approach')
    axs[1].semilogy(np.abs(b_f1), 'go', label = 'Integration Approach')
    fig.suptitle(r'Magnitudes of coefficients in log scale for e^x')
    
    axs[0].set(ylabel=r'log(|$a_n$|)$\longrightarrow$',xlabel=r'$n\longrightarrow$')
    axs[1].set(ylabel=r'log(|$b_n$|)$\longrightarrow$',xlabel=r'$n\longrightarrow$')
    axs[0].grid()
    axs[1].grid()
    axs[0].legend()
    axs[1].legend()
    savefig("Ques5_1.jpg")
    close()

    
    figure(8)
    fig, axs = plt.subplots(2)
    
    axs[0].semilogy(np.abs(coeff2_a), 'bo', label = 'Least Squares Approach')
    axs[0].semilogy(np.abs(a_f2), 'go', label = 'Integration Approach')
    axs[1].semilogy(np.abs(coeff2_b), 'bo', label = 'Least Squares Approach')
    axs[1].semilogy(np.abs(b_f2), 'go', label = 'Integration Approach')
    fig.suptitle(r'Magnitudes of coefficients in log scale for cos(cos(x))')
    
    axs[0].set(ylabel=r'log(|$a_n$|)$\longrightarrow$',xlabel=r'$n\longrightarrow$')
    axs[1].set(ylabel=r'log(|$b_n$|)$\longrightarrow$',xlabel=r'$n\longrightarrow$')
    axs[0].grid()
    axs[1].grid()
    axs[0].legend()
    axs[1].legend()
    savefig("Ques5_2.jpg")
    close()


    figure(9)
    fig, axs = plt.subplots(2)
    
    axs[0].loglog(np.abs(coeff1_a), 'bo', label = 'Least Squares Approach')
    axs[0].loglog(np.abs(a_f1), 'go', label = 'Integration Approach')
    axs[1].loglog(np.abs(coeff1_b), 'bo', label = 'Least Squares Approach')
    axs[1].loglog(np.abs(b_f1), 'go', label = 'Integration Approach')
    fig.suptitle(r'Magnitudes of coefficients in log scale for e^x')
    
    axs[0].set(ylabel=r'log(|$a_n$|)$\longrightarrow$',xlabel=r'$n\longrightarrow$')
    axs[1].set(ylabel=r'log(|$b_n$|)$\longrightarrow$',xlabel=r'$n\longrightarrow$')
    axs[0].grid()
    axs[1].grid()
    axs[0].legend()
    axs[1].legend()
    savefig("Ques5_3.jpg")
    close()


    figure(10)
    fig, axs = plt.subplots(2)
    
    axs[0].loglog(np.abs(coeff2_a), 'bo', label = 'Least Squares Approach')
    axs[0].loglog(np.abs(a_f2), 'go', label = 'Integration Approach')
    axs[1].loglog(np.abs(coeff2_b), 'bo', label = 'Least Squares Approach')
    axs[1].loglog(np.abs(b_f2), 'go', label = 'Integration Approach')
    fig.suptitle(r'Magnitudes of coefficients in log scale for cos(cos(x))')
    
    axs[0].set(ylabel=r'log(|$a_n$|)$\longrightarrow$',xlabel=r'$n\longrightarrow$')
    axs[1].set(ylabel=r'log(|$b_n$|)$\longrightarrow$',xlabel=r'$n\longrightarrow$')
    axs[0].grid()
    axs[1].grid()
    axs[0].legend()
    axs[1].legend()
    savefig("Ques5_4.jpg")
    close()

#Analysing the deviation of the two approaches 

error_f1 = np.abs(coeff_f1 - coeff1)
error_f2 = np.abs(coeff_f2 - coeff2)
max_error_f1 = np.max(error_f1)
max_error_f2 = np.max(error_f2)

print("Maximum error for coefficients of e^x is ",max_error_f1)
print("Maximum error for coefficients of cos(cos(x)) is ",max_error_f2)

#Function to plot the values of coefficients found by the 2 methods to see the deviations
def plotting_convergence():                 
    fourier_f1 = np.matmul(A,coeff1)
    fourier_f2 = np.matmul(A,coeff2)

    figure(11)
    semilogy(fourier_f1, 'm', label = 'Fourier representation')
    semilogy(f1_periodic(x), 'c', label = 'Original function')

    grid(True)
    title(r'Convergence of Fourier Series representation to actual function for e^x')
    ylabel(r'Value in log scale')
    xlabel(r'$x$')
    legend()
    savefig("Ques7_1.jpg")
    close()


    figure(12)
    semilogy(fourier_f2, 'm', label = 'Fourier representation')
    semilogy(f2(x), 'b', label = 'Original function')
    
    grid(True)
    title(r'Convergence of Fourier Series representation to actual function for cos(cos(x))', fontsize = 10)
    ylabel(r'Value in log scale', fontsize = 8)
    xlabel(r'$x$')
    legend()
    savefig("Ques7_2.jpg")
    close()

plotting_exponentials()     #Q1 plots function call
plotting_coscosx()

plotting_coefficients()     #Q3 plots function call

plotting_comparing_coeff()  #Q5 plots function call

plotting_convergence()      #Q7 plots function call
              
