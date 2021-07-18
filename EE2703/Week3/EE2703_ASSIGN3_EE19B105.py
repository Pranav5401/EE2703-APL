"""
                EE2703 Applied Programming Lab 2021
                          Assignment 3

Purpose  : To determine how noise affects the estimation of the coefficients of a function
Author   : Pranav Phatak (EE19B105)
Input    : Program reads data from fitting.dat as input
Output   : Plots of different graphs as shown in the report

"""

from pylab import *
import scipy.special as sp
import warnings


warnings.filterwarnings("ignore")                   #To ignore the warnings so that the terminal doesn't become crowded when program is compiled

"""
Defining functions to solve all the questions, will be called from main

"""

def g(t,A,B):                                               #Function for calculating g(t,A,B)
    return(A*sp.jn(2,t)+B*t)                                #sp.jn denotes the bessel function

def calc_MSE(predictions, true_values):                         #Function to calculate Mean Square Error
    return ((predictions - true_values) ** 2).mean()


def make_M(t):                                              #Function for making matrix M
    return c_[sp.jn(2,t),t]


def solve_for_AB(M,b):                                      #Estimating the value of A and B
    return linalg.lstsq(M,b)                                #Inbuilt function which solves matrices to give best possible result     


def plot_all_data(t,values,sigma,ground_truth):             #Function to plot ground_truth values and the data read from fitting.dat
    values=c_[values,ground_truth]
    figure(0)
    plot(t,values)
    
    title("Q4:Data to be fitted to theory")
    xlabel(r"t $\rightarrow$")
    ylabel(r"$f(t)+Noise \rightarrow$")
    sigma=[f"Standard Deviation={round(i,4)}" for i in sigma]
    legend(list(sigma)+["True Value"])
    grid(True)
    
    savefig("Ques3_4.jpg")                                  #Saving plot    
    close()


def plot_errorbars(t,values,ground_truth):                  #Function to plot error bars
    figure(1)
    std_dev = std(values[:,0]-ground_truth)
    
    title(f"Q5: Data Points for stddev={round(std_dev,3)} along with exact function")
    xlabel(r"t $\rightarrow$")
    ylabel(r"f(t) $\rightarrow$")
    errorbar(t[::5],values[:,0][::5],std_dev,fmt='ro')
    plot(t,ground_truth)
    legend(["f(t)","Error Bar"])
    
    savefig("Ques5.jpg")                                    #Saving plot
    close()


def plot_contours(values,t):                                #Function to plot contours
    figure(2)
    A_range=arange(0,2.1,0.1)
    B_range=arange(-0.2,0.01,0.01)
    epsilon_matrix = zeros((len(A_range),len(B_range)))
    
    for count_A,A in enumerate(A_range):
        for count_B,B in enumerate(B_range):
            epsilon_matrix[count_A][count_B] = calc_MSE(values[:,0],g(t,A,B))
    
    contour_obj = contour(A_range,B_range,epsilon_matrix,levels=np.arange(0,20*0.025,0.025),cmap="magma")
    clabel(contour_obj,contour_obj.levels[0:5],inline=True,fontsize=10)
    
    title('Q8: Contour plot of $\epsilon_{ij}$')
    xlabel(r'A $\rightarrow$',size=12)
    ylabel(r'B $\rightarrow$',size=12)
    plot(1.05, -0.105,'ro', label = 'Exact Value')
    annotate("Exact Value",xy = [0.8,-0.100])
    
    savefig("Ques8.jpg")                                    #Saving plot
    close()


def plot_variation_of_error(sigma,error_a,error_b):         #Function to plot variation of error in A,B with change in sigma in normal scale
    figure(3)
    plot(sigma,error_a,'r--')
    scatter(sigma,error_a)
    plot(sigma,error_b, 'b--')
    scatter(sigma,error_b)
    
    legend(["Aerr","Berr"])
    title("Q10: Variation of error with noise")
    xlabel(r'$\sigma_{n}\rightarrow$',size=10)
    ylabel(r'Absolute Error $\rightarrow$',size=10)
    
    savefig("Ques10.jpg")                                  #Saving plot     
    close()


def plot_logvslog_variation_of_error(sigma,error_a,error_b):     #Function to plot variation of error in A,B with change in sigma in log vs log scale
    figure(5)
    stem(sigma,error_a,'ro')
    loglog(sigma,error_a,'ro')
    loglog(sigma,error_b,'bo')
    stem(sigma,(error_b),'bo')
    
    title("Q11: loglog Variation of error with noise")
    xlabel(r'$\sigma_{n}\rightarrow$',fontsize=15)
    ylabel(r'Error$\rightarrow$',fontsize=15)   
    legend(["Aerr","Berr"])
    
    savefig("Ques11.jpg")                                 #Saving plot              
    close()


"""
Main function starts here

"""


try:                                    #Question 2 (loading the data)
    data=loadtxt("fitting.dat")
except:
    print("No such file exists. Please run generate_data.py to generate fitting.py")
    exit()

t_values=data[:,0]                   #Parsing the data
func_values=data[:,1:]

ground_truth = g(t_values,1.05,-0.105)   #Calculating ground truth values for Question 6
    
sigma=logspace(-1,-3,9)              #noise stdev
plot_all_data(t_values,func_values,sigma,ground_truth)          #Question 3 and 4

plot_errorbars(t_values,func_values,ground_truth)                   #Question 5

    
M=make_M(t_values)                                             #Question 6
AB=np.array([1.05,-0.105])
    
    
print("The mean square error between M*AB and g(t,A,B) is: ",calc_MSE(matmul(M,AB),g(t_values,1.05,-0.105)))         #Confirm that M*AB and the function g(t,A,B) are equal


plot_contours(func_values,t_values)                                 #Question 7 and 8

AB_gt_values=solve_for_AB(M,ground_truth)[0]        #The ground truth values of A and B 
print("The best estimate of A and B for ground truth values are: ", round(AB_gt_values[0],3),round(AB_gt_values[1],3))          #Question 9
    
    
a_error=zeros(9)        
b_error=zeros(9)
total_error=zeros(9)
    
for i in range(9):                              #Iterating through all the columns of fitting.dat file and calculating the absolute error between true values of A and B vs the predicted value
    temp = solve_for_AB(M,func_values[:,i])
    prediction=temp[0]
    diff=temp[1]
    a_error[i],b_error[i] = abs(prediction[0]-AB_gt_values[0]),abs(prediction[1]-AB_gt_values[1])
    total_error[i] = diff
    
plot_variation_of_error(sigma,a_error,b_error)              #Question 10
    
plot_logvslog_variation_of_error(sigma,a_error,b_error)          #Question 11


print("Plots are saved in the same directory as the code")

