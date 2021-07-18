"""
                                                     EE2703 Applied Programming Lab 2021
                                                                 Assignment 5

 Purpose  : Solving, Laplace's Equation for values of the potential and current
 Author   : Pranav Phatak (EE19B105)
 Input    : Values of nx,ny,radius of wire and iterations to perform can be specified using the command line
 Output   : Plots different graphs as specified in the report

"""

import numpy as np
import os,sys
import scipy
import scipy.linalg as s
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3

#get user inputs
try:
    if(len(sys.argv)==5):
        Nx=int(sys.argv[1])
        Ny=int(sys.argv[2])
        radius=int(sys.argv[3])  
        Niter=int(sys.argv[4])
    else:
        Nx=25 # length along x
        Ny=25 # length along y
        radius=8 #radius of central lead
        Niter=1500 #number of iterations
except:
    Nx=25 # size along x
    Ny=25 # size along y
    radius=8 #radius of central lead
    Niter=1500 #number of iterations

#initialize potential matrix
phi_matrix=np.zeros((Nx,Ny),dtype = float)
#Range of points on X and Y axis
x_range,y_range=np.linspace(-0.5,0.5,num=Nx,dtype=float),np.linspace(-0.5,0.5,num=Ny,dtype=float)
Y,X=np.meshgrid(y_range,x_range,sparse=False)   #Make grid from x and y values
phi_matrix[np.where(X**2+Y**2<(0.35)**2)]=1.0   #Make electrode voltage as 1V

#plot initial potential
plt.figure(3)
plt.title("Initial potential 2D contour")
plt.xlabel(r"X$\longrightarrow$")
plt.ylabel(r"Y$\longrightarrow$")
plt.contourf(X,Y,phi_matrix,cmap=plt.get_cmap('coolwarm'))
plt.colorbar()
plt.savefig("3.jpg")
plt.close()

#function to update voltage matrix for every iteration
def new_phi(phi_matrix,phi_old):
    phi_matrix[1:-1,1:-1]=0.25*(phi_old[1:-1,0:-2]+ phi_old[1:-1,2:]+ phi_old[0:-2,1:-1] + phi_old[2:,1:-1])
    return phi_matrix

#function to set the boundary conditions on voltage matrix for every iteration
def voltage_boundary_conditions(phi_matrix):
    phi_matrix[:,Nx-1]=phi_matrix[:,Nx-2] # Right Boundary
    phi_matrix[:,0]=phi_matrix[:,1] # Left Boundary
    phi_matrix[0,:]=phi_matrix[1,:] # Top Boundary
    phi_matrix[Ny-1,:]=0
    center = np.where(X**2+Y**2<(0.35)**2)
    phi_matrix[center]=1.0  #Make electrode voltage as 1V
    return phi_matrix


#function to get exponential fit for the error plot
def fit_exponential(y,Niter,iteration_start):
    log_error = np.log(error)[-iteration_start:]
    X = np.vstack([(np.arange(Niter)+1)[-iteration_start:],np.ones(log_error.shape)]).T
    log_error = np.reshape(log_error,(1,log_error.shape[0])).T
    return s.lstsq(X, log_error)[0]

#boundary conditions
def temperature_boundary_conditions(phi_matrix):
    phi_matrix[:,Nx-1]=phi_matrix[:,Nx-2] # Right Boundary
    phi_matrix[:,0]=phi_matrix[:,1] # Left Boundary
    phi_matrix[0,:]=phi_matrix[1,:] # Top Boundary
    phi_matrix[Ny-1,:]=300.0
    center = np.where(X**2+Y**2<(0.35)**2)
    phi_matrix[center]=300.0    #Make electrode temperature as 300K
    return phi_matrix

#laplaces equation
def new_temperature(temp,tempold,Jx,Jy):
    temp[1:-1,1:-1]=0.25*(tempold[1:-1,0:-2]+ tempold[1:-1,2:]+ tempold[0:-2,1:-1] + tempold[2:,1:-1]+(Jx)**2 +(Jy)**2)
    return temp


#Function to find net error
find_net_error = lambda a,b,Niter : -a/b*np.exp(b*(Niter+0.5))

#Initialize error matrix
error = np.zeros(Niter)
#Iterate and update the voltage matrix

for iteration in range(Niter):
    phi_old = phi_matrix.copy()
    phi_matrix = new_phi(phi_matrix,phi_old)  #Update matrix
    phi_matrix = voltage_boundary_conditions(phi_matrix)    #Set boundary condition
    error[iteration] = np.max(np.abs(phi_matrix-phi_old))    #Error between old and new voltage matrix
    if(error[iteration] == 0):  #Break if error reaches steady state
        break

B1,A1 = fit_exponential(error,Niter,0)          #Expoential fit with all entries
B2,A2 = fit_exponential(error,Niter,500)        #Expoential fit with entries only after 500 iterations    

    
#Plotting best fit and error
#log-log scale
plt.figure(1)
plt.title("Best fit for error on loglog scale")
plt.xlabel(r"Iterations$\longrightarrow$")
plt.ylabel(r"Error$\longrightarrow$")
x = np.asarray(range(Niter))+1
plt.loglog(x,error,label="Error")
plt.loglog(x[::50],np.exp(A1+B1*np.asarray(range(Niter)))[::50],'ro',label="fit1")      #plotting in jumps of 50 as is asked
plt.loglog(x[::50],np.exp(A2+B2*np.asarray(range(Niter)))[::50],'go',label="fit2")
plt.legend()
plt.savefig("1.jpg")
plt.close()
    
#now semilog
plt.figure(2)
plt.title("Best fit for error on semilog scale")
plt.xlabel(r"Iterations$\longrightarrow$")
plt.ylabel(r"Error$\longrightarrow$")
plt.semilogy(x,error,label="Error")
plt.semilogy(x[::50],np.exp(A1+B1*np.asarray(range(Niter)))[::50],'ro',label="fit1")    #plotting in jumps of 50 as is asked    
plt.semilogy(x[::50],np.exp(A2+B2*np.asarray(range(Niter)))[::50],'go',label="fit2")
plt.legend()
plt.savefig("2.jpg")
plt.close()



#plotting cumulative error
iteration=np.arange(100,1501,100)
plt.figure(4)
plt.grid(True)
plt.title(r'Cumulative Error on loglog scale')
plt.xlabel(r"Iterations$\longrightarrow$")
plt.ylabel(r"Net maximum error$\longrightarrow$")
plt.loglog(iteration,np.abs(find_net_error(A2,B2,iteration)),'ro')
plt.savefig("4.jpg")
plt.close()


#plotting 2d contour of final potential
plt.figure(5)
plt.title("Final potential 2D contour")
plt.xlabel(r"X$\longrightarrow$")
plt.ylabel(r"Y$\longrightarrow$")
plt.contourf(Y,X[::-1],phi_matrix)  #Contour plot
x_electrode,y_electrode=np.where(X**2+Y**2<(0.35)**2)   #Points with electrode
plt.plot((x_electrode-Nx/2)/Nx,(y_electrode-Ny/2)/Ny,'ro')  #Mark those points as red
plt.colorbar()
plt.savefig("5.jpg")
plt.close()


#Plotting 3d contour of final potential
fig6=plt.figure(6)     # open a new figure
ax=p3.Axes3D(fig6) # Axes3D is the means to do a surface plot
plt.title('The 3-D surface plot of the final potential')
surf = ax.plot_surface(Y, X, phi_matrix.T, rstride=1, cstride=1, cmap=plt.cm.jet)
plt.savefig("6.jpg")
plt.close()

#finding and plotting current density
Jx,Jy = (1/2*(phi_matrix[1:-1,0:-2]-phi_matrix[1:-1,2:]),1/2*(phi_matrix[:-2,1:-1]-phi_matrix[2:,1:-1]))
plt.figure(7)
plt.title("Vector plot of current flow")
plt.quiver(Y[1:-1,1:-1],-X[1:-1,1:-1],-Jx[:,::-1],-Jy)
x_electrode,y_electrode=np.where(X**2+Y**2<(0.35)**2)   #Points with electrode
plt.plot((x_electrode-Nx/2)/Nx,(y_electrode-Ny/2)/Ny,'ro')  #Mark those points as red
plt.savefig("7.jpg")
plt.close()


#initialize temperature matrix
temp=300 * np.ones((Nx,Ny),dtype = float)

#Iterate and update the temperature matrix
for k in range(Niter):
    tempold = temp.copy()
    temp = new_temperature(temp,tempold,Jx,Jy)   #Update step
    temp = temperature_boundary_conditions(temp)    #Set boundary condition after every iteration


#plotting 3d contour of final temp
fig10=plt.figure(8)
ax=p3.Axes3D(fig10)
plt.title('The 3-D surface plot of the temperature')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('Temperature')
ax.plot_surface(X, Y, temp , rstride=1, cstride=1, cmap=plt.cm.jet,linewidth=0, antialiased=False)
plt.savefig("8.jpg")
plt.close()
