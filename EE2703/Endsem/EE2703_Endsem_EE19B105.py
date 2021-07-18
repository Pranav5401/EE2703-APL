"""
                                                     EE2703 Applied Programming Lab 2021
                                                                 Endsem

 Purpose  : Solving for Magnetic Field along the axis of current carrying loop using Vector Potential Curl
 Author   : Pranav Phatak (EE19B105)
 Input    : Values of N,radius of wire can be specified using the command line
 Output   : Plots of graphs of B and current in wire, least error fit to given data comparing to Bz = cz^b

"""
import numpy as np
from pylab import *
import warnings


warnings.filterwarnings("ignore")                   #To ignore the warnings so that the terminal doesn't become crowded when program is compiled


#get user inputs for N, radius if no input given take N = 100, radius =10 as mentioned in qs
try:
    if(len(sys.argv)==3):
        N=int(sys.argv[1])
        radius=int(sys.argv[2])  
    else:
        N=100         #Number used to make sections of wire 
        radius=10     #radius of wire
except:
    N=100             #Number used to make sections of wire
    radius=10         #radius of wire
    
#Fixed constants
mu = 4*np.pi*pow(10,-7)         
k = 1/radius

#Making the vectors of r and current in x,y,z form
phi = linspace(0,2*pi,N)                   
rl = c_[radius*cos(phi),radius*sin(phi),zeros(N)]         #rl[l] gives expected rl
dphi = (2*np.pi/N)*c_[-sin(phi),cos(phi),zeros(N)]        #dphi vector
dphi_abs = 2*np.pi/N                                      #dphi magnitude
dl = radius*dphi                                          #using that dl = rdphi 
I = c_[(-abs(cos(phi))*sin(phi))*4*np.pi/mu,abs(cos(phi))*cos(phi)*4*np.pi/mu,zeros(N)]
I_plot = I/(4*np.pi/mu)     #scaling all I factors to get better quiver
 
#Plot current elements
title("Centres of elements of broken wire")
xlabel(r"$x\longrightarrow$")
ylabel(r"$y\longrightarrow$")
plot(rl[:,0],rl[:,1],'bo',markersize=3)
savefig("Current elements.jpg")
close()

#Plot the quiver for I vs r
title("Quiver plot for current")
xlabel(r"$x\longrightarrow$")
ylabel(r"$y\longrightarrow$")
quiver(rl[:,0],rl[:,1],I_plot[:,0],I_plot[:,1],color ='blue',scale=30)
grid()
savefig("Current Quiver Plot.jpg")
close()


#Make a meshgrid of 3,3,1000 with x,y centred at 0 whereas z ranging from 1 to 1000
xi = 3
yi = 3
zi = 1000
x=np.linspace(-pow(10,-5),pow(10,-5),xi)
y=np.linspace(-pow(10,-5),pow(10,-5),yi)
z=np.linspace(1,zi+1,zi)
#Indexing = 'ij' makes meshgrid of form x,y,z        
xx,yy,zz = np.meshgrid(x,y,z,indexing='ij')           
#General vector r for any position of size 3,3,1000      
r = np.stack((xx,yy,zz),axis=3)                                
#General vector r for any position, copied 100 times and reshaped
r_vectorization = np.tile(r,100).reshape(3,3,1000,100,3)

#Function calc(l) for given l to find R =  |r-rl[l]|
def calc_l(r,rl,l):                                         
    return norm(r-rl[l],axis=-1)

#Bz using for loop using calc(l) function only for Dynamic case
def Bz_by_for_loop():
    dAx = []
    dAy = []
    for l in range(N):
        R = calc_l(r,rl,l)
        dAx.append((mu/(4*np.pi))*I[l][0]*exp(-1j*k*R)*radius*dphi_abs/(R))
        dAy.append((mu/(4*np.pi))*I[l][1]*exp(-1j*k*R)*radius*dphi_abs/(R))
    
    #Ax,Ay given by sum of dAx,dAy respectively
    Ax = np.sum(dAx, axis=0).reshape(3,3,1000)
    Ay = np.sum(dAy, axis=0).reshape(3,3,1000)
    
    #Since Bz is curl of A, its given as follows
    Bz = (Ay[2,1,:]-Ay[0,1,:]-(Ax[1,2,:]-Ax[1,0,:]))/(2*pow(10,-5))       
    return Bz
    
#Function calc_all extended from calc(l) to find R = |r-rl| for all l     
def calc_all(r,rl):
    '''
    Extra dimension is added since this function stores the value of |r-rl| for
    all of the 100 rl's corresponding to the 3,3,1000 points in the meshgrid
    '''
    return norm(r-rl,axis=-1)    

#Bz using vectorization using calc_all extended from calc(l)    
def Bz_by_vectorization():
    R = calc_all(r_vectorization,rl)                                          

    #Taking sum of -1th elements in array since it contains dl element of sum
    Ax = np.sum((mu/(4*np.pi))*I[:,0]*exp(-1j*k*R)*radius*dphi_abs/(R),axis=-1)          
    Ay = np.sum((mu/(4*np.pi))*I[:,1]*exp(-1j*k*R)*radius*dphi_abs/(R),axis=-1)   

    #Finding Potential Fields for static case. No e^(-jkR) factor in this summation
    Ax_st = np.sum((mu/(4*np.pi))*I[:,0]*radius*dphi_abs/(R),axis=-1)         
    Ay_st = np.sum((mu/(4*np.pi))*I[:,1]*radius*dphi_abs/(R),axis=-1)   

    #Taking Bz as the curl = dAy/dx - dAx/dy
    Bz = (Ay[2,1,:]-Ay[0,1,:]-(Ax[1,2,:]-Ax[1,0,:]))/(2*pow(10,-5))      

    #Similar calculations of Bz for static current
    Bz_st = (Ay_st[2,1,:]-Ay_st[0,1,:]-(Ax_st[1,2,:]-Ax_st[1,0,:]))/(2*pow(10,-5))

    return Bz, Bz_st
    
Bz, Bz_static = Bz_by_vectorization()
Bz1 = Bz_by_for_loop()

#Plotting absolute value of Bz vs z in loglog scale for the for loop and vectorization method
fig,axs = subplots(2)
fig.suptitle("Bz using i) for loop ii) Vectorization")
axs[0].loglog()
axs[0].grid()
axs[0].set_xlabel("log(z)")
axs[0].set_ylabel("log(Bz)")
axs[0].plot(z,abs(Bz1))
axs[1].loglog()
axs[1].grid()
axs[1].set_xlabel("log(z)")
axs[1].set_ylabel("log(Bz)")
axs[1].plot(z,abs(Bz))
savefig("Different computations for B(z).jpg")
close()

#Finding least error fit and plotting fit vs actual curve in loglog scale
def least_error_fit(Bz,z):
    '''     
    Taking elements from z>10 onwards since it can be seen from graph of actual 
    curve that it becomes linear in log scale after certain value,
    This is done to get the best possible fit for the curve by considering as 
    much domain of z as possible. The c1 we get is actually log(c) if we 
    have fit to cz^b since we actually fit log(Bz) to b*log(z) + log(c), 
    so c = e^c1 so we return exp(c)         
    '''        
    y = log(abs(Bz))[10:]                       
    x = log(z[10:])                                 
    A = c_[x,ones(1000)[10:]]
    b,c = lstsq(A,y)[0]
    return b,exp(c)


#Plotting fits for Dynamic case with right side derivative only since absolute gives noise
b1,c1 = least_error_fit(Bz,z)          
print("c = "+str(c1)+" and b = "+str(b1)+" where Bz is fitted to cz^b for Dynamic case")

x = log(z[10:])
title("Actual Curve and Linear Fit vs z for Dynamic case in loglog")
xlabel("log(z)")
ylabel("log(|Magnetic Field|)")
loglog()
plot(z,abs(Bz), label="Actual Curve")
plot(exp(x),c1*exp(b1*x),label="Linear Fit")
grid()
legend()
savefig("Actual vs Fit for Dynamic.jpg")
close()        

#Plotting fits for Static case with right side derivative only since absolute gives noise
b2,c2 = least_error_fit(Bz_static,z)       
print("c = "+str(c2)+" and b = "+str(b2)+" where Bz is fitted to cz^b for Statics case\n")

title("Actual Curve and Linear Fit vs z for Static case in loglog")
xlabel("log(z)")
ylabel("log(|Magnetic Field|)")
loglog()
plot(z,abs(Bz_static), label="Actual Curve")
plot(exp(x),c2*exp(b2*x),label="Linear Fit")
grid()
legend()
savefig("Actual vs Fit for Static.jpg")
close()        

