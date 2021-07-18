'''
                            EE2703 Applied Programming Lab Week 9
Author :  Pranav Phatak
Roll No : EE19B105
Assignment 9 : DFT of Aperiodic Signals
Inputs: None
Outputs: Several graphs of DFT of various functions
'''


import numpy as np
from matplotlib import cm
from pylab import *
from mpl_toolkits.mplot3d import Axes3D


def log10(x):
    return np.log(x)/2.303
                                                            ### Examples ###

# Plotting Spectrum of sin(sqrt(2)t) :

t=linspace(-np.pi,np.pi,65)[:-1]
dt=t[1]-t[0];fmax=1/dt
y=np.sin(np.sqrt(2)*t)
y[0]=0                                      
y=fftshift(y)                               
Y=fftshift(fft(y))/64.0
w=linspace(-np.pi*fmax,np.pi*fmax,65)[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y))
xlim([-10,10])
ylabel(r"$|Y|\longrightarrow$")
title(r"Spectrum of $sin(\sqrt{2}t)$")
grid(True)
subplot(2,1,2)
plot(w,np.angle(Y),'ro')
xlim([-10,10])
ylabel(r"Phase of $Y\longrightarrow$")
xlabel(r"$\omega\longrightarrow$")
grid(True)
savefig("Spectrum of sin(sqrt(2)t).jpg")


# Plotting pure function f(t) = sin(sqrt(2)t)

t1=linspace(-pi,pi,65)[:-1]
t2=linspace(-3*pi,-pi,65)[:-1]
t3=linspace(pi,3*pi,65)[:-1]
figure(2)
plot(t1,np.sin(np.sqrt(2)*t1),'b')
plot(t2,np.sin(np.sqrt(2)*t2),'r')
plot(t3,np.sin(np.sqrt(2)*t3),'r')
ylabel(r"$y\longrightarrow$")
xlabel(r"$t\longrightarrow$")
title(r"$sin(\sqrt{2}t)$")
grid(True)
savefig("Pure function of sin(sqrt(2)t).jpg")


# Plotting sin(sqrt(2)t) periodically extended from [-pi,pi] to [-inf, inf] 


t1=linspace(-pi,pi,65)[:-1]
t2=linspace(-3*pi,-pi,65)[:-1]
t3=linspace(pi,3*pi,65)[:-1]
y=np.sin(np.sqrt(2)*t1)
figure(3)
plot(t1,y,'bo')
plot(t2,y,'ro')
plot(t3,y,'ro')
ylabel(r"$y\longrightarrow$")
xlabel(r"$t\longrightarrow$")
title(r"$sin(\sqrt{2}t)$ with $t$ wrapping every $2\pi$ ")
grid(True)
savefig("Periodically extended sin(sqrt(2)t).jpg")


# Plotting Spectrum of periodically extended Ramp function  


t=linspace(-np.pi,np.pi,65)[:-1]
dt=t[1]-t[0]
fmax=1/dt
y=t
y[0]=0                                  
y=fftshift(y)                           
Y=fftshift(fft(y))/64.0
w=linspace(-np.pi*fmax,np.pi*fmax,65)[:-1]
figure(4)
semilogx(abs(w),20*log10(abs(Y)))
xlim([1,10])
ylim([-20,0])
xticks([1,2,5,10],["1","2","5","10"])
ylabel(r"$|Y|(dB)\longrightarrow$")
title(r"Spectrum of a digital ramp")
xlabel(r"$\omega\longrightarrow$")
grid(True)
savefig("Spectrum of ramp.jpg")


# Plot of sin(sqrt(2)t)*w(t) periodically extended from [-np.pi,np.pi] to [-inf,inf]


t1=linspace(-pi,pi,65)[:-1]
t2=linspace(-3*pi,-pi,65)[:-1]
t3=linspace(pi,3*pi,65)[:-1]
n=np.arange(64)
wnd=fftshift(0.54+0.46*np.cos(2*pi*n/63))
y=np.sin(np.sqrt(2)*t1)*wnd
figure(5)
plot(t1,y,'bo')
plot(t2,y,'ro')
plot(t3,y,'ro')
ylabel(r"$y\longrightarrow$")
xlabel(r"$t\longrightarrow$")
title(r"$sin(\sqrt{2}t)\times w(t)$ with $t$ wrapping every $2\pi$ ")
grid(True)
savefig("Periodically extended sin(sqrt(2)t)w(t).jpg")


# Spectrum of sin(sqrt(2)t)*w(t) with period 2pi


t=linspace(-np.pi,np.pi,65)[:-1]
dt=t[1]-t[0];
fmax=1/dt
n=np.arange(64)
wnd=fftshift(0.54+0.46*np.cos(2*np.pi*n/63))
y=np.sin(np.sqrt(2)*t)*wnd
y[0]=0 
y=fftshift(y) 
Y=fftshift(fft(y))/64.0
w=linspace(-np.pi*fmax,np.pi*fmax,65)[:-1]
figure(6)
subplot(2,1,1)
plot(w,abs(Y))
xlim([-8,8])
ylabel(r"$|Y|\longrightarrow$")
title(r"Spectrum of $sin(\sqrt{2}t)\times w(t)$ with period $2\pi$")
grid(True)
subplot(2,1,2)
plot(w,np.angle(Y),'ro')
xlim([-8,8])
ylabel(r"Phase of $Y\longrightarrow$")
xlabel(r"$\omega\longrightarrow$")
grid(True)
savefig("Spectrum of sin(sqrt(2)t)w(t) with period 2pi.jpg")


# Spectrum of sin(sqrt(2)t)*w(t) with period 8pi


t=linspace(-4*np.pi,4*np.pi,257)[:-1]
dt=t[1]-t[0];
fmax=1/dt
n=np.arange(256)
wnd=fftshift(0.54+0.46*np.cos(2*np.pi*n/256))
y=np.sin(np.sqrt(2)*t)
y=y*wnd
y[0]=0 
y=fftshift(y) 
Y=fftshift(fft(y))/256.0
w=linspace(-np.pi*fmax,np.pi*fmax,257)[:-1]
figure()
subplot(2,1,1)
plot(w,abs(Y))
xlim([-8,8])
ylabel(r"$|Y|\longrightarrow$")
title(r"Spectrum of $sin(\sqrt{2}t)\times w(t)$ with period $8\pi$")
grid(True)
subplot(2,1,2)
plot(w,np.angle(Y),'ro',lw=2)
xlim([-8,8])
ylabel(r"Phase of $Y\longrightarrow$")
xlabel(r"$\omega\longrightarrow$")
grid(True)
savefig("Spectrum of sin(sqrt(2)t)w(t) with period 8pi.jpg")


                                                    ### Assignment Questions ###


# General function for finding spectrum and plotting it :
def FFT_and_plotting_spectrum(lim,n,f,t_=0,t_lims = False, windowing= False,xlim1=10,title1 = r"Spectrum of $sin(\sqrt{2}t)$",xlabel1 = r"$\omega\longrightarrow$",ylabel1= r"$|Y|\longrightarrow$", ylabel2 = r"Phase of $Y\longrightarrow$", display = True,savename=None):
    if(t_lims):
        t = t_
    else:
        t=linspace(-lim,lim,n+1)[:-1]
    dt=t[1]-t[0];
    fmax=1/dt
    y = f(t)
    if (windowing):
        m=np.arange(n)
        wnd=fftshift(0.54+0.46*np.cos(2*np.pi*m/n))
        y = y*wnd
    y[0]=0 
    y=fftshift(y) 
    Y=fftshift(fft(y))/float(n)
    w=linspace(-np.pi*fmax,np.pi*fmax,n+1)[:-1]
    
    mag = abs(Y)
    ph = np.angle(Y)
    
    if (display):

        figure()
        subplot(2,1,1)
        plot(w,mag)
        xlim([-xlim1,xlim1])
        ylabel(ylabel1)
        title(title1)
        grid(True)
        subplot(2,1,2)
        ph[np.where(mag<3e-3)] = 0
        plot(w,ph,'ro')
        xlim([-xlim1,xlim1])
        ylabel(ylabel2)
        xlabel(xlabel1)
        grid(True)
        savefig(savename)
    
    return w,Y


# Function f(t) = cos^3(w0*t) :
def cos3(t,w0=0.86):
    return (np.cos(w0*t))**3


# Function f(t) = cos(w0*t + delta) :
def cosine(t,w0=1.5,delta=0.5):
    return np.cos(w0*t + delta)

# Question 2 
a,b = FFT_and_plotting_spectrum(4*np.pi,64*4,cos3,xlim1= 3,windowing=False, title1 = r"Spectrum of $cos^3(w_0t)$  without windowing", display = True,savename='Spectrum of cos^3(wt) without windowing.jpg')
a,b = FFT_and_plotting_spectrum(4*np.pi,64*4,cos3,xlim1= 3,windowing=True, title1 = r"Spectrum of $cos^3(w_0t)$  with windowing",display = True,savename='Spectrum of cos^3(wt) with windowing.jpg')

# Question 3
#FFT_and_plotting_spectrum of cos(wt+delta) windowed windowing to estimate w, delta
w,Y = FFT_and_plotting_spectrum(np.pi,128,cosine,xlim1= 3,windowing=True, title1 = r"Spectrum of $cos(w_0t + \delta)$  without noise",display = True,savename='cos(wt+d) without noise.jpg')
print("Without noise : ")

def estimate_omega_delta(w,Y):
    # by finding the weighted average of all w>0, we find w0. Delta is found by calculating the phase at w closest to w0.
    ii = where(w>=0)
    w_from_spectrum = sum(abs(Y[ii])**2*w[ii])/sum(abs(Y[ii])**2)
    i = abs(w-w_from_spectrum).argmin()
    delta = angle(Y[i])
    print("Value of w0 without noise from the spectrum: ",w_from_spectrum)
    print("Value of delta without noise from the spectrum: ",delta)

estimate_omega_delta(w,Y)


# Question 4 : Added Gaussian Noise
print("With added Gaussian Noise : ")
def cosine_with_noise(t,w0=1.5,delta=0.5):
    return np.cos(w0*t + delta) + 0.1*np.random.randn(128)

w,Y = FFT_and_plotting_spectrum(np.pi,128,cosine_with_noise,xlim1= 3,windowing=True, title1 = r"Spectrum of $cos(w_0t + \delta)$  with noise",display = True,savename='cos(wt+d) with noise.jpg')

estimate_omega_delta(w,Y)

# Question 5
def chirp(t):
    return np.cos(16*(1.5 + t/(2*np.pi))*t) 

w,Y = FFT_and_plotting_spectrum(np.pi,1024,chirp,xlim1= 60,windowing=True, title1 = r"Spectrum of chirp function with windowing",display = True,savename='Chirp without windowing.jpg')
w,Y = FFT_and_plotting_spectrum(np.pi,1024,chirp,xlim1= 60,windowing=False, title1 = r"Spectrum of chirp function wihout windowing",display = True,savename='Chirp with windowing.jpg')

#  Question 6
t=np.linspace(-np.pi,np.pi,1025);t=t[:-1]
t_arrays=np.split(t,16)

Y1_mags=np.zeros((16,64))
Y1_angles=np.zeros((16,64))
Y2_mags=np.zeros((16,64))
Y2_angles=np.zeros((16,64))

for i in range(len(t_arrays)):
    w1,Y1 = FFT_and_plotting_spectrum(lim = 10,t_ = t_arrays[i],t_lims=True,n = 64,f = chirp,xlim1= 60,windowing=False, title1 = r"Spectrum of chirp function",display = False)
    Y1_mags[i] =  abs(Y1)
    Y1_angles[i] = np.angle(Y1)

for i in range(len(t_arrays)):
    w2,Y2 = FFT_and_plotting_spectrum(lim = 10,t_ = t_arrays[i],t_lims=True,n = 64,f = chirp,xlim1= 60,windowing=True, title1 = r"Spectrum of chirp function",display = False)
    Y2_mags[i] =  abs(Y2)
    Y2_angles[i] = np.angle(Y2)


# 3D Plot of frequency, magnitude and time for non-windowed Chirp function :
fig = figure()
ax = fig.add_subplot(111, projection='3d')
ax.view_init(45,15)
t=np.linspace(-np.pi,np.pi,1025);t=t[:-1]
fmax = 1/(t[1]-t[0])
t=t[::64]
w1=np.linspace(-fmax*np.pi,fmax*np.pi,64+1);w1=w1[:-1]
t,w1=np.meshgrid(t,w1)

surf=ax.plot_surface(t,w1,Y1_mags.T,cmap=cm.BrBG,linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.5, aspect=5)
ylabel("Frequency$\longrightarrow$")
xlabel("Time$\longrightarrow$")
title("3D Magnitude Plot for Non-Windowed Chirp Function")
savefig("Magnitude 3D plot non-windowed chirp.jpg")

fig = figure()
ax = fig.add_subplot(111, projection='3d')
ax.view_init(45,15)
surf=ax.plot_surface(t,w1,Y1_angles.T,cmap=cm.BrBG,linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.5, aspect=5)
ylabel("Frequency$\longrightarrow$")
xlabel("Time$\longrightarrow$")
title("3D Phase Plot for Non-Windowed Chirp Function")
savefig("Phase 3D plot non-windowed chirp.jpg")


# 3D Plot of frequency, magnitude and time for Windowed Chirp function :
fig = figure()
ax = fig.add_subplot(111, projection='3d')
ax.view_init(45,15)
t=np.linspace(-np.pi,np.pi,1025);t=t[:-1]
fmax = 1/(t[1]-t[0])
t=t[::64]
w2=np.linspace(-fmax*np.pi,fmax*np.pi,64+1);w2=w2[:-1]
t,w2=np.meshgrid(t,w2)

surf=ax.plot_surface(t,w2,Y2_mags.T,cmap=cm.BrBG,linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.5, aspect=5)
ylabel("Frequency$\longrightarrow$")
xlabel("Time$\longrightarrow$")
title("3D Magnitude Plot for Windowed Chirp Function")
savefig("Magnitude 3D plot windowed chirp.jpg")

fig = figure()
ax = fig.add_subplot(111, projection='3d')
ax.view_init(45,15)
surf=ax.plot_surface(t,w2,Y2_angles.T,cmap=cm.BrBG,linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.5, aspect=5)
ylabel("Frequency$\longrightarrow$")
xlabel("Time$\longrightarrow$")
title("3D Phase Plot for Windowed Chirp Function")
savefig("Phase 3D plot windowed chirp.jpg")
