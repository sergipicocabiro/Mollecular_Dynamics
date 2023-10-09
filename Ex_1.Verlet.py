import numpy as np
import matplotlib.pyplot as plt

#Important constants. 
    ### We will be using Astronomic Units for distances, Solar masses for the mass and years for time. All of this changes the value of the G constant
UA = 149597870700
M = 1.98847*10**(30)
any = 3600*24*365
G_old = 6.67428*10**(-11)
G = G_old*(any**2)*M/(UA**3) 

#Time parameters
N = 100000
deltat = 10**(-2)

#Initial parameters and variable initialization.
x0_metres = 5.28*10**(12)
y0_metres = 0
vx0_metressec = 0
vy0_metressec = 9.13*10**(2)
    ###We also need to change units:
x0 = x0_metres/UA
y0 = y0_metres/UA
vx0 = vx0_metressec*any/UA
vy0 = vy0_metressec*any/UA
    ###That is the way of initializing an array
X, Y, Vx, Vy, E, R = np.empty((6, N))
Vx2, Vy2, Ax, Ay = np.empty((4, N-1))

X[0]=x0
Y[0]=y0
Vx[0]=vx0
Vy[0]=vy0
E[0]=-G/(X[0]**2+Y[0]**2)**(1/2)+0.5*(Vx[0]**2+Vy[0]**2)
R[0]=np.sqrt(X[0]**2+Y[0]**2)


#Acceleration as function at time t
def ax(X,Y):
    return(-G*X/(X**2+Y**2)**(3/2))
def ay(X,Y):
    return(-G*Y/(X**2+Y**2)**(3/2))

#Verlet Algorithm
for i in range(N-1):
    #Step 1
    X[i+1]=X[i]+Vx[i]*deltat+0.5*ax(X[i],Y[i])*deltat**2
    Y[i+1]=Y[i]+Vy[i]*deltat+0.5*ay(X[i],Y[i])*deltat**2

    #Step 2
    Vx2[i]=Vx[i]+0.5*ax(X[i],Y[i])*deltat
    Vy2[i]=Vy[i]+0.5*ay(X[i],Y[i])*deltat

    #Step 3
    Ax[i]=ax(X[i+1],Y[i+1]) 
    Ay[i]=ay(X[i+1],Y[i+1])

    #Step 4
    Vx[i+1]=Vx2[i]+0.5*Ax[i]*deltat
    Vy[i+1]=Vy2[i]+0.5*Ay[i]*deltat

    #Energy calculation
    E[i+1]=-G/(X[i+1]**2+Y[i+1]**2)**(1/2)+0.5*(Vx[i+1]**2+Vy[i+1]**2)

    #Distance to the sun calculation
    R[i+1]=np.sqrt(X[i+1]**2+Y[i+1]**2)

plt.plot(X,Y)
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.title('Trajectory in the xy plane of the comet')
plt.show()
plt.plot(np.linspace(0,N*deltat,N),E)
plt.xlabel(r'$t$')
plt.ylabel(r'$E(t)=-\frac{GM}{r}+\frac{1}{2}v^{2}$')
plt.title('Energy as a function of time')
    ###There is a little error when computing the energy close to the origin. Probably a computation error of a very small quantity. Do not know yet how to fix it.
plt.show()
plt.plot(np.linspace(0,N*deltat,N),R)
plt.xlabel(r'$t$')
plt.ylabel(r'r(t)=\sqrt{x^2+y^2}')
plt.title('Distance to the sun as a function of time')
plt.show()
print(f'The minimal distance to the sun is {np.min(R)} AU')



