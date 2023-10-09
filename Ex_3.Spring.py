import numpy as np
import matplotlib.pyplot as plt

#Time parameters
N = 50000
deltat = 10**(-4)

#Important cosntants
g = 10
mindist = deltat*10

#Variable initialization and initial conditions
def spring_simulation(n):
    Vel, Pos = np.empty((2,n,N))
    Ekin,Egrav,Eel, Et = np.empty((4,N))
    for k in range(n):
        Vel[:,0] = 0
        Pos[k,0] = k+4
        Ekin[0] = np.sum(0.5*Vel[k,0]**2)
        Egrav[0] = np.sum(Pos[k,0]*g)
        Eel[0] = np.sum(0.5*np.diff(Pos[:,0])**2)
        Et[0] = Ekin[0]+Egrav[0]+Eel[0]

    #Acceleration at time t
    def atop(y_bot,y_top):
        return(-(y_top-y_bot)-1)
    def abot(y_bot,y_top):
        return(-(y_bot-y_top)-1)

    #Verlet algorithm function
    def Verlet(pos,vel):
        Pos_new, Vel_half, A_new, Vel_new = np.empty((4,n))
        for k in range(n):
            if(k!=0 and k!=n-1):
                #Step 1
                Pos_new[k]=pos[k]+vel[k]*deltat+0.5*(abot(pos[k],pos[k+1])+atop(pos[k],pos[k-1])-g)*deltat**2

                #Step 2
                Vel_half[k]=vel[k]+0.5*(abot(pos[k],pos[k+1])+atop(pos[k],pos[k-1])-g)*deltat

                #Step 3
                A_new[k]=abot(Pos_new[k],Pos_new[k+1])+atop(Pos_new[k],Pos_new[k-1])-g

                #Step 4
                Vel_new[k]=Vel_half[k]+0.5*A_new[k]*deltat
            if(k==0):
                #Step 1
                Pos_new[k]=pos[k]+vel[k]*deltat+0.5*(abot(pos[k],pos[k+1])-g)*deltat**2

                #Step 2
                Vel_half[k]=vel[k]+0.5*(abot(pos[k],pos[k+1])-g)*deltat

                #Step 3
                A_new[k]=abot(Pos_new[k],Pos_new[k+1])-g

                #Step 4
                Vel_new[k]=Vel_half[k]+0.5*A_new[k]*deltat
            if(k==n-1):
                #Step 1
                Pos_new[k]=pos[k]+vel[k]*deltat+0.5*(atop(pos[k],pos[k-1])-g)*deltat**2

                #Step 2
                Vel_half[k]=vel[k]+0.5*(atop(pos[k],pos[k-1])-g)*deltat

                #Step 3
                A_new[k]=atop(Pos_new[k],Pos_new[k-1])-g

                #Step 4
                Vel_new[k]=Vel_half[k]+0.5*A_new[k]*deltat

        #Energy 
        Egrav = g*(np.sum(Pos_new)) 
        Eel = np.sum(0.5*np.diff(Pos_new)**2)
        Ekin = np.sum(0.5*Vel_new**2)
        Et = Egrav+ Eel + Ekin

        return(Pos_new,Vel_new,Egrav,Eel,Ekin,Et)

    #Algorithm Implementation
    for i in range(N-1):
        if(abs(np.any(np.diff(Pos[:,i])))>=mindist and (np.any(Pos[:,i]))>mindist):
            Pos[:,i+1],Vel[:,i+1],Egrav[i+1],Eel[i+1],Ekin[i+1],Et[i+1] = Verlet(Pos[:,i],Vel[:,i])
        if(np.any(Pos[:,i])>0):
            indices = np.where(Pos[:,i]<0)[0]
            Vel[indices,i] = -Vel[indices,i]
            Pos[:,i+1],Vel[:,i+1],Egrav[i+1],Eel[i+1],Ekin[i+1],Et[i+1] = Verlet(Pos[:,i],Vel[:,i])
        if(np.any(np.diff(Pos[:,i]))>=mindist):
            indices = np.where(np.diff(Pos[:,i]) < mindist)[0]
            Vel[indices,i] = -Vel[indices,i]
            Vel[indices+1,i] = -Vel[indices+1,i]
            Pos[:,i+1],Vel[:,i+1],Egrav[i+1],Eel[i+1],Ekin[i+1],Et[i+1] = Verlet(Pos[:,i],Vel[:,i])




    #Plot
    time=np.linspace(0,N*deltat,N)
    for k in range(n):
        plt.plot(time, Pos[k,:])
    plt.show()
    plt.plot(time, Et, label = 'Total Energy')
    plt.plot(time, Ekin, label = 'Kinetic Energy')
    plt.plot(time,Egrav, label = 'Gravitational Energy')
    plt.plot(time,Eel, label='Elastic Energy')
    plt.legend()
    plt.show()

spring_simulation(20)



