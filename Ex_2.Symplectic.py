import numpy as np
import matplotlib.pyplot as plt

# Time parameters
deltat = 10**(-2)
N = int(10/deltat)

#Initialization and initial conditions
X_1, X_2, P_1, P_2, error_1, error_2, H_1, H_2, Hs_1,Hs_2 = np.empty((10,N))
X_1[0] = X_2[0] = 1
P_1[0] = P_2[0] = error_1[0] = error_2[0] = 0
H_1[0] = H_2[0] = 0.5*(X_1[0]**2+ P_1[0]**2)
Hs_1[0] = Hs_2[0] = 0.5*(X_1[0]**2+ P_1[0]**2-X_1[0]*P_1[0])

#Real solutions
X = np.cos(deltat * np.arange(N))
P = -np.sin(deltat * np.arange(N))

#Algorithms
for i in range(N-1):
    P_1[i+1]=P_1[i]-X_1[i]*deltat
    X_1[i+1]=X_1[i]+P_1[i+1]*deltat
    X_2[i+1]=X_2[i]+P_2[i]*deltat
    P_2[i+1]=P_2[i]-X_2[i]*deltat

    #Computing the error
    error_1[i+1]=(X_1[i+1]-X[i+1])**2
    error_2[i+1]=(X_2[i+1]-X[i+1])**2

    #Computing H
    H_1[i+1]=0.5*(X_1[i+1]**2+P_1[i+1]**2)
    H_2[i+1]=0.5*(X_2[i+1]**2+P_2[i+1]**2)
    Hs_1[i+1]=0.5*(X_1[i+1]**2+P_1[i+1]**2-X_1[i+1]*P_1[i+1]*deltat)
    Hs_2[i+1]=0.5*(X_2[i+1]**2+P_2[i+1]**2-X_2[i+1]*P_2[i+1]*deltat)


#Error with each method
Symplectic_rms = np.sqrt((1/len(X))*np.sum(error_1))
Non_symplectic_rms = np.sqrt((1/len(X))*np.sum(error_2))
print(f'The RMS of the symplectic algorithm is {Symplectic_rms}')
print(f'The RMS of the non-symplectic algorithm is {Non_symplectic_rms}')
print(f'The symplectic algorithm error is a {Symplectic_rms*100/(Symplectic_rms+Non_symplectic_rms)} % of the non-symplectic error')

#Plot
time = np.linspace(0,10,N)

plt.plot(time,X_1, label='Symplectic')
plt.plot(time,X_2, label='Non-Symplectic')
plt.plot(time, X, label ='Real solution')
plt.legend()
plt.title('Solution '+r'$x(t)$')
plt.xlabel(r'$t$')
plt.ylabel(r'x(t)')
plt.show()

plt.plot(time, H_1, label='Symplectic')
plt.plot(time,H_2, label='Non-Symplectic')
plt.legend()
plt.title('Hamiltonian as a function of time')
plt.xlabel(r'$t$')
plt.ylabel(r'$H(t)$')
plt.show()

plt.plot(time, Hs_1, label='Symplectic')
plt.plot(time,Hs_2, label='Non-Symplectic')
plt.legend()
plt.title('Shadow Hamiltonian as a function of time')
plt.xlabel(r'$t$')
plt.ylabel(r'$H\'(t)$')
plt.show()

