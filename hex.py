import numpy as np
import numpy.linalg as alg
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size': 16})

def A(n,k,nPhi):return 2*np.exp(1j*(n*np.pi*nPhi + k/2))*np.cos(n*np.pi*nPhi + k/2)
def Astar(n,k,nPhi):return 2*np.exp(-1j*(n*np.pi*nPhi + k/2))*np.cos(n*np.pi*nPhi + k/2)

def hexH(nPhi,nu,Enot=0,HNum = 50):#HNum must be even (2qx2q eigenproblem)
    outArr = np.zeros((HNum,HNum),dtype = complex)
    for i in range(0,HNum,2):
        outArr[i+1][i] = 1
        outArr[i][i+1] = 1           #m = i/2+1
        outArr[i+1][(i+2)%HNum] = Astar(i/2+1,nu,nPhi)
        outArr[(i+2)%HNum][i+1] = A(i/2+1,nu,nPhi)
    outArr[0][-1] = A(1,nu,nPhi)*np.exp(-1j*nu*HNum/2)
    outArr[-1][0] = Astar(1,nu,nPhi)*np.exp(1j*nu*HNum/2)
    return outArr

#nuRes,pRes,latticeRes = 10,600,150
nuRes,pRes,latticeRes = 1,100,100

fig1=plt.figure()
ax1=fig1.gca()
fig1.set_size_inches(11,8)

ax1.set_xlabel(r"$\phi / \phi_0$",size=24)
ax1.set_ylabel(r"$\epsilon$",size=24)

ax1.tick_params(which='both',direction="in")

for i in range(nuRes):
    nu = 2*np.pi*i/nuRes
    for j in range(pRes+1):
        nPhi = j/pRes
        H = hexH(nPhi,nu,HNum = latticeRes)
        eVals = alg.eig(H)
        ax1.plot(nPhi*np.ones(latticeRes),eVals[0].real,'o',color = 'tab:blue',ms = 5/pRes)

fig1.savefig("Butterfly_Hex.png")#,transparent = True)
