import numpy as np
import numpy.linalg as alg
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size': 16})

def squareH(nPhi,nu,Enot=0,HNum = 50):
    outArr = np.zeros((HNum,HNum))
    for m in range(HNum):
        outArr[m][m] = 2*np.cos(2*np.pi*(m+1)*nPhi - nu)
        outArr[m][(m+1)%HNum] = 1
        outArr[m][(m-1)%HNum] = 1
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
        H = squareH(nPhi,nu,HNum = latticeRes)
        eVals = alg.eigh(H)
        ax1.plot(nPhi*np.ones(latticeRes),eVals[0],'o',color = 'tab:blue',ms = 5/pRes)

fig1.savefig("Butterfly.png")#,transparent = True)
	
