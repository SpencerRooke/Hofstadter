import numpy as np
import numpy.linalg as alg
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation

matplotlib.rcParams.update({'font.size': 16})

#division by 1+Lambda for normalization (keeps us within [-4,4]) 
def K(n,Lambda,nPhi,nu=0):return (1+Lambda*np.cos(2*np.pi*nPhi*n+nu)) / (1+Lambda) 
def V(n,Lambda,nPhi,nu=0):return K(n+1,Lambda,nPhi,nu) +  K(n-1,Lambda,nPhi,nu)

def qChain(nPhi,Lambda = 1,HNum = 100,nu=0):
    outArr = np.zeros((HNum,HNum))
    #This can probably be done better
    for i in range(0,HNum):
        outArr[i][i] = V(i,Lambda,nPhi,nu)
        outArr[(i+1)%HNum][i] =  outArr[i][(i+1)%HNum] = -K(i,Lambda,nPhi,nu)
    #Makes this an Hnum length chain. Can remove for periodicity(Circular Chain)
    #outArr[0][-1] = 0
    #outArr[-1][0] = 0
    return outArr

pRes,latticeRes = 800,100 

def Update(frame):
    Lambda = frame**1.25/40 #nonlinear, as the later bit of the animation is slow otherwise
    
    newX = np.empty(0)
    newY = np.empty(0)

    for j in range(pRes+1):
        nPhi = j/pRes
        H = qChain(nPhi,Lambda,HNum = latticeRes,nu=0)
        newX = np.append(newX, nPhi*np.ones(latticeRes))
        newY = np.append(newY, alg.eigh(H)[0])
    
    ax.clear()
    line = ax.plot(newX,newY,'o',color = 'tab:blue',ms = 150/pRes)
    return line   
    
##############################################################################
fig=plt.figure()
ax=fig.gca()
fig.set_size_inches(11,8)

ax.set_xlabel(r"$\phi / \phi_0$",size=24)
ax.set_ylabel(r"$\epsilon = \omega^2$",size=24)

anim = animation.FuncAnimation(fig, Update,frames=176, interval=10, blit=True)

anim.save('qChain.mp4', fps=15, extra_args=['-vcodec', 'libx264'])

plt.show()
