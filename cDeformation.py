import numpy as np
import numpy.linalg as alg
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation

matplotlib.rcParams.update({'font.size': 16})


def A(j,nPhi, nu=0,t_a=1):return -2*np.cos(2*np.pi*nPhi*j + nu)
def B(j,nPhi,t_c,t_cP, nu=0):
    return 1 + t_c*np.exp(1j*(nPhi*np.pi*(2*j+1))+nu) + t_cP*np.exp(-1j*(nPhi*np.pi*(2*j+1))+nu)
def BStar(j,nPhi,t_c,t_cP, nu=0):
    return 1 + t_c*np.exp(-1j*(nPhi*np.pi*(2*j+1))+nu) + t_cP*np.exp(1j*(nPhi*np.pi*(2*j+1))+nu)

#Assuming k_y = k_x0 = 0, and t_a = t_b = 1
def Hamil(nPhi,t_c,t_cP, HNum=100, nu=0):
    outArr = np.zeros((HNum,HNum),dtype = complex)
    for i in range(0,HNum):
        outArr[i][i] = A(i,nPhi)
        outArr[(i+1)%HNum][i] = BStar(i,nPhi,t_c,t_cP)
        outArr[i][(i+1)%HNum] = B(i,nPhi,t_c,t_cP)
    #Can Remove periodicity
    #outArr[0][-1] = 0
    #outArr[-1][0] = 0
    return outArr
    
t_cInitial,t_cFinal = 0,1
t_cPInitial,t_cPFinal = 0,1
pRes,latticeRes = 800,100 
nFrames = 100

dtc = (t_cFinal-t_cInitial)/nFrames
dtcP = (t_cPFinal-t_cPInitial)/nFrames

def Update(frame):
    tc =  t_cInitial + frame*dtc
    tcP = t_cPInitial + frame*dtcP
    
    newX = np.empty(0)
    newY = np.empty(0)

    for j in range(pRes+1):
        nPhi = j/pRes
        H = Hamil(nPhi,tc,tcP, HNum = latticeRes)
        newX = np.append(newX, nPhi*np.ones(latticeRes))
        newY = np.append(newY, alg.eig(H)[0].real)
    
    ax.clear()
    line = ax.plot(newX,newY,'o',color = 'tab:blue',ms = 150/pRes)
    return line   
    
##############################################################################
fig=plt.figure()
ax=fig.gca()
fig.set_size_inches(11,8)

ax.set_xlabel(r"$\phi / \phi_0$",size=24)
ax.set_ylabel(r"$\epsilon = \omega^2$",size=24)

anim = animation.FuncAnimation(fig, Update,frames=nFrames, interval=10, blit=True)

anim.save('Deform2.mp4', fps=10, extra_args=['-vcodec', 'libx264'])

plt.show()