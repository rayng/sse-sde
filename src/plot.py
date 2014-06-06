
# This script plots stiffness and compressibility versus field value
# It can be easily modified to include other observables as well.
# Not for beta doubling.

import numpy as np
import matplotlib.pyplot as plt

L=8
nmax=10
nobs=13
sys=str(L)+"x"+str(L)
beta=2.*L
N=L*L
hbi=2.75
hbf=3.15
nsteps=8
hh=np.arange(hbi,hbf, (hbf-hbi)/nsteps )
samples = 1000
print hh


for h in np.arange(0,8):
    fn=np.loadtxt("../data/"+sys+"/DIST/TIM-"+sys+"h"+str(h)+"p0beta0")  #unpack sorts it according to columns
    
    rho_col=(0)*nobs+7
    rhosq_col=(0)*nobs+8
    m_col =(0)*nobs+3
    msq_col = (0)*nobs+4
    #plt.subplot(4,1,1)
    #plt.title('Stiffness')
    #plt.scatter([hh[h]], [fn[rho_col]], c='r')#
    #rhoerr=1/np.sqrt(samples) * ( pow(np.pi/beta,2)*fn[rhosq_col] - fn[rho_col]*fn[rho_col] ) 
    #plt.errorbar([hh[h]], [fn[rho_col]], yerr=rhoerr, c='r', marker='o')
    
    
    #plt.subplot(4,1,2)
    #comp = 1.*beta*N*(fn[msq_col] - fn[m_col]*fn[m_col] )
    #plt.title('Compressibility')
    #plt.scatter([hh[h]], [comp], c='b')
    
    
    #plt.subplot(4,1,3)
    #plt.title('magnetization')
    #plt.scatter([hh[h]], [fn[m_col]], c='g')


    #plt.subplot(4,1,4)
    plt.title('Energy susceptibility')
    Esus= (fn[12] - fn[11]*fn[11] - fn[11])/(beta*hh[h]*hh[h]*N)
    plt.scatter([hh[h]], [Esus], c='g')
    
    
    
    
    #plt.show()
plt.show()
    
    
    
