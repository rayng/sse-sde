
# Not for beta doubling.

import numpy as np
import matplotlib.pyplot as plt


L=8
nmax=10
nobs=13
sys=str(L)+"x"+str(L)
beta=L
N=L*L
hh=np.arange(1,3.25,0.25)
samples = 1000
print hh


for h in np.arange(0,8):
    fn=np.loadtxt("../data/"+sys+"/DIST/"+sys+"h"+str(h)+"p0beta0")  #unpack sorts it according to columns
    
    rho_col=1+(0)*nobs+7
    rhosq_col=1+(0)*nobs+8
    m_col =1+(0)*nobs+3
    msq_col = 1+(0)*nobs+4
    plt.subplot(1,2,1)
    plt.title('Stiffness')
    #plt.scatter([hh[h]], [fn[rho_col]], c='r')#
    rhoerr=1/np.sqrt(samples) * ( pow(np.pi/beta,2)*fn[rhosq_col] - fn[rho_col]*fn[rho_col] ) 
    plt.errorbar([hh[h]], [fn[rho_col]], yerr=rhoerr, c='r', marker='o')
    
    
    plt.subplot(1,2,2)
    comp = 1.*beta*N*(fn[msq_col] - fn[m_col]*fn[m_col] )
    print h, comp
    plt.title('Compressibility')
    plt.scatter([hh[h]], [comp], c='b')

    
        
    
    #plt.show()
plt.show()
    
    
    
