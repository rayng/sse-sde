
# Not for beta doubling.

import numpy as np
import matplotlib.pyplot as plt


L=4
nmax=6
nobs=13
sys=str(L)+"x"+str(L)
beta=L
N=L*L
hh=np.arange(1,3.25,0.25)
print hh

for h in np.arange(0,8):
    fn=np.loadtxt("../data/"+sys+"/DIST/"+sys+"h"+str(h)+"p0beta0")  #unpack sorts it according to columns

    for n in np.arange(0,nmax):
    
        rho_col=1+(n)*nobs+7
        m_col =1+(n)*nobs+3
        msq_col = 1+(n)*nobs+4
        plt.subplot(1,2,1)
        plt.title('Stiffness')
        plt.scatter([hh[h]], [fn[rho_col]], c='r')#
        plt.subplot(1,2,2)
        comp = 1.*beta*N*(fn[msq_col] - fn[m_col]*fn[m_col] )
        plt.title('Compressibility')
        plt.scatter([hh[h]], [comp], c='b')
        
        
    plt.show()
    
    #plt.show()
plt.show()
    
    
    
