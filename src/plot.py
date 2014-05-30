
# Not for beta doubling.

import numpy as np
import matplotlib.pyplot as plt


L=8
nmax=5
nobs=13
sys=str(L)+"x"+str(L)
beta=L
N=L*L
hh=np.arange(1,3.25,0.25)
print hh

for h in np.arange(0,8):
    fn=np.loadtxt("../data/"+sys+"/DIST/"+sys+"h"+str(h)+"p0beta0")  #unpack sorts it according to columns
    temp=[]
    for n in np.arange(0,nmax):
        temp.append(pow(2.,n))
        rho_col=1+(n)*nobs+7
        m_col =1+(n)*nobs+3
        msq_col = 1+(n)*nobs+4
        plt.subplot(1,2,1)
        plt.title('Stiffness')
        plt.xscale('log')
        plt.scatter([temp[n]], [fn[rho_col]], c='r')#
        plt.subplot(1,2,2)
        comp = 1.*temp[n]*N*(fn[msq_col] - fn[m_col]*fn[m_col] )
        plt.title('Compressibility')
        plt.xscale('log')
        plt.scatter([temp[n]], [comp], c='b')
        
        
    plt.show()
    plt.clf()
    #plt.show()

    
    
    
