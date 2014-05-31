
# Not for beta doubling.

import numpy as np
import matplotlib.pyplot as plt



L=16
nmax=10
nobs=13
sys=str(L)+"x"+str(L)
beta=L
N=L*L
hh=np.arange(1,3.25,0.25)
samples = 1000
nlags=100
print hh



def autocorr(x, nlags):
    N=len(x)
    xt=np.zeros(nlags)
    for lag in np.arange(0,nlags):  
        for i in np.arange(0,N-lag):
            xt[lag] += x[i]*x[i+lag]
        xt[lag] /= (N-lag)

    xsq_av = np.mean(x)*np.mean(x)
    
    return (xt-xsq_av)/ np.var(x)
        


nlags=100

for h in np.arange(0,8):
    fn=np.loadtxt("../data/"+sys+"/DIST/"+sys+"h"+str(h)+"p0beta0autocorr",unpack=True)  #unpack sorts it according to columns

    rhosq_col=(0)*nobs+8
    m_col =(0)*nobs+3
    msq_col = (0)*nobs+4
    plt.subplot(1,2,1)
    plt.title('energy autocorr')
    energ_col=1
    energy=fn[energ_col]
    plt.plot(np.arange(0,nlags), autocorr(energy,nlags), 'bo-')

    plt.subplot(1,2,2)
    plt.title('stiffness autocorr')
    rho_col=(0)*nobs+7
    stiff = fn[rho_col]
    plt.plot(np.arange(0,nlags), autocorr(stiff,nlags), 'bo-')

    plt.savefig(sys+'h'+str(h)+'autocorrNe16.png')
    #plt.show()
    plt.clf()
    
    
    
    
    
"""
    plt.errorbar([hh[h]], [fn[rho_col]], yerr=rhoerr, c='r', marker='o')
    plt.subplot(1,2,2)
    comp = 1.*beta*N*(fn[msq_col] - fn[m_col]*fn[m_col] )
    print h, comp
    plt.title('Compressibility')
    plt.scatter([hh[h]], [comp], c='b')
"""
    
    #plt.show()
#plt.show()
    
    
    
