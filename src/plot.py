
# This script plots stiffness and compressibility versus field value
# It can be easily modified to include other observables as well.
# Not for beta doubling.

import numpy as np
import matplotlib.pyplot as plt



for L in (8,8):
    if(L==8):
        col='r'
    else:
        col='b'
    nmax=10
    nobs=13
    sys=str(L)+"x"+str(L)
    beta=2.*L
    N=L*L
    hbi=1.375
    hbf=1.575
    nsteps=16
    hh=np.arange(hbi,hbf, (hbf-hbi)/nsteps )
    samples = 500
    
    for h in np.arange(0,len(hh)):
        fn=np.loadtxt("../data/"+sys+"/DIST/TIM-"+sys+"h"+str(h)+"p0beta0")  #unpack sorts it according to columns
    
        rho_col=(0)*nobs+7
        rhosq_col=(0)*nobs+8
        m_col =(0)*nobs+3
        msq_col = (0)*nobs+4

        plt.title('Energy susceptibility')
        Esus= (fn[12] - fn[11]*fn[11] - fn[11])/(beta*hh[h]*hh[h]*N)
        Esus2= (fn[18] - fn[17]*fn[17] - fn[17])/(beta*hh[h]*hh[h]*N)
        Esus3= (fn[20] - fn[19]*fn[19] - fn[19])/(beta*hh[h]*hh[h]*N)
        Esus4= (fn[16] - fn[15]*fn[15] - fn[15])/(beta*hh[h]*hh[h]*N)
        Esus5= (fn[14] - fn[13]*fn[13] - fn[13])/(beta*hh[h]*hh[h]*N)
        
        #Esus=fn[19]
        #Esus= (fn[10] - fn[9]*fn[9] - fn[9])/(beta*hh[h]*hh[h]*N)
        #plt.scatter([2.*hh[h]], [fn[m_col]], c='g')
        plt.scatter([2.*hh[h]], [Esus2], c=col)
        plt.scatter([2.*hh[h]], [Esus], c=col)
    #plt.scatter([2.*hh[h]], [Esus2], c='r')
    #plt.scatter([2.*hh[h]], [Esus3], c='r')
    #plt.scatter([2.*hh[h]], [Esus4], c='b')
    #plt.scatter([2.*hh[h]], [Esus5], c='b')
    #plt.show()
plt.show()
    
    
    
