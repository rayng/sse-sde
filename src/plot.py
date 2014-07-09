
# This script plots stiffness and compressibility versus field value
# It can be easily modified to include other observables as well.
# Not for beta doubling.

import numpy as np
import matplotlib.pyplot as plt



for L in (8,12):
    if(L==8):
        col='r'
    else:
        col='b'
    nmax=10
    nobs=24
    sys=str(L)+"x"+str(L)
    #sys=str(L)+"x1"
    beta=2.*L
    N=L*L
    #N=L
    hbi=1.35  #0.25
    hbf=1.6  #0.75
    nsteps=16
    hh=np.arange(hbi,hbf, (hbf-hbi)/nsteps )
    samples = 500
    
    for h in np.arange(0,len(hh)):
        fn=np.loadtxt("../data/"+sys+"/DIST/TIM-"+sys+"h"+str(h)+"p0beta"+str(int(beta)))  #unpack sorts it according to columns

        rho_col=(0)*nobs+7
        rhosq_col=(0)*nobs+8
        m_col =(0)*nobs+3
        msq_col = (0)*nobs+4
        
        plt.title('Energy susceptibility')
        Esus= (fn[12] - fn[11]*fn[11] - fn[11])/(beta*hh[h]*hh[h]*N)     # hop 0 
        Esus2= (fn[18] - fn[17]*fn[17] - fn[17])/(beta*hh[h]*hh[h]*N)    # hop +1
        Esus3= (fn[20] - fn[19]*fn[19] - fn[19])/(beta*hh[h]*hh[h]*N)    # hop -1
        Esus4= (fn[16] - fn[15]*fn[15] - fn[15])/(beta*hh[h]*hh[h]*N)    # dis +1
        Esus5= (fn[14] - fn[13]*fn[13] - fn[13])/(beta*hh[h]*hh[h]*N)    # dis -1
        
        fs0 = (fn[21] - fn[11]*fn[11]/(8.*hh[h]*hh[h]))/N             # fs-hop 0
        fsp = (fn[22] - fn[17]*fn[17]/(8.*hh[h]*hh[h]))/N             # fs-hop +
        fsm = (fn[23] - fn[19]*fn[19]/(8.*hh[h]*hh[h]))/N             # fs-hop -
        
        #Esus=fn[19]
        #Esus= (fn[10] - fn[9]*fn[9] - fn[9])/(beta*hh[h]*hh[h]*N)
        #plt.scatter([2.*hh[h]], [fn[m_col]], c='g')

        #plt.scatter([2.*hh[h]], [fs0/4], c='g')
        #plt.scatter([2.*hh[h]], [fsm], c='y')
        plt.scatter([2.*hh[h]], [fsp/4], c=col)
        #plt.scatter([2.*hh[h]], [fsm/2], c='r')

        #print 2.*hh[h], fsm/2
        
        #print fn[23], fn[22], fn[21]
        #plt.scatter([2.*hh[h]], [fn[22]], c=col)
        #plt.scatter([2.*hh[h]], [fn[23]], c=col)
        #plt.scatter([2.*hh[h]], [Esus2], c='b')
    #plt.scatter([2.*hh[h]], [Esus3], c='r')
    #plt.scatter([2.*hh[h]], [Esus4], c='b')
    #plt.scatter([2.*hh[h]], [Esus5], c='b')

    plt.yticks(np.arange(0,0.55,0.05))
    plt.xlim(2.70,3.15,0.1)
    plt.ylim(0,0.55,0.05)
plt.show()

#plt.axhline(y=0.475)
#plt.axhline(y=0.15)

#ff=np.loadtxt("fsm", unpack=True)
#plt.xticks(np.arange(0.5,1.5,0.1))
#plt.yticks(np.arange(0,0.55,0.05))
#plt.plot( ff[0], ff[1]/2, 'ro')
#plt.xlim(0.5,1.5,0.1)
#plt.ylim(0,0.55,0.05)
#plt.show()
    
    
    
