import numpy as np
import itertools as it
import matplotlib.pyplot as plt

out=np.loadtxt("../data/16x1/SDE/SDE-TIM-h0p0beta0", dtype=np.str)

i=iter(out)
nT=10000
nV=16
nn=2*2*16+1
Nsamples=1

trj1= np.shape(list(it.islice(i,0,nT)))


#zav=np.zeros(shape=(nn,nT))
zav=np.zeros(shape=(nn,nT))

Sz = np.zeros(shape=(nV,nT),dtype=np.complex_)
Sy = np.zeros(shape=(nV,nT),dtype=np.complex_)
Sx = np.zeros(shape=(nV,nT),dtype=np.complex_)

zcmp = np.zeros(shape=(nV,nT), dtype=np.complex_)
zcmppr = np.zeros(shape=(nV,nT), dtype=np.complex_)
R = np.zeros(shape=(nV,nT), dtype=np.complex_)



with open('../data/16x1/SDE/SDE-TIM-h0p0beta0') as t_in:
    
    for ntraj in np.arange(0,Nsamples):
        out = np.genfromtxt(it.islice(t_in, nT))
        #print np.shape(out)
        out2=np.transpose(out)
        print np.shape(out2)
        
        for i in np.arange(0,nV):
            zcmp[i]=out2[4*i+1] + 1j*out2[4*i+2]
            zcmppr[i]=out2[4*i+3] + 1j*out2[4*i+4]

        #print np.shape(zcmp)

        
        R = (zcmp+zcmppr)*0.5
        
        Sxtemp=0.5*(np.cosh(zcmp)- np.sinh(zcmp)*np.tanh(R) )
        Sytemp=0.5*1j*(np.sinh(zcmp) - np.cosh(zcmp)*np.tanh(R)  )
        Sztemp=0.5*np.tanh(R)
        
        Sx=Sx+Sxtemp
        Sy=Sy+Sytemp
        Sz=Sz+Sztemp

        #print 0.5*np.tanh(R)[0]
        #plt.plot(out2[0], Sytemp[0].real )
        
        #look at variable 0#
        #plt.plot(out2[0], out2[4*i+1],'ro-')
         #plt.plot(out2[0], out2[4*i+2], 'bo-')
              #plt.plot(out2[4*i+1], out2[4*i+2], 'bo-')
#plt.plot(out2[0], out2[4*i+3],'ro-')
    #plt.plot(out2[0], out2[4*i+4], 'bo-')
    #plt.show()
    
        zav = zav+ out2


    
#plt.show()

zav=zav/Nsamples
Sz=Sz/Nsamples
Sy=Sy/Nsamples
Sx=Sx/Nsamples




# plot some observables

#print Sz

#print Sz

#for i in np.arange(0,nV):
    #plt.plot(zav[0], Sx[i].real, label="site "+str(i) )
    #plt.plot(zav[0], Sy[i].real, label="site "+str(i) )
    #plt.plot(zav[0], Sz[i].real, label="site "+str(i) )
    
    #plt.plot(zav[0], Sx[i].real, 'b')
    #plt.plot(zav[0], Sy[i].real, 'g')
#plt.legend(loc='best')
#plt.show()


# Lattice averaged
plt.plot( zav[0], np.mean(Sx, axis=0), 'ro-' , markevery=200,label=r'$\langle S_x \rangle$' )
plt.plot( zav[0], np.mean(Sy, axis=0),  'bo-', markevery=200, label=r'$\langle S_y \rangle$' )
plt.plot( zav[0], np.mean(Sz, axis=0), 'yo-', markevery=200, label=r'$\langle S_z \rangle$' )
plt.title(str(Nsamples)+' samples')
plt.savefig('../figs/'+str(Nsamples)+'samples-exactdrift.png')
plt.legend(loc='best')

plt.show()

    
#plt.plot(zav[2*i+1], zav[2*i+2], 'ro-')
    #plt.plot(zav[0], zav[2*i+1], 'ro-')
    #plt.plot(zav[0], zav[2*i+2], 'bo-')
    #plt.show()
    
#print len(out2[2*i+1])
        
    


