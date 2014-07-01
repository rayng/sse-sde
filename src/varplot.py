import numpy as np
import itertools as it
import matplotlib.pyplot as plt


h=0
nV=16
out=np.loadtxt("../data/"+str(nV)+"x1/SDE/SDE-TIM-h"+str(h)+"p0beta0", dtype=np.str)
i=iter(out)
nT=3000

nn=2*2*nV+1
Nsamples=50

trj1= np.shape(list(it.islice(i,0,nT)))

zav=np.zeros(shape=(nn,nT))
Sz = np.zeros(shape=(nV,nT),dtype=np.complex_)
Sy = np.zeros(shape=(nV,nT),dtype=np.complex_)
Sx = np.zeros(shape=(nV,nT),dtype=np.complex_)

zcmp = np.zeros(shape=(nV,nT), dtype=np.complex_)
zcmppr = np.zeros(shape=(nV,nT), dtype=np.complex_)
R = np.zeros(shape=(nV,nT), dtype=np.complex_)



with open('../data/'+str(nV)+'x1/SDE/SDE-TIM-h0p0beta0') as t_in:
    
    for ntraj in np.arange(0,Nsamples):

        out = np.genfromtxt(it.islice(t_in, nT))
        print np.shape(out)
        out2=np.transpose(out)
        print np.shape(out2)


        tarr=out2[0]
        for i in np.arange(0,nV):
            zcmp[i]=out2[4*i+1] + 1j*out2[4*i+2]
            zcmppr[i]=out2[4*i+3] + 1j*out2[4*i+4]

        #print np.shape(zcmp)
        
        R= (zcmp+zcmppr)*0.5
        S= (zcmp-zcmppr)*(0.5*1j)
        
        Sxtemp=0.5*(np.cosh(zcmp)- np.sinh(zcmp)*np.tanh(R) )
        Sytemp=0.5*1j*(np.sinh(zcmp) - np.cosh(zcmp)*np.tanh(R)  )
        Sztemp=0.5*np.tanh(R)
        
        Sx=Sx+Sxtemp
        Sy=Sy+Sytemp
        Sz=Sz+Sztemp

        #print 0.5*np.tanh(R)[0]
        #plt.plot(out2[0], Sytemp[0].real )
        
        #look at variable 0#
        
        
        for i in np.arange(0,nV):
            """
            z-variables
            """
            plt.plot(out2[0], out2[4*i+1],'r--')
            plt.plot(out2[0], out2[4*i+2], 'b--')
            plt.plot(out2[0], out2[4*i+3],'r--')
            plt.plot(out2[0], out2[4*i+4], 'b--')
            """
            R,S-variables
            """
            #plt.plot(out2[0], R[i].real,'r--')
            #plt.plot(out2[0], R[i].imag, 'b--')
            #plt.plot(out2[0], S[i].real,'r--')
            #plt.plot(out2[0], S[i].imag, 'b--')
            

    #plt.show()
    #plt.xlim(-15,15)
    
    plt.legend(loc='best')
    
    
    zav = zav+ out2


    
#plt.show()

plt.clf()
plt.plot(tarr, 'ro-');
#plt.show()

zav=zav/Nsamples
Sz=Sz/Nsamples
Sy=Sy/Nsamples
Sx=Sx/Nsamples



print np.shape(Sx), np.shape(Sy), np.shape(Sz)

# plot some observables


"""
for i in np.arange(0,nV):
    plt.plot(zav[0], Sx[i].real, label="site "+str(i) )
plt.legend(loc='best')
plt.savefig('Sx-sites.png')
plt.show()

for i in np.arange(0,nV):
    plt.plot(zav[0], Sy[i].real, label="site "+str(i) )
plt.legend(loc='best')
plt.savefig('Sy-sites.png')
plt.show()


for i in np.arange(0,nV):
    plt.plot(zav[0], Sx[i].real, label="site "+str(i) )
plt.legend(loc='best')
plt.savefig('Sz-sites.png')
plt.show()
"""

    #plt.plot(zav[0], Sy[i].real, label="site "+str(i) )
    #plt.plot(zav[0], Sz[i].real, label="site "+str(i) )
    
    #plt.plot(zav[0], Sx[i].real, 'b')
    #plt.plot(zav[0], Sy[i].real, 'g')



# Lattice averaged
plt.clf()



Sxav=np.mean(Sx, axis=0).real
Syav=np.mean(Sy, axis=0).real
Szav=np.mean(Sz, axis=0).real


plt.plot( tarr, Sxav, 'ro-' , markevery=100,markersize=1,label=r'$\langle S_x \rangle$' )
plt.plot( tarr, Syav,  'bo-', markevery=100,markersize=1, label=r'$\langle S_y \rangle$' )
plt.plot( tarr, Szav, 'yo-', markevery=100, markersize=1,label=r'$\langle S_z \rangle$' )
plt.title(str(Nsamples)+' samples')
plt.savefig('../figs/'+str(Nsamples)+'samples-polar-x.png')
plt.legend(loc='best')


fp=open("../data/"+str(nV)+"h"+str(h)+"obs.dat","w")
for t in np.arange(0,len(Sxav)):
    fp.write('%f %f %f %f \n' % (tarr[t], Sxav[t], Syav[t], Szav[t]) )

fp.close()
#plt.show()

    
#plt.plot(zav[2*i+1], zav[2*i+2], 'ro-')
    #plt.plot(zav[0], zav[2*i+1], 'ro-')
    #plt.plot(zav[0], zav[2*i+2], 'bo-')
    #plt.show()
    
#print len(out2[2*i+1])
        
    


