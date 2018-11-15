import numpy as np
import matplotlib.pyplot as plt
import os.path

def curplot(device):
    
    filename = device+"_currents.dat"
    currents = np.loadtxt(filename,skiprows=1)

    R = currents[:,0]
    tot = currents[:,1]/1000.
    ext = currents[:,2]/1000.
    bs = currents[:,3]/1000.
    di = currents[:,4]/1000.
    ps = currents[:,5]/1000.
    nb = currents[:,6]/1000.
    ext2 = currents[:,7]/1000.

    plt.plot(R, tot, label='Total Current')
    plt.plot(R, ext, label='External')
    plt.plot(R, ext2, label='External  2')
   
    plt.plot(R, bs, label='Bootstrap')
    plt.plot(R, di, label='Diamagnetic')
    plt.plot(R, ps, label='Pfirsch-Schluter')
    plt.plot(R, nb, label='Neutral Beam Current')
    

    plt.xlabel('R (m)')
    plt.ylabel('Current (MA m$^{-2}$)')
    plt.legend(loc=4,fontsize=10)
    plt.grid()

    plt.savefig('currents.png')
    plt.show()



    plt.plot(R, nb, label='Neutral Beam Current', linewidth=2)
    #plt.plot(R, ext, label='Extra needed')
    plt.plot(R, tot-bs-di-ps, label='Total external needed', linewidth=2)
    plt.xlabel('R (m)', fontsize=20)
    plt.legend()
    plt.ylabel('Current (MA m$^{-2}$)', fontsize=20)
    plt.xticks(size=20)
    plt.yticks(size=20)
    plt.title('Neutral beam current',fontsize=24)
    plt.savefig('nbi.png')
    plt.show()
    plt.clf()

    file2=device+"_currentprof.dat"
    if ( os.path.isfile(file2)):
        data=np.loadtxt(file2)

        rho=data[:,0]
        eps=data[:,1]
        jtot=data[:,2]/1000000
        jnb=data[:,3]/1000000
        te=data[:,4]/1000

        print(np.mean(te))
        plt.plot(eps, jtot)
        plt.plot(eps, jnb)
        #plt.show()
        plt.clf()

        plt.plot(rho, te)
        #plt.show()
        plt.clf()


def nbi_radial(device):
    
    filename=device+"_flxav.dat"
    if ( os.path.isfile(filename)):
        rawdata=np.loadtxt(filename)

 
        rho=rawdata[:,0]
        nbc=rawdata[:,1]/1000000
        nb_e=rawdata[:,2]
        nb_i=rawdata[:,3]
        nb_den=rawdata[:,4]*1e-19
        nb_vel=rawdata[:,5]
        flxvol=rawdata[:,6]
        nb_fast=rawdata[:,7]/1000000
    
        plt.plot(rho,nb_den)
        plt.xlabel('$\\rho$')
        plt.ylabel('Density $m^{-3}$')
        #plt.show()
        plt.close()
        

        plt.plot(rho,nbc)
        plt.xlabel('$\\rho$')
        plt.ylabel('Current $MAm^{-2}$')
        #plt.show()
        plt.close()


        plt.plot(rho,nb_fast)
        plt.xlabel('$\\rho$')
        plt.ylabel('Current $MAm^{-2}$')
        plt.title('Current without back electron effect')
        #plt.show()
        plt.close()

        
        plt.plot(rho, nbc/(2*np.pi*0.85))
        plt.xlabel('$\\rho$')
        plt.ylabel('Current $MAm^{-3}$')
        # plt.show()
        plt.close()

        plt.plot(rho, nb_e, label='Power to Electrons')
        plt.plot(rho,nb_i, label='Power to Ions')
        plt.xlabel('$\\rho$')
        plt.ylabel('Power ($MWm^{-2}$)')
        plt.legend()
        #plt.show()
        plt.close()
