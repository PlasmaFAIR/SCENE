import numpy as np
import matplotlib.pyplot as plt
import os.path
import pickle

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

        psi=data[:,0]
        jtot=data[:,1]/1000000
        jnb=data[:,2]/1000000
        jbs=data[:,3]/1000000

        rho = np.sqrt(psi/max(psi))

        plt.plot(rho, jtot-jbs, label='Needed')
        plt.plot(rho, jnb, label='NBI')
        plt.plot(rho, jbs, label='BS')
        plt.plot(rho, jtot, label='Total')
        plt.xlabel('$ \\rho_{\\psi}$')
        plt.title('$\\frac{<J.B>}{<B^2>}$ profiles')
        plt.legend()
        plt.show()



def nbi_radial(device):
    
    filename=device+"_flxav.dat"
    if ( os.path.isfile(filename)):
        rawdata=np.loadtxt(filename, skiprows=1)

 
        psi=rawdata[:,0]
        nbc=rawdata[:,1]/1000000
        nb_e=rawdata[:,2]
        nb_i=rawdata[:,3]
        nb_den=rawdata[:,4]*1e-19
        nb_vel=rawdata[:,5]
        flxvol=rawdata[:,6]
        nb_fast=rawdata[:,7]/1000000

        pow_e=sum(nb_e*flxvol)
        pow_i=sum(nb_i*flxvol)

        rho = np.sqrt(psi/max(psi))

        
    filename=device+"_srcfast.dat"
    if ( os.path.isfile(filename)):
        rawdata=np.loadtxt(filename, skiprows=1)
        
        data = rawdata/1.6
        print(np.shape(rawdata), np.shape(rho))
        ax = plt.subplot(111)
        dep = np.sum(data,axis=1)
        plt.plot(rho,dep)
        plt.xlabel('$\\rho_\\psi$')
        plt.ylabel('FI source $\\# m^{-3}s^{-1}$')
        #pickle.dump(ax, open('nbidep.pickle', 'wb'))
        plt.show()
        #plt.close()

        '''
        plt.plot(rho,nb_den)
        plt.xlabel('$\\rho_\\psi$')
        plt.ylabel('Density $m^{-3}$')
        #plt.show()
        plt.close()

        plt.plot(rho,nb_vel)
        plt.xlabel('$\\rho_\\psi$')
        plt.ylabel('Vel $m^{-2}s^{-1}$')
        #plt.show()
        plt.close()
               

        plt.plot(rho**2,nbc)
        plt.xlabel('$\\psi$')
        plt.ylabel('Current $MAm^{-2}$')
        #plt.show()
        plt.close()


        plt.plot(rho,nb_fast)
        plt.xlabel('$\\rho_\\psi$')
        plt.ylabel('Current $MAm^{-2}$')
        plt.title('Current without back electron effect')
        #plt.show()
        plt.close()

        
        plt.plot(rho, nbc/(2*np.pi*0.85))
        plt.xlabel('$\\rho_\\psi$')
        plt.ylabel('Current $MAm^{-3}$')
        # plt.show()
        plt.close()
        '''
        plt.plot(rho, nb_e, label='Total power to electrons %5.2f MW'%pow_e)
        plt.plot(rho,nb_i, label='Total power to ions %5.2f MW'%pow_i)
        plt.plot(rho,nb_i+nb_e, label='Total power MW')
        plt.xlabel('$\\rho_\\psi$')
        plt.ylabel('Power ($MWm^{-2}$)')
        plt.legend()
        plt.show()
        #plt.close()

        print(sum(flxvol))

        
if __name__ == "__main__":

    device = input('Enter device name: ')
    curplot(device)
    nbi_radial(device)

    

    
