import numpy as np
import matplotlib.pyplot as plt
import os.path


def flxplot(device):

    filename=device+'_fluxsur.dat'
         
    fluxsur = np.loadtxt(filename)


    npts, ncon = fluxsur[0,:]

    fluxsur = np.delete(fluxsur, 0,axis=0)


    fluxsur = np.reshape(fluxsur, (int(ncon),  int(npts), 2 ) )


    for i in range(int(ncon)):

        plt.plot(fluxsur[i,:,0], fluxsur[i,:,1])

        

    if ( os.path.isfile('bdy.txt')):
         bdy_dat='bdy.txt'

         bdy=np.loadtxt(bdy_dat, skiprows=1)

         plt.plot(bdy[:,0], bdy[:,1], 'k', linewidth=4.0)

         
    ax = plt.gca()
    ax.set_aspect('equal')
    plt.xlabel('R (m)')
    plt.ylabel('Z (m)')
    plt.title(' Flux surfaces of ')
    plt.savefig('flxsur.png')
    plt.show()

    print('Flux Surface plot saved')


def jettodata(device):

    filename = device+'.jetto'

    data =np.loadtxt(filename, skiprows=1)

    phi = data[:,0]
    ne =data[:,1]*1.0e-19
    te =data[:,3]*1.0e-3
    rhotor = np.sqrt(phi/max(phi))
    q = data[:,5]



    
    plt.plot(rhotor,te, 'r--', lw=5)
    plt.title('Electron temperature profile', fontsize=24)
    plt.xlabel('$ \\rho_{\\phi}$', fontsize=20)
    plt.ylabel('$keV$', fontsize=20)
    plt.xticks(size=20)
    plt.yticks(size=20)
    plt.savefig('te_torflux.png')
    plt.show()
   
    plt.plot(rhotor,ne, 'r--', lw=5)
    plt.title('Electron density profile', fontsize=24)
    plt.xlabel('$ \\rho_{\\phi}$', fontsize=20)
    plt.ylabel('# $10^{19} m^{-3}$', fontsize=20)
    plt.xticks(size=20)
    plt.yticks(size=20)
    plt.savefig('ne_torflux.png')
    plt.show()
   
    plt.plot(rhotor,q, 'r--', lw=5)
    plt.title('Safety Factor', fontsize=24)
    plt.xlabel('$ \\rho_{\\phi}$', fontsize=20)
    plt.ylabel('q', fontsize=20)
    plt.xticks(size=20)
    plt.yticks(size=20)
    plt.savefig('Safetyfactor_torflux.png')
    plt.show()


if __name__ == "__main__":

    device = input('Enter device name: ')

    flxplot(device)
    #jettodata(device)

    
