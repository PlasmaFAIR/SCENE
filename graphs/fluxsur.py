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
