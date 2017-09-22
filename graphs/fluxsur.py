import numpy as np
import matplotlib.pyplot as plt


def flxplot(device):

    filename=device+'_fluxsur.dat'

    bdy_dat='bdy.txt'

    bdy=np.loadtxt(bdy_dat, skiprows=1)
    fluxsur = np.loadtxt(filename)


    npts, ncon = fluxsur[0,:]

    fluxsur = np.delete(fluxsur, 0,axis=0)


    fluxsur = np.reshape(fluxsur, (int(ncon),  int(npts), 2 ) )


    for i in range(int(ncon)):

        plt.plot(fluxsur[i,:,0], fluxsur[i,:,1])



    plt.plot(bdy[:,0], bdy[:,1], 'k', linewidth=4.0)
    ax = plt.gca()
    ax.set_aspect('equal')
    plt.xlabel('R (m)')
    plt.ylabel('Z (m)')
    plt.title(' Flux surfaces of ')
    #plt.show()
    plt.savefig('flxsur.png')
    plt.close()

    #print('Flux Surface plot saved')
