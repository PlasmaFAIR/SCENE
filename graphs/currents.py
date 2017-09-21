import numpy as np
import matplotlib.pyplot as plt


def curplot(device):
    
    filename = device+"_currents.dat"
    currents = np.loadtxt(filename,skiprows=1)

    R = currents[:,0]
    tot = currents[:,1]
    ext = currents[:,2]
    bs = currents[:,3]
    di = currents[:,4]
    ps = currents[:,5]
    nb = currents[:,6]


    plt.plot(R, tot, label='Total Current')
    plt.plot(R, ext, label='External')
    plt.plot(R, bs, label='Bootstrap')
    plt.plot(R, di, label='Diamagnetic')
    plt.plot(R, ps, label='Pfirsch-Schluter')
    plt.plot(R, nb, label='Neutral Beam Current')
    

    plt.xlabel('R (m)')
    plt.ylabel('Current (kA m^(-2))')
    plt.legend(loc=2)
    plt.grid()

    plt.savefig('currents.png')
    plt.show()
    plt.close()


    plt.plot(R, nb, label='Neutral Beam Current')

    plt.xlabel('R (m)')
    plt.ylabel('Current (kA m^(-2))')
    plt.title('Neutral beam current')
    plt.show()
    
