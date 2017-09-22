import numpy as np
import matplotlib.pylab as plt

def plots(device):


    filename = device+'_overall.dat'
    data = np.loadtxt(filename,skiprows=1)

    length = len(data)/3

    data = np.reshape(data, ( length,9))

    R = data[:,0]
    ne = data[:,1]
    ni = data[:,2]
    B_t = data[:,3]
    B_p = data[:,4]
    B = data[:,5]
    q = data[:,6]
    T = data[:,7]
    p = data[:,8]


    plt.plot(R, ne, label='$n_e$')
    plt.plot(R,ni, label='$n_i$')
    plt.legend(loc=1)
    plt.xlabel('R (m)')
    plt.ylabel(r'Density $(m^{-3})$')
    plt.savefig('densities.png')
    plt.clf()
    
    plt.plot(R, B_t, label=r'$B_\phi$')
    plt.plot(R, B_p, label=r'$B_\theta$')
    plt.plot(R, B, label='$B$')
    plt.xlabel('R (m)')
    plt.ylabel('Field (T)')
    plt.title('Magnetic fields')
    plt.legend(loc=1)
    plt.savefig('magfield.png')
    plt.clf()

    plt.plot(R, q)
    plt.title('Safety Factor')
    plt.xlabel('R (m) ')
    plt.ylabel('Safety factor')
    plt.savefig('Safety factor.png')
    plt.clf()


    fig, ax1 = plt.subplots()

    ax1.plot(R, T, 'r--', label='$Temp$')
    ax1.set_xlabel('R (m)')
    ax1.set_ylabel('Temp (keV)', color='r')
    ax1.tick_params('y', colors='r')

    ax2 = ax1.twinx()
    ax2.plot(R, p, 'b.', label='$Pressure$')
    ax2.set_ylabel('Pressure', color='b')
    ax2.tick_params('y', colors='b')
    plt.title('Electron Pressure and Temp')

    plt.savefig('Electron temp and press.png')
    plt.close()


