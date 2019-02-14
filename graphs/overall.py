import numpy as np
import matplotlib.pylab as plt


def plots(device):

    filename = device+'_overall.dat'
    data = np.loadtxt(filename, skiprows=1)

    #length = len(data)/3

    #data = np.reshape(data, ( length,9))

    R = data[:, 0]
    ne = data[:, 1]
    ni = data[:, 2]
    B_t = data[:, 3]
    B_p = data[:, 4]
    B = data[:, 5]
    q = data[:, 6]
    T = data[:, 7]
    p = data[:, 8]

    plt.plot(R, ne, 'r', label='$n_e$', linewidth=2)
    plt.plot(R, ni, 'b', label='$n_i$', linewidth=2)
    plt.legend(loc=1)
    plt.xlabel('R (m)', fontsize=20)
    plt.ylabel('Density $(m^{-3})$', fontsize=20)
    plt.title('Electron density', fontsize=24)
    plt.xticks(size=20)
    plt.yticks(size=20)
    plt.savefig('densities.png')
    plt.show()

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
    plt.show()

    fig, ax1 = plt.subplots()

    ax1.plot(R, T, 'r--', label='$Temp$', linewidth=2)
    ax1.set_xlabel('R (m)', fontsize=20)
    ax1.set_ylabel('Temp (keV)', fontsize=20)
    ax1.tick_params('y', colors='r')
    plt.xticks(size=20)
    plt.yticks(size=20)
    ax2 = ax1.twinx()
    ax2.plot(R, p, 'b.', label='$Pressure$')
    ax2.set_ylabel('Pressure', color='b')
    ax2.tick_params('y', colors='b')
    plt.title('Electron Temperature', fontsize=24)

    plt.savefig('Electron temp and press.png')
    plt.show()
    plt.close()
