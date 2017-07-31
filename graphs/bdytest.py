import numpy as np
import matplotlib.pyplot as plt


def fourier_test(): 

    f1 = open('bdy.txt', 'r')
    
    ndat = int(f1.readline())

    f1.close()


    bdy = np.loadtxt('bdy.txt',skiprows=1)






    f = open('bdytest.dat', 'r')

    npts = int(f.readline())


    theta = np.zeros(npts)
    ztmp = np.zeros(npts)
    rtmp = np.zeros(npts)

    thdim = np.zeros(ndat)

    for i in range(npts):

        data = f.readline().split()
        theta[i] = float(data[0])
        ztmp[i] = float(data[1])
        rtmp[i] = float(data[2])


    for j in range(ndat):

        thdim[j] = f.readline()

    f.close()


    diff = bdy[0,0] - rtmp[0]

    bdy[:,0] = bdy[:,0] - diff



    plt.plot(theta, ztmp)
    plt.xlabel('Poloidal angle (rad)')
    plt.ylabel('Z (m)')
    plt.plot(thdim,bdy[:,1])
    plt.title('Fourier components fit in Z')
    plt.savefig('fourier_z.png')
    plt.close()

    plt.plot(theta, rtmp)
    plt.plot(thdim, bdy[:,0])
    plt.xlabel('Poloidal angle (rad)')
    plt.ylabel('R (m)')
    plt.title('Fourier components fit in R')
    plt.savefig('fourier_r.png')
    plt.close()
