import numpy as np
import currents
import fluxsur
import overall
import bdytest

device = raw_input('Enter device name: ')

print(device, 'graphs will be produced')

filename = device+'_safety.dat'
data = np.loadtxt(filename)


R = data[:,0]
q = data[:,1]

    
plt.plot(R, q)
plt.title('Safety Factor')
plt.xlabel('R (m) ')
plt.ylabel('Safety factor')
plt.savefig('Safety factor.png')
plt.clf()







