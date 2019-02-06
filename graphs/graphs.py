import numpy as np
import currents
import fluxsur
import overall

device = input('Enter device name: ')

print(device, 'graphs will be produced')


#Produce fourier fits
#bdytest.fourier_test()


#Produce other plots
overall.plots(device)

#Produce currents plot
currents.curplot(device)

currents.nbi_radial(device)

#Produce flux plot
fluxsur.flxplot(device)







