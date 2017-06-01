import numpy as np
import currents
import fluxsur
import overall
import bdytest

device = raw_input('Enter device name: ')

print(device, 'graphs will be produced')


#Produce fourier fits
bdytest.fourier_test()

#Produce currents plot
currents.curplot(device)

#Produce flux plot
fluxsur.flxplot(device)


#Produce other plots
overall.plots(device)





