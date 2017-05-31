import numpy as np
import currents
import fluxsur
import overall

device = raw_input('Enter device name: ')

print(device, 'graphs will be produced')


#Produce currents plot
currents.curplot(device)

#Produce flux plot
fluxsur.flxplot(device)


#Produce other plots
overall.plots(device)



