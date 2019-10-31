from freegs import geqdsk
from freegs import machine
from freegs import critical
from freegs.plotting import plotEquilibrium, plotCoils
import numpy as np
import matplotlib.pyplot as plt



#Rough positions of coils in STPP
coils = [("Div", machine.Circuit([("DivL", machine.Coil(1,-5.7),1.0),
                                  ("DivU", machine.Coil(1,5.7),1.0)])),
        ("P1", machine.Circuit([("P1L", machine.Coil(4.9,-4.8),1.0),
                                ("P1U", machine.Coil(4.9,4.8),1.0)])),
        ("P2",machine.Circuit([("P2L", machine.Coil(6.4,-2.2),1.0),
                               ("P2U", machine.Coil(6.4,2.2),1.0)]))#,
         #("P3",machine.Circuit([("P3L", machine.Coil(6.9,-1.2),1.0),
         #                       ("P3U", machine.Coil(6.9,1.2),-1.0)])),
         ]

coils = [("Div", machine.Circuit([("DivL", machine.Coil(2.5,-11.0),1.0),
                                  ("DivU", machine.Coil(2.5,11.0),1.0)])),
        ("P1", machine.Circuit([("P1L", machine.Coil(6.9,-5.8),1.0),
                                ("P1U", machine.Coil(6.9,5.8),1.0)])),
        ("P2",machine.Circuit([("P2L", machine.Coil(8.4,-2.2),1.0),
                               ("P2U", machine.Coil(8.4,2.2),1.0)]))#,
         #("P3",machine.Circuit([("P3L", machine.Coil(6.9,-1.2),1.0),
         #                       ("P3U", machine.Coil(6.9,1.2),-1.0)])),
         ]

'''
#No. of coils + distance
coil_num = 10
dist = 12
coils = []

#Make a semi circle with evenly spaced points
angles = np.linspace(0,np.pi/2,coil_num)
zpos = dist*np.sin(angles)
rpos = dist*np.cos(angles)*0.75

plt.plot(rpos,zpos,'x')
plt.plot(rpos,-zpos,'x')
ax=plt.gca()
ax.set_aspect('equal')
plt.show()

for i in range(coil_num):
    #No coil on axis
    if i!=0:
        coils.append(["U"+str(i),machine.Coil(rpos[i],zpos[i])])
        coils.append(["L"+str(i),machine.Coil(rpos[i],-zpos[i])])
'''

step = machine.Machine(coils)


#with open("g014220.00200") as f:

device = input('Enter device name: ')

filename = device+'.geqdsk'

with open(filename) as f:
    eq = geqdsk.read(f, step, show=True,domain=(0.5,5.0,-5.0,5.0))


#q_before = 
q = critical.find_safety(eq,npsi=200)

plotEquilibrium(eq, show=False)

for i in range(len(step.coils)):
    R = step.coils[i][1].coils[0][1].R
    Z = step.coils[i][1].coils[0][1].Z
    cur = step.coils[i][1].current
    cur = "{:.1f}".format(cur/1e6)
    plt.plot(R,Z, 'ks',markersize=12)
    plt.annotate(cur+" MA", (R+0.5,Z))
    plt.plot(R,-Z, 'ks',markersize=12)
    plt.annotate(cur+" MA", (R+0.5,-Z))
   
plt.xlim([0,8])
plt.ylim([-6,6])


boundary = critical.find_separatrix(eq,ntheta=100)
bdy = np.array(boundary)
rbdy = bdy[:,0]
zbdy = bdy[:,1]

#Start from Max R
ind = np.argmax(rbdy)
rbdy = np.roll(rbdy,-ind)
zbdy = np.roll(zbdy,-ind)

zbdy = zbdy-zbdy[0]
rbdy = np.append(rbdy,rbdy[0])
zbdy = np.append(zbdy,zbdy[0])

#Reverse Arrays
rbdy = rbdy[::-1]
zbdy = zbdy[::-1]


ind = np.argmin(rbdy)
nbdy = len(rbdy)

f = open('freegs_bdy.txt','w')
f.write(str(nbdy)+'\n')

for i in range(ind):
    f.write('{:03.4f} {:03.4f}\n'.format(rbdy[i],zbdy[i]))


for i in range(nbdy-ind):
    f.write('{:03.4f} {:03.4f}\n'.format(rbdy[ind-i],-zbdy[ind-i]))

f.close()



fileOut = "step_freegs.geqdsk"



with open(fileOut,"w+") as fh:
    geqdsk.write(eq,fh)
     


plt.show()

plotCoils(coils)

'''
#tokamak = machine.MAST_sym()
tokamak = machine.TestTokamak()

#with open("g014220.00200") as f:
with open("lsn.geqdsk") as f:
    eq = geqdsk.read(f, tokamak, show=True)

plotEquilibrium(eq)
'''
