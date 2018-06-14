import sys
import os
try:
  pflotran_dir = os.environ['PFLOTRAN_DIR']
except KeyError:
  print('PFLOTRAN_DIR must point to PFLOTRAN installation directory and be defined in system environment variables.')
  sys.exit(1)
sys.path.append(pflotran_dir + '/src/python')
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
import pflotran as pft
from h5py import *

fig = plt.figure(figsize=(8,6))
plt.subplot(1,1,1)
fig.suptitle("Comparison: Decay",fontsize=16)
plt.xlabel('Time [y]')
plt.ylabel('A [M]')

linewidth = 2

# analytical
path = 'analytical/'
filename = 'decay_only_analytical-obs-0.tec'
data = pft.Dataset(path+filename,1,2)
plt.plot(data.get_array('x'),data.get_array('y'),
         label='ufd analytical',lw=linewidth)

# implicit
path = 'implicit/'
filename = 'decay_only_implicit-obs-0.tec'
data = pft.Dataset(path+filename,1,2)
plt.plot(data.get_array('x'),data.get_array('y'),
         label='ufd implicit',lw=linewidth)

# implicit
path = 'conventional/'
filename = 'decay_only_conventional-obs-0.tec'
data = pft.Dataset(path+filename,1,2)
plt.plot(data.get_array('x'),data.get_array('y'),
         label='conventional',ls=':',lw=linewidth)

num_values = 1001
dt = 12.5/(num_values-1)
k = 1.e-8
yr_to_sec = 3600.*24.*365.
x = np.arange(0,num_values*dt,dt)
y = np.zeros(num_values)
for i in range(num_values):
  y[i] = 1.e-5*math.exp(-k*i*dt*yr_to_sec)
plt.plot(x,y,label='analytical solution',ls='--',lw=linewidth)

'''
# conventionsl small time step
path = 'conventional_sm_ts/'
filename = 'decay_only_conventional-obs-0.tec'
data = pft.Dataset(path+filename,1,2)
plt.plot(data.get_array('x'),data.get_array('y'),
         label='conventional',ls=':',lw=4,color='black')
'''


#plt.xlim(0.,1.)
#plt.ylim(4.8,8.2)
#plt.grid(True)

#'best'         : 0, (only implemented for axis legends)
#'upper right'  : 1,
#'upper left'   : 2,
#'lower left'   : 3,
#'lower right'  : 4,
#'right'        : 5,
#'center left'  : 6,
#'center right' : 7,
#'lower center' : 8,
#'upper center' : 9,
#'center'       : 10,
plt.legend(loc=1,title='Scenario')
# xx-small, x-small, small, medium, large, x-large, xx-large, 12, 14
#plt.setp(plt.gca().get_legend().get_texts(),fontsize='small')
#plt.setp(plt.gca().get_legend().get_texts(),linespacing=0.)
plt.setp(plt.gca().get_legend().get_frame().set_fill(False))
plt.setp(plt.gca().get_legend().draw_frame(False))
#plt.gca().yaxis.get_major_formatter().set_powerlimits((-1,1))

fig.subplots_adjust(hspace=0.2,wspace=0.2,
                  bottom=.12,top=.9,
                  left=.12,right=.9)

plt.show()
