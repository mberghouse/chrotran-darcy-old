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
fig.suptitle("Comparison: 1D Reactive Transport",fontsize=16)
plt.xlabel('X [m]')
plt.ylabel('U234L [M]')

# 1D structured transport only
path = '1d_structured_uniform/'
filename = 'rt-014.tec'
data = pft.Dataset(path+filename,1,7)
plt.plot(data.get_array('x'),data.get_array('y'),
         label='structured - transport only')


# 1D structured with flow
path = '1d_structured/'
filename = 'rt-014.tec'
data = pft.Dataset(path+filename,1,13)
plt.plot(data.get_array('x'),data.get_array('y'),
         label='structured',ls='--')

# 1D unstructured uniform
grid_path = '../transport_only/1d_unstructured_uniform/'
path = '1d_unstructured_uniform/'
filename = 'grid.uge'
f = open(grid_path+filename,'r')
num_cells = int(f.readline().split()[1])
x = np.zeros(num_cells)
for i in range(num_cells):
  x[i] = float(f.readline().split()[1])
f.close()
filename = 'rt.h5'
h5file = File(path+filename,mode='r')
y = h5file['  14 Time  1.25000E+01 y/Total U234L [M]']
plt.plot(x,y,label='unstructured - uniform',ls='-.')
h5file.close()

# 1D unstructured non-uniform
grid_path = '../transport_only/1d_unstructured_irregular/'
path = '1d_unstructured_irregular/'
filename = 'grid.uge'
f = open(grid_path+filename,'r')
num_cells = int(f.readline().split()[1])
x = np.zeros(num_cells)
for i in range(num_cells):
  x[i] = float(f.readline().split()[1])
f.close()
filename = 'rt.h5'
h5file = File(path+filename,mode='r')
y = h5file['  14 Time  1.25000E+01 y/Total U234L [M]']
plt.plot(x,y,label='unstructured - irregular')
h5file.close()

# 1D unstructured implicit
grid_path = '../transport_only/1d_unstructured_implicit/'
path = '1d_unstructured_implicit/'
filename = 'grid.ugi'
f = open(grid_path+filename,'r')
strings = f.readline().split()
num_cells = int(strings[0])
num_vertices = int(strings[1])
x = np.zeros(num_cells)
for i in range(num_cells):
  f.readline()
x[0] += float(f.readline().split()[0])
for i in range(1,num_cells+1):
  value = float(f.readline().split()[0])
  x[i-1] += value
  if i < num_cells:
    x[i] += value
f.close()
x *= 0.5
filename = 'rt.h5'
h5file = File(path+filename,mode='r')
y = h5file['  14 Time  1.25000E+01 y/Total U234L [M]']
plt.plot(x,y,label='unstructured - implicit',ls=':')
h5file.close()

# 1D structured transport only
path = '1d_structured_uniform_sm_ts/'
filename = 'rt-014.tec'
data = pft.Dataset(path+filename,1,7)
plt.plot(data.get_array('x'),data.get_array('y'),
         label='structured - small time step',
         ls=':',lw=1,color='black')


#plt.xlim(40.,60.)
#plt.ylim(0,0.1)
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
plt.legend(loc=2,title='Scenario')
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
