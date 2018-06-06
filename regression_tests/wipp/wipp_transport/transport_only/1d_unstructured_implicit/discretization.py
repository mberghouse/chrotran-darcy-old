import sys
from h5py import * 
import numpy

filename = 'grid.ugi'
f = open(filename,'w')

length = 100.
nx = 1000
dx = length/float(nx)
dy = 1.
dz = 1.

ii = [0]*8
ii[0] = 1
ii[1] = 2
ii[2] = ii[1]+nx+1
ii[3] = ii[0]+nx+1
ii[4] = ii[0]+2*(nx+1)
ii[5] = ii[1]+2*(nx+1)
ii[6] = ii[2]+2*(nx+1)
ii[7] = ii[3]+2*(nx+1)
f.write('%d %d\n'%(nx,(nx+1)*4))
for i in range(nx):
  f.write('H')
  for j in range(8):
    f.write(' %d'%(ii[j]+i))
  f.write('\n')
for k in range(2):
  for j in range(2):
    for i in range(nx+1):
      f.write('%f %f %f\n'%(i*dx,j*dy-dy/2.,k*dz-dz/2.))
f.close()
