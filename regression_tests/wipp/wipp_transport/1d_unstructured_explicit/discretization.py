import sys
from h5py import *
import numpy

filename = 'grid.uge'
f = open(filename,'w')

length = 100.
nx = 1000
dx = length/float(nx)
dy = 1.
dz = 1.

dx1 = 0.4*dx
dx2 = 1.6*dx
dy1 = 2.5
dy2 = 0.625

dxx = [dx1,dx2]
dyy = [dy1,dy2]

# all cells have a volume of 1 m^3
# the first two and last two grid cells are 0.1 x 1 x 1 = 1 m^3 
# the cells in between alternate between:
# i % 2 == 0 :  0.04 x 2.500 x 1 = 1 m^3
# i % 2 == 1 :  0.16 x 0.625 x 1 = 1 m^3

f.write('CELLS %d\n'%nx)
x_offset = 0.
f.write('1 %f 0. 0. %f\n'%(x_offset+0.5*dx,dx*dy*dz))
x_offset += dx
f.write('2 %f 0. 0. %f\n'%(x_offset+0.5*dx,dx*dy*dz))
x_offset += dx

for i in range(1,nx/2-1):
  center1 = x_offset + 0.5 * dx1
  f.write('%d %f 0. 0. %f\n'%(i*2+1,center1,dx1*dy1*dz))
  x_offset +=  dx1
  center2 = x_offset + 0.5 * dx2
  f.write('%d %f 0. 0. %f\n'%(i*2+2,center2,dx2*dy2*dz))
  x_offset +=  dx2
  
f.write('%d %f 0. 0. %f\n'%(nx-1,x_offset+0.5*dx,dx*dy*dz))
x_offset += dx
f.write('%d %f 0. 0. %f\n'%(nx,x_offset+0.5*dx,dx*dy*dz))
x_offset += dx

f.write('CONNECTIONS %d\n'%(nx-1))
x_offset = dx
f.write('%d %d %f 0. 0. %f\n'%(1,2,x_offset,dz))
x_offset += dx
f.write('%d %d %f 0. 0. %f\n'%(2,3,x_offset,dz))
for i in range(2,nx-2):
  x_offset += dxx[i%2]
  f.write('%d %d %f 0. 0. %f\n'%(i+1,i+2,x_offset,dz))
x_offset += dx
f.write('%d %d %f 0. 0. %f\n'%(nx-1,nx,x_offset,dz))

f.close()
print('done')

filename = 'alpha.h5'
h5file = File(filename,mode='w')

iarray = numpy.arange(1,nx+1,1,'=i4')
dataset_name = 'Cell Ids'
h5dset = h5file.create_dataset(dataset_name, data=iarray)
rarray = numpy.ones(nx,'=f8')
for i in range(2,nx-2):
  rarray[i] = dyy[i%2]
dataset_name = 'Alpha'
h5dset = h5file.create_dataset(dataset_name, data=rarray)

h5file.close()

filename = 'west.ex'
f = open(filename,'w')
f.write('CONNECTIONS 1\n')
f.write('1 0. 0. 0. %f\n'%(dy*dz))
f.close()

filename = 'east.ex'
f = open(filename,'w')
f.write('CONNECTIONS 1\n')
f.write('%d %f 0. 0. %f\n'%(nx,length,dy*dz))
f.close()