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

dx1 = dx
dx2 = dx
dy1 = dy
dy2 = dy

print('dx1: %f\n'%dx1)
print('dx2: %f\n'%dx2)
print('dy1: %f\n'%dy1)
print('dy2: %f\n'%dy2)

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
f.write('%d %d %f 0. 0. %f\n'%(1,2,x_offset,dy*dz))
x_offset += dx
f.write('%d %d %f 0. 0. %f\n'%(2,3,x_offset,dy*dz))
for i in range(2,nx-2):
  x_offset += dxx[i%2]
  f.write('%d %d %f 0. 0. %f\n'%(i+1,i+2,x_offset,dy*dz))
x_offset += dx
f.write('%d %d %f 0. 0. %f\n'%(nx-1,nx,x_offset,dy*dz))

f.close()
print('done')

'''
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
'''

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

filename = '1d_uge.h5'
h5file = File(filename,mode='w')


vertices = numpy.zeros((nx*8,3),'=f8')
cells = numpy.zeros((nx*9),'=i4')
cell_ids = numpy.arange(1,nx+1,1,'=i4')
volumes = numpy.zeros(nx,'=f8')
xc = numpy.zeros(nx,'=f8')
yc = numpy.zeros(nx,'=f8')
zc = numpy.zeros(nx,'=f8')
zc[:] = dz/2.

a = numpy.arange(8)
for i in range(nx):
  cells[i*9] = 9
  cells[i*9+1:i*9+9] = a + i*8

icount = 0
for i in range(2):
  xc[i] = i*dx+dx/2.
  volumes[i] = dx*dy*dz
  vertices[icount][0] = i*dx
  vertices[icount][1] = -dy/2.
  vertices[icount][2] = 0.
  icount += 1
  vertices[icount][0] = (i+1)*dx
  vertices[icount][1] = -dy/2.
  vertices[icount][2] = 0.
  icount += 1
  vertices[icount][0] = (i+1)*dx
  vertices[icount][1] = dy/2
  vertices[icount][2] = 0.
  icount += 1
  vertices[icount][0] = i*dx
  vertices[icount][1] = dy/2
  vertices[icount][2] = 0.
  icount += 1
  vertices[icount][0] = i*dx
  vertices[icount][1] = -dy/2.
  vertices[icount][2] = dz
  icount += 1
  vertices[icount][0] = (i+1)*dx
  vertices[icount][1] = -dy/2.
  vertices[icount][2] = dz
  icount += 1
  vertices[icount][0] = (i+1)*dx
  vertices[icount][1] = dy/2.
  vertices[icount][2] = dz
  icount += 1
  vertices[icount][0] = i*dx
  vertices[icount][1] = dy/2.
  vertices[icount][2] = dz
  icount += 1

offset = dx*2
for i in range(2,nx-2):
  dxxx = dxx[i%2]
  dyyy = dyy[i%2]
  xc[i+2] = offset + dxxx/2.
  volumes[i] = dxxx*dyyy*dz
  vertices[icount][0] = offset
  vertices[icount][1] = -dyyy/2.
  vertices[icount][2] = 0
  icount += 1
  vertices[icount][0] = offset + dxxx
  vertices[icount][1] = -dyyy/2.
  vertices[icount][2] = 0.
  icount += 1
  vertices[icount][0] = offset + dxxx
  vertices[icount][1] = dyyy/2.
  vertices[icount][2] = 0.
  icount += 1
  vertices[icount][0] = offset
  vertices[icount][1] = dyyy/2.
  vertices[icount][2] = 0.
  icount += 1
  vertices[icount][0] = offset
  vertices[icount][1] = -dyyy/2.
  vertices[icount][2] = dz
  icount += 1
  vertices[icount][0] = offset + dxxx
  vertices[icount][1] = -dyyy/2.
  vertices[icount][2] = dz
  icount += 1
  vertices[icount][0] = offset + dxxx
  vertices[icount][1] = dyyy/2.
  vertices[icount][2] = dz
  icount += 1
  vertices[icount][0] = offset
  vertices[icount][1] = dyyy/2
  vertices[icount][2] = dz
  icount += 1
  offset += dxx[i%2]

for i in range(2):
  xc[i+nx-2] = offset + i*dx+dx/2.
  volumes[i+nx-2] = dx*dy*dz
  vertices[icount][0] = offset + i*dx
  vertices[icount][1] = -dy/2.
  vertices[icount][2] = 0.
  icount += 1
  vertices[icount][0] = offset + (i+1)*dx
  vertices[icount][1] = -dy/2.
  vertices[icount][2] = 0.
  icount += 1
  vertices[icount][0] = offset + (i+1)*dx
  vertices[icount][1] = dy/2
  vertices[icount][2] = 0.
  icount += 1
  vertices[icount][0] = offset + i*dx
  vertices[icount][1] = dy/2
  vertices[icount][2] = 0.
  icount += 1
  vertices[icount][0] = offset + i*dx
  vertices[icount][1] = -dy/2.
  vertices[icount][2] = dz
  icount += 1
  vertices[icount][0] = offset + (i+1)*dx
  vertices[icount][1] = -dy/2.
  vertices[icount][2] = dz
  icount += 1
  vertices[icount][0] = offset + (i+1)*dx
  vertices[icount][1] = dy/2.
  vertices[icount][2] = dz
  icount += 1
  vertices[icount][0] = offset + i*dx
  vertices[icount][1] = dy/2.
  vertices[icount][2] = dz
  icount += 1

dataset_name = 'Domain/Cell_Ids'
h5dset = h5file.create_dataset(dataset_name, data=cell_ids)
dataset_name = 'Domain/Cells'
h5dset = h5file.create_dataset(dataset_name, data=cells)
dataset_name = 'Domain/Vertices'
h5dset = h5file.create_dataset(dataset_name, data=vertices)
dataset_name = 'Domain/Volumes'
h5dset = h5file.create_dataset(dataset_name, data=volumes)
dataset_name = 'Domain/XC'
h5dset = h5file.create_dataset(dataset_name, data=xc)
dataset_name = 'Domain/YC'
h5dset = h5file.create_dataset(dataset_name, data=yc)
dataset_name = 'Domain/ZC'
h5dset = h5file.create_dataset(dataset_name, data=zc)

h5file.close()
