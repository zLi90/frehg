""" Create input data for FREHG """
import numpy as np
import numpy.matlib

L = 100
H = 0.2
dim = [20, 200]

rBC = 5.5e-6
# rBC = 3e-6
tBC = 0.1
iBC = 0.0
fname = ['rain.dat','tideP.dat','inflow.dat']
fvalue = [rBC, tBC, iBC]

# time number for 2018-01-01
tBase = 1506837600.0
ndays = 100
dt = 3600

# bathymetry (P1, P2)
topo = 0.0 * np.ones((dim[1]))
for ii in range(L):
    topo[ii] = round(H - ii * H/L, 4) + 0.2
bath = np.matlib.repmat(topo, dim[0], 1)
bath1D = np.reshape(bath, (dim[0]*dim[1],1), order='F')
# bathymetry (P3)
# bath = np.zeros((162,200))
# bath[0:80,100:] = 65.0
# bath[82:,100:] = 65.0
# slope = np.linspace(0,40,81)
# topo = np.hstack((np.fliplr(np.reshape(slope,(1,81))),np.reshape(slope,(1,81))))
# for ii in range(100):
#     bath[:,ii] = topo + 20.0 - 0.2*ii
# bath1D = np.reshape(bath, (162*200,1), order='F')

# create time vector for rainfall, tide and inflow
tVec = np.zeros((int(ndays * 86400 / dt)))
for ii in range(len(tVec)):
    tVec[ii] = tBase + ii * dt

# create rainfall, tide and inflow
# write to file
for ii in range(len(fname)):
    BC = fvalue[ii] * np.ones((len(tVec)))
    data = np.zeros((len(tVec),2))
    data[:,0] = tVec
    data[:,1] = BC
    np.savetxt(fname[ii], data, delimiter=' ')

np.savetxt('bath.dat', bath1D, delimiter=' ')
