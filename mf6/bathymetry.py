"""
    Create bathymetry ('top' and 'botm' variables of Modflow DIS package) for
    Modflow6-flopy.
    ZhiLi 20191007
"""
import numpy as np
import numpy.matlib
import math
import matplotlib.pyplot as plt

def createBath(Nx, Ny):
    # Create 1D bottom elevation
    slope1 = np.linspace(0.0, -2.0, Ny/2.0)
    slope1[0:round(Ny/7)] = 0.0
    slope2 = np.linspace(-2.0, 0.0, Ny/2.0)
    bath1d = np.zeros((1,Ny))
    bath1d[0,0:round(Ny/2.0)] = slope1
    bath1d[0,round(Ny/2.0):Ny] = slope2
    bath1d[0,round(Ny*6/7):Ny] = 0.0
    bath1d[0,round(Ny*4/7):round(Ny*5/7)] = -0.5
    # Convert to 2D bottom elevation
    bath = np.matlib.repmat(bath1d[0,:], Nx, 1)

    return bath

def FrehdBath(bath, Nx, Ny):
    bmax = np.amax(bath)
    bath[0,:] = bmax + 0.01
    bath[-1,:] = bmax + 0.01
    bath[:,0] = bmax + 0.01
    bath[:,-1] = bmax + 0.01
    aa = bath == 0.0
    bath[aa] = 0.01
    fid = open('bath.dat','w')
    for ii in range(0,Ny):
        for jj in range(0,Nx):
            fid.write(str(bath[jj,ii])+"\n")
    return 0

def computeDepth(bath, h):
    dim = np.shape(bath)
    depth = np.zeros((dim))
    for ii in range(0,dim[0]):
        for jj in range(0,dim[1]):
            if h > bath[ii,jj]:
                depth[ii,jj] = h - bath[ii,jj]
            else:
                depth[ii,jj] = 0.0
    return depth

def createDis(bath, depth, Nx, Ny, Nlay, H):
    # Calculate bottom of each cell
    dz0 = H / (Nlay - 1)
    zbot = np.zeros((Nlay, Nx, Ny))
    for ii in range(0,Nx):
        for jj in range(0,Ny):
            Z = np.amax(bath) + dz0
            for kk in range(0,Nlay):
                Z -= dz0
                zbot[kk,ii,jj] = Z
                if Z < bath[ii,jj] and (Z+dz0) > bath[ii,jj]:
                    zbot[kk-1,ii,jj] = bath[ii,jj]
    # Check active status of grid cells
    actv = np.zeros((Nlay, Nx, Ny), dtype=int)
    for ii in range(0,Nx):
        for jj in range(0,Ny):
            for kk in range(0,Nlay):
                if zbot[kk,ii,jj] >= bath[ii,jj]:
                    actv[kk,ii,jj] = 0
#                elif zbot[kk,ii,jj] == bath[ii,jj]:
#                    if depth[ii,jj] > 0:
#                        actv[kk,ii,jj] = 1
#                    else:
#                        actv[kk,ii,jj] = 0
#                elif ii == 0 or ii == Nx-1 or jj == 0 or jj == Ny-1 or kk == Nlay-1:
#                    actv[kk,ii,jj] = 0
                else:
                    actv[kk,ii,jj] = 1
    # Calculate top of each column
    ztop = (np.amax(bath) + dz0) * np.ones((Nx, Ny))
    # Calculate index of surface cells
    isurf = np.zeros((Nx, Ny), dtype=int)
    for ii in range(0,Nx):
        for jj in range(0,Ny):
            for kk in range(0,Nlay):
                if zbot[kk,ii,jj] == bath[ii,jj]:
                    isurf[ii,jj] = int(kk+1)

    return zbot, ztop, actv, isurf
