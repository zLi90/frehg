"""
    Plot Frehg output
"""
import numpy as np
import matplotlib.pyplot as plt

Nx = 8
Ny = 150
Nz = 112

dx = 5.0
dy = 5.0
dz = 0.05

dt = 5.0
itvl = 36
Tend = 3024

wcs = 0.4
wcr = 0.08

domain = 'subsurf'
fdir = 'outputP5'

fhead = fdir+'/head_'
fsurf = fdir+'/surfaceZ_'
fdept = fdir+'/depth_'
fsatu = fdir+'/satu_'
fvelo = fdir+'/vVelocity_'
ind = [0, 2340, 3024]

def extractData3(fdir, ind):
    Nt = len(ind)
    data = np.zeros((Nx,Ny,Nz,Nt))
    for tt in range(0,Nt):
        fname = fdir+str(ind[tt])+'.dat'
        fid = open(fname,"r")
        for ii in range(0,Ny):
            for jj in range(0,Nx):
                for kk in range(0,Nz):
                    h = float(fid.readline())
                    data[jj,ii,kk,tt] = h
        fid.close()
    return data

def extractData2(fdir, ind):
    Nt = len(ind)
    data = np.zeros((Nx,Ny,Nt))
    for tt in range(0,Nt):
        fname = fdir+str(ind[tt])+'.dat'
        fid = open(fname,"r")
        for ii in range(0,Ny):
            for jj in range(0,Nx):
                h = float(fid.readline())
                data[jj,ii,tt] = h
        fid.close()
    return data

def extractVelo(fdir, fdep, itvl, Tend):
    tVec = np.linspace(0, Tend, round(Tend/itvl)+1)
    data1D = np.zeros((len(tVec)))
    for tt in range(len(tVec)):
        data = np.zeros((Nx,Ny))
        fname = fdir+str(int(tVec[tt]))+'.dat'
        fdept = fdep+str(int(tVec[tt]))+'.dat'
        fid1 = open(fname, 'r')
        fid2 = open(fdept, 'r')
        for ii in range(0,Ny):
            for jj in range(0,Nx):
                v = float(fid1.readline())
                d = float(fid2.readline())
                data[jj,ii] = v * d
        data1D[tt] = np.mean(data[:,79])
        fid1.close()
        fid2.close()
    return data1D

def extractWT(fdir, itvl, Tend):
    eps = 1e-2
    tVec = np.linspace(0, Tend, round(Tend/itvl)+1)
    data1D = np.zeros((len(tVec)))
    for tt in range(len(tVec)):
        data3D = np.zeros((Nx,Ny,Nz))
        fname = fdir+str(int(tVec[tt]))+'.dat'
        fid = open(fname, 'r')
        # read 3D data
        for ii in range(0,Ny):
            for jj in range(0,Nx):
                for kk in range(0,Nz):
                    data3D[jj,ii,kk] = float(fid.readline())
        fid.close()
        # extract a 2D slice
        slice = data3D[2,:,:]
        # search for wc < wcs
        flag = 0
        for ii in range(0,Ny):
            for kk in range(0,Nz):
                if kk == 0:
                    if slice[ii,kk] > wcs-eps:
                        data1D[tt] = (ii+1)*dy
                        flag = 1
                        break
                else:
                    if slice[ii,kk] > wcs-eps and slice[ii,kk-1] == wcr:
                        data1D[tt] = (ii+1)*dy
                        flag = 1
                        break
            if ii > 80:
                data1D[tt] = 400
                break
            if flag == 1:
                break
        # print('T = ',tVec[tt],', --- data1D = ',data1D[tt])
    return data1D

if domain == 'subsurf':
    head = extractData3(fhead, ind)
    satu = extractData3(fsatu, ind)
surf = extractData2(fsurf, ind)
dept = extractData2(fdept, ind)
velo = extractData2(fvelo, ind)
bath = surf - dept
flow = extractVelo(fvelo, fdept, itvl, Tend) * dx * Nx * 60.0
dist = extractWT(fsatu, itvl, Tend)

# plot surface domain
if domain == 'subsurf':
    ifig = 1
    for ff in range(0,len(ind)):
        plt.subplot(2,len(ind),ifig)
        slice = head[2,0:80,:,ff]
        fig = plt.imshow(np.transpose(slice),cmap='jet',vmin=-5,vmax=5)
        # plt.colorbar()
        ifig += 1
    for ff in range(0,len(ind)):
        plt.subplot(2,len(ind),ifig)
        slice = satu[2,0:80,:,ff]
        fig = plt.imshow(np.transpose(slice),cmap='jet',vmin=0.3,vmax=0.401)
        ifig += 1

# plot subsurface domain
ifig = 1
plt.figure(2)
for ff in range(0,len(ind)):
    plt.subplot(len(ind)*2,1,ifig)
    slice = dept[:,:,ff]
    fig = plt.imshow(slice,cmap='jet',vmin=0.0,vmax=0.04)
    plt.colorbar()
    ifig += 1

    plt.subplot(len(ind)*2,1,ifig)
    slice = velo[:,:,ff]
    fig = plt.imshow(slice,cmap='jet',vmin=-0.1,vmax=0.4)
    plt.colorbar()
    ifig += 1

# check mass conservation
outflow = sum(flow) * dt * itvl / 60
inflow = 5.5e-6 * 200 * 60.0 * Nx*dx * (Ny/2.0)*dy
print('Rain fall = ',inflow, 'Runoff = ',outflow)

tVec = np.linspace(0,int(Tend*dt),int(Tend/itvl+1))
plt.figure(3)
plt.plot(tVec/60.0, flow)
plt.ylim(0.0,10.0)
plt.xlim(0.0,Tend*dt/60.0)
plt.grid()

plt.figure(4)
plt.plot(tVec/60.0, dist)
# plt.ylim(0.0,10.0)
plt.xlim(0.0,Tend*dt/60.0)
plt.grid()

plt.show()
