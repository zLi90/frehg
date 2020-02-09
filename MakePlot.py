"""
    Plot Frehg output
"""
import numpy as np
import matplotlib.pyplot as plt

Nx = 20
Ny = 200
Nz = 26

dx = 4.0
dy = 4.0
dz = 0.2

dt = 4.0
itvl = 45
Tend = 4500

domain = 'subsurf'
fdir = 'outputP2'

fhead = fdir+'/head_'
fsurf = fdir+'/surfaceZ_'
fdept = fdir+'/depth_'
fsatu = fdir+'/satu_'
fvelo = fdir+'/vVelocity_'
ind = [180, 2970, 4500]

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
        data1D[tt] = np.mean(data[:,96])
        fid1.close()
        fid2.close()
    return data1D

if domain == 'subsurf':
    head = extractData3(fhead, ind)
    satu = extractData3(fsatu, ind)
surf = extractData2(fsurf, ind)
dept = extractData2(fdept, ind)
velo = extractData2(fvelo, ind)
bath = surf - dept
flow = extractVelo(fvelo, fdept, itvl, Tend) * dx * Nx * 60.0

# plot surface domain
if domain == 'subsurf':
    ifig = 1
    for ff in range(0,len(ind)):
        if Ny < Nz:
            plt.subplot(1,len(ind),ifig)
        else:
            plt.subplot(len(ind),2,ifig)
        slice = head[1,0:100,:,ff]
        fig = plt.imshow(np.transpose(slice),cmap='jet',vmin=-5,vmax=5)
        # plt.colorbar()
        ifig += 1

        plt.subplot(len(ind),2,ifig)
        slice = satu[1,0:100,:,ff]
        fig = plt.imshow(np.transpose(slice),cmap='jet',vmin=0,vmax=0.45)
        ifig += 1

# plot subsurface domain
ifig = 1
plt.figure(2)
for ff in range(0,len(ind)):
    plt.subplot(len(ind),2,ifig)
    slice = dept[:,:,ff]
    fig = plt.imshow(slice,cmap='jet',vmin=0.0,vmax=0.04)
    plt.colorbar()
    ifig += 1

    plt.subplot(len(ind),2,ifig)
    slice = velo[:,:,ff]
    fig = plt.imshow(slice,cmap='jet',vmin=-0.1,vmax=0.4)
    plt.colorbar()
    ifig += 1

# check mass conservation
outflow = sum(flow) * dt * itvl / 60
inflow = 5.5e-6 * 200 * 60.0 * Nx*dx * (Ny/2.0)*dy
# print(inflow, outflow)

tVec = np.linspace(0,int(Tend*dt),int(Tend/itvl+1))
plt.figure(3)
plt.plot(tVec/60.0, flow)
plt.ylim(0.0,12.0)
plt.xlim(0.0,Tend*dt/60.0)
plt.grid()
#slice = head[5,:,:,len(ind)-1]
#fig = plt.imshow(np.transpose(slice),cmap='jet',vmin=-1.0,vmax=2.0)
#plt.colorbar()
plt.show()
