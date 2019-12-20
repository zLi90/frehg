"""
    Plot Frehg output
"""
import numpy as np
import matplotlib.pyplot as plt

Nx = 11
Ny = 80
Nz = 41

dx = 10.0
dy = 10.0
dz = 1.0

domain = 'subsurf'
fdir = 'outputFlopy'

fhead = fdir+'/head_'
fsurf = fdir+'/surfaceZ_'
fdept = fdir+'/depth_'
fsatu = fdir+'/satu_'
ind = [0,720]

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


surf = extractData2(fsurf, ind)
bath = surf - extractData2(fdept, ind)
if domain == 'subsurf':
    head = extractData3(fhead, ind)
    satu = extractData3(fsatu, ind)

ifig = 1
if domain == 'subsurf':
    # slice = head[5,:,:,0]
    # mdfw = np.genfromtxt('flopy_out_sat.csv',delimiter=',')
#    diff = (mdfw - np.transpose(slice))
    # diff = np.transpose(slice)
    # fig = plt.imshow(diff,cmap='jet',vmin=-1.0,vmax=2.0)
    # plt.colorbar()
    for ff in range(0,len(ind)):
        plt.subplot(1,len(ind),ifig)
        slice = head[1,:,:,ff]
        fig = plt.imshow(np.transpose(slice),cmap='jet',vmin=-1.0,vmax=1.0)
        plt.colorbar()
        ifig += 1

elif domain == 'surf':
    for ff in range(0,len(ind)):
        plt.subplot(len(ind),1,ifig)
        slice = surf[:,:,ff]
        fig = plt.imshow(slice,cmap='jet',vmin=-0.1,vmax=0.1)
        plt.colorbar()
        ifig += 1


#slice = head[5,:,:,len(ind)-1]
#fig = plt.imshow(np.transpose(slice),cmap='jet',vmin=-1.0,vmax=2.0)
#plt.colorbar()
plt.show()
