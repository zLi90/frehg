"""
Plot surface runoff modeled by Frehg (Maxwell_2014)
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

sim_groundwater = True
sim_shallowwater = True
validation = True
savefig = False


dim = [3, 10, 33]
delta = [80.0, 80.0, 0.1]
dt = 2.0
Tend = 18000
Titvl = 180
tVec = np.linspace(0,Tend,Tend/Titvl+1)

# directory of the model outputs
fdir = ['ex0_output/K-5T2/']
fsub = ['swe/','dif/']

# active domain dimension
Nx = 1
Ny = 5

# monitor location
xslice = 1
yslice = 5
tslice = 33

fdp = 'depth_'
fvv = 'vv_'

xVec = np.linspace(0,dim[0]*delta[0],dim[0])
yVec = np.linspace(0,dim[1]*delta[1],dim[1])
zVec = np.linspace(-5.0,0.1,dim[2])
co = ['k','r','b']
ls = ['-','--','-']
mk = ['o','d']
fs = 7
lw = 0.7
ind = ['()','(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)']
X,Y = np.meshgrid(xVec,yVec)

def extractModel(fname, dimension, tVec, domain):
    if domain == '3D':
        N = dim[0]*dim[1]*dim[2]
        Nt = len(tVec)
        data3d = np.zeros((dim[2],dim[0],dim[1],Nt))
        for tt in range(Nt):
            data = []
            fullname = fname+str(int(tVec[tt]))
            fid = open(fullname,'r')
            for ii in range(N):
                line = fid.readline()
                data.append(float(line))
            fid.close()
            data3d[:,:,:,tt] = np.transpose(np.reshape(np.array(data),(dim[1],dim[0],dim[2])))
    elif domain == '2D':
        N = dim[0]*dim[1]
        Nt = len(tVec)
        data3d = np.zeros((1,dim[0],dim[1],Nt))
        for tt in range(Nt):
            data = []
            fullname = fname+str(int(tVec[tt]))
            fid = open(fullname,'r')
            for ii in range(N):
                line = fid.readline()
                data.append(float(line))
            fid.close()
            data3d[:,:,:,tt] = np.transpose(np.reshape(np.array(data),(dim[1],dim[0],1)))
    return data3d

#
# Execution
#

# extract data
simu = {}
for ii in range(len(fdir)):
    simu[ii] = {}
    for jj in range(len(fsub)):
        simu[ii][jj] = {}
        simu[ii][jj]['dept'] = extractModel(fdir[ii]+fsub[jj]+fdp, dim, tVec, '2D')
        simu[ii][jj]['velo'] = extractModel(fdir[ii]+fsub[jj]+fvv, dim, tVec, '2D')
        simu[ii][jj]['flow'] = simu[ii][jj]['velo'] * simu[ii][jj]['dept'] * delta[0]
        simu[ii][jj]['runoff'] = np.squeeze((simu[ii][jj]['flow'][:,xslice,yslice,:]))

#
#   Make plots
#
font = {'family' : 'Arial',
        'size'   : fs}
matplotlib.rc('font', **font)
cm = 1.0 / 2.54

ii = 0
if sim_shallowwater:
    # runoff
    fig = plt.figure(1, figsize=[8*cm,6*cm])
    ax = fig.gca()
    pos1 = ax.get_position()
    pos2 = [pos1.x0 + 0.01, pos1.y0 + 0.03, pos1.width, pos1.height]
    ax.set_position(pos2)
    for jj in range(len(fsub)):
        plt.plot(tVec/60.0, simu[ii][jj]['runoff']*60.0, color=co[ii], linestyle=ls[jj], linewidth=lw)
    plt.xlabel(r'Time [min]', fontsize=fs)
    plt.ylim([0,12])
    plt.ylabel(r'Runoff [$\mathrm{m^3/min}$]', fontsize=fs)
    plt.legend(['SWE','DW'], fontsize=fs)

    if savefig:
        plt.savefig('Figure_2.eps',format='eps')

plt.show()
