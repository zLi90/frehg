"""
    An example illustrates the basic usage of Modflow6(flopy)
     - ZhiLi20191004
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import flopy
from bathymetry import *

'''
    FLAGS
'''
runSim = 1
makePlot = 1
saveData = 0
savename = 'flopy_out_unsat.csv'

'''
    ======== Problem definition ========
'''
workspace = os.path.join('data')
if not os.path.exists(workspace):
    os.makedirs(workspace)
name = 'mf6base'
# Dimension of domain
Lx = 110.0
Ly = 800.0
H = 10.0
Nx = 11
Ny = 80
Nlay = 21
# Hydraulic properties
Kx = 1e-5
Ky = 1e-5
Kz = 1e-6
porosity = 0.4
compW = 5e-10
compR = 1e-7
Ss = 9.8e3 * (compR + porosity*compW)
# Initial conditions
H0 = 0.0
# Boundary conditions
hTop = 0.0
hSide = 0.0
# Simulation Control (in sec)
Nt = 288
dt = 300.0
# Use evaporation, extinction depth
useEvap = 1
dext = 10.0
q_et = 0.0

'''
    ======== Create mf6 objects ========
'''
if runSim == 1:
    # Create flopy object
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name='mf6', version='mf6', sim_ws=workspace)
    # Create time discretization (just one step equilibrium)
    tdis = flopy.mf6.modflow.mftdis.ModflowTdis(sim, pname='tdis', time_units='seconds',
                                    nper=1, perioddata=[(Nt*dt, Nt, 1.0)])
    # Create groundwater flow object
    gwf = flopy.mf6.modflow.ModflowGwf(sim, modelname=name, model_nam_file='{}.nam'.format(name))
    # Create iterative solver
    ims = flopy.mf6.modflow.mfims.ModflowIms(sim, pname='ims', complexity='MODERATE')
    # Create spatial discretization
    dx = Lx/(Nx-1)
    dy = Ly/(Ny-1)
    bath = createBath(Nx, Ny)
    depthH = computeDepth(bath, hTop)
    zbot, ztop, actv, isurf = createDis(bath, depthH, Nx, Ny, Nlay, H)
    dis = flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(gwf, pname='dis', length_units='METERS',
                                    nlay=Nlay, nrow=Nx, ncol=Ny, delr=dx, delc=dy,
                                    top=ztop, botm=zbot, idomain = actv)
    # Create storage package
    sto = flopy.mf6.modflow.mfgwfsto.ModflowGwfsto(gwf, pname='sto', iconvert=1, transient=True,
                                    storagecoefficient=False, ss=Ss)
    # Create initial condition
    hic = H0 * np.ones((Nlay, Nx, Ny))
    ic = flopy.mf6.modflow.mfgwfic.ModflowGwfic(gwf, pname='ic', strt=hic)
    # Create nodal properties
    npf = flopy.mf6.modflow.mfgwfnpf.ModflowGwfnpf(gwf, pname='npf', icelltype=1,
                                    k=Kx, k22=Ky, k33=Kz, wetdry=1, save_flows=True)
    # Create boundary conditions
    chd_rec = []
    for ii in range(0,Nx):
        for jj in range(0,Ny):
            if actv[isurf[ii,jj],ii,jj] == 1:
                if depthH[ii,jj] > 0.001:
                    chd_rec.append(((isurf[ii,jj], ii, jj), depthH[ii,jj]))
#                    chd_rec.append(((isurf[ii,jj], ii, jj), 1.0))

    chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(gwf, pname='chd', maxbound=len(chd_rec),
                                    stress_period_data=chd_rec, save_flows=True)

    # Create unsaturated zone
    uzf_rec = []
    uzf_pkg = []
    iuzno = 0
    for ii in range(0,Nx):
        for jj in range(0,Ny):
            for kk in range(0,Nlay):
                if zbot[kk,ii,jj] > -dext and actv[kk,ii,jj] == 1:
                    if kk >= 2 and actv[kk-2,ii,jj] == 0 and depthH[ii,jj] <= 0.0:
                        land = 1
                    else:
                        land = 0
                    uzf_pkg.append((iuzno, (kk,ii,jj), land, 1, 0.1, Kz, 0.1, porosity, porosity, 5.0))
                    uzf_rec.append((iuzno, 0.0, q_et, dext, 0.1, 0.1, 1, 1))
                    iuzno += 1
    if useEvap == 1:
        uzf = flopy.mf6.modflow.mfgwfuzf.ModflowGwfuzf(gwf, pname='uzf', save_flows=True, nuzfcells=iuzno,                                            simulate_et=False, square_gwet=False, unsat_etwc=False, ntrailwaves=20, nwavesets=80, simulate_gwseep=False, packagedata=uzf_pkg, perioddata=uzf_rec)
        
    # Create output control
    head_filerecord = ['{}.hds'.format(name)]
    budget_filerecord = ['{}.cbb'.format(name)]
    saverecord = [('HEAD', 'ALL'),
                  ('BUDGET', 'ALL')]
    printrecord = [('HEAD', 'LAST')]
    oc = flopy.mf6.modflow.mfgwfoc.ModflowGwfoc(gwf, pname='oc', saverecord=saverecord,
                                                head_filerecord=head_filerecord,
                                                budget_filerecord=budget_filerecord,
                                                printrecord=printrecord)
'''
    ======== Execute simulation ========
'''
if runSim == 1:
    sim.write_simulation()
    success, buff = sim.run_simulation()
    print('\nSuccess is: ', success)
'''
    ======== Post-processing ========
'''
if makePlot == 1:
#    bath = createBath(Nx, Ny)
#    FrehdBath(bath, Nx, Ny)
    fname = os.path.join(workspace, '{}.hds'.format(name))
    hds = flopy.utils.binaryfile.HeadFile(fname)
    h = hds.get_data(kstpkper=(Nt-1,0))
    x = np.linspace(0, Lx, Nx)
    y = np.linspace(0, Ly, Ny)
    z = np.linspace(-H/Nlay/2, -H+H/Nlay/2, Nlay)
    print(np.shape(h))
    
    if saveData == 1:
        np.savetxt(savename,h[:,round((Nx-1)/2),:],delimiter=',')

    rg = np.linspace(0,2,10)
    # fig = plt.contourf(x,z, np.transpose(h[:,round((Nx-1)/2),:]),rg)
    fig = plt.imshow(h[:,round((Nx-1)/2),:],vmin=-1.0,vmax=2.0,cmap='jet')
    # fig = plt.imshow(h[:,round((Nx-1)/2),:])
    plt.colorbar()
    plt.show()
