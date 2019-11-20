#!/usr/bin/env python

#Author: Craig Lage, UC Davis;
#Date: 04-Nov-19
#This program plots the Poisson equation solutions from the C++ Poisson solver

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os, sys, time, h5py
sys.path.append(os.path.realpath('./pysrc'))
from pysubs import *  # These are the plotting subroutines

#****************MAIN PROGRAM*****************

# First, read the .cfg file

configfile = sys.argv[1]
run = int(sys.argv[2])
ConfigData = ReadConfigFile(configfile)

outputfilebase = ConfigData["outputfilebase"]
outputfiledir = ConfigData["outputfiledir"]
ScaleFactor = ConfigData["ScaleFactor"]
GridsPerPixelX = ConfigData["GridsPerPixelX"]
# Create the output directory if it doesn't exist
if not os.path.isdir(outputfiledir+"/plots"):
    os.mkdir(outputfiledir+"/plots")


print("Making Potential and Charge Density Convergence plots\n")
plt.figure()

plt.suptitle("Convergence Tests")
plt.subplot(1,2,1)
plt.title("Phi-Collect Gate", fontsize=12)
plt.subplot(1,2,2)
plt.title("Phi-Collect Gate", fontsize=12)
zmax = 64 * ScaleFactor
colors = ['black','red','green','blue','magenta','cyan']

NumMultis = int(3 + np.log2(ScaleFactor))

for i in range(NumMultis):

    # This holds all of the data
    dat = Array3dHDF5(outputfiledir, outputfilebase, run, multi=i)
    nxx = dat.nx - 1
    nyy = dat.ny - 1
    nzz = dat.nz - 1
    nxcenter = int(nxx/2)
    nycenter = int(nyy/2)
    plt.subplot(2,2,1)
    plt.plot(dat.z, dat.phi[nxcenter, nycenter,:], label = "Grid %d"%i, color = colors[i], marker = 'x')
    plt.ylim(0.0, 25.0)
    plt.xlim(0.0,5.0)
    plt.subplot(2,2,2)
    zm = int(zmax / pow(2,i))
    plt.plot(dat.z[0:zm], dat.phi[nxcenter, nycenter,0:zm], label = "Grid %d"%i, color = colors[i], marker = 'x')
    plt.ylim(-15.0, 25.0)
    nxcenter = int(nxcenter + ScaleFactor * GridsPerPixelX / 2 / pow(2,i))
    plt.subplot(2,2,3)
    plt.plot(dat.z, dat.phi[nxcenter, nycenter,:], label = "Grid %d"%i, color = colors[i], marker = 'x')
    plt.ylim(-60.0, 20.0)
    plt.subplot(2,2,4)
    plt.plot(dat.z[0:zm], dat.phi[nxcenter, nycenter,0:zm], label = "Grid %d"%i, color = colors[i], marker = 'x')
    plt.ylim(-6.0, 9.0)
plt.subplot(2,2,3)
plt.legend(loc = 'upper right')
plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_Grid_Comparison.pdf")


print("Making 1D Electron Density plots\n")
plt.figure()

plt.suptitle("Convergence Tests")
plt.subplot(1,1,1)
plt.title("Electron Density - Collect Gate", fontsize=12)
colors = ['black','red','green','blue','magenta','cyan']
for i in range(NumMultis):

    # This holds all of the data
    dat = Array3dHDF5(outputfiledir, outputfilebase, run, multi=i)
    nxx = dat.nx - 1
    nyy = dat.ny - 1
    nzz = dat.nz - 1
    zmax = dat.Elec.shape[2]
    nxcenter = int(nxx/2)
    nycenter = int(nyy/2)
    plt.plot(dat.z[0:zmax], dat.Elec[nxcenter, nycenter,:] / (dat.dx[0] * dat.dy[0]) / dat.dz[0:zmax], label = "Grid %d"%i, color = colors[i], marker = '*')
plt.legend(loc = 'upper right')
plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_Elec_Comparison.pdf")
