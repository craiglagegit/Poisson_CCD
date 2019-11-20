#!/usr/bin/env python

#Author: Craig Lage, UC Davis;
#Date: 24-Oct-19
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

# This holds all of the data
dat = Array3dHDF5(outputfiledir, outputfilebase, run)

ScaleFactor = ConfigData["ScaleFactor"]
GridsPerPixelX = ConfigData["GridsPerPixelX"]
GridsPerPixelY = ConfigData["GridsPerPixelY"]
ChargeFactor = 1.6E-19 * 1.0E6 / (11.7 * 8.85E-12)/((dat.x[1]-dat.x[0])*(dat.y[1]-dat.y[0])) #(QE*MICRON_PER_M/(EPSILON_0*EPSILON_SI))/(dx*dy)


nxx = dat.nx - 1
nyy = dat.ny - 1
nzz = dat.nz - 1

nxcenter = int(nxx/2)
nycenter = int(nyy/2)
nxcenter2 = nxcenter

NumPixelsPlottedX = 3
NumPixelsPlottedY = 4    
nycenter2 = nycenter
nymin = int(nycenter - (NumPixelsPlottedY * ScaleFactor * GridsPerPixelY)/2)
nymax = int(nycenter + (NumPixelsPlottedY * ScaleFactor * GridsPerPixelY)/2)

nxmin = int(nxcenter - (NumPixelsPlottedX * ScaleFactor * GridsPerPixelX)/2)
nxmax = int(nxcenter + (NumPixelsPlottedX * ScaleFactor * GridsPerPixelX)/2)

nzmin = 0
nzmax = 16 * ScaleFactor
"""
print "Channel CG"
for i in range(dat.Elec.shape[2]):
    print "k = %d, z=%.3f, phi = %.3f, holes = %.3f, elec = %.3f, rho = %.3f, Ez = %.3f"%(i,dat.z[i],dat.phi[nxcenter, nycenter,i],dat.Hole[nxcenter, nycenter,i],dat.Elec[nxcenter, nycenter,i],dat.rho[nxcenter, nycenter,i],dat.Ez[nxcenter, nycenter,i])

print "Channel BG"
for i in range(dat.Elec.shape[2]):
    print "k = %d, z=%.3f, phi = %.3f, holes = %.3f, elec = %.3f, rho = %.3f, Ez = %.3f"%(i,dat.z[i],dat.phi[nxcenter, nycenter+ScaleFactor*GridsPerPixelY/2,i],dat.Hole[nxcenter, nycenter+ScaleFactor*GridsPerPixelY/2,i],dat.Elec[nxcenter, nycenter+ScaleFactor*GridsPerPixelY/2,i],dat.rho[nxcenter, nycenter+ScaleFactor*GridsPerPixelY/2,i],dat.Ez[nxcenter, nycenter+ScaleFactor*GridsPerPixelY/2,i])
"""

# Create the output directory if it doesn't exist
if not os.path.isdir(outputfiledir+"/plots"):
    os.mkdir(outputfiledir+"/plots")

mpl.rcParams['contour.negative_linestyle'] = 'solid'
mpl.rcParams.update({'font.size': 6})

print("Making 1D Potential and Charge Density plots\n")
plt.figure()

plt.suptitle("1D Potential and Charge Density Slices. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
plotcounter = 1
plt.subplots_adjust(hspace=0.3, wspace=0.3)
phinumzs = 160
numzs = 160
elecnumzs = ConfigData["Nzelec"] * ConfigData["ScaleFactor"]

centers = [(0,0), (-3,3), (0,3), (3,3), (-3,-3), (0,-3), (3,-3)]
fluxes = []
fluxes.append(0)
for kk in range(6):
    kkk = 3 * kk
    fluxes.append(ConfigData["CollectedCharge_0_%d"%kkk])

plt.subplot(2,2,1)
plt.title("Phi-Collect Gate", fontsize=12)

for m,(dnx, dny) in enumerate(centers):
    nxcenter2 = int(nxcenter + dnx * ScaleFactor * GridsPerPixelX)
    nycenter2 = int(nycenter + dny * ScaleFactor * GridsPerPixelY)
    plt.plot(dat.z[0:phinumzs],(dat.phi[nxcenter2,nycenter2,0:phinumzs]+dat.phi[nxcenter2-1,nycenter2,0:phinumzs]+dat.phi[nxcenter2,nycenter2-1,0:phinumzs]+dat.phi[nxcenter2-1,nycenter2-1,0:phinumzs])/4.0, label = "%d"%fluxes[m])
plt.legend(loc = "upper right")
#plt.xlabel("Z-Dimension (microns)")
plt.ylabel('$\phi(x,y,z)$ [V]',fontsize=12)
plt.ylim(-10.0, 30.0)
plt.xlim(0.0,4.0)

plt.subplot(2,2,2)
plt.title("Phi-Barrier Gate", fontsize=12)

for m,(dnx, dny) in enumerate(centers):
    nxcenter2 = int(nxcenter + dnx * ScaleFactor * GridsPerPixelX)
    nycenter2 = int(nycenter + dny * ScaleFactor * GridsPerPixelY + ScaleFactor * GridsPerPixelY / 2)
    plt.plot(dat.z[0:phinumzs],(dat.phi[nxcenter2,nycenter2,0:phinumzs]+dat.phi[nxcenter2-1,nycenter2,0:phinumzs]+dat.phi[nxcenter2,nycenter2-1,0:phinumzs]+dat.phi[nxcenter2-1,nycenter2-1,0:phinumzs])/4.0, label = "%d"%fluxes[m])
plt.legend(loc = "upper right")
#plt.xlabel("Z-Dimension (microns)")
plt.ylabel('$\phi(x,y,z)$ [V]',fontsize=12)
plt.ylim(-10.0, 30.0)
plt.xlim(0.0,4.0)

plt.subplot(2,2,3)
plt.title("Rho-Collect Gate", fontsize=12)
nxcenter2 = nxcenter
nycenter2 = nycenter
plt.plot(dat.z[0:numzs], (dat.rho[nxcenter2,nycenter2,0:numzs]+dat.rho[nxcenter2-1,nycenter2,0:numzs]+dat.rho[nxcenter2,nycenter2-1,0:numzs]+dat.rho[nxcenter2-1,nycenter2-1,0:numzs])/4.0, label = "Fixed charge", color='green')
for m,(dnx, dny) in enumerate(centers):
    nxcenter2 = int(nxcenter + dnx * ScaleFactor * GridsPerPixelX)
    nycenter2 = int(nycenter + dny * ScaleFactor * GridsPerPixelY)
    plt.plot(dat.z[0:elecnumzs], -ChargeFactor * (dat.Elec[nxcenter2,nycenter2,0:elecnumzs]+dat.Elec[nxcenter2-1,nycenter2,0:elecnumzs]+dat.Elec[nxcenter2,nycenter2-1,0:elecnumzs]+dat.Elec[nxcenter2-1,nycenter2-1,0:elecnumzs])/4.0 / dat.dz[0:elecnumzs], label = "%d"%fluxes[m])

plt.legend(loc = "upper right")
plt.xlabel("Z-Dimension (microns)", fontsize=12)
plt.ylabel('$\\rho(x,y,z)/\epsilon_{Si}$ [V/um$^2$]',fontsize=12)
plt.ylim(-80.0, 80.0)
plt.xlim(0.0,4.0)

slicez = 8 * ScaleFactor
ddny = int(3 * ScaleFactor * GridsPerPixelY / 2)
plt.subplot(2,2,4)
plt.title("Potential Barrier, z = %.2f"%dat.z[slicez])
for m,(dnx, dny) in enumerate(centers):
    nxcenter2 = int(nxcenter + dnx * ScaleFactor * GridsPerPixelX)
    nycenter2 = int(nycenter + dny * ScaleFactor * GridsPerPixelY)
    nymin = int(nycenter2 - ddny)
    nymax = int(nycenter2 + ddny)
    deltay = dat.y[nycenter2]
    plt.plot(dat.y[nymin:nymax]-deltay,dat.phi[nxcenter2,nymin:nymax, slicez], label = "%d"%fluxes[m])
plt.ylim(-5.0, 15.0)
plt.xlim(dat.y[nycenter-ddny]-dat.y[nycenter], dat.y[nycenter+ddny]-dat.y[nycenter])
plt.xlabel("Y-Dimension (microns)")
plt.ylabel("Potential(Volts)")
plt.legend(loc = "lower right")


plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_1D_Multi_%d.pdf"%run)

