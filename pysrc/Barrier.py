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
from scipy import stats

#****************MAIN PROGRAM*****************
# First, read the .cfg file
configfile = sys.argv[1]
run = int(sys.argv[2])
ConfigData = ReadConfigFile(configfile)
outputfilebase = ConfigData["outputfilebase"]
outputfiledir = ConfigData["outputfiledir"]
recomb_factor = 1.0
file = open(outputfiledir+'/barrier.txt', 'w')
file.write("   Flux      Barrier      \n")

# This holds all of the data
dat = Array3dHDF5(outputfiledir, outputfilebase, run)

ScaleFactor = ConfigData["ScaleFactor"]
GridsPerPixel = ConfigData["GridsPerPixelX"]
Vph = float(ConfigData["Vparallel_hi"])    
Vpl = float(ConfigData["Vparallel_lo"])    
qfh = float(ConfigData["qfh"])    
dn = GridsPerPixel*ScaleFactor
nxx = dat.nx - 1
nyy = dat.ny - 1
nzz = dat.nz - 1

nxc = int(nxx/2)
nyc = int(nyy/2)

kmin = 4

centers = [(nxc-3*dn,nyc+3*dn),(nxc,nyc+3*dn),(nxc+3*dn,nyc+3*dn),(nxc-3*dn,nyc-3*dn),(nxc,nyc-3*dn),(nxc+3*dn,nyc-3*dn)]
fluxes = []
for kk in range(6):
    fluxes.append(ChargeDepth(dat,None,centers[kk][0], centers[kk][1], int(dn/2), int(dn/2), dat.Elec.shape[2], recomb_factor=recomb_factor))
bheights = []
linefluxes = []
linebheights = []

for m, (nxcenter,nycenter) in enumerate(centers):
    elecmax = 0
    for i in range(dat.Elec.shape[2]):
        if dat.Elec[nxcenter,nycenter,i] > 0.0 and i > elecmax: elecmax = i

    bheight = 20.0
    for k in range(4, elecmax):
        if dat.Elec[nxcenter, nycenter,k] < 1.0:
            continue

        beginphi = 0.0
        minphi = 20.0
        endphi = 0.0
        ReachedEnd = False
        for i in range(nycenter,nycenter+dn):
            if dat.Elec[nxcenter,i,k] < 1.0 and ReachedEnd:
                if dat.phi[nxcenter,i,k] < minphi:
                    minphi = dat.phi[nxcenter,i,k]
                    imin = i
                if dat.Elec[nxcenter,i+1,k] > 1.0:
                    endphi = dat.phi[nxcenter,i+1,k]
                    iend = i
                    break
            if dat.Elec[nxcenter,i,k] > 1.0 and not ReachedEnd:
                if dat.Elec[nxcenter,i+1,k] < 1.0:
                    beginphi = dat.phi[nxcenter,i,k]
                    ibegin = i
                    ReachedEnd = True
        if ReachedEnd == False: bheight = 0.0
        bheight = min(bheight,min(abs(beginphi - minphi),abs(endphi-minphi)))
        #print nxcenter, nycenter, beginphi, minphi,endphi, ibegin,imin,iend
        beginphi = 0.0
        minphi = 20.0
        endphi = 0.0
        ReachedEnd = False
        for i in range(nycenter-dn,nycenter):
            if dat.Elec[nxcenter,i,k] < 1.0 and ReachedEnd:
                if dat.phi[nxcenter,i,k] < minphi:
                    minphi = dat.phi[nxcenter,i,k]
                if dat.Elec[nxcenter,i+1,k] > 1.0:
                    endphi = dat.phi[nxcenter,i+1,k]
                    break
            if dat.Elec[nxcenter,i,k] > 1.0 and not ReachedEnd:
                if dat.Elec[nxcenter,i+1,k] < 1.0:
                    beginphi = dat.phi[nxcenter,i,k]
                    ReachedEnd = True
        if ReachedEnd == False: bheight = 0.0
        bheight = min(bheight,min(abs(beginphi - minphi),abs(endphi-minphi)))
        #print nxcenter, nycenter, beginphi, minphi,endphi, ibegin,imin,iend
    bheights.append(bheight)
    if bheight > 0:
        linebheights.append(bheight * 1000.0)
        linefluxes.append(fluxes[m])
slope, intercept, r_value, p_value, std_err = stats.linregress(linefluxes[-3:], linebheights[-3:])
value = (1100.0 - intercept) / slope

for kk in range(6):
    file.write(" %.1f     %.4f\n"%(fluxes[kk], bheights[kk]))
file.write("Flux at 1.1V = %f\n"%value)
file.close()

# Create the output directory if it doesn't exist
if not os.path.isdir(outputfiledir+"/plots"):
    os.mkdir(outputfiledir+"/plots")

plt.figure()
plt.title("Interwell Barrier Height vs Pixel Charge")
plt.plot(fluxes, np.array(bheights)*1000.0, marker = '*', label = "Vph = %.1f"%Vph)
plt.plot([0,300000],[1100.0, 1100.0], color = 'blue', ls = '--')
plt.text(200000, 1100.0,"Limit")
plt.xlabel("Pixel Charge (e-)")
plt.ylabel("Barrier Height (mV)")
plt.legend()
plt.savefig(outputfiledir+"/plots/Barrier_Height_Vpl_%.1f_Vph_%.1f.pdf"%(Vpl,Vph))

