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
from matplotlib.colors import from_levels_and_colors

#****************MAIN PROGRAM*****************

# First, read the .cfg file

configfile = sys.argv[1]
run = int(sys.argv[2])
    
ConfigData = ReadConfigFile(configfile)
outputfilebase = ConfigData["outputfilebase"]
outputfiledir = ConfigData["outputfiledir"]
Vpl = ConfigData["Vparallel_lo"]
Vph = ConfigData["Vparallel_hi"]
# This holds all of the data
dat = Array3dHDF5(outputfiledir, outputfilebase, run)

ScaleFactor = ConfigData["ScaleFactor"]
GridsPerPixelX = ConfigData["GridsPerPixelX"]
GridsPerPixelY = ConfigData["GridsPerPixelY"]
PixelSizeX = ConfigData["PixelSizeX"]
PixelSizeY = ConfigData["PixelSizeY"]
ZMult = 4.0
#kmax = int(dat.nz*0.9)
kmax = 75 * ScaleFactor

nxx = dat.nx - 1
nyy = dat.ny - 1
nzz = dat.nz - 1
#nxcenter = nxx/2
#nycenter = nyy/2
#nxmin = 0
#nxmax = nxx
#nymin = 0
#nymax = nyy
nxcenter = 42 * ScaleFactor
nycenter = 36 * ScaleFactor
nxmin = 18 * ScaleFactor
nxmax = 66 * ScaleFactor
nymin = 12 * ScaleFactor
nymax = 60 * ScaleFactor
dnx = 2 * ScaleFactor
dny = 2 * ScaleFactor

dgatenx = int(ScaleFactor*GridsPerPixelX/PixelSizeX * 4.9 / 2)
dgateny = int(ScaleFactor*GridsPerPixelY/PixelSizeY * 27.0 / 2)

ChargeDepth(dat, outputfiledir+'/charge.txt', 42*ScaleFactor, 36*ScaleFactor, dgatenx, dgateny, 44)

carriers = ['Electron', 'Hole', 'Mobile', 'Fixed']
plotdatas = [dat.Elec, dat.Hole, dat.Hole-dat.Elec, dat.rho]
cmap0 = plt.cm.get_cmap("jet")
cmap1 = plt.cm.get_cmap("seismic")
cmaps = [cmap0, cmap0, cmap1, cmap1]
ForceZeros = [False, False, True, True]
xslicemins = [nxcenter-dnx, nxcenter-dnx, nxcenter-dnx, nxcenter-dnx]
xslicemaxs = [nxcenter+dnx, nxcenter+dnx, nxcenter+dnx, nxcenter+dnx]
# Create the output directory if it doesn't exist
if not os.path.isdir(outputfiledir+"/plots"):
    os.mkdir(outputfiledir+"/plots")

for i, plotdata in enumerate(plotdatas):

    fig = plt.figure(figsize = (12,12))
    plt.suptitle("%s Charge Distribution"%carriers[i], fontsize = 36)

    ax1=plt.axes([0.10,0.40,0.50,0.50],aspect=1)
    ax1.set_title("X-Y Slice")
    ax1.set_xticks([])
    ax1.set_yticks([])
    [plotarray, dxx, dyy, levels, my_cmap] = BuildPlotArray(dat, plotdata, 2, nxmin, nxmax, nymin, nymax, 0, kmax, ZMult, ForceZeros[i], Vph, Vpl, cmaps[i])
    ax1.contourf(dxx, dyy, plotarray, levels = levels, cmap = my_cmap, extend='both')

    ax2=plt.axes([0.10,0.20,0.50,0.20], aspect=ZMult)
    ax2.set_title("X-Z Slice")
    ax2.set_xticks([])
    ax2.set_yticks([0.0,1.0,2.0])
    ax2.set_ylabel("Z (Microns)")
    [plotarray, dxx, dyy, levels, my_cmap] = BuildPlotArray(dat, plotdata, 1, nxmin, nxmax, nymin, nymax, 0, kmax, ZMult, ForceZeros[i], Vph, Vpl, cmaps[i])
    ax2.contourf(dxx, dyy, plotarray, levels = levels, cmap = my_cmap, extend='both')

    ax3=plt.axes([0.10,0.10,0.50,0.10])
    ax3.set_title("X-Cut, Y = %.2f"%dat.y[nycenter])
    ax3.set_xlabel("X (Microns)")
    ax3.set_ylabel("Charge Density \n(arb. units)")
    [plotslice, xpoints] = BuildPlotSlice(dat, plotdata, 0, nxmin, nxmax, nycenter-1, nycenter+2, 0, kmax)
    ax3.plot(xpoints, plotslice)
    
    ax4=plt.axes([0.60,0.40,0.20,0.50], aspect=1.0/ZMult)
    ax4.set_title("Y-Z Slice")
    ax4.set_xticks([0.0,1.0,2.0])
    ax4.set_yticks([])
    ax4.set_xlabel("Z (Microns)")
    [plotarray, dxx, dyy, levels, my_cmap] = BuildPlotArray(dat, plotdata, 0, 0, kmax, nymin, nymax, xslicemins[i], xslicemaxs[i], ZMult, ForceZeros[i], Vph, Vpl, cmaps[i])
    ax4.contourf(dxx, dyy, plotarray, levels = levels, cmap = my_cmap, extend='both')        
    ax4.invert_xaxis()

    ax5=plt.axes([0.80,0.40,0.10,0.50])
    ax5.set_title("Y-Cut, X = %.2f"%dat.x[nxcenter])
    ax5.yaxis.tick_right()
    ax5.yaxis.set_label_position("right")
    for tick in ax5.get_xticklabels():
        tick.set_rotation(90)
    ax5.set_ylabel("Y (Microns)")
    ax5.set_xlabel("Charge Density \n(arb. units)")
    [plotslice, xpoints] = BuildPlotSlice(dat, plotdata, 1, nymin, nymax, nxcenter-1, nxcenter+2, 0, kmax)
    ax5.plot(plotslice, xpoints)
    ax5.set_xlim(ax5.get_xlim()[::-1])

    ax6=plt.axes([0.65,0.20,0.10,0.10])
    ax6.set_title("Z-Cut")
    ax6.set_ylabel("Log Charge Density \n(arb. units)")
    ax6.set_xticks([0.0,1.0,2.0])
    ax6.set_xlabel("Z ( Microns)")
    ax6.yaxis.set_label_position("right")
    ax6.set_ylim(-4.0, 1.0)
    ax6.text(3.0,-8.0,"Z-Axis has a %.1fX Scale Multiplier"%ZMult, fontsize = 16)
    if i == 2:
        [plotslicen, xpoints] = BuildPlotSlice(dat, plotdata, 2, 0, kmax, nxcenter-1, nxcenter+2, nymin, nymax)
        [plotslicep, xpoints] = BuildPlotSlice(dat, plotdata, 2, 0, kmax, nxcenter-1, nxcenter+2, nymin, nymax)        
        ax6.plot(xpoints, np.log10(plotslicen/plotslicen.min()+1.0E-5), color='blue')
        ax6.plot(xpoints, np.log10(plotslicep/plotslicep.max()+1.0E-5), color='red')        
        ax6.text(3.0,-9.0,"+ Charge in Red, - Charge in Blue", fontsize = 16)        
    elif i == 3:
        [plotslicep, xpoints] = BuildPlotSlice(dat, plotdata, 2, 0, kmax, nxcenter-1, nxcenter+2, nymin, nymax)
        [plotslicen, xpoints] = BuildPlotSlice(dat, plotdata, 2, 0, kmax, nxcenter-1, nxcenter+2, nymin, nymax)        
        ax6.plot(xpoints, np.log10(plotslicen/plotslicen.min()+1.0E-5), color='blue')
        ax6.plot(xpoints, np.log10(plotslicep/plotslicep.max()+1.0E-5), color='red')        
        ax6.text(3.0,-9.0,"+ Charge in Red, - Charge in Blue", fontsize = 16)        
    else:
        [plotslice, xpoints] = BuildPlotSlice(dat, plotdata, 2, 0, kmax, nxmin, nxmax, nymin, nymax)
        ax6.plot(xpoints, np.log10(plotslice/plotslice.max()+1.0E-5))
    ax6.set_xlim(ax6.get_xlim()[::-1])
    ax6.text(3.0,-10.0,"Oxide in yellow, Gates in green", fontsize = 16)        
    plt.savefig(outputfiledir+"/plots/%sDistribution_30Apr17_%d.pdf"%(carriers[i],run))
    plt.close(fig)
