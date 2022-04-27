#!/usr/bin/env python

#Author: Craig Lage, UC Davis;
#Date: 24-Oct-19
#This program plots the Poisson equation solutions from the C++ Poisson solver

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os, sys, time, h5py
sys.path.append(os.path.realpath('./pysrc'))
from pysubs import *  # These are the plotting subroutines
#***************SUBROUTINES*****************
def colorbar(mappable):
    last_axes = plt.gca()
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(mappable, cax=cax)
    plt.sca(last_axes)
    return cbar

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
PixelSizeX = ConfigData["PixelSizeX"]
PixelSizeY = ConfigData["PixelSizeY"]

SensorThickness = ConfigData["SensorThickness"]
NumElec = ConfigData["NumElec"]
ContactVoltage = ConfigData["Vcontact"]
DeltaV = ConfigData["DeltaV_0_0"]
CenterVoltage = ContactVoltage + DeltaV

nxx = dat.nx - 1
nyy = dat.ny - 1
nzz = dat.nz - 1

nxcenter = int(nxx/2)
nycenter = int(nyy/2)
nxcenter2 = int(nxcenter + GridsPerPixelX * ScaleFactor / 2)
nxcenter3 = int(nxcenter + GridsPerPixelX * ScaleFactor)

NumPixelsPlotted = 5
nycenter2 = int(nycenter + GridsPerPixelY * ScaleFactor / 2)
nycenter3 = int(nycenter + GridsPerPixelY * ScaleFactor)

nymin = int(nycenter - (NumPixelsPlotted * ScaleFactor * GridsPerPixelY)/2)
nymax = int(nycenter + (NumPixelsPlotted * ScaleFactor * GridsPerPixelY)/2)

nxmin = int(nxcenter - (NumPixelsPlotted * ScaleFactor * GridsPerPixelX)/2)
nxmax = int(nxcenter + (NumPixelsPlotted * ScaleFactor * GridsPerPixelX)/2)

nzmin = 0
nzmax = nzz

# Create the output directory if it doesn't exist
if not os.path.isdir(outputfiledir+"/plots"):
    os.mkdir(outputfiledir+"/plots")

mpl.rcParams['contour.negative_linestyle'] = 'solid'
mpl.rcParams.update({'font.size': 6})

# Next, plots of the pixel boundaries
print("Making pixel plots\n")
plt.figure()
mpl.rcParams['contour.negative_linestyle'] = 'solid'
#mpl.rcParams.update({'font.size': 18})

plt.suptitle("MIRI Charge Collection Plot Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)

filename = outputfiledir+"/"+outputfilebase+'_'+str(run)+"_CC.dat"
file = open(filename,"r")
lines = file.readlines()
file.close()
if len(lines) < 2:
    print("No data in CC file.  Quitting")
    sys.exit()

Nx = ConfigData["PixelBoundaryNx"]
Ny = ConfigData["PixelBoundaryNy"]
    
CC_array = np.zeros([Nx, Ny])
for i, line in enumerate(lines):
    if i == 0:
        continue
    items = line.split()
    CC_array[int(items[0]), int(items[1])] = items[2]

TotalElectrons = CC_array.sum()
LostElectrons = NumElec * (run + 1) - TotalElectrons
print("A total of %d electrons were collected, so %d electrons recombined"%(TotalElectrons,LostElectrons))    
img1 = plt.imshow(CC_array, interpolation='nearest')
colorbar(img1)    
    
plt.xlabel("X(Pixels)",fontsize = 18)
plt.ylabel("Y(Pixels)",fontsize = 18)
plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_CC_%d.pdf"%run)


plt.suptitle("CCD Pixel Plots. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 24)
plotcounter = 1
plt.subplots_adjust(hspace=0.3, wspace=0.1)

filename = outputfiledir+"/"+outputfilebase+'_'+str(run)+"_Pts.dat"
file = open(filename,"r")
lines = file.readlines()
file.close()
if len(lines) < 2:
    print("No data in Pts file.  Quitting")
    sys.exit()
redsx=[]
redsy=[]
blacksx=[]
blacksy=[]
plottedxin = -1000.0
plottedyin = -1000.0
lines.remove(lines[0])
for line in lines:
    values = line.split()
    phase = int(values[2])
    if phase == 0:
        xin = float(values[3])
        yin = float(values[4])
    elif phase == 4:
        xout = float(values[3])
        yout = float(values[4])
        if np.isnan(xout) or np.isnan(yout):
            print("xin = %.3f, yin = %.3f is a nan")
            continue
        pixxout = int((xout - ConfigData["PixelBoundaryLowerLeft"][0]) / PixelSizeX)
        pixyout = int((yout - ConfigData["PixelBoundaryLowerLeft"][1]) / PixelSizeY)
        if (pixxout + pixyout) % 2 == 0:
            redsx.append(xin)
            redsy.append(yin)
        else:
            blacksx.append(xin)
            blacksy.append(yin)
        continue
    else:
        continue

plt.subplot(1,1,1,aspect=1)
plt.title("Pixel Boundaries",fontsize = 12)
if ConfigData["PixelBoundaryTestType"] == 0:
    spotsize = 10.0 * ConfigData["PixelBoundaryStepSize"][0] * ConfigData["PixelBoundaryStepSize"][1]
else:
    spotsize = 0.1
plt.scatter(redsx,redsy,s=spotsize,color="red")
plt.scatter(blacksx,blacksy,s=spotsize,color="black")

plt.xlabel("X(microns)",fontsize = 18)
plt.ylabel("Y(microns)",fontsize = 18)
plt.xlim(ConfigData["PixelBoundaryLowerLeft"][0], ConfigData["PixelBoundaryUpperRight"][0])
plt.ylim(ConfigData["PixelBoundaryLowerLeft"][1], ConfigData["PixelBoundaryUpperRight"][1])


plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_Pixels_%d.pdf"%run)


if ConfigData["LogPixelPaths"] != 0 and run % ConfigData["LogPixelPaths"] == 0:
    # Last, plots of the electron paths
    print("Making array electron path plots\n")
    # Plotting the paths along a line through the center)
    yline = (ConfigData["PixelBoundaryLowerLeft"][1] + ConfigData["PixelBoundaryUpperRight"][1]) / 2.0
    xline = (ConfigData["PixelBoundaryLowerLeft"][0] + ConfigData["PixelBoundaryUpperRight"][0]) / 2.0

    vertical_zoom = 4
    plt.figure()
    plt.suptitle("Electron Path Plot - Vertical Zoom = %d"%vertical_zoom, fontsize = 24)
    plt.subplots_adjust(wspace=0.2)
    #print(xline, yline)
    for line in lines:
        values = line.split()
        phase = int(values[2])
        if phase == 0:
            xin = float(values[3])
            yin = float(values[4])
            #print(xin, yin)
            if (yin > yline - 0.2) and (yin < yline + 0.2):
                #print("YPlotting thisID")
                YPlotThisID = True
                xpaths=[]
                zxpaths=[]
            else:
                YPlotThisID = False
            if (xin > xline - 0.2) and (xin < xline + 0.2):
                #print("XPlotting this ID")
                XPlotThisID = True
                ypaths=[]
                zypaths=[]
            else:
                XPlotThisID = False
            continue
        if XPlotThisID or YPlotThisID:
            xout = float(values[3])
            yout = float(values[4])
            zout = float(values[5])
            if np.isnan(xout) or np.isnan(yout) or np.isnan(zout):
                continue
            if YPlotThisID:
                xpaths.append(xout)
                zxpaths.append(zout)
                if phase == 4:
                    pixxin = int((xin - ConfigData["PixelBoundaryLowerLeft"][0]) / PixelSizeX)                
                    if pixxin % 2 == 0:
                        color = "red"
                    else:
                        color = "black"
                    plt.subplot(1,2,1,aspect=vertical_zoom)
                    plt.plot(xpaths, zxpaths, color = color, linewidth = 0.1)

            if XPlotThisID:
                ypaths.append(yout)
                zypaths.append(zout)
                if phase == 4:
                    pixyin = int((yin - ConfigData["PixelBoundaryLowerLeft"][1]) / PixelSizeY)
                    if pixyin % 2 == 0:
                        color = "red"
                    else:
                        color = "black"
                    plt.subplot(1,2,2,aspect=vertical_zoom)                        
                    plt.plot(ypaths, zypaths, color = color, linewidth = 0.1)

    plt.subplot(1,2,1,aspect=vertical_zoom)
    plt.ylabel("Z(microns)")
    plt.xlabel("X (microns)")
    plt.ylim(0.0, SensorThickness)
    plt.xlim(ConfigData["PixelBoundaryLowerLeft"][0], ConfigData["PixelBoundaryUpperRight"][0])
    plt.subplot(1,2,2,aspect=vertical_zoom)
    plt.ylabel("Z(microns)")
    plt.xlabel("Y (microns)")
    plt.ylim(0.0,SensorThickness)
    plt.xlim(ConfigData["PixelBoundaryLowerLeft"][1], ConfigData["PixelBoundaryUpperRight"][1])
    plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_Paths_%d.pdf"%run)

