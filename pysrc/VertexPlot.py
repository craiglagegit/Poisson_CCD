#!/usr/bin/env python

#Author: Craig Lage, UC Davis;
#Date: 24-Oct-19
#This program plots the Poisson equation solutions from the C++ Poisson solver

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os, sys, time
sys.path.append(os.path.realpath('./pysrc'))
from pysubs import *  # These are the plotting subroutines

#****************MAIN PROGRAM*****************

# First, read the .cfg file

configfile = sys.argv[1]
run = int(sys.argv[2])
ConfigData = ReadConfigFile(configfile)
outputfilebase = ConfigData["outputfilebase"]
outputfiledir = ConfigData["outputfiledir"]
Nx = ConfigData["PixelBoundaryNx"]
Ny = ConfigData["PixelBoundaryNy"]
XCenter = ConfigData["FilledPixelCoords_0_0"][0]
YCenter = ConfigData["FilledPixelCoords_0_0"][1]
PixelSizeX = ConfigData["PixelSizeX"]
PixelSizeY = ConfigData["PixelSizeY"]
GridsPerPixelY = ConfigData["GridsPerPixelY"]
ScaleFactor = ConfigData["ScaleFactor"]
NxCenter = int((XCenter - ConfigData["PixelBoundaryLowerLeft"][0]) / PixelSizeX)
NyCenter = int((YCenter - ConfigData["PixelBoundaryLowerLeft"][1]) / PixelSizeY)
Area_0 = 100.0
NumAngles = 4 * ConfigData["NumVertices"] + 4
NumElec = ConfigData["CollectedCharge_0_0"]

PlotDelta = int(sys.argv[3])

filename = outputfiledir + '/' + outputfilebase +'_%d_Area'%run + '.dat'

[area,sim] = ReadAreaFile(filename, Nx, Ny, NxCenter, NyCenter, Area_0)

# Create the output directory if it doesn't exist
if not os.path.isdir(outputfiledir+"/plots"):
    os.mkdir(outputfiledir+"/plots")

plt.figure()
plt.subplot(1,1,1,aspect = 1)
plt.title("Pixel Area: %d e-"%NumElec)

filename = outputfiledir + '/' + outputfilebase +'_%d_Vertices'%run + '.dat'

(vx, vy) = ReadVertexFile(filename, Nx, Ny, NumAngles)

LineXMin = XCenter - (PlotDelta + 0.5) * PixelSizeX
LineXMax = XCenter + (PlotDelta + 0.5) * PixelSizeX
LineYMin = YCenter - (PlotDelta + 0.5) * PixelSizeY
LineYMax = YCenter + (PlotDelta + 0.5) * PixelSizeY

# In the case of an even number of collecting phases,
# The pixel boundaries are shifted down by 1/2 of a grid cell.
if ConfigData["CollectingPhases"]%2 == 1:
    DeltaY = 0.0
else:
    DeltaY = PixelSizeY / float(GridsPerPixelY * ScaleFactor * 2) 

plt.title("Pixel Vertices: %d e-"%NumElec)

for i in range(NxCenter-PlotDelta, NxCenter+PlotDelta+1):
    plt.plot([XCenter+PixelSizeX*(i-NxCenter-0.5), XCenter+PixelSizeX*(i-NxCenter-0.5)], [LineYMin, LineYMax], color = 'black', ls = '--')
for j in range(NyCenter-PlotDelta, NyCenter+PlotDelta+1):
    plt.plot([LineXMin, LineXMax], [YCenter+PixelSizeY*(j-NyCenter-0.5)-DeltaY, YCenter+PixelSizeY*(j-NyCenter-0.5)-DeltaY], color = 'black', ls = '--')


for i in range(NxCenter-PlotDelta, NxCenter+PlotDelta+1):
    for j in range(NyCenter-PlotDelta, NyCenter+PlotDelta+1):
        if i == NxCenter and j == NyCenter:
            textcolor = 'red'
        else:
            textcolor = 'black'
        plt.text(XCenter+PixelSizeX*(i-NxCenter-0.2), YCenter+PixelSizeY*(j-NyCenter-0.1), "%.4f"%area[i,j], color = textcolor, fontsize = 12/PlotDelta, fontweight = 'bold')
        x = []
        y = []
        for k in range(NumAngles):
            x.append(vx[i,j,k])
            y.append(vy[i,j,k])
        x.append(vx[i,j,0])
        y.append(vy[i,j,0])
        plt.plot(x, y, lw = 0.5)

plt.savefig(outputfiledir+"/plots/PixelVertices_%d.pdf"%run)
