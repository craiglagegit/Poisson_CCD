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
ConfigData = ReadConfigFile(configfile)
outputfilebase = ConfigData["outputfilebase"]
outputfiledir = ConfigData["outputfiledir"]
Nx = ConfigData["PixelBoundaryNx"]
Ny = ConfigData["PixelBoundaryNy"]
XCenter = ConfigData["DeltaVPixelCoords_0_0"][0]
YCenter = ConfigData["DeltaVPixelCoords_0_0"][1]
PixelSizeX = ConfigData["PixelSizeX"]
PixelSizeY = ConfigData["PixelSizeY"]
NxCenter = int((XCenter - ConfigData["PixelBoundaryLowerLeft"][0]) / PixelSizeX)
NyCenter = int((YCenter - ConfigData["PixelBoundaryLowerLeft"][1]) / PixelSizeY)
Area_0 = PixelSizeX * PixelSizeY
AreaStep = ConfigData["PixelAreas"]
NumSteps = ConfigData["NumSteps"]
NumElec = ConfigData["NumElec"]

pixelX = int(sys.argv[2])
pixelY = int(sys.argv[3])

elecs = []
areas=[]
for run in range(0, NumSteps, AreaStep):
    filename = outputfiledir + '/' + outputfilebase +'_%d_Area'%run + '.dat'
    [area,sim] = ReadAreaFile(filename, Nx, Ny, NxCenter, NyCenter, Area_0)
    filename = outputfiledir + '/' + outputfilebase +'_%d_CC'%run + '.dat'
    elec = ReadCCFile(filename, Nx, Ny)
    elecs.append(elec[pixelX, pixelY])
    areas.append(area[pixelX, pixelY])
# Create the output directory if it doesn't exist
if not os.path.isdir(outputfiledir+"/plots"):
    os.mkdir(outputfiledir+"/plots")

plt.figure()
plt.subplot(1,1,1)
plt.title("Pixel Area: Pixel(%d, %d)"%(pixelX, pixelY))

filename = outputfiledir + '/' + outputfilebase +'_PixelArea_%d_%d'%(pixelX, pixelY) + '.dat'
plt.plot(elecs, areas, marker='x')
plt.xlabel("Stored electrons")
plt.ylabel("PixelArea")

plt.savefig(outputfiledir+"/plots/PixelAreas_%d_%d.pdf"%(pixelX, pixelY))
