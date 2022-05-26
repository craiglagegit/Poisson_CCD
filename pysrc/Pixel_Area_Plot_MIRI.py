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

run = 0
Area_0 = 625.0
DeltaVs = []
Areas = []
DeltaVs.append(0.0)
Areas.append(0.0)
for i in range(7):
    try:
        configfile = f"data/pixel_area_0p{i+1}/pixel.cfg"
        ConfigData = ReadConfigFile(configfile)
        outputfilebase = ConfigData["outputfilebase"]
        outputfiledir = ConfigData["outputfiledir"]
        Nx = ConfigData["PixelBoundaryNx"]
        Ny = ConfigData["PixelBoundaryNy"]
        XCenter = ConfigData["DeltaVPixelCoords_0_0"][0]
        YCenter = ConfigData["DeltaVPixelCoords_0_0"][1]
        PixelSizeX = ConfigData["PixelSizeX"]
        PixelSizeY = ConfigData["PixelSizeY"]
        GridsPerPixelY = ConfigData["GridsPerPixelY"]
        ScaleFactor = ConfigData["ScaleFactor"]
        NxCenter = int((XCenter - ConfigData["PixelBoundaryLowerLeft"][0]) / PixelSizeX)
        NyCenter = int((YCenter - ConfigData["PixelBoundaryLowerLeft"][1]) / PixelSizeY)
        DeltaV = ConfigData["DeltaV_0_0"]
        filename = outputfiledir + '/' + outputfilebase +'_%d_Area'%run + '.dat'
        [area,sim] = ReadAreaFile(filename, Nx, Ny, NxCenter, NyCenter, Area_0)
        DeltaVs.append(DeltaV)
        Areas.append((Area_0 - area[2,2]) / Area_0 * 100.0)
    except:
        continue
print(DeltaVs)
print(Areas)
    
# Create the output directory if it doesn't exist
if not os.path.isdir(outputfiledir+"/plots"):
    os.mkdir(outputfiledir+"/plots")

plt.figure()
plt.subplot(1,1,1)
plt.title("Pixel Area Loss vs DeltaV")

plt.plot(DeltaVs, Areas, marker = 'x')
plt.xlabel("Pixel DeltaV (V)")
plt.ylabel("Pixel Area Loss (%)")

print(outputfiledir)
plt.savefig(outputfiledir+"/plots/PixelAreas.pdf")
