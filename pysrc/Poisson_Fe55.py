#!/usr/bin/env python

#Author: Craig Lage, UC Davis;
#Date: 24-Oct-19

#This program plots the Poisson equation solutions from the C++ Poisson solver
import numpy as np
import matplotlib.pyplot as plt
import os, sys, time
sys.path.append(os.path.realpath('./pysrc'))
from pysubs import *  # These are the plotting subroutines

#****************MAIN PROGRAM*****************

# First, read the .cfg file

configfile = sys.argv[1]
run = int(sys.argv[2])
skip = int(sys.argv[3])
ConfigData = ReadConfigFile(configfile)

outputfilebase = ConfigData["outputfilebase"]
outputfiledir = ConfigData["outputfiledir"]
NumElec = ConfigData["NumElec"]
PixelSizeX = ConfigData["PixelSizeX"]
PixelSizeY = ConfigData["PixelSizeY"]

# Create the output directory if it doesn't exist
if not os.path.isdir(outputfiledir+"/plots"):
    os.mkdir(outputfiledir+"/plots")

filename = outputfiledir+"/"+outputfilebase+'_'+str(run)+"_Pts.dat"
file = open(filename,"r")
lines = file.readlines()
file.close()
if len(lines) < 2:
    print("No data in Pts file.  Quitting")
    sys.exit()
lines.remove(lines[0])

if ConfigData["LogPixelPaths"] != 0 and run % ConfigData["LogPixelPaths"] == 0:
    # Last, plots of the electron paths
    print("Making array electron path plots\n")
    # Plotting the paths along a line through the center

    NumSteps = 2000
    xpaths = np.zeros([NumElec, NumSteps])
    zxpaths = np.zeros([NumElec, NumSteps])
    ypaths = np.zeros([NumElec, NumSteps])
    zypaths = np.zeros([NumElec, NumSteps])
    maxsteps = np.zeros([NumElec], dtype=int)
    
    vertical_zoom = 0.5
    plt.figure()
    plt.suptitle("Fe55 Path Plot - Vertical Zoom = %.1f"%vertical_zoom, fontsize = 24)
    plt.subplots_adjust(wspace=0.2)

    for line in lines:
        values = line.split()
        id = int(values[0])
        step = int(values[1])
        phase = int(values[2])
        if step > NumSteps - 1:
            continue
        
        x = float(values[3])
        y = float(values[4])
        z = float(values[5])
        if np.isnan(x) or np.isnan(y) or np.isnan(z):
            continue
        xpaths[id, step] = x
        zxpaths[id, step] = z            
        ypaths[id, step] = y
        zypaths[id, step] = z
        maxsteps[id] = step
    fe_counter = 0
    for id in range(0, NumElec, skip):
        if maxsteps[id] == 1999:        
            fe_counter += 1
        else:
            plt.subplot(1,2,1,aspect=vertical_zoom)
            plt.plot(xpaths[id,0:maxsteps[id]], zxpaths[id,0:maxsteps[id]], color='blue', linewidth = 0.01)
            plt.subplot(1,2,2,aspect=vertical_zoom)
            plt.plot(ypaths[id,0:maxsteps[id]], zypaths[id,0:maxsteps[id]], color='blue', linewidth = 0.01)


    print("%d electrons did not finish"%fe_counter)

# Now plot the holes
        
filename = outputfiledir+"/"+outputfilebase+'_'+str(run)+"_HPts.dat"
file = open(filename,"r")
lines = file.readlines()
file.close()
if len(lines) < 2:
    print("No data in Pts file.  Quitting")
    sys.exit()
lines.remove(lines[0])

if ConfigData["LogPixelPaths"] != 0 and run % ConfigData["LogPixelPaths"] == 0:
    # Last, plots of the hole paths
    print("Making array hole path plots\n")
    # Plotting the paths along a line through the center

    NumSteps = 2000
    xpaths = np.zeros([NumElec, NumSteps])
    zxpaths = np.zeros([NumElec, NumSteps])
    ypaths = np.zeros([NumElec, NumSteps])
    zypaths = np.zeros([NumElec, NumSteps])
    maxsteps = np.zeros([NumElec], dtype=int)

    for line in lines:
        values = line.split()
        id = int(values[0])
        step = int(values[1])
        phase = int(values[2])
        if step > NumSteps - 1:
            continue
        
        x = float(values[3])
        y = float(values[4])
        z = float(values[5])
        if np.isnan(x) or np.isnan(y) or np.isnan(z):
            continue
        xpaths[id, step] = x
        zxpaths[id, step] = z            
        ypaths[id, step] = y
        zypaths[id, step] = z
        maxsteps[id] = step
    fh_counter = 1
    for id in range(0, NumElec, skip):
        if maxsteps[id] == 1999:
            fh_counter += 1
        else:
            plt.subplot(1,2,1,aspect=vertical_zoom)
            plt.plot(xpaths[id,0:maxsteps[id]], zxpaths[id,0:maxsteps[id]], color='red', linewidth = 0.01)
            plt.subplot(1,2,2,aspect=vertical_zoom)
            plt.plot(ypaths[id,0:maxsteps[id]], zypaths[id,0:maxsteps[id]], color='red', linewidth = 0.01)


    print("%d holes did not finish"%fh_counter)

        
    plt.subplot(1,2,1,aspect=vertical_zoom)
    plt.ylabel("Z(microns)")
    plt.xlabel("X (microns)")
    plt.xticks([20,30,40,50])
    plt.ylim(0,110)
    plt.xlim(ConfigData["PixelBoundaryLowerLeft"][0]+PixelSizeX, ConfigData["PixelBoundaryUpperRight"][0]-PixelSizeX)
    plt.subplot(1,2,2,aspect=vertical_zoom)
    plt.ylabel("Z(microns)")
    plt.xlabel("Y (microns)")
    plt.xticks([20,30,40,50])
    plt.ylim(0,110)
    plt.xlim(ConfigData["PixelBoundaryLowerLeft"][1]+PixelSizeY, ConfigData["PixelBoundaryUpperRight"][1]-PixelSizeY)




    plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_Paths_%d.pdf"%run)

