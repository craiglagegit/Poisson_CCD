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

#****************SUBROUTINES*****************

def PlotCenterPotentials(ConfigData, run):
    outputfilebase = ConfigData["outputfilebase"]    
    outputfiledir = ConfigData["outputfiledir"]
    #dat = Array3dHDF5(outputfiledir, outputfilebase, ConfigData["LogEField"], run) 

    # This holds all of the data
    nxx = dat.nx - 1
    nyy = dat.ny - 1
    nzz = dat.nz - 1

    nxcenter = int(nxx/2)
    nycenter = int(nyy/2)
    nxcenter2 = nxcenter

    mpl.rcParams['contour.negative_linestyle'] = 'solid'
    mpl.rcParams.update({'font.size': 6})

    print("Making center potential plots\n")

    plt.plot(dat.y[:],dat.phi[nxcenter,:,0], label="z=%.2f"%dat.z[0])
    plt.plot(dat.y[:],dat.phi[nxcenter,:,10], label="z=%.2f"%dat.z[10])
    plt.plot(dat.y[:],dat.phi[nxcenter,:,int(dat.nz/8)], label="z=%.2f"%dat.z[int(dat.nz/8)])
    plt.plot(dat.y[:],dat.phi[nxcenter,:,int(dat.nz/4)], label="z=%.2f"%dat.z[int(dat.nz/4)])
    plt.plot(dat.y[:],dat.phi[nxcenter,:,int(dat.nz/2)], label="z=%.2f"%dat.z[int(dat.nz/2)])
    plt.plot(dat.y[:],dat.phi[nxcenter,:,int(3*dat.nz/4)], label="z=%.2f"%dat.z[int(3*dat.nz/4)])
    zmaxx = min(100.0, dat.z[dat.nz-1])
    plt.plot(dat.y[:],dat.phi[nxcenter,:,dat.nz-1], label="z=%.2f"%zmaxx)
    plt.ylabel('Potential(V)')
    plt.xlabel('Y (micron)')
    plt.legend()

    plt.savefig(outputfiledir+"/plots/Center_Potential.pdf")
    return

def PlotElectronPaths(ConfigData, run):
    outputfilebase = ConfigData["outputfilebase"]    
    outputfiledir = ConfigData["outputfiledir"]
    xmin = ConfigData["PixelBoundaryLowerLeft"][0]
    ymin = ConfigData["PixelBoundaryLowerLeft"][1]
    xmax = ConfigData["PixelBoundaryUpperRight"][0]
    ymax = ConfigData["PixelBoundaryUpperRight"][1]
    PixelSizeX = ConfigData["PixelSizeX"]
    PixelSizeY = ConfigData["PixelSizeY"]
    nx = ConfigData["PixelBoundaryNx"]
    ny = ConfigData["PixelBoundaryNy"]
    NumPixels = ny - 1
    # PixelBoundary needs to go one pixel beyond the edge of the
    # array, so we can track electrons which start outside the array

    array = Array2dSet(xmin, xmax, nx, ymin, ymax-PixelSizeY, NumPixels, 1)

    print("Making array electron path plots\n")

    filename = outputfiledir+"/"+outputfilebase+"_" + str(run) + "_Pts.dat"
    file = open(filename,"r")
    lines = file.readlines()
    file.close()
    lines.remove(lines[0])    

    pixel_edge_paths = np.zeros(NumPixels)

    xline = (xmin + xmax) / 2.0
    plt.figure()
    plt.suptitle("Electron Paths", fontsize = 18)

    plt.title("Electron Paths")
    oldyin = 1000000.0
    #lastpixyout = 10000


    for line in lines:
        values = line.split()
        phase = int(values[2])
        if phase == 0:
            xin = float(values[3])
            yin = float(values[4])
            if (xin > xline - .10) and (xin < xline + .10):
                XPlotThisID = True
                ypaths=[]
                zypaths=[]
            else:
                XPlotThisID = False
            continue
        if XPlotThisID:
            xout = float(values[3])
            yout = float(values[4])
            zout = float(values[5])
            if np.isnan(xout) or np.isnan(yout) or np.isnan(zout):
                continue

            ypaths.append(yout)
            zypaths.append(zout)
            if phase == 2:
                pixyout = int(ypaths[-1]/10.0) - int(array.y[0]/10.0)
                if pixyout % 2 == 0:
                    color = "red"
                else:
                    color = "black"
                plt.plot(ypaths, zypaths, color = color, linewidth = 0.1)
                XPlotThisID = False

                if pixyout >= 0 and pixyout < NumPixels:
                    if yin > pixel_edge_paths[pixyout]:
                        pixel_edge_paths[pixyout] = yin
                
    plt.ylabel("Z(microns")
    plt.xlabel("Y (microns)")
    plt.ylim(0, 100.0)
    plt.xlim(ymin, ymax)

    print("Pixel edges as determined from electron paths = ", pixel_edge_paths, "\n")

    plt.savefig(outputfiledir + "/plots/" + outputfilebase + "_FullPaths.pdf")

    return [array, pixel_edge_paths]
    
def PlotPixelEdges(ConfigData, run):
    outputfilebase = ConfigData["outputfilebase"]
    outputfiledir = ConfigData["outputfiledir"]
    nx = ConfigData["PixelBoundaryNx"]
    ny = ConfigData["PixelBoundaryNy"]
    PixelSizeX = ConfigData["PixelSizeX"]
    PixelSizeY = ConfigData["PixelSizeY"]
    Area_0 = 100.0
    NumAngles = 4 * ConfigData["NumVertices"] + 4

    filename = outputfiledir + '/' + outputfilebase +'_%d_Area.dat'%run
    [area,sim] = ReadAreaFile(filename, nx, ny, 0, ny/2, Area_0)
    plt.figure()
    plt.subplot(1,1,1,aspect = 1)
    plt.title("Pixel Area")

    filename = outputfiledir + '/' + outputfilebase +'_%d_Vertices.dat'%run

    (vx, vy) = ReadVertexFile(filename, nx, ny, NumAngles)

    LineXMin = ConfigData["PixelBoundaryLowerLeft"][0]
    LineXMax = ConfigData["PixelBoundaryUpperRight"][0]
    LineYMin = ConfigData["PixelBoundaryLowerLeft"][1]
    LineYMax = ConfigData["PixelBoundaryUpperRight"][1]


    plt.title("Pixel Vertices")

    for i in range(nx+1):
        xline = ConfigData["PixelBoundaryLowerLeft"][0] + i * PixelSizeX
        plt.plot([xline, xline], [LineYMin, LineYMax], color = 'black', ls = '--')
    for j in range(ny+1):
        yline = ConfigData["PixelBoundaryLowerLeft"][1] + j * PixelSizeY
        plt.plot([LineXMin, LineXMax],[yline, yline], color = 'black', ls = '--')

    for i in range(nx):
        for j in range(ny):
            x = []
            y = []
            for k in range(NumAngles):
                x.append(vx[i,j,k])
                y.append(vy[i,j,k])
            x.append(vx[i,j,0])
            y.append(vy[i,j,0])
            plt.plot(x, y, lw = 0.5)

    NumPixels = ny - 1
    # Need to track one pixel beyond ymax because of shift
    pixel_edge_vertices = np.zeros(NumPixels)
    for i in range(NumPixels):
        # Edge is the average of vertices 8 and 9, which are the center of the top edge
        pixel_edge_vertices[i] = (vy[0,i,8] + vy[0,i,9]) / 2.0

    print("Pixel edges as determined from vertex finding = ", pixel_edge_vertices, "\n")

    plt.savefig(outputfiledir + "/plots/PixelVertices.pdf")

    return pixel_edge_vertices


def PlotEdgeShift(ConfigData, array, pixel_edge_paths, pixel_edge_vertices):
    outputfilebase = ConfigData["outputfilebase"]
    outputfiledir = ConfigData["outputfiledir"]
    xmin = ConfigData["PixelBoundaryLowerLeft"][0]
    ymin = ConfigData["PixelBoundaryLowerLeft"][1]
    xmax = ConfigData["PixelBoundaryUpperRight"][0]
    ymax = ConfigData["PixelBoundaryUpperRight"][1]
    nx = ConfigData["PixelBoundaryNx"]
    ny = ConfigData["PixelBoundaryNy"]
    NumPixels = ny - 1
    # Need to track one pixel beyond ymax because of shift

    pixel_edge_best = np.zeros(NumPixels)
    pixel_edge_shift = np.zeros(NumPixels)    
    for i in range(NumPixels):
        if abs(pixel_edge_paths[i] - pixel_edge_vertices[i]) < 0.50:
            # If the two agree within 0.5 micron, use the vertex finding value
            # otherwise, use the electron path value
            pixel_edge_best[i] = pixel_edge_vertices[i]
        else:
            pixel_edge_best[i] = pixel_edge_paths[i]
        pixel_edge_shift[i] = pixel_edge_best[i] - (array.y[i] + array.dy / 2.0)
    plt.figure()

    plt.plot(array.y, pixel_edge_shift, lw=2)
    plt.plot((ymin, ymax),(0.0,0.0),ls="-", lw=4, color="black")
    plt.xlim(ymin, ymax)
    plt.ylim(-25.0, 25.0)
    plt.xlabel("Pixel Center(microns)", fontsize = 18)
    plt.ylabel("Pixel Shift(microns)", fontsize = 18)
    #plt.legend(loc='upper right', handlelength=5)
    plt.savefig(outputfiledir + "/plots/Pixel_Shift.pdf")
    return pixel_edge_shift
    
def PlotSpotRolloff(ConfigData, array, pixel_edge_shift):
    outputfilebase = ConfigData["outputfilebase"]
    outputfiledir = ConfigData["outputfiledir"]
    xmin = ConfigData["PixelBoundaryLowerLeft"][0]
    ymin = ConfigData["PixelBoundaryLowerLeft"][1]
    xmax = ConfigData["PixelBoundaryUpperRight"][0]
    ymax = ConfigData["PixelBoundaryUpperRight"][1]
    nx = ConfigData["PixelBoundaryNx"]
    ny = ConfigData["PixelBoundaryNy"]
    NumPixels = ny - 1
    # Need to track one pixel beyond ymax because of shift
    PixelSizeX = ConfigData["PixelSizeX"]
    PixelSizeY = ConfigData["PixelSizeY"]
    
    sigma = 25.0 / 2.355 # Assumes 25 micron FWHM
    # This matches closely to what we measure, which is a sigma of a little over 10 microns.
    plt.figure()
    plt.title("Spot Rolloff", fontsize=24)
    # The array contains the ideal pixel edges.
    # pixel_edge_shift contains the edge shifts
    xspot = (xmin + xmax) / 2.0

    deltas = []
    step_size = 2.0 # Stepping the spot in 2 micron steps.
    ystart = ymin + 3.0  * PixelSizeY # Start 3 pixels into the array
    numsteps = NumPixels * int(PixelSizeY / step_size)
    yspots = [ystart + step_size * i for i in range(numsteps)] 
    #First, fill the array with data for the spot
    for yspot in yspots:
        for i in range(array.nx):
            for j in range(array.ny):
                xl = array.x[i] - array.dx/2.0 - xspot
                xh = array.x[i] + array.dx/2.0 - xspot
                if j == 0:
                    yl = array.y[j] - array.dy / 2.0 - yspot
                else:
                    yl = array.y[j] - array.dy / 2.0 + pixel_edge_shift[j-1] - yspot
                yh = array.y[j] + array.dy / 2.0 + pixel_edge_shift[j] - yspot
                array.data[i,j] = Area(xl, xh, yl, yh, sigma, sigma, 1.0)
        # Now find the spot centroid
        sumx = 0.0
        sumy = 0.0
        for i in range(array.nx):
            for j in range(array.ny):
                sumx += array.data[i,j] * array.x[i]
                sumy += array.data[i,j] * array.y[j]
        sumx = sumx / array.data.sum()
        sumy = sumy / array.data.sum()
        delta = sumy - yspot
        deltas.append(delta)

        #print "Actual xspot = %.3f, Actual yspot = %.3f, Measured xspot = %.3f,  Measured yspot = %.3f, delta = %.3f"%(xspot,yspot, sumx,sumy, delta)

    plt.plot(yspots, deltas)
    plt.plot([ymax - PixelSizeY, ymax - PixelSizeY],[-20.0,20.0]) # Array Edge
    plt.text(ymax - PixelSizeY * 0.98, 5.0, "Array Edge")
    plt.plot([ymin, ymax+2.0*PixelSizeY],[0.0,0.0]) # Zero line
    plt.ylim(-20.0,10.0)
    plt.xlim(ymin, ymax+2.0*PixelSizeY)
    
    plt.xlabel("Spot Centroid(microns)")
    plt.ylabel("Centroid Shift(microns)")
    #plt.legend(loc='upper right', handlelength=5)
    plt.savefig(outputfiledir + "/plots/SpotRolloff.pdf")


#****************MAIN PROGRAM*****************

# First, read the .cfg file

configfile = sys.argv[1]
run = int(sys.argv[2])
ConfigData = ReadConfigFile(configfile)

outputfilebase = ConfigData["outputfilebase"]
outputfiledir = ConfigData["outputfiledir"]

# This holds all of the data
dat = Array3dHDF5(outputfiledir, outputfilebase, run)

# Create the output directory if it doesn't exist
if not os.path.isdir(outputfiledir+"/plots"):
    os.mkdir(outputfiledir+"/plots")

PlotCenterPotentials(ConfigData, run)
[array, pixel_edge_paths] = PlotElectronPaths(ConfigData, run)
pixel_edge_vertices = PlotPixelEdges(ConfigData, run)

pixel_edge_shift = PlotEdgeShift(ConfigData, array, pixel_edge_paths, pixel_edge_vertices)

PlotSpotRolloff(ConfigData, array, pixel_edge_shift)

