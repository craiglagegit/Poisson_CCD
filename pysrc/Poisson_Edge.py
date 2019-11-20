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
nxcenter2 = int(((ConfigData["PixelRegionLowerLeft_0"][0] + ConfigData["PixelRegionUpperRight_0"][0]) / 2.0 - ConfigData["SimulationRegionLowerLeft"][0]) / (GridsPerPixelX * ScaleFactor))
nycenter2 = int(((ConfigData["PixelRegionLowerLeft_0"][1] + ConfigData["PixelRegionUpperRight_0"][1]) / 2.0 - ConfigData["SimulationRegionLowerLeft"][1]) / (GridsPerPixelY * ScaleFactor))

nxmin = 0
nxmax = nxx
nymin = 0
nymax = nyy

nzmin = 0
nzmax = 16 * ScaleFactor

# Create the output directory if it doesn't exist
if not os.path.isdir(outputfiledir+"/plots"):
    os.mkdir(outputfiledir+"/plots")

mpl.rcParams['contour.negative_linestyle'] = 'solid'
mpl.rcParams.update({'font.size': 6})

print("Making array edge potential plots\n")
plt.figure()
plt.suptitle("Array Edge Potentials. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
plt.subplot(2,2,1)
plt.title("Front Edge")
plt.ylim(-20.0, 20.0)
for slicez in [0,1,2,3,10]:
    plt.plot(dat.x[:],dat.phi[:,0,slicez], label = '$z_0=%.1f$'%dat.z[slicez])
plt.xlabel('$x$ [um]')
plt.ylabel('$\phi(x,y_F,z_0)$ [V]')
plt.legend()

plt.subplot(2,2,2)
plt.title("Back Edge")
plt.ylim(-80.0, 20.0)
for slicez in [0,1,2,3,10]:
    plt.plot(dat.x[:],dat.phi[:,dat.ny-1,slicez], label = '$z_0=%.1f$'%dat.z[slicez])
plt.xlabel('$x$ [um]')
plt.ylabel('$\phi(x,y_B,z_0)$ [V]')
plt.legend()

plt.subplot(2,2,3)
plt.title("Left Edge")
plt.ylim(-75.0, 25.0)
for slicez in [0,1,2,3,10]:
    plt.plot(dat.y[:],dat.phi[0,:,slicez], label = '$z_0=%.1f$'%dat.z[slicez])
plt.xlabel('$y$ [um]')
plt.ylabel('$\phi(x_L,y,z_0)$ [V]')
plt.legend()

plt.subplot(2,2,4)
plt.title("Right Edge")
plt.ylim(-75.0, 25.0)
for slicez in [0,1,2,3,10]:
    plt.plot(dat.y[:],dat.phi[dat.nx-1,:,slicez], label = '$z_0=%.1f$'%dat.z[slicez])
plt.xlabel('$y$ [um]')
plt.ylabel('$\phi(x_R,y,z_0)$ [V]')
plt.legend()

plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_Edge_Potentials_%d.pdf"%run)

print("Making 1D potential and Charge Density plots\n")
plt.figure()

plt.suptitle("1D Potential and Charge Density Slices. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
plotcounter = 1
plt.subplots_adjust(hspace=0.3, wspace=0.3)
phinumzs = 160
numzs = 160
elecnumzs = ConfigData["Nzelec"] * ConfigData["ScaleFactor"]

plt.subplot(2,3,1)
plt.title("Phi-Collect Gate", fontsize=12)

nxcenter3 = int(nxcenter2 + 2 * GridsPerPixelX*ScaleFactor)

plt.plot(dat.z[0:phinumzs],(dat.phi[nxcenter3,nycenter2,0:phinumzs]+dat.phi[nxcenter3-1,nycenter2,0:phinumzs]+dat.phi[nxcenter3,nycenter2-1,0:phinumzs]+dat.phi[nxcenter3-1,nycenter2-1,0:phinumzs])/4.0, label = "Empty Well")
plt.legend(loc = "lower left")
#plt.xlabel("Z-Dimension (microns)")
plt.ylabel('$\phi(x,y,z)$ [V]',fontsize=12)
plt.ylim(-5.0, 25.0)
plt.xlim(0.0,4.0)

plt.subplot(2,3,4)
plt.title("Rho-Collect Gate", fontsize=12)
plt.plot(dat.z[0:numzs], dat.rho[nxcenter2,nycenter2,0:numzs], label = "Fixed charge", color='green')
plt.plot(dat.z[0:elecnumzs], -ChargeFactor * dat.Elec[nxcenter3,nycenter2,0:elecnumzs] / dat.dz[0:elecnumzs], label = "Empty well", color='magenta')

plt.legend(loc = "upper right")
plt.xlabel("Z-Dimension (microns)", fontsize=12)
plt.ylabel('$\\rho(x,y,z)/\epsilon_{Si}$ [V/um$^2$]',fontsize=12)
plt.ylim(-80.0, 250.0)
plt.xlim(0.0,4.0)
nxcenter3 = int(nxcenter2 + GridsPerPixelX * ScaleFactor / 2)
nycenter3 = nycenter2
nycenter4 = int(nycenter2 + GridsPerPixelY * ScaleFactor / 2)
plt.subplot(2,3,2)
plt.title("Phi-ChanStop", fontsize=12)
plt.plot(dat.z[0:phinumzs],dat.phi[nxcenter3,nycenter3,0:phinumzs], label = "CollectGate")
plt.plot(dat.z[0:phinumzs],dat.phi[nxcenter3,nycenter4,0:phinumzs], label = "Barrier Gate")
plt.legend(loc = "lower left")
#plt.xlabel("Z-Dimension (microns)")
#plt.ylabel('$\phi(x,y,z)$ [V]',fontsize=9)
plt.ylim(-10.0, 12.0)
plt.xlim(0.0,10.0)
plt.subplot(2,3,5)
plt.title("Rho-ChanStop", fontsize=12)

plt.plot(dat.z[0:numzs], dat.rho[nxcenter3,nycenter3,0:numzs], color = 'green', label = "Fixed charge")

plt.plot(dat.z[0:elecnumzs], ChargeFactor / dat.dz[0:elecnumzs] * dat.Hole[nxcenter3,nycenter3,0:elecnumzs], label = "Holes Collect Gate", color = 'red')

plt.plot(dat.z[0:elecnumzs], ChargeFactor / dat.dz[0:elecnumzs] * dat.Hole[nxcenter3,nycenter4,0:elecnumzs], label = "Holes Barrier Gate", color = 'orange', )

plt.legend(loc = "lower right")
plt.xlabel("Z-Dimension (microns)", fontsize=12)
#plt.ylabel('$\\rho(x,y,z)/\epsilon_{Si}$ [V/um$^2$]',fontsize=9)
plt.ylim(-100.0, 100.0)
plt.xlim(0.0,10.0)

nxcenter3 = nxcenter2
nycenter3 = int(nycenter2 + GridsPerPixelY * ScaleFactor / 2)
plt.subplot(2,3,3)
plt.title("Phi-Barrier Gate", fontsize=12)
plt.plot(dat.z[0:phinumzs],dat.phi[nxcenter3,nycenter3,0:phinumzs])
#plt.xlabel("Z-Dimension (microns)")
#plt.ylabel('$\phi(x,y,z)$ [V]',fontsize=9)
plt.ylim(-10.0, 12.0)
plt.xlim(0.0,4.0)
plt.subplot(2,3,6)
plt.title("Rho-Barrier Gate", fontsize=12)

plt.plot(dat.z[0:numzs], dat.rho[nxcenter3,nycenter3,0:numzs], color = 'green', label = "Fixed charge")

plt.plot(dat.z[0:elecnumzs], ChargeFactor / dat.dz[0:elecnumzs] * dat.Hole[nxcenter3,nycenter3,0:elecnumzs], color = 'red', label = "Holes")

plt.legend(loc = "lower left")
plt.xlabel("Z-Dimension (microns)", fontsize=12)
#plt.ylabel('$\\rho(x,y,z)/\epsilon_{Si}$ [V/um$^2$]',fontsize=9)
plt.ylim(-80.0, 250.0)
plt.xlim(0.0,4.0)
plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_1D_Potentials_%d.pdf"%run)

print("Making summary plots\n")
plt.figure()
plt.suptitle("Potentials. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
[yy,xx] = np.meshgrid(dat.y[nymin:nymax],dat.x[nxmin:nxmax])
plt.subplots_adjust(hspace=0)

plt.subplot(2,1,1)
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

slicez = 0
plt.subplot(2,1,2, aspect=1)
plt.title("Potentials at z=0", fontsize = 12)
levels = np.linspace(-60.0, 20.0, 161)
plt.contour(yy,xx,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,linewidths=0.1,cmap='jet')
plt.contourf(yy,xx,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,cmap='jet')
plt.ylabel("X-Dimension (microns)")
plt.yticks([40.0, 70.0])
plt.xlabel("Y-Dimension (microns)")
cb=plt.colorbar(orientation='horizontal')
cb.set_ticks([-60.0,-45.0,-30.0,-15.0,0.0,15.0])
cb.set_label('$\phi(x,y)$ [V]')
plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_Bottom_Potential_Vbb_%.0f.pdf"%ConfigData["Vbb"])


# Next, plots of the pixel boundaries
print("Making pixel plots\n")
plt.figure()
mpl.rcParams['contour.negative_linestyle'] = 'solid'
#mpl.rcParams.update({'font.size': 18})

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
    #zout = float(values[5])        
    if phase == 0:
        xin = float(values[3])
        yin = float(values[4])
    elif phase == 2:
    #elif zout < 1.20:
        xout = float(values[3])
        yout = float(values[4])
        if np.isnan(xout) or np.isnan(yout):
            print("xin = %.3f, yin = %.3f is a nan")
            continue
        pixxout = int(xout/10.0)
        pixyout = int(yout/10.0)
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

for linex in np.linspace(120.0,200.0,9):
    plt.plot((linex,linex),(20.0,70.0),linewidth=1.0, color='blue')

plt.xlabel("X(microns)",fontsize = 18)
plt.xticks([40.0, 70.0])
plt.ylabel("Y(microns)",fontsize = 18)
plt.xlim(ConfigData["PixelBoundaryLowerLeft"][0], ConfigData["PixelBoundaryUpperRight"][0])
plt.ylim(ConfigData["PixelBoundaryLowerLeft"][1], ConfigData["PixelBoundaryUpperRight"][1])


plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_Pixels_%d.pdf"%run)
