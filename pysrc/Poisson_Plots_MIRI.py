#!/usr/bin/env python

#Author: Craig Lage, UC Davis;
#Date: 25-Oct-19

#This program plots the Poisson equation solutions from the C++ Poisson solver
import os, sys, time, h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
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


SensorThickness = ConfigData["SensorThickness"]
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

"""
print("Center")
for k in range(dat.nz):
    print(f"k={k}, z={dat.z[k]}, Rho={dat.rho[nxcenter,nycenter,k]}, elec = {dat.Elec[nxcenter,nycenter,k]}, hole = {dat.Hole[nxcenter,nycenter,k]}, Phi={dat.phi[nxcenter,nycenter,k]}, Ez={dat.Ez[nxcenter,nycenter,k]}")

print("Between contacts")
for k in range(dat.nz):
    print(f"k={k}, z={dat.z[k]}, Rho={dat.rho[nxcenter,nycenter3,k]}, elec = {dat.Elec[nxcenter,nycenter3,k]}, hole = {dat.Hole[nxcenter,nycenter3,k]}, Phi={dat.phi[nxcenter,nycenter3,k]}")

print("Corner")
for k in range(dat.nz):
    print(f"k={k}, z={dat.z[k]}, Rho={dat.rho[nxcenter3,nycenter3,k]}, elec = {dat.Elec[nxcenter3,nycenter3,k]}, hole = {dat.Hole[nxcenter3,nycenter3,k]}, Phi={dat.phi[nxcenter3,nycenter3,k]}")
"""


# Create the output directory if it doesn't exist
if not os.path.isdir(outputfiledir+"/plots"):
    os.mkdir(outputfiledir+"/plots")

mpl.rcParams['contour.negative_linestyle'] = 'solid'
mpl.rcParams.update({'font.size': 6})


print("Making summary plots\n")

plt.figure()
plt.suptitle("MIRI Simulation. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 12)
[yy,xx] = np.meshgrid(dat.y[nymin:nymax],dat.x[nxmin:nxmax])

slicez = 0
ax1=plt.axes([0.07,0.50,0.40,0.40],aspect=1)
ax1.set_title("Phi, z = %.1f"%dat.z[slicez], fontsize=12)
levels = np.linspace(-3.0, -1.0, 21)
ax1.contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,linewidths=0.1,cmap='jet')
img1 = ax1.contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,cmap='jet')
ax1.set_xlabel("X-Dimension (microns)", fontsize=12)
ax1.set_ylabel("Y-Dimension (microns)", fontsize=12)
colorbar(img1)
ax1.set_label('$\phi(x,y,z)$ [V]')


ax2=plt.axes([0.57,0.50,0.40,0.40])
ax2.set_title("Phi", fontsize=12)
ax2.plot(dat.z[0:nzmax], dat.phi[nxcenter,nycenter,0:nzmax], label='Vcontact=%.1f'%CenterVoltage)
ax2.plot(dat.z[0:nzmax], dat.phi[nxcenter3,nycenter3,0:nzmax], label='Vcontact=%.1f'%ContactVoltage)
ax2.set_xlabel("Z-Dimension (microns)", fontsize=12)
ax2.set_ylabel('$\phi(x,y,z)$ [V]',fontsize=9)
ax2.set_ylim(-5.0, 0.0)
ax2.set_xlim(0.0,35.0)
ax2.legend()

logHoleData = np.log10(dat.Hole+1.0E-12)
cmap0 = plt.cm.get_cmap("jet")
zMin = 25 # minimum in microns
zTicks = range(zMin, int(SensorThickness) + 1)
for k in range(nzz):
    if dat.z[k] >= zMin:
        kmin = k
        break

dny = int(ScaleFactor*GridsPerPixelY/16)
nymin = nycenter - dny
nymax = nycenter + dny

ax3=plt.axes([0.10,0.05,0.80,0.32])
ax3.set_title("Log10 Holes: X-Z Slice", fontsize=12)
#ax3.set_xticks([])
ax3.set_yticks(zTicks)
[plotarray, dxx, dyy, levels, my_cmap] = BuildPlotArray(dat, logHoleData, 1, nxmin, nxmax,kmin, nzz, nymin, nymax, 1.0, False, cmap0, specialColors=False)
ax3.contourf(dxx, dyy, plotarray, levels = levels, cmap = my_cmap, extend='both')

plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_Summary_%d.pdf"%run)

"""
# Commenting out the dopant plot.
# Can always uncomment if needed
print("Making dopant plot\n")

numzs = 25 # controls how deep the plot goes
# Converts doping in code units back into cm^-3
# ChargeFactor = (QE*MICRON_PER_M/(EPSILON_0*EPSILON_SI)) / pow(MICRON_PER_CM, 3);
ChargeFactor = 1.6E-19 * 1.0E6 / (11.7 * 8.85E-12) / 1.0E12
PixelVolume = (dat.x[1] - dat.x[0]) * (dat.y[1] - dat.y[0]) * dat.dz[0:numzs] / 1.0E12 
plt.figure()
plt.suptitle("MIRI Simulation. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 12)

plt.title("Contact Doping", fontsize=12)

plt.plot(dat.z[0:numzs], dat.rho[nxcenter,nycenter,0:numzs] / ChargeFactor, color = 'green', label = "Fixed charge")

plt.plot(dat.z[0:numzs], dat.Elec[nxcenter,nycenter,0:numzs] / PixelVolume, label = "Electrons", color = 'red', ls = ':')

plt.legend(loc = "upper right")
plt.xlabel("Z-Dimension (microns)", fontsize=12)
plt.ylabel('$\\rho(x,y,z)$ cm$^3$')

#plt.ylim(-100.0, 100.0)
plt.xlim(0.0,1.0)

plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_ContactDoping_%d.pdf"%run)
"""
