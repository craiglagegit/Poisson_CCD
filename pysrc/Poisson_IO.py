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
nxcenter2 = nxcenter
nycenter2 = nycenter

nxmin = 0
nxmax = nxx
nymin = 0
nymax = nyy

nzmin = 0
nzmax = 16 * ScaleFactor

"""
print "OG"
xplot = 88 * ScaleFactor
yplot = 53 * ScaleFactor
for i in range(dat.Elec.shape[2]):
    print "k = %d, z=%.3f, phi = %.3f, holes = %.3f, elec = %.3f, rho = %.3f, Ez = %.3f"%(i,dat.z[i],dat.phi[xplot,yplot,i],dat.hole[xplot,yplot,i],dat.Elec[xplot,yplot,i],dat.rho[xplot,yplot,i],dat.Ez[xplot,yplot,i])

print "Last Serial"
xplot = 88 * ScaleFactor
yplot = 45 * ScaleFactor
for i in range(dat.Elec.shape[2]):
    print "k = %d, z=%.3f, phi = %.3f, holes = %.3f, elec = %.3f, rho = %.3f, Ez = %.3f"%(i,dat.z[i],dat.phi[xplot,yplot,i],dat.hole[xplot,yplot,i],dat.Elec[xplot,yplot,i],dat.rho[xplot,yplot,i],dat.Ez[xplot,yplot,i])

print "Hole Region"
xplot = 34 * ScaleFactor
for i in range(dat.Elec.shape[2]):
    print "k = %d, z=%.3f, phi = %.3f, holes = %.3f, rho = %.3f, Ez = %.3f"%(i,dat.z[i],dat.phi[xplot,nycenter,i],dat.hole[xplot,nycenter,i],dat.rho[xplot,nycenter,i],dat.Ez[xplot,nycenter,i])

print "Center"
xplot = nxcenter
for i in range(dat.Elec.shape[2]):
    print "k = %d, z=%.3f, phi = %.3f, holes = %.3f, rho = %.3f, Ez = %.3f"%(i,dat.z[i],dat.phi[xplot,nycenter,i],dat.hole[xplot,nycenter,i],dat.rho[xplot,nycenter,i],dat.Ez[xplot,nycenter,i])

print "Right Edge"
#xplot = 560 # Top and Bottom
xplot = 150 # Left and Right
for i in range(dat.Elec.shape[2]):
    print "k = %d, z=%.3f, phi = %.3f, elec = %.3f, rho = %.3f, Ez = %.3f"%(i,dat.z[i],dat.phi[xplot,nycenter,i],dat.Elec[xplot,nycenter,i],dat.rho[xplot,nycenter,i],dat.Ez[xplot,nycenter,i])
"""
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
plt.ylim(-25.0, 25.0)
for slicez in [0,1,2,3,10]:
    plt.plot(dat.x[:],dat.phi[:,0,slicez], label = '$z_0=%.1f$'%dat.z[slicez])
plt.xlabel('$x$ [um]')
plt.ylabel('$\phi(x,y_F,z_0)$ [V]')
plt.legend(loc = 'lower left')

plt.subplot(2,2,2)
plt.title("Back Edge")
plt.ylim(-25.0, 25.0)
for slicez in [0,1,2,3,10]:
    plt.plot(dat.x[:],dat.phi[:,dat.ny-1,slicez], label = '$z_0=%.1f$'%dat.z[slicez])
plt.xlabel('$x$ [um]')
plt.ylabel('$\phi(x,y_B,z_0)$ [V]')
plt.legend(loc = 'lower left')

plt.subplot(2,2,3)
plt.title("Left Edge")
plt.ylim(-25.0, 25.0)
for slicez in [0,1,2,3,10]:
    plt.plot(dat.y[:],dat.phi[0,:,slicez], label = '$z_0=%.1f$'%dat.z[slicez])
plt.xlabel('$y$ [um]')
plt.ylabel('$\phi(x_L,y,z_0)$ [V]')
plt.legend(loc = 'lower left')

plt.subplot(2,2,4)
plt.title("Right Edge")
plt.ylim(-25.0, 25.0)
for slicez in [0,1,2,3,10]:
    plt.plot(dat.y[:],dat.phi[dat.nx-1,:,slicez], label = '$z_0=%.1f$'%dat.z[slicez])
plt.xlabel('$y$ [um]')
plt.ylabel('$\phi(x_R,y,z_0)$ [V]')
plt.legend(loc = 'lower left')

plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_Edge_Potentials_%d.pdf"%run)

print("Making 1D Potential and Charge Density plots\n")
plt.figure()

plt.suptitle("1D Potential and Charge Density Slices. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
plotcounter = 1
plt.subplots_adjust(hspace=0.3, wspace=0.3)
phinumzs = 160
numzs = 160
elecnumzs = ConfigData["Nzelec"] * ConfigData["ScaleFactor"]

plt.subplot(1,2,1)
plt.title("Vertical Isolation", fontsize=12)
nxcenter2 = 32 * ScaleFactor
nycenter2 = 72 * ScaleFactor
plt.plot(dat.z[0:phinumzs],dat.phi[nxcenter2,nycenter2,0:phinumzs], label = "Ground Contact")
nxcenter2 = 72 * ScaleFactor
nycenter2 = 64 * ScaleFactor
plt.plot(dat.z[0:phinumzs],dat.phi[nxcenter2,nycenter2,0:phinumzs], label = "Between Trans.")
plt.xlabel("Z-Dimension (microns)")
plt.ylabel('$\phi(z)$ [V]',fontsize=12)
plt.ylim(-5.0, 5.0)
plt.xlim(0.0,20.0)
plt.legend(loc = 'upper right')

plt.subplot(1,2,2)
nxcenter2 = 88 * ScaleFactor
nymin4 = 0
nymax4 = 64 * ScaleFactor

plt.title("Phi-Serial Chain", fontsize=12)
for sz in [4,6,7,8,9,10]:
    slicez = sz * ScaleFactor
    plt.plot(dat.y[nymin4:nymax4],dat.phi[nxcenter2,nymin4:nymax4,slicez], label = '$z_0=%.1f$'%dat.z[slicez])
plt.xlabel("Y-Dimension (microns)")
plt.ylabel('$\phi(y)$ [V]',fontsize=12)
plt.ylim(-5.0, 30.0)
plt.xlim(0.0,40.0)
plt.legend(loc = 'upper left')

plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_1D_Potentials_%d.pdf"%run)


print("Making summary plots\n")
plt.figure()
plt.suptitle("Potentials on Bottom. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
[yy,xx] = np.meshgrid(dat.y[nymin:nymax],dat.x[nxmin:nxmax])

slicez = 0
plt.subplot(1,2,1, aspect = 1)
levels = np.linspace(-10.0, 26.0, 37)
plt.title ("Z = %.2f"%dat.z[slicez])
plt.contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,linewidths=0.1,cmap='jet')
plt.contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,cmap='jet')
plt.xlabel("X-Dimension (microns)")
plt.ylabel("Y-Dimension (microns)")
cb=plt.colorbar(orientation='horizontal')
cb.set_ticks([-10.0,0.0,10.0,20.0,30.0])
cb.set_label('$\phi(x,y)$ [V]')

slicez = 32
plt.subplot(1,2,2, aspect = 1)
levels = np.linspace(-5.0, 20.0, 51)
plt.title ("Z = %.2f"%dat.z[slicez])
plt.contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,linewidths=0.1,cmap='jet')
plt.contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,cmap='jet')
plt.xlabel("X-Dimension (microns)")
plt.ylabel("Y-Dimension (microns)")
cb=plt.colorbar(orientation='horizontal')
cb.set_ticks([-5.0,0.0,5.0,10.0,15.0,20.0])
cb.set_label('$\phi(x,y)$ [V]')

plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_Summary_1_%d.pdf"%run)

print("Making summary plots\n")
plt.figure()
plt.suptitle("Potentials on Bottom. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
[yy,xx] = np.meshgrid(dat.y[nymin:nymax],dat.x[nxmin:nxmax])

slicez = 0
plt.subplot(1,2,1, aspect = 1)
levels = np.linspace(-10.0, 26.0, 37)
plt.title ("Z = %.2f"%dat.z[slicez])
plt.contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,linewidths=0.1,cmap='jet')
plt.contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,cmap='jet')
plt.xlabel("X-Dimension (microns)")
plt.ylabel("Y-Dimension (microns)")
cb=plt.colorbar(orientation='horizontal')
cb.set_ticks([-10.0,0.0,10.0,20.0,30.0])
cb.set_label('$\phi(x,y)$ [V]')

slicez = 9
plt.subplot(1,2,2, aspect = 1)
levels = np.linspace(-10.0, 26.0, 37)
plt.title ("Z = %.2f"%dat.z[slicez])
plt.contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,linewidths=0.1,cmap='jet')
plt.contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,cmap='jet')
plt.xlabel("X-Dimension (microns)")
plt.ylabel("Y-Dimension (microns)")
cb=plt.colorbar(orientation='horizontal')
cb.set_ticks([-5.0,0.0,5.0,10.0,15.0,20.0])
cb.set_label('$\phi(x,y)$ [V]')

plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_Summary_2_%d.pdf"%run)
