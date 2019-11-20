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
dat = Array3dHDF5(outputfiledir, outputfilebase,  run)

ScaleFactor = ConfigData["ScaleFactor"]
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
nzmax = nzz



"""
print "Vdd contact"
xplot = 41 * ScaleFactor
yplot = 90 * ScaleFactor
for i in range(dat.Elec.shape[2]):
    print "k = %d, z=%.3f, phi = %.3f, holes = %.3f, elec = %.3f, rho = %.3f, Ez = %.3f"%(i,dat.z[i],dat.phi[xplot,yplot,i],dat.Hole[xplot,yplot,i]*ChargeFactor/dat.dz[i],dat.Elec[xplot,yplot,i]*ChargeFactor/dat.dz[i],dat.rho[xplot,yplot,i],dat.Ez[xplot,yplot,i])

print "Vdd Region"
xplot = 41 * ScaleFactor
yplot = 72 * ScaleFactor
for i in range(dat.Elec.shape[2]):
    print "k = %d, z=%.3f, phi = %.3f, holes = %.3f, elec = %.3f, rho = %.3f, Ez = %.3f"%(i,dat.z[i],dat.phi[xplot,yplot,i],dat.Hole[xplot,yplot,i]*ChargeFactor/dat.dz[i],dat.Elec[xplot,yplot,i]*ChargeFactor/dat.dz[i],dat.rho[xplot,yplot,i],dat.Ez[xplot,yplot,i])

print "OG"
xplot = 41 * ScaleFactor
yplot = 48 * ScaleFactor
for i in range(dat.Elec.shape[2]):
    print "k = %d, z=%.3f, phi = %.3f, holes = %.3f, elec = %.3f, rho = %.3f, Ez = %.3f"%(i,dat.z[i],dat.phi[xplot,yplot,i],dat.Hole[xplot,yplot,i]*ChargeFactor/dat.dz[i],dat.Elec[xplot,yplot,i]*ChargeFactor/dat.dz[i],dat.rho[xplot,yplot,i],dat.Ez[xplot,yplot,i])

print "OD"
xplot = 30 * ScaleFactor
yplot = 48 * ScaleFactor
for i in range(dat.Elec.shape[2]):
    print "k = %d, z=%.3f, phi = %.3f, holes = %.3f, elec = %.3f, rho = %.3f, Ez = %.3f"%(i,dat.z[i],dat.phi[xplot,yplot,i],dat.Hole[xplot,yplot,i]*ChargeFactor/dat.dz[i],dat.Elec[xplot,yplot,i]*ChargeFactor/dat.dz[i],dat.rho[xplot,yplot,i],dat.Ez[xplot,yplot,i])

print "Field Region"
xplot = 47 * ScaleFactor
yplot = 17 * ScaleFactor
for i in range(dat.Elec.shape[2]):
    print "k = %d, z=%.3f, phi = %.3f, holes = %.3f, elec = %.3f, rho = %.3f, Ez = %.3f"%(i,dat.z[i],dat.phi[xplot,yplot,i],dat.Hole[xplot,yplot,i]*ChargeFactor/dat.dz[i],dat.Elec[xplot,yplot,i]*ChargeFactor/dat.dz[i],dat.rho[xplot,yplot,i],dat.Ez[xplot,yplot,i])

print "Hole Region"
xplot = 20
for i in range(dat.Elec.shape[2]):
    print "k = %d, z=%.3f, phi = %.3f, holes = %.3f, rho = %.3f, Ez = %.3f"%(i,dat.z[i],dat.phi[xplot,nycenter,i],dat.Hole[xplot,nycenter,i],dat.rho[xplot,nycenter,i],dat.Ez[xplot,nycenter,i])

print "Center"
xplot = nxcenter
for i in range(dat.Elec.shape[2]):
    print "k = %d, z=%.3f, phi = %.3f, holes = %.3f, rho = %.3f, Ez = %.3f"%(i,dat.z[i],dat.phi[xplot,nycenter,i],dat.Hole[xplot,nycenter,i],dat.rho[xplot,nycenter,i],dat.Ez[xplot,nycenter,i])

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


print("Making 1D Potential and Charge Density plots\n")
plt.figure()

plt.suptitle("1D Potential and Charge Density Slices. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
plotcounter = 1
plt.subplots_adjust(hspace=0.3, wspace=0.3)
phinumzs = 160 * ScaleFactor
numzs = 160 * ScaleFactor
elecnumzs = ConfigData["Nzelec"] * ConfigData["ScaleFactor"]

plt.subplot(2,2,1)
plt.title("Phi - S/D Regions", fontsize=12)
nxcenter2 = 30 * ScaleFactor
nycenter2 = 48 * ScaleFactor
plt.plot(dat.z[0:phinumzs],dat.phi[nxcenter2,nycenter2,0:phinumzs], label = "OD")
nxcenter2 = 54 * ScaleFactor
nycenter2 = 48 * ScaleFactor
plt.plot(dat.z[0:phinumzs],dat.phi[nxcenter2,nycenter2,0:phinumzs], label = "OS")
plt.xlabel("Z-Dimension (microns)")
plt.ylabel('$\phi(z)$ [V]',fontsize=12)
plt.ylim(-1.0, 1.0)
plt.xlim(0.0,10.0)
plt.legend(loc = 'upper right')

plt.subplot(2,2,2)
plt.title("Phi - Gate Region", fontsize=12)
nxcenter2 = 41 * ScaleFactor
nycenter2 = 48 * ScaleFactor
plt.plot(dat.z[0:phinumzs],dat.phi[nxcenter2,nycenter2,0:phinumzs], label = "Gate")
plt.xlabel("Z-Dimension (microns)")
plt.ylabel('$\phi(z)$ [V]',fontsize=12)
plt.ylim(-2.0, 2.0)
plt.xlim(0.0,10.0)
plt.legend(loc = 'upper right')

plt.subplot(2,2,3)
plt.title("Rho-S/D Regions", fontsize=12)
nxcenter2 = 30 * ScaleFactor
nycenter2 = 48 * ScaleFactor
plt.plot(dat.z[0:numzs], dat.rho[nxcenter2,nycenter2,0:numzs], label = "Fixed charge", color='green')
plt.plot(dat.z[0:elecnumzs], -ChargeFactor * dat.Elec[nxcenter2,nycenter2,0:elecnumzs] / dat.dz[0:elecnumzs], label = "Electrons", color='blue')
plt.plot(dat.z[0:elecnumzs], ChargeFactor * dat.Hole[nxcenter2,nycenter2,0:elecnumzs] / dat.dz[0:elecnumzs], label = "Holes", color='red')
plt.legend(loc = "lower right")
plt.xlabel("Z-Dimension (microns)", fontsize=12)
plt.ylabel('$\\rho(x,y,z)/\epsilon_{Si}$ [V/um$^2$]',fontsize=12)
plt.ylim(-250.0, 250.0)
plt.xlim(0.0,4.0)

plt.subplot(2,2,4)
plt.title("Rho-Gate Region", fontsize=12)
nxcenter2 = 41 * ScaleFactor
nycenter2 = 48 * ScaleFactor
plt.plot(dat.z[0:numzs], dat.rho[nxcenter2,nycenter2,0:numzs], label = "Fixed charge", color='green')
plt.plot(dat.z[0:elecnumzs], -ChargeFactor * dat.Elec[nxcenter2,nycenter2,0:elecnumzs] / dat.dz[0:elecnumzs], label = "Electrons", color='blue')
plt.plot(dat.z[0:elecnumzs], ChargeFactor * dat.Hole[nxcenter2,nycenter2,0:elecnumzs] / dat.dz[0:elecnumzs], label = "Holes", color='red')
plt.legend(loc = "lower right")
plt.xlabel("Z-Dimension (microns)", fontsize=12)
plt.ylabel('$\\rho(x,y,z)/\epsilon_{Si}$ [V/um$^2$]',fontsize=12)
plt.ylim(-250.0, 250.0)
plt.xlim(0.0,4.0)
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
levels = np.linspace(-5.0, 15.0, 41)
plt.title ("Z = %.2f"%dat.z[slicez])
plt.contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,linewidths=0.1,cmap='jet')
plt.contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,cmap='jet')
plt.xlabel("X-Dimension (microns)")
plt.ylabel("Y-Dimension (microns)")
cb=plt.colorbar(orientation='horizontal')
cb.set_ticks([-5.0,0.0,5.0,10.0,15.0])
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
levels = np.linspace(-10.0, 20.0, 61)
plt.title ("Z = %.2f"%dat.z[slicez])
plt.contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,linewidths=0.1,cmap='jet')
plt.contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,cmap='jet')
plt.xlabel("X-Dimension (microns)")
plt.ylabel("Y-Dimension (microns)")
cb=plt.colorbar(orientation='horizontal')
cb.set_ticks([-10.0,0.0,10.0,20.0])
cb.set_label('$\phi(x,y)$ [V]')

plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_Summary_2_%d.pdf"%run)

plt.figure()
plt.suptitle("Vertical Cross-Section. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
plt.subplots_adjust(hspace=0.3, wspace=0.3)
nxcenter = 38 * ScaleFactor
nymin = 5 * ScaleFactor
nymax = 120 * ScaleFactor
nycenter = 42 * ScaleFactor
nxmin = 16 * ScaleFactor
nxmax = 64 * ScaleFactor
plt.subplot(1,2,1)
plt.title("Phi(x,z) y = %.2f"%dat.y[nycenter])
plt.xlabel("X-Dimension (microns)")
plt.ylabel("Z-Dimension (microns)")
nzmin = 0
nzmax = 160 * ScaleFactor
[zz,xx] = np.meshgrid(dat.z[nzmin:nzmax],dat.x[nxmin:nxmax])
plt.contour(xx,zz,dat.phi[nxmin:nxmax,nycenter,nzmin:nzmax],levels,linewidths=0.1,cmap='jet')
plt.contourf(xx,zz,dat.phi[nxmin:nxmax,nycenter,nzmin:nzmax],levels,cmap='jet')
cb=plt.colorbar(orientation='horizontal')
cb.set_ticks([-60.0,-45.0,-30.0,-15.0,0.0,15.0])
cb.set_label('$\phi(x,z)$ [V]')
[zz,xx] = np.meshgrid(dat.z[nzmin:nzmax],dat.x[nxmin:nxmax])
plt.ylim(zz[0,0], zz[-1,-1])
plt.xlim(xx[0,0], xx[-1,-1])

plt.subplot(1,2,2)
plt.title("Phi(y,z) x = %.2f"%dat.x[nxcenter])
plt.xlabel("Y-Dimension (microns)")
plt.ylabel("Z-Dimension (microns)")
nzmin = 0
nzmax = 160 * ScaleFactor

[zz,yy] = np.meshgrid(dat.z[nzmin:nzmax],dat.y[nymin:nymax])
plt.contour(yy,zz,dat.phi[nxcenter,nymin:nymax,nzmin:nzmax],levels,linewidths=0.1,cmap='jet')
plt.contourf(yy,zz,dat.phi[nxcenter,nymin:nymax,nzmin:nzmax],levels,cmap='jet')
cb=plt.colorbar(orientation='horizontal')
cb.set_ticks([-60.0,-45.0,-30.0,-15.0,0.0,15.0])
cb.set_label('$\phi(y,z)$ [V]')
[zz,yy] = np.meshgrid(dat.z[nzmin:nzmax],dat.y[nymin:nymax])
plt.ylim(zz[0,0], zz[-1,-1])
plt.xlim(yy[0,0], yy[-1,-1])

plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_Summary_3_%d.pdf"%run)


nzmax = 40 * ScaleFactor

plt.figure()
plt.suptitle("Vertical Cross-Section. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
plt.subplots_adjust(hspace=0.3, wspace=0.3)
nxcenter = 38 * ScaleFactor
nymin = 5 * ScaleFactor
nymax = 120 * ScaleFactor
nycenter = 42 * ScaleFactor
nxmin = 16 * ScaleFactor
nxmax = 64 * ScaleFactor
plt.subplot(1,2,1)
plt.title("Elec(x,z) y = %.2f"%dat.y[nycenter])
plt.xlabel("X-Dimension (microns)")
plt.ylabel("Z-Dimension (microns)")
nzmin = 0
nzmax = 160 * ScaleFactor
[zz,xx] = np.meshgrid(dat.z[nzmin:nzmax],dat.x[nxmin:nxmax])
plt.contour(xx,zz,dat.Elec[nxmin:nxmax,nycenter,nzmin:nzmax],linewidths=0.1,cmap='jet')
plt.contourf(xx,zz,dat.Elec[nxmin:nxmax,nycenter,nzmin:nzmax],cmap='jet')
cb=plt.colorbar(orientation='horizontal')
cb.set_ticks([-60.0,-45.0,-30.0,-15.0,0.0,15.0])
cb.set_label('elec(x,z) [V]')
[zz,xx] = np.meshgrid(dat.z[nzmin:nzmax],dat.x[nxmin:nxmax])
plt.ylim(zz[0,0], zz[-1,-1])
plt.xlim(xx[0,0], xx[-1,-1])

plt.subplot(1,2,2)
plt.title("Elec(y,z) x = %.2f"%dat.x[nxcenter])
plt.xlabel("Y-Dimension (microns)")
plt.ylabel("Z-Dimension (microns)")
nzmin = 0
nzmax = 160 * ScaleFactor

[zz,yy] = np.meshgrid(dat.z[nzmin:nzmax],dat.y[nymin:nymax])
plt.contour(yy,zz,dat.Elec[nxcenter,nymin:nymax,nzmin:nzmax],linewidths=0.1,cmap='jet')
plt.contourf(yy,zz,dat.Elec[nxcenter,nymin:nymax,nzmin:nzmax],cmap='jet')
cb=plt.colorbar(orientation='horizontal')
cb.set_ticks([-60.0,-45.0,-30.0,-15.0,0.0,15.0])
cb.set_label('elec(y,z) [V]')
[zz,yy] = np.meshgrid(dat.z[nzmin:nzmax],dat.y[nymin:nymax])
plt.ylim(zz[0,0], zz[-1,-1])
plt.xlim(yy[0,0], yy[-1,-1])

plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_Summary_4_%d.pdf"%run)

plt.figure()
plt.suptitle("Vertical Cross-Section. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
plt.subplots_adjust(hspace=0.3, wspace=0.3)
nxcenter = 38 * ScaleFactor
nymin = 5 * ScaleFactor
nymax = 120 * ScaleFactor
nycenter = 42 * ScaleFactor
nxmin = 16 * ScaleFactor
nxmax = 64 * ScaleFactor
plt.subplot(1,2,1)
plt.title("Hole(x,z) y = %.2f"%dat.y[nycenter])
plt.xlabel("X-Dimension (microns)")
plt.ylabel("Z-Dimension (microns)")
nzmin = 0
nzmax = 160 * ScaleFactor
[zz,xx] = np.meshgrid(dat.z[nzmin:nzmax],dat.x[nxmin:nxmax])
plt.contour(xx,zz,dat.Hole[nxmin:nxmax,nycenter,nzmin:nzmax],linewidths=0.1,cmap='jet')
plt.contourf(xx,zz,dat.Hole[nxmin:nxmax,nycenter,nzmin:nzmax],cmap='jet')
cb=plt.colorbar(orientation='horizontal')
cb.set_ticks([-60.0,-45.0,-30.0,-15.0,0.0,15.0])
cb.set_label('hole(x,z) [V]')
[zz,xx] = np.meshgrid(dat.z[nzmin:nzmax],dat.x[nxmin:nxmax])
plt.ylim(zz[0,0], zz[-1,-1])
plt.xlim(xx[0,0], xx[-1,-1])

plt.subplot(1,2,2)
plt.title("Hole(y,z) x = %.2f"%dat.x[nxcenter])
plt.xlabel("Y-Dimension (microns)")
plt.ylabel("Z-Dimension (microns)")
nzmin = 0
nzmax = 160 * ScaleFactor

[zz,yy] = np.meshgrid(dat.z[nzmin:nzmax],dat.y[nymin:nymax])
plt.contour(yy,zz,dat.Hole[nxcenter,nymin:nymax,nzmin:nzmax],linewidths=0.1,cmap='jet')
plt.contourf(yy,zz,dat.Hole[nxcenter,nymin:nymax,nzmin:nzmax],cmap='jet')
cb=plt.colorbar(orientation='horizontal')
cb.set_ticks([-60.0,-45.0,-30.0,-15.0,0.0,15.0])
cb.set_label('hole(y,z) [V]')
[zz,yy] = np.meshgrid(dat.z[nzmin:nzmax],dat.y[nymin:nymax])
plt.ylim(zz[0,0], zz[-1,-1])
plt.xlim(yy[0,0], yy[-1,-1])

plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_Summary_5_%d.pdf"%run)
