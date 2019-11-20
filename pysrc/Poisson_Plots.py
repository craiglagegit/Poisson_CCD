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

NumPixelsPlotted = 4
nycenter2 = nycenter
nxmin = int(nxcenter - (NumPixelsPlotted * ScaleFactor * GridsPerPixelX)/2)
nxmax = int(nxcenter + (NumPixelsPlotted * ScaleFactor * GridsPerPixelX)/2)
nymin = int(nycenter - (NumPixelsPlotted * ScaleFactor * GridsPerPixelY)/2)
nymax = int(nycenter + (NumPixelsPlotted * ScaleFactor * GridsPerPixelY)/2)

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
plt.ylim(-20.0, 20.0)
for slicez in [0,1,2,3,10]:
    plt.plot(dat.x[:],dat.phi[:,dat.ny-1,slicez], label = '$z_0=%.1f$'%dat.z[slicez])
plt.xlabel('$x$ [um]')
plt.ylabel('$\phi(x,y_B,z_0)$ [V]')
plt.legend()

plt.subplot(2,2,3)
plt.title("Left Edge")
plt.ylim(-20.0, 20.0)
for slicez in [0,1,2,3,10]:
    plt.plot(dat.y[:],dat.phi[0,:,slicez], label = '$z_0=%.1f$'%dat.z[slicez])
plt.xlabel('$y$ [um]')
plt.ylabel('$\phi(x_L,y,z_0)$ [V]')
plt.legend()

plt.subplot(2,2,4)
plt.title("Right Edge")
plt.ylim(-20.0, 20.0)
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
plt.plot(dat.z[0:phinumzs],(dat.phi[nxcenter2,nycenter2,0:phinumzs]+dat.phi[nxcenter2-1,nycenter2,0:phinumzs]+dat.phi[nxcenter2,nycenter2-1,0:phinumzs]+dat.phi[nxcenter2-1,nycenter2-1,0:phinumzs])/4.0, label = "100K electrons")

nxcenter3 = nxcenter2 + 2 * GridsPerPixelX*ScaleFactor

plt.plot(dat.z[0:phinumzs],(dat.phi[nxcenter3,nycenter2,0:phinumzs]+dat.phi[nxcenter3-1,nycenter2,0:phinumzs]+dat.phi[nxcenter3,nycenter2-1,0:phinumzs]+dat.phi[nxcenter3-1,nycenter2-1,0:phinumzs])/4.0, label = "Empty Well")
plt.legend(loc = "lower left")
#plt.xlabel("Z-Dimension (microns)")
plt.ylabel('$\phi(x,y,z)$ [V]',fontsize=12)
plt.ylim(-5.0, 25.0)
plt.xlim(0.0,4.0)

plt.subplot(2,3,4)
plt.title("Rho-Collect Gate", fontsize=12)
plt.plot(dat.z[0:numzs], dat.rho[nxcenter2,nycenter2,0:numzs], label = "Fixed charge", color='green')
plt.plot(dat.z[0:elecnumzs], -ChargeFactor * dat.Elec[nxcenter2,nycenter2,0:elecnumzs] / dat.dz[0:elecnumzs], label = "100K electrons", color='blue')
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

print("Making 1D potential Plots #2 \n")
plt.figure()
plt.suptitle("1D Potentials in Storage Region. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
plt.subplots_adjust(hspace=0.3, wspace=0.3)
slicez = 16 * ScaleFactor
plt.subplot(1,2,1, aspect = 1)
plt.title("Phi, z = %.2f"%dat.z[slicez])
levels = np.linspace(-20.0, 20.0, 21)
[yy,xx] = np.meshgrid(dat.y[nymin:nymax],dat.x[nxmin:nxmax])
plt.contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels, linewidths=0.1,cmap='jet')
plt.contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,cmap='jet')
plt.xlabel("X-Dimension (microns)")
plt.ylabel("Y-Dimension (microns)")
plt.plot([dat.x[nxmin+1],dat.x[nxmax-1]],[dat.y[nycenter2],dat.y[nycenter2]],ls = "-", color="k")
plt.plot([dat.x[nxcenter2],dat.x[nxcenter2]],[dat.y[nymin+1],dat.y[nymax-1]],ls = "-", color="k")
plt.colorbar(orientation='horizontal').set_label('$\phi(x,y,z)$ [V]')

plt.subplot(1,2,2)
plt.title("Phi-Collect Gate, z = %.2f"%dat.z[slicez])
plt.plot(dat.x[nxmin:nxmax],dat.phi[nxmin:nxmax, nycenter2, slicez], label = "XSlice, y = %.2f"%dat.y[nycenter2])
plt.plot(dat.y[nymin:nymax],dat.phi[nxcenter2,nymin:nymax, slicez], label = "YSlice, x = %.2f"%dat.x[nxcenter2])
plt.ylim(-10.0, 20.0)
plt.xlim(dat.x[nxmin],dat.x[nxmax])
plt.xlabel("X,Y-Dimension (microns)")
plt.ylabel("Potential(Volts)")
plt.legend()
plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_1D_Potentials_2_%d.pdf"%run)

print("Making 1D potential Plots #3 \n")
slicezs = [4,8,12,14,16,18,20,24]
plt.figure()
plt.suptitle("1D Potentials in Storage Region. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
plt.subplots_adjust(hspace=0.3, wspace=0.3)

for i,slicez in enumerate(slicezs):
    plt.subplot(2,4,i+1)
    plt.title("Phi, z = %.2f"%dat.z[slicez])
    plt.plot(dat.x[nxmin:nxmax],dat.phi[nxmin:nxmax, nycenter2, slicez], label = "XSlice, y = %.2f"%dat.y[nycenter2])
    plt.plot(dat.y[nymin:nymax],dat.phi[nxcenter2,nymin:nymax, slicez], label = "YSlice, x = %.2f"%dat.x[nxcenter2])
    plt.ylim(-10.0, 20.0)
    plt.xlim(dat.x[nxmin],dat.x[nxmax])
    plt.xlabel("X,Y-Dimension (microns)")
    plt.ylabel("Potential(Volts)")
    plt.legend()
plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_1D_Potentials_3_%d.pdf"%run)

print("Making 1D potential Plots #5 \n")
slicezs = np.array([8,9,10,12]) * ScaleFactor
plt.figure()
plt.suptitle("1D Potentials in Storage Region. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
plt.subplots_adjust(hspace=0.3, wspace=0.3)

for i,slicez in enumerate(slicezs):
    plt.subplot(2,4,i+1)
    plt.title("Phi, z = %.2f"%dat.z[slicez])
    plt.plot(dat.x[nxmin:nxmax],dat.phi[nxmin:nxmax, nycenter2, slicez], label = "XSlice, y = %.2f"%dat.y[nycenter2])
    plt.plot(dat.y[nymin:nymax],dat.phi[nxcenter2,nymin:nymax, slicez], label = "YSlice, x = %.2f"%dat.x[nxcenter2])
    plt.ylim(-10.0, 20.0)
    plt.xlim(dat.x[nxmin],dat.x[nxmax])
    plt.xlabel("X,Y-Dimension (microns)")
    plt.ylabel("Potential(Volts)")
    plt.legend()
for i,slicez in enumerate(slicezs):
    plt.subplot(2,4,i+5)
    plt.title("Holes, z = %.2f"%dat.z[slicez])
    plt.plot(dat.x[nxmin:nxmax],dat.Hole[nxmin:nxmax, nycenter2, slicez], label = "XSlice, y = %.2f"%dat.y[nycenter2])
    plt.plot(dat.y[nymin:nymax],dat.Hole[nxcenter2,nymin:nymax, slicez], label = "YSlice, x = %.2f"%dat.x[nxcenter2])
    #plt.ylim(-5.0, 0.0)
    plt.xlim(dat.x[nxmin],dat.x[nxmax])
    plt.xlabel("X,Y-Dimension (microns)")
    plt.ylabel("Holes")
    plt.legend()
plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_1D_Potentials_5_%d.pdf"%run)

print("Making 1D potential Plots #4 \n")
plt.figure()

plt.suptitle("1D Potentials in Isolation Regions. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
plt.subplots_adjust(hspace=0.3, wspace=0.3)
slicez = 16 * ScaleFactor
nxcenter3 = int(nxcenter2 + GridsPerPixelX * ScaleFactor / 2)
nycenter3 = int(nycenter2 + GridsPerPixelY * ScaleFactor / 2)
plt.subplot(1,2,1, aspect = 1)
plt.title("Phi, z = %.2f"%dat.z[slicez])
levels = np.linspace(-20.0, 20.0, 21)
[yy,xx] = np.meshgrid(dat.y[nymin:nymax],dat.x[nxmin:nxmax])
plt.contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels, linewidths=0.1,cmap='jet')
plt.contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,cmap='jet')
plt.xlabel("X-Dimension (microns)")
plt.ylabel("Y-Dimension (microns)")
plt.plot([dat.x[nxmin+1],dat.x[nxmax-1]],[dat.y[nycenter3],dat.y[nycenter3]],ls = "-", color="k")
plt.plot([dat.x[nxcenter3],dat.x[nxcenter3]],[dat.y[nymin+1],dat.y[nymax-1]],ls = "-", color="k")
#plt.colorbar()

plt.subplot(1,2,2)
plt.title("Phi-Collect Gate, z = %.2f"%dat.z[slicez])
plt.plot(dat.x[nxmin:nxmax],dat.phi[nxmin:nxmax, nycenter3, slicez], label = "XSlice, y = %.2f"%dat.y[nycenter3])
plt.plot(dat.y[nymin:nymax],dat.phi[nxcenter3,nymin:nymax, slicez], label = "YSlice, x = %.2f"%dat.x[nxcenter3])
plt.ylim(-10.0, 10.0)
plt.xlim(dat.x[nxmin],dat.x[nxmax])
plt.xlabel("X,Y-Dimension (microns)")
plt.ylabel("Potential(Volts)")
plt.legend ()
plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_1D_Potentials_4_%d.pdf"%run)


print("Making summary plots\n")
plt.figure()
plt.suptitle("CCD Charge Collection. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
plotcounter = 1
plt.subplots_adjust(hspace=0.3, wspace=0.3)
[yy,xx] = np.meshgrid(dat.y[nymin:nymax],dat.x[nxmin:nxmax])

slicez = 0
plt.subplot(2,2,1, aspect = 1)
plt.title("Phi, z = 0.0")
levels = np.linspace(-10.0, 10.0, 31)
plt.contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels, linewidths=0.1,cmap='jet')
plt.contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,cmap='jet')
plt.xlabel("X-Dimension (microns)")
plt.ylabel("Y-Dimension (microns)")
plt.colorbar().set_label('$\phi(x,y,z)$ [V]')

plt.subplot(2,2,2, aspect = 1)
nxmin2 = int(nxmin - 9 * GridsPerPixelX * ScaleFactor / 2)
nymin2 = int(nymin - 9 * GridsPerPixelY * ScaleFactor / 2)
rho0 = dat.rho[nxmin2,nymin2,slicez+1]

plt.title("Rho, z = %.2f - %.2f"%(dat.z[nzmin],dat.z[nzmax]))
levels = np.linspace(-40.0,20.0,61)
plotarray = np.array(dat.rho[nxmin:nxmax,nymin:nymax,nzmin:nzmax].sum(axis=2)/(nzmax-nzmin))
plt.contour(xx,yy,plotarray, levels, linewidths=0.1,cmap='jet')
plt.contourf(xx,yy,plotarray, levels, cmap='jet')
plt.xlabel("X-Dimension (microns)")
plt.ylabel("Y-Dimension (microns)")
plt.colorbar().set_label('$\\rho(x,y,z) / \epsilon_{Si}$ [V/um$^2$]')

slicez = 8 * ScaleFactor
plt.subplot(2,2,3, aspect = 1)
plt.title("Phi, z = %.2f"%dat.z[slicez])
levels = np.linspace(-20.0, 20.0, 21)
plt.contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels, linewidths=0.1,cmap='jet')
plt.contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,cmap='jet')
plt.xlabel("X-Dimension (microns)")
plt.ylabel("Y-Dimension (microns)")
plt.plot([dat.x[nxmin+1],dat.x[nxmax-1]],[dat.y[nycenter2],dat.y[nycenter2]],ls = "-", color="k")
plt.plot([dat.x[nxcenter2],dat.x[nxcenter2]],[dat.y[nymin+1],dat.y[nymax-1]],ls = "-", color="k")
plt.colorbar().set_label('$\phi(x,y,z)$ [V]')

slicez = 16 * ScaleFactor
plt.subplot(2,2,4, aspect = 1)
plt.title("Phi, z = %.2f"%dat.z[slicez])
plt.contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels, linewidths=0.1,cmap='jet')
plt.contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,cmap='jet')
plt.xlabel("X-Dimension (microns)")
plt.ylabel("Y-Dimension (microns)")
plt.colorbar().set_label('$\phi(x,y,z)$ [V]')

plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_Summary_1_%d.pdf"%run)

plt.figure()
plt.suptitle("CCD Charge Collection. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
plt.subplots_adjust(hspace=0.3, wspace=0.3)
levels = np.linspace(-20.0, 20.0, 41)

plt.subplot(1,2,1)
plt.title("Phi and (-)E in Gate Region. y = %.2f"%dat.y[nycenter2])
plt.xlabel("X-Dimension (microns)")
plt.ylabel("Z-Dimension (microns)")
nzmin = 0
nzmax = 16 * ScaleFactor
[zz,xx] = np.meshgrid(dat.z[nzmin:nzmax],dat.x[nxmin:nxmax])
plt.contour(xx,zz,dat.phi[nxmin:nxmax,nycenter2,nzmin:nzmax],levels, linewidths=0.1,cmap='jet')
plt.contourf(xx,zz,dat.phi[nxmin:nxmax,nycenter2,nzmin:nzmax],levels,cmap='jet')
plt.colorbar().set_label('$\phi(x,y,z)$ [V]')
if ConfigData["LogEField"] == 1:
    nzmin = 1
    [zz,xx] = np.meshgrid(dat.z[nzmin:nzmax],dat.x[nxmin:nxmax])
    plt.quiver(xx, zz, dat.Ex[nxmin:nxmax,nycenter2,nzmin:nzmax], dat.Ez[nxmin:nxmax,nycenter2,nzmin:nzmax], color='b', scale = 150.0)#, scale_units="width")
[zz,xx] = np.meshgrid(dat.z[nzmin:nzmax],dat.x[nxmin:nxmax])
plt.ylim(zz[0,0], zz[-1,-1])
plt.xlim(xx[0,0], xx[-1,-1])


plt.subplot(1,2,2)
plt.title("Phi and (-)E in Gate Region. x = %.2f"%dat.x[nxcenter2])
plt.xlabel("Y-Dimension (microns)")
plt.ylabel("Z-Dimension (microns)")
nzmin = 0
[zz,yy] = np.meshgrid(dat.z[nzmin:nzmax],dat.y[nymin:nymax])
plt.contour(yy,zz,dat.phi[nxcenter2,nymin:nymax,nzmin:nzmax],levels, linewidths=0.1,cmap='jet')
plt.contourf(yy,zz,dat.phi[nxcenter2,nymin:nymax,nzmin:nzmax],levels,cmap='jet')
plt.colorbar().set_label('$\phi(x,y,z)$ [V]')
if ConfigData["LogEField"] == 1:
    nzmin = 1
    [zz,yy] = np.meshgrid(dat.z[nzmin:nzmax],dat.y[nymin:nymax])
    plt.quiver(yy, zz, dat.Ey[nxcenter2,nymin:nymax,nzmin:nzmax], dat.Ez[nxcenter,nymin:nymax,nzmin:nzmax], color='b', scale = 150.0)#, scale_units="width")
nzmin = 0
[zz,yy] = np.meshgrid(dat.z[nzmin:nzmax],dat.y[nymin:nymax])
plt.ylim(zz[0,0], zz[-1,-1])
plt.xlim(yy[0,0], yy[-1,-1])

plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_Summary_2_%d.pdf"%run)

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
    if phase == 0:
        xin = float(values[3])
        yin = float(values[4])
    elif phase == 4:
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

    vertical_zoom = 1
    plt.figure()
    plt.suptitle("Electron Path Plot - Vertical Zoom = %d"%vertical_zoom, fontsize = 24)
    plt.subplots_adjust(wspace=0.2)

    for line in lines:
        values = line.split()
        phase = int(values[2])
        if phase == 0:
            xin = float(values[3])
            yin = float(values[4])
            if (yin > yline - .10) and (yin < yline + .10):
                YPlotThisID = True
                xpaths=[]
                zxpaths=[]
            else:
                YPlotThisID = False
            if (xin > xline - .10) and (xin < xline + .10):
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
                    pixxin = int(xin/10.0)
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
                    pixyin = int(yin/10.0)
                    if pixyin % 2 == 0:
                        color = "red"
                    else:
                        color = "black"
                    plt.subplot(1,2,2,aspect=vertical_zoom)                        
                    plt.plot(ypaths, zypaths, color = color, linewidth = 0.1)

    plt.subplot(1,2,1,aspect=vertical_zoom)
    plt.ylabel("Z(microns)")
    plt.xlabel("X (microns)")
    plt.ylim(0.0,110.0)
    plt.xlim(ConfigData["PixelBoundaryLowerLeft"][0], ConfigData["PixelBoundaryUpperRight"][0])
    plt.subplot(1,2,2,aspect=vertical_zoom)
    plt.ylabel("Z(microns)")
    plt.xlabel("Y (microns)")
    plt.ylim(0.0,110.0)
    plt.xlim(ConfigData["PixelBoundaryLowerLeft"][1], ConfigData["PixelBoundaryUpperRight"][1])
    plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_Paths_%d.pdf"%run)

