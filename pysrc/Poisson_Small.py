#!/usr/bin/env python

#Author: Craig Lage, UC Davis;
#Date: 25-Oct-19

#This program plots the Poisson equation solutions from the C++ Poisson solver
import os, sys, time, h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
sys.path.append(os.path.realpath('./pysrc'))
from pysubs import *  # These are the plotting subroutines

#****************MAIN PROGRAM*****************

# First, read the .cfg file

configfile = sys.argv[1]
run = int(sys.argv[2])
multi = int(sys.argv[3])

ConfigData = ReadConfigFile(configfile)

outputfilebase = ConfigData["outputfilebase"]
outputfiledir = ConfigData["outputfiledir"]

# This holds all of the data
dat = Array3dHDF5(outputfiledir, outputfilebase, run, multi=multi)

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
nxcenter3 = int(nxcenter2 + GridsPerPixelX * ScaleFactor / 2 / 2**(multi))

NumPixelsPlotted = 1
nycenter2 = nycenter
nycenter3 = int(nycenter2 + GridsPerPixelY * ScaleFactor / 2/ 2**(multi))

nymin = int(nycenter - (NumPixelsPlotted * ScaleFactor * GridsPerPixelY)/2)
nymax = int(nycenter + (NumPixelsPlotted * ScaleFactor * GridsPerPixelY)/2)

nxmin = int(nxcenter - (NumPixelsPlotted * ScaleFactor * GridsPerPixelX)/2)
nxmax = int(nxcenter + (NumPixelsPlotted * ScaleFactor * GridsPerPixelX)/2)

nzmin = 0
nzmax = 16 * ScaleFactor

print(dat.Elec.sum(), dat.Hole.sum())

print("Center")
for k in range(dat.nz):
    print(f"k={k}, z={dat.z[k]}, Rho={dat.rho[nxcenter,nycenter,k]}, elec = {dat.Elec[nxcenter,nycenter,k]}, hole = {dat.Hole[nxcenter,nycenter,k]}, Phi={dat.phi[nxcenter,nycenter,k]}")

print("Between contacts")
for k in range(dat.nz):
    print(f"k={k}, z={dat.z[k]}, Rho={dat.rho[nxcenter3,nycenter3,k]}, elec = {dat.Elec[nxcenter,nycenter,k]}, hole = {dat.Hole[nxcenter,nycenter,k]}, Phi={dat.phi[nxcenter3,nycenter3,k]}")

#sys.exit()


# Create the output directory if it doesn't exist
if not os.path.isdir(outputfiledir+"/plots"):
    os.mkdir(outputfiledir+"/plots")

mpl.rcParams['contour.negative_linestyle'] = 'solid'
mpl.rcParams.update({'font.size': 6})


print("Making 1D Potential and Charge Density plots\n")
plt.figure()

plt.suptitle("1D Potential and Charge Density Slices. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 12)
plotcounter = 1
plt.subplots_adjust(hspace=0.5, wspace=0.3)
phinumzs = dat.nz
numzs = dat.nz
elecnumzs = ConfigData["Nzelec"] * ConfigData["ScaleFactor"]

plt.subplot(2,2,1)
plt.title("Phi-Contact", fontsize=12)
nxcenter3 = nxcenter2

plt.plot(dat.z[0:phinumzs],(dat.phi[nxcenter3,nycenter2,0:phinumzs]+dat.phi[nxcenter3-1,nycenter2,0:phinumzs]+dat.phi[nxcenter3,nycenter2-1,0:phinumzs]+dat.phi[nxcenter3-1,nycenter2-1,0:phinumzs])/4.0)
#plt.legend(loc = "lower left")
plt.xlabel("Z-Dimension (microns)", fontsize=12)
plt.ylabel('$\phi(x,y,z)$ [V]',fontsize=9)
plt.ylim(-5.0, 0.0)
plt.xlim(0.0,35.0)

plt.subplot(2,2,2)
plt.title("Rho-Contact", fontsize=12)
plt.plot(dat.z[0:numzs], (dat.rho[nxcenter2,nycenter2,0:numzs]+dat.rho[nxcenter2-1,nycenter2,0:numzs]+dat.rho[nxcenter2,nycenter2-1,0:numzs]+dat.rho[nxcenter2-1,nycenter2-1,0:numzs])/4.0)

#plt.legend(loc = "upper right")
plt.xlabel("Z-Dimension (microns)", fontsize=12)
plt.ylabel('$log10  \\rho(x,y,z)/\epsilon_{Si}$ [V/um$^2$]',fontsize=9)
plt.yscale('symlog')
plt.ylim(-10.0, 1.0E4)
plt.xlim(0.0,35.0)
nxcenter3 = int(nxcenter2 + GridsPerPixelX * ScaleFactor / 2)
nycenter3 = nycenter2
nycenter4 = int(nycenter2 + GridsPerPixelY * ScaleFactor / 2)
plt.subplot(2,2,3)
plt.title("Phi-Between Contacts", fontsize=12)
plt.plot(dat.z[0:phinumzs],(dat.phi[nxcenter3,nycenter3,0:phinumzs]+dat.phi[nxcenter3-1,nycenter3,0:phinumzs]+dat.phi[nxcenter3,nycenter3-1,0:phinumzs]+dat.phi[nxcenter3-1,nycenter3-1,0:phinumzs])/4.0)
#plt.legend(loc = "lower left")
plt.xlabel("Z-Dimension (microns)", fontsize=12)
plt.ylabel('$\phi(x,y,z)$ [V]',fontsize=9)
plt.ylim(-5.0, 0.0)
plt.xlim(0.0,35.0)
plt.subplot(2,2,4)
plt.title("Rho-Between Contacts", fontsize=12)

plt.plot(dat.z[0:numzs], ((dat.rho[nxcenter3,nycenter3,0:numzs]+dat.rho[nxcenter3-1,nycenter3,0:numzs]+dat.rho[nxcenter3,nycenter3-1,0:numzs]+dat.rho[nxcenter3-1,nycenter3-1,0:numzs])/4.0))

plt.legend(loc = "lower right")
plt.xlabel("Z-Dimension (microns)", fontsize=12)
plt.ylabel('$log10  \\rho(x,y,z)/\epsilon_{Si}$ [V/um$^2$]',fontsize=9)
plt.ylim(-10.0, 1E4)
plt.yscale('symlog')
plt.xlim(0.0,35.0)
plt.savefig(outputfiledir+"/plots/"+outputfilebase+"_1D_Potentials_%d.pdf"%run)

print("Making 1D potential Plots #2 \n")
plt.figure()
plt.suptitle("1D Potentials in Storage Region. Grid = %d*%d*%d."%(nxx,nyy,nzz),fontsize = 18)
plt.subplots_adjust(hspace=0.3, wspace=0.3)
slicez = 8 * ScaleFactor
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
slicezs = [4,6,8,10,12,14,16,18]
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
slicezs = [4,8,12,16]
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
plt.legend()
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
plt.contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,linewidths=0.1,cmap='jet')
plt.contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,cmap='jet')
plt.xlabel("X-Dimension (microns)")
plt.ylabel("Y-Dimension (microns)")
plt.colorbar().set_label('$\phi(x,y,z)$ [V]')

plt.subplot(2,2,2, aspect = 1)
nxmin2 = int(nxmin -  GridsPerPixelX * ScaleFactor / 2)
nymin2 = int(nymin -  GridsPerPixelY * ScaleFactor / 2)
rho0 = dat.rho[nxmin2,nymin2,slicez+1]

plt.title("Rho, z = %.2f - %.2f"%(dat.z[nzmin],dat.z[nzmax]))
levels = np.linspace(-30.0,20.0,51)
plotarray = np.array(dat.rho[nxmin:nxmax,nymin:nymax,nzmin:nzmax].sum(axis=2)/(nzmax-nzmin))
plt.contour(xx,yy,plotarray, levels,linewidths=0.1,cmap='jet')
plt.contourf(xx,yy,plotarray, levels,cmap='jet')
plt.xlabel("X-Dimension (microns)")
plt.ylabel("Y-Dimension (microns)")
plt.colorbar().set_label('$\\rho(x,y,z) / \epsilon_{Si}$ [V/um$^2$]')

slicez = 8 * ScaleFactor
plt.subplot(2,2,3, aspect = 1)
plt.title("Phi, z = %.2f"%dat.z[slicez])
levels = np.linspace(-20.0, 20.0, 21)
plt.contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,linewidths=0.1,cmap='jet')
plt.contourf(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,cmap='jet')
plt.xlabel("X-Dimension (microns)")
plt.ylabel("Y-Dimension (microns)")
plt.plot([dat.x[nxmin+1],dat.x[nxmax-1]],[dat.y[nycenter2],dat.y[nycenter2]],ls = "-", color="k")
plt.plot([dat.x[nxcenter2],dat.x[nxcenter2]],[dat.y[nymin+1],dat.y[nymax-1]],ls = "-", color="k")
plt.colorbar().set_label('$\phi(x,y,z)$ [V]')

slicez = 16 * ScaleFactor
plt.subplot(2,2,4, aspect = 1)
plt.title("Phi, z = %.2f"%dat.z[slicez])
plt.contour(xx,yy,dat.phi[nxmin:nxmax,nymin:nymax,slicez],levels,linewidths=0.1,cmap='jet')
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

