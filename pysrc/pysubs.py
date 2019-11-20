#!/usr/bin/env python

#Author: Craig Lage, UC Davis;
#Date: 24-Oct-19

#These are a set of Python subroutines used in the various plotting routines for
#plotting the output from the Poisson_CCD simulator
import numpy as np
import matplotlib.pyplot as plt
from random import random as rand
import os, sys, time, h5py, xlrd
from scipy.special import erf

#****************SUBROUTINES*****************
class Array3dHDF5(object):
    # This structure holds the data from the C++ run
    # The multi option is used if the multigrids have been written
    def __init__(self, directory, filebase, run, multi=None):
        # First get the grids
        coords = ['x', 'y', 'z']
        for coord in coords:
            if multi is None:
                exec("grid = np.loadtxt('%s/grid_%s.dat', skiprows=1)"%(directory,coord))
            else:
                exec("grid = np.loadtxt('%s/multigrid_%s_%s.dat', skiprows=1)"%(directory,str(multi),coord))
            exec('self.n%s=grid.shape[0]'%coord)
            exec('self.%smin=grid[0,1]'%coord)
            exec('self.%smax=grid[self.n%s-1,3]'%(coord,coord))            
            exec('self.%s=grid[:,2]'%coord)
            exec('self.d%s=(grid[:,3] - grid[:,1])'%coord)                        
        # Then get the data
        datanames = ['phi', 'rho', 'Ex', 'Ey', 'Ez', 'Hole', 'Elec']
        for name in datanames:
            try:
                if multi is None:                
                    exec("hdf = h5py.File('%s/%s_%s_%s.hdf5', 'r')"%(directory,filebase,str(run),name))
                else:
                    exec("hdf = h5py.File('%s/%s_Multi_%s_%s_%s.hdf5', 'r')"%(directory,filebase,str(multi),str(run),name))
                exec('self.%s=np.array(hdf[next(hdf.__iter__())])'%name)
            except:
                continue
        return


class Array2dSet:
    def __init__(self,xmin,xmax,nx,ymin,ymax,ny,nstamps):
        # This packages up a set of nstamps postage stamp images,
        # each image of which is nx * ny pixels
        self.nx=nx
        self.ny=ny
        self.nstamps=nstamps

        self.xmin=xmin
        self.ymin=ymin
        
        self.xmax=xmax
        self.ymax=ymax
        
        self.dx=(xmax-xmin)/nx
        self.dy=(ymax-ymin)/ny
        
        self.x=np.linspace(xmin+self.dx/2,xmax-self.dx/2,nx)
        self.y=np.linspace(ymin+self.dy/2,ymax-self.dy/2,ny)

        self.data=np.zeros([nx,ny,nstamps])
        self.xoffset=np.zeros([nstamps])
        self.yoffset=np.zeros([nstamps])
        self.imax=np.zeros([nstamps])

    
def ReadConfigFile(filename):
    # This reads the Poisson simulator config file for
    # the settings that were run
    # and returns a dictionary with the values

    with open(filename,'r') as file:
        lines=file.readlines()
    lines = [ l.strip() for l in lines ]
    lines = [ l.split() for l in lines if len(l) > 0 and l[0] != '#' ]
    for line in lines:
        if line[1] != '=':
            print("Check line: ",line)
            raise IOError("Error reading config file %s"%filename)
    config = {}
    for line in lines:
        try:
            # First get the ordered pairs
            config.update({line[0]:[eval(line[2]), eval(line[3])]})
        except:
            try:
                # Then get the numeric values
                config.update({line[0]:eval(line[2])})
            except:
                try:
                    # Last, get the strings
                    config.update({line[0]:str(line[2])})
                except:
                    pass
    return config


def BuildPlotArray(dat, plotdata, axis, nxmin, nxmax, nymin, nymax, nzmin, nzmax, ZMult, ForceZero, Vph, Vpl, cmap):
    # This builds a 2D array for use with the charge plotting
    plotarray = np.zeros([nxmax-nxmin+1,nymax-nymin+1])
    dxx = np.zeros([nxmax-nxmin+1])
    dyy = np.zeros([nymax-nymin+1])
    for i in range(nxmax-nxmin+1):
        if axis == 0:
            dxx[i] = dat.z[nxmin+i] * ZMult
        elif axis == 1:
            dxx[i] = dat.x[nxmin+i]            
        elif axis == 2:             
            dxx[i] = dat.x[nxmin+i]            
    for j in range(nymax-nymin+1):
        if axis == 0:
            dyy[j] = dat.y[nymin+j]
        elif axis == 1:
            dyy[j] = dat.z[nymin+j] * ZMult
        elif axis == 2:             
            dyy[j] = dat.y[nymin+j]            

    for i in range(nxmax-nxmin+1):
        for j in range(nymax-nymin+1):
            for k in range(nzmax-nzmin+1):
                if axis == 0: 
                    # Y-Z slice
                    plotarray[i,j] += plotdata[nzmin+k,nymin+j,nxmin+i]
                elif axis == 1: 
                    # X-Z slice
                    plotarray[i,j] += plotdata[nxmin+i,nzmin+k,nymin+j]                
                elif axis == 2: 
                    # X-Y slice
                    plotarray[i,j] += plotdata[nxmin+i,nymin+j,nzmin+k]                

    num_levels = 100
    if ForceZero:
        pdata = plotarray.clip(0.0,1.0E12)
        ndata = plotarray.clip(-1.0E12,0.0)
        plotarray = (pdata/pdata.max()-ndata/ndata.min())
    vmax = plotarray.max()*1.01 + 0.01
    vmin = plotarray.min()*1.01 - 0.01
    levels = np.linspace(vmin, vmax, num_levels)    

    for i in range(nxmax-nxmin):
        for j in range(nymax-nymin):
            if axis == 0:
                if abs(dat.rho[int((nzmin+nzmax)/2),nymin+j,nxmin+i]) < 1.0E-12:
                    plotarray[i,j] = vmax+1.0E6
                if abs(Vph - dat.phi[int((nzmin+nzmax)/2),nymin+j,nxmin+i]) < 1.0E-12 or abs(Vpl - dat.phi[int((nzmin+nzmax)/2),nymin+j,nxmin+i]) < 1.0E-12:
                    plotarray[i,j] = vmin-1.0E6
            elif axis == 1:
                if abs(dat.rho[nxmin+i, int((nzmin+nzmax)/2),nymin+j]) < 1.0E-12:
                    plotarray[i,j] = vmax+100000.0
                if abs(Vph - dat.phi[nxmin+i,int((nzmin+nzmax)/2),nymin+j]) < 1.0E-12 or abs(Vpl - dat.phi[nxmin+i,int((nzmin+nzmax)/2),nymin+j]) < 1.0E-12:
                    plotarray[i,j] = vmin-1.0E6
    if axis == 0:
        for j in range(nymax-nymin+1):
            plotarray[1,j] =  vmin-1.0E6
    if axis == 1:
        for i in range(nxmax-nxmin+1):
            plotarray[i,1] =  vmin-1.0E6

    cmap.set_under("green")
    cmap.set_over("yellow")
                
    [pyy,pxx] = np.meshgrid(dyy, dxx)# Data grid for plots
    return [plotarray, pxx, pyy, levels, cmap]


def BuildPlotSlice(dat, plotdata, axis, nxmin, nxmax, nymin, nymax, nzmin, nzmax):
    # This builds a 1D plot slice for use with the Charge plotting routines
    plotslice = np.zeros([nxmax-nxmin+1])
    xpoints = np.zeros([nxmax-nxmin+1])
    for i in range(nxmax-nxmin+1):
        if axis == 0: 
            xpoints[i] = dat.x[nxmin+i]            
        elif axis == 1:
            xpoints[i] = dat.y[nxmin+i]            
        elif axis == 2:             
            xpoints[i] = dat.z[nxmin+i]

    for i in range(nxmax-nxmin+1):
        for j in range(nymax-nymin+1):
            for k in range(nzmax-nzmin+1):
                if axis == 0: 
                    # X Slice
                    plotslice[i] += plotdata[nxmin+i,nymin+j,nzmin+k]                
                elif axis == 1: 
                    # Y Slice
                    plotslice[i] += plotdata[nymin+j,nxmin+i,nzmin+k]                
                elif axis == 2: 
                    # Z Slice
                    plotslice[i] += plotdata[nymin+j,nzmin+k,nxmin+i]                
    return [plotslice, xpoints]

def ChargeDepth(dat, filename, nxcenter, nycenter, dnx, dny, nzmax, recomb_factor = 1.0):
    # Calculates the charges in a region and the average depth of the electron cloud
    # The recomb_factor removes a fraction of electrons which have reached the interface
    if filename is not  None:
        file = open(filename,"w")
    nzxmin = nxcenter - dnx
    nzxmax = nxcenter + dnx
    nzymin = nycenter - dny
    nzymax = nycenter + dny
    ncenter = 0.0
    nzcenter = 0.0
    ncenter_int = 0.0
    # Find Ckmin
    for nz in range(nzmax):
        if dat.rho[nxcenter, nycenter, nz] > 1.0E-12:
            Ckmin = nz
            break
    Holes_in_region = 0.0    
    Below_kmin = 0.0
    for nx in range(nzxmin, nzxmax+1):
        for ny in range(nzymin, nzymax+1):
            for nz in range(nzmax):
                if dat.Elec[nx,ny,nz] > 0.1:
                    if nz == Ckmin:
                        ncenter_int += dat.Elec[nx,ny,nz]
                        ncenter += dat.Elec[nx,ny,nz] * recomb_factor
                        nzcenter += dat.z[nz] * dat.Elec[nx,ny,nz] * recomb_factor
                    else:
                        ncenter += dat.Elec[nx,ny,nz]
                        nzcenter += dat.z[nz] * dat.Elec[nx,ny,nz]
                if dat.Hole[nx,ny,nz] > 0.1:
                    Holes_in_region += dat.Hole[nx,ny,nz]                    
    if ncenter > 0:
        meanz = nzcenter / ncenter
    else:
        meanz = 0.0
    print("Electrons at Interface = %.1f\n"%ncenter_int)
    print("Electrons In Region = %.1f, Mean z = %.3f microns\n"%(ncenter, meanz))
    print("Holes In Region = %.1f\n"%Holes_in_region)
    print("Total Electrons = %.1f\n"%dat.Elec.sum())
    if filename is not  None:
        file.write("Electrons in Region = %.1f, Mean z = %.3f microns\n"%(ncenter, meanz))
        file.close()
    return ncenter

def ReadCorrData(filename):
    # This reads the correlation data file
    # and returns arrays with the values and sigmas
    systematic_error = 0.0001
    data = np.zeros([6,6])
    sigma = np.zeros([6,6])
    file = open(filename,'r')
    lines=file.readlines()
    file.close()
    lines.remove(lines[0]) # Strip the title line
    try:
        for line in lines:
            items = line.split()
            if items[0] == 'Slope':
                break
            i = int(items[0])
            j = int(items[1])        
            data[i,j] = float(items[2])
            sigma[i,j] = float(items[3]) + systematic_error        
    except:
        print("Error reading data file")
        sys.exit()
    return [data, sigma]

def ReadAreaFile(filename, nx, ny, nxcenter, nycenter, Area_0):
    # This reads the correlation data file
    # and returns an array with the expected correlations
    area = np.zeros([nx, ny])
    file = open(filename, 'r')
    lines = file.readlines()
    file.close()
    lines.remove(lines[0]) # Strip the title line    
    for line in lines:
        items = line.split()
        i = int(items[0])
        j = int(items[1])
        area[i,j] = float(items[2])

    sim = np.zeros([6,6])
    num = np.zeros([6,6], dtype = int)    
    try:
        for i in range(nx):
            for j in range(ny):
                ii = abs(i - nxcenter)
                jj = abs(j - nycenter)
                sim[jj,ii] += (area[i,j] - Area_0) / Area_0
                num[jj,ii] += 1
        for i in range(6):
            for j in range(6):
                if num[i,j] >0:
                    sim[i,j] /= float(num[i,j])
    except:
        pass
    return [area,sim]

def ReadVertexFile(filename, nx, ny, NumAngles):
    vx = np.zeros([nx, ny, NumAngles])
    vy = np.zeros([nx, ny, NumAngles])
    file = open(filename, 'r')
    lines = file.readlines()
    file.close()
    FirstPass = True
    for line in lines:
        items = line.split()
        if items[0] == "X0":
            continue
        x = int(float(items[0]))
        y = int(float(items[1]))
        if FirstPass:
            lastx = x; lasty = y; i = 0; j = 0; k = 0; FirstPass = False
        if y != lasty:
            k = 0
            j += 1
            lasty = y
        if x != lastx:
            k = 0
            j = 0
            i += 1
            lastx = x
            lasty = y

        #print(i,j,k)
        vx[i,j,k] = float(items[3])
        vy[i,j,k] = float(items[4])
        k += 1
    return (vx, vy)

def ReadCCFile(filename, Nx, Ny):
    file = open(filename,'r')
    lines=file.readlines()
    file.close()
    lines.remove(lines[0]) # Strip the title line
    elec = np.zeros([Nx,Ny])
    count = 0
    for line in lines:
        items = line.split()
        i = int(items[0]) - 1
        j = int(items[1]) - 1        
        elec[i,j] = int(items[2])

    return elec


def New_BF_Cfg_File(incfgfile, outcfgfile, newrun):
    # This increments the run number to start a new spot
    # and assigns a random offset value within the central pixel
    ConfigData = ReadConfigFile(incfgfile)
    lines = OpenFile(incfgfile)
    dirbase = ConfigData['outputfiledir'].split('bfrun')
    xoff = -5.0 + 10.0 * rand()
    yoff = -5.0 + 10.0 * rand()
    for i, line in enumerate(lines):
        try:
            if line.split()[0] == 'Xoffset':
                lines[i] = 'Xoffset = '+str(xoff)+'\n'
        except:
            continue
        try:
            if line.split()[0] == 'Yoffset':
                lines[i] = 'Yoffset = '+str(yoff)+'\n'
        except:
            continue
        try:
            if line.split()[0] == 'outputfiledir':
                lines[i] = 'outputfiledir = '+dirbase[0]+'bfrun_%d\n'%newrun
        except:
            continue
    newcfgfile = open(outcfgfile, 'w')
    for line in lines:
        newcfgfile.write(line)
    newcfgfile.close()
    return 

def New_Fe55_Cfg_File(incfgfile, outcfgfile, newrun):
    # This increments the run number to start a new spot
    # and assigns a random offset value within the central pixel
    ConfigData = ReadConfigFile(incfgfile)
    lines = OpenFile(incfgfile)
    dirbase = ConfigData['outputfiledir'].split('run')
    newseed = int(2**32 * rand())
    for i, line in enumerate(lines):
        try:
            if line.split()[0] == 'Seed':
                lines[i] = 'Seed = %d'%newseed
        except:
            continue
        try:
            if line.split()[0] == 'outputfiledir':
                lines[i] = 'outputfiledir = '+dirbase[0]+'run_%d\n'%newrun
        except:
            continue
    newcfgfile = open(outcfgfile, 'w')
    for line in lines:
        newcfgfile.write(line)
    newcfgfile.close()
    return 

def New_Trans_Cfg_File(incfgfile, outcfgfile, newrun, Vg):
    # This increments the run number and the gate voltage
    # to simulate a transistor I-V curve
    ConfigData = ReadConfigFile(incfgfile)
    lines = OpenFile(incfgfile)
    dirbase = ConfigData['outputfiledir'].split('transrun')
    for i, line in enumerate(lines):
        try:
            if line.split()[0] == 'FixedRegionVoltage_11':
                lines[i] = 'FixedRegionVoltage_11 = '+str(Vg)+'\n'
        except:
            continue
        try:
            if line.split()[0] == 'FixedRegionVoltage_36':
                lines[i] = 'FixedRegionVoltage_36 = '+str(Vg)+'\n'
        except:
            continue
        try:
            if line.split()[0] == 'outputfiledir':
                lines[i] = 'outputfiledir = '+dirbase[0]+'transrun_%d\n'%newrun
        except:
            continue
    newcfgfile = open(outcfgfile, 'w')
    for line in lines:
        newcfgfile.write(line)
    newcfgfile.close()
    return 

def New_Sat_Cfg_File(incfgfile, outcfgfile, newrun, Vph, Vpl):
    # This increments the run number, the parallel voltages
    # and the charges to simulate a saturationcurve
    qmax = int(15000 * (Vph - Vpl))
    qmin = qmax / 4
    qstep = qmax * 15 / 100

    ConfigData = ReadConfigFile(incfgfile)
    lines = OpenFile(incfgfile)
    dirbase = ConfigData['outputfiledir'].split('satrun')
    for i, line in enumerate(lines):
        try:
            if line.split()[0] == 'Vparallel_lo':
                lines[i] = 'Vparallel_lo = '+str(Vpl)+'\n'
        except:
            continue
        try:
            if line.split()[0] == 'Vparallel_hi':
                lines[i] = 'Vparallel_hi = '+str(Vph)+'\n'
        except:
            continue
        try:
            if line.split()[0] == 'outputfiledir':
                lines[i] = 'outputfiledir = '+dirbase[0]+'satrun_%d\n'%newrun
        except:
            continue
        # Now set the 6 charges to give a reasonable BarrierHeight curve
        for m in range(6):
            for n in range(3):
                index = 3 * m + n
                charge = qmin + m * qstep
                try:
                    if line.split()[0] == 'CollectedCharge_0_%d'%index:
                        lines[i] = 'CollectedCharge_0_%d = %d\n'%(index, charge)
                except:
                    continue

    newcfgfile = open(outcfgfile, 'w')
    for line in lines:
        newcfgfile.write(line)
    newcfgfile.close()
    return 


def OpenFile(filename):
    Numtries = 10
    tries = 0
    lines = []
    while tries < Numtries:
        try:
            file = open(filename, 'r')
            lines = file.readlines()
            file.close()
            break
        except:
            time.sleep(1.0)
            tries += 1

    return lines

def ReadSIMSData():

    epsilon_si = 11.7 * 8.85E-14
    qe = 1.6E-19
    boron_wb = xlrd.open_workbook('measurements/C0HVL528x02_B site.xls')
    boron_data = boron_wb.sheet_by_name('Processed data')
    boron_conc = []
    boron_depth = []
    for i in range(boron_data.nrows):
        try:
            if type(boron_data.row(i)[0].value) is float:
                if boron_data.row(i)[1].value < 1.0E17: # One data point is garbage
                    boron_depth.append(boron_data.row(i)[0].value)
                    boron_conc.append(boron_data.row(i)[1].value)
        except:
            continue

    phos_wb = xlrd.open_workbook('measurements/C0HVL528L05_P Site.xls')
    phos_data = phos_wb.sheet_by_name('Processed data')
    phos_conc = []
    phos_depth = []
    for i in range(phos_data.nrows):
        try:
            if type(phos_data.row(i)[0].value) is float:
                phos_depth.append(phos_data.row(i)[0].value)
                phos_conc.append(phos_data.row(i)[1].value)
        except:
            continue

    boron_conc = -np.array(boron_conc) * qe / epsilon_si * 1.0E-8 # Convert to V/um/2
    phos_conc = np.array(phos_conc) * qe / epsilon_si * 1.0E-8 # Convert to V/um/2
    boron_depth = np.array(boron_depth)
    phos_depth = np.array(phos_depth)    
    return [boron_depth, boron_conc, phos_depth, phos_conc]


def mu_si(E, T):
    # Electron mobility
    vm = 1.53E9 * pow(T, -0.87)
    Ec = 1.01 * pow(T, 1.55)
    beta = 2.57E-2 * pow(T, 0.66)
    return ((vm/Ec) / pow(1.0 + pow(abs(E)/Ec,beta), 1/beta))

def Read_STA3800_IV_Data(filename):
    data_wb = xlrd.open_workbook(filename)
    idvg1_data = data_wb.sheet_by_name('8-1011.TXT')
    idvg3_data = data_wb.sheet_by_name('8-1014.TXT')
    Vgs = []
    Ids = []
    for i in range(idvg1_data.nrows):
        try:
            if type(idvg1_data.row(i)[0].value) is float:
                #if idvg1_data.row(i)[1].value > 0.0:
                Vgs.append(idvg1_data.row(i)[1].value)
                Ids.append(idvg1_data.row(i)[3].value)
        except:
            continue
    return [Vgs, Ids]


def Area(xl, xh, yl, yh, sigmax, sigmay, Imax):
    # Calculates how much of a 2D Gaussian falls within a rectangular box
    ssigx = np.sqrt(2) * sigmax
    ssigy = np.sqrt(2) * sigmay    
    I = (erf(xh/ssigx)-erf(xl/ssigx))*(erf(yh/ssigy)-erf(yl/ssigy))
    return Imax * I / 4.0

