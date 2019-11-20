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

def ReadCorrDataBig(filename):
    # This reads the correlation data file
    # and returns arrays with the values and sigmas
    systematic_error = 0.0
    data = np.zeros([8,8])
    sigma = np.zeros([8,8])
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


#****************MAIN PROGRAM*****************

# First, read the .cfg file

configfile = sys.argv[1]
run = int(sys.argv[2])
ConfigData = ReadConfigFile(configfile)

outputfilebase = ConfigData["outputfilebase"]
outputfiledir = ConfigData["outputfiledir"]
Nx = ConfigData["PixelBoundaryNx"]
Ny = ConfigData["PixelBoundaryNy"]
NumElec = ConfigData["CollectedCharge_0_0"]
NumPhases = ConfigData["NumPhases"]

NxCenter = 4
NyCenter = 4
Area_0 = 99.9974

areafilename = outputfiledir + '/' + outputfilebase +'_%d_Area.dat'%run
if NumPhases == 3:
    datafilename = 'measurements/corr_itl_08nov19.txt'
elif NumPhases == 4:
    datafilename = 'measurements/corr_e2v_08nov19.txt'
[area,sim] = ReadAreaFile(areafilename, Nx, Ny, NxCenter, NyCenter, Area_0)
[data,sigma] = ReadCorrDataBig(datafilename)


x = []
y = []
xfit = []
yfit = []
xneg = []
yneg = []
xdata = []
ydata = []
yerr = []
ysigma = []
xdataneg = []
ydataneg = []
yerrneg = []
ysigmaneg = []

FOM = 0.0
Num_FOM = 0.0

for i in range(6):
    for j in range(6):
        yvalue = sim[i,j] / NumElec
        if i + j == 0:
            rsquared = 0.85
        else:
            rsquared = float(i**2 + j**2)
        #if rsquared > 18:
        #    continue
        FOM_contrib = ((data[i,j] - yvalue) / sigma[i,j])**2
        print("i = %d, j = %d, data = %.6g, sim = %.6g, FOM_contrib = %.4f"%(i,j,data[i,j],yvalue,FOM_contrib))
        FOM += FOM_contrib
        Num_FOM += 1.0
        if yvalue < 0:
            xneg.append(rsquared)
            yneg.append(-yvalue)
        else:
            x.append(rsquared)
            y.append(yvalue)
        yvalue = data[i,j]
        if yvalue < 0:
            xdataneg.append(rsquared)
            ydataneg.append(-yvalue)
            yerrneg.append(sigma[i,j])
        else:
            xdata.append(rsquared)
            ydata.append(yvalue)
            yerr.append(sigma[i,j])
        if rsquared > 1.1 and i < 3 and j < 3 and yvalue > 0:
            xfit.append(rsquared)
            yfit.append(yvalue)


chisquared = FOM / Num_FOM
print("Chi-squared = %.2f"%chisquared) 
plt.figure()
plt.title("Covariance Matrix", fontsize=18)
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.8, 100.0)
plt.ylim(1E-10, 1E-5)
yerr = np.array(yerr)
ylower = np.maximum(1.1E-10, ydata - yerr)
yerr_lower = ydata - ylower
yerrneg = np.array(yerrneg)
ylowerneg = np.maximum(1.1E-10, ydataneg - yerrneg)
yerr_lowerneg = ydataneg - ylowerneg
plt.errorbar(xdata,ydata, yerr = [yerr_lower, 2.0*yerr] , ls = 'None',marker = '.', ms = 10, color = 'green', label = 'Data')
plt.scatter(np.array(x), np.array(y), marker = 'x', s = 50, color = 'blue', label = 'Sim')
if len(xdataneg) > 0:
    plt.errorbar(xdataneg,ydataneg, yerr = [yerr_lowerneg, 2.0*yerrneg] , ls = 'None',marker = '.', ms = 10, color = 'cyan', label = 'Data-Neg')
if len(xneg) > 0:
    plt.scatter(np.array(xneg), np.array(yneg), marker = 'x', s = 50, color = 'magenta', label = 'Sim-Neg')

from scipy import stats
slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(xfit),np.log10(yfit))
xplot=np.linspace(0.0, 2.0, 100)
yplot = slope * xplot + intercept
plt.plot(10**xplot, 10**yplot, color='red', lw = 2, ls = '--')
plt.text(1.2, 5E-6, "C00: Meas = %.4g, Sim = %.4g"%(data[0,0],sim[0,0]/NumElec),fontsize = 10)
plt.text(1.2, 2.5E-6, "C01: Meas = %.4g, Sim = %.4g"%(data[0,1],sim[0,1]/NumElec),fontsize = 10)
plt.text(1.2, 1.25E-6, "C10: Meas = %.4g, Sim = %.4g"%(data[1,0],sim[1,0]/NumElec),fontsize = 10)
plt.text(1.2, 6.25E-7, "C11: Meas = %.4g, Sim = %.4g"%(data[1,1],sim[1,1]/NumElec),fontsize = 10)
plt.text(1.2, 3.10E-7, "$\chi^2$/DOF = %.2f"%chisquared, fontsize = 10)
plt.xlabel("$i^2 + j^2$")
plt.ylabel("$\delta$ Area / Area")
plt.xlim(0.8,100)
plt.legend()
plt.savefig(outputfiledir+"/plots/Area_Covariance_%d.pdf"%run)
