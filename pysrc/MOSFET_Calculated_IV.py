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

configfile = sys.argv[1]
run = int(sys.argv[2])
ConfigData = ReadConfigFile(configfile)
outputfilebase = ConfigData["outputfilebase"]
outputfiledir = ConfigData['outputfiledir']
dirbase = ConfigData['outputfiledir'].split('transrun')[0]


q = 1.6E-19   # MKS
me = 9.11E-31   # MKS
k = 1.38E-23   # MKS
T = 273.0
E = 20000.0
Tox = ConfigData["GateOxide"] * 1.0E-4
Eps_ox = 8.85E-14 * 4.2
mu = mu_si(E, T)
print(("Mobility = %f"%mu))
Vds = 0.5
L = 5.0E-4
Vgates = []
Is = []
for run in range(13):
    newcfgfile = dirbase+'transrun_%d/trans.cfg'%run
    ConfigData = ReadConfigFile(newcfgfile)

    Vgate = ConfigData["FixedRegionVoltage_11"]
    file = open(dirbase+'transrun_%d/charge.txt'%run,'r')
    lines = file.readlines()
    file.close
    Ne = float(lines[0].split()[4].strip(','))
    print(run, Vgate, Ne)
    I = Vds * q * mu * Ne / L**2
    Vgates.append(Vgate)
    Is.append(I)

Qss = ConfigData["ChannelSurfaceCharge"]
[Vgs, Ids] = Read_STA3800_IV_Data("measurements/STA3800_meas.xls")

#Delta_V = -0.7
#Delta_Q = Eps_ox / Tox * Delta_V / q

# Create the output directory if it doesn't exist
if not os.path.isdir(outputfiledir+"/plots"):
    os.mkdir(outputfiledir+"/plots")

plt.figure()
ax1=plt.axes([0.2,0.1,0.6,0.6])
ax1.set_title("STA3800C Id-Vg")
ax1.plot(Vgs, np.array(Ids)*1000.0, color = 'green', lw = 2, label='Measured')
ax1.scatter(np.array(Vgates), np.array(Is)*1000.0, marker = 'x', color='red', label='Sim - Qss = %g'%Qss)
#ax1.scatter(np.array(Vgates)+Delta_V, np.array(Is)*1000.0, marker = 'x', color='red', label='Sim - QS=1.4E12')
ax1.set_xlabel("Vgs (volts)")
ax1.set_ylabel("Ids(mA)")
ax1.set_ylim(0,1.0)
ax1.legend(loc='upper left')
ax1.text(-23, 0.3, 'Vds = 0.5V')
#ax1.text(-23, 0.5, 'DeltaV = %.1f, DeltaQ = %g'%(Delta_V, Delta_Q))
plt.savefig(outputfiledir+"/plots/IdVg_QSS_%s.pdf"%str(Qss))
