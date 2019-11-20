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
import pickle

#****************MAIN PROGRAM*****************

# First, read the .cfg file

configfile = sys.argv[1]
run = int(sys.argv[2])
ConfigData = ReadConfigFile(configfile)
mainfiledir = ConfigData["outputfiledir"]

# Now get the measured data

pkl_file_20170509 = open('measurements/sat_data_20170509.pkl', 'rb')
data = pickle.load(pkl_file_20170509, encoding='latin1')
pkl_file_20170509.close()
[[VphVplm10,ElecVplm10],[VphVplm9,ElecVplm9],[VphVplm8,ElecVplm8],[VphVplm6,ElecVplm6],[VphVplm4,ElecVplm4],[VphVplm2,ElecVplm2]] = data

#****************MAIN PROGRAM*****************
NumSats = int(sys.argv[3])
SatSims = []
colors = {'-10.0':'red', '-9.0':'blue', '-8.0':'green', '-6.0': 'black', '-4.0':'magenta', '-2.0':'cyan'}
for satrun in range(NumSats):
    configfile = 'data/satrun_%d/sat.cfg'%satrun
    run = 0
    ConfigData = ReadConfigFile(configfile)

    Vph = ConfigData["Vparallel_hi"]
    Vpl = ConfigData["Vparallel_lo"]
    outputfiledir = ConfigData["outputfiledir"]
    barrierfile = open(outputfiledir+'/barrier.txt', 'r')
    lines = barrierfile.readlines()
    barrierfile.close()
    barrier_height = float(lines[-1].split()[4])
    SatSims.append([Vph, Vpl, barrier_height])

# Create the output directory if it doesn't exist
if not os.path.isdir(mainfiledir+"/plots"):
    os.mkdir(mainfiledir+"/plots")

plt.figure()
plt.suptitle("Saturation Curve - ITL STA3800 - 029")
plt.subplots_adjust(wspace = 2.0)
plt.subplot(1,1,1)
plt.xlabel("Vparallel-hi (V)")
plt.ylabel("Peak Signal(Electrons)")
plt.ylim(0,600000)
plt.xlim(-1.0,6.0)

plt.plot(VphVplm2, ElecVplm2, marker = 'x', label = 'Meas:Vpl = -2V', color = colors['-2.0'])
plt.plot([entry[0] for entry in SatSims[0:6]], [entry[2] for entry in SatSims[0:6]], marker = '*', ls = 'None', color = colors['-2.0'], markersize = 20)

plt.plot(VphVplm4, ElecVplm4, marker = 'x', label = 'Meas:Vpl = -4V', color = colors['-4.0'])
plt.plot([entry[0] for entry in SatSims[6:12]], [entry[2] for entry in SatSims[6:12]], marker = '*', ls = 'None', color = colors['-4.0'], markersize = 20)

plt.plot(VphVplm6, ElecVplm6, marker = 'x', label = 'Meas:Vpl = -6V', color = colors['-6.0'])
plt.plot([entry[0] for entry in SatSims[12:18]], [entry[2] for entry in SatSims[12:18]], marker = '*', ls = 'None', color = colors['-6.0'], markersize = 20)

plt.plot(VphVplm8, ElecVplm8, marker = 'x', label = 'Meas:Vpl = -8V', color = colors['-8.0'])
plt.plot([entry[0] for entry in SatSims[18:24]], [entry[2] for entry in SatSims[18:24]], marker = '*', ls = 'None', color = colors['-8.0'], markersize = 20)

plt.legend(loc = 'upper right')
plt.text(-0.8,500000, "Simulations are large stars")
plt.savefig(mainfiledir+"/plots/Saturation_Plot.pdf")


#************END MAIN PROGRAM*************************


