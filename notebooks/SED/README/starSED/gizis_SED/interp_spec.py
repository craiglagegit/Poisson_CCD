import scipy.interpolate as si
import numpy
import sys

fh = open(sys.argv[1]+".dat", 'r')
convFlux = bool(int(sys.argv[2]))
fhout = open(sys.argv[1]+"_interp.dat", 'w')
waves = []
fluxs = []
for l in fh:
    flds = l.rstrip().split()
    waves.append(float(flds[0])/10.)
    if convFlux:
        print "doing flux conversion"
        fluxs.append(float(flds[1])/waves[-1])
    else:
        fluxs.append(float(flds[1]))
fh.close()
waves = numpy.asarray(waves)
fluxs = numpy.asarray(fluxs)

intrp = si.interp1d(waves, fluxs)
norm = intrp(500.)
iwaves = numpy.array(range(int(waves.max() - waves.min()))) + int(waves.min())
ifluxs = intrp(iwaves)/norm
for x,y in zip(iwaves,ifluxs):
    fhout.write("%f %f\n"%(x,y))
fhout.close()
