from __future__ import with_statement
import numpy
import os
import gzip
import eups

rootFile = os.path.join(eups.productDir('sims_sed_library'),'starSED','gizis_SED','PNastrom.dat.gz')
outputFile = 'PNastrom_interp.dat.gz'

dtype = numpy.dtype([('wavelen',numpy.float),
                     ('flux',numpy.float)
                     ])


data = numpy.loadtxt(rootFile,dtype=dtype)
data['wavelen'] *= 0.1

dwav=10000.0 #in nm
redLimit = 30000.0 #30 microns in nm
redWavelen = numpy.arange(data['wavelen'][-3]+5.0,redLimit,dwav)

coeffs = numpy.polyfit(numpy.log10(data['wavelen'][-20:-3]), numpy.log10(data['flux'][-20:-3]), deg=1)
redFlux = numpy.power(10.0, numpy.polyval(coeffs, numpy.log10(redWavelen)))

wavelen = numpy.append(data['wavelen'][:-3], redWavelen)
flux = numpy.append(data['flux'][:-3], redFlux)

norm = numpy.interp(500.0, wavelen, flux)
flux=flux/norm

with gzip.open(outputFile, 'wb') as output:
    for w, f in zip(wavelen, flux):
        output.write('%e %e\n' % (w,f))

