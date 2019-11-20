import os
import eups
import numpy
import gzip

fileNames = ['RedE1astrom.dat.gz', 'RedE2astrom.dat.gz']
rootDir = os.path.join(eups.productDir('sims_sed_library'),'starSED','gizis_SED')

dtype = numpy.dtype([('wavelen', numpy.float),
                     ('flux', numpy.float)])

for name in fileNames:
    fullName = os.path.join(rootDir,name)
    data = numpy.loadtxt(fullName, dtype=dtype)
    data['wavelen'] *= 0.1
    
    wavelen = numpy.arange(9.0,data['wavelen'][0]-0.5,0.1*(data['wavelen'][0]-0.5-9.0))
    flux = numpy.ones(len(wavelen))*data['flux'][0]
    
    wavelen = numpy.append(wavelen, data['wavelen'])
    flux = numpy.append(flux, data['flux'])
    
    blueWavelen = numpy.arange(data['wavelen'][-1]+0.5, 30000.0, 0.1*(30000.0 - data['wavelen'][-1] -0.5))
    blueFlux = numpy.ones(len(blueWavelen))*data['flux'][-1]
    
    wavelen = numpy.append(wavelen, blueWavelen)
    flux = numpy.append(flux, blueFlux)
    
    flux = flux/wavelen
    
    norm = numpy.interp(500.0, wavelen, flux)
    flux=flux/norm
    
    
    with gzip.open(name.replace('.dat','_interp.dat'),'wb') as output:
        for w, f in zip(wavelen, flux):
            output.write('%e %e\n' % (w,f))
    
    
    
    
