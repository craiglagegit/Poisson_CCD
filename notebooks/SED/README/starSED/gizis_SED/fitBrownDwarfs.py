from __future__ import with_statement
import os
import numpy
import gzip
import eups
from lsst.sims.photUtils import Bandpass, Sed, PhotometryBase

lsstDir = os.path.join(eups.productDir('sims_sed_library'),'starSED','gizis_SED')
burrowsDir = '/Users/danielsf/physics/burrowsBrownDwarfs'
btSettlDir = '/Users/danielsf/physics/newSedLibrary/starSED/mlt/'

sdssPhotometry = PhotometryBase()
sdssPhotometry.loadTotalBandpassesFromFiles(bandpassNames = ['u','g','r','i','z'],
                                            bandpassDir = os.path.join(eups.productDir('throughputs'),'sdss'),
                                            bandpassRoot = 'sdss_')

controlBandpass = Bandpass()
controlBandpass.imsimBandpass()

linesDict = {}
with open('linesLookup.txt', 'r') as linesInput:
    for line in linesInput:
        vv = line.split()
        if vv[0] =='>':
            linesDict[vv[1]] = numpy.float(vv[2])

listOfSubDirs = os.listdir(burrowsDir)

dtype = numpy.dtype([
                    ('id', numpy.int),
                    ('freq', numpy.float),
                    ('wavelen', numpy.float),
                    ('fnu', numpy.float),
                    ('flambda', numpy.float),
                    ('fdet', numpy.float)
                    ])


_fitNames = []
_fitWavelen = []
_fitFlux = []
_wavStart = None

for dirName in listOfSubDirs:
    files = os.listdir(os.path.join(burrowsDir,dirName))
    for name in files:
        if 'tar' not in name:
            fullName = os.path.join(burrowsDir,dirName,name)
            _fitNames.append(os.path.join(dirName,name))
        
            data = numpy.loadtxt(fullName, dtype=dtype, skiprows=1)
            data['wavelen'] *= 1000.0
            
            _fitWavelen.append(data['wavelen'])
            _fitFlux.append(data['flambda'])
        
        
            if _wavStart is None:
                _wavStart = data['wavelen'][0]
            else:
                if data['wavelen'][0] != _wavStart:
                    print 'WARNING wavstart ',_wavStart,data['wavelen'][0]
                    exit()


print 'read in burrows spectra'

dtype=numpy.dtype([
                  ('wavelen', numpy.float),
                  ('flux', numpy.float)
                  ])


for name in os.listdir(btSettlDir):
    if name[:3]=='lte':
        fullName = os.path.join(btSettlDir,name)
        _fitNames.append(fullName)
        data = numpy.loadtxt(fullName, dtype=dtype)
        _fitWavelen.append(data['wavelen'])
        _fitFlux.append(data['flux'])
 

print 'read in btSettl'

lsstNamesRaw = os.listdir(lsstDir)
lsstNames = []
lsstWavelen = []
lsstFlux = []
for name in lsstNamesRaw:
    if name[:2] == 'BD' and 'interp' not in name:
        lsstNames.append(name)
        data = numpy.loadtxt(os.path.join(lsstDir,name), dtype=dtype)
        data['wavelen'] *= 0.1

        lsstWavelen.append(data['wavelen'])
        lsstFlux.append(data['flux'])

lsstWavelen = lsstWavelen
lsstFlux = lsstFlux

print 'read in lsst spectra'

def bestFit(wavelen, flux):
    """
    wavelen and flux are the arrays to be fit (the test spectrum)
    
    fitWavelen and fitFlux are 2d numpy arrays of the spectra you are
    choosing from
    """
    
    dexes = numpy.where(wavelen>_wavStart)
    wavs = wavelen[dexes]
    ff = flux[dexes]
    
    distance = []
    normStorage = []
    normedDistance = []

    for fitwavs, fitff, name in zip(_fitWavelen, _fitFlux, _fitNames):
        fff = numpy.interp(wavs, fitwavs, fitff)
        norm = (fff*ff).sum()/(fff*fff).sum()
        if numpy.isnan(norm):
            print 'WARNING norm is nan ',(fff*ff).sum(), (fff*fff).sum(), name
            for x,y,z in zip(wavs, ff, fff):
                print x, y, z
            exit()
        normStorage.append(norm)
        normedDistance.append(numpy.power(ff-fff*norm,2).sum())
        dd = numpy.power(ff-fff*norm,2).sum()
        distance.append(dd)
    distance = numpy.array(distance)
    
    dex = numpy.argmin(distance)
    return dex, normStorage[dex]

with open('burrowsMapping.txt','w') as output:
    for name, wavelen, flux in zip(lsstNames, lsstWavelen, lsstFlux):
        if 'BD1000e' not in name:
            dex, norm = bestFit(wavelen, flux)

            outName = name.replace('.dat','_interp.dat')

            wavelen = numpy.copy(_fitWavelen[dex])
            flux = numpy.copy(_fitFlux[dex])
            flux *= norm
            
            if wavelen[0]>9.0:
                coeffs = numpy.polyfit(numpy.log10(wavelen[:100]), numpy.log10(flux[:100]), deg=1)
                redWavelen = numpy.arange(9.0, wavelen[0]-0.5, 0.4*(wavelen[0]-0.5-9.0))
                redFlux = numpy.ones(len(redWavelen))*10e-5
                wavelen = numpy.append(redWavelen, wavelen)
                flux = numpy.append(redFlux, flux)

            with gzip.open(outName,'wb') as spectralFile:
                for w,f in zip(wavelen, flux):
                    spectralFile.write('%e %e\n' % (w,f))


            dummySed = Sed(wavelen=wavelen, flambda=flux)
            sdssMags = sdssPhotometry.manyMagCalc_list(dummySed)
            dummySed = Sed(wavelen=wavelen, flambda=flux)
            magNorm = dummySed.calcMag(controlBandpass)
            output.write('%s is %s -- norm %e -- %e %e %e %e %e %e\n' % 
                        (name, _fitNames[dex], norm,
                         sdssMags[0], sdssMags[1], sdssMags[2], sdssMags[3], sdssMags[4],
                         magNorm))

            if 'BD1000' in outName:
                outName = outName.replace('BD1000','BD1000e')
                for wavelenStr in linesDict:
                    dex = numpy.argmin(numpy.abs(wavelen-0.1*numpy.float(wavelenStr)))
                    flux[dex] = linesDict[wavelenStr]
                with gzip.open(outName,'wb') as spectralFile:
                    for w, f in zip(wavelen, flux):
                        spectralFile.write('%e %e\n' % (w,f))
