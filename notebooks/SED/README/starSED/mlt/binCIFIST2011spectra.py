import os
import bz2
import gzip
import numpy
from lsst.sims.photUtils import Bandpass, Sed

inputDir = '/astro/store/pogo3/danielsf/CIFIST2011raw'
outputDir = '/astro/store/pogo3/danielsf/CIFIST2011binned'

normBandpass = Bandpass()
normBandpass.imsimBandpass()

inNames = []
for name in os.listdir(inputDir):
    if name[:3]=='lte':
        inNames.append(name)

dwav = 0.5 #bin width in nm
wavMax = 3.1e4 #maximum wavelength in nm

for filename in inNames:
    outname = filename.strip('bz2')
    outname = outname.strip('.7')
    outname = outname + '.gz'
    outname = os.path.join(outputDir, outname)
    if not os.path.exists(outname):
        name = os.path.join(inputDir,filename)

        rawArray = []

        inSpectrum = bz2.BZ2File(name,'r')
        for line in inSpectrum:
            vv = line.split()
            rawArray.append((numpy.float(vv[0].replace('D','E'))*0.1, numpy.float(vv[1].replace('D','E'))))
        inSpectrum.close()

        rawArray = numpy.array(rawArray, dtype=numpy.dtype([('wavelen',numpy.float),('flux',numpy.float)]))

        rawArray = numpy.sort(rawArray, order='wavelen')

        rawWav = rawArray['wavelen']
        rawFlux = numpy.power(10.0, rawArray['flux'])

        if numpy.diff(rawWav).min()<0.0:
            print 'WARNING diff min in wav ',numpy.diff(rawWav).min()
            exit()

        remainder1 = numpy.abs((rawWav % dwav))
        remainder = numpy.where(remainder1<0.5*dwav, remainder1, dwav-remainder1)

        dexes = numpy.array([ii for ii in range(len(remainder)) if rawWav[ii] <= wavMax and (ii==0 or ii==len(remainder)-1 or \
                                                 (remainder[ii]<remainder[ii-1] and remainder[ii]<remainder[ii+1]) \
						 or rawWav[ii]-rawWav[ii-1]>dwav)])

        minDexes = dexes[:len(dexes)-1]
        maxDexes = dexes[1:]+1
        wavOut = rawWav[minDexes]

        fluxOut = numpy.array([
                          numpy.trapz(rawFlux[minDexes[ii]:maxDexes[ii]], rawWav[minDexes[ii]:maxDexes[ii]])/numpy.diff(rawWav[minDexes[ii]:maxDexes[ii]]).sum()
                          for ii in range(len(wavOut))
                          ])

	
	tempSed = Sed(wavelen=wavOut, flambda=fluxOut)
        ff = tempSed.calcFluxNorm(15.0, normBandpass)
	fluxOut = fluxOut*ff

        output = gzip.open(outname, 'w')
        for w, f in zip(wavOut, fluxOut):
            output.write('%.7e %.7e\n' % (w,f))
        output.close()
        print 'wrote ',outname,fluxOut.min(),fluxOut.max()





