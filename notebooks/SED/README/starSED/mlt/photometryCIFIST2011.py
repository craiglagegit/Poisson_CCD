import eups
import os
import numpy
from lsst.sims.photUtils import PhotometryBase
from lsst.sims.photUtils import Sed, Bandpass

bandpassDir = os.path.join(eups.productDir('throughputs'),'sdss')

#directory where the resampled SEDs are stored
directory = '/astro/store/pogo3/danielsf/CIFIST2011binned'
fileNames = os.listdir(directory)

controlBandpass = Bandpass()
controlBandpass.imsimBandpass()

sdssPhot = PhotometryBase()
sdssPhot.loadTotalBandpassesFromFiles(bandpassNames=['u','g','r','i','z'],
                                      bandpassDir=bandpassDir, bandpassRoot='sdss_')

output = open('MLTlookup.txt','w')

for name in fileNames:
    if name[:3] == 'lte':
        sedName = os.path.join(directory,name)
        sdssSED = Sed()
        sdssSED.readSED_flambda(sedName)
        if sdssSED.wavelen[-1]<30000.:
            print 'WARNING ',name,' max wave ',sdssSED.wavelen[-1]
        else:
            lsstSED = Sed()
            lsstSED.readSED_flambda(sedName)
            normSED = Sed()
            normSED.readSED_flambda(sedName)
            sdssMags = sdssPhot.manyMagCalc_list(sdssSED)
            magNorm = normSED.calcMag(controlBandpass)
    
            output.write(name+', ')
            for m in sdssMags:
                output.write('%.12le, ' % m)
            output.write('%.12le ' % magNorm)
            output.write('\n')
    
output.close()
