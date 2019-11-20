import os
import numpy
import gzip


def parseBergeron(logg, desiredTemperatures, chemistry='H',
                  outputDir = 'extendedBergeron1503'):

    bergeronDirectory = 'baseBergeronModels'
    
    if chemistry == 'He':
        suffix = '_300_ML125_IR'
        inName = os.path.join(bergeronDirectory, logg+'0'+suffix)
    else:
        inName = os.path.join(bergeronDirectory, logg+'0')
    
    inFile = open(inName, 'r')
    lines = inFile.readlines()
    inFile.close()
    nwav = -1
    wavelen = []
    
    iterator = lines.__iter__()
    nwav = int(iterator.next().strip())
    wavelen = []
    while len(wavelen)<nwav:
        line = iterator.next()
        line = line.split()
        for vv in line:
            wavelen.append(float(vv))
    if len(wavelen) != nwav:
        print 'WARNING logg %s len(wav) %d nwav %d ' % (logg, len(wavelen), nwav)
        exit()

    wavelen = numpy.array(wavelen)
    
    flux = []
    Teff = []
    for line in iterator:
        line = line.replace('Effective temperature', 'Teff')
        line = line.replace('=','')
        line = line.split()
        temperature = float(line[1])
        Teff.append(temperature)
        fluxDummy = []
        while len(fluxDummy) < nwav:
            line = iterator.next()
            line = line.split()
            for vv in line:
            
                #Because sometimes the formatting is broken
                #and the 'E' is missing from scientific notation
                if '-' in vv and 'E' not in vv:
                    vv = vv.replace('-','E-')
                if '+' in vv and 'E' not in vv:
                    vv = vv.replace('+', 'E+')
                fluxDummy.append(float(vv))
        
        if len(fluxDummy) != nwav:
            print 'WARNING logg %s Teff %e len(f) %d nwav %d' % \
                  (logg, temperature, len(fluxDummy), nwav)
    
            exit()
    
        flux.append(fluxDummy)
        
    flux = numpy.array(flux)
    wavelen = wavelen*0.1
    flux = flux/(wavelen*wavelen)
    
    fluxT = numpy.transpose(flux)
    Teff = numpy.array(Teff)
    
    if chemistry == 'He':
        prefix = 'bergeron_He_'
    else:
        prefix = 'bergeron_'
    
    for desiredT in desiredTemperatures:
        
        outName = prefix + str(desiredT[1]) + '_' + logg +'.dat_' + str(desiredT[0]) + '.gz'
        outName = os.path.join(outputDir, outName)
        
        if desiredT[0] == desiredT[1]:
            ix = numpy.where(numpy.abs(float(desiredT[0]) - Teff)<1.0)
            if len(ix) != 1:
                print 'WARNING found %d indices ' % len(ix)
                exit()
            else:
                ff = flux[ix[0][0]]
        else:
            iInterp = numpy.where(Teff<float(desiredT[0]))
            if desiredT[1] != int(Teff[iInterp[0][-1]]):
                print 'using interp ', (desiredT[1] == int(Teff[iInterp[0][-1]])), '%d %d -- %d' % \
                (desiredT[1], int(Teff[iInterp[0][-1]]), desiredT[0])
            
            ff = []
            for i in range(nwav):
                ff.append(numpy.interp(float(desiredT[0]), Teff, fluxT[i]))
            ff = numpy.array(ff)
            if len(ff) != nwav:
                print 'WARNING numpy.interp returne len %d ' % len(ff)

        normVal = numpy.interp(500., wavelen, ff)

        output = gzip.open(outName, 'wb')
        output.write('# Version 2.0 - nm, norm@500 to 1')
        output.write('# Wavelength(nm)  Flambda(ergs/cm^s/s/nm)')
        
        if len(wavelen)==0 or len(ff)==0:
            print 'WARNING writing empty file to ',outName
            exit()
        
        for (ll, f) in zip(wavelen, ff):
            output.write('%e %e\n' % (ll, f/normVal))
        output.close()

if __name__ == "__main__":

    fileFile = open('bergeron_models_needed.sav','r')
    lines = fileFile.readlines()
    fileFile.close()
    
    #will be a dict keyed on log(g)
    #values will be tuples: first entry will be actual temperature
    #second entry will be minimum temperature for interpolation
    desiredHFiles = {}
    desiredHeFiles = {}
    
    for name in lines:
        name = name.replace('_He','He')
        name = name.replace('.', '_')
        name = name.split('_')
        
        if name[0] == 'bergeron':
            if name[2] not in desiredHFiles:
                desiredHFiles[name[2]] = []
            
            desiredHFiles[name[2]].append((int(name[4]), int(name[1])))
        elif name[0] == 'bergeronHe':
            if name[2] not in desiredHeFiles:
                desiredHeFiles[name[2]] = []
            
            desiredHeFiles[name[2]].append((int(name[4]), int(name[1])))
        else:
            print 'WARNING cannot interpret: ',name
            exit()

    for logg in desiredHFiles:
        parseBergeron(logg, desiredHFiles[logg], chemistry='H')
        
    for logg in desiredHeFiles:
        parseBergeron(logg, desiredHeFiles[logg], chemistry='He')
    
 
       
