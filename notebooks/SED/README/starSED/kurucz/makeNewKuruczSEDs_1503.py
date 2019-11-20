import pyfits
import numpy
import math
import gzip
import os

def write_spectra(metal, desiredFiles=[], directory='.'):

    ct = 0
    desiredGravity = {}
    for name in desiredFiles:
        name = name.replace('.','_')
	name = name.split('_')
	if name[0] == metal:
	   gg = name[3]
	   tt = name[4]
	   
	   if gg not in desiredGravity:
	       desiredGravity[gg] = []
	   desiredGravity[gg].append(tt)

    if len(desiredGravity)==0:
        return
    
    print desiredGravity
    
    gravity = ['g00', 'g05', 'g10', 'g15', 'g20', 'g25', 'g30', 'g35', 'g40', 'g45', 'g50']

    temps = [3500, 10000, 13000, 35000, 50001]

    if metal != 'kp10' and metal != 'kp05':
        temps = [3500, 10000, 13000, 35000, 50001]
    #Need to treat kp10 differently
    elif metal == 'kp10':
        temps = [3500, 10000, 13000, 35000, 40001]
    elif metal == 'kp05':
        temps = [3500, 10000, 13000, 35000, 45001]

    tempSpacing = [250, 500, 1000, 2500]
    tempInterpSpacing = [20, 50, 50, 250]
    tempNum = numpy.arange(len(temps)-1)
    temperature = ([])
    interpolatedTemps = ([])

    for temp in tempNum:
        temperature = numpy.append(temperature, numpy.arange(temps[temp], temps[temp+1], tempSpacing[temp]))
    for temp in range(0, len(temperature)-1):
        if temperature[temp] < temps[1]:
            interpolatedTemps = numpy.append(interpolatedTemps, numpy.arange(temperature[temp], temperature[temp+1], tempInterpSpacing[0]))
        elif temperature[temp] < temps[2]:
            interpolatedTemps = numpy.append(interpolatedTemps, numpy.arange(temperature[temp], temperature[temp+1], tempInterpSpacing[1]))
        elif temperature[temp] < temps[3]:
            interpolatedTemps = numpy.append(interpolatedTemps, numpy.arange(temperature[temp], temperature[temp+1], tempInterpSpacing[2]))
        else:
            interpolatedTemps = numpy.append(interpolatedTemps, numpy.arange(temperature[temp], temperature[temp+1], tempInterpSpacing[3]))
    interpolatedTemps = numpy.append(interpolatedTemps, 50000)
    
    for grav in desiredGravity:
        flux = numpy.zeros((len(temperature), 1221))

        newFlux = numpy.zeros((1221, len(interpolatedTemps)))
        newFluxLog = numpy.zeros((1221, len(interpolatedTemps)))
        waveNumber = numpy.arange(1221)
        tempNumber = numpy.arange(len(temperature))
        x=0
        for temp in temperature:

            modelFile = pyfits.open('/astro/net/pogo3/jbkalmbach/k93models/' + metal + '/' + metal + '_' + str(int(temp)) + '.fits')

            modelData = modelFile[1].data
            fl = modelData.field(grav)
            flux[x] = fl
            x+=1
	
	#need to transpose flux so that we can interpolate each wavelength
	#bin over temperature
        fluxW = numpy.transpose(flux)
        fluxWLog = numpy.zeros((1221, len(temperature)))

        wavelength = modelData.field('wavelength')
        numLambda = numpy.arange(len(wavelength))
        for y in range(0, len(fluxW)):
            for z in range(0, len(fluxW[y])):
                fluxW[y,z] = float(fluxW[y,z])
                if fluxW[y,z] != 0.0:
                    fluxWLog[y,z] = math.log10(fluxW[y,z])
        interpolatedTempsLog = numpy.zeros(len(interpolatedTemps))
        for y in range(0, len(interpolatedTemps)):
            interpolatedTemps[y] = float(interpolatedTemps[y])
            interpolatedTempsLog[y] = math.log10(interpolatedTemps[y])
        temperatureLog = numpy.zeros(len(temperature))
        for y in range(0, len(temperature)):
            temperature[y] = float(temperature[y])
            temperatureLog[y] = math.log10(temperature[y])
        
	#now, for each wavelength bin, you are interpolating every spectrum over temperature
	#(y data is an array of fluxes, one for each model; x data is an array of temperatures)
        for num in numLambda:
            newFluxLog[num] = numpy.interp(interpolatedTempsLog, temperatureLog, fluxWLog[num])
        
        for y in range(0, len(newFluxLog)):
            for z in range(0, len(newFluxLog[y])):
                if newFluxLog[y,z] != 0.0:
                    newFlux[y,z] = math.pow(10.0, newFluxLog[y,z])

        newFluxFinal = numpy.transpose(newFlux)

        wavelength = wavelength/10.0

        for t in range(0, len(interpolatedTemps)):
            tempValue = interpolatedTemps[t]
	    if str(int(tempValue)) in desiredGravity[grav]:
	    
	        #old resampling code
                #wavelen_grid = numpy.arange(9.0, 30000.1, 0.1, dtype = 'float')
                #flux_grid = numpy.interp(wavelen_grid, wavelength, newFluxFinal[t])
                
		#just print the raw Kurucz grids
		wavelen_grid = wavelength
		flux_grid=newFluxFinal[t]
		
		#print wavelen_grid, flux_grid
                #wavelen_grid value is 500.0 nm at the 2000th entry, normalize flux_grid at 500 nm
                normValue = numpy.interp(500.0, wavelen_grid, flux_grid)
		
		if numpy.sum(flux_grid > 0):
                    for w in range(0, len(wavelen_grid)):
                        if normValue != 0.0:
                            flux_grid[w] = flux_grid[w]/normValue
                    nameTemp = temperature[numpy.where(temperature <= tempValue)]

                    #f = gzip.open('/astro/net/pogo3/jbkalmbach/newKurucz/' + str(metal) + '_' + str(int(nameTemp[-1])) + '.fits_' + str(grav) + '_' + str(int(tempValue)) + '.gz', 'wb')
                    outputName = str(metal) + '_' + str(int(nameTemp[-1])) + '.fits_' + str(grav) + '_' + str(int(tempValue)) + '.gz'
	            print outputName, outputName in desiredFiles
		    print wavelen_grid
		    print flux_grid
		    f = gzip.open(os.path.join(directory,outputName), 'wb')
                
		    f.write('# Version 3.0 - nm, norm@500 to 1\n')
                    f.write('#Wavelength(nm) Flambda(ergs/cm^s/s/nm)\n')
                    for line in range(0, len(wavelen_grid)):
                        stringLambda = str(wavelen_grid[line])
                        stringFlux = str(flux_grid[line])
                        f.write(stringLambda + ' ' + stringFlux + '\n')
                    f.close()
                    print 'wrote: ' + str(metal) + '_' + str(int(nameTemp[-1])) + '.fits_' + str(grav) + '_' + str(int(tempValue)) + '.gz'

if __name__ == "__main__":

    #read in the list of files we want to generate
    contents = open('kurucz_contents.sav','r')
    lines = contents.readlines()
    contents.close()
    desiredFiles = []
    for name in lines:
        desiredFiles.append(name.strip())

    metallicities = ['km01', 'km02', 'km03', 'km05', 'km10', 'km15', 'km20', 'km25', 'km30', 'km35', 'km40', 'km45', 'km50', 'kp00', 'kp01', 'kp02', 'kp03', 'kp05', 'kp10']

    for metal in metallicities:
        write_spectra(metal, desiredFiles=desiredFiles, directory='extendedKurucz1503/')
