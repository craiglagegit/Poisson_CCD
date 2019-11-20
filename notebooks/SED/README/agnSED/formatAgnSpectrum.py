import gzip
import numpy

redPower = -1.5 #power law index of F_lambda (lambda < 89 nm)
bluePower = -0.5 #power law index of F_lambda (lambda > 800 nm)

redTransition = 89.0 #wavelength in nm at which to transition to red power law
blueTransition = 800.0 #wavelength in nm at which to transition to blue power law

rawFile = 'agn_spectrum_raw.txt'
dtype = numpy.dtype([('wavelen',numpy.float),
                     ('flux',numpy.float),
                     ('sigma',numpy.float)])

data = numpy.genfromtxt(rawFile,dtype=dtype)
data['wavelen'] *= 0.1

wavelen = numpy.arange(9.0, redTransition, 1.0)

normFlux = numpy.interp(redTransition, data['wavelen'], data['flux'])
norm = normFlux*numpy.power(redTransition, -1.0*redPower)
flux = norm*numpy.power(wavelen,redPower)

dexes = numpy.where(numpy.logical_and(data['wavelen'] > redTransition,
                    data['wavelen'] < blueTransition))
wavelen = numpy.append(wavelen, data['wavelen'][dexes])
flux = numpy.append(flux, data['flux'][dexes])

normFlux = numpy.interp(blueTransition, wavelen, flux)
norm = normFlux*numpy.power(blueTransition, -1.0*bluePower)

blue = numpy.arange(wavelen[-1]+1.0,30000.0,10.0)
blueFlux = norm*numpy.power(blue,bluePower)
wavelen = numpy.append(wavelen, blue)
flux = numpy.append(flux, blueFlux)

norm = numpy.interp(500.0, wavelen, flux)
flux = flux/norm

with gzip.open('agn.spec.gz','wb') as output:
    for w, f in zip(wavelen, flux):
        output.write('%e %e\n' % (w,f))

print numpy.diff(wavelen).min()
