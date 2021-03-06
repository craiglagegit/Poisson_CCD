March 2015

The models in this directory come from

http://www.stsci.edu/hst/observatory/crds/k93models.html

The specific models we use are generated by interpreting the baseline Kurucz
grid in effective temperature.

The naming convention of the spectra is

- The letter and number after the 'k' denotes the logarithmic metallicity
relative to solar. km30 means that log(Z/Z_solar) = -3.0.  kp25 means that
log(Z/Z_solar) = 2.5

- The next number is the effective temperature of the baseline Kurucz model in
the low-Teff bin used to interpolate the effective temperature, i.e.
km30_5250.fits_g00_5270.gz is a model with an effective temperature of 5270 K
which was interpolated using baseline models with Teff=5250 K and Teff=5500 K

- The number after the g denotes the log of the surface gravity

- The final number is the effective temperature of the model

The python script makeNewKuruczSEDs_1503.py (courtesy of Bryce Kalmbach; edited
by Scott Daniel) was used to do the interpolation.

kurucz_contents.sav lists the SED files we wanted to end up with and is used
as an input ot makeNewKuruczSEDs_1503.py.

Note that if you run makeNewKuruczSEDs_1503.py, a few of the spectra will come
out with the wrong name due to an ambiguity of how to treat spectra that require
no interpoloation (e.g. km30_5250.fits_g00_5250.gz)
