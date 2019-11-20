April 2015

This spectrum is taken from table 1 of Vanden Berk et al. 2001 (AJ 122:549)

It is extended as follows:

## composite quasar rest-frame spectrum
## Fnu, for wavelength in the range 300-12000 Ang
## Based on the SDSS composite QSO spectrum from
## Van den Berk et al., which is extended from
## its 890-8500 Ang wavelength range to
## 300-12,000 Ang range which ensures complete
## coverage in LSST wavelength range for 0 < z < 9.
## The range is extended by assuming power laws:
##  for lambda<890 Ang: F_lambda propto 1/lambda^1.5
## for lambda>8000 Ang: F_lambda propto 1/lambda^0.5
## The wavelength grid is equidistant with a step
## of 5 Ang (similar to intrinsic resolution),
## and Fnu is normalized to 1 at 5000 Ang.
##                       (Zeljko Ivezic, Oct 2009)

As of this writing, the power laws have been used to further extend the
spectrum from 9 nm to 30 microns to allow for coverage of other, non-LSST
bandpasses.

The file raw_agn_spectrum.txt in this docs directory is the raw table
as downloaded from the electronic journal.

The script formatAgnSpectrumpy converts raw_agn_spectrum.txt to the provided
spectrum (note that you will need to remove all of the header information
from raw_agn_spectrum.txt)
