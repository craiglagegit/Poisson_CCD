The .dat versions are the originals from John Gizis in Mar 2011.  K. Krughoff
modified these on Mar 28 2011 as follows:
1. Converted wavelength from angstrom to nanometers
2. Interpolated and rebinned to 1nm from 300nm to 1200nm inclusive
3. In the cases of RedE1 and RedE2 converted from lambda*f(lambda) to
f(lambda)
4. Normalized spectra to unity at 500nm wavelength bin.

The interpolated spectra are in *_interp.dat

-------------------------------

11 May 2015

The spectra were extended to cover the wavelength range
9nm < lambda < 30,000nm.  In the case of PNastrom_interp.dat,
RedE1astrom_interp.dat and RedE2astrom_interp.dat, this
was a simple matter of extending the existing spectra
through interpolation (see the scripts extendPNastrom.py
and extendRedAstrom.py).

For the brown dwarf spectra, the original spectra provided
by Gizis were matched with spectra from the Bt-Settl
grid of spectra (which can be found either here

https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011/SPECTRA/

or in the mlt/ subdirectory of this library) and the 
non-equilibrium Brown Dwarf models provided by Adam Burrows here

http://www.astro.princeton.edu/~burrows/non/non.html

Where necessary, the fit spectra were padded with
very small flux values (see the script fitBrownDwarfs.py).
The following indicates the mapping between the provided
Brown Dwarf spectra and the spectral grids

BD1000_interp.dat.gz is clr_neq/T900_g50_d6f0.clr from Burrows
BD1000e_interp.dat is the same with emission lines added by hand
(see linesLookup.txt for the specific lines)
BD1500_interp.dat.gz is c100_neq/T1600_g50_d4f0.c100 from Burrows
BD1800_interp.dat.gz is lte023-5.0+0.5a+0.0.BT-Settl.spec.gz from Bt-Settl
BD2000.dat.gz is lte024-5.0+0.5a+0.0.BT-Settl.spec.gz from Bt-Settl
BD325.dat.gz is lte004.5-4.5-0.0a+0.0.BT-Settl.spec.gz from Bt-Settl
BD555.dat.gz is c100_neq/T700_g55_d2f2.c100 from Burrows
