The SEDs in this directory come from the BT-Settl CIFIST 2011 grid provided here

https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011/SPECTRA/

Below we reproduce the README.txt file that comes with those spectra.

In this directory we provide the scripts used to process and format those
files (since they are quite large).

CIFISTcurlScript.sh will download the raw SED files

binCIFIST2011spectra.py will resample the SEDs on a more reasonable wavelength grid

photometryCIFIST2011.py will integrate the resampled SEDs over the sdss bandpasses,
resulting in the lookup table used to match SEDs to objects in teh database (MLTlookup.txt)


------------- original README below --------------

Models are provided for different solar abundances via the following
subdirectories: 
GNS93: using the Grevesse, Noels & Sauval (1993) solar abundances
AGSS2009: using the Asplund et al. (2009) solar abundances
CIFIST2011: using the Caffau et al. (2011) solar abundances
CIFIST2011b: accounting also for a calibration of the mixing length based on
RHD simulations by Freytag et al. (2010).
CIFIST2011bc: accounting for the calibration of the mixing length based on
RHD simulations by Freytag et al. (2010). Additonal adjustments to the MLT
equations are taken into account.
CIFIST2011_2015: published version of the BT-Settl grid (Baraffe et al. 2015,
Allard et al. 2015 in preparation).  This grid will be the most complete
of the CIFIST2011 grids above, but currently: Teff= 1200 - 7000K, logg=2.5 - 5.5, 
[M/H]= 0.0.  To be put online when Baraffe et al. (2015) is accepted.  
Since the model atmospheres are using spherical radiative transfer the 
interior+atmosphere problem becomes iterative, and futur versions of the 
synthetic spectra grid will use the radius and the lithium abundance of 
the current evolution models.  

And different results of the model atmosphere computations are distributed in
different sub-directories:
   STRUCTURES  the standard output (fort.6) of Phoenix
   RESTARTS    part of the same information with a more concise format (fort.20)
   SPECTRA     the synthetic spectra (integrated spectra: fort.7 or intensity spectra: fort.9 )
   BIN_SPECTRA same spectra but in Hierarchical Data Format (HDF) binary format
   FITS        same spectra but in FITS format
   COLORS      corresponding  synthetic colors tables 
   ISOCHRONES  the Baraffe et al. (1997, 1998, 2003) (soon to be updated with Baraffe et al. 2015) 
               theoretical isochrones interpolated into the synthetic color tables
   OLD         those directories located in the respective ISOCHRONES and COLORS
               directories contain the results computed all from the same early 
               version of the BT-Settl model atmosphere grid.

Alpha element enhancement is taken into account:
For [M/H]= -0.0, +0.3, +0.5 no anhancement
    [M/H]= -0.5 we use an value of [alpha/H]=+0.2
    [M/H]= -1.0, -1.5, -2.0, -2.5, -3.0, -3.0, -3.5, -4.0 with [alpha/H]=+0.4
But some incomplete grids of spectra are available for other values of the alpha
enrichment. In the _hot packages the files for Teff =2600K to 70000K mention
the alpha enhancement value in the file name: ex. a+0.2 means [alpha/H]=+0.2
relative to solar values.

The spectral resolution used to compute the *.spec.7.bz2 spectra is following or better:
0.02A in the optical, and 0.05A in the infrared for the most recent files. Otherwise:
0.1A between 10 and 3600A, 0.05A between 10600 and 25000A, 0.2A between 52000 and 270000A
The simulator computes the sampling spacing according to lam1/500000 (where
lam1 is the short wavelength bound of the range specified by the user) so that
the resulting spacing is small in the UV and visual (0.01A at 5000A) and
increase towards the infrared (0.1A at 5mu) as is usually the case of the
observations.  However if the user asks for an even smaller spacing, the user
value is used in the limit of the size of the Phoenix code vectors (2500000
wavelengths points allowed currently). 

Why 2600K? This is because of the phase change due to cristalisation of
silicates in atmospheres below Teff= 2600K. Models are provided for different
treatments of this phase change. In this directory, a cloud model is used.
Associated papers can be found here:
http://phoenix.ens-lyon.fr/Grids/Papers/
Baraffe15.pdf latest publication of the consistent interior and evolution BT-Settl models.
BT-Settl.pdf is the only refereed paper (Allard et al. 2012) on the BT-Settl model atmospheres thus far.
Allard13.pdf latest news about the development of the BT-Settl models.
Allard12.pdf detailed and in depth description of the BT-Settl model atmosphere
Allard03.pdf explains in a conference proceedings the basic approach
Allard07.pdf is the publication of line broadening opacity results

The names and format of the files is explained under this link:
http://phoenix.ens-lyon.fr/Grids/FORMAT

[Scott Daniel, 2015 Juy 27: this file says that the filenames are formatted according to
      lte{Teff/10}....
upon inspection, it appears that this is supposed to say lte{Teff/100}...]

If you need models between the grid points, at a higher spectral resolution or
for different elemental abundances you have the possibility of running the
model atmosphere code PHOENIX via this job submission page:
http://phoenix.ens-lyon.fr/simulator 
by clicking the Run Phoenix button (after selecting the model grid 
opacity setup you want to use from the scrolling menu to the left).  

Please address questions and remarks to fallard@ens-lyon.fr


 
