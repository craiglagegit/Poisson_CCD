April 2015

The SEDs in this directory correspond to asteroid types in the Bus-Demeo
taxonomy.  To generate them, download the machine-readable asteroid
reflectances available here

http://sbn.psi.edu/pds/resource/busdemeotax.html

and run the script astcolors.py on the provided file meanspectra.tab

astcolors.py will take the provided asteroid reflectance spectra and convolve
them with the solar spectrum in kurucz_sun.gz to produce the observed SED of the
given asteroid types.
