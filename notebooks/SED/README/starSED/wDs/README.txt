March 2015

The models in this directory were provided by Pierre Bergeron via private
communication.  They were parsed using the script extendBergeron.py. 
bergeron_models_needed.sav is a list of the models we needed to provide for our
SED library.

see here for citation

http://www.astro.umontreal.ca/~bergeron/CoolingModels/

The models in this directory are named something like

bergeron_6000_65.dat_6100

The first number represents the lower of two Teffs used to interpolate the model
spectrum (i.e. in the example above, we have a model with Teff=6100 arrived at
by interpolating between a model with Teff=6000 and a model with Teff=6250). 
The midell number is the log of the surface gravity (the example above has
log(g)=6.5).  The final number is the actual Teff of the model.

Helium white dwarf models are named

bergeron_He_####_##.dat_####
