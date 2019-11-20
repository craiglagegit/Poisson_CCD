---
author:
- Craig Lage
title: '**Poisson Solver Code**'
---

This code is a simple grid-based Poisson's equation solver intended to
simulate pixel distortion effects in thick fully-depleted CCD's. The
code builds a 3D rectilinear grid to represent a portion of the CCD,
assigns the appropriate charge densities and applied potentials, then
solves Poisson's equation using multi-grid methods. The code
self-consistently solves for the electron distribution and the potential
inside the CCD. A 360^3 grid, which is adequate for most purposes,
solves in less than one minute on a typical laptop. The code also
includes prescriptions to propagate electrons from a given point of
creation by an incoming photon down to the point of collection,
including both drift and diffusion. Most data is saved as hdf files. The
current code is configured to model both the ITL STA3800C CCD and the
E2V CCD250, but other CCDs can be modeled by editing the configuration
file. Plotting routines are available to plot the potentials, E-Fields,
pixel shapes, and electron paths. A description of the code, the
measurements which were used to validate the code, and some samples of
the output are in the file docs/Poisson_Paper_19Nov19.pdf. Below is a
basic description of how to install the code and a number of examples.

The code contains many options, and not all combinations have been
tested together. If you find a set of options that does not work as you
expect, please let me know. However, all of the example configuration
files described in the Examples Section below have been tested.

This work is supported by DOE HEP Grant DE-SC0009999.

Installing: Read the Installation Section below.

Running: The basic syntax is:

src/Poisson   config_file.cfg

More details are provided in the Examples Section.

Hopefully you find the code useful. Comments and questions are
encouraged and should be addressed to: cslage@ucdavis.edu

Dependencies:

There is one dependency that need to be installed before you can compile
the Poisson code. The former dependency on the C++ Boost library has
been removed.

1.  HDF5 libraries. There are several options for installing these:

    1.  Ubuntu: Install the hdf5 libraries using: sudo apt-get install
        hdf5-tools

    2.  Mac OSX: Assuming you are using homebrew, install using: brew
        install hdf5

    3.  Build them from source. They can be downloaded from:
        www.hdfgroup.org/HDF5/release/obtain5.html

2.  After installing the hdf5 libraries, edit the src/Makefile line
    HDF5_DIR to point to the correct location.

3.  In the src directory, type "make". This should build the Poisson
    code, and create an executable called src/Poisson. Depending on
    where you have installed the above libraries, you may need to edit
    your LD_LIBRARY_PATH environment variable so the system can find
    the appropriate files for linking.

4.  Several sample Makefiles have been included in the src directory,
    including a Makefile.nersc that works for me on NERSC Cori.

Running the python plotting routines also requires that you install h5py
so that Python can read the HDF5 files. All of the Python plot routines
have been edited to be Python3 compatible.

If you run the forward modeling code in order to generate
brighter-fatter plots as described in the bfrun example below, you will
also need to build the forward.so Python extension. Instructions for
this are in the pysrc/forward_model_varying_i directory.

The code has been cleaned up quite a bit so that all of the
configuration parameters in the various examples work together. The
previous most up-to-date version was at
<https://github.com/craiglagegit/Poisson_CCD22/tree/hole20>, but all
work should move to the master branch at
<https://github.com/craiglagegit/Poisson_CCD>.

There are a total of 14 examples included with the code. Each example is
in a separate directory in the data directory, and has a configuration
file of the form file.cfg. The parameters in the file.cfg files are
commented to explain(hopefully) the purpose of each parameter, and a
detailed listing of all configuration parameters is in Appendix A of the
docs/Poisson_Paper_19Nov19.pdf paper. Python plotting routines are
included with instructions below on how to run the plotting routines and
the expected output. The plot outputs are placed in the data/*/plots
files, so you can see the expected plots without having to run the code.
If you edit the .cfg files, it is likely that you will need to customize
the Python plotting routines as well.

-   Example 1: data/smallpixel/smallpixel.cfg

    1.  Purpose: A single pixel and surroundings. The central pixel
        contains 100,000 electrons. No electron tracking or pixel
        boundary plotting is done. The subgrids are saved so one can
        look at convergence. This is useful for getting things set up
        rapidly

    2.  Syntax: src/Poisson data/smallpixel/smallpixel.cfg

    3.  Expected run time:    < 1 minute.

    4.  Plot Syntax: python pysrc/Poisson_Small.py
        data/smallpixel/smallpixel.cfg 0

    5.  Plot Syntax: python pysrc/Poisson_Convergence.py
        data/smallpixel/smallpixel.cfg 0

-   Example 2: data/pixel0/pixel.cfg

    1.  Purpose: A 9x9 grid of pixels at low resolution (ScaleFactor=1).
        The central pixel contains 100,000 electrons. No electron
        tracking or pixel boundary plotting is done.

    2.  Syntax: src/Poisson data/pixel0/pixel.cfg

    3.  Expected run time:    ~ 1 minute.

    4.  Plot Syntax: python pysrc/Poisson_Plots.py
        data/pixel0/pixel.cfg 0

    5.  Plot Syntax: python pysrc/ChargePlots.py data/pixel0/pixel.cfg 0
        2

-   Example 3: data/pixel-itl/pixel.cfg

    1.  Purpose: A 9x9 grid of pixels at higher resolution
        (ScaleFactor=2). The central pixel contains 100,000 electrons.
        The parameters are set up for the ITL STA3800C CCD. This should
        give physically meaningful results, and is what was used in the
        published papers. After solving Poisson's equation, electron
        tracking is done to determine the pixel distortions.

    2.  Syntax: src/Poisson data/pixel-itl/pixel.cfg

    3.  Expected run time:    ~ 30 minutes.

    4.  Plot Syntax: python pysrc/Poisson_Plots.py
        data/pixel-itl/pixel.cfg 0

    5.  Plot Syntax: python pysrc/ChargePlots.py
        data/pixel-itl/pixel.cfg 0 2

    6.  Plot Syntax: python pysrc/VertexPlot.py data/pixel-itl/pixel.cfg
        0 2

    7.  Plot Syntax: python pysrc/Area_Covariance_Plot.py
        data/pixel-itl/pixel.cfg 0

-   Example 4: data/pixel-e2v/pixel.cfg

    1.  Purpose: A 9x9 grid of pixels at higher resolution
        (ScaleFactor=2). The central pixel contains 100,000 electrons.
        The parameters are set up for the E2V CCD250 CCD. This should
        give physically meaningful results, and is what was used in the
        published papers. After solving Poisson's equation, electron
        tracking is done to determine the pixel distortions.

    2.  Syntax: src/Poisson data/pixel-e2v/pixel.cfg

    3.  Expected run time:    ~ 30 minutes.

    4.  Plot Syntax: python pysrc/Poisson_Plots.py
        data/pixel-e2v/pixel.cfg 0

    5.  Plot Syntax: python pysrc/ChargePlots.py
        data/pixel-e2v/pixel.cfg 0 2

    6.  Plot Syntax: python pysrc/VertexPlot.py data/pixel-e2v/pixel.cfg
        0 2

    7.  Plot Syntax: python pysrc/Area_Covariance_Plot.py
        data/pixel-e2v/pixel.cfg 0

-   Example 5: data/pixel1/pixel.cfg

    1.  Purpose: A 9x9 grid of pixels at higher resolution
        (ScaleFactor=2). The central pixel contains 100,000 electrons.
        The parameters are set up for the ITL STA3800C CCD. This is
        included as an example of ElectronMethod=1, which iterates to
        find the value of the electron Quasi_Fermi level. After solving
        Poisson's equation, electron tracking is done to determine the
        pixel distortions.

    2.  Syntax: src/Poisson data/pixel1/pixel.cfg

    3.  Expected run time:    ~ 45 minutes.

    4.  Plot Syntax: python pysrc/Poisson_Plots.py
        data/pixel1/pixel.cfg 0

    5.  Plot Syntax: python pysrc/ChargePlots.py data/pixel1/pixel.cfg 0
        2

    6.  Plot Syntax: python pysrc/VertexPlot.py data/pixel1/pixel.cfg 0
        2

    7.  Plot Syntax: python pysrc/Area_Covariance_Plot.py
        data/pixel1/pixel.cfg 0

-   Example 6: data/smallcapf/smallcap.cfg

    1.  Purpose: A simple field oxide capacitor for testing convergence.

    2.  Syntax: src/Poisson data/smallcapf/smallcap.cfg

    3.  Expected run time:    ~ 1 minute.

    4.  Plot Syntax: python pysrc/Poisson_CapConvergence.py
        data/smallcapf/smallcap.cfg 0

-   Example 7: data/smallcapg/smallcap.cfg

    1.  Purpose: A simple gate oxide capacitor for testing convergence.

    2.  Syntax: src/Poisson data/smallcapg/smallcap.cfg

    3.  Expected run time:    ~ 1 minute.

    4.  Plot Syntax: python pysrc/Poisson_CapConvergence.py
        data/smallcapg/smallcap.cfg 0

-   Example 8: data/edgerun/edge.cfg

    1.  Purpose: A simulation of the pixel distortion at the top and
        bottom edges of the ITL STA3800C CCD.

    2.  Syntax: src/Poisson data/edgerun/edge.cfg

    3.  Expected run time:    ~ 15 minutes.

    4.  Plot Syntax: python pysrc/Poisson_Edge.py data/edgerun/edge.cfg
        0

    5.  Plot Syntax: python pysrc/Pixel_Shift.py data/edgerun/edge.cfg
        0

-   Example 9: data/iorun/io.cfg

    1.  Purpose: A simulation of the output transistor region of the ITL
        STA3800C CCD, including the last few serial stages.

    2.  Syntax: src/Poisson data/iorun/io.cfg

    3.  Expected run time:    ~ 5 minutes.

    4.  Plot Syntax: python pysrc/Poisson_IO.py data/iorun/io.cfg 0

    5.  Plot Syntax: python pysrc/ChargePlots_IO.py data/iorun/io.cfg 0

-   Example 10: data/treering/treering.cfg

    1.  Purpose: A simulation of a 9x9 pixel region with a sinusoidal
        treering dopant variation introduced, so one can see the
        resulting pixel distortions.

    2.  Syntax: src/Poisson data/treering/treering.cfg

    3.  Expected run time:    ~ 20 minutes.

    4.  Plot Syntax: python pysrc/VertexPlot.py
        data/treering/treering.cfg 0 4

-   Example 11: data/fe55/pixel.cfg

    1.  Purpose: A simulation of a 5x5 pixel region with a single Fe55
        event above the center pixel. To determine the distribution of
        pixel counts, one needs to run many of these simulations with a
        random center point.

    2.  Syntax: src/Poisson data/fe55/pixel.cfg

    3.  Expected run time:    ~ 5 minutes.

    4.  Plot Syntax: python pysrc/Poisson_Fe55.py data/fe55/pixel.cfg 0
        40

-   Example 12: data/transrun/trans.cfg

    1.  Purpose: A simulation of the output transistor of the ITL
        STA3800C CCD. In order to plot the I-V characteristic, this
        needs to be re-run with different gate voltage values.

    2.  Syntax: src/Poisson data/transrun/trans.cfg

    3.  Expected run time:    ~ 10 minutes.

    4.  Plot Syntax: python pysrc/Poisson_Trans.py
        data/transrun/trans.cfg 0

    5.  Plot Syntax: python pysrc/ChargePlots_Trans.py
        data/transrun/trans.cfg 0

    6.  Plot Syntax: python pysrc/MOSFET_Calculated_IV.py
        data/transrun/trans.cfg 0

-   Example 13: data/satrun/sat.cfg

    1.  Purpose: A simulation of the saturation level, as determined by
        the barrier lowering. Again, this needs to be run multiple times
        with varying parallel voltages.

    2.  Syntax: src/Poisson data/satrun/sat.cfg

    3.  Expected run time:    ~ 20 minutes.

    4.  Plot Syntax: python pysrc/Poisson_Sat.py data/satrun/sat.cfg 0

    5.  Plot Syntax: python pysrc/ChargePlots.py data/satrun/sat.cfg 0 9

    6.  Plot Syntax: python pysrc/Barrier.py data/satrun/sat.cfg 0

    7.  Plot Syntax: python pysrc/Plot_SatLevel.py data/satrun/sat.cfg
        0

-   Example 14: data/bfrun/bf.cfg

    1.  Purpose: This builds up a Gaussian spot in the center of a 9x9
        pixel region by sequentially solving Poisson's equation, adding
        electrons, and repeating. As written, this is done 80 times,
        building up the central pixel charge to about 80,000 electrons.

    2.  Syntax: src/Poisson data/bfrun/bf.cfg

    3.  Expected run time:   ~ 8 hours.

    4.  Plot Syntax: python Poisson_Plots.py data/bfrun/bf.cfg 80

    5.  Plot Syntax: python ChargePlots.py data/bfrun/bf.cfg 80 9

    6.  Plot Syntax: python Plot_BF_Spots.py data/bfrun/bf.cfg 0
        NumSpots 80 (Assumes this was run NumSpots times with random
        center locations)
