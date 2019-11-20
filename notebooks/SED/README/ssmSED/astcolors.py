import os
import numpy
import pylab
import gzip
from lsst.sims.photUtils.Bandpass import Bandpass
from lsst.sims.photUtils.Sed import Sed


filterlist= ('u', 'g', 'r', 'i', 'z', 'y')
colornames = ('ug', 'gr', 'ri', 'iz', 'zy')

def read_asteroids_reflectance(dataDir='.'):
    # Read the sun's spectrum.
    sun = Sed()
    sun.readSED_flambda('kurucz_sun')
    # Read the asteroid reflectance spectra.
    allfiles = os.listdir(dataDir)
    
    asteroidDtype = numpy.dtype([
                                ('wavelength', numpy.float),
                                ('A', numpy.float),
                                ('A_sig', numpy.float),
                                ('B', numpy.float),
                                ('B_sig', numpy.float),
                                ('C', numpy.float),
                                ('C_sig', numpy.float),
                                ('Cb', numpy.float),
                                ('Cb_sig', numpy.float),
                                ('Cg', numpy.float),
                                ('Cg_sig', numpy.float),
                                ('Cgh', numpy.float),
                                ('Cgh_sig', numpy.float),
                                ('Ch', numpy.float),
                                ('Ch_sig', numpy.float),
                                ('D', numpy.float),
                                ('D_sig', numpy.float),
                                ('K', numpy.float),
                                ('K_sig', numpy.float),
                                ('L', numpy.float),
                                ('L_sig', numpy.float),
                                ('O', numpy.float),
                                ('O_sig', numpy.float),
                                ('Q', numpy.float),
                                ('Q_sig', numpy.float),
                                ('R', numpy.float),
                                ('R_sig', numpy.float),
                                ('S', numpy.float),
                                ('S_sig', numpy.float),
                                ('Sa', numpy.float),
                                ('Sa_sig', numpy.float),
                                ('Sq', numpy.float),
                                ('Sq_sig', numpy.float),
                                ('Sr', numpy.float),
                                ('Sr_sig', numpy.float),
                                ('Sv', numpy.float),
                                ('Sv_sig', numpy.float),
                                ('T', numpy.float),
                                ('T_sig', numpy.float),
                                ('V', numpy.float),
                                ('V_sig', numpy.float),
                                ('X', numpy.float),
                                ('X_sig', numpy.float),
                                ('Xc', numpy.float),
                                ('Xc_sig', numpy.float),
                                ('Xe', numpy.float),
                                ('Xe_sig', numpy.float),
                                ('Xk', numpy.float),
                                ('Xk_sig', numpy.float),
                                ])

    data = numpy.loadtxt(os.path.join(dataDir,'meanspectra.tab'),
                            dtype = asteroidDtype)

    data['wavelength'] *= 1000.0 #because spectra are in microns
 
    wavelen_step = min(numpy.diff(data['wavelength']).min(), numpy.diff(sun.wavelen).min())
    wavelen = numpy.arange(sun.wavelen[0], data['wavelength'][-1], wavelen_step)

    ast_reflect = {}
    for a in data.dtype.names:
        if a == 'wavelength' or a[-3:] == 'sig':
            continue

        # Read the basic reflectance data
        ast_reflect[a] = Sed(wavelen=data['wavelength'], flambda=data[a])

        # And now add an extrapolation to the blue end.
        # Binzel cuts off at 450nm.
        condition = ((ast_reflect[a].wavelen >= 450) & (ast_reflect[a].wavelen < 700))
        x = ast_reflect[a].wavelen[condition]
        y = ast_reflect[a].flambda[condition]
        p = numpy.polyfit(x, y, deg=2)
        condition = (wavelen < 450)
        flambda = numpy.zeros(len(wavelen), 'float')
        interpolated = numpy.polyval(p, wavelen[condition])
        flambda[condition] = [ii if ii > 0.0 else 0.0 for ii in interpolated]
        condition = (wavelen >= 450)
        flambda[condition] = numpy.interp(wavelen[condition], ast_reflect[a].wavelen, ast_reflect[a].flambda)
        ast_reflect[a] = Sed(wavelen, flambda)

    ast = {}
    for a in ast_reflect:
        ast[a] = sun.multiplySED(ast_reflect[a], wavelen_step=wavelen_step)

    for a in ast:
        name = a + '.dat'
        normalizedSed = Sed(wavelen=ast[a].wavelen, flambda=ast[a].flambda)
        norm = numpy.interp(500.0, normalizedSed.wavelen, normalizedSed.flambda)
        normalizedSed.multiplyFluxNorm(1.0/norm)
        normalizedSed.writeSED(name, print_header='21 April 2015; normalized to flambda=1 at 500nm')

    return ast_reflect, sun, ast


def plot_seds(tast):
    pylab.figure()
    for a in tast.keys():
        pylab.plot(tast[a].wavelen, tast[a].flambda)
    pylab.xlim(300, 1100)
    pylab.xlabel('Wavelength (nm)')
    pylab.ylabel('F_lambda')
    return

def calc_colors(ast):
    fdir = os.getenv('LSST_THROUGHPUTS_DEFAULT')
    lsst = {}
    for f in filterlist:
        lsst[f] = Bandpass()
        lsst[f].readThroughput(os.path.join(fdir, 'total_'+f+'.dat'))
    mags = {}
    for f in filterlist:
        mags[f] = {}
        for a in ast.keys():
            mags[f][a] = ast[a].calcMag(lsst[f])
    colormags = {}
    for i in range(0, len(filterlist)-1):
        c = colornames[i]
        c1 = c[0]
        c2 = c[1]
        colormags[c] = {}
        for a in ast.keys():
            colormags[c][a] = mags[c1][a] - mags[c2][a]
    writestring = '#Name '
    for c in colornames:
        writestring = writestring + ' %s ' %(c)
    print writestring
    for a in ast.keys():
        writestring= '%s ' %(a.strip('.txt'))
        for c in colornames:
            writestring =  writestring+ '%.3f ' %(colormags[c][a])
        print writestring
    return mags, colormags

def plot_colors(colormags, add_KBO = True):
    plotpoints = {}
    for a in colormags[colornames[0]].keys():
        plotpoints[a] = numpy.random.random(3)
    kbocolor = {'ug' : 1.52 , 'gr':.815 , 'ri':.368 , 'iz':.196, 'zy':0.2  }
    pylab.figure()
    pylab.prism()
    pylab.subplots_adjust(top=0.92, bottom=0.08, wspace=0.3, hspace=0.3)
    subplots = 4
    for i in range(0, subplots):
        pylab.subplot(2,2,i+1)
        pcolor1 = colornames[i]
        pcolor2 = colornames[i+1]
        for a in colormags[pcolor1].keys():
            pylab.plot(colormags[pcolor2][a], colormags[pcolor1][a], color=plotpoints[a],
                       marker='o', linestyle='')
            namelabel = a.strip('.txt')
            pylab.annotate(namelabel, (colormags[pcolor2][a]+0.005, colormags[pcolor1][a]+0.005),
                           color=plotpoints[a])
            if add_KBO:
                pylab.plot(kbocolor[pcolor2], kbocolor[pcolor1], marker='*', linestyle='', color='k')
        pylab.xlabel(pcolor2)
        pylab.ylabel(pcolor1)
        xmin, xmax = pylab.xlim()
        step = (xmax - xmin) / 5.0
        locs = numpy.arange(xmin, xmax+step, step)
        pylab.xticks(locs)
        ymin, ymax = pylab.ylim()
        step = (ymax - ymin) / 5.0
        locs = numpy.arange(ymin, ymax+step, step)
        pylab.yticks(locs)
        pylab.grid()
        """
    pylab.subplot(2,2,subplots)
    pcolor1 = 'gr'
    pcolor2 = 'ri'
    for a in colormags[pcolor1].keys():
        pylab.plot(colormags[pcolor2][a], colormags[pcolor1][a], marker='', linestyle='')
        namelabel = a.strip('.txt')
        pylab.annotate(namelabel, (colormags[pcolor2][a], colormags[pcolor1][a]), color=plotpoints[a])
        pylab.tick_params(labelbottom='off', labelleft='off')
        """
    pylab.suptitle('Binzel Asteroids - Colors in LSST')
    return

def plot_colors2(colormags, add_KBO = True):
    plotpoints = {}
    for a in colormags[colornames[0]].keys():
        if a[0] == 'C':
            plotpoints[a] = 'r'
        elif a[0] == 'S':
            plotpoints[a] = 'b'
        elif a[0] == 'X':
            plotpoints[a] = 'g'
        else:
            plotpoints[a] = 'y'
    kbocolor = {'ug' : 1.52 , 'gr':.815 , 'ri':.368 , 'iz':.196, 'zy':0.2  }
    for i in range(1, 3):
        pylab.figure()
        pylab.prism()
        pcolor1 = colornames[i]
        pcolor2 = colornames[i+1]
        for a in colormags[pcolor1].keys():
            pylab.plot(colormags[pcolor2][a], colormags[pcolor1][a], color=plotpoints[a],
                       marker='o', linestyle='')
            namelabel = a.strip('.txt')
            pylab.annotate(namelabel, (colormags[pcolor2][a]+0.005, colormags[pcolor1][a]+0.005),
                           color=plotpoints[a])
            if add_KBO:
                pylab.plot(kbocolor[pcolor2], kbocolor[pcolor1], marker='*', linestyle='', color='k')
        pylab.xlabel(pcolor2)
        pylab.ylabel(pcolor1)
        xmin, xmax = pylab.xlim()
        #step = (xmax - xmin) / 5.0
        #locs = numpy.arange(xmin, xmax+step, step)
        #pylab.xticks(locs)
        ymin, ymax = pylab.ylim()
        #step = (ymax - ymin) / 5.0
        #locs = numpy.arange(ymin, ymax+step, step)
        #pylab.yticks(locs)
        pylab.grid()
        pylab.title('Binzel Taxonomy - Colors in LSST')
    return



def plot_families(mags):
    pcolor1 = 'iz'
    izcolor = {}
    for a in mags['i'].keys():
        izcolor[a] = mags['i'][a] - mags['z'][a]
    pcolor2 = 'a*'
    acolor = {}
    for a in mags['i'].keys():
        acolor[a] = 0.89*(mags['g'][a]-mags['r'][a]) + 0.45*(mags['r'][a] - mags['i'][a]) - 0.57
    pylab.figure()
    for a in mags['i'].keys():
        pylab.plot(acolor[a], izcolor[a], marker='o', linestyle='')
        namelabel = a.strip('.txt')
        pylab.annotate(namelabel, (acolor[a]+0.005, izcolor[a]+0.005))
    pylab.xlabel('a*')
    pylab.ylabel('i-z')
    pylab.title('Binzel colors compared to SDSS')
    return


if __name__ == '__main__':
    ast_reflect, sun, ast = read_asteroids_reflectance()
    mags, colormags = calc_colors(ast)
    plot_seds(ast)
    plot_colors(colormags, add_KBO=False)
    plot_colors2(colormags, add_KBO=False)
    plot_families(mags)

    pylab.show()
