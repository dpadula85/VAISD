#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import rc, rcParams
from matplotlib import ticker, gridspec

# rcParams.update({'figure.autolayout': True})
rcParams.update({'figure.subplot.left': 0.175})
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

eV2wn = 8065.544005

# Lorentzian gamma : FWHM
# Gaussian sigma : Std Dev
def lorentzian(SpecRange, freq, Gamma):
    '''
    Function to compute a Lorentzian shape extending over SpecRange, centred
    at freq, and with Gamma width.

    Parameters
    ----------
    SpecRange: np.array (N).
        spectral range.
    freq: float.
        centre of the Lorentzian.
    Gamma: float.
        width of the Lorentzian.

    Returns
    -------
    lsfunc: np.array (N).
        line shape function extending in the desired region.
    '''
    
    lsfunc = 1 / np.pi * Gamma / (Gamma**2 + (SpecRange - freq)**2)

    return lsfunc


def damped_lorentzian(SpecRange, freq, Gamma):
    '''
    Function to compute a Damped Lorentzian shape extending over SpecRange,
    centred at freq, and with Gamma width.

    Parameters
    ----------
    SpecRange: np.array (N).
        spectral range.
    freq: float.
        centre of the Damped Lorentzian.
    Gamma: float.
        width of the Damped Lorentzian.

    Returns
    -------
    lsfunc: np.array (N).
        line shape function extending in the desired region.
    '''
    
    lsfunc = 2 / np.pi * SpecRange**2 * Gamma / ((Gamma * SpecRange)**2 + (SpecRange**2 - freq**2)**2)

    return lsfunc


def gaussian(SpecRange, freq, Sigma):
    '''
    Function to compute a Gaussian shape extending over SpecRange, centred
    at freq, and with Sigma width.

    Parameters
    ----------
    SpecRange: np.array (N).
        spectral range.
    freq: float.
        centre of the Lorentzian.
    Sigma: float.
        width of the Lorentzian.

    Returns
    -------
    lsfunc: np.array (N).
        line shape function extending in the desired region.
    '''

    expfac = -0.5 * ((SpecRange - freq) / Sigma)**2
    lsfunc = 1 / Sigma * np.sqrt(2 * np.pi) * np.exp(expfac)

    return lsfunc


def calcspecden(freqs, lambdas, **kwargs):
    '''
    Function to compute the Spectral Density of a set of oscillators with
    frequencies freqs and reorganisation energies lambdas.

    Parameters
    ----------
    freqs: np.array (N).
        oscillation frequencies (in wavenumbers).
    lambdas: np.array (N)
        reorganisation energies (in wavenumbers).

    Returns
    -------
    result: np.array (N,3).
        Spectral Density in the desired Spectral Range. The first two columns
        describe the Spectral range in wavenumbers or eV, respectively.
        The last column contains the Spectral Density in wavenumbers.
    '''

    widths = kwargs.pop("widths", None)
    ls = kwargs.pop("widths", "lor")
    xl = kwargs.pop("llim", 0)
    xu = kwargs.pop("ulim", 2000)
    step = kwargs.pop("step", 1)

    if not widths:
        widths = np.ones_like(freqs) * 5.0

    if ls == "gau":
        ls = gaussian

    elif ls == "lor":
        ls = lorentzian

    elif ls == "dlor":
        ls = damped_lorentzian

    SpecRange = np.arange(xl, xu, step)
    SpecDen = np.zeros(SpecRange.shape[0])

    for i in range(len(freqs)):
        SpecDen += np.pi * freqs[i] * lambdas[i] * ls(SpecRange, freqs[i], widths[i])

    # Convert SpecRange to eV
    SpecRange_eV = SpecRange / eV2wn

    # Prepare results
    result = np.c_[SpecRange, SpecRange_eV, SpecDen]

    return result


def plot(data, unit="wn", figname="specden.svg"):


    if unit == "eV":
        x = data[:,1]
        xlabel = 'E / eV'
        xmaj = ticker.MultipleLocator(0.05)

    elif unit == "wn":
        x = data[:,0]
        xlabel = r'$\tilde{\nu}$ / cm$^{-1}$'
        xmaj = ticker.MultipleLocator(500)

    specden = data[:,-1]

    # Set up a plot
    fig = plt.figure() 
    gs = gridspec.GridSpec(1, 1)
    gs.update(wspace=0.1, hspace=0.1)

    ax = plt.subplot(gs[0])
    ax.plot(x, specden)

    # Plot options
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r'$J(\omega)$ / cm$^{-1}$')

    ax.set_xlim(x.min(), x.max() + 1)
    ax.set_ylim(0)

    xminloc = ticker.AutoMinorLocator(5)
    yminloc = ticker.AutoMinorLocator(5)

    ax.xaxis.set_major_locator(xmaj)
    ax.xaxis.set_minor_locator(xminloc)
    ax.yaxis.set_minor_locator(yminloc)

#    ax.yaxis.set_label_coords(-0.1, 0.5)

    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.tick_params(axis='both', which='major', direction='in', pad=10, length=5)
    ax.tick_params(axis='both', which='minor', direction='in', pad=10, length=2)
    # ax.relim()
    # ax.autoscale_view(True,True,True)
    # ax.set_aspect(0.65/ax.get_data_ratio())
    rcParams.update({'font.size': 18})
    # plt.gcf().subplots_adjust(bottom=0.15)
    plt.savefig(figname, dpi=600)
    plt.close()
    # plt.show()

    return


if __name__ == '__main__':

    import textwrap
    import argparse as arg

    parser = arg.ArgumentParser(description='Convolutes Stick Spectra',
                                formatter_class=arg.ArgumentDefaultsHelpFormatter)
    
    #
    # Input files
    #
    inp = parser.add_argument_group("Input Data")
    
    inp.add_argument('-i', '--input',
                     default="results.txt", type=str, dest="InputFile",
                     help='''Input file''')
    
    #
    # Spectra Options
    #
    spec = parser.add_argument_group("Spectra Convolution Options")
    
    spec.add_argument('--ls',
                     default="lor", type=str, choices=["gau", "lor", "dlor"],
                     dest="LineShape", help='''Spectral LineShape.''')
    
    spec.add_argument('--lw',
                     default=[5.0], type=float, nargs='+', dest="LineWidth",
                     help='''Spectral LineWidth in wavenumbers (gamma for Lorentzian,
                     sigma for Gaussian LineShape.''')
    
    spec.add_argument('--unit',
                     default="wn", type=str, choices=["eV", "wn"], dest="SpecUnit",
                     help='''X axis unit for plotting Spectra.''')

    #
    # Output Options
    #
    out = parser.add_argument_group("Output Options")

    out.add_argument('-o', '--out',
                     default="Spec", type=str, dest="OutPref",
                     help='''Output Prefix for files''')

    out.add_argument('--figext',
                     default=None, type=str, choices=["svg", "png", "eps", "pdf"],
                     dest="FigExt", help='''Format for image output''')
    
    # Parse and create the Options Dictionary
    args = parser.parse_args()
    Opts = vars(args)

    if Opts['OutPref']:
        Opts['OutPref'] = Opts['OutPref'] + "."

    # Convolute Spectra
    data = np.loadtxt(Opts['InputFile'])
    freqs = data[:,1]
    lambdas = data[:,-1]

    widths = Opts['LineWidth']
    widths = Opts['LineWidth']

    # Check LineWidths array length in comparison to ExcEnergies
    try:
        assert len(widths) == len(freqs)

    except AssertionError:
        widths = Opts['LineWidth'] * len(freqs)

    finally:
        assert len(widths) == len(freqs)

    specden = calcspecden(freqs, lambdas, **Opts)

    # Convoluted Spectral Density
    titles = ("E (cm^-1)", "E (eV)", "Spectral Density (cm^-1)")
    header = ("\nCalculated Spectral Density\n"
              "Lineshape: %s\n"
              "Linewidths (cm^-1): %s \n\n" % (Opts['LineShape'], widths[0]))

    header1 = "%9s %10s %25s\n" % (titles)
    fmt = "%11.4f %10.8f %17.8e"
    np.savetxt(Opts['OutPref'] + "specden.txt", specden, fmt=fmt, header=header+header1)

    if Opts['FigExt']:

        figname = Opts['OutPref'] + "specden." + Opts['FigExt'] 
        plot(specden, unit=Opts['SpecUnit'], figname=figname)

    pass
