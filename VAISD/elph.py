#!/usr/bin/env python

import numpy as np

# From https://docs.scipy.org/doc/scipy/reference/constants.html
au2J = 4.35974465e-18
wn2Hz = 2.99792458e10
ang2m = 1.0e-10
au2ang = 0.5291771
AMU2Kg = 1.66053904e-27
h = 6.62607004e-34

# For the theory, refer to PCCP, 2012, 14, 13549-13563.
def adiabshift(masses, gsfreqs, gsmodes, gscoor, escoor):
    '''
    Function to compute dimensionless shifts, reorganisation energies, and
    Huang-Rhys factors between two geometries of two different electronic
    states within the Adiabatic Shift approximation.

    Parameters
    ----------
    masses: np.array (N).
        masses.
    gsfreqs: np.array (3N-6).
        ground state vibrational frequencies.
    gsmodes: np.array (3N-6,N,3).
        ground state normal modes.
    gscoor: np.array (N,3).
        ground state equilibrium geometry.
    escoor: np.array (N,3).
        excited state equilibrium geometry.

    Returns
    -------
    ks: np.array (3N-6).
        dimensionless shifts.
    lambdas: np.array (3N-6).
        reorganisation energies.
    HRs: np.array (3N-6).
        Huang-Rhys factors.
    '''

    # This is a matrix product between M^0.5 and each element of gsmodes
    mwc_gsmodes = np.einsum('j,ijk->ijk', np.sqrt(masses), gsmodes)

    # Calculate the norm of mwc_gsmodes array, keeping the i axis (axis=0)
    # (dim = 3N-6)
    norm = np.linalg.norm(mwc_gsmodes, axis=(1,2))

    # This divides each element of mwc_gsmodes by the corresponding element of norm
    mwc_gsmodes = np.einsum('i,ijk->ijk', 1 / norm, mwc_gsmodes)

    # Equilibrium Geometries in Mass-Weighted Cartesian Coords
    # This is a matrix product between M^0.5 and gscoor
    mwc_gscoor = np.einsum('j,jk->jk', np.sqrt(masses), gscoor)
    mwc_escoor = np.einsum('j,jk->jk', np.sqrt(masses), escoor)

    # Unit conversion to SI to obtain Shifts and Dimensionless Shifts
    freqsHz = 2 * np.pi * gsfreqs * wn2Hz
    mwc_gscoor = mwc_gscoor * ang2m * np.sqrt(AMU2Kg)
    mwc_escoor = mwc_escoor * ang2m * np.sqrt(AMU2Kg)
    
    # mwc_ks : MWC shift with respect to the ES equilibrium geometry (SI units m * sqrt(Kg))
    # projmwc_ks : MWC shift projected onto MWC normal modes
    # ks : dimensionless shift with respect to the ES equilibrium geometry (checked against FCClasses)
    # HRs : Huang-Rhys factors (dimensionless), total by np.sum()
    # lambdas : reorganisation energies (wn), total by np.sum (checked against FCClasses)
    # All these arrays have only the i axis (dim = 3N-6)
    mwc_ks = mwc_gscoor - mwc_escoor
    projmwc_ks = np.einsum('jk,ijk->i', mwc_ks, mwc_gsmodes)
    ks = projmwc_ks * np.sqrt(freqsHz * 2 * np.pi/h)
    HRs = 0.5 * ks**2
    lambdas = 0.5 * gsfreqs * ks**2
    
    return ks, lambdas, HRs


# For the theory, refer to PCCP, 2012, 14, 13549-13563.
def vg(masses, gsfreqs, gsmodes, esgrad):
    '''
    Function to compute dimensionless shifts, reorganisation energies, and
    Huang-Rhys factors between two geometries of two different electronic
    states within the Vertical Gradient approximation.

    Parameters
    ----------
    masses: np.array (N).
        masses.
    gsfreqs: np.array (3N-6).
        ground state vibrational frequencies.
    gsmodes: np.array (3N-6,N,3).
        ground state normal modes.
    esgrad: np.array (N,3).
        excited state gradients at the Franck-Condon point.

    Returns
    -------
    ks: np.array (3N-6).
        dimensionless shifts.
    lambdas: np.array (3N-6).
        reorganisation energies.
    HRs: np.array (3N-6).
        Huang-Rhys factors.
    '''

    # Define M^-1 1d array, with axis j (dim = N)
    invmasses = 1 / masses

    # This is a matrix product between M^0.5 and each element of gsmodes
    mwc_gsmodes = np.einsum('j,ijk->ijk', np.sqrt(masses), gsmodes)

    # Calculate the norm of mwc_gsmodes array, keeping the i axis (axis=0)
    # (dim = 3N-6)
    norm = np.linalg.norm(mwc_gsmodes, axis=(1,2))

    # This divides each element of mwc_gsmodes by the corresponding element of norm
    mwc_gsmodes = np.einsum('i,ijk->ijk', 1 / norm, mwc_gsmodes)

    # Cartesian Gradient in Mass-Weighted Cartesian Coords
    # This is a matrix product between M^-0.5 and esgrad
    mwc_esgrad = np.einsum('j,jk->jk', np.sqrt(invmasses), esgrad)
    
    # Project MWC Gradient onto MWC modes
    # This multiplies mwc_esgrad with each element of mwc_gsmodes and sums along
    # j and k axes, leavinh the i axis (dim = 3N-6)
    projgrad = np.einsum('jk,ijk->i', mwc_esgrad, mwc_gsmodes)
    # np.savetxt(sys.stdout, projgrad, fmt="%18.6e")
    
    # Unit conversion to SI to obtain Shifts and Dimensionless Shifts
    freqsHz = 2 * np.pi * gsfreqs * wn2Hz
    projgrad = projgrad * au2J / (ang2m * np.sqrt(AMU2Kg))
    
    # mwc_ks : MWC shift of the ES at the GS geometry (SI units m * sqrt(Kg))
    # ks : dimensionless shift of the ES at the GS geometry (checked against FCClasses)
    # HRs : Huang-Rhys factors (dimensionless), total by np.sum()
    # lambdas : reorganisation energies (wn), total by np.sum (checked against FCClasses)
    # All these arrays have only the i axis (dim = 3N-6)
    mwc_ks = -projgrad / freqsHz**2
    ks = mwc_ks * np.sqrt(freqsHz * 2 * np.pi/h)
    HRs = 0.5 * ks**2
    lambdas = 0.5 * gsfreqs * ks**2
    
    return ks, lambdas, HRs


if __name__ == '__main__':
    pass
