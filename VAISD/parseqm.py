#!/usr/bin/env python

from __future__ import division
import sys
import numpy as np
from rdkit.Chem import GetPeriodicTable
pt = GetPeriodicTable()


def skiplines(openfile, nlines=0):
    '''
    Function to skip nlines + 1 lines in openfile. In other words, if nlines=0 it will
    go to the next line.

    Parameters
    ----------
    openfile: object.
        File object to process.
    nlines: int.
        Number of lines to skip.

    Returns
    -------
    line: string.
        Line after skipping nlines + 1 lines.
    '''

    for i in range(nlines + 1):
        line = next(openfile)

    return line


def parsefreqs_QChem(filename):
    '''Parses Q-Chem frequencies logfile for geometric and vibrational properties.'''

    with open(filename) as f:

        for line in f:

            #
            # Geometry
            #
            if "Standard Nuclear Orientation" in line:

                # This next line guarantees we only retrieve the last structure
                # in case a double job opt+freq is contained in the logfile
                structure = []
                line = skiplines(f, 2)
                data = line.split()

                while len(data) == 5:

                    atom = data[1]
                    atom_x = float(data[2])
                    atom_y = float(data[3])
                    atom_z = float(data[4])
                    structure.append([atom, atom_x, atom_y, atom_z])

                    data = next(f).split()

            #
            # Vibrational Properties
            #
            if "VIBRATIONAL ANALYSIS" in line:

                freqs = []
                forcecns = []
                redmasses = []
                modes = []
                line = skiplines(f, 10)

                while line:

                    #
                    # Frequencies
                    #
                    if "Frequency" in line:
                        tmpfreqs = list(map(float, line.split()[1:]))
                        freqs.extend(tmpfreqs)

                    #
                    # Force Constants
                    #
                    if "Force Cnst" in line:
                        tmpcns = list(map(float, line.split()[2:]))
                        forcecns.extend(tmpcns)

                    #
                    # Reduced Masses
                    #
                    if "Red. Mass" in line:
                        tmpredms = list(map(float, line.split()[2:]))
                        redmasses.extend(tmpredms)

                    #
                    # Cartesian Displacements
                    #
                    if "Raman Active" in line:
                        line = skiplines(f, 1)
                        disps = []

                        while "TransDip" not in line:
                            data = list(map(float, line.split()[1:]))
                            N = len(data) // 3

                            if not disps:
                                for n in range(N):
                                    disps.append([])

                            for n in range(N):
                                disps[n].append(data[3*n:(3*n)+3])

                            line = skiplines(f)

                        modes.extend(disps)

                    try:
                        line = skiplines(f)
                    
                    except StopIteration:
                        break


        Z_atoms = [ pt.GetAtomicNumber(x[0]) for x in structure ]
        masses = np.array([ pt.GetAtomicWeight(x[0]) for x in structure ])
        atoms = [ x[0] for x in structure ]
        coords = np.array([ x[1:] for x in structure ])
        freqs = np.array(freqs)
        forcecns = np.array(forcecns)
        redmasses = np.array(redmasses)
        modes = np.array(modes)

        return Z_atoms, masses, atoms, coords, freqs, forcecns, redmasses, modes


def parseopt_QChem(filename):
    '''Parses Q-Chem optimisations logfile for geometric properties.'''

    with open(filename) as f:

        for line in f:

            #
            # Geometry
            #
            if "Standard Nuclear Orientation" in line:

                # This next line guarantees we only retrieve the last structure
                # in case a double job opt+freq is contained in the logfile
                structure = []
                line = skiplines(f, 2)
                data = line.split()

                while len(data) == 5:

                    atom = data[1]
                    atom_x = float(data[2])
                    atom_y = float(data[3])
                    atom_z = float(data[4])
                    structure.append([atom, atom_x, atom_y, atom_z])

                    data = next(f).split()

        Z_atoms = [ pt.GetAtomicNumber(x[0]) for x in structure ]
        masses = np.array([ pt.GetAtomicWeight(x[0]) for x in structure ])
        atoms = [ x[0] for x in structure ]
        coords = np.array([ x[1:] for x in structure ])

        return Z_atoms, masses, atoms, coords


def parsefreqs_G09(filename):
    '''Parses G09 frequencies logfile for geometric and vibrational properties.'''

    with open(filename) as f:

        structures = []
        freqs = []
        forcecns = []
        redmasses = []
        masses = []
        modes = []
        HPmodes = False

        struct_done = False

        for line in f:

            #
            # Atomic Masses
            #
            if "AtmWgt" in line and not struct_done:
                data = list(map(float, line.split()[1:]))
                masses.extend(data)

            #
            # Non reoriented geometry
            #
            if "Input orientation" in line:
                structure = []
                line = skiplines(f, 4)
                data = line.split()

                while len(data) == 6:

                    Z_atom = int(data[1])
                    atom_x = float(data[3])
                    atom_y = float(data[4])
                    atom_z = float(data[5])
                    structure.append([Z_atom, atom_x, atom_y, atom_z])

                    data = next(f).split()

                structures.append(structure)
                struct_done = True

            #
            # Reoriented geometry, overwrite previous
            #
            if "Standard orientation" in line:
                structure = []
                line = skiplines(f, 4)
                data = line.split()

                while len(data) == 6:

                    Z_atom = int(data[1])
                    atom_x = float(data[3])
                    atom_y = float(data[4])
                    atom_z = float(data[5])
                    structure.append([Z_atom, atom_x, atom_y, atom_z])

                    data = next(f).split()

                structures.append(structure)
                struct_done = True

            #
            # Vibrational Properties
            #

            #
            # HPmodes
            #
            if "Frequencies ---" in line:
                HPmodes = True
                NAtoms = len(structure)
                tmpfreqs = list(map(float, line.split()[2:]))
                freqs.extend(tmpfreqs)

            if "Reduced masses ---" in line:
                tmpredms = list(map(float, line.split()[3:]))
                redmasses.extend(tmpredms)

            if "Force constants ---" in line:
                tmpcns = list(map(float, line.split()[3:]))
                forcecns.extend(tmpcns)

            if "Coord Atom Element" in line:
                disps = []

                for n in range(3 * NAtoms):
                    line = skiplines(f)
                    data = line.split()
                    atomindex = int(data[1]) - 1
                    numbers = list(map(float, data[3:]))
                    numbermodes = len(numbers)

                    if not disps:
                        for mode in range(numbermodes):
                            # For each mode, make list of list [atom][coord_index]
                            disps.append([[] for x in range(0, NAtoms)])

                    for mode in range(numbermodes):
                        disps[mode][atomindex].append(numbers[mode])

                modes.extend(disps)

            #
            # No HPmodes
            #
            if not HPmodes:

                if "Frequencies --" in line:
                    NAtoms = len(structure)
                    tmpfreqs = list(map(float, line.split()[2:]))
                    freqs.extend(tmpfreqs)

                if "Red. masses --" in line:
                    tmpredms = list(map(float, line.split()[3:]))
                    redmasses.extend(tmpredms)

                if "Frc consts --" in line:
                    tmpcns = list(map(float, line.split()[3:]))
                    forcecns.extend(tmpcns)

                if "X      Y      Z" in line:
                    line = skiplines(f)
                    disps = []

                    i = 1
                    while i <= NAtoms:
                        data = list(map(float, line.split()[2:]))
                        N = len(data) // 3

                        if not disps:
                            for n in range(N):
                                disps.append([])

                        for n in range(N):
                            disps[n].append(data[3*n:(3*n)+3])

                        line = skiplines(f)
                        i += 1

                    modes.extend(disps)

        if not masses:
            masses = np.array([ pt.GetAtomicWeight(x[0]) for x in structure ])

        Z_atoms = [ int(x[0]) for x in structure ]
        masses = np.array(masses)
        atoms = [ pt.GetElementSymbol(x[0]) for x in structure ]
        coords = np.array([ x[1:] for x in structures[-1] ])
        freqs = np.array(freqs)
        forcecns = np.array(forcecns)
        redmasses = np.array(redmasses)
        modes = np.array(modes)

        return Z_atoms, masses, atoms, coords, freqs, forcecns, redmasses, modes


def parseopt_G09(filename):
    '''Parses G09 optimisation logfile for geometric properties.'''

    with open(filename) as f:

        structures = []
        masses = []
        struct_done = False

        for line in f:

            #
            # Atomic Masses
            #
            if "AtmWgt" in line and not struct_done:
                data = list(map(float, line.split()[1:]))
                masses.extend(data)

            #
            # Non reoriented geometry
            #
            if "Input orientation" in line: # and not struct_done:
                structure = []
                line = skiplines(f, 4)
                data = line.split()

                while len(data) == 6:

                    Z_atom = int(data[1])
                    atom_x = float(data[3])
                    atom_y = float(data[4])
                    atom_z = float(data[5])
                    structure.append([Z_atom, atom_x, atom_y, atom_z])

                    data = next(f).split()

                structures.append(structure)
                struct_done = True

            #
            # Reoriented geometry, overwrite previous
            #
            if "Standard orientation" in line:
                structure = []
                line = skiplines(f, 4)
                data = line.split()

                while len(data) == 6:

                    Z_atom = int(data[1])
                    atom_x = float(data[3])
                    atom_y = float(data[4])
                    atom_z = float(data[5])
                    structure.append([Z_atom, atom_x, atom_y, atom_z])

                    data = next(f).split()

                structures.append(structure)
                struct_done = True


        Z_atoms = [ int(x[0]) for x in structure ]
        masses = np.array(masses)
        atoms = [ pt.GetElementSymbol(x[0]) for x in structure ]
        coords = np.array([ x[1:] for x in structures[-1] ])

        return Z_atoms, masses, atoms, coords


def parseforces_QChem(filename):
    '''Parses QChem forces logfile to retrieve only forces.'''

    with open(filename) as f:

        for line in f:

            #
            # Geometry
            #
            if "Standard Nuclear Orientation" in line:

                # This next line guarantees we only retrieve the last structure
                # in case a double job opt+freq is contained in the logfile
                structure = []
                line = skiplines(f, 2)
                data = line.split()

                while len(data) == 5:

                    atom = data[1]
                    atom_x = float(data[2])
                    atom_y = float(data[3])
                    atom_z = float(data[4])
                    structure.append([atom, atom_x, atom_y, atom_z])

                    data = next(f).split()

                NAtoms = len(structure)

            #
            # Forces
            #
            if "Gradient of the state energy" in line:
                forces = np.array([])
                tmpforces = []
                line = skiplines(f)

                while forces.size < 3 * NAtoms:
                    data = list(map(float, line.split()))
                    N = len(data)

                    if np.array_equal(np.diff(data), np.ones(len(data) - 1)) or np.diff(data) is None:
                        line = skiplines(f)
                        continue

                    cooridx = int(line.split()[0]) - 1

                    tmpforces.append(data[1:])

                    if cooridx == 2:
                        if forces.size == 0:
                            forces = np.array(tmpforces).T

                        else:
                            forces = np.concatenate((forces, np.array(tmpforces).T), axis=0)

                        tmpforces = []

                    line = skiplines(f)

            #
            # Forces
            #
            if "Updating gradient with" in line:
                forces = []
                line = skiplines(f, 4)

                while '--------------------' not in line:
                    data = list(map(float, line.split()[1:]))
                    forces.append(data)
                    line = skiplines(f)

    return np.array(forces)


def parseforces_G09(filename):
    '''Parses G09 forces logfile to retrieve only forces.'''

    with open(filename) as f:

        forces = []

        for line in f:

            #
            # Forces
            #
            if "Forces (Hartrees/Bohr)" in line:
                line = skiplines(f, 2)

                while len(line.split()) == 5:
                    data = list(map(float, line.split()[2:]))
                    forces.append(data)

                    line = skiplines(f)


    return np.array(forces)


def guess(filename):
    '''Returns the correct class needed to parse filename, if it exists.'''

    #
    # Dictionary of unique sentences in QM packages output files to guess
    # the correct parser to use
    #
    filetypes = {}
    
    filetypes["Gaussian(R) 09 "] = "G09"
    filetypes["A Quantum Leap Into The Future Of Chemistry"] = "QChem"

    filetype = None
    done = False
    with open(filename) as f:

        for line in f:
            for sentence in filetypes.keys():

                if sentence in line:
                    filetype = filetypes[sentence]
                    done = True
                    break

            # once the type has been identified, exit the cycles
            if done:
                break

    if not filetype:
        print(" %s" % filename)
        print(" File type not known")
        sys.exit()

    return filetype


if __name__ == '__main__':
    pass
