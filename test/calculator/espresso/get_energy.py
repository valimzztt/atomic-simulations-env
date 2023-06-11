import re
import os 
import ase.io.espresso
curr_directory  = os.getcwd()
filepath = curr_directory + "/espresso.pwo"

import os
import operator as op
import re
import warnings
from collections import OrderedDict
from os import path

import numpy as np

from ase.atoms import Atoms
from ase.calculators.singlepoint import (SinglePointDFTCalculator,
                                         SinglePointKPoint)
from ase.calculators.calculator import kpts2ndarray, kpts2sizeandoffsets
from ase.dft.kpoints import kpoint_convert
from ase.constraints import FixAtoms, FixCartesian
from ase.data import chemical_symbols, atomic_numbers
from ase.units import create_units
from ase.utils import iofunction

filename = "espresso_nonhubbard.pwo"
# Quantum ESPRESSO uses CODATA 2006 internally
units = create_units('2006')


# Section identifiers
_PW_START = 'Program PWSCF'
_PW_END = 'End of self-consistent calculation'
_PW_CELL = 'CELL_PARAMETERS'
_PW_POS = 'ATOMIC_POSITIONS'
_PW_MAGMOM = 'Magnetic moment per site'
_PW_FORCE = 'Forces acting on atoms'
_PW_TOTEN = '!    total energy'
_PW_STRESS = 'total   stress'
_PW_FERMI = 'the Fermi energy is'
_PW_HIGHEST_OCCUPIED = 'highest occupied level'
_PW_HIGHEST_OCCUPIED_LOWEST_FREE = 'highest occupied, lowest unoccupied level'
_PW_KPTS = 'number of k points='
_PW_BANDS = _PW_END
_PW_BANDSTRUCTURE = 'End of band structure calculation'


class Namelist(OrderedDict):
    """Case insensitive dict that emulates Fortran Namelists."""
    def __contains__(self, key):
        return super(Namelist, self).__contains__(key.lower())

    def __delitem__(self, key):
        return super(Namelist, self).__delitem__(key.lower())

    def __getitem__(self, key):
        return super(Namelist, self).__getitem__(key.lower())

    def __setitem__(self, key, value):
        super(Namelist, self).__setitem__(key.lower(), value)

    def get(self, key, default=None):
        return super(Namelist, self).get(key.lower(), default)
# TODO: index -1 special case?
# Index all the interesting points
indexes = {
    _PW_START: [],
    _PW_END: [],
    _PW_CELL: [],
    _PW_POS: [],
    _PW_MAGMOM: [],
    _PW_FORCE: [],
    _PW_TOTEN: [],
    _PW_STRESS: [],
    _PW_FERMI: [],
    _PW_HIGHEST_OCCUPIED: [],
    _PW_HIGHEST_OCCUPIED_LOWEST_FREE: [],
    _PW_KPTS: [],
    _PW_BANDS: [],
    _PW_BANDSTRUCTURE: [],
}


def parse_position_line(line):
    """Parse a single line from a pw.x output file.

    The line must contain information about the atomic symbol and the position,
    e.g.

    995           Sb  tau( 995) = (   1.4212023   0.7037863   0.1242640  )

    Parameters
    ----------
    line : str
        Line to be parsed.

    Returns
    -------
    sym : str
        Atomic symbol.
    x : float
        x-position.
    y : float
        y-position.
    z : float
        z-position.
    """
    pat = re.compile(r'\s*\d+\s*(\S+)\s*tau\(\s*\d+\)\s*='
                     r'\s*\(\s*(\S+)\s+(\S+)\s+(\S+)\s*\)')
    match = pat.match(line)
    assert match is not None
    sym, x, y, z = match.group(1, 2, 3, 4)
    return sym, float(x), float(y), float(z)
def label_to_symbol(label):
    """Convert a valid espresso ATOMIC_SPECIES label to a
    chemical symbol.

    Parameters
    ----------
    label : str
        chemical symbol X (1 or 2 characters, case-insensitive)
        or chemical symbol plus a number or a letter, as in
        "Xn" (e.g. Fe1) or "X_*" or "X-*" (e.g. C1, C_h;
        max total length cannot exceed 3 characters).

    Returns
    -------
    symbol : str
        The best matching species from ase.utils.chemical_symbols

    Raises
    ------
    KeyError
        Couldn't find an appropriate species.

    Notes
    -----
        It's impossible to tell whether e.g. He is helium
        or hydrogen labelled 'e'.
    """

    # possibly a two character species
    # ase Atoms need proper case of chemical symbols.
    if len(label) >= 2:
        test_symbol = label[0].upper() + label[1].lower()
        if test_symbol in chemical_symbols:
            return test_symbol
    # finally try with one character
    test_symbol = label[0].upper()
    if test_symbol in chemical_symbols:
        return test_symbol
    else:
        raise KeyError('Could not parse species from label {0}.'
                       ''.format(label))


def parse_pwo_start(lines, index=0):
    """Parse Quantum ESPRESSO calculation info from lines,
    starting from index. Return a dictionary containing extracted
    information.

    - `celldm(1)`: lattice parameters (alat)
    - `cell`: unit cell in Angstrom
    - `symbols`: element symbols for the structure
    - `positions`: cartesian coordinates of atoms in Angstrom
    - `atoms`: an `ase.Atoms` object constructed from the extracted data

    Parameters
    ----------
    lines : list[str]
        Contents of PWSCF output file.
    index : int
        Line number to begin parsing. Only first calculation will
        be read.

    Returns
    -------
    info : dict
        Dictionary of calculation parameters, including `celldm(1)`, `cell`,
        `symbols`, `positions`, `atoms`.

    Raises
    ------
    KeyError
        If interdependent values cannot be found (especially celldm(1))
        an error will be raised as other quantities cannot then be
        calculated (e.g. cell and positions).
    """
    # TODO: extend with extra DFT info?

    info = {}

    for idx, line in enumerate(lines[index:], start=index):
        if 'celldm(1)' in line:
            # celldm(1) has more digits than alat!!
            info['celldm(1)'] = float(line.split()[1]) * units['Bohr']
            info['alat'] = info['celldm(1)']
        elif 'number of atoms/cell' in line:
            info['nat'] = int(line.split()[-1])
        elif 'number of atomic types' in line:
            info['ntyp'] = int(line.split()[-1])
        elif 'crystal axes:' in line:
            info['cell'] = info['celldm(1)'] * np.array([
                [float(x) for x in lines[idx + 1].split()[3:6]],
                [float(x) for x in lines[idx + 2].split()[3:6]],
                [float(x) for x in lines[idx + 3].split()[3:6]]])
        elif 'positions (alat units)' in line:
            info['symbols'], info['positions'] = [], []

            for at_line in lines[idx + 1:idx + 1 + info['nat']]:
                sym, x, y, z = parse_position_line(at_line)
                info['symbols'].append(label_to_symbol(sym))
                info['positions'].append([x * info['celldm(1)'],
                                          y * info['celldm(1)'],
                                          z * info['celldm(1)']])
            # This should be the end of interesting info.
            # Break here to avoid dealing with large lists of kpoints.
            # Will need to be extended for DFTCalculator info.
            break

    # Make atoms for convenience
    info['atoms'] = Atoms(symbols=info['symbols'],
                          positions=info['positions'],
                          cell=info['cell'], pbc=True)

    return info


with open(filename,"r") as file:
    # For each iteration in the DFT calculations we get new updated energy values
    energies = []
    fil = file    
    pwo_lines = fil.readlines()
    print(pwo_lines)
    results_required = True
    # TODO: index -1 special case?
    # Index all the interesting points

    indexes = {
        _PW_START: [],
        _PW_END: [],
        _PW_CELL: [],
        _PW_POS: [],
        _PW_MAGMOM: [],
        _PW_FORCE: [],
        _PW_TOTEN: [],
        _PW_STRESS: [],
        _PW_FERMI: [],
        _PW_HIGHEST_OCCUPIED: [],
        _PW_HIGHEST_OCCUPIED_LOWEST_FREE: [],
        _PW_KPTS: [],
        _PW_BANDS: [],
        _PW_BANDSTRUCTURE: [],
    }
    index = -1
    for idx, line in enumerate(pwo_lines):
        for identifier in indexes:
            if identifier in line:
                indexes[identifier].append(idx)

    all_config_indexes = sorted(indexes[_PW_START] +
                                indexes[_PW_POS])

    if results_required:
        print("HELLO")
        results_indexes = sorted(indexes[_PW_TOTEN] + indexes[_PW_FORCE] +
                                 indexes[_PW_STRESS] + indexes[_PW_MAGMOM] +
                                 indexes[_PW_BANDS] +
                                 indexes[_PW_BANDSTRUCTURE])
        print(results_indexes)
        # Prune to only configurations with results data before the next
        # configuration
        results_config_indexes = []
        for config_index, config_index_next in zip(
                all_config_indexes,
                all_config_indexes[1:] + [len(pwo_lines)]):
            if any([config_index < results_index < config_index_next
                    for results_index in results_indexes]):
                results_config_indexes.append(config_index)


        # slice from the subset
        print(results_config_indexes)
        image_indexes = results_config_indexes
    else:
        image_indexes = all_config_indexes[index]

    # Extract initialisation information each time PWSCF starts
    # to add to subsequent configurations. Use None so slices know
    # when to fill in the blanks.
    pwscf_start_info = dict((idx, None) for idx in indexes[_PW_START])
    for image_index in image_indexes:
        # Find the nearest calculation start to parse info. Needed in,
        # for example, relaxation where cell is only printed at the
        # start.
        if image_index in indexes[_PW_START]:
            prev_start_index = image_index
        else:
            # The greatest start index before this structure
            prev_start_index = [idx for idx in indexes[_PW_START]
                                if idx < image_index][-1]

        # add structure to reference if not there
        if pwscf_start_info[prev_start_index] is None:
            pwscf_start_info[prev_start_index] = parse_pwo_start(
                pwo_lines, prev_start_index)

        # Get the bounds for information for this structure. Any associated
        # values will be between the image_index and the following one,
        # EXCEPT for cell, which will be 4 lines before if it exists.
        for next_index in all_config_indexes:
            if next_index > image_index:
                break
        else:
            # right to the end of the file
            next_index = len(pwo_lines)
        

    for line in file:
        pattern="!"
        if re.search(pattern, line):
            energies.append(re.findall(r'-?\d+', line))
    # Get the last calculated energy value before the calculation converges
    a = energies[0][0]
    b = energies[0][1]
    d = float(f'{a}.{b}')  
    print("The final energies is ", d)

    calc = SinglePointDFTCalculator(structure, energy=energy,
                                    free_energy=energy,
                                    forces=forces, stress=stress,
                                    magmoms=magmoms, efermi=efermi,
                                    ibzkpts=None)
    calc.kpts = kpts
    structure.calc = calc

    yield structure

