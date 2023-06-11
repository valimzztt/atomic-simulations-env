import pytest
import os
from ..my_io.espresso import *
import numpy as np
from ase import Atoms
from ase.constraints import FixAtoms
from ase.calculators.espresso import Espresso
from ase.constraints import FixCartesian
from ase.calculators.espresso import Espresso
# Importing the StringIO module.
from io import StringIO

# Test an instance where Hubbard correction is not inserted 
def test_write_espresso_in_nohubbard():
    # Create a file-like object to capture the output
    fd = StringIO()

    # Create an example Atoms object
    atoms = Atoms(symbols=['Ti', 'Zr', 'O'], positions=[[9.44200000000000,0.00000000000000,0.00000000000000],
                                                       [0.00000000000000, 5.74540000000000, 5.01400000000000],
                                                       [0.00000000000000, 0.00000000000000, 5.01400000000000]])
    curr_directory = os.getcwd()
    pseudo_dir = os.path.join(curr_directory, 'pseudos')
    # Set up other necessary inputs
    pseudopotentials = {'Ti': 'Ti.pbesol-spn-kjpaw_psl.1.0.0.UPF',
                        'Zr': 'Zr.pbesol-spn-kjpaw_psl.1.0.0.UPF',
                        'O': 'O.pbesol-n-kjpaw_psl.1.0.0.UPF'}
    input_data = {
        'control': {
           'calculation': 'relax',
           'restart_mode': 'from_scratch',
           'pseudo_dir': pseudo_dir,
           'prefix': 'tutorial_hubbard',        
            },
        'system': {
           'ecutwfc': 40.424222147,
           'occupations': 'smearing',
           'degauss': 0.0146997171,
           'lda_plus_u' : '.FALSE.',
          }
    }
    koffset = (0, 0, 0)

    # Call the method being tested
    write_espresso_in(fd, atoms, input_data=input_data, pseudopotentials=pseudopotentials, koffset=koffset)

    # Get the captured output
    output = fd.getvalue()
    print(output)

    # Assert expected output
    assert 'CONTROL' in output
    # Assert statements
    assert atoms.get_chemical_symbols() == ['Ti', 'Zr', 'O']
    assert len(atoms) == 3
    assert atoms.get_positions()[0].tolist() == [9.44200000000000, 0.00000000000000, 0.00000000000000]
    assert atoms.get_positions()[1].tolist() == [0.00000000000000, 5.74540000000000, 5.01400000000000]
    assert atoms.get_positions()[2].tolist() == [0.00000000000000, 0.00000000000000, 5.01400000000000]

    assert curr_directory == os.getcwd()
    assert pseudo_dir == os.path.join(curr_directory, 'pseudos')

    assert pseudopotentials['Ti'] == 'Ti.pbesol-spn-kjpaw_psl.1.0.0.UPF'
    assert pseudopotentials['Zr'] == 'Zr.pbesol-spn-kjpaw_psl.1.0.0.UPF'
    assert pseudopotentials['O'] == 'O.pbesol-n-kjpaw_psl.1.0.0.UPF'

    assert input_data['control']['calculation'] == 'relax'
    assert input_data['control']['restart_mode'] == 'from_scratch'
    assert input_data['control']['pseudo_dir'] == pseudo_dir
    assert input_data['control']['prefix'] == 'tutorial_hubbard'

    assert input_data['system']['ecutwfc'] == 40.424222147
    assert input_data['system']['occupations'] == 'smearing'
    assert input_data['system']['degauss'] == 0.0146997171
    assert input_data['system']['lda_plus_u'] == '.FALSE.'

    assert koffset == (0, 0, 0)
    # Clean up the file-like object
    fd.close()

# Run the test
test_write_espresso_in_nohubbard()

# Test an instance where Hubbard correction is inserted 
def test_write_espresso_in_hubbard():
    # Create a file-like object to capture the output
    fd = StringIO()

    # Create an example Atoms object
    atoms = Atoms(symbols=['Ti', 'Zr', 'O'], positions=[[9.44200000000000,0.00000000000000,0.00000000000000],
                                                       [0.00000000000000, 5.74540000000000, 5.01400000000000],
                                                       [0.00000000000000, 0.00000000000000, 5.01400000000000]])
    curr_directory = os.getcwd()
    pseudo_dir = os.path.join(curr_directory, 'pseudos')
    # Set up other necessary inputs
    pseudopotentials = {'Ti': 'Ti.pbesol-spn-kjpaw_psl.1.0.0.UPF',
                        'Zr': 'Zr.pbesol-spn-kjpaw_psl.1.0.0.UPF',
                        'O': 'O.pbesol-n-kjpaw_psl.1.0.0.UPF'}
    input_data = {
        'control': {
           'calculation': 'relax',
           'restart_mode': 'from_scratch',
           'pseudo_dir': pseudo_dir,
           'prefix': 'tutorial_hubbard',        
            },
        'system': {
           'ecutwfc': 40.424222147,
           'occupations': 'smearing',
           'degauss': 0.0146997171,
           'lda_plus_u' : '.TRUE.',
          },
        'HUBBARD': {
            'hubbard_corrections': [('Ti', "3d", 5.0), ('Zr',"4d", 3.0)]
        },
    }
    koffset = (0, 0, 0)

    # Call the method being tested
    write_espresso_in(fd, atoms, input_data=input_data, pseudopotentials=pseudopotentials, koffset=koffset)

    # Get the captured output
    output = fd.getvalue()
    print(output)

    # Assert expected output
    assert 'CONTROL' in output
    # Assert statements
    assert atoms.get_chemical_symbols() == ['Ti', 'Zr', 'O']
    assert len(atoms) == 3
    assert atoms.get_positions()[0].tolist() == [9.44200000000000, 0.00000000000000, 0.00000000000000]
    assert atoms.get_positions()[1].tolist() == [0.00000000000000, 5.74540000000000, 5.01400000000000]
    assert atoms.get_positions()[2].tolist() == [0.00000000000000, 0.00000000000000, 5.01400000000000]

    assert curr_directory == os.getcwd()
    assert pseudo_dir == os.path.join(curr_directory, 'pseudos')

    assert koffset == (0, 0, 0)
    # Perform assertions to validate the written input
    assert 'U Ti-3d 5.0' in output
    assert 'U Zr-4d 3.0' in output
    # Clean up the file-like object
    fd.close()

# Run the test
test_write_espresso_in_hubbard()