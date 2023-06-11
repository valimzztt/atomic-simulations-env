# atomic-simulations-env
Updated version of the ASE (Atomic Simulation Environment) Python package. 
The atomic simulations environment provides modules for performing many standard simulation tasks such as structure optimization, molecular dynamics, handling of constraints and performing nudged elastic band calculations.
## Changes from original ASE package
### 1: Hubbard Parameter: 
The user can now provide the Hubbard correction for single atomic species. Each element in the array that serves as the key of the hubbard_corrections key represents a single Hubbard correction. The first element of each triple represents the atomic species to which the Hubbard correction is applied, the second element represents the orbital, and the third element represents the energy of the correction in electron volts (eV).
The new input syntax for the input_data variable when using DFT+U (Dudarev’s formulation) for (Ti, Zr)O will, for example, be: 
```
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
```
