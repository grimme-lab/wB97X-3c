# ORCA4wB97X-3c
This is a Fortran script for setting up a ωB97X-3c calculation with ORCA.

### Release version (recommended)
The use of the release binary [`o4wb3c`] is recommended. The binary has to be added to a location belonging to your `$PATH` variable.

### Building with Fortran package Manager
You can use the Fortran Package Manager (https://github.com/fortran-lang/fpm) to build the project.
To install the project in your prefered path, just use 
```
fpm install -profile release -prefix [path]
```
More information about FPM can be found in the respective documentation.

### Usage
To set up an input file for ωB97X-3c, you have to execute the binary in a directory with a molecular structure file (can be either .xyz, coord, or common formats (see `mctc-lib` (https://github.com/grimme-lab/mctc-lib) for possible formats).

You need the files:
- `basis_vDZP` and `ecp` in your `$HOME` or you specify an individual location of the files (see example below).

Example use cases:

```
o4wb3c --struc coord.benzene
o4wb3c --struc ch3.xyz --basisfile /home/$USER/basissets/basis_vDZP_TM.txt --chrg -1
```
If no `--struc` file is explicitly given, `o4wb3c` assumes a `coord` file.

See the "-help" flag for further input possibilities.
