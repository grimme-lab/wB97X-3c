# ORCA4wB97X-3c
This is a Fortran script for setting up a wB97X-3c calculation with ORCA.

### Building with Fortran package Manager
You can use the Fortran Package Manager (https://github.com/fortran-lang/fpm) to build the project.
To install the project in your prefered path, just use 
```
fpm install -profile release -prefix [path]
```
More information about FPM can be found in the respective documentation.

To set up an input file for wB97X-3c, you have to execute the binary in a directory with a molecular structure file (can be either .xyz, coord, or other usual formats (see mctc-lib for possible formats).

You need the files:
- ".basis_vDZP" and ".ecp" in your $HOME

Example usage:

```
o4wb3c --struc ch4.xyz
```

See the "-help" flag for further possibilities.
