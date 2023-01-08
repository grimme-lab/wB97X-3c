# ORCA4wB97X-3c
Fortran script for setting up a wB97X-3c calculation with ORCA.

To set up an input file for wB97X-3c, you have to execute the binary in a directory with a molecular structure file (can be either .xyz, coord, or other usual formats (see mctc-lib for possible formats).

You need the files:
- ".basis_vDZP" and ".ecp" in your $HOME

Example usage:

o4wb3c --struc ch4.xyz

See the "-help" flag for further possibilities.

WARNING: Due to the use of a deprecated I/O routine, the use of 'coord' files with scientific notation is not possible.
