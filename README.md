> [!WARNING]  
> ωB97X-3c is now natively available in `ORCA-6.0.0` (see [here](https://www.faccts.de/docs/orca/6.0/manual/contents/detailed/model.html#omegab97x-3c-a-composite-range-separated-hybrid-dft-method-with-a-molecule-optimized-polarized-valence-double-zeta-basis-set) for details and download ORCA [here](https://orcaforum.kofo.mpg.de/app.php/dlext/?cat=23)).
> This program is therefore considered obsolete and **_should not be used_** anymore for standard purposes.


## ωB97X-3c for ORCA
This is a Fortran script for setting up a ωB97X-3c calculation with `ORCA>=5.0.3`. When using the script or the original implementation, please cite the original [publication](https://pubs.aip.org/aip/jcp/article-abstract/158/1/014103/2867476/B97X-3c-A-composite-range-separated-hybrid-DFT):

> M. Müller; A. Hansen; S. Grimme; _J. Chem. Phys._ **158**, 014103 (2023)

```
@article{muller_b97x-3c_2023,
	title = {ω{B97X}-3c: {A} composite range-separated hybrid {DFT} method with a molecule-optimized polarized valence double- zeta basis set},
	volume = {158},
	issn = {10897690},
	url = {https://aip.scitation.org/doi/abs/10.1063/5.0133026},
	doi = {10.1063/5.0133026},
	abstract = {A new composite density functional theory (DFT) method is presented. It is based on ωB97X-V as one of the best-performing density functionals for the GMTKN55 thermochemistry database and completes the family of "3c"methods toward range-separated hybrid DFT. This method is consistently available for all elements up to Rn (Z = 1-86). Its further key ingredients are a polarized valence double-ζ (vDZP) Gaussian basis set, which was fully optimized in molecular DFT calculations, in combination with large-core effective core potentials and a specially adapted D4 dispersion correction. Unlike most existing double-ζ atomic orbital sets, vDZP shows only small basis set superposition errors (BSSEs) and can compete with standard sets of triple-ζ quality. Small residual BSSE effects are efficiently absorbed by the D4 damping scheme, which overall eliminates the need for an explicit treatment or empirical corrections for BSSE. Thorough tests on a variety of thermochemistry benchmark sets show that the new composite method, dubbed ωB97X-3c, is on par with or even outperforms standard hybrid DFT methods in a quadruple-zeta basis set at a small fraction of the computational cost. Particular strengths of this method are the description of non-covalent interactions and barrier heights, for which it is among the best-performing density functionals overall.},
	number = {1},
	urldate = {2023-03-19},
	journal = {Journal of Chemical Physics},
	author = {Müller, Marcel and Hansen, Andreas and Grimme, Stefan},
	month = jan,
	year = {2023},
	pmid = {36610980},
	note = {Publisher: AIP Publishing LLCAIP Publishing},
	pages = {014103},
}
```

### Release version (recommended)
The use of the release binary [`o4wb3c`] is recommended. The binary has to be added to a location belonging to your `$PATH` variable.

### Building with Fortran package Manager
You can use the Fortran Package Manager (https://github.com/fortran-lang/fpm) in version 0.9.0 or higher to build the project.
To install the project in your prefered path, just use 
```
fpm install -profile release -prefix [path]
```
More information about FPM can be found in the respective documentation.

### Usage
To set up an input file for ωB97X-3c, you have to execute the binary in a directory with a molecular structure file (can be either .xyz, coord, or common formats (see `mctc-lib` (https://github.com/grimme-lab/mctc-lib) for possible formats).

You need the files:
- `.basis_vDZP` and `.ecp_vDZP` (can be downloaded [here](https://github.com/grimme-lab/ORCA4wB97X-3c/blob/main/basis_vDZP) and here [here](https://github.com/grimme-lab/ORCA4wB97X-3c/blob/main/ecp_vDZP)) in your `$HOME` (ATTENTION: the file names of the basis set and ECP files do not yet contain the `.`) or you specify an individual location of the files (see example below or press `--help`).

Example use cases:

```
o4wb3c --struc coord.benzene
o4wb3c --struc ch3.xyz --basisfile /home/$USER/basissets/basis_vDZP_TM.txt --chrg -1
```
If no `--struc` file is explicitly given, `o4wb3c` assumes a `coord` file.

See the "-help" flag for further input possibilities.
If you should observe instabilities with the `PModel` guess in ORCA, try to use `o4wb3c` together with the `--suggestedguess` flag or provide an individual guess option with `--guess`.
