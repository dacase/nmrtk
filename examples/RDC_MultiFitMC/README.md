This is a simple example to run `multiFitMC` command for fitting
experimental RDC dataset to an ensemble structures using Monte Carlo
algorithm. Similar to `multiFitGA` We assume that a unique ensemble of structures
has been obtained from conformational search or molecular dynamics simulation, 
and the RDC values for each unique structure have been calculated by other programs 
such as PALES. Generally, `multiFitMC` is less efficient to converge to a good solution. 

## Prepare input files
We need to prepare three input files to find the population weights of the
best ensemble structures as below:

1) Put the experimental RDC dataset in the file `RDCs_exp.dat`, similar to run `multiFitGA`.

2) similar to run `multiFitGA`, sort and put the calculated RDC data from all unique 
conformers into the file, `RDCs_sorted.dat`.

3) Specify the fitting parameters in the file, `multiFitMC.inp`. The
detailed meaning of each parameter can be found after the `//` symbol. This
file can be used as a template for other systems. Only necessary values
after the colon (:) symbols need to be updated.

## Fit to multiple structures
After the above input files are ready, we can run the fitting command as below:

   `path_to_executable/multiFitMC multiFitMC.inp  >& screen.out &`

When the fitting procedure is finished, we will find several output files
with the names of `*.out`.  The last column of `probBest.out` shows the
population weights of the best structure ensemble, with the first three
columns are the same as `RDCs_sorted.dat`. The ensemble averaged RDC values
are output as the second column of `nmrBest.out`.
