This is a simple example to run `multiFitGA` command for fitting
experimental RDC dataset to an ensemble structures using genetic
algorithm. We assume that a unique ensemble of structures has been obtained
from conformational search or molecular dynamics simulation, and the RDC
values for each unique structure have been calculated by other programs such
as PALES.

## Prepare input files
We need to prepare three input files to find the population weights of the
best ensemble structures as below:

1) Put the experimental RDC dataset in the file, `RDCs_exp.dat`, with one
value per row.

2) Sort and put the calculated RDC data from all unique conformers into
the file, `RDCs_sorted.dat`, with one conformer per row. The first column
denotes the row ID, the second column represents the conformer ID in the
original conformational trajectory, the third column is the R factor of
calculated RDCs from corresponding conformer (used to sort all unique
conformers), and the rest values are calculated RDCs.

3) Specify the fitting parameters in the file, `multiFitGA.inp`. The
detailed meaning of each parameter can be found after the `//` symbol. This
file can be used as a template for other systems. Only necessary values
after the colon (:) symbols need to be updated.

## Fit to multiple structures
After the above input files are ready, we can run the fitting command as below:

   multiFitGA multiFitGA.inp  > screen.out &`

When the fitting procedure is finished, we will find several output files
with the names of `*.out`.  The last column of `probBest.out` shows the
population weights of the best structure ensemble, with the first three
columns are the same as `RDCs_sorted.dat`. The ensemble averaged RDC values
are output as the third column of `nmrBest.out`, and the synthesized
experimental RDC values with random errors (added to the experimental values
in the first column) are listed in the second column.

 
