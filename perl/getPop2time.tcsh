#!/bin/tcsh
set ifileS = $1
set ifileE = $2
set inpfile = $3 
set outfile = $4 

set ifileCur=$ifileS

while ( $ifileCur <= $ifileE )

  probIntRect.pl 2 $inpfile fes.dat.$ifileCur 2.479
  cat pops.out >> $outfile  
  echo finishing $ifileCur
  @ ifileCur = $ifileCur + 1

end

