# Introduction 

The first step for installing NMRTK is cloning or obtaining the source code
from github.  NMRTK also requires GNU Scientific Library (GSL) and Genetic
Algorithms Library (GAlib) from MIT.  Both of them can be installed into
Linux system-widely (such as /usr/local). For users without administration
permission, it is more convenient to install both of them under the NMRTK
package as shown below. For other customized installation, users need to go
through the instructions under corresponding packages.

# Install GSL lib 

On Ubuntu, simply type "sudo apt-get install libgsl-dev".  Similar commands
probably work on other Linux distributions.

Or, the GSL lib can be downloaded from the                                      
[webpage](https://www.gnu.org/software/gsl/). Although the tested versions  
are 2.6 and 1.14, other versions should work well. To install GSL, please   
follow the steps as listed below:                                           

1) move the GSL source code and untar it under the NMRTK folder such as
`nmrtk_git`.

   `tar -xzvf XXX.tar.gz`

2) Change to the folder of GSL source code such as (gsl-2.6), compile       
and install it under NMRTK folder such as `$HOME/software/nmrtk` using      
following commands.                                                         

    `./configure --prefix=$HOME/software/nmrtk_git`

    `make` 

    `make install` 

# Install GAlib

On Ubuntu, simply type "sudo apt-get install libga-dev".  Similar commands
probably work on other Linux distributions.

Or, the 2.4.7 version of [GAlib](http://web.mit.edu/galib/www/GAlib.html) has
been included in this NMRTK package. Brief instruction on compiling and
installing is listed as below:

1) Change to the `galib247` under the NMRTK folder such as
`$HOME/software/nmrtk_git`.

2) Modify the `DESTDIR` variable in the `makevars` file for the destination
folder such as
   
   `DESTDIR=$HOME/software/nmrtk_git`
    
3) Add `-fpermissive` option to the `CXXFLAGS` variable to allow old C++
code, such as
   
    `CXXFLAGS    = -O3 -Wall -fpermissive` 
    
4) Compile and install GAlib by two commands as below:
   
   `make all` 
   
   `make install`

# Install NMRTK

1) Change to the `src` folder of NMRTK and set the `BASEDIR`, `GALIBS` and
 `GSLLIBS`  variables in `Makefile` (for sudo install) or `Makefile_local` 
 (for customized install) to match your system
   
2) if you installed GSL lib and GAlib using sudo command, run `make all` to compile and install NMRTK. 

3) Or if installed GSL lib and GAlib in a customized way, run `make all -f Makefile_local` 
   to compile and install NMRTK. 

The related executable commands can be found under the `bin` folder of NMR
package.
