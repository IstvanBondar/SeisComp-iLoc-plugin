Installing sciLoc library
=========================

The Makefile checks for the existence of dependencies (lapack);
installs the RSTT libraries;
installs the sciLoc library;
compiles the test programs.

Just run make.

Tested on:
  Mac Os 10.12.6,
  Ubuntu 16.04.4
  CentOS 6.9

On Linux, add $HOME/lib to your LD_LIBRARY_PATH:

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/lib

For the SeisComp4 interface see sciLocInterface.h in the include directory.

Istvan Bondar
ibondar2014@gmail.com



