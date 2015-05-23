# Ruby necpp module

This module allows you to do antenna simulations in Ruby using the nec2++ antenna
simulation package. This is a wrapper using SWIG of the C interface, so the syntax
is quite simple. Have a look at the file test.rb, for an example of how this 
library can be used.

### Author

Tim Molteno. tim@physics.otago.ac.nz

## Instructions

To use this ruby module, you must have the necpp library installed on your system. This can
be installed in the main part of the necpp code distribution.

### To generate this module

You should install SWIG (on Debian 'aptitude install swig ruby-dev'), and then
issue the following commands

    cd ext/necpp
    swig -v -I../../../src -c++ -ruby necpp.i
    ruby extconf.rb
    make
    sudo make install

Alternatively use the build.sh script.

    cd ext/necpp
    ./build.sh
      
Then test with 

    ruby ../example/test.rb

## Genetic Optimization of Antennas

Have a look in the directory genetic_optimizer for some Ruby code that optimizes antenna designs.
