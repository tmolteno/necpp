# Python necpp module

This module allows you to do antenna simulations in Python using the nec2++ antenna
simulation package. This is a wrapper using SWIG of the C interface, so the syntax
is quite simple. Have a look at the file test.py, for an example of how this 
library can be used.

### Author

Tim Molteno. tim@physics.otago.ac.nz

## Instructions

To use this python module, you must have the necpp library installed on your system. This can
be installed in the main part of the necpp code distribution.

### To generate this module

You should install SWIG (on Debian 'aptitude install swig python-dev'), and then
issue the following commands


    ./build.sh
      
Then test with 

    python ../example/test.py

