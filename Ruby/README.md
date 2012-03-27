#Ruby necpp module

This module allows you to do antenna simulations in Ruby using the nec2++ antenna
simulation package. This is a wrapper using SWIG of the C interface, so the syntax
is quite simple. Have a look at the file test.rb, for an example of how this 
library can be used.

###Author

Tim Molteno. tim@physics.otago.ac.nz

##Instructions

To use this library, you must have the necpp library installed on your system:

On Debian based systems:

	aptitude install ruby-dev swig

You should have built the nec2++ distribution and installed it.

###To generate this module

You should install SWIG (on Debian 'aptitude install swig ruby-dev'), and then
issue the following commands

	swig -v -c++ -ruby necpp.i
	ruby extconf.rb
	make
	sudo make install

Then test with 

	ruby test.rb
