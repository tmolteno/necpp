#!/bin/sh
# Script to build the nec2++ ruby module.
# Change the RUBY environment variable 
# as appropriate for your system
#RUBY=ruby1.9.1
RUBY=ruby
swig -v -I../../../src -c++ -ruby necpp.i
${RUBY} extconf.rb
make V=1
sudo make install
