require 'mkmf'
# Add autoitx.lib

$libs = append_library($libs,"necpp")
#$libs = append_library($libs,"stdc++")

create_makefile('necpp')