require 'mkmf'

# NECPP source directory
NECPP_SRC = "../../../src"
$CFLAGS << " -I#{NECPP_SRC}"
$CPPFLAGS << " -I#{NECPP_SRC} -I#{NECPP_SRC}/eigen3"

# Eigen is bundled — no external LAPACK/BLAS needed
dir_config("necpp", "/usr/local/include", "/usr/local/lib")
unless have_header("libnecpp.h") and find_library("necpp", "nec_create")
   abort "necpp library is missing. Build with: cd ../../.. && make"
end

create_makefile('necpp')
