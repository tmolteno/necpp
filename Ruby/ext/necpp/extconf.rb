require 'mkmf'

#NECPP_SRC = "/home/tim/github/necpp/src"
NECPP_SRC = "../../../src"
$CFLAGS << " -I#{NECPP_SRC}"
$CPPFLAGS << " -I#{NECPP_SRC}"

dir_config("lapack", "/usr/include/atlas", "/usr/lib/atlas-base/atlas")
#unless have_header("clapack.h") and have_library("lapack", "clapack_zgetrf", %w(clapack.h))
unless have_header("clapack.h") 
  abort "lapack is missing. please install ATLAS"
end
unless have_library("lapack", "clapack_zgetrf", "clapack.h")
  abort "lapack is missing. please install ATLAS or 
  issue the command update-alternatives liblapack.so"
end

dir_config("necpp", "/usr/local/include", "/usr/local/lib")
unless have_header("libnecpp.h") and find_library("necpp", "nec_create")
   abort "necpp library is missing."
end

#$LDFLAGS << " -Wl,-rpath,/usr/lib/lapack"

create_makefile('necpp')
