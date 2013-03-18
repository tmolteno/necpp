require 'mkmf'


dir_config("lapack", "/usr/include/atlas", "/usr/lib/lapack")
unless have_header("clapack.h") and find_library("lapack", "clapack_zgetrf")
  abort "lapack is missing. please install ATLAS"
end

dir_config("necpp", "/usr/local/include", "/usr/local/lib")
unless have_header("libnecpp.h") and find_library("necpp", "nec_create")
   abort "necpp library is missing."
end

$LDFLAGS << " -Wl,-rpath,/usr/lib/lapack"

create_makefile('necpp')
