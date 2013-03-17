require 'mkmf'

$CFLAGS << " -O3"

$libs = append_library($libs,"necpp")
$libs = append_library($libs,"lapack")

dir_config("necpp", "/usr/local/include", "/usr/local/lib")

# unless find_library("necpp")
#   abort "necpp library is missing."
# end

dir_config("lapack", "/usr/include/atlas", "/usr/local/lib")
unless have_header("clapack.h") and find_library("lapack", "clapack_dgetrf", "/usr/lib/lapack")
  abort "lapack is missing. please install ATLAS"
end

$LDFLAGS << " -Wl,-rpath,/usr/lib/lapack"

create_makefile('necpp')
