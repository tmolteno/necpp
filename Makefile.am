# not a GNU package. You can remove this line, if
# have all needed files, that a GNU package needs

AUTOMAKE_OPTIONS = foreign 1.4
ACLOCAL_AMFLAGS = -I config/m4
SUBDIRS = src

EXTRA_DIST = \
	example \
	html \
	win32

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = necpp.pc

dist-hook:
	rm -rf `find $(distdir)testharness/data -name *.out*`
	rm -rf `find $(distdir)testharness/data -name *.diff*`
