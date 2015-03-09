# not a GNU package. You can remove this line, if
# have all needed files, that a GNU package needs

AUTOMAKE_OPTIONS = foreign 1.4
ACLOCAL_AMFLAGS = -I m4
SUBDIRS = src

EXTRA_DIST = \
	example \
	html \
	test_data \
	win32 \
	Makefile.test

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = necpp.pc

dist-hook:
	rm -rf `find $(distdir)testharness/data -name *.out*`
	rm -rf `find $(distdir)testharness/data -name *.diff*`
