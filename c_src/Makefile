#Makefile for nec2c     21 Aug 2003

SHELL = /bin/sh
PROJECT = nec2c
BINDIR  = /usr/local/bin
CC = gcc -Wall -g3

objects = calculations.o fields.o geometry.o ground.o input.o \
	  main.o matrix.o misc.o network.o radiation.o somnec.o

$(PROJECT) : $(objects)
	    $(CC) -lm -o $(PROJECT) $(objects)

$(objects) : nec2c.h

nec2dx : 
	g77 -o nec2dx nec2dx.f
	install -m 755 --strip nec2dx $(BINDIR)

install : $(PROJECT)
	  install -m 755 --strip $(PROJECT) $(BINDIR)

.PHONY : distclean
distclean  :
	-rm -f *.o *~ $(PROJECT) nec2dx

