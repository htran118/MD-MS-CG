# This sample Makefile allows you to make a gsl application
#
# Do not use this Makefile for your projects in this class.  
# Rather, use it as a starting point for writing a Makefile specifically
#  that project.
#
# To use this Makefile, you type:
#
#        make xxxx
#                  
# where
#       xxxx.c or xxxx.cpp is the name of the file you wish to compile 
#       
# A binary named xxxx will be produced
#
# Additional sources for the make can passed as the value of SRCS on 
#  the command line:
#      make SRCS="yyyy.o zzzz.o" xxxx
#
# Proper location of libraries depends on the ENV variable OSTYPE
# To insure this is picked up, you should export OSTYPE prior to running make
#                  

CC = gcc 

CPP = g++

CFLAGS = -ggdb -O3

CPPFLAGS = $(CFLAGS) -std=c++11

# default
ifndef OSTYPE
OSTYPE = linux-gnu
INCLUDES = -L/usr/include
endif
#

# OSX
ifeq ($(findstring darwin,$(OSTYPE)), darwin)
INCLUDES = -I/opt/local/include -I/usr/X11/include 
LDLIBS = -L/opt/local/lib -lgsl -lgslcblas -lm 
DEFS =
endif
#


# FreeBSD  These are out of date
ifeq ($(findstring freebsd,$(OSTYPE)), freebsd)
INCLUDES = -I/usr/local/include
LDLIBS = -L/usr/local/lib -lgsl -lgslcblas -lm 
DEFS =
endif
#

# linux
ifeq ($(findstring linux,$(OSTYPE)), linux)
INCLUDES = 
LDLIBS = -lgsl -lgslcblas -lm
DEFS =
endif
#

## --------------------------------------------------------------  ##

driver/simula.o: driver/simula.cpp libs/tlist.h libs/ptype.h libs/btype.h libs/atype.h libs/dtype.h libs/rtype.h libs/world.h libs/partk.h libs/bound.h libs/angle.h libs/dihed.h libs/resid.h libs/memor.h libs/const.h

driver/memExt.o: driver/memExt.cpp libs/tlist.h libs/ptype.h libs/btype.h libs/atype.h libs/dtype.h libs/rtype.h libs/world.h libs/partk.h libs/bound.h libs/angle.h libs/dihed.h libs/resid.h libs/memor.h libs/const.h

driver/convrt.o: driver/convrt.cpp libs/tlist.h libs/ptype.h libs/btype.h libs/atype.h libs/dtype.h libs/rtype.h libs/world.h libs/partk.h libs/bound.h libs/angle.h libs/dihed.h libs/resid.h libs/memor.h libs/const.h

driver/convrt.o: driver/convrt.cpp libs/tlist.h libs/ptype.h libs/btype.h libs/atype.h libs/dtype.h libs/rtype.h libs/world.h libs/partk.h libs/bound.h libs/angle.h libs/dihed.h libs/resid.h libs/memor.h libs/const.h

#driver/coarse.o: driver/coarse.cpp libs/tlist.h libs/ptype.h libs/btype.h libs/atype.h libs/dtype.h libs/rtype.h libs/world.h libs/partk.h libs/bound.h libs/angle.h libs/dihed.h libs/resid.h libs/memor.h libs/const.h

driver/dtAnal.o: driver/dtAnal.cpp libs/world.h libs/tlist.h libs/const.h

driver/latGen.o: driver/latGen.cpp libs/const.h

test/testCrap.o: test/testCrap.cpp libs/const.h

test/testIntegr.o: test/testIntegr.cpp libs/partk.h libs/const.h

test/testMemory.o: test/testMemory.cpp libs/memor.h libs/partk.h libs/const.h

test/testAuxMem.o: test/testAuxMem.cpp

test/testMZForc.o: test/testMZForc.cpp libs/const.h

test/testMakeCG.o: test/testMakeCG.cpp libs/const.h

test/testFDT.o: test/testFDT.cpp libs/const.h

test/testVecCrs.o: test/testVecCrs.cpp libs/const.h

libs/world.o: libs/world.cpp libs/world.h libs/memor.h libs/resid.h libs/partk.h libs/bound.h libs/angle.h libs/dihed.h libs/tlist.h libs/const.h
libs/resid.o: libs/resid.cpp libs/resid.h libs/partk.h libs/bound.h libs/angle.h libs/dihed.h libs/rtype.h libs/tlist.h libs/const.h
libs/partk.o: libs/partk.cpp libs/partk.h libs/bound.h libs/angle.h libs/dihed.h libs/resid.h libs/ptype.h libs/world.h libs/memor.h libs/const.h
libs/bound.o: libs/bound.cpp libs/bound.h libs/btype.h libs/partk.h libs/world.h
libs/angle.o: libs/angle.cpp libs/angle.h libs/atype.h libs/partk.h libs/bound.h libs/world.h
libs/dihed.o: libs/dihed.cpp libs/dihed.h libs/dtype.h libs/partk.h libs/bound.h libs/angle.h libs/world.h
libs/memor.o: libs/memor.cpp libs/memor.h libs/partk.h libs/const.h
libs/auxMm.o: libs/auxMm.cpp libs/auxMm.h libs/const.h

libs/tlist.o: libs/tlist.cpp libs/tlist.h libs/rtype.h libs/ptype.h libs/btype.h libs/atype.h libs/dtype.h libs/const.h
libs/rtype.o: libs/rtype.cpp libs/rtype.h libs/ptype.h libs/btype.h libs/atype.h libs/dtype.h libs/tlist.h libs/const.h
libs/ptype.o: libs/ptype.cpp libs/ptype.h libs/const.h
libs/btype.o: libs/btype.cpp libs/btype.h libs/ptype.h libs/tlist.h
libs/atype.o: libs/atype.cpp libs/atype.h libs/ptype.h libs/tlist.h
libs/dtype.o: libs/dtype.cpp libs/dtype.h libs/ptype.h libs/tlist.h

libs/const.o: libs/const.cpp libs/const.h

simula: driver/simula.o libs/tlist.o libs/ptype.o libs/btype.o libs/atype.o libs/dtype.o libs/rtype.o libs/world.o libs/partk.o libs/bound.o libs/angle.o libs/dihed.o libs/resid.o libs/memor.o libs/const.o
	$(CPP) $(CPPFLAGS) driver/simula.o $(TYPES) $(OBJECTS) $(SRCS) $(LDLIBS) -o simula

memExt: driver/memExt.o libs/tlist.o libs/ptype.o libs/btype.o libs/atype.o libs/dtype.o libs/rtype.o libs/world.o libs/partk.o libs/bound.o libs/angle.o libs/dihed.o libs/resid.o libs/memor.o libs/const.o
	$(CPP) $(CPPFLAGS) driver/memExt.o $(TYPES) $(OBJECTS) $(SRCS) $(LDLIBS) -o memExt

#coarse: driver/coarse.o libs/tlist.o libs/ptype.o libs/btype.o libs/atype.o libs/dtype.o libs/rtype.o libs/world.o libs/partk.o libs/bound.o libs/angle.o libs/dihed.o libs/resid.o libs/memor.o libs/const.o
#	$(CPP) $(CPPFLAGS) driver/coarse.o $(TYPES) $(OBJECTS) $(SRCS) $(LDLIBS) -o coarse

dtAnal: driver/dtAnal.o libs/tlist.o libs/ptype.o libs/btype.o libs/atype.o libs/dtype.o libs/rtype.o libs/world.o libs/partk.o libs/bound.o libs/angle.o libs/dihed.o libs/resid.o libs/memor.o libs/const.o
	$(CPP) $(CPPFLAGS) driver/dtAnal.o $(TYPES) $(OBJECTS) $(SRCS) $(LDLIBS) -o dtAnal

convrt: driver/convrt.o libs/tlist.o libs/ptype.o libs/btype.o libs/atype.o libs/dtype.o libs/rtype.o libs/world.o libs/partk.o libs/bound.o libs/angle.o libs/dihed.o libs/resid.o libs/const.o
	$(CPP) $(CPPFLAGS) driver/convrt.o $(TYPES) $(OBJECTS) $(SRCS) $(LDLIBS) -o convrt

latGen: driver/latGen.o libs/const.o
	$(CPP) $(CPPFLAGS) driver/latGen.o libs/const.o $(LDLIBS) -o latGen

testCrap: test/testCrap.o libs/const.o
	$(CPP) $(CPPFLAGS) test/testCrap.o libs/const.o $(LDLIBS) -o testCrap

testPBC: test/testPBC.cpp
	$(CPP) $(CPPFLAGS) test/testPBC.cpp -o testPBC

testMod: test/testMod.cpp
	$(CPP) $(CPPFLAGS) test/testMod.cpp -o testMod

testRandom: test/testRandom.cpp
	$(CPP) $(CPPFLAGS) test/testRandom.cpp -o testRandom

testSizeOf: test/testSizeOf.cpp
	$(CPP) $(CPPFLAGS) test/testSizeOf.cpp -o testSizeOf

testIntegr: test/testIntegr.o libs/tlist.o libs/ptype.o libs/btype.o libs/atype.o libs/dtype.o libs/rtype.o libs/world.o libs/partk.o libs/bound.o libs/angle.o libs/dihed.o libs/resid.o libs/memor.o libs/const.o
	$(CPP) $(CPPFLAGS) test/testIntegr.o $(TYPES) $(OBJECTS) $(SRCS) $(LDLIBS) -o testIntegr

testMemory: test/testMemory.o libs/tlist.o libs/ptype.o libs/btype.o libs/atype.o libs/dtype.o libs/rtype.o libs/world.o libs/partk.o libs/bound.o libs/angle.o libs/dihed.o libs/resid.o libs/memor.o libs/const.o
	$(CPP) $(CPPFLAGS) test/testMemory.o $(TYPES) $(OBJECTS) $(SRCS) $(LDLIBS) -o testMemory

testAuxMem: test/testAuxMem.o
	$(CPP) $(CPPFLAGS) test/testAuxMem.o $(LDLIBS) -o testAuxMem

testFDT: test/testFDT.o libs/const.o
	$(CPP) $(CPPFLAGS) test/testFDT.o libs/const.o  $(LDLIBS) -o testFDT

testVecCrs: test/testVecCrs.o libs/const.o
	$(CPP) $(CPPFLAGS) test/testVecCrs.o libs/const.o  $(LDLIBS) -o testVecCrs

testMZForc: test/testMZForc.o libs/const.o
	$(CPP) $(CPPFLAGS) test/testMZForc.o libs/const.o  $(LDLIBS) -o testMZForc

testMakeCG: test/testMakeCG.o libs/const.o
	$(CPP) $(CPPFLAGS) test/testMakeCG.o libs/const.o  $(LDLIBS) -o testMakeCG


## --------------------------------------------------------------  ##

TYPES = libs/ptype.o libs/btype.o libs/atype.o libs/dtype.o libs/rtype.o libs/tlist.o
OBJECTS = libs/partk.o libs/bound.o libs/angle.o libs/dihed.o libs/resid.o libs/memor.o libs/world.o libs/const.o 
LINKS = driver/simula.o driver/dtAnal.o driver/convrt.o driver/latGen.o driver/memExt.o test/testIntegr.o test/testMemory.o test/testFDT.o test/testAuxMem.o test/testVecCrs.o test/testMZForc.o test/testMZForc.o
EXPENDABLES = *~ *stackdump $(LINKS) $(OBJECTS) $(TYPES)
EXECUTABLES = simula latGen memExt dtAnal convrt $(TESTS)
TESTS = testPBC testMod testRandom testSizeOf testIntegr testMemory testFDT testAuxMem testVecCrs testMZForc testMakeCG

.o:
	$(CPP) $@.o $(SRCS) $(LDLIBS) -o $@
.cpp:
	$(CPP) $(CPPFLAGS) $(INCLUDES) $(DEFS) $@.cpp $(SRCS) $(LDLIBS) -o $@

clear:
	rm -f $(EXPENDABLES) $(TESTS)
clean: 
	rm -f $(EXPENDABLES) $(EXECUTABLES)
all:
	rm -f $(EXPENDABLES)
	@$(MAKE) $(EXECUTABLES)
links:
	rm -f $(LINKS)
	@$(MAKE) $(LINKS)
objects:
	rm -f $(OBJECTS)
	@$(MAKE) $(OBJECTS)
types:
	rm -f $(TYPES)
	@$(MAKE) $(TYPES)
tests:
	rm -f $(EXPENDABLES)
	@$(MAKE) $(TESTS)

