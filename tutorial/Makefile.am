noinst_PROGRAMS = singlemode
singlemode_SOURCES = singlemode.cpp 
singlemode_LDADD   = ../model/libMODEL.a  -lm
INCLUDES = -I../ -I../..

EXTRA_DIST = singlemode.cpp plotresults ssa.gp statplot.gp dynplot.gp 

test: singlemode
	./singlemode
	./plotresults

CLEANFILES = *.dat *.*~ *.ps

AUTOMAKE_OPTIONS = foreign
