#-------------------------------------------------------------------------------
# UMFPACK makefile (for GNU make or original make)
#-------------------------------------------------------------------------------

# UMFPACK requires the AMD package to be in ../AMD

default: library

include ../UFconfig/UFconfig.mk

# compile just the C-callable libraries and demo programs (not mexFunctions)
library:
	( cd ../AMD && $(MAKE) library )
	( cd Lib    && $(MAKE) )
	cat Doc/License

# remove object files, but keep the compiled programs and library archives
clean:
	( cd ../AMD && $(MAKE) clean )
	( cd Lib    && $(MAKE) clean )
	( cd Doc    && $(MAKE) clean )

# clean, and then remove compiled programs and library archives
purge:
	( cd ../AMD && $(MAKE) purge )
	( cd Lib    && $(MAKE) purge )
	( cd Doc    && $(MAKE) purge )

# create PDF documents for the original distribution
doc:
	( cd ../AMD && $(MAKE) doc )
	( cd Doc    && $(MAKE) )

# get ready for distribution
dist: purge
	( cd ../AMD && $(MAKE) dist )
	( cd Doc    && $(MAKE) )

distclean: purge

ccode: library
