#------------------------------------------------------------------------------
# AMD Makefile (for GNU Make or original make)
#------------------------------------------------------------------------------

default: library

include ../UFconfig/UFconfig.mk

# compile just the C-callable libraries (not mexFunctions or Demos)
library:
	( cd Lib    ; $(MAKE) )

# remove object files, but keep the compiled programs and library archives
clean:
	( cd Lib    ; $(MAKE) clean )

# clean, and then remove compiled programs and library archives
purge:
	( cd Lib    ; $(MAKE) purge )
	( cd Doc    ; $(MAKE) purge )

distclean: purge

# create PDF documents for the original distribution
doc:
	( cd Doc    ; $(MAKE) )

# get ready for distribution
dist: purge
	( cd Doc    ; $(MAKE) )

ccode: library

lib: library
