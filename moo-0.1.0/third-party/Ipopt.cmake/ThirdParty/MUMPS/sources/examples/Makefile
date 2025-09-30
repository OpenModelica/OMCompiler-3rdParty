#
#  This file is part of MUMPS 5.8.1, released
#  on Wed Jul 30 16:49:18 UTC 2025
#
topdir = ..
libdir = $(topdir)/lib

default: d

.PHONY: default all s d c z multi clean
.SECONDEXPANSION:

all:	c z s d multi

c:	csimpletest csimpletest_save_restore
z:	zsimpletest zsimpletest_save_restore
s:	ssimpletest ssimpletest_save_restore
d:	dsimpletest dsimpletest_save_restore c_example_save_restore c_example
multi:	multiple_arithmetics_example


include $(topdir)/Makefile.inc



LIBSMUMPS = -L$(libdir) -lsmumps$(PLAT) -lmumps_common$(PLAT)
LIBCMUMPS = -L$(libdir) -lcmumps$(PLAT) -lmumps_common$(PLAT)
LIBSDMUMPS = -L$(libdir) -ldmumps$(PLAT) -lsmumps$(PLAT) -lmumps_common$(PLAT)
LIBCZMUMPS = -L$(libdir) -lzmumps$(PLAT) -lcmumps$(PLAT) -lmumps_common$(PLAT)

ssimpletest:   $$@.o
	$(FL) -o $@ $(OPTL) ssimpletest.o  $(LIBSMUMPS) $(LORDERINGS) $(LIBS) $(RPATH_OPT) $(LIBBLAS) $(LIBOTHERS)

dsimpletest:  $$@.o 
	$(FL) -o $@ $(OPTL) dsimpletest.o  $(LIBSDMUMPS) $(LORDERINGS) $(LIBS) $(RPATH_OPT) $(LIBBLAS) $(LIBOTHERS)

csimpletest:  $$@.o
	$(FL) -o $@ $(OPTL) csimpletest.o  $(LIBCMUMPS) $(LORDERINGS) $(LIBS) $(RPATH_OPT) $(LIBBLAS) $(LIBOTHERS)

zsimpletest:  $$@.o
	$(FL) -o $@ $(OPTL) zsimpletest.o  $(LIBCZMUMPS) $(LORDERINGS) $(LIBS) $(RPATH_OPT) $(LIBBLAS) $(LIBOTHERS)

c_example: $$@.o
	$(FL) -o $@ $(OPTL) $@.o $(LIBSDMUMPS) $(LORDERINGS) $(LIBS) $(RPATH_OPT) $(LIBBLAS) $(LIBOTHERS)


multiple_arithmetics_example: $$@.o
	$(FL) -o $@ $(OPTL) $@.o $(LIBSDMUMPS) $(LIBCZMUMPS) $(LORDERINGS) $(LIBS) $(RPATH_OPT) $(LIBBLAS) $(LIBOTHERS)

ssimpletest_save_restore:   $$@.o
	$(FL) -o $@ $(OPTL) ssimpletest_save_restore.o  $(LIBSMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

dsimpletest_save_restore:  $$@.o 
	$(FL) -o $@ $(OPTL) dsimpletest_save_restore.o  $(LIBSDMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

csimpletest_save_restore:  $$@.o
	$(FL) -o $@ $(OPTL) csimpletest_save_restore.o  $(LIBCMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

zsimpletest_save_restore:  $$@.o
	$(FL) -o $@ $(OPTL) zsimpletest_save_restore.o  $(LIBCZMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

c_example_save_restore: $$@.o
	$(FL) -o $@ $(OPTL) $@.o $(LIBSDMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

.SUFFIXES: .c .F .o
.F.o:
	$(FC) $(OPTF) -I. -I$(topdir)/include -I$(topdir)/src $(INCS) -c $*.F $(OUTF)$*.o
.c.o:
	$(CC) $(OPTC) $(CDEFS) -I. -I$(topdir)/include -I$(topdir)/src $(INCS) -c $*.c $(OUTC)$*.o

clean:
	$(RM) *.o [sdcz]simpletest c_example multiple_arithmetics_example ssimpletest_save_restore dsimpletest_save_restore csimpletest_save_restore zsimpletest_save_restore c_example_save_restore
