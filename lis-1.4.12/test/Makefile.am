LDADD = -L$(top_builddir)/src -llis@LIBSUFFIX@
AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir)/include
AM_FFLAGS = -I$(top_builddir) -I$(top_srcdir)/include
AM_LDFLAGS = $(AMDEFS)

EXTRA_DIST = defs test.sh testmat.mtx testmat2.mtx

testdir = ${prefix}/share/examples/lis

test_DATA = testmat.mtx testmat2.mtx

test_SCRIPTS = defs test.sh

test_PROGRAMS = test1 test2 test3 test4 test5 test6 etest1 etest2 etest3 etest4 etest5 etest6 spmvtest1 spmvtest2 spmvtest3 spmvtest4 spmvtest5
if ENABLE_FORTRAN
  test_PROGRAMS += test1f test4f etest1f etest4f
endif

test1_SOURCES  = test1.c
test2_SOURCES  = test2.c
test3_SOURCES  = test3.c
test4_SOURCES  = test4.c
test5_SOURCES  = test5.c
test6_SOURCES  = test6.c
etest1_SOURCES  = etest1.c
etest2_SOURCES  = etest2.c
etest3_SOURCES  = etest3.c
etest4_SOURCES  = etest4.c
etest5_SOURCES  = etest5.c
etest6_SOURCES  = etest6.c
spmvtest1_SOURCES  = spmvtest1.c
spmvtest2_SOURCES  = spmvtest2.c
spmvtest3_SOURCES  = spmvtest3.c
spmvtest4_SOURCES  = spmvtest4.c
spmvtest5_SOURCES  = spmvtest5.c
if ENABLE_FORTRAN
  test1f_SOURCES  = test1f.F
  test4f_SOURCES  = test4f.F
  etest1f_SOURCES  = etest1f.F
  etest4f_SOURCES  = etest4f.F
endif

bin_PROGRAMS = lsolve esolve spmvtest1 spmvtest2 spmvtest3 spmvtest4 spmvtest5

lsolve_SOURCES  = test1.c
esolve_SOURCES  = etest5.c

PPF77COMPILE = $(F77) $(AMDEFS) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FFLAGS) $(FFLAGS)
LTPPF77COMPILE = $(LIBTOOL) --mode=compile --tag=F77 $(F77) $(AMDEFS) $(DEFS) \
	$(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
	$(AM_FFLAGS) $(FFLAGS)
COMPILE = $(CC) $(AMDEFS) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) \
	$(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS)
LTCOMPILE = $(LIBTOOL) --mode=compile --tag=CC $(CC) $(AMDEFS) $(DEFS) \
	$(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
	$(AM_CFLAGS) $(CFLAGS)

if ENABLE_SAAMG
  LINK = $(LIBTOOL) --mode=link $(FC) \
	$(AM_LDFLAGS) $(LDFLAGS) $(FCLDFLAGS) -o $@ $(LIBS) $(CLIBS)
  F77LINK = $(LIBTOOL) --mode=link $(FC) \
	$(AM_LDFLAGS) $(LDFLAGS) -o $@ $(LIBS) $(CLIBS)
else 
if ENABLE_FORTRAN
  LINK = $(LIBTOOL) --mode=link $(FC) \
	$(AM_LDFLAGS) $(LDFLAGS) $(FCLDFLAGS) -o $@ $(LIBS) $(CLIBS)
  F77LINK = $(LIBTOOL) --mode=link $(F77) \
	$(AM_LDFLAGS) $(LDFLAGS) -o $@ $(LIBS) $(CLIBS)
else
  LINK = $(LIBTOOL) --mode=link $(CC) $(AM_LDFLAGS) $(LDFLAGS) -o $@ $(LIBS)
endif
endif

CLEANFILES = *.il


TESTS_ENVIRONMENT = top_builddir=$(top_builddir) testdir=$(testdir) \
                    enable_mpi=$(enable_mpi) \
                    enable_omp=$(enable_omp) enable_saamg=$(enable_saamg) \
                    enable_fortran=$(enable_fortran) enable_quad=$(enable_quad) \
                    MPIRUN=$(MPIRUN) MPINP=$(MPINP)

TESTS = test.sh

installcheck: 
	@if test ! -f $(testdir)/test.sh; then \
	    echo '=== Test programs not installed'; \
	    exit; \
	fi; \
	$(TESTS_ENVIRONMENT) $(testdir)/test.sh

