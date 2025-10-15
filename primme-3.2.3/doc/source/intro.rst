.. highlight:: bash

PRIMME: PReconditioned Iterative MultiMethod Eigensolver
--------------------------------------------------------

PRIMME, pronounced as *prime*, is a high-performance library for computing a few eigenvalues/eigenvectors, and singular values/vectors.
PRIMME is especially optimized for large, difficult problems.
Real symmetric and complex Hermitian problems, standard :math:`A x = \lambda x` and generalized :math:`A x = \lambda B x`, are supported.
Besides standard eigenvalue problems with a normal matrix are supported.
It can find largest, smallest, or interior singular/eigenvalues, and can use preconditioning to accelerate convergence. 
PRIMME is written in C99, but complete interfaces are provided for Fortran, MATLAB, Python, and R.

Incompatibilities
^^^^^^^^^^^^^^^^^

From PRIMME 2.2 to 3.0:

* Removed constants ``primme_thick`` and ``primme_dtr``, and the member ``scheme`` from ``restarting_params``.

* Added members |numBroadcast|, |volumeBroadcast|, ``flopsDense``, |timeBroadcast|, ``timeDense``, ``estimateBNorm``, ``estimateInvBNorm``, and |lockingIssue| to ``primme_stats``.

* Added members |matrixMatvec_type|, |applyPreconditioner_type|, |massMatrixMatvec_type|, |globalSumReal_type|, |broadcastReal|, |broadcastReal_type|, |BNorm|, |invBNorm|, |orth|, |internalPrecision|, |massMatrix|, |convTestFun_type|, |monitorFun_type|, |queue|, and ``profile`` to :c:type:`primme_params`.

* Added members |SmatrixMatvec_type|, |SapplyPreconditioner_type|, |SglobalSumReal_type|, |SbroadcastReal|, |SbroadcastReal_type|, |SinternalPrecision|, |SconvTestFun_type|, |SmonitorFun_type|, |Squeue|, and ``profile`` to :c:type:`primme_svds_params`.

* Changed callbacks |monitorFun| and |SmonitorFun|.

* Changed the value of all constants; see :c:func:`primme_get_member_f77`, :c:func:`primme_set_member_f77`, :c:func:`primme_svds_get_member_f77`, and :c:func:`primme_svds_set_member_f77`.

* Removed ``intWorkSize``, ``realWorkSize``, ``intWork``, ``realWork`` from :c:type:`primme_params` and :c:type:`primme_svds_params`.

From PRIMME 2.0 to 2.1:

* Added members |monitorFun| and ``monitor`` to :c:type:`primme_params`.

* Added members |SmonitorFun| and ``monitor`` to :c:type:`primme_svds_params`.

* Renamed ``PRIMME_SUBSPACE_ITERATION`` as |STEEPEST_DESCENT|.

From PRIMME 1.x to 2.0:

* Prototype of callbacks has changed: |matrixMatvec|, |applyPreconditioner|, |massMatrixMatvec| and |globalSumReal|.

* The next parameters are :c:type:`PRIMME_INT`: |n|, |nLocal|, |maxMatvecs|, |iseed|, |numOuterIterations|, |numRestarts|, |numMatvecs| and |numMatvecs|; use the macro ``PRIMME_INT_P`` to print the values.

* Rename the values of the enum :c:type:`primme_preset_method`.

* Rename ``primme_Free`` to :c:func:`primme_free`.

* Integer parameters in Fortran functions are of the same size as :c:type:`PRIMME_INT`, which is ``integer*8`` by default.

* Extra parameter in many Fortran functions to return the error code.

* Removed ``primme_display_stats_f77``.

Changelog
^^^^^^^^^

Changes in PRIMME 3.2.3 (released on Nov 8, 2024):

* Fixed bug in updating the convergence flags when blockSize>1.

* Added a new GPU interface that does not need MAGMA, supporting CUDA & ROCM
  (see :c:func:`cublas_dprimme`, :c:func:`cublas_zprimme_normal`, and :c:func:`cublas_dprimme_svds`).

* Fixed compilation issues with several GPU libraries: MAGMA, cu/hipBLAS, cu/hipSPARSE.

* Updated the Python package to support Python 3.6+.

Changes in PRIMME 3.2 (released on Jan 29, 2021):

* Fixed Intel 2021 compiler error ``"Unsupported combination of types for <tgmath.h>."``

* Fixed compiling issues with PGI compiler also about ``tgmath.h``.

* Fixed :c:func:`dprimme` and other variants not returning error code `PRIMME_MAIN_ITER_FAILURE` when it should do in some corner cases.

* Fixed warnings from gcc/clang undefined behavior sanitizers.

* Matlab: renamed `disp` to `reportLevel`.

* Matlab: add flag `returnUnconverged` to return unconverged pairs optionally.

* Matlab: return primme_params/primme_svds_params. 

Changes in PRIMME 3.1 (released on May 2, 2020):

* Fixed compilation issues in F90 interface and examples.

* Fixed bug in block orthogonalization.

* Updated Python interface to Python version 3.8.

Changes in PRIMME 3.0 (released on December 14, 2019):

* Added support for the generalized Hermitian eigenvalue problem (see |massMatrixMatvec|) and the standard normal eigenvalue problem (see :c:func:`zprimme_normal`).

* Added support for GPU (see :c:func:`magma_dprimme`, :c:func:`magma_zprimme_normal`, and :c:func:`magma_dprimme_svds`).

* Added support for half precision (see :c:func:`hprimme` and :c:func:`kprimme`, and other variants for normal eigenproblems and singular value problems).

* Added block orthogonalization (see |orth|).

* Resolution of all linear system of equations simultaneously in Jacobi-Davidson.

* Added interface for Fortran 90.

* Added an optional callback for broadcasting (see |broadcastReal| and |SbroadcastReal|).

* The callbacks can work with different precision than the main call (see for instance |matrixMatvec_type| and |globalSumReal_type|).

* Added new counters: |numGlobalSum|, |volumeGlobalSum|, |numBroadcast|, |volumeGlobalSum|, |timeOrtho|, |timeGlobalSum|, |timeBroadcast|.

* Added :c:func:`primme_params_create`, :c:func:`primme_params_destroy`, :c:func:`primme_svds_params_create`, and :c:func:`primme_svds_params_destroy`.

Changes in PRIMME 2.2 (released on October 26, 2018):

* Improved stability for single precision.

* Improved support for the shared library.

* Updated PETSc examples in Fortran; added new examples in single precision.

* Improved support for FreeBSD and MacOS.

* New install and uninstall actions.

* MATLAB interface support for user-defined stopping criterion (see |convTestFun| and |SconvTestFun|) and stopping with ctr+c.

* Optional suffix on BLAS_/LAPACK_ function names (for OpenBLAS, see :ref:`PRIMME_BLAS_SUFFIX <making>`).

* Replaced XHEGV by XHEGVX (to support ESSL).

* Fixed bugs in the library and in the Matlab interface.

Changes in PRIMME 2.1 (released on April 4, 2017):

* Improved robustness by broadcasting the result of critical LAPACK_
  operations instead of replicating them on every process; this is
  useful when using a threaded BLAS_/LAPACK_ or when some parallel
  processes may run on different architectures or libraries.

* New stopping criteria in QMR that improve performance for interior
  problems.

* MATLAB interface reimplementation with support for singular value
  problems, :mat:func:`primme_svds()`, with double and single precision, and 
  compatible with Octave.

* R interface

* Proper reporting of convergence history for singular value solvers.

Changes in PRIMME 2.0 (released on September 19, 2016):

* Changed license to BSD 3-clause.

* New support for singular value problems; see :c:func:`dprimme_svds`.

* New support for ``float`` and ``complex float`` arithmetic.

* Support for problem dimensions larger than 2^31, without requiring
  BLAS_ and LAPACK_ compiled with 64-bits integers.

* Improve robustness and performance for interior problems; implemented advanced refined
  and harmonic-Ritz extractions.

* Python interface compatible with NumPy_ and `SciPy Library`_.

* Added parameter to indicate the leading dimension of the input/output matrices and to
  return an error code in callbacks |matrixMatvec|, |applyPreconditioner|,
  |massMatrixMatvec| and |globalSumReal|.

* Changed to type :c:type:`PRIMME_INT` the options |n|, |nLocal|, |maxMatvecs|
  and |iseed|, and the stats counters |numOuterIterations|, |numRestarts|, |numMatvecs|,
  |numPreconds|. Also changed |realWorkSize| to ``size_t``. Fortran interface functions
  will expect an ``interger`` of size compatible with :c:type:`PRIMME_INT` for
  all parameters with integer type: ``int``, :c:type:`PRIMME_INT` and ``size_t``;
  see also parameter ``value`` in functions
  :c:func:`primmetop_set_member_f77`,
  :c:func:`primmetop_get_member_f77`,
  :c:func:`primme_set_member_f77` and
  :c:func:`primme_get_member_f77`.

* Added parameter to return an error code in Fortran interface functions:
  :c:func:`primmetop_set_member_f77`,
  :c:func:`primmetop_get_member_f77`,
  :c:func:`primme_set_member_f77` and
  :c:func:`primme_get_member_f77`.

* Added leading dimension for ``evecs`` |ldevecs| and preferred leading dimension
  for the operators |ldOPs|, such as |matrixMatvec|.

* Optional user-defined convergence function, |convTestFun|.

* Prefixed methods with ``PRIMME_``. Rename Fortran constants from ``PRIMMEF77_``
  to ``PRIMME_``.

* Removed ``primme_display_stats_f77``.

Changes in PRIMME 1.2.2 (released on October 13, 2015):

* Fixed wrong symbols in :file:`libdprimme.a` and :file:`libzprimme.a`.

* :c:func:`primme_set_method` sets |JDQMR| instead of |JDQMR_ETol| for preset methods
  |DEFAULT_MIN_TIME| and |DYNAMIC| when seeking interior values.

* Fixed compilation of driver with a PETSc_ installation without HYPRE.

* Included the content of the environment variable ``INCLUDE`` for compiling the driver.


Changes in PRIMME 1.2.1 (released on September 7, 2015):

* Added MATLAB interface to full PRIMME functionality.

* Support for BLAS_/LAPACK_ with 64bits integers (``-DPRIMME_BLASINT_SIZE=64``).

* Simplified configuration of Make_flags and Make_links (removed ``TOP`` variable 
  and replaced defines ``NUM_SUM`` and ``NUM_IBM`` by ``F77UNDERSCORE``).

* Replaced directories :file:`DTEST` and :file:`ZTEST` by :file:`TEST`, that has:

  * :file:`driver.c`: read matrices in MatrixMarket format and PETSc_ binary and
    call PRIMME with the parameters specified in a file; support
    complex arithmetic and MPI and can use PETSc_ preconditioners.
  * :file:`ex*.c` and :file:`ex*.f`: small, didactic examples of usage in C and Fortran
    and in parallel (with PETSc_).

* Fixed a few minor bugs and improved documentation (especially the F77 interface).

* Using Sphinx_ to manage documentation.

Changes in PRIMME 1.2 (released on December 21, 2014):

* A Fortran compiler is no longer required for building the PRIMME library.
  Fortran programs can still be linked to PRIMME's F77 interface.

* Fixed some uncommon issues with the F77 interface.

* PRIMME can be called now multiple times from the same program.

* Performance improvements in the QMR inner solver, especially for 
  complex arithmetic.

* Fixed a couple of bugs with the locking functionality. 

  * In certain extreme cases where all eigenvalues of a matrix were needed.
  * The order of selecting interior eigenvalues.

  The above fixes have improved robustness and performance. 

* PRIMME now assigns unique random seeds per parallel process 
  for up to 4096^3  (140 trillion) processes.

* For the |DYNAMIC| method, fixed issues with initialization and 
  synchronization decisions across multiple processes.

* Fixed uncommon library interface bugs, coordinated better
  setting the method and the user setting of parameters, and improved 
  the interface in the sample programs and makefiles.

* Other performance and documentation improvements.

License Information
^^^^^^^^^^^^^^^^^^^

PRIMME is licensed under the 3-clause license BSD.
Python and MATLAB interfaces have BSD-compatible licenses.
Source code under :file:`tests` is compatible with LGPLv3.
Details can be taken from :file:`COPYING.txt`::

   Copyright (c) 2018, College of William & Mary
   All rights reserved.


Citing the code 
^^^^^^^^^^^^^^^ 

.. only:: latex

   Please cite [r1]_ and [r6]_. Find the BibTeX in the following and also in :file:`doc/primme.bib`:

   .. literalinclude:: ../primme.bib
      :language: bibtex

.. only:: text

   Please cite (find the BibTeX in :file:`doc/primme.doc`):

.. only:: not latex and not text

   Please cite (`BibTeX`_):

.. [r1] A. Stathopoulos and J. R. McCombs PRIMME: *PReconditioned Iterative
   MultiMethod Eigensolver: Methods and software description*, ACM
   Transaction on Mathematical Software Vol. 37, No. 2, (2010),
   21:1-21:30.

.. [r6] L. Wu, E. Romero and A. Stathopoulos, *PRIMME_SVDS: A High-Performance
  Preconditioned SVD Solver for Accurate Large-Scale Computations*,
  J. Sci. Comput., Vol. 39, No. 5, (2017), S248--S271.

.. only:: latex

   More information on the algorithms and research that led to this
   software can be found in the rest of the papers [r2]_, [r3]_, [r4]_, [r5]_, [r7]_.
   The work has been supported by a number of grants from the
   National Science Foundation.

.. only:: not latex

   More information on the algorithms and research that led to this
   software can be found in the rest of the papers. The work has been
   supported by a number of grants from the National Science Foundation.

.. [r2] A. Stathopoulos, *Nearly optimal preconditioned methods for Hermitian
   eigenproblems under limited memory. Part I: Seeking one eigenvalue*, SIAM
   J. Sci. Comput., Vol. 29, No. 2, (2007), 481--514.

.. [r3] A. Stathopoulos and J. R. McCombs, *Nearly optimal preconditioned
   methods for Hermitian eigenproblems under limited memory. Part II:
   Seeking many eigenvalues*, SIAM J. Sci. Comput., Vol. 29, No. 5, (2007),
   2162-2188.

.. [r4] J. R. McCombs and A. Stathopoulos, *Iterative Validation of
   Eigensolvers: A Scheme for Improving the Reliability of Hermitian
   Eigenvalue Solvers*, SIAM J. Sci. Comput., Vol. 28, No. 6, (2006),
   2337-2358.

.. [r5] A. Stathopoulos, *Locking issues for finding a large number of eigenvectors
   of Hermitian matrices*, Tech Report: WM-CS-2005-03, July, 2005.

.. [r7] L. Wu and A. Stathopoulos, *A Preconditioned Hybrid SVD Method for Computing
  Accurately Singular Triplets of Large Matrices*, SIAM J. Sci. Comput. 37-5(2015),
  pp. S365-S388.


Contact Information 
^^^^^^^^^^^^^^^^^^^

For reporting bugs or questions about functionality contact `Andreas Stathopoulos`_ by
email, `andreas` at `cs.wm.edu`. See further information in
the webpage http://www.cs.wm.edu/~andreas/software and on github_.

Support
^^^^^^^

- National Science Foundation through grants CCF 1218349, ACI SI2-SSE 1440700, and NSCI 1835821
- Department of Energy through grant Exascale Computing Project 17-SC-20-SC



Directory Structure
^^^^^^^^^^^^^^^^^^^

The next directories and files should be available:

* :file:`COPYING.txt`, license;
* :file:`Make_flags`,  flags to be used by makefiles to compile library and tests;
* :file:`Link_flags`,  flags needed in making and linking the test programs;
* :file:`include/`,    directory with headers files;
* :file:`src/`,        directory with the source code for :file:`libprimme`:

   * :file:`include/`,   common headers;
   * :file:`eigs/`,      eigenvalue interface and implementation;
   * :file:`svds/`,      singular value interface and implementation;
   * :file:`tools/`,     tools used to generated some headers;

* :file:`Matlab/`,       MATLAB interface;
* :file:`Python/`,       Python interface;
* :file:`examples/`,     sample programs in C, C++ and F77, both sequential and parallel;
* :file:`tests/`,        drivers for testing purpose and test cases;
* :file:`lib/libprimme.a`,   the PRIMME library (to be made);
* :file:`makefile`       main make file;
* :file:`readme.txt`     text version of the documentation;
* :file:`doc/`           directory with the HTML and PDF versions of the documentation.

.. _making :

Making and Linking
^^^^^^^^^^^^^^^^^^

:file:`Make_flags` has the flags and compilers used to make :file:`libprimme.a`:

* `CC`, compiler program such as ``gcc``, ``clang`` or ``icc``.
* `CFLAGS`, compiler options such as ``-g`` or ``-O3`` and macro definitions
   like the ones described next.

Compiler flags for the BLAS_ and LAPACK_ libraries:

* ``-DF77UNDERSCORE``, if Fortran appends an underscore to function names
  (usually it does).
* ``-DPRIMME_BLASINT_SIZE=64``, if the library integers are 64-bit integer (``kind=8``) type,
  aka ILP64 interface; usually integers are 32-bits even in 64-bit architectures (aka LP64 interface).
* ``-DPRIMME_BLAS_SUFFIX=<suffix>``, set a suffix to BLAS/LAPACK function names; for instance,
  OpenBlas compiled with ILP64 may append ``_64`` to the function names.

By default PRIMME sets the integer type for matrix dimensions and counters (:c:type:`PRIMME_INT`)
to 64 bits integer ``int64_t``. This can be changed by setting the macro ``PRIMME_INT_SIZE``
to one of the following values:

- ``0``: use the regular ``int`` of your compiler.
- ``32``: use C99 ``int32_t``. 
- ``64``: use C99 ``int64_t``. 

.. note::

   When ``-DPRIMME_BLASINT_SIZE=64`` is set the code uses the type ``int64_t``
   supported by the C99 standard. In case the compiler doesn't honor the
   standard, you can set the corresponding type name supported, for instance
   ``-DPRIMME_BLASINT_SIZE=__int64``.

After customizing :file:`Make_flags`, type this to generate :file:`libprimme.a`::

    make lib

Making can be also done at the command line::

    make lib CC=clang CFLAGS='-O3'

:file:`Link_flags` has the flags for linking with external libraries and making the executables
located in :file:`examples` and :file:`tests`:

* `LDFLAGS`, linker flags such as ``-framework Accelerate``.
* `LIBS`, flags to link with libraries (BLAS_ and LAPACK_ are required), such as ``-lprimme -llapack -lblas -lgfortran -lm``.

After that, type this to compile and execute a simple test::

    $ make test
    ...
    Test passed!
    ...
    Test passed! 

In case of linking problems check flags in `LDFLAGS` and `LIBS` and consider
to add/remove ``-DF77UNDERSCORE`` from `CFLAGS`. If the execution fails consider
to add/remove ``-DPRIMME_BLASINT_SIZE=64`` from `CFLAGS`.

Full description of actions that `make` can take:

* `make lib`, builds the static library :file:`libprimme.a`.
* `make solib`, builds the shared library :file:`libprimme.so`.
* `make install`, installs header files and the static and dynamic libraries.
* `make uninstall`, uninstalls header files and the static and dynamic libraries.
* `make matlab`, builds `libprimme.a` compatible with MATLAB and the MATLAB module.
* `make octave`, builds `libprimme.a` and the Octave module.
* `make python`, builds `libprimme.a` and the Python module.
* `make python_install`, install the Python module.
* `make R_install`, builds and installs the R package.
* `make test`, build and execute simple examples.
* `make clean`, removes all :file:`*.o`, :file:`a.out`, and core files from :file:`src`.

Considerations using an IDE
"""""""""""""""""""""""""""

PRIMME can be built in other environments such as Anjuta, Eclipse, KDevelop, Qt Creator,
Visual Studio and XCode. To build the PRIMME library do the following:

#. Create a new project and include the source files under the directory :file:`src`.
#. Add the directories :file:`include` and :file:`src/include` as include directories.

To build an example code using PRIMME make sure:

- to add a reference for PRIMME, BLAS_ and LAPACK_ libraries;
- to add the directory :file:`include` as an include directory.

Main Contributors
^^^^^^^^^^^^^^^^^

* James R. McCombs
* Eloy Romero Alcalde
* Andreas Stathopoulos
* Lingfei Wu

.. _`github`: https://github.com/primme/primme
.. _`BibTeX`: https://raw.githubusercontent.com/primme/primme/master/doc/primme.bib

.. include:: epilog.inc
