
PRIMME: PReconditioned Iterative MultiMethod Eigensolver
========================================================

PRIMME, pronounced as *prime*, is a high-performance library for computing a few eigenvalues/eigenvectors, and singular values/vectors.
PRIMME is especially optimized for large, difficult problems.
Real symmetric and complex Hermitian problems, standard :math:`A x = \lambda x` and generalized :math:`A x = \lambda B x`, are supported.
Besides, standard eigenvalue problems with a normal matrix are supported.
It can find largest, smallest, or interior singular/eigenvalues, and can use preconditioning to accelerate convergence. 
PRIMME is written in C99, but complete interfaces are provided for Fortran, MATLAB, Python, and R.

Making and Linking
------------------

To generate the static and the shared library type::

    make lib     #  builds lib/libprimme.a
    make solib   #  builds lib/libprimme.so (or lib/libprimme.dylib)

The shared library is generated with the action `solib` instead. Usual flags are supported

* `CC`, compiler program such as ``gcc``, ``clang`` or ``icc``
* `CFLAGS`, compiler options such as ``-g`` or ``-O3``
* `CUDADIR`, directory of CUDA installation (optional)
* `MAGMADIR`, directory of MAGMA_ installation (optional)
* `PRIMME_WITH_HALF`, activates support for half precision if it set to `yes`;
  compiler supporting `__fp16` is required, e.g., clang.

The flags can be indicated by customizing `Make_flags` or directly introduced at the command line::

    make lib CC=clang CFLAGS='-O3'

For building the external interfaces just do::

    make matlab       # Set MATLAB=/path/Matlab/bin/matlab MEX=/path/Matlab/bin/mex if needed
    make matlab-cuda  # Requires to set CUDADIR and MAGMADIR
    make octave
    make python
    make R_install

We provide packages of the released version for R (see `R PRIMME`_)::

    install.packages('PRIMME')

and Python (see `Python primme`_)::

    pip install numpy   # if numpy is not installed yet
    pip install scipy   # if scipy is not installed yet
    pip install future  # if using python 2
    conda install mkl-devel # if using Anaconda Python distribution
    pip install primme

C Library Interface
-------------------

To compute few eigenvalues and eigenvectors from a real symmetric matrix call::

    int dprimme(double *evals, double *evecs, double *resNorms, 
                primme_params *primme);

The call arguments are:

* `evals`, array to return the found eigenvalues;
* `evecs`, array to return the found eigenvectors;
* `resNorms`, array to return the residual norms of the found eigenpairs; and
* `primme`, structure that specify the matrix problem, which eigenvalues are wanted and several method options.

To compute few singular values and vectors from a matrix call::

    int dprimme_svds(double *svals, double *svecs, double *resNorms, 
                primme_svds_params *primme_svds);

The call arguments are:

* `svals`, array to return the found singular values;
* `svecs`, array to return the found vectors;
* `resNorms`, array to return the residual norms of the triplets; and
* `primme_svds`, structure that specify the matrix problem, which values are wanted and several method options.

There are available versions for `half` and `float` and complex variants.
See documentation in `readme.txt` file and in ``doc`` directory; also it is online at doc_.
The `examples` directory has several self-contained examples in C, C++ and F77, some of them using PETSc_ and MAGMA_.

Citing this code 
----------------

Please cite (bibtex_):

* A. Stathopoulos and J. R. McCombs *PRIMME: PReconditioned Iterative
  MultiMethod Eigensolver: Methods and software description*, ACM
  Transaction on Mathematical Software Vol. 37, No. 2, (2010),
  21:1-21:30.

* L. Wu, E. Romero and A. Stathopoulos, *PRIMME_SVDS: A High-Performance
  Preconditioned SVD Solver for Accurate Large-Scale Computations*,
  SIAM J. Sci. Comput., Vol. 39, No. 5, (2017), S248--S271.

More information on the algorithms and research that led to this
software can be found in the rest of the papers. The work has been
supported by a number of grants from the National Science Foundation.

* A. Stathopoulos, *Nearly optimal preconditioned methods for Hermitian
  eigenproblems under limited memory. Part I: Seeking one eigenvalue*, SIAM
  J. Sci. Comput., Vol. 29, No. 2, (2007), 481--514.

* A. Stathopoulos and J. R. McCombs, *Nearly optimal preconditioned
  methods for Hermitian eigenproblems under limited memory. Part II:
  Seeking many eigenvalues*, SIAM J. Sci. Comput., Vol. 29, No. 5, (2007),
  2162-2188.

* J. R. McCombs and A. Stathopoulos, *Iterative Validation of
  Eigensolvers: A Scheme for Improving the Reliability of Hermitian
  Eigenvalue Solvers*, SIAM J. Sci. Comput., Vol. 28, No. 6, (2006),
  2337-2358.

* A. Stathopoulos, *Locking issues for finding a large number of eigenvectors
  of Hermitian matrices*, Tech Report: WM-CS-2005-03, July, 2005.

* L. Wu and A. Stathopoulos, *A Preconditioned Hybrid SVD Method for Computing
  Accurately Singular Triplets of Large Matrices*, SIAM J. Sci. Comput. 37-5(2015),
  pp. S365-S388.

License Information
-------------------

PRIMME is licensed under the 3-clause license BSD.
Python and Matlab interfaces have BSD-compatible licenses.
Source code under `tests` is compatible with LGPLv3.
Details can be taken from COPYING.txt.

Contact Information 
-------------------

For reporting bugs or questions about functionality contact `Andreas Stathopoulos`_ by
email, `andreas` at `cs.wm.edu`. See further information in
the webpage http://www.cs.wm.edu/~andreas/software.

Support
-------

- National Science Foundation through grants CCF 1218349, ACI SI2-SSE 1440700, and NSCI 1835821
- Department of Energy through grant Exascale Computing Project 17-SC-20-SC

.. _`Andreas Stathopoulos`: http://www.cs.wm.edu/~andreas/software
.. _`github`: https://github.com/primme/primme
.. _`doc`: http://www.cs.wm.edu/~andreas/software/doc/readme.html
.. _`R PRIMME`: https://cran.r-project.org/web/packages/PRIMME/index.html
.. _`Python primme`: https://pypi.org/project/primme/
.. _PETSc : http://www.mcs.anl.gov/petsc/
.. _`bibtex`: https://raw.githubusercontent.com/primme/primme/master/doc/primme.bib
.. _MAGMA: http://icl.cs.utk.edu/magma/
