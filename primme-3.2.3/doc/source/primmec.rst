
.. role:: ccode(code) 
   :language: c

.. highlight:: c

C Library Interface
-------------------

The PRIMME interface is composed of the following functions.
To solve real symmetric and complex Hermitian problems, standard :math:`A x = \lambda x` and generalized :math:`A x = \lambda B x`, call:

.. only:: not text

   .. parsed-literal::

      int :c:func:`dprimme <dprimme>` (double \*evals, double \*evecs, double \*resNorms,
                              primme_params \*primme)
      int :c:func:`zprimme <zprimme>` (double \*evals, PRIMME_COMPLEX_DOUBLE \*evecs,
                       double \*resNorms, primme_params \*primme)

.. only:: text

   ::

      int dprimme(double *evals, double *evecs, double *resNorms, 
                  primme_params *primme);

      int zprimme(double *evals, PRIMME_COMPLEX_DOUBLE *evecs, double *resNorms, 
                  primme_params *primme);

There are more versions for matrix problems working in other precisions:

+-----------+--------------------+--------------------+
| Precision |        Real        |       Complex      |
+===========+====================+====================+
| half      | :c:func:`hprimme`  | :c:func:`kprimme`  |
|           | :c:func:`hsprimme` | :c:func:`ksprimme` |
+-----------+--------------------+--------------------+
| single    | :c:func:`sprimme`  | :c:func:`cprimme`  |
+-----------+--------------------+--------------------+
| double    | :c:func:`dprimme`  | :c:func:`zprimme`  |
+-----------+--------------------+--------------------+

To solve standard eigenproblems with normal but not necessarily Hermitian matrices call:

.. only:: not text

   .. parsed-literal::

      int :c:func:`zprimme_normal <zprimme_normal>` (PRIMME_COMPLEX_DOUBLE \*evals,
                       PRIMME_COMPLEX_DOUBLE \*evecs,
                       double \*resNorms, primme_params \*primme)

.. only:: text

   ::

      int zprimme_normal(PRIMME_COMPLEX_DOUBLE *evals, PRIMME_COMPLEX_DOUBLE *evecs,
                       double *resNorms, primme_params *primme);

There are more versions for matrix problems working in other precisions:

+-----------+---------------------------+
| Precision |         Complex           |
+===========+===========================+
| half      | :c:func:`kprimme_normal`  |
|           | :c:func:`kcprimme_normal` |
+-----------+---------------------------+
| single    | :c:func:`cprimme_normal`  |
+-----------+---------------------------+
| double    | :c:func:`zprimme_normal`  |
+-----------+---------------------------+

Other useful functions:

.. only:: not text

   .. parsed-literal::

      void :c:func:`primme_initialize <primme_initialize>` (primme_params \*primme)
      int :c:func:`primme_set_method <primme_set_method>` (primme_preset_method method,
                                                           primme_params \*params)
      void :c:func:`primme_display_params <primme_display_params>` (primme_params primme)
      void :c:func:`primme_free <primme_Free>` (primme_params \*primme)

.. only:: text

   ::

      void primme_initialize(primme_params *primme);
      int primme_set_method(primme_preset_method method,
                                           primme_params *params);
      void primme_display_params(primme_params primme);
      void primme_free(primme_params primme);

PRIMME stores its data on the structure :c:type:`primme_params`.
See :ref:`guide-params` for an introduction about its fields.


Running
^^^^^^^

To use PRIMME, follow these basic steps.

#. Include::

      #include "primme.h"   /* header file is required to run primme */

#. Initialize a PRIMME parameters structure for default settings:

   .. only:: not text
   
       .. parsed-literal::

          :c:type:`primme_params` primme;
          :c:func:`primme_initialize <primme_initialize>` (&primme);

   .. only:: text
   
      ::
   
         primme_params primme;
         
         primme_initialize(&primme);
   
#. Set problem parameters (see also :ref:`guide-params`), and,
   optionally, set one of the :c:type:`preset methods <primme_preset_method>`:

   .. only:: not text

      .. parsed-literal::

         primme.\ |matrixMatvec| = LaplacianMatrixMatvec; /\* MV product \*/
         primme.\ |n| = 100;                   /\* set problem dimension \*/
         primme.\ |numEvals| = 10;       /\* Number of wanted eigenpairs \*/
         ret = :c:func:`primme_set_method <primme_set_method>` (method, &primme);
         ...

   .. only:: text

      ::

         primme.matrixMatvec = LaplacianMatrixMatvec; /* MV product */
         primme.n = 100;                   /* set problem dimension */
         primme.numEvals = 10;       /* Number of wanted eigenpairs */
         ret = primme_set_method(method, &primme);
         ...

#. Then call the solver:

   .. only:: not text
  
      .. parsed-literal::
 
         ret = :c:func:`dprimme <dprimme>` (evals, evecs, resNorms, &primme);
   
   .. only:: text
   
      ::
   
         ret = dprimme(evals, evecs, resNorms, &primme);

   The call arguments are:

   * `evals`, array to return the found eigenvalues;
   * `evecs`, array to return the found eigenvectors;
   * `resNorms`, array to return the residual norms of the found eigenpairs; and
   * `ret`, returned error code.

#. To free the work arrays in PRIMME:

   .. only:: not text
  
      .. parsed-literal::
 
         :c:func:`primme_free <primme_Free>` (&primme);
   
   .. only:: text
   
      ::
   
         primme_free(&primme);

See usage examples at the directory `examples`.

.. _guide-params:

Parameters Guide
^^^^^^^^^^^^^^^^

PRIMME stores the data on the structure :c:type:`primme_params`, which has the next fields:
   
.. only:: not text

      | *Basic*
      | ``PRIMME_INT`` |n|,  matrix dimension.
      | ``void (*`` |matrixMatvec| ``)(...)``, matrix-vector product.
      | ``void (*`` |massMatrixMatvec| ``)(...)``, mass matrix-vector product (null for standard problems).
      | ``int`` |numEvals|, how many eigenpairs to find.
      | ``primme_target`` |target|, which eigenvalues to find.
      | ``int`` |numTargetShifts|, for targeting interior eigenpairs.
      | ``double *`` |targetShifts|
      | ``double`` |eps|, tolerance of the residual norm of converged eigenpairs.
      |
      | *For parallel programs*
      | ``int`` |numProcs|, number of processes
      | ``int`` |procID|,  rank of this process
      | ``PRIMME_INT`` |nLocal|,  number of rows stored in this process
      | ``void (*`` |globalSumReal| ``)(...)``, sum reduction among processes
      | ``void (*`` |broadcastReal| ``)(...)``, broadcast array among processes
      |
      | *Accelerate the convergence*
      | ``void (*`` |applyPreconditioner| ``)(...)``, preconditioner-vector product.
      | ``int`` |initSize|, initial vectors as approximate solutions.
      | ``int`` |maxBasisSize|
      | ``int`` |minRestartSize|
      | ``int`` |maxBlockSize|
      |
      | *User data*
      | ``void *`` |commInfo|
      | ``void *`` |matrix|
      | ``void *`` |massMatrix|
      | ``void *`` |preconditioner|
      | ``void *`` |convtest|
      | ``void *`` |monitor|
      |
      | *Advanced options*
      | ``PRIMME_INT`` |ldevecs|, leading dimension of the evecs.
      | ``int`` |numOrthoConst|, orthogonal constrains to the eigenvectors.
      | ``int`` |dynamicMethodSwitch|
      | ``int`` |locking|
      | ``PRIMME_INT`` |maxMatvecs|
      | ``PRIMME_INT`` |maxOuterIterations|
      | ``PRIMME_INT`` |iseed| ``[4]``
      | ``double`` |aNorm|
      | ``double`` |BNorm|
      | ``double`` |invBNorm|
      | ``int`` |printLevel|
      | ``FILE *`` |outputFile|
      | ``double *`` |ShiftsForPreconditioner|
      | ``primme_init`` |initBasisMode|
      | ``struct projection_params`` :c:member:`projectionParams <primme_params.projectionParams.projection>`
      | ``struct restarting_params`` :c:member:`restartingParams <primme_params.restartingParams.scheme>`
      | ``struct correction_params`` :c:member:`correctionParams <primme_params.correctionParams.precondition>`
      | ``struct primme_stats`` :c:member:`stats <primme_params.stats.numOuterIterations>`
      | ``void (*`` |convTestFun| ``)(...)``, custom convergence criterion.
      | ``PRIMME_INT`` |ldOPS|, leading dimension to use in |matrixMatvec|.
      | ``void (*`` |monitorFun| ``)(...)``, custom convergence history.
      | ``primme_op_datatype`` |matrixMatvec_type|
      | ``primme_op_datatype`` |massMatrixMatvec_type|
      | ``primme_op_datatype`` |applyPreconditioner_type|
      | ``primme_op_datatype`` |globalSumReal_type|
      | ``primme_op_datatype`` |broadcastReal_type|
      | ``primme_op_datatype`` |internalPrecision|
      | ``primme_orth`` |orth|

.. only:: text

   ::

      /* Basic */
      PRIMME_INT n;                               // matrix dimension
      void (*matrixMatvec)(...);             // matrix-vector product
      void (*massMatrixMatvec)(...);    // mass matrix-vector product
      int numEvals;                    // how many eigenpairs to find
      primme_target target;              // which eigenvalues to find
      int numTargetShifts;       // for targeting interior eigenpairs
      double *targetShifts;
      double eps;            // tolerance of the converged eigenpairs
      
      /* For parallel programs */
      int numProcs;           // number of processes
      int procID;             // rank of this process 
      PRIMME_INT nLocal;      // number of rows stored in this process
      void (*globalSumReal)(...); // sum reduction among processes
      void (*broadcastReal)(...); // broadcast array among processes
      
      /* Accelerate the convergence */
      void (*applyPreconditioner)(...);     // precond-vector product
      int initSize;       // initial vectors as approximate solutions
      int maxBasisSize;
      int minRestartSize;
      int maxBlockSize;
      
      /* User data */
      void *commInfo;
      void *matrix;
      void *massMatrix;
      void *preconditioner;
      void *convtest;
      void *monitor;
      
      /* Advanced options */
      PRIMME_INT ldevecs; // leading dimension of the evecs
      int numOrthoConst; // orthogonal constrains to the eigenvectors
      int dynamicMethodSwitch;
      int locking;
      PRIMME_INT maxMatvecs;
      PRIMME_INT maxOuterIterations;
      PRIMME_INT iseed[4];
      double aNorm;
      double BNorm;
      double invBNorm;
      int printLevel;
      FILE *outputFile;
      double *ShiftsForPreconditioner;
      primme_init initBasisMode;
      struct projection_params projectionParams;
      struct restarting_params restartingParams;
      struct correction_params correctionParams;
      struct primme_stats stats;
      void (*convTestFun)(...); // custom convergence criterion
      PRIMME_INT ldOPS;   // leading dimension to use in matrixMatvec
      void (*monitorFun)(...); // custom convergence history
      primme_op_datatype matrixMatvec_type;
      primme_op_datatype massMatrixMatvec_type;
      primme_op_datatype applyPreconditioner_type;
      primme_op_datatype globalSumReal_type;
      primme_op_datatype broadcastReal_type;
      primme_op_datatype internalPrecision;
      primme_orth orth;
 
PRIMME requires the user to set at least the dimension of the matrix (|n|) and
the matrix-vector product (|matrixMatvec|), as they define the problem to be solved.
For parallel programs, |nLocal|, |procID| and |globalSumReal| are also required.

In addition, most users would want to specify how many eigenpairs to find,
|numEvals|, and provide a preconditioner |applyPreconditioner| (if available).

It is useful to have set all these before calling :c:func:`primme_set_method`.
Also, if users have a preference on |maxBasisSize|, |maxBlockSize|, etc, they
should also provide them into :c:type:`primme_params` prior to the
:c:func:`primme_set_method` call. This helps :c:func:`primme_set_method` make
the right choice on other parameters. It is sometimes useful to check the actual
parameters that PRIMME is going to use (before calling it) or used (on return)
by printing them with :c:func:`primme_display_params`.

Interface Description
^^^^^^^^^^^^^^^^^^^^^

The next enumerations and functions are declared in ``primme.h``.

?primme
"""""""

.. c:function:: int hprimme(PRIMME_HALF *evals, PRIMME_HALF *evecs, PRIMME_HALF *resNorms, primme_params *primme)
.. c:function:: int hsprimme(float *evals, PRIMME_HALF *evecs, float *resNorms, primme_params *primme)
.. c:function:: int kprimme(PRIMME_HALF *evals, PRIMME_COMPLEX_HALF *evecs, PRIMME_HALF *resNorms, primme_params *primme)
.. c:function:: int ksprimme(float *evals, PRIMME_COMPLEX_HALF *evecs, float *resNorms, primme_params *primme)

   .. versionadded:: 3.0
.. c:function:: int sprimme(float *evals, float *evecs, float *resNorms, primme_params *primme)
.. c:function:: int cprimme(float *evals, PRIMME_COMPLEX_FLOAT *evecs, float *resNorms, primme_params *primme)

   .. versionadded:: 2.0
.. c:function:: int dprimme(double *evals, double *evecs, double *resNorms, primme_params *primme)
.. c:function:: int zprimme(double *evals, PRIMME_COMPLEX_DOUBLE *evecs, double *resNorms, primme_params *primme)

   Solve a real symmetric/Hermitian standard or generalized eigenproblem.

   All arrays should be hosted on CPU. The computations are performed on CPU (see :c:func:`magma_dprimme` for using GPUs).

   :param evals: array at least of size |numEvals| to store the
      computed eigenvalues; all processes in a parallel run return this local array with the same values.

   :param evecs: array at least of size |nLocal| times (|numOrthoConst| + |numEvals|) with leading dimension |ldevecs|
      to store column-wise the (local part for this process of the) computed eigenvectors.

   :param resNorms: array at least of size |numEvals| to store the
      residual norms of the computed eigenpairs; all processes in parallel run return this local array with
      the same values.

   :param primme: parameters structure.

   :return: error indicator; see :ref:`error-codes`.

   On input, ``evecs`` should start with the content of the |numOrthoConst| vectors,
   followed by the |initSize| vectors.
 
   On return, the i-th eigenvector starts at evecs[( |numOrthoConst| + i)\* |ldevecs| ], with value `evals[i]` and associated residual 2-norm `resNorms[i]`.
   The first vector has index i=0. The number of eigenpairs marked as converged (see |eps|) is returned on |initSize|. Since version 4.0, if the returned error code is `PRIMME_MAIN_ITER_FAILURE`, PRIMME may return also unconverged eigenpairs and its residual norms in `evecs`, `evals`, and `resNorms` starting at i=|initSize| and going up to either |numEvals|-1 or the last `resNorms[i]` with non-negative value.
 
   All internal operations are performed at the same precision than ``evecs`` unless the user sets |internalPrecision| otherwise.
   The functions :c:func:`hsprimme` and :c:func:`ksprimme` perform all computations in half precision by default and report the eigenvalues and the residual norms in single precision. These functions may help in applications that may be not built with a compiler supporting half precision.

   The type and precision of the callbacks is also the same as `evecs`. Although this can be changed. See details for |matrixMatvec|, |massMatrixMatvec|, |applyPreconditioner|, |globalSumReal|, |broadcastReal|, and |convTestFun|.

`cublas_?primme` & `magma_?primme` 
""""""""""""""""""""""""""""""""""

.. c:function:: int cublas_hprimme(PRIMME_HALF *evals, PRIMME_HALF *evecs, PRIMME_HALF *resNorms, primme_params *primme)
.. c:function:: int cublas_hsprimme(float *evals, PRIMME_HALF *evecs, float *resNorms, primme_params *primme)
.. c:function:: int cublas_kprimme(PRIMME_HALF *evals, PRIMME_COMPLEX_HALF *evecs, PRIMME_HALF *resNorms, primme_params *primme)
.. c:function:: int cublas_sprimme(float *evals, float *evecs, float *resNorms, primme_params *primme)
.. c:function:: int cublas_ksprimme(float *evals, PRIMME_COMPLEX_HALF *evecs, float *resNorms, primme_params *primme)
.. c:function:: int cublas_cprimme(float *evals, PRIMME_COMPLEX_FLOAT *evecs, float *resNorms, primme_params *primme)
.. c:function:: int cublas_dprimme(double *evals, double *evecs, double *resNorms, primme_params *primme)
.. c:function:: int cublas_zprimme(double *evals, PRIMME_COMPLEX_DOUBLE *evecs, double *resNorms, primme_params *primme)
.. c:function:: int magma_hprimme(PRIMME_HALF *evals, PRIMME_HALF *evecs, PRIMME_HALF *resNorms, primme_params *primme)
.. c:function:: int magma_hsprimme(float *evals, PRIMME_HALF *evecs, float *resNorms, primme_params *primme)
.. c:function:: int magma_kprimme(PRIMME_HALF *evals, PRIMME_COMPLEX_HALF *evecs, PRIMME_HALF *resNorms, primme_params *primme)
.. c:function:: int magma_sprimme(float *evals, float *evecs, float *resNorms, primme_params *primme)
.. c:function:: int magma_ksprimme(float *evals, PRIMME_COMPLEX_HALF *evecs, float *resNorms, primme_params *primme)
.. c:function:: int magma_cprimme(float *evals, PRIMME_COMPLEX_FLOAT *evecs, float *resNorms, primme_params *primme)
.. c:function:: int magma_dprimme(double *evals, double *evecs, double *resNorms, primme_params *primme)
.. c:function:: int magma_zprimme(double *evals, PRIMME_COMPLEX_DOUBLE *evecs, double *resNorms, primme_params *primme)

   Solve a real symmetric/Hermitian standard or generalized eigenproblem.

   Most of the computations are performed on GPU (see :c:func:`dprimme` for using only the CPU).

   :param evals: CPU array at least of size |numEvals| to store the
      computed eigenvalues; all processes in a parallel run return this local array with the same values.

   :param evecs: GPU array at least of size |nLocal| times (|numOrthoConst| + |numEvals|) with leading dimension |ldevecs|
      to store column-wise the (local part for this process of the) computed eigenvectors.

   :param resNorms: CPU array at least of size |numEvals| to store the
      residual norms of the computed eigenpairs; all processes in parallel run return this local array with
      the same values.

   :param primme: parameters structure.

   :return: error indicator; see :ref:`error-codes`.

   On input, ``evecs`` should start with the content of the |numOrthoConst| vectors,
   followed by the |initSize| vectors.
 
   On return, the i-th eigenvector starts at evecs[( |numOrthoConst| + i)\* |ldevecs| ], with value `evals[i]` and associated residual 2-norm `resNorms[i]`.
   The first vector has index i=0. The number of eigenpairs marked as converged (see |eps|) is returned on |initSize|. Since version 4.0, if the returned error code is `PRIMME_MAIN_ITER_FAILURE`, PRIMME may return also unconverged eigenpairs and its residual norms in `evecs`, `evals`, and `resNorms` starting at i=|initSize| and going up to either |numEvals|-1 or the last `resNorms[i]` with non-negative value.
 
   All internal operations are performed at the same precision than ``evecs`` unless the user sets |internalPrecision| otherwise.
   The functions :c:func:`hsprimme` and :c:func:`ksprimme` perform all computations in half precision by default and report the eigenvalues and the residual norms in single precision. These functions may help in applications that may be not built with a compiler supporting half precision.

   The type and precision of the callbacks is also the same as `evecs`. Although this can be changed. See details for |matrixMatvec|, |massMatrixMatvec|, |applyPreconditioner|, |globalSumReal|, |broadcastReal|, and |convTestFun|.

   .. versionadded:: 3.0

?primme_normal
""""""""""""""

.. c:function:: int kprimme_normal(PRIMME_COMPLEX_HALF *evals, PRIMME_COMPLEX_HALF *evecs, PRIMME_HALF *resNorms, primme_params *primme)
.. c:function:: int kcprimme_normal(PRIMME_COMPLEX_FLOAT *evals, PRIMME_COMPLEX_HALF *evecs, float *resNorms, primme_params *primme)
.. c:function:: int cprimme_normal(PRIMME_COMPLEX_FLOAT *evals, PRIMME_COMPLEX_FLOAT *evecs, float *resNorms, primme_params *primme)
.. c:function:: int zprimme_normal(PRIMME_COMPLEX_DOUBLE *evals, PRIMME_COMPLEX_DOUBLE *evecs, double *resNorms, primme_params *primme)

   Solve a normal standard eigenproblem, which may not be Hermitian.

   All arrays should be hosted on CPU. The computations are performed on CPU (see :c:func:`magma_zprimme_normal` for using GPUs).

   :param evals: array at least of size |numEvals| to store the
      computed eigenvalues; all processes in a parallel run return this local array with the same values.

   :param evecs: array at least of size |nLocal| times (|numOrthoConst| + |numEvals|) with leading dimension |ldevecs|
      to store column-wise the (local part for this process of the) computed eigenvectors.

   :param resNorms: array at least of size |numEvals| to store the
      residual norms of the computed eigenpairs; all processes in parallel run return this local array with
      the same values.

   :param primme: parameters structure.

   :return: error indicator; see :ref:`error-codes`.

   On input, ``evecs`` should start with the content of the |numOrthoConst| vectors,
   followed by the |initSize| vectors.
 
   On return, the i-th eigenvector starts at evecs[( |numOrthoConst| + i)\* |ldevecs| ], with value `evals[i]` and associated residual 2-norm `resNorms[i]`.
   The first vector has index i=0. The number of eigenpairs marked as converged (see |eps|) is returned on |initSize|. Since version 4.0, if the returned error code is `PRIMME_MAIN_ITER_FAILURE`, PRIMME may return also unconverged eigenpairs and its residual norms in `evecs`, `evals`, and `resNorms` starting at i=|initSize| and going up to either |numEvals|-1 or the last `resNorms[i]` with non-negative value.
 
   All internal operations are performed at the same precision than ``evecs`` unless the user sets |internalPrecision| otherwise.
   The functions :c:func:`hsprimme` and :c:func:`ksprimme` perform all computations in half precision by default and report the eigenvalues and the residual norms in single precision. These functions may help in applications that may be not built with a compiler supporting half precision.

   The type and precision of the callbacks is also the same as `evecs`. Although this can be changed. See details for |matrixMatvec|, |applyPreconditioner|, |globalSumReal|, |broadcastReal|, and |convTestFun|.

   .. versionadded:: 3.0

`cublas_?primme_normal` & `magma_?primme_normal`
""""""""""""""""""""""""""""""""""""""""""""""""

.. c:function:: int magma_kprimme_normal(PRIMME_COMPLEX_HALF *evals, PRIMME_COMPLEX_HALF *evecs, PRIMME_HALF *resNorms, primme_params *primme)
.. c:function:: int magma_kcprimme_normal(PRIMME_COMPLEX_FLOAT *evals, PRIMME_COMPLEX_HALF *evecs, float *resNorms, primme_params *primme)
.. c:function:: int magma_cprimme_normal(PRIMME_COMPLEX_FLOAT *evals, PRIMME_COMPLEX_FLOAT *evecs, float *resNorms, primme_params *primme)
.. c:function:: int magma_zprimme_normal(PRIMME_COMPLEX_DOUBLE *evals, PRIMME_COMPLEX_DOUBLE *evecs, double *resNorms, primme_params *primme)

   Solve a normal standard eigenproblem, which may not be Hermitian.

   Most of the computations are performed on GPU (see :c:func:`zprimme_normal` for using only the CPU).

   :param evals: CPU array at least of size |numEvals| to store the
      computed eigenvalues; all processes in a parallel run return this local array with the same values.

   :param evecs: GPU array at least of size |nLocal| times (|numOrthoConst| + |numEvals|) with leading dimension |ldevecs|
      to store column-wise the (local part for this process of the) computed eigenvectors.

   :param resNorms: CPU array at least of size |numEvals| to store the
      residual norms of the computed eigenpairs; all processes in parallel run return this local array with
      the same values.

   :param primme: parameters structure.

   :return: error indicator; see :ref:`error-codes`.

   On input, ``evecs`` should start with the content of the |numOrthoConst| vectors,
   followed by the |initSize| vectors.
 
   On return, the i-th eigenvector starts at evecs[( |numOrthoConst| + i)\* |ldevecs| ], with value `evals[i]` and associated residual 2-norm `resNorms[i]`.
   The first vector has index i=0. The number of eigenpairs marked as converged (see |eps|) is returned on |initSize|. Since version 4.0, if the returned error code is `PRIMME_MAIN_ITER_FAILURE`, PRIMME may return also unconverged eigenpairs and its residual norms in `evecs`, `evals`, and `resNorms` starting at i=|initSize| and going up to either |numEvals|-1 or the last `resNorms[i]` with non-negative value.
 
   All internal operations are performed at the same precision than ``evecs`` unless the user sets |internalPrecision| otherwise.
   The functions :c:func:`hsprimme` and :c:func:`ksprimme` perform all computations in half precision by default and report the eigenvalues and the residual norms in single precision. These functions may help in applications that may be not built with a compiler supporting half precision.

   The type and precision of the callbacks is also the same as `evecs`. Although this can be changed. See details for |matrixMatvec|, |applyPreconditioner|, |globalSumReal|, |broadcastReal|, and |convTestFun|.

   .. versionadded:: 3.0

primme_initialize
"""""""""""""""""

.. c:function:: void primme_initialize(primme_params *primme)

   Initialize a PRIMME parameters structure to the default values.

   After calling :c:func:`dprimme` (or a variant), call :c:func:`primme_free` to release allocated resources by PRIMME.

   :param primme: parameters structure.

   Example::

      primme_params primme;
      primme_initialize(&primme);

      primme.n = 100;
      ...
      dprimme(evals, evecs, rnorms, &primme);
      ...

      primme_free(&primme);
 
   See the alternative function :c:func:`primme_params_create` that also allocates the `primme_params` structure.

primme_params_create
""""""""""""""""""""

.. c:function:: primme_params* primme_params_create(void)

   Allocate and initialize a parameters structure to the default values.

   After calling :c:func:`dprimme` (or a variant), call :c:func:`primme_params_destroy` to release allocated resources by PRIMME.

   :return: pointer to a parameters structure.

   Example::

      primme_params *primme = primme_params_create();

      primme->n = 100;
      ...
      dprimme(evals, evecs, rnorms, primme);
      ...

      primme_params_destroy(primme);

   See the alternative function :c:func:`primme_initialize` that only initializes the structure.

   .. versionadded:: 3.0

primme_set_method
"""""""""""""""""

.. c:function:: int primme_set_method(primme_preset_method method, primme_params *primme)

   Set PRIMME parameters to one of the preset configurations.

   :param method: preset configuration; one of

      | |DYNAMIC|
      | |DEFAULT_MIN_TIME|
      | |DEFAULT_MIN_MATVECS|
      | |Arnoldi|
      | |GD|
      | |GD_plusK|
      | |GD_Olsen_plusK|
      | |JD_Olsen_plusK|
      | |RQI|
      | |JDQR|
      | |JDQMR|
      | |JDQMR_ETol|
      | |STEEPEST_DESCENT|
      | |LOBPCG_OrthoBasis|
      | |LOBPCG_OrthoBasis_Window|

   :param primme: parameters structure.

   See also :ref:`methods`.

primme_display_params
"""""""""""""""""""""

.. c:function:: void primme_display_params(primme_params primme)

   Display all printable settings of ``primme`` into the file descriptor |outputFile|.

   :param primme: parameters structure.

primme_free
"""""""""""

.. c:function:: void primme_free(primme_params *primme)

   Free memory allocated by PRIMME.

   :param primme: parameters structure.

primme_params_destroy
"""""""""""""""""""""

.. c:function:: int primme_params_destroy(primme_params *primme)

   Free memory allocated by PRIMME associated to a parameters structure created
   with :c:func:`primme_params_create`.

   :param primme: parameters structure.

   :return: nonzero value if the call is not successful.

   .. versionadded:: 3.0

.. include:: epilog.inc
