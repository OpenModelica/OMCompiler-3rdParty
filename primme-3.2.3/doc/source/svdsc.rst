
.. role:: ccode(code) 
   :language: c

.. highlight:: c

C Library Interface
-------------------

.. versionadded:: 2.0

The PRIMME SVDS interface is composed of the following functions.
To solve real and complex singular value problems call respectively:

.. only:: not text

   .. parsed-literal::

      int :c:func:`dprimme_svds <dprimme_svds>` (double \*svals, double \*svecs, double \*resNorms,
                             primme_svds_params \*primme_svds)
      int :c:func:`zprimme_svds <zprimme_svds>` (double \*svals, :c:type:`PRIMME_COMPLEX_DOUBLE` \*svecs,
                       double \*resNorms, primme_svds_params \*primme_svds)

.. only:: text

   ::

      int dprimme_svds(double *svals, double *svecs, double *resNorms, 
                  primme_svds_params *primme);

      int zprimme_svds(double *svals, PRIMME_COMPLEX_DOUBLE *svecs, double *resNorms, 
                  primme_svds_params *primme);

There are more versions for matrix problems working in other precisions:

+-----------+-------------------------+-------------------------+
| Precision |        Real             |       Complex           |
+===========+=========================+=========================+
| half      | :c:func:`hprimme_svds`  | :c:func:`kprimme_svds`  |
|           | :c:func:`hsprimme_svds` | :c:func:`ksprimme_svds` |
+-----------+-------------------------+-------------------------+
| single    | :c:func:`sprimme_svds`  | :c:func:`cprimme_svds`  |
+-----------+-------------------------+-------------------------+
| double    | :c:func:`dprimme_svds`  | :c:func:`zprimme_svds`  |
+-----------+-------------------------+-------------------------+

Other useful functions:

.. only:: not text

   .. parsed-literal::

      void :c:func:`primme_svds_initialize <primme_svds_initialize>` (primme_svds_params \*primme_svds)
      int :c:func:`primme_svds_set_method <primme_svds_set_method>` (primme_svds_preset_method method,
         primme_preset_method methodStage1,
         primme_preset_method methodStage2, primme_svds_params \*primme_svds)
      void :c:func:`primme_svds_display_params <primme_svds_display_params>` (primme_svds_params primme_svds)
      void :c:func:`primme_svds_free <primme_svds_free>` (primme_svds_params \*primme_svds)

.. only:: text

   ::

      void primme_svds_initialize(primme_svds_params *primme_svds);
      int primme_svds_set_method(primme_svds_preset_method method,
            primme_preset_method methodStage1, primme_preset_method methodStage2,
            primme_svds_params *primme_svds);
      void primme_svds_display_params(primme_svds_params primme_svds);
      void primme_svds_Free(primme_svds_params *primme_svds);

PRIMME SVDS stores its data on the structure :c:type:`primme_svds_params`.
See :ref:`svds-guide-params` for an introduction about its fields.


Running
^^^^^^^

To use PRIMME SVDS, follow these basic steps.

#. Include::

      #include "primme.h"   /* header file is required to run primme */

#. Initialize a PRIMME SVDS parameters structure for default settings:

   .. only:: not text
   
       .. parsed-literal::

          :c:type:`primme_svds_params` primme_svds;
          :c:func:`primme_svds_initialize <primme_svds_initialize>` (&primme_svds);

   .. only:: text
   
      ::
   
         primme_svds_params primme_svds;
         
         primme_svds_initialize(&primme_svds);
   
#. Set problem parameters (see also :ref:`svds-guide-params`), and,
   optionally, set one of the :c:type:`preset methods <primme_svds_preset_method>`:

   .. only:: not text

      .. parsed-literal::

         primme_svds.\ |SmatrixMatvec| = matrixMatvec; /\* MV product \*/
         primme_svds.\ |Sm| = 1000;     /\* set the matrix dimensions \*/
         primme_svds.\ |Sn| = 100;
         primme_svds.\ |SnumSvals| = 10; /\* Number of singular values \*/
         :c:func:`primmesvds_set_method <primme_svds_set_method>` (primme_svds_hybrid, PRIMME_DEFAULT_METHOD,
                                  PRIMME_DEFAULT_METHOD, &primme_svds);
         ...

   .. only:: text

      ::

         primme_svds.matrixMatvec = matrixMatvec; /* MV product */
         primme_svds.m = 1000;                    /* set problem dimension */
         primme_svds.n = 100;
         primme_svds.numSvals = 10;    /* Number of wanted singular values */
         primme_svds_set_method(primme_svds_hybrid, PRIMME_DEFAULT_METHOD,
                                   PRIMME_DEFAULT_METHOD, &primme_svds);
         ...

#. Then to solve a real singular value problem call:

   .. only:: not text
  
      .. parsed-literal::
 
         ret = :c:func:`dprimme_svds <dprimme_svds>` (svals, svecs, resNorms, &primme_svds);
   
   .. only:: text
   
      ::
   
         ret = dprimme_svds(svals, svecs, resNorms, &primme_svds);

   The previous is the double precision call. There is available calls for complex
   double, single and complex single; check :c:func:`zprimme_svds`, :c:func:`sprimme_svds`
   and :c:func:`cprimme_svds`.

   To solve complex singular value problems call:

   .. only:: not text
   
      .. parsed-literal::

         ret = :c:func:`zprimme_svds <zprimme_svds>` (svals, svecs, resNorms, &primme_svds);
   
   .. only:: text
   
      ::
   
         ret = zprimme_svds(svals, svecs, resNorms, &primme_svds);

   The call arguments are:

   * `svals`, array to return the found singular values;
   * `svecs`, array to return the found left and right singular vectors;
   * `resNorms`, array to return the residual norms of the found triplets; and
   * `ret`, returned error code.

#. To free the work arrays in PRIMME SVDS:

   .. only:: not text
  
      .. parsed-literal::
 
         :c:func:`primme_svds_free <primme_svds_Free>` (&primme_svds);
   
   .. only:: text
   
      ::
   
         primme_svds_free(&primme_svds);

.. _svds-guide-params:

Parameters Guide
^^^^^^^^^^^^^^^^

PRIMME SVDS stores the data on the structure :c:type:`primme_svds_params`, which has the next fields:
   
.. only:: not text

      | *Basic*
      | ``PRIMME_INT`` |Sm|,  number of rows of the matrix.
      | ``PRIMME_INT`` |Sn|,  number of columns of the matrix.
      | ``void (*`` |SmatrixMatvec| ``)(...)``, matrix-vector product.
      | ``int`` |SnumSvals|, how many singular triplets to find.
      | ``primme_svds_target`` |Starget|, which singular values to find.
      | ``double`` |Seps|, tolerance of the residual norm of converged triplets.
      |
      | *For parallel programs*
      | ``int`` |SnumProcs|, number of processes
      | ``int`` |SprocID|,  rank of this process
      | ``PRIMME_INT`` |SmLocal|, number of rows stored in this process
      | ``PRIMME_INT`` |SnLocal|, number of columns stored in this process
      | ``void (*`` |SglobalSumReal| ``)(...)``, sum reduction among processes
      |
      | *Accelerate the convergence*
      | ``void (*`` |SapplyPreconditioner| ``)(...)``, preconditioner-vector product.
      | ``int`` |SinitSize|, initial vectors as approximate solutions.
      | ``int`` |SmaxBasisSize|
      | ``int`` |SminRestartSize|
      | ``int`` |SmaxBlockSize|
      |
      | *User data*
      | ``void *`` |ScommInfo|
      | ``void *`` |Smatrix|
      | ``void *`` |Spreconditioner|
      | ``void *`` |Sconvtest|
      | ``void *`` |Smonitor|
      |
      | *Advanced options*
      | ``int`` |SnumTargetShifts|, for targeting interior singular values.
      | ``double *`` |StargetShifts|
      | ``int`` |SnumOrthoConst|, orthogonal constrains to the singular vectors.
      | ``int`` |Slocking|
      | ``PRIMME_INT`` |SmaxMatvecs|
      | ``PRIMME_INT`` |Siseed| ``[4]``
      | ``double`` |SaNorm|
      | ``int`` |SprintLevel|
      | ``FILE *`` |SoutputFile|
      | ``primme_svds_operator`` |Smethod|
      | ``primme_svds_operator`` |SmethodStage2|
      | |primme_params| |Sprimme|
      | |primme_params| |SprimmeStage2|
      | ``void (*`` |SconvTestFun| ``)(...)``, custom convergence criterion.
      | ``void (*`` |SmonitorFun| ``)(...)``, custom convergence history.
      | ``primme_op_datatype`` |SmatrixMatvec_type|
      | ``primme_op_datatype`` |SapplyPreconditioner_type|
      | ``primme_op_datatype`` |SglobalSumReal_type|
      | ``primme_op_datatype`` |SbroadcastReal_type|
      | ``primme_op_datatype`` |SinternalPrecision|

.. only:: text

   ::

      /* Basic */
      PRIMME_INT m;                    // number of rows of the matrix
      PRIMME_INT n;                 // number of columns of the matrix
      void (*matrixMatvec)(...);              // matrix-vector product
      int numSvals;              // how many singular triplets to find
      primme_svds_target target;      // which singular values to find
      double eps;               // tolerance of the converged triplets
      
      /* For parallel programs */
      int numProcs;          // number of processes
      int procID;            // rank of this process
      PRIMME_INT mLocal;     // number of rows stored in this process
      PRIMME_INT nLocal;     // number of columns stored in this process
      void (*globalSumReal)(...); // sum reduction among processes
      
      /* Accelerate the convergence */
      void (*applyPreconditioner)(...); // preconditioner-vector product
      int initSize;        // initial vectors as approximate solutions
      int maxBasisSize;
      int minRestartSize;
      int maxBlockSize;
      
      /* User data */
      void *commInfo;
      void *matrix;
      void *preconditioner;
      void *convtest;
      void *monitor;
      
      /* Advanced options */
      int numTargetShifts;        // for targeting interior values
      double *targetShifts;
      int numOrthoConst;   // orthogonal constrains to the vectors
      int locking;
      PRIMME_INT maxMatvecs;
      PRIMME_INT iseed[4];
      double aNorm;
      int printLevel;
      FILE * outputFile;
      primme_svds_operator method;
      primme_svds_operator methodStage2;
      primme_params primme;
      primme_params primmeStage2;
      void (*convTestFun)(...); // custom convergence criterion
      void (*monitorFun)(...); // custom convergence history
      primme_op_datatype matrixMatvec_type;
      primme_op_datatype applyPreconditioner_type;
      primme_op_datatype globalSumReal_type;
      primme_op_datatype broadcastReal_type;
      primme_op_datatype internalPrecision;


PRIMME SVDS requires the user to set at least the matrix dimensions (|Sm| x |Sn|) and
the matrix-vector product (|SmatrixMatvec|), as they define the problem to be solved.
For parallel programs, |SmLocal|, |SnLocal|, |SprocID| and |SglobalSumReal| are also required.

In addition, most users would want to specify how many singular triplets to find,
and provide a preconditioner (if available).

It is useful to have set all these before calling :c:func:`primme_svds_set_method`.
Also, if users have a preference on |SmaxBasisSize|, |SmaxBlockSize|, etc, they
should also provide them into :c:type:`primme_svds_params` prior to the
:c:func:`primme_svds_set_method` call. This helps :c:func:`primme_svds_set_method` make
the right choice on other parameters. It is sometimes useful to check the actual
parameters that PRIMME SVDS is going to use (before calling it) or used (on return)
by printing them with :c:func:`primme_svds_display_params`.

Interface Description
^^^^^^^^^^^^^^^^^^^^^

The next enumerations and functions are declared in ``primme.h``.

?primme_svds
""""""""""""

.. c:function:: int hprimme_svds(PRIMME_HALF *svals, PRIMME_HALF *svecs, PRIMME_HALF *resNorms, primme_svds_params *primme_svds)
.. c:function:: int hsprimme_svds(float *svals, PRIMME_HALF *svecs, float *resNorms, primme_svds_params *primme_svds)
.. c:function:: int kprimme_svds(PRIMME_HALF *svals, PRIMME_COMPLEX_HALF *svecs, PRIMME_HALF *resNorms, primme_svds_params *primme_svds)
.. c:function:: int ksprimme_svds(float *svals, PRIMME_COMPLEX_HALF *svecs, float *resNorms, primme_svds_params *primme_svds)

   .. versionadded:: 3.0
.. c:function:: int sprimme_svds(float *svals, float *svecs, float *resNorms, primme_svds_params *primme_svds)
.. c:function:: int cprimme_svds(float *svals, PRIMME_COMPLEX_FLOAT *svecs, float *resNorms, primme_svds_params *primme_svds)

   .. versionadded:: 2.0
.. c:function:: int dprimme_svds(double *svals, double *svecs, double *resNorms, primme_svds_params *primme_svds)
.. c:function:: int zprimme_svds(double *svals, PRIMME_COMPLEX_DOUBLE *svecs, double *resNorms, primme_svds_params *primme_svds)

   Solve a real singular value problem.

   All arrays should be hosted on CPU. The computations are performed on CPU (see :c:func:`magma_dprimme_svds` for using GPUs).

   :param svals: array at least of size |SnumSvals| to store the
      computed singular values; all processes in a parallel run return this local array with the same values.

   :param svecs: array at least of size (|SmLocal| + |SnLocal|) times (|SnumOrthoConst| + |SnumSvals|)
      to store column-wise the (local part for this process of the) computed left singular vectors
      and the right singular vectors.

   :param resNorms: array at least of size |SnumSvals| to store the
      residual norms of the computed triplets; all processes in parallel run return this local array with
      the same values.

   :param primme_svds: parameters structure.

   :return: error indicator; see :ref:`error-codes-svds`.

   On input, ``svecs`` should start with the content of the |SnumOrthoConst| left vectors,
   followed by the |SinitSize| left vectors, followed by the |SnumOrthoConst| right vectors, and
   followed by the |SinitSize| right vectors.
 
   On return, the i-th left singular vector starts at svecs[( |SnumOrthoConst| +i)\* |SmLocal| ].
   The i-th right singular vector starts at svecs[( |SnumOrthoConst| + |SinitSize| )\* |SmLocal| + ( |SnumOrthoConst| +i)\* |SnLocal| ].
   The first vector has i=0.
 
   All internal operations are performed at the same precision than ``svecs`` unless the user sets |SinternalPrecision| otherwise.
   The functions :c:func:`hsprimme_svds` and :c:func:`ksprimme_svds` perform all computations in half precision by default and report the eigenvalues and the residual norms in single precision. These functions may help in applications that may be not built with a compiler supporting half precision.

   The type and precision of the callbacks depends on the type and precision of ``svecs``. Although this can be changed. See details for |SmatrixMatvec|, |SapplyPreconditioner|, |SglobalSumReal|, |SbroadcastReal|, and |SconvTestFun|.


`cublas_?primme_svds` & `magma_?primme_svds`
""""""""""""""""""""""""""""""""""""""""""""

.. c:function:: int cublas_hprimme_svds(PRIMME_HALF *svals, PRIMME_HALF *svecs, PRIMME_HALF *resNorms, primme_svds_params *primme_svds)
.. c:function:: int cublas_hsprimme_svds(float *svals, PRIMME_HALF *svecs, float *resNorms, primme_svds_params *primme_svds)
.. c:function:: int cublas_kprimme_svds(PRIMME_HALF *svals, PRIMME_COMPLEX_HALF *svecs, PRIMME_HALF *resNorms, primme_svds_params *primme_svds)
.. c:function:: int cublas_ksprimme_svds(float *svals, PRIMME_COMPLEX_HALF *svecs, float *resNorms, primme_svds_params *primme_svds)
.. c:function:: int cublas_sprimme_svds(float *svals, float *svecs, float *resNorms, primme_svds_params *primme_svds)
.. c:function:: int cublas_cprimme_svds(float *svals, PRIMME_COMPLEX_FLOAT *svecs, float *resNorms, primme_svds_params *primme_svds)
.. c:function:: int cublas_dprimme_svds(double *svals, double *svecs, double *resNorms, primme_svds_params *primme_svds)
.. c:function:: int cublas_zprimme_svds(double *svals, PRIMME_COMPLEX_DOUBLE *svecs, double *resNorms, primme_svds_params *primme_svds)
.. c:function:: int magma_hprimme_svds(PRIMME_HALF *svals, PRIMME_HALF *svecs, PRIMME_HALF *resNorms, primme_svds_params *primme_svds)
.. c:function:: int magma_hsprimme_svds(float *svals, PRIMME_HALF *svecs, float *resNorms, primme_svds_params *primme_svds)
.. c:function:: int magma_kprimme_svds(PRIMME_HALF *svals, PRIMME_COMPLEX_HALF *svecs, PRIMME_HALF *resNorms, primme_svds_params *primme_svds)
.. c:function:: int magma_ksprimme_svds(float *svals, PRIMME_COMPLEX_HALF *svecs, float *resNorms, primme_svds_params *primme_svds)
.. c:function:: int magma_sprimme_svds(float *svals, float *svecs, float *resNorms, primme_svds_params *primme_svds)
.. c:function:: int magma_cprimme_svds(float *svals, PRIMME_COMPLEX_FLOAT *svecs, float *resNorms, primme_svds_params *primme_svds)
.. c:function:: int magma_dprimme_svds(double *svals, double *svecs, double *resNorms, primme_svds_params *primme_svds)
.. c:function:: int magma_zprimme_svds(double *svals, PRIMME_COMPLEX_DOUBLE *svecs, double *resNorms, primme_svds_params *primme_svds)

   Solve a real singular value problem.

   Most of the computations are performed on GPU (see :c:func:`dprimme_svds` for using only the CPU).

   :param svals: CPU array at least of size |SnumSvals| to store the
      computed singular values; all processes in a parallel run return this local array with the same values.

   :param svecs: GPU array at least of size (|SmLocal| + |SnLocal|) times (|SnumOrthoConst| + |SnumSvals|)
      to store column-wise the (local part for this process of the) computed left singular vectors
      and the right singular vectors.

   :param resNorms: CPU array at least of size |SnumSvals| to store the
      residual norms of the computed triplets; all processes in parallel run return this local array with
      the same values.

   :param primme_svds: parameters structure.

   :return: error indicator; see :ref:`error-codes-svds`.

   On input, ``svecs`` should start with the content of the |SnumOrthoConst| left vectors,
   followed by the |SinitSize| left vectors, followed by the |SnumOrthoConst| right vectors, and
   followed by the |SinitSize| right vectors.

   On return, the i-th left singular vector starts at svecs[( |SnumOrthoConst| +i)\* |SmLocal| ].
   The i-th right singular vector starts at svecs[( |SnumOrthoConst| + |SinitSize| )\* |SmLocal| + ( |SnumOrthoConst| +i)\* |SnLocal| ].
   The first vector has i=0.
 
   All internal operations are performed at the same precision than ``svecs`` unless the user sets |SinternalPrecision| otherwise.
   The functions :c:func:`magma_hsprimme_svds` and :c:func:`magma_ksprimme_svds` perform all computations in half precision by default and report the eigenvalues and the residual norms in single precision. These functions may help in applications that may be not built with a compiler supporting half precision.

   The type and precision of the callbacks depends on the type and precision of ``svecs``. Although this can be changed. See details for |SmatrixMatvec|, |SapplyPreconditioner|, |SglobalSumReal|, |SbroadcastReal|, and |SconvTestFun|.

   .. versionadded:: 3.0

primme_svds_initialize
""""""""""""""""""""""

.. c:function:: void primme_svds_initialize(primme_svds_params *primme_svds)

   Initialize PRIMME SVDS parameters structure to the default values.

   After calling :c:func:`dprimme_svds` (or a variant), call :c:func:`primme_svds_free` to release allocated resources by PRIMME.

   :param primme_svds: parameters structure.

   Example::

      primme_svds_params primme_svds;
      primme_svds_initialize(&primme_svds);

      primme_svds.n = 100;
      ...
      dprimme_svds(svals, svecs, rnorms, &primme_svds);
      ...

      primme_svds_free(&primme_svds);

   See the alternative function :c:func:`primme_svds_params_create` that also allocates the structure.


primme_svds_create
""""""""""""""""""

.. c:function:: primme_svds_params* primme_svds_params_create(void)

   Allocate and initialize a parameters structure to the default values.

   After calling :c:func:`dprimme_svds` (or a variant), call :c:func:`primme_svds_params_destroy` to release allocated resources by PRIMME.

   :param primme_sv: parameters structure.

   Example::

      primme_svds_params *primme_svds = primme_svds_params_create();

      primme_svds->n = 100;
      ...
      dprimme_svds(svals, svecs, rnorms, primme_svds);
      ...

      primme_svds_params_destroy(primme_svds);

   See the alternative function :c:func:`primme_svds_initialize` that only initializes the structure.

   .. versionadded:: 3.0

primme_svds_set_method
""""""""""""""""""""""

.. c:function:: int primme_svds_set_method(primme_svds_preset_method method, primme_preset_method methodStage1, primme_preset_method methodStage2, primme_svds_params *primme_svds)

   Set PRIMME SVDS parameters to one of the preset configurations.

   :param method: preset method to compute the singular triplets; one of

      * |primme_svds_default|, currently set as |primme_svds_hybrid|.
      * |primme_svds_normalequations|, compute the eigenvectors of :math:`A^*A` or :math:`A A^*`.
      * |primme_svds_augmented|, compute the eigenvectors of the augmented matrix, :math:`\left(\begin{array}{cc} 0 & A^* \\ A & 0 \end{array}\right)`.
      * |primme_svds_hybrid|, start with |primme_svds_normalequations|; use the
        resulting approximate singular vectors as initial vectors for
        |primme_svds_augmented| if the required accuracy was not achieved.

   :param methodStage1: preset method to compute the eigenpairs at the first stage; see available values at :c:func:`primme_set_method`.

   :param methodStage2: preset method to compute the eigenpairs with
      the second stage of ``primme_svds_hybrid``; see available values at :c:func:`primme_set_method`.

   :param primme_svds: parameters structure.

   See also :ref:`methods_svds`.

primme_svds_display_params
""""""""""""""""""""""""""

.. c:function:: void primme_svds_display_params(primme_svds_params primme_svds)

   Display all printable settings of ``primme_svds`` into the file descriptor |SoutputFile|.

   :param primme_svds: parameters structure.

primme_svds_free
""""""""""""""""

.. c:function:: void primme_svds_free(primme_svds_params *primme_svds)

   Free memory allocated by PRIMME SVDS.

   :param primme_svds: parameters structure.

primme_svds_params_destroy
""""""""""""""""""""""""""

.. c:function:: int primme_svds_params_destroy(primme_svds_params *primme)

   Free memory allocated by PRIMME associated to a parameters structure created
   with :c:func:`primme_svds_params_create`.

   :param primme_svds: parameters structure.

   :return: nonzero value if the call is not successful.

   .. versionadded:: 3.0

.. include:: epilog.inc
