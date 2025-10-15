
.. highlight:: c

Parameters Description
----------------------

Types
"""""

The following data types are macros used in PRIMME as followed.

.. c:type:: PRIMME_INT

   Integer type used in matrix dimensions (such as |n| and |nLocal|) and counters (such as |numMatvecs|).

   The integer size is controlled by the compilation flag  ``PRIMME_INT_SIZE``, see :ref:`making`.

.. c:type:: PRIMME_HALF

   Macro that is ``__fp16`` if half precision is supported by the compiler.
   Otherwise it is a struct with the same size as ``int short``.

   .. versionadded:: 3.0

.. c:type:: PRIMME_COMPLEX_HALF

   Macro that is a struct with fields ``r`` and ``i`` with type ``PRIMME_HALF``.

   .. versionadded:: 3.0

.. c:type:: PRIMME_COMPLEX_FLOAT

   Macro that is ``complex float`` in C and ``std::complex<float>`` in C++.

   .. versionadded:: 2.0

.. c:type:: PRIMME_COMPLEX_DOUBLE

   Macro that is ``complex double`` in C and ``std::complex<double>`` in C++.

   .. versionadded:: 2.0

Other macros
""""""""""""

.. c:macro:: PRIMME_VERSION_MAJOR

   Constant ``int`` with the major version number.

   For instance, the value of the macro is ``3`` for version 3.0.

   .. versionadded:: 3.0

.. c:macro:: PRIMME_VERSION_MINOR

   Constant ``int`` with the minor version number.

   For instance, the value of the macro is ``0`` for version 3.0.

   .. versionadded:: 3.0


primme_params
"""""""""""""

.. c:type:: primme_params

   Structure to set the problem matrices and eigensolver options.

   .. c:member:: PRIMME_INT n

      Dimension of the matrix.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | this field is read by :c:func:`dprimme`.

   .. c:member:: void (*matrixMatvec) (void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr)

      Block matrix-multivector multiplication, :math:`y = A x` in solving :math:`A x = \lambda x` or :math:`A x = \lambda B x`.
   
      :param x: matrix of size |nLocal| x ``blockSize`` in column-major_ order with leading dimension ``ldx``.
      :param ldx: the leading dimension of the array ``x``.
      :param y: matrix of size |nLocal| x ``blockSize`` in column-major_ order with leading dimension ``ldy``.
      :param ldy: the leading dimension of the array ``y``.
      :param blockSize: number of columns in ``x`` and ``y``.
      :param primme: parameters structure.
      :param ierr: output error code; if it is set to non-zero, the current call to PRIMME will stop.

      The actual type of ``x`` and ``y`` matches the type of ``evecs`` of the
      calling  :c:func:`dprimme` (or a variant), unless the user sets |matrixMatvec_type| to
      another precision.

      Input/output:

         | :c:func:`primme_initialize` sets this field to NULL;
         | this field is read by :c:func:`dprimme`.

      .. note::

            If you have performance issues with leading dimension different from |nLocal|,
            set |ldOPs| to |nLocal|.

      .. versionchanged:: 2.0

   .. c:member:: primme_op_datatype matrixMatvec_type

      Precision of the vectors ``x`` and ``y`` passed to |matrixMatvec|.

      If it is ``primme_op_default``, the vectors' type matches the calling
      :c:func:`dprimme` (or a variant). Otherwise, the user can force the precision of the vectors ``x`` and ``y`` to be a particular precision regardless of the calling :c:func:`dprimme` (or a variant) function: half,
      single, or double, if |matrixMatvec_type| is ``primme_half``, ``primme_float``
      or ``primme_double`` respectively.

      It is not recommended to set a lower precision than the one required to converge.
      An example of this is calling :c:func:`dprimme` setting |eps| to 1e-10 and |matrixMatvec_type| to ``primme_op_half``.

      Input/output:

         | :c:func:`primme_initialize` sets this field to ``primme_op_default``;
         | this field is read by :c:func:`dprimme`, and if it is
           ``primme_op_default`` it is set to the value that matches the precision of
           calling function.

      .. versionadded:: 3.0

   .. c:member:: void (*applyPreconditioner)(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr)

      Block preconditioner-multivector application, :math:`y = M^{-1}x` where :math:`M` is usually an approximation of :math:`A - \sigma I` or :math:`A - \sigma B` for finding eigenvalues close to :math:`\sigma`.

      :param x: matrix of size |nLocal| x ``blockSize`` in column-major_ order with leading dimension ``ldx``.
      :param ldx: the leading dimension of the array ``x``.
      :param y: matrix of size |nLocal| x ``blockSize`` in column-major_ order with leading dimension ``ldy``.
      :param ldy: the leading dimension of the array ``y``.
      :param blockSize: number of columns in ``x`` and ``y``.
      :param primme: parameters structure.
      :param ierr: output error code; if it is set to non-zero, the current call to PRIMME will stop.

      The actual type of ``x`` and ``y`` matches the type of ``evecs`` of the
      calling  :c:func:`dprimme` (or a variant), unless the user sets |applyPreconditioner_type| to
      another precision.

      Input/output:

         | :c:func:`primme_initialize` sets this field to NULL;
         | this field is read by :c:func:`dprimme`.
 
      .. versionchanged:: 2.0

   .. c:member:: primme_op_datatype applyPreconditioner_type

      Precision of the vectors ``x`` and ``y`` passed to |applyPreconditioner|.

      If it is ``primme_op_default``, the vectors' type matches the calling
      :c:func:`dprimme` (or a variant). Otherwise, the user can force the precision of the vectors ``x`` and ``y`` to be a particular precision regardless of the calling :c:func:`dprimme` (or a variant) function: half,
      single, or double, if |matrixMatvec_type| is ``primme_half``, ``primme_float``
      or ``primme_double`` respectively.

      Input/output:

         | :c:func:`primme_initialize` sets this field to ``primme_op_default``;
         | this field is read by :c:func:`dprimme`, and if it is
           ``primme_op_default`` it is set to the value that matches the precision of
           calling function.

      .. versionadded:: 3.0

   .. c:member:: void (*massMatrixMatvec) (void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr)

      Block matrix-multivector multiplication, :math:`y = B x` in solving :math:`A x = \lambda B x`.
      If it is NULL, the standard eigenvalue problem :math:`A x = \lambda x` is solved.

      :param x: matrix of size |nLocal| x ``blockSize`` in column-major_ order with leading dimension ``ldx``.
      :param ldx: the leading dimension of the array ``x``.
      :param y: matrix of size |nLocal| x ``blockSize`` in column-major_ order with leading dimension ``ldy``.
      :param ldy: the leading dimension of the array ``y``.
      :param blockSize: number of columns in ``x`` and ``y``.
      :param primme: parameters structure.
      :param ierr: output error code; if it is set to non-zero, the current call to PRIMME will stop.

      The actual type of ``x`` and ``y`` matches the type of ``evecs`` of the
      calling  :c:func:`dprimme` (or a variant), unless the user sets |matrixMatvec_type| to
      another precision.

      Input/output:

         | :c:func:`primme_initialize` sets this field to NULL;
         | this field is read by :c:func:`dprimme`.

      .. versionchanged:: 2.0

   .. c:member:: primme_op_datatype massMatrixMatvec_type

      Precision of the vectors ``x`` and ``y`` passed to |massMatrixMatvec|.

      If it is ``primme_op_default``, the vectors' type matches the calling
      :c:func:`dprimme` (or a variant). Otherwise, the user can force the precision of the vectors ``x`` and ``y`` to be a particular precision regardless of the calling :c:func:`dprimme` (or a variant) function: half,
      single, or double, if |massMatrixMatvec_type| is ``primme_half``, ``primme_float``
      or ``primme_double`` respectively.

      It is not recommended to set a lower precision than the one required to converge.
      An example of this is calling :c:func:`dprimme` setting |eps| to 1e-10 and |massMatrixMatvec_type| to ``primme_op_half``.

      Input/output:

         | :c:func:`primme_initialize` sets this field to ``primme_op_default``;
         | this field is read by :c:func:`dprimme`, and if it is
           ``primme_op_default`` it is set to the value that matches the precision of
           calling function.

      .. versionadded:: 3.0

   .. c:member:: int numProcs

      Number of processes calling :c:func:`dprimme` or :c:func:`zprimme` in parallel.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 1;
         | this field is read by :c:func:`dprimme`.

   .. c:member:: int procID

      The identity of the local process within a parallel execution calling :c:func:`dprimme` or
      :c:func:`zprimme`.
      Only the process with id 0 prints information.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | :c:func:`dprimme` sets this field to 0 if |numProcs| is 1;
         | this field is read by :c:func:`dprimme`.

   .. c:member:: int nLocal

      Number of local rows on this process.

      Input/output:

         | :c:func:`primme_initialize` sets this field to -1;
         | :c:func:`dprimme` sets this field to |n| if |numProcs| is 1;
         | this field is read by :c:func:`dprimme`.

   .. c:member:: void *commInfo

      A pointer to whatever parallel environment structures needed.
      For example, with MPI, it could be a pointer to the MPI communicator.
      PRIMME does not use this. It is available for possible use in 
      user functions defined in |matrixMatvec|,
      |applyPreconditioner|, |massMatrixMatvec|,
      |globalSumReal|, and |broadcastReal|.

      Input/output:

         | :c:func:`primme_initialize` sets this field to NULL;

   .. c:member:: void (*globalSumReal)(void *sendBuf, void *recvBuf, int *count, primme_params *primme, int *ierr)

      Global sum reduction function. No need to set for sequential programs.

      :param sendBuf: array of size ``count`` with the local input values.
      :param recvBuf: array of size ``count`` with the global output values
         so that the i-th element of recvBuf is the sum over all processes of the i-th element
         of ``sendBuf``.
      :param count: array size of ``sendBuf`` and ``recvBuf``.
      :param primme: parameters structure.
      :param ierr: output error code; if it is set to non-zero, the current call to PRIMME will stop.

      The actual type of ``sendBuf`` and ``recvBuf`` matches the type of ``evecs`` of the
      calling  :c:func:`dprimme` (or a variant), unless the user sets |globalSumReal_type| to
      another precision. See the recomendation about precision in |globalSumReal_type|.

      Note that ``count`` is the number of values of the real type.
 
      Input/output:

         | :c:func:`primme_initialize` sets this field to an internal function;
         | :c:func:`dprimme` sets this field to an internal function if |numProcs| is 1 and |globalSumReal| is NULL;
         | this field is read by :c:func:`dprimme`.

      When MPI is used, this can be a simply wrapper to MPI_Allreduce() as shown below:

      .. code:: c

         void par_GlobalSumForDouble(void *sendBuf, void *recvBuf, int *count, 
                                  primme_params *primme, int *ierr) {
            MPI_Comm communicator = *(MPI_Comm *) primme->commInfo;
            if (sendBuf == recvBuf) {
              *ierr = MPI_Allreduce(MPI_IN_PLACE, recvBuf, *count, MPIU_REAL, MPI_SUM, communicator) != MPI_SUCCESS;
            } else {
              *ierr = MPI_Allreduce(sendBuf, recvBuf, *count, MPIU_REAL, MPI_SUM, communicator) != MPI_SUCCESS;
            }
         }

      When calling :c:func:`sprimme` and :c:func:`cprimme` replace ``MPI_DOUBLE`` by ```MPI_FLOAT``.

      .. versionchanged:: 2.0

   .. c:member:: primme_op_datatype globalSumReal_type

      Precision of the vectors ``sendBuf`` and ``recvBuf`` passed to |globalSumReal|.

      If it is ``primme_op_default``, the vectors' type matches the calling
      :c:func:`dprimme` (or a variant). Otherwise, the user can force the precision of the vectors ``sendBuf`` and ``recvBuf`` to be a particular precision regardless of the calling :c:func:`dprimme` (or a variant) function: half,
      single, or double, if |globalSumReal_type| is ``primme_half``, ``primme_float``
      or ``primme_double`` respectively.

      It is recommended to set a precision so that the machine precision times log2(|numProcs|) is smaller than the precision required to converge.
      An example of this is calling :c:func:`hprimme` setting |eps| to 0.01 and |globalSumReal_type| to ``primme_op_single`` for 1000 processes.


      Input/output:

         | :c:func:`primme_initialize` sets this field to ``primme_op_default``;
         | this field is read by :c:func:`dprimme`, and if it is
           ``primme_op_default`` it is set to the value that matches the precision of
           calling function.

      .. versionadded:: 3.0

   .. c:member:: void (*broadcastReal)(void *buffer, int *count, primme_params *primme, int *ierr)

      Broadcast function from process with ID zero. It is optional in parallel executions, and not needed for sequential programs.

      :param buffer: array of size ``count`` with the local input values.
      :param count: array size of ``sendBuf`` and ``recvBuf``.
      :param primme: parameters structure.
      :param ierr: output error code; if it is set to non-zero, the current call to PRIMME will stop.

      The actual type of ``buffer`` matches the type of ``evecs`` of the
      calling  :c:func:`dprimme` (or a variant), unless the user sets |broadcastReal_type| to
      another precision.

      If |broadcastReal| is not provided, PRIMME uses |globalSumReal| for broadcasting, which is usually a bit more expensive.

      Input/output:

         | :c:func:`primme_initialize` sets this field to NULL;
         | this field is read by :c:func:`dprimme`.

      When MPI is used, this can be a simply wrapper to MPI_Bcast() as shown below:

      .. code:: c

         void broadcastForDouble(void *buffer, int *count, 
                                 primme_params *primme, int *ierr) {
            MPI_Comm communicator = *(MPI_Comm *) primme->commInfo;
            if(MPI_Bcast(buffer, *count, MPI_DOUBLE, 0 /* root */,
                          communicator) == MPI_SUCCESS) {
               *ierr = 0;
            } else {
               *ierr = 1;
            }
         }

      When calling :c:func:`sprimme` and :c:func:`cprimme` replace ``MPI_DOUBLE`` by ```MPI_FLOAT``.

      .. versionadded:: 3.0

   .. c:member:: primme_op_datatype broadcastReal_type

      Precision of the vector ``buffer``` passed to |broadcastReal|.

      If it is ``primme_op_default``, the vectors' type matches the calling
      :c:func:`dprimme` (or a variant). Otherwise, the user can force the precision of the vectors ``x`` and ``y`` to be a particular precision regardless of the calling :c:func:`dprimme` (or a variant) function: half,
      single, or double, if |matrixMatvec_type| is ``primme_half``, ``primme_float``
      or ``primme_double`` respectively.

      It is not recommended to set a lower precision than the one required to converge.
      An example of this is calling :c:func:`dprimme` setting |eps| to 1e-10 and |broadcastReal_type| to ``primme_op_half``.

      Input/output:

         | :c:func:`primme_initialize` sets this field to ``primme_op_default``;
         | this field is read by :c:func:`dprimme`, and if it is
           ``primme_op_default`` it is set to the value that matches the precision of
           calling function.

      .. versionadded:: 3.0

   .. c:member:: primme_op_datatype internalPrecision

      Internal working precision.

      If it is ``primme_op_default``, most of the vectors are stored with the
      same precision as the calling :c:func:`dprimme` (or a variant), and most of the
      computations are done in that precision too. Otherwise, the working precision
      is changed to half, single, or double, if the user sets |internalPrecision| to
      ``primme_half``, ``primme_float`` or ``primme_double`` respectively.

      Input/output:

         | :c:func:`primme_initialize` sets this field to ``primme_op_default``;
         | this field is read by :c:func:`dprimme`.

      .. versionadded:: 3.0

   .. c:member:: int numEvals

      Number of eigenvalues wanted.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 1;
         | this field is read by :c:func:`primme_set_method` (see :ref:`methods`) and :c:func:`dprimme`.

   .. index:: interior problem

   .. c:member:: primme_target target

      Which eigenpairs to find:

      ``primme_smallest``
         Smallest algebraic eigenvalues; |targetShifts| is ignored.

      ``primme_largest``
         Largest algebraic eigenvalues; |targetShifts| is ignored.

      ``primme_closest_geq``
         Closest to, but greater or equal than the shifts in |targetShifts|.

      ``primme_closest_leq``
         Closest to, but less or equal than the shifts in |targetShifts|.

      ``primme_closest_abs``
         Closest in absolute value to the shifts in |targetShifts|.

      ``primme_largest_abs``
         Furthest in absolute value to the shifts in |targetShifts|.

      Input/output:

         | :c:func:`primme_initialize` sets this field to |primme_smallest|;
         | this field is read by :c:func:`dprimme`.

   .. index:: interior problem

   .. c:member:: int numTargetShifts
 
      Size of the array |targetShifts|.
      Used only when |target| is |primme_closest_geq|,
      |primme_closest_leq|, |primme_closest_abs| or |primme_largest_abs|.
      The default values is 0.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | this field is read by :c:func:`dprimme`.

   .. index:: interior problem

   .. c:member:: double *targetShifts

      Array of shifts, at least of size |numTargetShifts|.
      Used only when |target| is |primme_closest_geq|,
      |primme_closest_leq|, |primme_closest_abs| or |primme_largest_abs|.

      Eigenvalues are computed in order so that the
      i-th eigenvalue is the closest (or closest but left or closest but right, see |target|)
      to the i-th shift. If |numTargetShifts| < |numEvals|, the last shift given
      is used for all the remaining i's.

      Input/output:

         | :c:func:`primme_initialize` sets this field to NULL;
         | this field is read by :c:func:`dprimme`.

      .. note::

         Considerations for interior problems:

         * PRIMME will try to compute the eigenvalues in the order given in the
           |targetShifts|. However, for code efficiency and robustness, the shifts
           should be ordered. Order them in ascending (descending) order for shifts
           closer to the lower (higher) end of the spectrum.
         * If some shift is close to the lower (higher) end of the spectrum,
           use either |primme_closest_geq| (|primme_closest_leq|) or
           |primme_closest_abs|.
         * |primme_closest_leq| and |primme_closest_geq| are more efficient
           than |primme_closest_abs|.
         * For interior eigenvalues larger |maxBasisSize| is usually more robust.
         * To find the largest magnitude eigenvalues set |target| to |primme_largest_abs|,
           |numTargetShifts| to 1 and |targetShifts| to an array with a zero value.

   .. c:member:: int printLevel

      The level of message reporting from the code. All output is written in |outputFile|.

      One of:
 
      * 0: silent.
      * 1: print some error messages when these occur.
      * 2: as in 1, and info about targeted eigenpairs when they are marked as converged::
      
            #Converged $1 eval[ $2 ]= $3 norm $4 Mvecs $5 Time $7

        or locked::

            #Lock epair[ $1 ]= $3 norm $4 Mvecs $5 Time $7

      * 3: in as 2, and info about targeted eigenpairs every outer iteration::
      
            OUT $6 conv $1 blk $8 MV $5 Sec $7 EV $3 |r| $4

        Also, if it is used the dynamic method, show JDQMR/GDk performance ratio and
        the current method in use.
      * 4: in as 3, and info about targeted eigenpairs every inner iteration::
      
            INN MV $5 Sec $7 Eval $3 Lin|r| $9 EV|r| $4
      
      * 5: in as 4, and verbose info about certain choices of the algorithm.
      
      Output key:

      | $1: Number of converged pairs up to now.
      | $2: The index of the pair currently converged.
      | $3: The eigenvalue.
      | $4: Its residual norm.
      | $5: The current number of matrix-vector products.
      | $6: The current number of outer iterations.
      | $7: The current elapsed time.
      | $8: Index within the block of the targeted pair .
      | $9: QMR norm of the linear system residual.

      In parallel programs, output is produced in call with
      |procID| 0 when |printLevel|
      is from 0 to 4.
      If |printLevel| is 5 output can be produced in any of
      the parallel calls.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 1;
         | this field is read by :c:func:`dprimme`.

   .. note::

      Convergence history for plotting may be produced simply by:

      .. code-block:: bash

         grep OUT outpufile | awk '{print $8" "$14}' > out
         grep INN outpufile | awk '{print $3" "$11}' > inn

      Then in gnuplot:

      .. code-block:: gnuplot

         plot 'out' w lp, 'inn' w lp

   .. c:member:: double aNorm

      An estimate of the norm of :math:`A`, which is used in the default convergence
      criterion (see |eps|).

      If |aNorm| is less than or equal to 0, the code uses the largest absolute
      Ritz value seen divided by |invBNorm|. On return, |aNorm| is then replaced with that value.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0.0;
         | this field is read and written by :c:func:`dprimme`.

   .. c:member:: double BNorm

      An estimate of the norm of :math:`B`, which is used to estimate the
      conditioning number of the matrix :math:`B`.

      If |BNorm| is less than or equal to 0, the code uses the largest inner-product with :math:`B`
      seen. On return, |BNorm| is then replaced with that value.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0.0;
         | this field is read and written by :c:func:`dprimme`.

      .. versionadded:: 3.0

   .. c:member:: double invBNorm

      An estimate of the norm of the inverse of :math:`B`, which is used in the default convergence
      criterion (see |eps|).

      If |invBNorm| is less than or equal to 0, the code uses the inverse of the smallest inner-product with :math:`B`
      seen. On return, |invBNorm| is then replaced with that value.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0.0;
         | this field is read and written by :c:func:`dprimme`.

      .. versionadded:: 3.0

   .. c:member:: primme_orth orth

      Selects the orthogonalization method used by PRIMME.

      If the value is ``primme_orth_implicit_I``, the bases are orthogonalized with
      classical Gram-Schmidt with reorthogonalization stopping when the new vector's
      norm is not reduced more than :math:`1/sqrt{2}` (Daniel's test) from the previous
      iteration. If several vectors are going to be orthogonalized, the algorithm is
      applied vector by vector.

      If the value is ``primme_orth_explicit_I``, the bases are orthogonalized with
      iterative Cholesky QR (or SVQB if Cholesky factorization fails), stopping
      when the conditioning of the basis is around :math:`sqrt(3)`. That
      deviation of the orthogonality level is taken into account in the Galerkin
      method.
      
      The option ``primme_orth_explicit_I`` is usually more expensive in FLOPS,
      but it may be faster in time than ``primme_orth_implicit_I`` when |maxBlockSize|
      is large.

      ``primme_orth_implicit_I`` is set by default if the precision is higher than
      single precision and |maxBlockSize| is 1. Otherwise, ``primme_orth_explicit_I``
      is set by default.

      Input/output:

         | :c:func:`primme_initialize` sets this field to ``primme_orth_default``;
         | this field is read and written by :c:func:`dprimme`.

      .. versionadded:: 3.0

   .. c:member:: double eps

      If |convTestFun| is NULL, an eigenpairs is marked as converged when the 2-norm
      of the residual vector is less than |eps| \* |aNorm| \* |invBNorm|.
      The residual vector is :math:`A x - \lambda x` or :math:`A x - \lambda B x`.

      The default value is machine precision times :math:`10^4`.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0.0;
         | this field is read and written by :c:func:`dprimme`.
 
   .. c:member:: FILE *outputFile

      Opened file to write down the output.

      Input/output:

         | :c:func:`primme_initialize` sets this field to the standard output;
         | this field is read by :c:func:`dprimme` and :c:func:`primme_display_params`.

   .. c:member:: int dynamicMethodSwitch

      If this value is 1, it alternates dynamically between |DEFAULT_MIN_TIME|
      and |DEFAULT_MIN_MATVECS|, trying to identify the fastest method.

      On exit, it holds a recommended method for future runs on this problem:

         | -1: use |DEFAULT_MIN_MATVECS| next time.
         | -2: use |DEFAULT_MIN_TIME| next time.
         | -3: close call, use |DYNAMIC| next time again.
      
      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`primme_set_method` (see :ref:`methods`);
         | this field is read by :c:func:`dprimme`.

      .. note::

         Even for expert users we do not recommend setting |dynamicMethodSwitch|
         directly, but through :c:func:`primme_set_method`.

      .. note::

         The code obtains timings by the ``gettimeofday`` Unix utility. If a cheaper, more
         accurate timer is available, modify the ``PRIMMESRC/COMMONSRC/wtime.c``

   .. c:member:: int locking

      If set to 1, hard locking will be used (locking converged eigenvectors
      out of the search basis). If set to 0, the code will try to use soft
      locking (à la ARPACK), when large enough |minRestartSize| is available.

      Input/output:

         | :c:func:`primme_initialize` sets this field to -1;
         | written by :c:func:`primme_set_method` (see :ref:`methods`);
         | this field is read by :c:func:`dprimme`.

   .. c:member:: int initSize
 
      On input, the number of initial vector guesses provided in ``evecs`` argument in
      :c:func:`dprimme` or :c:func:`zprimme`.

      On output, |initSize| holds the number of converged eigenpairs. Without |locking| all
      |numEvals| approximations are in ``evecs`` but only the |initSize| ones are
      converged.

      During execution, it holds the current number of converged eigenpairs.
      In addition, if locking is used, these are accessible in ``evals`` and ``evecs``.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | this field is read and written by :c:func:`dprimme`.

   .. c:member:: PRIMME_INT ldevecs

      The leading dimension of ``evecs``. The default is |nLocal|.

      Input/output:

         | :c:func:`primme_initialize` sets this field to -1;
         | this field is read by :c:func:`dprimme`.

      .. versionadded:: 2.0

   .. c:member:: int numOrthoConst

      Number of vectors to be used as external orthogonalization constraints.
      These vectors are provided in the first |numOrthoConst| positions of
      the ``evecs`` argument in :c:func:`dprimme` or :c:func:`zprimme` and must be orthonormal.

      PRIMME finds new eigenvectors orthogonal to these constraints (equivalent to solving
      the problem with :math:`(I-YY^*)A(I-YY^*)` and :math:`(I-YY^*)B(I-YY^*)` matrices
      where :math:`Y` are the given constraint vectors).
      This is a handy feature if some eigenvectors are already known, or 
      for finding more eigenvalues after a call to :c:func:`dprimme` or :c:func:`zprimme`,
      possibly with different parameters (see an example in :file:`TEST/ex_zseq.c`).

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | this field is read by :c:func:`dprimme`.

   .. c:member:: int maxBasisSize

      The maximum basis size allowed in the main iteration. This has memory
      implications.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | this field is read and written by :c:func:`primme_set_method` (see :ref:`methods`);
         | this field is read by :c:func:`dprimme`.

   .. c:member:: int minRestartSize

      Maximum Ritz vectors kept after restarting the basis.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | this field is read and written by :c:func:`primme_set_method` (see :ref:`methods`);
         | this field is read by :c:func:`dprimme`.

   .. c:member:: int maxBlockSize
 
      The maximum block size the code will try to use.

      The user should set
      this based on the architecture specifics of the target computer, 
      as well as any a priori knowledge of multiplicities. The code does 
      *not* require that |maxBlockSize| > 1 to find multiple eigenvalues. For some 
      methods, keeping to 1 yields the best overall performance.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 1;
         | this field is read and written by :c:func:`primme_set_method` (see :ref:`methods`);
         | this field is read by :c:func:`dprimme`.

      .. note::

         Inner iterations of QMR are not performed in a block fashion.
         Every correction equation from a block is solved independently.

   .. index:: stopping criterion

   .. c:member:: PRIMME_INT maxMatvecs

      Maximum number of matrix vector multiplications (approximately equal to 
      the number of preconditioning operations) that the code is allowed to 
      perform before it exits.

      Input/output:

         | :c:func:`primme_initialize` sets this field to ``INT_MAX``;
         | this field is read by :c:func:`dprimme`.

   .. index:: stopping criterion

   .. c:member:: PRIMME_INT maxOuterIterations

      Maximum number of outer iterations that the code is allowed to perform 
      before it exits.

      Input/output:

         | :c:func:`primme_initialize` sets this field to ``INT_MAX``;
         | this field is read by :c:func:`dprimme`.

   .. c:member:: PRIMME_INT iseed

      The ``PRIMME_INT iseed[4]`` is an array with the seeds needed by the LAPACK_ dlarnv and zlarnv.

      The default value is an array with values -1, -1, -1 and -1. In that case, ``iseed``
      is set based on the value of |procID| to avoid every parallel process generating
      the same sequence of pseudorandom numbers.

      Input/output:

         | :c:func:`primme_initialize` sets this field to ``[-1, -1, -1, -1]``;
         | this field is read and written by :c:func:`dprimme`.

   .. c:member:: void *matrix

      This field may be used to pass any required information 
      in the matrix-vector product |matrixMatvec|.

      Input/output:

         | :c:func:`primme_initialize` sets this field to NULL;

   .. c:member:: void *massMatrix

      This field may be used to pass any required information 
      in the matrix-vector product |massMatrixMatvec|.

      Input/output:

         | :c:func:`primme_initialize` sets this field to NULL;
       
      .. versionadded:: 3.0

   .. c:member:: void *preconditioner

      This field may be used to pass any required information 
      in the preconditioner function |applyPreconditioner|.

      Input/output:

         | :c:func:`primme_initialize` sets this field to NULL;

   .. c:member:: double *ShiftsForPreconditioner

      Array of size ``blockSize`` provided during execution of :c:func:`dprimme` and
      :c:func:`zprimme` holding
      the shifts to be used (if needed) in the preconditioning operation.

      For example if the block size is 3,
      there will be an array of three shifts in |ShiftsForPreconditioner|.
      Then the user can invert a shifted preconditioner for each of the 
      block vectors :math:`(M-ShiftsForPreconditioner_i)^{-1} x_i`.
      Classical Davidson (diagonal) preconditioning is an example of this.
   
      | this field is read and written by :c:func:`dprimme`.

   .. c:member:: primme_init initBasisMode

      Select how the search subspace basis is initialized up to |minRestartSize| vectors
      if not enough initial vectors are provided (see |initSize|):

      * ``primme_init_krylov``, with a block Krylov subspace generated by the matrix problem and
        the last initial vectors if given or a random vector otherwise; the size of the block is |maxBlockSize|.
      * ``primme_init_random``, with random vectors.
      * ``primme_init_user``, the initial basis will have only initial vectors if given,
        or a single random vector.

      Input/output:

         | :c:func:`primme_initialize` sets this field to |primme_init_krylov|;
         | this field is read by :c:func:`dprimme`.

      .. versionadded:: 2.0

   .. c:member:: primme_projection projectionParams.projection

      Select the extraction technique, i.e., how the approximate eigenvectors :math:`x_i` and
      eigenvalues :math:`\lambda_i` are computed from the search subspace :math:`\mathcal V`:

      * ``primme_proj_RR``, Rayleigh-Ritz, :math:`Ax_i - Bx_i\lambda_i \perp \mathcal V`.
      * ``primme_proj_harmonic``, Harmonic Rayleigh-Ritz,
        :math:`Ax_i - Bx_i\lambda_i \perp (A-\tau B)\mathcal V`, where :math:`\tau` is the current
        target shift (see |targetShifts|).
      * ``primme_proj_refined``, refined extraction, compute :math:`x_i` with :math:`||x_i||=1` that
        minimizes :math:`||(A-\tau B)x_i||`; the eigenvalues are computed as the
        Rayleigh quotients, :math:`\lambda_i=\frac{x_i^*Ax_i}{x_i^*Bx_i}`.

      Input/output:

         | :c:func:`primme_initialize` sets this field to |primme_proj_default|;
         | :c:func:`primme_set_method` and :c:func:`dprimme` sets it to |primme_proj_RR| if it is |primme_proj_default|.
 
   .. c:member:: int restartingParams.maxPrevRetain

      Number of approximations from previous iteration to be retained
      after restart (this is the locally optimal restarting, see [r2]_).
      The restart size is |minRestartSize| plus |maxPrevRetain|.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | this field is read and written by :c:func:`primme_set_method` (see :ref:`methods`);
         | this field is read by :c:func:`dprimme`.

   .. c:member:: int correctionParams.precondition

      Set to 1 to use preconditioning.
      Make sure |applyPreconditioner| is not NULL then!

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | this field is read and written by :c:func:`primme_set_method` (see :ref:`methods`);
         | this field is read by :c:func:`dprimme`.

   .. c:member:: int correctionParams.robustShifts

      Set to 1 to use robust shifting. It tries to avoid stagnation and 
      misconvergence by providing as shifts in |ShiftsForPreconditioner|
      the Ritz values displaced by an approximation of the eigenvalue error.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`primme_set_method` (see :ref:`methods`);
         | this field is read by :c:func:`dprimme`.

   .. c:member:: int correctionParams.maxInnerIterations

      Control the maximum number of inner QMR iterations:

      * 0:  no inner iterations;
      * >0: perform at most that number of inner iterations per outer step;
      * <0: perform at most the rest of the remaining matrix-vector products
        up to reach |maxMatvecs|.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | this field is read and written by :c:func:`primme_set_method` (see :ref:`methods`);
         | this field is read by :c:func:`dprimme`.

      See also |convTest|.

   .. c:member:: double correctionParams.relTolBase

      Parameter used when |convTest|
      is |primme_decreasing_LTolerance|.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`primme_set_method` (see :ref:`methods`);
         | this field is read by :c:func:`dprimme`.

   .. c:member:: primme_convergencetest correctionParams.convTest

      Set how to stop the inner QMR method:

      * ``primme_full_LTolerance``: stop by iterations only;

      * ``primme_decreasing_LTolerance``, stop when
        :math:`\text{relTolBase}^{-\text{outIts}}` where outIts
        is the number of outer iterations and retTolBase is set in
        |relTolBase|;
        This is a legacy option from classical JDQR and we recommend
        **strongly** against its use.

      * ``primme_adaptive``, stop when the estimated eigenvalue residual
        has reached the required tolerance (based on Notay's JDCG).

      * ``primme_adaptive_ETolerance``, as |primme_adaptive| but also
        stopping when the estimated eigenvalue residual has reduced 10
        times.

      Input/output:

         | :c:func:`primme_initialize` sets this field to ``primme_adaptive_ETolerance``;
         | written by :c:func:`primme_set_method` (see :ref:`methods`);
         | this field is read by :c:func:`dprimme`.

      .. note::

         Avoid to set |maxInnerIterations| to -1 and |convTest| to |primme_full_LTolerance|.

      See also |maxInnerIterations|.

   .. c:member:: int correctionParams.projectors.LeftQ
   .. c:member:: int correctionParams.projectors.LeftX
   .. c:member:: int correctionParams.projectors.RightQ
   .. c:member:: int correctionParams.projectors.RightX
   .. c:member:: int correctionParams.projectors.SkewQ
   .. c:member:: int correctionParams.projectors.SkewX

      Control the projectors involved in the computation of the correction
      appended to the basis every (outer) iteration.

      Consider the current selected Ritz value :math:`\Lambda` and vectors :math:`X`,
      the residual associated vectors :math:`R=AX-X\Lambda`, the previous locked vectors
      :math:`Q`, and the preconditioner :math:`M^{-1}`.

      When |maxInnerIterations| is 0, the correction :math:`D` appended to the basis
      in GD is:

      .. only:: not text

         +--------+-------+------------------------------------------------------------------+
         | RightX | SkewX |        :math:`D`                                                 |
         +========+=======+==================================================================+
         |    0   |   0   | :math:`M^{-1}R` (Classic GD)                                     |
         +--------+-------+------------------------------------------------------------------+
         |    1   |   0   | :math:`M^{-1}(R-\Delta X)` (cheap Olsen's Method)                |
         +--------+-------+------------------------------------------------------------------+
         |    1   |   1   | :math:`(I- M^{-1}X(X^*M^{-1}X)^{-1}X^*)M^{-1}R` (Olsen's Method) |
         +--------+-------+------------------------------------------------------------------+
         |    0   |   1   | error                                                            |
         +--------+-------+------------------------------------------------------------------+

      .. only:: text

         +--------+-------+----------------------------------------------------------+
         | RightX | SkewX | D                                                        |
         +========+=======+==========================================================+
         |    0   |   0   | M^{-1}R (Classic GD)                                     |
         +--------+-------+----------------------------------------------------------+
         |    1   |   0   | M^{-1}(R-\Delta X) (cheap Olsen's Method)                |
         +--------+-------+----------------------------------------------------------+
         |    1   |   1   | (I- M^{-1}X(X^*M^{-1}X)^{-1}X^*)M^{-1}R (Olsen's Method) |
         +--------+-------+----------------------------------------------------------+
         |    0   |   1   | error                                                    |
         +--------+-------+----------------------------------------------------------+


      Where :math:`\Delta` is a diagonal matrix that :math:`\Delta_{i,i}` holds an estimation
      of the error of the approximate eigenvalue :math:`\Lambda_{i,i}`.
 
      The values of ``RightQ``, ``SkewQ``, ``LeftX`` and ``LeftQ`` are ignored.

      When |maxInnerIterations| is not 0, the correction :math:`D` in Jacobi-Davidson results
      from solving:

      .. math::

         P_Q^l P_X^l (A-\sigma I) P_X^r P_Q^r M^{-1} D' = -R, \ \ \  D = P_X^r P_Q^l M^{-1}D'.

      For ``LeftQ``:

         | 0: :math:`P_Q^l = I`;
         | 1: :math:`P_Q^l = I - QQ^*`.

      For ``LeftX``:

         | 0: :math:`P_X^l = I`;
         | 1: :math:`P_X^l = I - XX^*`.

      For ``RightQ`` and ``SkewQ``:

      +--------+-------+-------------------------------+
      | RightQ | SkewQ |        :math:`P_Q^r`          |
      +========+=======+===============================+
      |    0   |   0   | :math:`I`                     |
      +--------+-------+-------------------------------+
      |    1   |   0   | :math:`I - QQ^*`              |
      +--------+-------+-------------------------------+
      |    1   |   1   | :math:`I - KQ(Q^*KQ)^{-1}Q^*` |
      +--------+-------+-------------------------------+
      |    0   |   1   | error                         |
      +--------+-------+-------------------------------+

      For ``RightX`` and ``SkewX``:

      +--------+-------+-------------------------------+
      | RightX | SkewX |        :math:`P_X^r`          |
      +========+=======+===============================+
      |    0   |   0   | :math:`I`                     |
      +--------+-------+-------------------------------+
      |    1   |   0   | :math:`I - XX^*`              |
      +--------+-------+-------------------------------+
      |    1   |   1   | :math:`I - KX(X^*KX)^{-1}X^*` |
      +--------+-------+-------------------------------+
      |    0   |   1   | error                         |
      +--------+-------+-------------------------------+

      Input/output:

         | :c:func:`primme_initialize` sets all of them to 0;
         | this field is written by :c:func:`primme_set_method` (see :ref:`methods`);
         | this field is read by :c:func:`dprimme`.

      See [r3]_ for a study about different projector configurations in JD.

   .. c:member:: PRIMME_INT ldOPs

      Recommended leading dimension to be used in |matrixMatvec|, |applyPreconditioner| and |massMatrixMatvec|.
      The default value is zero, which means no user recommendation. In that case,
      PRIMME computes ldOPs internally to get better memory performance.

      Input/output:

         | :c:func:`primme_initialize` sets this field to -1;
         | this field is read by :c:func:`dprimme`.

      .. versionadded:: 2.0

   .. c:member:: void (*monitorFun)(void *basisEvals, int *basisSize, int *basisFlags, int *iblock, int *blockSize, void *basisNorms, int *numConverged, void *lockedEvals, int *numLocked, int *lockedFlags, void *lockedNorms, int *inner_its, void *LSRes, const char *msg, double *time, primme_event *event, struct primme_params *primme, int *ierr)


      Convergence monitor. Used to customize how to report solver 
      information during execution (iteration number, matvecs, time, 
      unconverged and converged eigenvalues, residual norms, targets, etc).

      :param basisEvals:   array with approximate eigenvalues of the basis.
      :param basisSize:    size of the arrays, ``basisEvals``, ``basisFlags`` and ``basisNorms``.
      :param basisFlags:   state of every approximate pair in the basis.
      :param iblock:       indices of the approximate pairs in the block targeted during current iteration.
      :param blockSize:    size of array ``iblock``.
      :param basisNorms:   array with residual norms of the pairs in the basis.
      :param numConverged: number of pairs converged in the basis plus the number of the locked pairs (note that this value isn't monotonic).
      :param lockedEvals:  array with the locked eigenvalues.
      :param numLocked:    size of the arrays ``lockedEvals``, ``lockedFlags`` and ``lockedNorms``.
      :param lockedFlags:  state of each locked eigenpair.
      :param lockedNorms:  array with the residual norms of the locked pairs.
      :param inner_its:    number of performed QMR iterations in the current correction equation. It resets for each block vector.
      :param LSRes:        residual norm of the linear system at the current QMR iteration.
      :param msg:          output message or function name.
      :param time:         time duration.
      :param event:        event reported.
      :param primme:       parameters structure; the counter in ``stats`` are updated with the current number of matrix-vector products, iterations, elapsed time, etc., since start.
      :param ierr:         output error code; if it is set to non-zero, the current call to PRIMME will stop.

      This function is called at the following events:

      * ``*event == primme_event_outer_iteration``: every outer iterations.

        For this event the following inputs are provided:
        ``basisEvals``, ``basisNorms``, ``basisSize``, ``basisFlags``, ``iblock`` and ``blockSize``.

        ``basisNorms[iblock[i]]`` has the residual norm for the selected pair in the block.
        PRIMME avoids computing the residual of soft-locked pairs, ``basisNorms[i]`` for ``i<iblock[0]``.
        So those values may correspond to previous iterations. The values ``basisNorms[i]`` for ``i>iblock[blockSize-1]``
        are not valid.

        If |locking| is enabled, ``lockedEvals``, ``numLocked``, ``lockedFlags`` and ``lockedNorms`` are also provided.

        ``inner_its`` and  ``LSRes`` are not provided.

      * ``*event == primme_event_inner_iteration``: every QMR iteration.

        ``basisEvals[0]`` and ``basisNorms[0]`` provides the approximate eigenvalue and the residual norm
        of the pair which is improved in the current correction equation. If |convTest| is |primme_adaptive|
        or |primme_adaptive_ETolerance|, ``basisEvals[0]`` and ``basisNorms[0]`` are updated every QMR iteration.

        ``inner_its`` and  ``LSRes`` are also provided.

        ``lockedEvals``, ``numLocked``, ``lockedFlags`` and ``lockedNorms`` may not be provided.

      * ``*event == primme_event_converged``: a new eigenpair in the basis passed the convergence criterion.

        ``iblock[0]`` is the index of the newly converged pair in the basis which will be locked or soft-locked.
        The following are provided: ``basisEvals``, ``basisNorms``, ``basisSize``, ``basisFlags`` and ``blockSize[0]==1``.

        ``lockedEvals``, ``numLocked``, ``lockedFlags`` and ``lockedNorms`` may not be provided.

        ``inner_its`` and  ``LSRes`` are not provided.

      * ``*event == primme_event_locked``: new pair was added to the locked eigenvectors.

        ``lockedEvals``, ``numLocked``, ``lockedFlags`` and ``lockedNorms`` are provided.
        The last element of ``lockedEvals``, ``lockedFlags`` and ``lockedNorms`` corresponds
        to the recent locked pair.
 
        ``basisEvals``, ``numConverged``, ``basisFlags`` and ``basisNorms`` may not be provided.

        ``inner_its`` and ``LSRes`` are not provided.

      * ``*event == primme_event_message``: output message

        ``msg`` is the message to print.

        The rest of the arguments are not provided.

      The values of ``basisFlags`` and ``lockedFlags`` are:

      * ``0``: unconverged.
      * ``1``: internal use; only in ``basisFlags``.
      * ``2``: passed convergence test |convTestFun|.
      * ``3``: *practically converged* because the solver may not be able 
        to reduce the residual norm further without recombining 
        the locked eigenvectors.

      The actual type of ``basisEvals``, ``basisNorms``, ``lockedEvals``, ``lockedNorms`` and ``LSRes`` matches the type of ``evecs`` of the
      calling  :c:func:`dprimme` (or a variant), unless the user sets |monitorFun_type| to
      another precision.

      Input/output:

         | :c:func:`primme_initialize` sets this field to NULL;
         | :c:func:`dprimme` sets this field to an internal function if it is NULL;
         | this field is read by :c:func:`dprimme`.

      .. versionchanged:: 3.0

   .. c:member:: primme_op_datatype monitorFun_type

      Precision of the vectors ``basisEvals``, ``basisNorms``, ``lockedEvals``, ``lockedNorms`` and ``LSRes`` passed to |monitorFun|.

      If it is ``primme_op_default``, the vectors' type matches the calling
      :c:func:`dprimme` (or a variant). Otherwise, the precision is half,
      single, or double, if |monitorFun_type| is ``primme_half``, ``primme_float``
      or ``primme_double`` respectively.

      Input/output:

         | :c:func:`primme_initialize` sets this field to ``primme_op_default``;
         | this field is read by :c:func:`dprimme`, and if it is
           ``primme_op_default`` it is set to the value that matches the precision of
           calling function.

      .. versionadded:: 3.0

   .. c:member:: void *monitor

      This field may be used to pass any required information 
      to the function |monitorFun|.

      Input/output:

         | :c:func:`primme_initialize` sets this field to NULL;

      .. versionadded:: 2.0

   .. c:member:: PRIMME_INT stats.numOuterIterations

      Hold the number of outer iterations. The value is available during execution and at the end.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`dprimme`.

   .. c:member:: PRIMME_INT stats.numRestarts

      Hold the number of restarts during execution and at the end.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`dprimme`.

   .. c:member:: PRIMME_INT stats.numMatvecs

      Hold how many vectors the operator in |matrixMatvec| has been applied on.
      The value is available during execution and at the end.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`dprimme`.

   .. c:member:: PRIMME_INT stats.numPreconds

      Hold how many vectors the operator in |applyPreconditioner| has been applied on.
      The value is available during execution and at the end.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`dprimme`.

   .. c:member:: PRIMME_INT stats.numGlobalSum

      Hold how many times |globalSumReal| has been called.
      The value is available during execution and at the end.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`dprimme`.

   .. c:member:: double stats.volumeGlobalSum

      Hold how many :c:type:`REAL` have been reduced by |globalSumReal|.
      The value is available during execution and at the end.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`dprimme`.

   .. c:member:: PRIMME_INT stats.numBroadcast

      Hold how many times |broadcastReal| has been called.
      The value is available during execution and at the end.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`dprimme`.

      .. versionadded:: 3.0

   .. c:member:: double stats.volumeBroadcast

      Hold how many :c:type:`REAL` have been broadcast by |broadcastReal|.
      The value is available during execution and at the end.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`dprimme`.

      .. versionadded:: 3.0

   .. c:member:: PRIMME_INT stats.numOrthoInnerProds

      Hold how many inner products with vectors of length |nLocal| have been computed during orthogonalization.
      The value is available during execution and at the end.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`dprimme`.

      .. versionadded:: 3.0

   .. c:member:: double stats.elapsedTime

      Hold the wall clock time spent by the call to :c:func:`dprimme` or :c:func:`zprimme`.
      The value is available at the end of the execution.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`dprimme`.

   .. c:member:: double stats.timeMatvec

      Hold the wall clock time spent by |matrixMatvec|.
      The value is available at the end of the execution.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`dprimme`.

   .. c:member:: double stats.timePrecond

      Hold the wall clock time spent by |applyPreconditioner|.
      The value is available at the end of the execution.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`dprimme`.

   .. c:member:: double stats.timeOrtho

      Hold the wall clock time spent by orthogonalization.
      The value is available at the end of the execution.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`dprimme`.

   .. c:member:: double stats.timeGlobalSum

      Hold the wall clock time spent by |globalSumReal|.
      The value is available at the end of the execution.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`dprimme`.

   .. c:member:: double stats.timeBroadcast

      Hold the wall clock time spent by |broadcastReal|.
      The value is available at the end of the execution.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`dprimme`.

      .. versionadded:: 3.0

   .. c:member:: double stats.estimateMinEVal

      Hold the estimation of the smallest eigenvalue for the current eigenproblem.
      The value is available during execution and at the end.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`dprimme`.

   .. c:member:: double stats.estimateMaxEVal

      Hold the estimation of the largest eigenvalue for the current eigenproblem.
      The value is available during execution and at the end.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`dprimme`.

   .. c:member:: double stats.estimateLargestSVal

      Hold the estimation of the largest singular value (i.e., the absolute value of
      the eigenvalue with largest absolute value) for the current eigenproblem.
      The value is available during execution and at the end.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`dprimme`.

   .. c:member:: double stats.maxConvTol

      Hold the maximum residual norm of the converged eigenvectors.
      The value is available during execution and at the end.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`dprimme`.

   .. c:member:: PRIMME_INT stats.lockingIssue

      It is set to a nonzero value if some of the returned eigenpairs do not pass the convergence criterion.
      See |convTestFun| and |eps|.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`dprimme`.

      .. versionadded:: 3.0

   .. c:member:: void (*convTestFun) (double *eval, void *evec, double *resNorm, int *isconv, primme_params *primme, int *ierr)

      Function that evaluates if the approximate eigenpair has converged.
      If NULL, it is used the default convergence criteria (see |eps|).
   
      :param eval: the approximate value to evaluate.
      :param evec: one dimensional array of size |nLocal| containing the approximate vector; it can be NULL.
      :param resNorm: the norm of the residual vector.
      :param isconv: (output) the function sets zero if the pair is not converged and non zero otherwise.
      :param primme: parameters structure.
      :param ierr: output error code; if it is set to non-zero, the current call to PRIMME will stop.

      The actual type of ``evec`` matches the type of ``evecs`` of the
      calling  :c:func:`dprimme` (or a variant), unless the user sets |convTestFun_type| to
      another precision.

      Input/output:

         | :c:func:`primme_initialize` sets this field to NULL;
         | this field is read by :c:func:`dprimme`.

      .. versionadded:: 2.0

   .. c:member:: primme_op_datatype convTestFun_type

      Precision of the vectors ``evec`` passed to |convTestFun|.

      If it is ``primme_op_default``, ``evec``'s type matches the calling
      :c:func:`dprimme` (or a variant). Otherwise, the precision is half,
      single, or double, if |convTestFun_type| is ``primme_half``, ``primme_float``
      or ``primme_double`` respectively.

      Input/output:

         | :c:func:`primme_initialize` sets this field to ``primme_op_default``;
         | this field is read by :c:func:`dprimme`, and if it is
           ``primme_op_default`` it is set to the value that matches the precision of
           calling function.

      .. versionadded:: 3.0

   .. c:member:: void *convtest

      This field may be used to pass any required information 
      to the function |convTestFun|.

      Input/output:

         | :c:func:`primme_initialize` sets this field to NULL;

      .. versionadded:: 2.0

   .. c:member:: void *queue

      Pointer to the accelerator's data structure.

      If the main call is :c:func:`dprimme_magma` or a variant, this field
      should have the pointer to an initialized ``magma_queue_t``.

      See example :file:`examples/ex_eigs_dmagma.c`.

      Input/output:

         | :c:func:`primme_initialize` sets this field to NULL;
         | this field is read by :c:func:`dprimme_magma`.

      .. versionadded:: 3.0

.. _methods:

Preset Methods
--------------

.. c:enum:: primme_preset_method

   .. c:enumerator:: PRIMME_DEFAULT_MIN_TIME

      Set as |JDQMR_ETol| when |target| is either ``primme_smallest`` or
      ``primme_largest``, and as |JDQMR| otherwise. This method is usually
      the fastest if the cost of the matrix vector product is inexpensive.

   .. c:enumerator:: PRIMME_DEFAULT_MIN_MATVECS

      Currently set as |GD_Olsen_plusK|; this method usually performs
      fewer matrix vector products than other methods, so it's a good
      choice when this operation is expensive.

   .. c:enumerator:: PRIMME_DYNAMIC

      Switches to the best method dynamically; currently, between
      methods |DEFAULT_MIN_TIME| and |DEFAULT_MIN_MATVECS|.

      With |DYNAMIC| :c:func:`primme_set_method` sets
      |dynamicMethodSwitch| = 1 and makes the same changes as
      for method |DEFAULT_MIN_TIME|.

   .. c:enumerator:: PRIMME_Arnoldi

      Arnoldi implemented à la Generalized Davidson.

      With |Arnoldi| :c:func:`primme_set_method` sets:

      .. hlist::

         * |locking| = 0;
         * |maxPrevRetain| = 0;
         * |precondition| = 0;
         * |maxInnerIterations| = 0.

   .. c:enumerator:: PRIMME_GD

      Generalized Davidson.

      With |GD| :c:func:`primme_set_method` sets:

      .. hlist::

         * |locking| = 0;
         * |maxPrevRetain| = 0;
         * |robustShifts| = 1;
         * |maxInnerIterations| = 0;
         * |RightX| = 0;
         * |SkewX| = 0.

   .. c:enumerator:: PRIMME_GD_plusK

      GD with locally optimal restarting. 

      With |GD_plusK| :c:func:`primme_set_method` sets
      |maxPrevRetain| = 2 if |maxBlockSize| is 1 and |numEvals| > 1;
      otherwise it sets |maxPrevRetain| to |maxBlockSize|. Also:

      .. hlist::

         * |locking| = 0;
         * |maxInnerIterations| = 0;
         * |RightX| = 0;
         * |SkewX| = 0.

   .. c:enumerator:: PRIMME_GD_Olsen_plusK

      GD+k and the cheap Olsen's Method.

      With |GD_Olsen_plusK| :c:func:`primme_set_method` makes the
      same changes as for method |GD_plusK| and sets |RightX| = 1.

   .. c:enumerator:: PRIMME_JD_Olsen_plusK

      GD+k and Olsen's Method.

      With |JD_Olsen_plusK| :c:func:`primme_set_method` makes the
      same changes as for method |GD_plusK| and also sets
      |robustShifts| = 1, |RightX| to 1, and |SkewX| to 1.

   .. c:enumerator:: PRIMME_RQI

      (Accelerated) Rayleigh Quotient Iteration.

      With |RQI| :c:func:`primme_set_method` sets:

      .. hlist::

         * |locking| = 1;
         * |maxPrevRetain| = 0;
         * |robustShifts|  = 1;
         * |maxInnerIterations| = -1;
         * |LeftQ|   = 1;
         * |LeftX|   = 1;
         * |RightQ|  = 0;
         * |RightX|  = 1;
         * |SkewQ|   = 0;
         * |SkewX|   = 0;
         * |convTest| = |primme_full_LTolerance|.

      .. note::

         If |numTargetShifts| > 0 and |targetShifts| are provided, the interior problem
         solved uses these shifts in the correction equation. Therefore RQI becomes
         INVIT (inverse iteration) in that case.

   .. c:enumerator:: PRIMME_JDQR

      Jacobi-Davidson with fixed number of inner steps.

      With |JDQR| :c:func:`primme_set_method` sets:

      .. hlist::

         * |locking|     = 1;
         * |maxPrevRetain|      = 1;
         * |robustShifts|       = 0;
         * |maxInnerIterations| = 10 if it is 0;
         * |LeftQ|   = 0;
         * |LeftX|   = 1;
         * |RightQ|  = 1;
         * |RightX|  = 1;
         * |SkewQ|   = 1;
         * |SkewX|   = 1;
         * |relTolBase| = 1.5;
         * |convTest| = |primme_full_LTolerance|.

   .. c:enumerator:: PRIMME_JDQMR

      Jacobi-Davidson with adaptive stopping criterion for inner Quasi Minimum Residual (QMR).

      With |JDQMR| :c:func:`primme_set_method` sets:

      .. hlist::

         * |locking| = 0;
         * |maxPrevRetain| = 1 if it is 0
         * |maxInnerIterations| = -1;
         * |LeftQ|   = |precondition|;
         * |LeftX|   = 1;
         * |RightQ|  = 0;
         * |RightX|  = 0;
         * |SkewQ|   = 0;
         * |SkewX|   = 1;
         * |convTest|  = |primme_adaptive|.

   .. c:enumerator:: PRIMME_JDQMR_ETol

      JDQMR but QMR stops after residual norm reduces by a 0.1 factor.

      With |JDQMR_ETol| :c:func:`primme_set_method` makes the same
      changes as for the method |JDQMR| and sets
      |convTest| = |primme_adaptive_ETolerance|.

   .. c:enumerator:: PRIMME_STEEPEST_DESCENT

      Steepest descent.

      With |STEEPEST_DESCENT| :c:func:`primme_set_method` sets:

      .. hlist::

         * |locking|    = 1;
         * |maxBasisSize| = |numEvals| \* 2;
         * |minRestartSize| = |numEvals|;
         * |maxBlockSize| = |numEvals|;
         * |maxPrevRetain|      = 0;
         * |robustShifts|       = 0;
         * |maxInnerIterations| = 0;
         * |RightX|  = 1;
         * |SkewX|   = 0.

   .. c:enumerator:: PRIMME_LOBPCG_OrthoBasis

      LOBPCG with orthogonal basis.

      With |LOBPCG_OrthoBasis| :c:func:`primme_set_method` sets:

      .. hlist::

         * |locking|    = 0;
         * |maxBasisSize| = |numEvals| \* 3;
         * |minRestartSize| = |numEvals|;
         * |maxBlockSize| = |numEvals|;
         * |maxPrevRetain|      = |numEvals|;
         * |robustShifts|       = 0;
         * |maxInnerIterations| = 0;
         * |RightX|  = 1;
         * |SkewX|   = 0.

   .. c:enumerator:: PRIMME_LOBPCG_OrthoBasis_Window

      LOBPCG with sliding window of |maxBlockSize| < 3 \* |numEvals|.

      With |LOBPCG_OrthoBasis_Window| :c:func:`primme_set_method` sets:

      .. hlist::

         * |locking|    = 0;
         * |maxBasisSize| = |maxBlockSize| \* 3;
         * |minRestartSize| = |maxBlockSize|;
         * |maxPrevRetain|      = |maxBlockSize|;
         * |robustShifts|       = 0;
         * |maxInnerIterations| = 0;
         * |RightX|  = 1;
         * |SkewX|   = 0.

.. _error-codes:

Error Codes
-----------

The functions :c:func:`dprimme` and :c:func:`zprimme` return one of the following error codes.
Some of the error codes have a macro associated which is indicated in brackets.

*  0: success; usually all requested eigenpairs have converged.
* -1: (``PRIMME_UNEXPECTED_FAILURE``) unexpected internal error; please consider to set |printLevel| to a value larger than 0 to see the call stack and to report these errors because they may be bugs.
* -2: (``PRIMME_MALLOC_FAILURE``) failure in allocating memory; it can be either CPU or GPU.
* -3: (``PRIMME_MAIN_ITER_FAILURE``) maximum number of outer iterations |maxOuterIterations| or matvecs |maxMatvecs| reached.
* -4: if argument ``primme`` is NULL.
* -5: if |n| < 0 or |nLocal| < 0 or |nLocal| > |n|.
* -6: if |numProcs| < 1.
* -7: if |matrixMatvec| is NULL.
* -8: if |applyPreconditioner| is NULL and |precondition| > 0.
* -10: if |numEvals| > |n|.
* -11: if |numEvals| < 0.
* -12: if |convTestFun| is not NULL and |eps| > 0 and |eps| < machine precision given by |internalPrecision| and the precision of PRIMME call (:c:func:`sprimme`, :c:func:`dprimme`...).
* -13: if |target| is not properly defined.
* -14: if |target| is one of |primme_closest_geq|,
  |primme_closest_leq|, |primme_closest_abs| or |primme_largest_abs| but
  |numTargetShifts| <= 0 (no shifts).
* -15: if |target| is one of |primme_closest_geq|,
  |primme_closest_leq|, |primme_closest_abs| or |primme_largest_abs| but
  |targetShifts| is NULL  (no shifts array).
* -16: if |numOrthoConst| < 0 or |numOrthoConst| > |n|.
  (no free dimensions left).
* -17: if |maxBasisSize| < 2 and |n| > 2.
* -18: if |minRestartSize| < 0, or |minRestartSize| is zero but |n| > 2 and |numEvals| > 0.
* -19: if |maxBlockSize| < 0, or |maxBlockSize| is zero but |numEvals| > 0.
* -20: if |maxPrevRetain| < 0.
* -22: if |initSize| < 0.
* -23: if |locking| == 0 and |initSize| > |maxBasisSize|.
* -24: if |locking| and |initSize| > |numEvals|.
* -25: if |maxPrevRetain| + |minRestartSize| >= |maxBasisSize|, and |n| > |maxBasisSize|.
* -26: if |minRestartSize| >= |n|, and |n| > 2.
* -27: if |printLevel| < 0 or |printLevel| > 5.
* -28: if |convTest| is not one of
  |primme_full_LTolerance|, |primme_decreasing_LTolerance|,
  |primme_adaptive_ETolerance| or |primme_adaptive|.
* -29: if |convTest| == |primme_decreasing_LTolerance| and |relTolBase| <= 1.
* -30: if ``evals`` is NULL.
* -31: if ``evecs`` is NULL, or is not a GPU pointer when calling a GPU variant (for instance ``magma_dprimme``).
* -32: if ``resNorms`` is NULL.
* -33: if |locking| == 0 and |minRestartSize| < |numEvals| and |n| > 2.
* -34: if |ldevecs| < |nLocal|.
* -35: if |ldOPs| is not zero and less than |nLocal|.
* -38: if |locking| == 0 and |target| is |primme_closest_leq| or |primme_closest_geq|.
* -40: (``PRIMME_LAPACK_FAILURE``) some LAPACK function performing a factorization returned an error code; set |printLevel| > 0 to see the error code and the call stack.
* -41: (``PRIMME_USER_FAILURE``) some of the user-defined functions (|matrixMatvec|, |applyPreconditioner|, ...) returned a non-zero error code; set |printLevel| > 0 to see the call stack that produced the error.
* -42: (``PRIMME_ORTHO_CONST_FAILURE``) the provided orthogonal constraints (see |numOrthoConst|) are not full rank.
* -43: (``PRIMME_PARALLEL_FAILURE``) some process has a different value in an input option than the process zero, or it is not acting coherently; set |printLevel| > 0 to see the call stack that produced the error.
* -44: (``PRIMME_FUNCTION_UNAVAILABLE``) PRIMME was not compiled with support for the requesting precision or for GPUs.

.. include:: epilog.inc
