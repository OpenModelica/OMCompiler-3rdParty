/* Copyright (C) 2007 VRTech Industrial Technologies - www.vrtech.com.br.
 * Copyright (C) 2007 Tong Kewei, Beihang University, - www.buaa.edu.cn.
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 */

package org.coinor;

import java.io.File;

/** A Java Native Interface for the Ipopt optimization solver.
 *
 * Ipopt is a solver for large scale nonlinear optimization problems (NLP).
 *
 * The Java Native Interface (JNI) is a programming framework that allows
 * Java code running in the Java Virtual Machine (JVM) to call and be
 * called by native applications (programs specific to a hardware and
 * operating system platform) and libraries written in other languages,
 * such as C and C++.
 *
 * This class is a JNI hook around the C++ interface of Ipopt, as a consequence
 * it will need a nativelly compiled DLL to run.
 * For more details about Ipopt [click here](https://github.com/coin-or/Ipopt).
 *
 * The user should subclass this class and implement the abstract methods.
 * At some point before solving the problem the
 * {@link #create(int, int, int, int, int)}
 * function should be called.
 * For simple cases you can call this function in the constructor of your class.
 *
 * Once the problem was created, {@link #OptimizeNLP()} will solve the problem.
 * Objects of this class can be reused to solve different problems, in other words,
 * {@link #create(int, int, int, int, int)}
 * and {@link #OptimizeNLP()} can be called multiple times.
 *
 * Programmers must call {@link #dispose()} when finished using a
 * Ipopt object, otherwise the nativelly allocated memory will be disposed of only
 * when the JVM call {@link #finalize()} on it.
 *
 * @author Rafael de Pelegrini Soares
 * @author Edson C. do Valle
 * @author Tong Kewei, BeiHang University
 */
public abstract class Ipopt
{
   /* Native function should not be used directly */
   private native boolean AddIpoptIntOption(
      long   ipopt,
      String keyword,
      int    val
   );

   /* Native function should not be used directly */
   private native boolean AddIpoptNumOption(
      long   ipopt,
      String keyword,
      double val
   );

   /* Native function should not be used directly */
   private native boolean AddIpoptStrOption(
      long   ipopt,
      String keyword,
      String val
   );

   /* Native function should not be used directly */
   private native long CreateIpoptProblem(
      int n,
      int m,
      int nele_jac,
      int nele_hess,
      int index_style
   );

   /* Native function should not be used directly */
   private native void FreeIpoptProblem(
      long ipopt
   );

   /* Native function should not be used directly */
   private native int OptimizeTNLP(
      long   ipopt,
      double x[],
      double g[],
      double obj_val[],
      double mult_g[],
      double mult_x_L[],
      double mult_x_U[],
      double callback_grad_f[],
      double callback_jac_g[],
      double callback_hess[]
   );

   /** Use C index style for iRow and jCol vectors */
   public final static int C_STYLE = 0;

   /** Use FORTRAN index style for iRow and jCol vectors */
   public final static int FORTRAN_STYLE = 1;

   /** The possible Ipopt status return codes: should be kept in sync with Ipopt return codes */
   public final static int SOLVE_SUCCEEDED = 0;
   public final static int ACCEPTABLE_LEVEL = 1;
   public final static int INFEASIBLE_PROBLEM = 2;
   public final static int SEARCH_DIRECTION_TOO_SMALL = 3;
   public final static int DIVERGING_ITERATES = 4;
   public final static int USER_REQUESTED_STOP = 5;
   public final static int ITERATION_EXCEEDED = -1;
   public final static int RESTORATION_FAILED = -2;
   public final static int ERROR_IN_STEP_COMPUTATION = -3;
   public final static int CPUTIME_EXCEEDED = -4;
   public final static int NOT_ENOUGH_DEGREES_OF_FRE = -10;
   public final static int INVALID_PROBLEM_DEFINITION = -11;
   public final static int INVALID_OPTION = -12;
   public final static int INVALID_NUMBER_DETECTED = -13;
   public final static int UNRECOVERABLE_EXCEPTION = -100;
   public final static int NON_IPOPT_EXCEPTION = -101;
   public final static int INSUFFICIENT_MEMORY = -102;
   public final static int INTERNAL_ERROR = -199;

   /** Pointer to the native optimization object */
   private long ipopt;

   /// Callback arguments
   private double callback_grad_f[];
   private double callback_jac_g[];
   private double callback_hess[];

   /** Final value of variable values */
   private double x[];

   /** Final value of objective function */
   private double obj_val[] = {0.0};

   /** Values of constraint at final point */
   private double g[];

   /** Final multipliers for lower variable bounds */
   private double mult_x_L[];

   /** Final multipliers for upper variable bounds */
   private double mult_x_U[];

   /** Final multipliers for constraints */
   private double mult_g[];

   /** Status returned by the solver */
   private int status = INVALID_PROBLEM_DEFINITION;

   /** Creates a new NLP Solver using a default as the DLL name.
    *
    * This expects the the Ipopt DLL can somehow be found
    * and that it has the canoncial name "ipopt" (on Unix, et.al.)
    * or "ipopt-3" or "ipopt-0" (on Windows).
    *
    * @see #Ipopt()
    */
   public Ipopt()
   {
      if( System.getProperty("os.name").toLowerCase().indexOf("win") >= 0 )
      {
         /* for Ipopt releases, it should be ipopt-3.dll
          * for other intermediate versions, it should be ipopt-0.dll
          * with MinGW, libtool adds a "lib" prefix
          * finally, try also without version info
          */
         final String[] candidates = { "ipopt-3", "ipopt-0", "libipopt-3", "libipopt-0", "ipopt", "libipopt" };
         boolean loadedlib = false;
         for( String c : candidates )
         {
            try
            {
               System.loadLibrary(c);
               loadedlib = true;
               break;
            }
            catch( UnsatisfiedLinkError e )
            { }
         }
         if( !loadedlib )
         {
            throw new UnsatisfiedLinkError("Could not load Ipopt library. Check your java.library.path.");
         }
      }
      else
      {
         System.loadLibrary("ipopt");
      }
   }

   /** Creates a NLP Solver for the given DLL file.
    * The given file must implement the native interface required by this class.
    * The given file must be located in some library search path.
    *
    * @param DLL the name of the DLL (without the extension or any platform dependent prefix).
    *
    * @see #Ipopt()
    */
   public Ipopt(
      String DLL)
   {
      // Loads the library
      System.loadLibrary(DLL);
   }

   /** Creates a NLP Solver for the given DLL file and path.
    * The given file must implement the native interface required by this class.
    *
    * @param path the path where the DLL is found.
    * @param DLL the name of the DLL (without the extension or any platform dependent prefix).
    *
    * @see #Ipopt()
    */
   public Ipopt(
      String path,
      String DLL)
   {
      // Loads the library
      File file = new File(path, System.mapLibraryName(DLL));
      System.load(file.getAbsolutePath());
   }


   /** Method to request bounds on the variables and constraints.
    *
    * The values of n and m that were specified in create() and are passed
    * here for debug checking. Setting a lower bound to a value less than or
    * equal to the value of the option "nlp_lower_bound_inf"
    * will cause Ipopt to assume no lower bound. Likewise, specifying the upper bound above or
    * equal to the value of the option nlp_upper_bound_inf
    * will cause Ipopt to assume no upper bound. These options are set to -10<sup>19</sup> and
    * 10<sup>19</sup>, respectively, by default, but may be modified by changing these
    * options.
    *
    *  @param n   (in) the number of variablesin the problem
    *  @param x_l (out) the lower bounds for the variables
    *  @param x_u (out) the upper bounds for the variables
    *  @param m   (in) the number of constraints in the problem
    *  @param g_l (out) the lower bounds for the constraints
    *  @param g_u (out) the upper bounds for the constraints
    *
    * @return true on success, otherwise false
    */
   abstract protected boolean get_bounds_info(
      int      n,
      double[] x_l,
      double[] x_u,
      int      m,
      double[] g_l,
      double[] g_u
   );

   /** Method to request the starting point before iterating.
    *
    *  The boolean variables indicate whether the algorithm requires to
    *  have x, z_L/z_u, and lambda initialized, respectively.  If, for some
    *  reason, the algorithm requires initializations that cannot be
    *  provided, false should be returned and Ipopt will stop.
    *  The default options only require initial values for the primal
    *  variables.
    *
    *  Note, that the initial values for bound multiplier components for
    *  absent bounds are ignored.
    *
    *  @param n      (in) the number of variables in the problem; it will have the same value that was specified in create()
    *  @param init_x (in) if true, this method must provide an initial value for the primal variables
    *  @param x      (out) the initial values for the primal variables
    *  @param init_z (in) if true, this method must provide an initial value for the bound multipliers
    *  @param z_L    (out) the initial values for the lower bound multipliers
    *  @param z_U    (out) the initial values for the upper bound multipliers
    *  @param m      (in) the number of constraints in the problem; it will have the same value that was specified in create()
    *  @param init_lambda (in) if true, this method must provide an initial value for the constraint multipliers
    *  @param lambda (out) the initial values for the constraint multipliers
    *
    * @return true on success, otherwise false
    */
   abstract protected boolean get_starting_point(
      int      n,
      boolean  init_x,
      double[] x,
      boolean  init_z,
      double[] z_L,
      double[] z_U,
      int      m,
      boolean  init_lambda,
      double[] lambda
   );

   /** Method to request the value of the objective function.
    *
    *  @param n     (in) the number of variables in the problem; it will have the same value that was specified in create()
    *  @param x     (in) the values for the primal variables at which the objective function is to be evaluated
    *  @param new_x (in) false if any evaluation method (`eval_*`) was previously called with the same values in x, true otherwise.
    *                    This can be helpful when users have efficient implementations that calculate multiple outputs at once.
    *                    Ipopt internally caches results from the TNLP and generally, this flag can be ignored.
    *  @param obj_value (out) storage for the value of the objective function
    *
    * @return true on success, otherwise false
    */
   abstract protected boolean eval_f(
      int      n,
      double[] x,
      boolean  new_x,
      double[] obj_value
   );

   /** Method to request the gradient of the objective function.
    *
    *  @param n     (in) the number of variables in the problem; it will have the same value that was specified in create()
    *  @param x     (in) the values for the primal variables at which the gradient is to be evaluated
    *  @param new_x (in) false if any evaluation method (`eval_*`) was previously called with the same values in x, true otherwise; see also eval_f()
    *  @param grad_f (out) array to store values of the gradient of the objective function.
    *                      The gradient array is in the same order as the variables
    *                      (i.e., the gradient of the objective with respect to `x[2]` should be put in `grad_f[2]`).
    *
    * @return true on success, otherwise false
    */
   abstract protected boolean eval_grad_f(
      int      n,
      double[] x,
      boolean  new_x,
      double[] grad_f
   );

   /** Method to request the constraint values.
    *
    *  @param n     (in) the number of variables in the problem; it will have the same value that was specified in create()
    *  @param x     (in) the values for the primal variables at which the constraint functions are to be evaluated
    *  @param new_x (in) false if any evaluation method (`eval_*`) was previously called with the same values in x, true otherwise; see also eval_f()
    *  @param m     (in) the number of constraints in the problem; it will have the same value that was specified in create()
    *  @param g     (out) array to store constraint function values, do not add or subtract the bound values.
    *
    * @return true on success, otherwise false
    */
   abstract protected boolean eval_g(
      int      n,
      double[] x,
      boolean  new_x,
      int      m,
      double[] g
   );

   /** Method to request either the sparsity structure or the values of the Jacobian of the constraints.
    *
    * The Jacobian is the matrix of derivatives where the derivative of
    * the i-th constraint function with respect to the j-th variable is placed in row
    * i and column j.
    *
    * The arrays iRow and jCol only need to be filled once.
    * If the iRow and jCol arguments are not NULL (first call to this function),
    * then Ipopt expects that the sparsity structure of the Jacobian
    * (the row and column indices only) are written into iRow and jCol.
    * At this call, the arguments x and values will be NULL.
    * If the arguments x and values are not NULL, then Ipopt
    * expects that the value of the Jacobian as calculated from array x
    * is stored in array values (using the same order as used when
    * specifying the sparsity structure).
    * At this call, the arguments iRow and jCol will be NULL.
    *
    *  @param n     (in) the number of variables in the problem; it will have the same value that was specified in create()
    *  @param x     (in) first call: NULL; later calls: the values for the primal variables at which the constraint Jacobian is to be evaluated
    *  @param new_x (in) false if any evaluation method (`eval_*`) was previously called with the same values in x, true otherwise; see also eval_f()
    *  @param m     (in) the number of constraints in the problem; it will have the same value that was specified in create()
    *  @param nele_jac (in) the number of nonzero elements in the Jacobian; it will have the same value that was specified in create()
    *  @param iRow  (out) first call: array of length nele_jac to store the row indices of entries in the Jacobian of the constraints; later calls: NULL
    *  @param jCol  (out) first call: array of length nele_jac to store the column indices of entries in the Jacobian of the constraints; later calls: NULL
    *  @param values (out) first call: NULL; later calls: array of length nele_jac to store the values of the entries in the Jacobian of the constraints
    *
    * @return true on success, otherwise false
    */
   abstract protected boolean eval_jac_g(
      int      n,
      double[] x,
      boolean  new_x,
      int      m,
      int      nele_jac,
      int[]    iRow,
      int[]    jCol,
      double[] values
   );


   /** Method to request either the sparsity structure or the values of the Hessian of the Lagrangian.
    *
    * The Hessian matrix that Ipopt uses is the sum of the Hessian matrices of objective function (multiplied by obj_factor)
    * and each constraint function (multiplied by lambda).
    *
    * The arrays iRow and jCol only need to be filled once.
    * If the iRow and jCol arguments are not NULL (first call to this function),
    * then Ipopt expects that the sparsity structure of the Hessian
    * (the row and column indices only) are written into iRow and jCol.
    * At this call, the arguments x, lambda, and values will be NULL.
    * If the arguments x, lambda, and values are not NULL, then Ipopt
    * expects that the value of the Hessian as calculated from arrays x
    * and lambda are stored in array values (using the same order as
    * used when specifying the sparsity structure).
    * At this call, the arguments iRow and jCol will be NULL.
    *
    * As this matrix is symmetric, Ipopt expects that only the lower diagonal entries are specified.
    *
    *  @param n     (in) the number of variables in the problem; it will have the same value that was specified in create()
    *  @param x     (in) first call: NULL; later calls: the values for the primal variables at which the Hessian is to be evaluated
    *  @param new_x (in) false if any evaluation method (`eval_*`) was previously called with the same values in x, true otherwise; see also eval_f()
    *  @param obj_factor (in) factor in front of the objective term in the Hessian
    *  @param m     (in) the number of constraints in the problem; it will have the same value that was specified in create()
    *  @param lambda (in) the values for the constraint multipliers at which the Hessian is to be evaluated
    *  @param new_lambda (in) false if any evaluation method was previously called with the same values in lambda, true otherwise
    *  @param nele_hess (in) the number of nonzero elements in the Hessian; it will have the same value that was specified in create()
    *  @param iRow  (out) first call: array of length nele_hess to store the row indices of entries in the Hessian; later calls: NULL
    *  @param jCol  (out) first call: array of length nele_hess to store the column indices of entries in the Hessian; later calls: NULL
    *  @param values (out) first call: NULL; later calls: array of length nele_hess to store the values of the entries in the Hessian
    *
    * @return true on success, otherwise false
    */
   abstract protected boolean eval_h(
      int      n,
      double[] x,
      boolean  new_x,
      double   obj_factor,
      int      m,
      double[] lambda,
      boolean  new_lambda,
      int      nele_hess,
      int[]    iRow,
      int[]    jCol,
      double[] values
   );

   /** Dispose of the natively allocated memory.
    *
    * Programmers must call the dispose method when finished
    * using a Ipopt object.
    *
    * An JIpopt object can be reused to solve different problems by calling again
    * {@link #create(int, int, int, int, int)}.
    * In this case, you should call the dispose method only when you
    * finished with the object and it is not needed anymore.
    */
   public void dispose()
   {
      // dispose the native implementation
      if( ipopt != 0 )
      {
         FreeIpoptProblem(ipopt);
         ipopt = 0;
      }
   }

   @Deprecated
   protected void finalize() throws Throwable
   {
      dispose();
   }

   /** Create a new problem.
    *
    * This is get_nlp_info in the C++ interface.
    *
    * @param n the number of variables in the problem.
    * @param m the number of constraints in the problem.
    * @param nele_jac the number of nonzero entries in the Jacobian.
    * @param nele_hess the number of nonzero entries in the Hessian.
    * @param index_style the numbering style used for row/col entries in the sparse matrix format (C_STYLE or FORTRAN_STYLE).
    *
    * @return true on success, otherwise false
    */
   public boolean create(
      int n,
      int m,
      int nele_jac,
      int nele_hess,
      int index_style)
   {
      // delete any previously created native memory
      dispose();

      x = new double[n];
      g = new double[m];

      // allocate the callback arguments
      callback_grad_f = new double[n];
      callback_jac_g  = new double[nele_jac];
      callback_hess   = new double[nele_hess];

      // the multiplier
      mult_x_U = new double[n];
      mult_x_L = new double[n];
      mult_g   = new double[m];

      // create the optimization problem and return a pointer to it
      ipopt = CreateIpoptProblem(n, m,  nele_jac, nele_hess, index_style);

      //System.out.println("Finish Java Obj");
      return ipopt == 0 ? false : true;
   }

   /** Function for setting an integer option.
    *
    * For a list of valid keywords check the Ipopt documentation.
    *
    * @param keyword the option keyword
    * @param val the value
    * @return false if the option could not be set (e.g., if keyword is unknown)
    */
   public boolean setIntegerOption(
      String keyword,
      int    val)
   {
      if( ipopt == 0 )
      {
         return false;
      }

      return AddIpoptIntOption(ipopt, keyword, val);
   }

   /** Function for setting a number option.
    *
    * For a list of valid keywords check the Ipopt documentation.
    *
    * @param keyword the option keyword
    * @param val the value
    * @return false if the option could not be set (e.g., if keyword is unknown)
    */
   public boolean setNumericOption(
      String keyword,
      double val)
   {
      if( ipopt == 0 )
      {
         return false;
      }

      return AddIpoptNumOption(ipopt, keyword, val);
   }

   /** Function for setting a string option.
    *
    * For a list of valid keywords check the Ipopt documentation.
    *
    * @param keyword the option keyword
    * @param val the value
    * @return false if the option could not be set (e.g., if keyword is unknown)
    */
   public boolean setStringOption(
      String keyword,
      String val)
   {
      if( ipopt == 0 )
      {
         return false;
      }

      return AddIpoptStrOption(ipopt, keyword, val.toLowerCase());
   }

   /** This function actually solve the problem.
    *
    * The solve status returned is one of the constant fields of this class,
    * e.g. SOLVE_SUCCEEDED. For more details about the valid solve status
    * check the Ipopt documentation.
    *
    * @return the solve status
    *
    * @see #getStatus()
    */
   public int OptimizeNLP()
   {
      this.status = this.OptimizeTNLP(ipopt,
                                      x, g, obj_val, mult_g, mult_x_L, mult_x_U,
                                      callback_grad_f, callback_jac_g, callback_hess);

      return this.status;
   }

   /** Gives primal variable values at final point.
    * @return the primal variable values at the final point.
    */
   public double[] getVariableValues()
   {
      return x;
   }

   /** Gives objective function value at final point.
    * @return the final value of the objective function.
    */
   public double getObjectiveValue()
   {
      return obj_val[0];
   }

   /** Gives Ipopt status of last OptimizeNLP call.
    * @return the status of the solver.
    *
    * @see #OptimizeNLP()
    */
   public int getStatus()
   {
      return status;
   }

   /** Gives constraint function values at final point.
    * @return Returns the final values for the constraints functions.
    */
   public double[] getConstraintValues()
   {
      return g;
   }

   /** Gives constraint dual multipliers in final point.
    * @return Returns the final multipliers for the constraints.
    */
   public double[] getConstraintMultipliers()
   {
      return mult_g;
   }

   /** Gives dual multipliers for variable lower bounds in final point.
    * @return Returns the final multipliers for the variable lower bounds.
    */
   public double[] getLowerBoundMultipliers()
   {
      return mult_x_L;
   }

   /** Gives dual multipliers for variable upper bounds in final point.
    * @return Returns the final multipliers for the variable upper bounds.
    */
   public double[] getUpperBoundMultipliers()
   {
      return mult_x_U;
   }

   /** If you using_scaling_parameters = true, this method should be overloaded.
    *
    * To instruct IPOPT to use scaling values for variables, the first element of use_x_g_scaling should be set.
    * To instruct IPOPT to use scaling values for constraints, the second element of use_x_g_scaling should be set.
    *
    * @param obj_scaling  double[1] to store a scaling factor for the objective (negative value leads to maximizing the objective function)
    * @param n  the number of variables in the problem
    * @param x_scaling  array to store the scaling factors for the variables
    * @param m  the number of constraints in the problem
    * @param g_scaling  array to store the scaling factors for the constraints
    * @param use_x_g_scaling boolean[2] to store whether scaling factors for variables (1st entry) and constraints (2nd entry) should be used
    *
    * @return true on success, otherwise false
    */
   public boolean get_scaling_parameters(
      double[]  obj_scaling,
      int       n,
      double[]  x_scaling,
      int       m,
      double[]  g_scaling,
      boolean[] use_x_g_scaling)
   {
      return false;
   }

   /** When LBFGS hessian approximation is used, this method should be overloaded.
    *
    * @return number of nonlinear variables, a negative value indicates that all variables are negative
    */
   public int get_number_of_nonlinear_variables()
   {
      return -1;
   }

   /** When LBFGS hessian approximation is used, this method should be overloaded.
    *
    * @param num_nonlin_vars number of nonlinear variables and length of pos_nonlin_vars array
    * @param pos_nonlin_vars the indices of all nonlinear variables
    *
    * @return true on success, otherwise false
    */
   public boolean get_list_of_nonlinear_variables(
      int   num_nonlin_vars,
      int[] pos_nonlin_vars)
   {
      return false;
   }
}
