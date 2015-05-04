/* dheatilu.f -- translated by f2c (version 20100827).
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../solver/ddaskr_types.h"

/* Table of constant values */

static integer c__2594 = 2594;
static integer c__4468 = 4468;
static integer c__1 = 1;

/* ***BEGIN PROLOGUE  DHEATILU */
/* ***REFER TO  DDASKR */
/* ***DATE WRITTEN   020813   (YYMMDD) */
/* ***REVISION DATE */

/* ***DESCRIPTON */

/* ----------------------------------------------------------------------- */
/* Example program for DDASKR. */
/* DAE system derived from the discretized heat equation on a square. */

/* This is the double precision version. */
/* ----------------------------------------------------------------------- */

/* This program solves a DAE system that arises from the heat equation, */
/*   du/dt = u   + u */
/*            xx    yy */
/* posed on the 2-D unit square with zero Dirichlet boundary conditions. */
/* An M+2 by M+2 mesh is set on the square, with uniform spacing 1/(M+1). */
/* The spatial deriviatives are represented by standard central finite */
/* difference approximations.  At each interior point of the mesh, */
/* the discretized PDE becomes an ODE for the discrete value of u. */
/* At each point on the boundary, we pose the equation u = 0.  The */
/* discrete values of u form a vector U, ordered first by x, then by y. */
/* The result is a DAE system G(t,U,U') = 0 of size NEQ = (M+2)*(M+2). */

/* Initial conditions are posed as u = 16x(1-x)y(1-y) at t = 0. */
/* The problem is solved by DDASKR on the time interval t .le. 10.24. */

/* The root functions are R1(U) = max(u) - 0.1, R2(U) = max(u) - 0.01. */

/* The Krylov linear system solution method, with preconditioning, is */
/* selected.  The preconditioner is a sparse matrix with half-bandwidths */
/* equal to 1, i.e. a tridiagonal matrix.  (The true half-bandwidths */
/* are equal to M+2.)  This corresponds to ignoring the y-direction */
/* coupling in the ODEs, for purposes of preconditioning.  The extra */
/* iterations resulting from this approximation are offset by the lower */
/* storage and linear system solution costs for a tridiagonal matrix. */

/* The routines DJACILU and DPSOLILU that generate and solve the sparse */
/* preconditioner are provided in a separate file for general use. */

/* The output times are t = .01 * 2**n (n = 0,...,10).  The maximum of */
/* abs(u) over the mesh, and various performance statistics, are printed. */

/* For details and test results on this problem, see the reference: */
/*   Peter N. Brown, Alan C. Hindmarsh, and Linda R. Petzold, */
/*   Using Krylov Methods in the Solution of Large-Scale Differential- */
/*   Algebraic Systems, SIAM J. Sci. Comput., 15 (1994), pp. 1467-1488. */
/* ----------------------------------------------------------------------- */

/* ***ROUTINES CALLED */
/*   UINIT, DDASKR */

/* ***END PROLOGUE  DHEATILU */

/* Here are necessary declarations.  The dimension statements use a */
/* maximum value for the mesh parameter M. */

/* Main program */ int main(void)
{
    /* Format strings */
    static char fmt_15[] = " Error return from DSPSETUP: IERR = %5d";
    static char fmt_30[] = " DHEATILU: Heat Equation Example Program for DDASKR\n\n"
    	"    M+2 by M+2 mesh, M =%3d,  System size NEQ =%3d\n\n"
    	"    Root functions are: R1 = max(u) - 0.1, and R2 = max(u) - 0.01\n\n"
    	"    Linear solver method flag INFO(12) =%3d    (0 = direct, 1 = Krylov)\n"
	    "    Preconditioner is a banded approximation with ML =%3d  MU =%3d\n\n"
    	"    Incomplete factorization option =%3d    (1 = ILUT, 2 = ILUTP)\n"
    	"    Tolerances are RTOL =%10.1E   ATOL =%10.1E\n\n";
    static char fmt_40[] = "     t           UMAX\t        NQ      H        "
    		"  STEPS   NNI     NLI\n";
    static char fmt_60[] = "    %10.4E\t%10.3E     %d    %10.2E\t%5d\t %5d\t %5d\n";
    static char fmt_61[] = "\t\t    *****   Root found, JROOT = %d  %d\n";
    static char fmt_65[] = "\n   Final time reached =  %12.4E\n";
    static char fmt_90[] = "\n Final statistics for this run..\n"
	    "   RWORK size =%6d\tIWORK size =%6d\n"
        "   Number of time steps ................ =%8d\n"
        "   Number of residual evaluations ...... =%8d\n"
    	"   Number of res. evals. for precond.    =%8d\n"
    	"   Number of root function evaluations . =%8d\n"
    	"   Number of preconditioner evaluations  =%8d\n"
    	"   Number of preconditioner solves ..... =%8d\n"
	    "   Number of nonlinear iterations ...... =%8d\n"
    	"   Number of linear iterations ......... =%8d\n"
    	"   Average Krylov subspace dimension ... =%8.4f\n"
    	"  %5d nonlinear conv. failures,  %5d linear conv. failures\n";
    static char fmt_100[] = " Minimum lengths for work arrays WP and IWP: %7d %7d\n";


    /* System generated locals */
    integer i__1, i__2;
    real_number d__1, d__2, d__3;

    /* Local variables */
    extern /* Subroutine */ int _daskr_dpsolilu_();
    extern /* Subroutine */ int _daskr_dspsetup_(integer *, integer *, integer *,
	    real_number *, integer *, integer *, integer *, integer *);
    static integer i__, m;
    static real_number t, u[144];
    static integer ml;
    static real_number dx, hu;
    static integer mu, nli, neq, nni, npe, nre, liw, nps, nrt, lrw, nqu, nst, 
	    idid, ncfl, ncfn, info[20], ipar[34];
    static real_number atol;
    extern /* Subroutine */ int resh_();
    static integer ierr;
    static real_number rpar[4];
    static integer nrte;
    static real_number umax, rtol;
    static integer iout, lout, nout;
    static real_number tout;
    static integer mband;
    static real_number coeff, avdim;
    extern /* Subroutine */ int uinit_(real_number *, real_number *, real_number
	    *, integer *);
    static integer iwork[4508], jroot[2];
    static real_number rwork[5293];
    extern /* Subroutine */ int _daskr_ddaskr_(Unknown_fp, integer *, real_number *,
	    real_number *, real_number *, real_number *, integer *, real_number *,
	     real_number *, integer *, real_number *, integer *, integer *,
	    integer *, real_number *, integer *, Unknown_fp, Unknown_fp, Unknown_fp, integer *,
	    integer *);
    extern /* Subroutine */ int rtheat_();
    static real_number uprime[144];
    static integer lwpmin;
    extern /* Subroutine */ int _daskr_djacilu_();
    static integer liwpmin;


/* Set LOUT, the unit number of the output device. */
    lout = 6;

/* Here set parameters for the problem being solved.  Use RPAR and IPAR */
/* to communicate these to the other routines. */

    m = 10;
    dx = 1. / (m + 1);
    neq = (m + 2) * (m + 2);
    coeff = 1. / (dx * dx);
    ipar[32] = neq;
    ipar[33] = m;
    rpar[2] = dx;
    rpar[3] = coeff;

/* Set NRT = number of root functions. */
    nrt = 2;

/* Here set the lengths of the preconditioner work arrays WP and IWP, */
/* load them into IWORK, and set the total lengths of WORK and IWORK. */
    iwork[26] = 2594;
    iwork[27] = 4468;
    lrw = 5293;
    liw = 4508;
/* Load values into IPAR and RPAR for sparse preconditioner. */
    ml = 1;
    mu = 1;
    ipar[0] = ml;
    ipar[1] = mu;
    ipar[2] = 5;
    ipar[3] = 5;
    ipar[4] = 1;
    ipar[5] = 5;
    ipar[6] = 1;
    ipar[7] = 1;
    ipar[8] = 2;
    ipar[9] = 0;
    ipar[10] = 1;
    ipar[29] = 0;
    rpar[0] = .001;
    rpar[1] = .01;
/* Check IPAR, RPAR, LENWP and LENIWP for illegal entries and long */
/* enough work array lengths. */
    _daskr_dspsetup_(&neq, &c__2594, &c__4468, rpar, ipar, &ierr, &lwpmin, &liwpmin);
    if (ierr != 0) {
    	printf(fmt_15, ierr);
		if (lwpmin > 2594) {
			puts(" More WP work array length needed");
		}
		if (liwpmin > 4468) {
			puts(" More IWP work array length needed");
		}
		exit(0);
    }

/* Call subroutine UINIT to initialize U and UPRIME. */

    uinit_(u, uprime, rpar, ipar);

/* ----------------------------------------------------------------------- */
/* Here we set up the INFO array, which describes the various options */
/* in the way we want DDASKR to solve the problem. */
/* In this case, we select the iterative preconditioned Krylov method, */
/* and we supply the sparse preconditioner routines DJACILU/DPSOLILU. */

/* We first initialize the entire INFO array to zero, then set select */
/* entries to nonzero values for desired solution options. */

/* To select the Krylov iterative method for the linear systems, */
/* we set INFO(12) = 1. */

/* Since we are using a preconditioner that involves approximate */
/* Jacobian elements requiring preprocessing, we have a JAC routine, */
/* namely subroutine DJACILU, and we must set INFO(15) = 1 to indicate */
/* this to DDASKR. */

/* No other entries of INFO need to be changed for this example. */
/* ----------------------------------------------------------------------- */

    for (i__ = 1; i__ <= 20; ++i__) {
/* L20: */
	info[i__ - 1] = 0;
    }

    info[11] = 1;
    info[14] = 1;

/* Here we set tolerances for DDASKR to indicate how much accuracy */
/* we want in the solution, in the sense of local error control. */
/* For this example, we ask for pure absolute error control with a */
/* tolerance of 1.0D-5. */
    rtol = 0.;
    atol = 1e-5;

/* Here we generate a heading with important parameter values. */

    printf(fmt_30, m, neq, info[11], ml, mu, c__1, rtol, atol);
    printf("%s", fmt_40);

/* ----------------------------------------------------------------------- */
/* Now we solve the problem. */

/* DDASKR will be called to compute 11 intermediate solutions from */
/* tout = 0.01 to tout = 10.24 by powers of 2. */

/* We pass to DDASKR the names DJACILU and DPSOLILU for the JAC and PSOL */
/* routines to do the preconditioning. */

/* At each output time, we compute and print the max-norm of the */
/* solution (which should decay exponentially in t).  We also print */
/* some relevant statistics -- the current method order and step size, */
/* the number of time steps so far, and the numbers of nonlinear and */
/* linear iterations so far. */

/* If a root was found, we flag this, and return to the DDASKR call. */

/* If DDASKR failed in any way (IDID .lt. 0) we print a message and */
/* stop the integration. */
/* ----------------------------------------------------------------------- */

    nout = 11;
    t = 0.;
    tout = .01;
    i__1 = nout;
    for (iout = 1; iout <= i__1; ++iout) {

L45:
	_daskr_ddaskr_((Unknown_fp)resh_, &neq, &t, u, uprime, &tout, info, &rtol, &atol, &
		idid, rwork, &lrw, iwork, &liw, rpar, ipar, (Unknown_fp)_daskr_djacilu_, (
		Unknown_fp)_daskr_dpsolilu_, (Unknown_fp)rtheat_, &nrt, jroot);

	umax = 0.;
	i__2 = neq;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L50: */
/* Computing MAX */
	    d__2 = umax, d__3 = (d__1 = u[i__ - 1], fabs(d__1));
	    umax = MAX(d__2,d__3);
	}

	hu = rwork[6];
	nqu = iwork[7];
	nst = iwork[10];
	nni = iwork[18];
	nli = iwork[19];
	printf(fmt_60, t, umax, nqu, hu, nst, nni, nli);

	if (idid == 5) {
		printf(fmt_61, jroot[0],jroot[1]);
	    goto L45;
	}

	if (idid < 0) {
		printf(fmt_65, t);
	    goto L80;
	}

	tout *= 2.;
/* L70: */
    }

/* Here we display some final statistics for the problem. */
/* The ratio of NLI to NNI is the average dimension of the Krylov */
/* subspace involved in the Krylov linear iterative method. */
L80:
    nst = iwork[10];
    npe = iwork[12];
    nre = iwork[11] + npe * mband;
    liw = iwork[16];
    lrw = iwork[17];
    nni = iwork[18];
    nli = iwork[19];
    nps = iwork[20];
    if (nni != 0) {
	avdim = (real_number) nli / (real_number) nni;
    }
    ncfn = iwork[14];
    ncfl = iwork[15];
    nrte = iwork[35];

    printf(fmt_90, lrw, liw, nst, nre, ipar[29], nrte, npe, nps, nni, nli, avdim, ncfn, ncfl);
    printf(fmt_100, lwpmin, liwpmin);

/* ------  End of main program for DHEATILU example program -------------- */
    exit(0);
    return 0;
} /* MAIN__ */

/* Subroutine */ int uinit_(real_number *u, real_number *uprime, real_number *
	rpar, integer *ipar)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, m;
    static real_number dx, xj, yk;
    static integer neq, ioff;


/* This routine computes and loads the vector of initial values. */
/* The initial U values are given by the polynomial u = 16x(1-x)y(1-y). */
/* The initial UPRIME values are set to zero.  (DDASKR corrects these */
/* during the first time step.) */


    /* Parameter adjustments */
    --ipar;
    --rpar;
    --uprime;
    --u;

    /* Function Body */
    neq = ipar[33];
    m = ipar[34];
    dx = rpar[3];

    i__1 = m + 1;
    for (k = 0; k <= i__1; ++k) {
	yk = k * dx;
	ioff = (m + 2) * k;
	i__2 = m + 1;
	for (j = 0; j <= i__2; ++j) {
	    xj = j * dx;
	    i__ = ioff + j + 1;
	    u[i__] = xj * 16. * (1. - xj) * yk * (1. - yk);
/* L10: */
	}
/* L20: */
    }
    i__1 = neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L30: */
	uprime[i__] = 0.;
    }
    return 0;
/* ------------  End of Subroutine UINIT  -------------------------------- */
} /* uinit_ */

/* Subroutine */ int resh_(real_number *t, real_number *u, real_number *uprime,
	real_number *cj, real_number *delta, integer *ires, real_number *rpar,
	integer *ipar)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, m, m2, neq, ioff;
    static real_number temx, temy, coeff;


/* This is the user-supplied RES subroutine for this example. */
/* It computes the residuals for the 2-D discretized heat equation, */
/* with zero boundary values. */


/* Set problem constants using IPAR and RPAR. */
    /* Parameter adjustments */
    --ipar;
    --rpar;
    --delta;
    --uprime;
    --u;

    /* Function Body */
    neq = ipar[33];
    m = ipar[34];
    coeff = rpar[4];
    m2 = m + 2;

/* Load U into DELTA, in order to set boundary values. */
    i__1 = neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	delta[i__] = u[i__];
    }

/* Loop over interior points, and load residual values. */
    i__1 = m;
    for (k = 1; k <= i__1; ++k) {
	ioff = m2 * k;
	i__2 = m;
	for (j = 1; j <= i__2; ++j) {
	    i__ = ioff + j + 1;
	    temx = u[i__ - 1] + u[i__ + 1];
	    temy = u[i__ - m2] + u[i__ + m2];
	    delta[i__] = uprime[i__] - (temx + temy - u[i__] * 4.) * coeff;
/* L20: */
	}
/* L30: */
    }

    return 0;
/* ------------  End of Subroutine RESH  --------------------------------- */
} /* resh_ */

/* Subroutine */ int rtheat_(integer *neq, real_number *t, real_number *u,
	real_number *up, integer *nrt, real_number *rval, real_number *rpar,
	integer *ipar)
{
    /* System generated locals */
    integer i__1;
    real_number d__1, d__2;

    /* Local variables */
    static integer i__;
    static real_number umax;


/* This routine finds the max of U, and sets RVAL(1) = max(u) - 0.1, */
/* RVAL(2) = max(u) - 0.01. */


    /* Parameter adjustments */
    --up;
    --u;
    --rval;

    /* Function Body */
    umax = 0.;
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
/* Computing MAX */
	d__1 = umax, d__2 = u[i__];
	umax = MAX(d__1,d__2);
    }
    rval[1] = umax - .1;
    rval[2] = umax - .01;

    return 0;
/* ------------  End of Subroutine RTHEAT  ------------------------------- */
} /* rtheat_ */

