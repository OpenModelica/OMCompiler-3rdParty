/* dbanpre.f -- translated by f2c (version 20100827).
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../solver/ddaskr_types.h"

/* ----------------------------------------------------------------------- */

/*             Preconditioner Routines for Banded Problems */
/*                          14 September 1995 */

/* The following pair of subroutines -- DBANJA and DBANPS -- provides a */
/* general-purpose banded preconditioner matrix for use with the DDASPK */
/* solver, with the Krylov linear system method.  When using DDASPK to */
/* solve a problem G(t,y,y') = 0, whose iteration matrix (Jacobian) */
/*    J = dG/dy + c * dG/dy'  (c = scalar) */
/* is either banded or approximately equal to a banded matrix, these */
/* routines can be used to generate a banded approximation to J as the */
/* preconditioner and to solve the resulting banded linear system, in */
/* conjunction with the Krylov method option (INFO(12) = 1) in DDASPK. */

/* Other than the user-supplied residual routine RES defining G(t,y,y'), */
/* the only other inputs required by these routines are the */
/* half-bandwidth parameters ML and MU of the approximate banded */
/* Jacobian.  If the system size is NEQ, the half-bandwidths are */
/* defined as integers between 0 and NEQ - 1 such that only elements */
/* with indices (i,j) satisfying */
/*    -ML .le. j - i .le. MU */
/* are to be retained in the preconditioner.  E.g., if ML = MU = 0, a */
/* diagonal matrix will be generated as the preconditioner.  The banded */
/* preconditioner is obtained by difference quotient approximations.  If */
/* the true problem Jacobian is not banded but is approximately equal to */
/* a matrix that is banded, the procedure used here will have the effect */
/* of lumping the elements outside of the band onto the elements within */
/* the band. */

/* To use these routines in conjunction with DDASPK, the user's calling */
/* program should include the following, in addition to setting the other */
/* DDASPK input parameters. */

/* (a) Dimension the array IPAR to have length at least 2, and load the */
/*     half-bandwidths into IPAR as */
/*       IPAR(1) = ML   and   IPAR(2) = MU */
/*     IPAR is used to communicate these parameters to DBANJA and DBANPS. */
/*     If the user program also uses IPAR for communication with RES, */
/*     that data should be located beyond the first 2 words of IPAR. */

/* (b) Include the names DBANJA and DBANPS in an EXTERNAL statement. */
/*     Set INFO(15) = 1 to indicate that a JAC routine exists. */
/*     Then in the call to DDASPK, pass the names DBANJA and DBANPS as */
/*     the arguments JAC and PSOL, respectively. */

/* (c) The DDASPK work arrays RWORK and IWORK must include segments WP */
/*     and IWP for use by DBANJA/DBANPS.  The lengths of these depend on */
/*     the problem size and half-bandwidths, as follows: */
/*       LWP =  length of RWORK segment WP = */
/*                     (2*ML + MU + 1)*NEQ + 2*( (NEQ/(ML+MU+1)) + 1) */
/*       LIWP = length of IWORK segment IWP = NEQ */
/*     (Note the integer divide in LWP.)  Load these lengths in IWORK as */
/*       IWORK(27) = LWP */
/*       IWORK(28) = LIWP */
/*     and include these values in the declared size of RWORK and IWORK. */


/* The DBANJA and DBANPS routines generate and solve the banded */
/* preconditioner matrix P within the preconditioned Krylov algorithm */
/* used by DDASPK when INFO(12) = 1.  P is generated and LU-factored */
/* periodically during the integration, and the factors are used to */
/* solve systems Px = b as needed. */
/* ----------------------------------------------------------------------- */
/* Subroutine */ int _daskr_dbanja_(Unknown_fp res, integer *ires, integer *neq,
	real_number *t, real_number *y, real_number *yprime, real_number *rewt,
	real_number *savr, real_number *wk, real_number *h__, real_number *cj,
	real_number *wp, integer *iwp, integer *ier, real_number *rpar, integer
	*ipar)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    real_number d__1, d__2, d__3, d__4, d__5;

    /* Builtin functions */
    real_number _daskr_real_sign(real_number *, real_number *);

    /* Local variables */
    static integer i__, j, k, n, i1, i2, ii, ml, mu, mba;
    static real_number del;
    static integer meb1, lenp;
    static real_number squr;
    extern /* Subroutine */ int _daskr_dgbfa_(real_number *, integer *, integer *,
	    integer *, integer *, integer *, integer *);
    static integer mband, isave, msave;
    extern real_number _daskr_d1mach_(integer *);
    static integer meband;
    static real_number delinv;
    static integer ipsave;
    static real_number uround;

    /* Table of constant values */
    static integer c__4 = 4;


/* ***BEGIN PROLOGUE  DBANJA */
/* ***DATE WRITTEN   891204   (YYMMDD) */
/* ***REVISION DATE  900122 */
/* ***REVISION DATE  920929   CJ in RES call sequence */
/* ***REVISION DATE  950914   Name change, minor revisions throughout */
/* ***AUTHORS  L. R. Petzold, P. N. Brown, A. C. Hindmarsh, C. W. Ulrich */
/*            Numerical Mathematics Group */
/*            Lawrence Livermore National Laboratory */
/*            Livermore, CA 94551 */

/* ***DESCRIPTION */

/* Subroutine DBANJA generates a banded preconditioner matrix P that */
/* approximates the DDASPK iteration matrix  J = dG/dy + CJ*dG/dy', */
/* where the DAE system is  G(t,y,y') = 0.  The band matrix P has */
/* half-bandwidths ML and MU.  It is computed by making (ML + MU + 1) */
/* calls to the user's RES routine and forming difference quotients, */
/* exactly as in the banded direct method option of DDASPK. */
/* DBANJA calls the LINPACK routine DGBFA to do an LU factorization of */
/* this matrix. */

/* The call sequence parameters have the following meanings. */

/*     RES      = External user-supplied subroutine to evaluate the */
/*                residuals.  See RES description in DDASPK prologue. */
/*     IRES     = Output flag set by RES.  See RES description in DDASPK. */
/*     NEQ      = Problem size. */
/*     T        = Independent variable t. */
/*     Y        = Array containing current dependent variables y. */
/*     YPRIME   = Array containing current derivative y'. */
/*     REWT     = Vector of reciprocal error weights, used here for */
/*                computing increments. */
/*     SAVR     = Current residual evaluated at (T,Y,YPRIME). */
/*     WK       = Real work space of length NEQ. */
/*     H        = Current step size. */
/*     CJ       = Scalar proportional to 1/H. */
/*     WP       = Real work array for P etc.  On output, it contains */
/*                the LU decomposition of the banded approximation P. */
/*     IWP      = Integer work space for matrix pivot information. */
/*     IER      = Output flag, > 0 if P is singular, and 0 otherwise. */
/*     RPAR,IPAR= Real and integer arrays used for communication between */
/*                the calling program and external user routines. */
/*                IPAR(1) and IPAR(2) must contain ML and MU, resp. */
/*                RPAR is not used here. */

/* ***ROUTINES CALLED */
/*   D1MACH, DGBFA, RES */

/* ***END PROLOGUE  DBANJA */


/* Set band parameters. */
    /* Parameter adjustments */
    --ipar;
    --rpar;
    --iwp;
    --wp;
    --wk;
    --savr;
    --rewt;
    --yprime;
    --y;

    /* Function Body */
    ml = ipar[1];
    mu = ipar[2];
    mband = ml + mu + 1;
    mba = MIN(mband,*neq);
    meband = mband + ml;
    meb1 = meband - 1;

/* Set the machine unit roundoff UROUND and SQRT(UROUND), used to */
/* set increments in the difference quotient procedure. */
    uround = _daskr_d1mach_(&c__4);
    squr = sqrt(uround);

/* Set pointers into WP.  LENP is the length of the segment for P. */
/* Following that are two segments of size (NEQ/MBAND), with offsets */
/* ISAVE and IPSAVE, for temporary storage of Y and YPRIME elements. */
    lenp = ((ml << 1) + mu + 1) * *neq;
    msave = *neq / mband + 1;
    isave = lenp;
    ipsave = isave + msave;

/* Initialize error flags. */
    *ier = 0;
    *ires = 0;

/* Generate the banded approximate iteration matrix P using */
/* difference quotients on the results of calls to RES. */

    i__1 = mba;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *neq;
	i__3 = mband;
	for (n = j; i__3 < 0 ? n >= i__2 : n <= i__2; n += i__3) {
	    k = (n - j) / mband + 1;
	    wp[isave + k] = y[n];
	    wp[ipsave + k] = yprime[n];
/* Computing MAX */
	    d__4 = (d__1 = y[n], fabs(d__1)), d__5 = (d__2 = *h__ * yprime[n],
		    fabs(d__2)), d__4 = MAX(d__4,d__5), d__5 = (d__3 = 1. /
		    rewt[n], fabs(d__3));
	    del = squr * MAX(d__4,d__5);
	    d__1 = *h__ * yprime[n];
	    del = _daskr_real_sign(&del, &d__1);
	    del = y[n] + del - y[n];
	    y[n] += del;
	    yprime[n] += *cj * del;
/* L10: */
	}
	(*res)(t, &y[1], &yprime[1], cj, &wk[1], ires, &rpar[1], &ipar[1]);
	if (*ires < 0) {
	    return 0;
	}
	i__3 = *neq;
	i__2 = mband;
	for (n = j; i__2 < 0 ? n >= i__3 : n <= i__3; n += i__2) {
	    k = (n - j) / mband + 1;
	    y[n] = wp[isave + k];
	    yprime[n] = wp[ipsave + k];
/* Computing MAX */
	    d__4 = (d__1 = y[n], fabs(d__1)), d__5 = (d__2 = *h__ * yprime[n],
		    fabs(d__2)), d__4 = MAX(d__4,d__5), d__5 = (d__3 = 1. /
		    rewt[n], fabs(d__3));
	    del = squr * MAX(d__4,d__5);
	    d__1 = *h__ * yprime[n];
	    del = _daskr_real_sign(&del, &d__1);
	    del = y[n] + del - y[n];
	    delinv = 1. / del;
/* Computing MAX */
	    i__4 = 1, i__5 = n - mu;
	    i1 = MAX(i__4,i__5);
/* Computing MIN */
	    i__4 = *neq, i__5 = n + ml;
	    i2 = MIN(i__4,i__5);
	    ii = n * meb1 - ml;
	    i__4 = i2;
	    for (i__ = i1; i__ <= i__4; ++i__) {
/* L20: */
		wp[ii + i__] = (wk[i__] - savr[i__]) * delinv;
	    }
/* L30: */
	}
/* L40: */
    }

/* Do LU decomposition of the band matrix P. */

    _daskr_dgbfa_(&wp[1], &meband, neq, &ml, &mu, &iwp[1], ier);
    return 0;

/* ------------  End of Subroutine DBANJA  ------------------------------- */
} /* dbanja_ */

/* Subroutine */ int _daskr_dbanps_(integer *neq, real_number *t, real_number *y,
	real_number *yprime, real_number *savr, real_number *wk, real_number *cj,
	real_number *wght, real_number *wp, integer *iwp, real_number *b,
	real_number *eplin, integer *ier, real_number *rpar, integer *ipar)
{
    static integer ml, mu;
    extern /* Subroutine */ int _daskr_dgbsl_(real_number *, integer *, integer *,
	    integer *, integer *, integer *, real_number *, integer *);
    static integer meband;

    /* Table of constant values */
    static integer c__0 = 0;


/* ***BEGIN PROLOGUE  DBANPS */
/* ***DATE WRITTEN   891204   (YYMMDD) */
/* ***REVISION DATE  900110   (YYMMDD) */
/* ***REVISION DATE  950914   Name change, minor revisions throughout */
/* ***AUTHORS  L. R. Petzold, P. N. Brown, A. C. Hindmarsh, C. W. Ulrich */
/*            Numerical Mathematics Group */
/*            Lawrence Livermore National Laboratory */
/*            Livermore, CA 94551 */

/* ***DESCRIPTION */

/* Subroutine DBANPS uses the factors produced by DBANJA to solve linear */
/* systems P x = b for the banded preconditioner P and a given vector b. */
/* It calls the LINPACK routine SGBSL for this. */

/* The call sequence parameters have the following meanings. */

/*     NEQ      = Problem size. */
/*     T        = Independent variable t (not used). */
/*     Y        = Array containing current dependent vars. (not used). */
/*     YPRIME   = Array containing current derivative (not used). */
/*     SAVR     = Current residual evaluated at (T,Y,YPRIME) (not used). */
/*     WK       = Real work space of length NEQ (not used). */
/*     CJ       = Scalar proportional to 1/H (H = step size) (not used). */
/*     WGHT     = Vector of error weights for computing norms (not used). */
/*     WP       = Real work array containing the LU decomposition of P. */
/*     IWP      = Integer array containing matrix pivot information. */
/*     B        = Right-hand side vector on input; solution on output. */
/*     EPLIN    = tolerance on linear system solver (not used). */
/*     IER      = Output error flag (not used; assumed 0 on input). */
/*     RPAR,IPAR= Real and integer arrays used for communication between */
/*                the calling program and external user routines. */
/*                IPAR(1) and IPAR(2) must contain ML and MU, resp. */
/*                RPAR is not used here. */

/* ***ROUTINES CALLED */
/*   DGBSL */

/* ***END PROLOGUE  DBANPS */


    /* Parameter adjustments */
    --ipar;
    --rpar;
    --b;
    --iwp;
    --wp;

    /* Function Body */
    ml = ipar[1];
    mu = ipar[2];
    meband = (ml << 1) + mu + 1;
    _daskr_dgbsl_(&wp[1], &meband, neq, &ml, &mu, &iwp[1], &b[1], &c__0);
    return 0;

/* ------------  End of Subroutine DBANPS  ------------------------------- */
} /* dbanps_ */

