/* drbdpre.f -- translated by f2c (version 20100827).
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../solver/ddaskr_types.h"

/* Common Block Declarations */

struct {
    real_number srur;
    integer mp, mpd, mpsq, meshx, meshy, mxmp;
} drpre1_;

#define drpre1_1 drpre1_

/* Table of constant values */

static integer c__4 = 4;
static integer c__0 = 0;

/* ----------------------------------------------------------------------- */

/*    Preconditioner Tools for Reaction-Transport Problems */
/*    Part I: Block-Diagonal Reaction-Based Factor without Grouping */
/*                        14 September 1995 */

/* The following three subroutines -- DMSET2, DRBDJA, DRBDPS -- */
/* are provided to assist in the generation and solution of */
/* preconditioner matrices for problems arising from reaction-transport */
/* systems, as solved with DASPK.  More specifically, they are intended */
/* as tools for preconditioners that include a contribution from the */
/* reaction terms of the system.  These are intended as auxiliary */
/* routines for the user-supplied routines JAC and PSOL called by */
/* DASPK when the Krylov method is selected. */

/* These routines are intended for a DAE system obtained from a system */
/* of reaction-transport PDEs, in which some of the PDE variables obey */
/* evolution equations, and the rest obey algebraic (time-independent) */
/* equations.  See Ref. 2, Section 4.  It is assumed that the */
/* right-hand sides of all the equations have the form of a sum of a */
/* reaction term R and a transport term S, that the transport term */
/* is discretized by finite differences, and that in the spatial */
/* discretization the PDE variables at each spatial point are kept */
/* together.  Thus the DAE system function, in terms of a dependent */
/* variable vector u, has the form */
/*     G(t,u,u') = I_d u' - R(t,u) - S(t,u) ,  where */
/*     I_d = identity matrix with zeros in the positions corresponding */
/*           to the algebraic components, ones in those for the */
/*           evolution (differential) components. */
/*     R(t,u) = the reaction terms (spatial coupling absent), and */
/*     S(t,u) = the spatial transport terms. */

/* As shown in [2], two possible preconditioners for such a system are: */
/* (a) P_R = c I_d - dR/du, based on the reaction term R alone, and */
/* (b) P_SR = (I - (1/c) dS/du) (c I_d - dR/du), the product of two */
/*     factors (in either order), one being P_R and the other being */
/*     based on the spatial term S alone. */
/* Here c is the scalar CJ that is input to the JAC and PSOL routines */
/* provided by the user (1/c is proportional to the step size H). */

/* The routines given here can be used to provide the reaction-based */
/* factor P_R.  More precisely, they provide an approximation A_R to */
/* P_R.  The matrix P_R is block-diagonal, with each block corresponding */
/* to one spatial point.  In A_R, we compute each block by difference */
/* quotient approximations, by way of calls to a user-supplied routine, */
/* subroutine RBLOCK, that evaluates the reaction terms at a single */
/* spatial point.  A_R has one such block for each spatial point in */
/* the mesh.  (For a more economical approximation, see Part II, */
/* on block-grouping in A_R.) */

/* The routines given here are specialized to the case of a 2-D problem */
/* on a rectangular mesh in the x-y plane.  However, they can be easily */
/* modified for a different problem geometry.  It is also assumed */
/* that the PDE variables are ordered so that the differential */
/* variables appear first, followed by the algebraic variables. */

/* To make use of these routines in a DASPK solution, the user must */
/* provide: */
/* (a) a calling program that sets the DASPK input parameters, and calls */
/*     DMSET2 to set mesh data and mesh-related DASPK inputs; */
/* (b) a JAC routine, as prescribed by the DASPK instructions, which */
/*     calls DRBDJA, and does any other Jacobian-related preprocessing */
/*     needed for preconditioning; and */
/* (c) a PSOL routine, as prescribed by the DASPK instructions, which */
/*     calls DRBDPS for the solution of systems A_R x = b, and does */
/*     any other linear system solving required by the preconditioner. */
/* Detailed descriptions and instructions are given below. */

/* In addition, the use of these routines requires: */
/*  * the LINPACK routines DGEFA and DGESL for dense linear sytems, and */
/*  * the machine constant routine D1MACH for the machine unit roundoff. */

/* (a) The calling program. */
/* The calling program sets the DASPK inputs and makes calls to DASPK. */
/* Here the DASPK inputs include */
/*   INFO(12) = 1 [to signal the Krylov method] */
/*   INFO(15) = 1 [to signal the presence of a JAC routine] */

/* Also, the use of the DRBDJA/DRBDPS routines in conjunction with */
/* DASPK requires that certain mesh-related data be set.  This can be */
/* done with the call */
/*     CALL DMSET2 (MX, MY, NS, NSD, LID, IWORK) */
/* The input arguments to DMSET2 are: */
/*   MX and MY = the mesh dimensions. */
/*   NS  = number of PDE variables. */
/*   NSD = number of differential PDE variables. */
/*   LID = offset in IWORK for array showing the differential and */
/*         algebraic components on input to DASPK, required if either */
/*         INFO(11) = 1 or INFO(16) = 1.  Set LID = 0 otherwise. */
/*         If this array is required, set LID = 40 or 40 + NEQ, */
/*         depending on the value of the constraint option INFO(10). */
/* DMSET2 loads mesh data in a COMMON block /DRPRE1/ used by the */
/* DRBDJA/DRBDPS routines. */

/* DMSET2 also loads the preconditioner work lengths into */
/* IWORK(27) and IWORK(28), and if LID > 0 it sets the ID array */
/* in IWORK showing the differential and algebraic components. */

/* (b) The JAC routine. */
/* The user-supplied JAC routine called by DASPK with the Krylov */
/* method specified, is to generate and preprocess Jacobian-related */
/* data as needed for later solution of the preconditioner system */
/* P x = b.  Assuming that P is to be an approximation of either P_R */
/* or P_SR, the JAC routine should call DRBDJA for the approximation */
/* A_R to P_R.  Subroutine DRBDJA generates A_R using difference */
/* quotients.  It then performs an LU decomposition of each block, */
/* using the LINPACK routine DGEFA. */

/* In terms of the arguments passed to JAC by DASPK, the call to */
/* DRBDJA should have the form */
/*     CALL DRBDJA (T, U, R0, RBLOCK, WK, REWT, CJ, WP, IWP, IER) */
/* where we use U instead of Y for the dependent variable array. */
/* The argument R0 is an array assumed to contain the current value */
/* of the R vector, at the current values (T,U).  This can be done, for */
/* example, by taking R0 to be RPAR, and loading RPAR with the */
/* vector R in the last call to the RES routine; in that case, the */
/* calling program must declare RPAR to have length at least NEQ. */
/* Alternatively, insert a call to RBLOCK (see below) within the */
/* loop over mesh points in DRBDJA. */

/* To use DRBDJA, the user must provide the following subroutine, */
/* which DRBDJA calls to obtain individual blocks of R: */
/*      SUBROUTINE RBLOCK (T, JX, JY, UXY, RXY) */
/* The input arguments to RBLOCK are: */
/*   T     = current time. */
/*   JX,JY = spatial indices in x- and y-directions. */
/*   UXY   = block of NS dependent variables at spatial point (JX,JY). */
/* RBLOCK is to load block (JX,JY) of R(t,u) into the array RXY. */

/* (c) The PSOL routine. */
/* The user-supplied PSOL routine must solve the linear system P x = b, */
/* where P is the preconditioner matrix.  For this, the PSOL routine */
/* should call DRBDPS for the solution of A_R.  Subroutine DRBDPS */
/* solves a linear system A_R x = b, using the LINPACK backsolve */
/* routine DGESL.  In terms of the arguments passed to PSOL by DASPK, */
/* the call to DRBDPS should have the form */
/*     CALL DRBDPS (B, WP, IWP) */
/* DRBDPS overwrites the B array (containing b) with the solution x. */

/* ----------------------------------------------------------------------- */

/* References */
/* [1] Peter N. Brown and Alan C. Hindmarsh, */
/*     Reduced Storage Matrix Methods in Stiff ODE Systems, */
/*     J. Appl. Math. & Comp., 31 (1989), pp. 40-91. */
/* [2] Peter N. Brown, Alan C. Hindmarsh, and Linda R. Petzold, */
/*     Using Krylov Methods in the Solution of Large-Scale Differential- */
/*     Algebraic Systems, SIAM J. Sci. Comput., 15 (1994), pp. 1467-1488. */
/* [3] Peter N. Brown, Alan C. Hindmarsh, and Linda R. Petzold, */
/*     Consistent Initial Condition Calculation for Differential- */
/*     Algebraic Systems, LLNL Report UCRL-JC-122175, August 1995; */
/*     submitted to SIAM J. Sci. Comp. */
/* ----------------------------------------------------------------------- */
/* Subroutine */ int _daskr_dmset2_(integer *mx, integer *my, integer *ns, integer *
	nsd, integer *lid, integer *iwork)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(real_number);

    /* Local variables */
    static integer i__, i0, jx, jy;
    extern real_number _daskr_d1mach_(integer *);
    static real_number uround;

/* ***BEGIN PROLOGUE  DMSET2 */
/* ***DATE WRITTEN   950830   (YYMMDD) */

/* ***AUTHORS  A. C. Hindmarsh */
/*            Lawrence Livermore National Laboratory */
/*            L-316, P.O. Box 808 */
/*            Livermore, CA 94551 */

/* ***DESCRIPTION */

/* ----------------------------------------------------------------------- */
/* This routine sets mesh parameters needed to use the routines */
/* DRBDJA and DRBDPS, assuming a 2-D rectangular problem. */
/* Given the mesh parameters, it loads the COMMON block /DRPRE1/, */
/* and the lengths LENWP and LENIWP in IWORK. */
/* Then if LID > 0, it also sets the ID array in IWORK, indicating */
/* which components are differential and which are algebraic. */

/* The variables in the COMMON block are defined as follows: */
/*   SRUR   = SQRT(unit roundoff), used in difference quotients. */
/*            UROUND = D1MACH(4) generates the unit roundoff. */
/*   MP     = NS = number of PDE variables, the size of each block in */
/*            the block-diagonal preconditioner matrix P_R. */
/*   MPD    = NSD = number of differential PDE variables.  In the DAE */
/*            system, the first MPD variables at each spatial point have */
/*            time  derivatives, and the remaining (MP - MPD) do not. */
/*   MPSQ   = MP*MP. */
/*   MESHX  = MX = x mesh size. */
/*   MESHY  = MY = y mesh size (the mesh is MESHX by MESHY). */
/* ----------------------------------------------------------------------- */
/* ***ROUTINES CALLED */
/*   D1MACH */

/* ***END PROLOGUE  DMSET2 */


/* Load the COMMON block. */
    /* Parameter adjustments */
    --iwork;

    /* Function Body */
    uround = _daskr_d1mach_(&c__4);
    drpre1_1.srur = sqrt(uround);
    drpre1_1.mp = *ns;
    drpre1_1.mpd = *nsd;
    drpre1_1.mpsq = *ns * *ns;
    drpre1_1.meshx = *mx;
    drpre1_1.meshy = *my;
    drpre1_1.mxmp = drpre1_1.meshx * drpre1_1.mp;

/* Here set the sizes of the preconditioning storage space segments */
/* in RWORK and IWORK. */
    iwork[27] = drpre1_1.mpsq * drpre1_1.meshx * drpre1_1.meshy;
    iwork[28] = drpre1_1.mp * drpre1_1.meshx * drpre1_1.meshy;

/* If LID .GT. 0, set the ID array in IWORK. */
    if (*lid == 0) {
	return 0;
    }
    i0 = *lid;
    i__1 = *my;
    for (jy = 1; jy <= i__1; ++jy) {
	i__2 = *mx;
	for (jx = 1; jx <= i__2; ++jx) {
	    i__3 = drpre1_1.mpd;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/* L10: */
		iwork[i0 + i__] = 1;
	    }
	    i__3 = drpre1_1.mp;
	    for (i__ = drpre1_1.mpd + 1; i__ <= i__3; ++i__) {
/* L20: */
		iwork[i0 + i__] = -1;
	    }
	    i0 += drpre1_1.mp;
/* L30: */
	}
/* L40: */
    }

    return 0;
/* ------------  End of Subroutine DMSET2  ------------------------------- */
} /* dmset2_ */

/* Subroutine */ int _daskr_drbdja_(real_number *t, real_number *u, real_number *r0,
	Unknown_fp rblock, real_number *r1, real_number *rewt, real_number *cj,
	real_number *bd, integer *ipbd, integer *ier)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    real_number d__1, d__2;

    /* Local variables */
    static integer i__, j, j0, js;
    static real_number uj;
    static integer jx, jy;
    static real_number fac;
    static integer ibd;
    static real_number del;
    static integer iip;
    static real_number dfac;
    extern /* Subroutine */ int _daskr_dgefa_(real_number *, integer *, integer *,
	    integer *, integer *);
    static integer idiag;

/* ***BEGIN PROLOGUE  DRBDJA */
/* ***DATE WRITTEN   950914   (YYMMDD) */

/* ***AUTHORS  A. C. Hindmarsh */
/*            Lawrence Livermore National Laboratory */
/*            L-316, P.O. Box 808 */
/*            Livermore, CA 94551 */

/* ***DESCRIPTION */

/* ----------------------------------------------------------------------- */
/* This routine generates and preprocesses a block-diagonal */
/* preconditioner matrix, based on the part of the Jacobian corresponding */
/* to the reaction terms R of the problem. */
/* It generates a matrix of the form CJ * I_d - dR/du. */
/* It calls DGEFA to do LU decomposition of each diagonal block. */
/* The computation of the diagonal blocks uses the mesh information in */
/* the COMMON block /DRPRE1/.  One block per spatial point is computed. */
/* The Jacobian elements are generated by difference quotients. */
/* This routine calls a user routine of the form */
/*      SUBROUTINE RBLOCK (T, JX, JY, UXY, RXY) */
/* which is to set RXY to block (JX,JY) of R, as a function of the */
/* current time T and block UXY of current dependent variable vector U. */
/* The array R0 is assumed to contain the current value of R at (T,U). */
/* ----------------------------------------------------------------------- */
/* On input: */
/*   T      = current value of independent variable. */
/*   U      = current dependent variable array. */
/*   R0  = array of current values of the vector R at (T,U) */
/*   RBLOCK = name of external routine that computes a single block of R. */
/*   R1     = array of length NEQ for work space. */
/*   REWT   = reciprocal error weights. */
/*   CJ     = scalar used in forming the system Jacobian. */

/* On output: */
/*   BD     = array containing the LU factors of the diagonal blocks. */
/*   IPBD   = integer array of pivots for the LU factorizations. */
/*   IER    = integer error flag.  If no error occurred, IER = 0. */
/*            If a zero pivot was found at stage k in one of the LU */
/*            factorizations, this routine returns IER = k > 0. */
/* Here BD is the RWORK segment WP, and IPBD is the IWORK segment IWP. */
/* ----------------------------------------------------------------------- */
/* ***ROUTINES CALLED */
/*   RBLOCK, DGEFA */

/* ***END PROLOGUE  DRBDJA */



/* Make MP calls to RBLOCK to approximate each diagonal block of dR/du. */
    /* Parameter adjustments */
    --ipbd;
    --bd;
    --rewt;
    --r1;
    --r0;
    --u;

    /* Function Body */
    dfac = .01;
    ibd = 0;
    j0 = 0;
    i__1 = drpre1_1.meshy;
    for (jy = 1; jy <= i__1; ++jy) {
	i__2 = drpre1_1.meshx;
	for (jx = 1; jx <= i__2; ++jx) {
/* If R0 has not been set previously as an array of length NEQ, it can */
/* be set here, as an array of length MP, with the call */
/*         CALL RBLOCK (T, JX, JY, U(J0+1), R0) */
/* In this case, change R0(J0+I) below to R0(I). */
	    i__3 = drpre1_1.mp;
	    for (js = 1; js <= i__3; ++js) {
		j = j0 + js;
		uj = u[j];
/* Computing MAX */
		d__1 = drpre1_1.srur * fabs(uj), d__2 = dfac / rewt[j];
		del = MAX(d__1,d__2);
		u[j] += del;
		fac = -1. / del;
		(*rblock)(t, &jx, &jy, &u[j0 + 1], &r1[1]);
		i__4 = drpre1_1.mp;
		for (i__ = 1; i__ <= i__4; ++i__) {
/* L10: */
		    bd[ibd + i__] = (r1[i__] - r0[j0 + i__]) * fac;
		}
		u[j] = uj;
		ibd += drpre1_1.mp;
/* L20: */
	    }
	    j0 += drpre1_1.mp;
/* L30: */
	}
/* L40: */
    }

/* Add matrix CJ * I_d, and do LU decomposition on blocks. -------------- */
    ibd = 1;
    iip = 1;
    i__1 = drpre1_1.meshx * drpre1_1.meshy;
    for (j = 1; j <= i__1; ++j) {
	idiag = ibd;
	i__2 = drpre1_1.mp;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (i__ <= drpre1_1.mpd) {
		bd[idiag] += *cj;
	    }
/* L70: */
	    idiag += drpre1_1.mp + 1;
	}
	_daskr_dgefa_(&bd[ibd], &drpre1_1.mp, &drpre1_1.mp, &ipbd[iip], ier);
	if (*ier != 0) {
	    goto L90;
	}
	ibd += drpre1_1.mpsq;
	iip += drpre1_1.mp;
/* L80: */
    }
L90:
    return 0;
/* ------------  End of Subroutine DRBDJA  ------------------------------- */
} /* drbdja_ */

/* Subroutine */ int _daskr_drbdps_(real_number *b, real_number *bd, integer *ipbd)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer ib, jx, jy, ibd, ier;
    extern /* Subroutine */ int _daskr_dgesl_(real_number *, integer *, integer *,
	    integer *, real_number *, integer *);

/* ***BEGIN PROLOGUE  DRBDPS */
/* ***DATE WRITTEN   950914   (YYMMDD) */

/* ***AUTHORS  A. C. Hindmarsh */
/*            Lawrence Livermore National Laboratory */
/*            L-316, P.O. Box 808 */
/*            Livermore, CA 94551 */

/* ***DESCRIPTION */

/* ----------------------------------------------------------------------- */
/* This routine solves a linear system A_R x = b, using the LU factors */
/* of the diagonal blocks computed in DRBDJA, and mesh parameters */
/* in the COMMON block /DRPRE1/. */
/* Here BD is the RWORK segment WP, and IPBD is the IWORK segment IWP. */
/* The right-hand side vector b, contained in B on entry, is overwritten */
/* with the solution vector x on return. */
/* ----------------------------------------------------------------------- */
/* ***ROUTINES CALLED */
/*   DGESL */

/* ***END PROLOGUE  DRBDPS */



    /* Parameter adjustments */
    --ipbd;
    --bd;
    --b;

    /* Function Body */
    ier = 0;
    ib = 1;
    ibd = 1;
    i__1 = drpre1_1.meshy;
    for (jy = 1; jy <= i__1; ++jy) {
	i__2 = drpre1_1.meshx;
	for (jx = 1; jx <= i__2; ++jx) {
	    _daskr_dgesl_(&bd[ibd], &drpre1_1.mp, &drpre1_1.mp, &ipbd[ib], &b[ib], &
		    c__0);
	    ib += drpre1_1.mp;
	    ibd += drpre1_1.mpsq;
/* L10: */
	}
/* L20: */
    }

    return 0;
/* ------------  End of Subroutine DRBDPS  ------------------------------- */
} /* drbdps_ */

