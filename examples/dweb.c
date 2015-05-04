/* dweb.f -- translated by f2c (version 20090411).
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../solver/ddaskr_types.h"

/* Common Block Declarations */

struct {
    real_number aa, ee, gg, bb, dprey, dpred;
} ppar1_;

#define ppar1_1 ppar1_

struct {
    integer np, ns;
    real_number ax, ay, acoef[4]	/* was [2][2] */, bcoef[2], dx, dy, alph,
	    beta, fpi, diff[2], cox[2], coy[2];
    integer mx, my, mxns;
} ppar2_;

#define ppar2_1 ppar2_

/* Table of constant values */

static integer c__1 = 1;
static integer c__40 = 40;
static integer c_b58 = 104883;
static integer c__1640 = 1640;

/* ***BEGIN PROLOGUE  DWEB */
/* ***REFER TO  DDASKR */
/* ***DATE WRITTEN   020813   (YYMMDD) */
/* ***REVISION DATE  021217   Added JROOT output value. */

/* ***AUTHORS  A. C. Hindmarsh, P. N. Brown */
/*            Lawrence Livermore National Laboratory */
/*            Livermore, CA 94551, USA */

/* ***DESCRIPTION */

/* ----------------------------------------------------------------------- */
/* Example program for DDASKR. */
/* DAE system derived from ns-species interaction PDE in 2 dimensions. */

/* This is the double precision version. */
/* ----------------------------------------------------------------------- */

/* This program solves a DAE system that arises from a system */
/* of partial differential equations.  The PDE system is a food web */
/* population model, with predator-prey interaction and diffusion on */
/* the unit square in two dimensions.  The dependent variable vector is */

/*         1   2        ns */
/*   c = (c , c , ..., c  ) */

/* and the PDEs are as follows.. */

/*     i               i      i */
/*   dc /dt  =  d(i)*(c    + c   )  +  R (x,y,c)  (i=1,...,ns/2) */
/*                     xx     yy        i */

/*                     i      i */
/*   0       =  d(i)*(c    + c   )  +  R (x,y,c)  (i=(ns/2)+1,...,ns) */
/*                     xx     yy        i */

/* where */
/*                  i          ns         j */
/*   R (x,y,c)  =  c *(b(i) + sum a(i,j)*c ) */
/*    i                       j=1 */

/* The number of species is ns = 2*np, with the first np being prey and */
/* the last np being predators.  The coefficients a(i,j), b(i), d(i) are */

/*   a(i,i) = -a  (all i) */
/*   a(i,j) = -g  (i .le. np, j .gt. np) */
/*   a(i,j) =  e  (i .gt. np, j .le. np) */
/*   b(i) =  b*(1 + alpha*x*y + beta*sin(4*pi*x)*sin(4*pi*y))  (i .le. np) */
/*   b(i) = -b*(1 + alpha*x*y + beta*sin(4*pi*x)*sin(4*pi*y))  (i .gt. np) */
/*   d(i) = dprey  (i .le. np) */
/*   d(i) = dpred  (i .gt. np) */

/* The various scalar parameters are set in subroutine SETPAR. */

/* The boundary conditions are.. normal derivative = 0. */
/* A polynomial in x and y is used to set the initial conditions. */

/* The PDEs are discretized by central differencing on a MX by MY mesh. */

/* The root function is R(T,Y,Y') = average(c1) - 20. */

/* The DAE system is solved by DDASKR with three different method options: */
/* (1) direct band method for the linear systems (internal Jacobian), */
/* (2) preconditioned Krylov method for the linear systems, without */
/*     block-grouping in the reaction-based factor, and */
/* (3) preconditioned Krylov method for the linear systems, with */
/*     block-grouping in the reaction-based factor. */

/* In the Krylov cases, the preconditioner is the product of two factors: */
/* (a) The spatial factor uses a fixed number of Gauss-Seidel iterations */
/* based on the diffusion terms only. */
/* (b) The reaction-based factor is a block-diagonal matrix based on */
/* the partial derivatives of the interaction terms R only. */
/* With block-grouping, only a subset of the ns by ns blocks are computed. */
/* An integer flag, JPRE, is set in the main program to specify whether */
/* the preconditioner is to use only one of the two factors or both, */
/* and in which order. */

/* The reaction-based preconditioner factor is set up and solved in */
/* seven subroutines -- */
/*   DMSET2, DRBDJA, DRBDPS  in the case of no block-grouping, and */
/*   DGSET2, GSET1, DRBGJA, DRBGPS  in the case of block-grouping. */
/* These routines are provided separately for general use on problems */
/* arising from reaction-transport systems. */

/* Two output files are written.. one with the problem description and */
/* performance statistics on unit LOUT = 9, and one with solution */
/* profiles at selected output times on unit LCOUT = 10. */
/* The solution file is written only in the case of the direct method. */
/* ----------------------------------------------------------------------- */
/* Note.. in addition to the main program and subroutines given below, */
/* this program requires the BLAS routine DAXPY. */
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
/* ***ROUTINES CALLED */
/*   SETPAR, DGSET2, CINIT, DDASKR, OUTWEB */

/* ***END PROLOGUE  DWEB */

/* Main program */ int main(void)
{
    /* Initialized data */

    static integer lcout = 10;

    /* Format strings */
    static char fmt_20[] = " DWEB: Example program for DDASKR package\n\n"
    	" Food web problem with NS species, NS =%4d\n"
    	" Predator-prey interaction and diffusion on a 2-D square\n";
    static char fmt_25[] =
    	" Matrix parameters..  a =%e12.4E   e =%e12.4E   g =%e12.4E\n"
    	" 					   b parameter =%12.4E\n"
	    " Diffusion coefficients.. dprey =%e12.4E   dpred =%e12.4E\n"
    	" Rate parameters alpha =%e12.4E  and beta =%e12.4E\n";
    static char fmt_30[] = " Mesh dimensions (MX,MY) =%4d %4d     "
    		"Total system size is NEQ =%7d\n";
    static char fmt_35[] = " Root function is R(Y) = average(c1) - 20\n";
    static char fmt_70[] =
    	" Tolerance parameters.. RTOL =%10.2E   ATOL =%10.2E\n"
    	" Internal I.C. calculation flag INFO(11) =%2d  (0 = off, 1 = on)\n"
    	" Predator I.C. guess =%10.2E\n"
    	" Alternate error test flag INFO(16) =%2d  (0 = off, 1 = on)\n";
    static char fmt_80[] =
    	"\n--------------------------------------------------------------------------------\n"
    	" Linear solver method flag INFO(12) =%2d   (0 = direct, 1 = Krylov)\n";
    static char fmt_90[] = " Difference-quotient banded Jacobian, half-bandwidths =%4d\n";
    static char fmt_100[] = " Preconditioner flag is JPRE =%3d\n"
    	"  (1 = reaction factor A_R, 2 = spatial factor A_S, 3 ="
	    " A_S*A_R, 4 = A_R*A_S )\n";
    static char fmt_110[] = " No block-grouping in reaction factor\n";
    static char fmt_120[] = " Block-grouping in reaction factor\n"
	    " Number of groups =%5d   (NGX by NGY, NGX =%3d,  NGY =%3d)\n";
    static char fmt_140[] = "\n   t             Ave.c1   NSTEP   NRE   NNI   NLI "
    		"  NPE   NQ     H          AVLIN\n";
    static char fmt_160[] = "%13.5E %10.5f %5d %6d %5d %5d %5d %4d %11.2E %9.4f\n";
    static char fmt_165[] = "\t\t*****   Root found, JROOT =%3d\n";
    static char fmt_170[] = "\n\n Final time reached =%12.4E\n\n";
    static char fmt_220[] = "\n\n Final statistics for this run..\n"
	    "\tRWORK size =%8d   IWORK size =%6d\n"
    	"\tNumber of time steps              =%5d\n"
    	"\tNumber of residual evaluations    =%5d\n"
    	"\tNumber of root fn. evaluations    =%5d\n"
    	"\tNumber of Jac. or prec. evals.    =%5d\n"
    	"\tNumber of preconditioner solves   =%5d\n"
    	"\tNumber of nonlinear iterations    =%5d\n"
    	"\tNumber of linear iterations       =%5d\n"
    	"\tAverage Krylov subspace dimension =%8.4F\n"
    	"\t%3d nonlinear conv. failures, %5d linear conv. failures\n";


    /* System generated locals */
    integer i__1;
    real_number d__1;

    /* Local variables */
    static integer i__;
    static real_number t, cc[800];
    static integer ng;
    static real_number hu;
    static integer jbg, nli, neq, nni, nre, npe, nxg, nyg, nps, nrt, nst, nqu;
    extern /* Subroutine */ int avc1_(real_number *, real_number *);
    static integer idid, ncfl, ncfn, info[20], ipar[2], meth;
    static real_number atol;
    static integer jpre;
    static real_number rpar[800];
    static integer nrte;
    static real_number rtol;
    static integer iout, nout;
    static real_number tout, c1ave;
    static integer imod3;
    extern /* Subroutine */ int jacrs_();
    extern /* Subroutine */ int cinit_(real_number *, real_number *, real_number
	    *, real_number *), setid_(integer *, integer *, integer *, integer
	    *, integer *, integer *);
    static real_number avlin;
    static integer leniw;
    extern /* Subroutine */ int rtweb_();
    static integer lenrw, iwork[1640], jroot;
    static real_number rwork[104883];
    extern /* Subroutine */ int _daskr_dgset2_(integer *, integer *, integer *,
	    integer *, integer *, integer *, integer *, integer *), _daskr_dmset2_(
	    integer *, integer *, integer *, integer *, integer *, integer *);
    static integer nlidif;
    static real_number predic;
    static integer nnidif;
    extern /* Subroutine */ int _daskr_ddaskr_(Unknown_fp, integer *, real_number *,
	    real_number *, real_number *, real_number *, integer *, real_number *,
	     real_number *, integer *, real_number *, integer *, integer *,
	    integer *, real_number *, integer *, Unknown_fp, Unknown_fp, Unknown_fp, integer *,
	    integer *);
    extern /* Subroutine */ int resweb_();
    extern /* Subroutine */ int setpar_(void), outweb_(real_number *,
	    real_number *, integer *, integer *, integer *, integer *, FILE*);
    extern /* Subroutine */ int psolrs_();
    static real_number ccprime[800];

    FILE* outFilewc, *outFilewd;


/* Dimension solution arrays and work arrays. */

/* When INFO(12) = 0, with INFO(5) = 0, INFO(6) = 1: */
/*    The length required for RWORK is */
/*        60 + 3*NRT + (2*ML+MU+11)*NEQ + 2*(NEQ/(ML+MU+1) + 1) . */
/*    For MX = MY = (even number) and ML = MU = NS*MX, this length is */
/*        60 + 3*NRT + (3*NS*MX + 11)*NEQ + MY . */
/*    The length required for IWORK is  40 + 2*NEQ . */

/* When INFO(12) = 1: */
/*    The length required for RWORK is */
/*        101 + 3*NRT + 19*NEQ + LENWP = 104 + 19*NEQ + NS*NS*NGRP . */
/*    The length required for IWORK is */
/*         40 + NEQ + LENIWP = 40 + NEQ + NS*NGRP . */

/* The dimensions for the various arrays are set below using parameters */
/*   MAXN    which must be .ge. NEQ = NS*MX*MY, */
/*   MAXS    which must be .ge. NS, */
/*   MAXM    which must be .ge. MAX(MX,MY). */



/* The COMMON blocks /PPAR1/ and /PPAR2/ contain problem parameters. */


/* Set output unit numbers for main output and tabulated solution. */

/* Open output files. */
	outFilewc = fopen("wccout", "w");
	outFilewd = fopen("wdout", "w");

/* Call SETPAR to set basic problem parameters. */
    setpar_();

/* Set remaining problem parameters. */
    neq = ppar2_1.ns * ppar2_1.mx * ppar2_1.my;
    ppar2_1.mxns = ppar2_1.mx * ppar2_1.ns;
    ppar2_1.dx = ppar2_1.ax / (real_number) (ppar2_1.mx - 1);
    ppar2_1.dy = ppar2_1.ay / (real_number) (ppar2_1.my - 1);
    i__1 = ppar2_1.ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = ppar2_1.dx;
	ppar2_1.cox[i__ - 1] = ppar2_1.diff[i__ - 1] / (d__1 * d__1);
/* L10: */
/* Computing 2nd power */
	d__1 = ppar2_1.dy;
	ppar2_1.coy[i__ - 1] = ppar2_1.diff[i__ - 1] / (d__1 * d__1);
    }

/* Set NRT = number of root functions. */
    nrt = 1;

    fprintf(outFilewd, fmt_20, ppar2_1.ns);

    fprintf(outFilewd, fmt_25, ppar1_1.aa, ppar1_1.ee, ppar1_1.gg, ppar1_1.bb, ppar1_1.dprey,
    		ppar1_1.dpred, ppar2_1.alph, ppar2_1.beta);

    fprintf(outFilewd, fmt_30, ppar2_1.mx, ppar2_1.my, neq);

    fputs(fmt_35, outFilewd);

/* Here set the flat initial guess for the predators. */
    predic = 1e5;

/* Set remaining method parameters for DDASKR. */
/* These include the INFO array and tolerances. */

    for (i__ = 1; i__ <= 20; ++i__) {
/* L50: */
	info[i__ - 1] = 0;
    }

/* Here set INFO(11) = 1, indicating I.C. calculation requested. */
    info[10] = 1;

/* Here set INFO(14) = 1 to get the computed initial values. */
    info[13] = 1;

/* Here set INFO(15) = 1 to signal that a preconditioner setup routine */
/* is to be called in the Krylov case. */
    info[14] = 1;

/* Here set INFO(16) = 1 to get alternative error test (on the */
/* differential variables only). */
    info[15] = 1;

/* Here set the tolerances. */
    rtol = 1e-5;
    atol = rtol;

    fprintf(outFilewd, fmt_70, rtol, atol, info[10], predic, info[15]);

/* Set NOUT = number of output times. */
    nout = 18;

/* Loop over method options: */
/* METH = 0 means use INFO(12) = 0 (direct) */
/* METH = 1 means use INFO(12) = 1 (Krylov) without block-grouping in */
/*          the reaction-based factor in the preconditioner. */
/* METH = 2 means use INFO(12) = 1 (Krylov) with block-grouping in */
/*          the reaction-based factor in the preconditioner. */
/* A block-grouping flag JBG, communicated through IPAR, is set to */
/* 0 (no block-grouping) or 1 (use block-grouping) with METH = 1 or 2. */
/* Reset INFO(1) = 0 and INFO(11) = 1. */

    for (meth = 0; meth <= 2; ++meth) {
		info[11] = MIN(meth,1);
		info[0] = 0;
		info[10] = 1;
		jbg = meth - 1;
		ipar[1] = jbg;

		fprintf(outFilewd, fmt_80, info[11]);

/* In the case of the direct method, set INFO(6) = 1 to signal a banded */
/* Jacobian, set IWORK(1) = IWORK(2) = MX*NS, the half-bandwidth, and */
/* call SETID to set the IWORK segment ID indicating the differential */
/* and algebraic components. */
		if (info[11] == 0) {
			info[5] = 1;
			iwork[0] = ppar2_1.mxns;
			iwork[1] = ppar2_1.mxns;
			setid_(&ppar2_1.mx, &ppar2_1.my, &ppar2_1.ns, &ppar2_1.np, &c__40,
				 iwork);
			fprintf(outFilewd, fmt_90, ppar2_1.mxns);
		}

/* In the case of the Krylov method, set and print various */
/* preconditioner parameters. */
		if (info[11] == 1) {
/* First set the preconditioner choice JPRE. */
/*  JPRE = 1 means reaction-only (block-diagonal) factor A_R */
/*  JPRE = 2 means spatial factor (Gauss-Seidel) A_S */
/*  JPRE = 3 means A_S * A_R */
/*  JPRE = 4 means A_R * A_S */
/* Use IPAR to communicate JPRE to the preconditioner solve routine. */
			jpre = 3;
			ipar[0] = jpre;
			fprintf(outFilewd, fmt_100, jpre);

/* Here call DMSET2 if JBG = 0, or DGSET2 if JBG = 1, to set the 2-D */
/* mesh parameters and block-grouping data, and the IWORK segment ID */
/* indicating the differential and algebraic components. */
			if (jbg == 0) {
				_daskr_dmset2_(&ppar2_1.mx, &ppar2_1.my, &ppar2_1.ns, &ppar2_1.np, &
						c__40, iwork);
				fputs(fmt_110, outFilewd);
			}
			if (jbg == 1) {
				nxg = 5;
				nyg = 5;
				ng = nxg * nyg;
				_daskr_dgset2_(&ppar2_1.mx, &ppar2_1.my, &ppar2_1.ns, &ppar2_1.np, &
					nxg, &nyg, &c__40, iwork);
				fprintf(outFilewd, fmt_120, ng, nxg, nyg);
			}
		}

/* Set the initial T and TOUT, and call CINIT to set initial values. */
		t = 0.;
		tout = 1e-8;
		cinit_(cc, ccprime, &predic, rpar);

		nli = 0;
		nni = 0;

		fprintf(outFilewd, fmt_140);

/* Loop over output times and call DDASKR.  At each output time, */
/* print average c1 value and performance data. */
/* The first call, with IOUT = 0, is to calculate initial values only. */
/* After the first call, reset INFO(11) = 0 and the initial TOUT. */
/* If a root was found, we flag this, and return to the DDASKR call. */

		i__1 = nout;
		for (iout = 0; iout <= i__1; ++iout) {

	L150:
			_daskr_ddaskr_((Unknown_fp)resweb_, &neq, &t, cc, ccprime, &tout, info, &rtol,
				&atol, &idid, rwork, &c_b58, iwork, &c__1640, rpar, ipar,
				(Unknown_fp)jacrs_, (Unknown_fp)psolrs_, (Unknown_fp)rtweb_, &nrt, &jroot);

			nst = iwork[10];
			nre = iwork[11];
			npe = iwork[12];
			nnidif = iwork[18] - nni;
			nni = iwork[18];
			nlidif = iwork[19] - nli;
			nli = iwork[19];
			nqu = iwork[7];
			hu = rwork[6];
			avlin = 0.;
			if (nnidif > 0) {
			avlin = (real_number) nlidif / (real_number) nnidif;
			}

			if (meth == 0) {
			imod3 = iout - iout / 3 * 3;
			if (imod3 == 0) {
				outweb_(&t, cc, &ppar2_1.ns, &ppar2_1.mx, &ppar2_1.my, &
					lcout, outFilewc);
			}
			}

			avc1_(cc, &c1ave);
			fprintf(outFilewd, fmt_160, t, c1ave, nst, nre, nni, nli, npe, nqu, hu, avlin);

			if (idid == 5) {
				fprintf(outFilewd, fmt_165, jroot);
			goto L150;
			}
			if (idid < 0) {
				fprintf(outFilewd, fmt_170, t);
			goto L210;
			}

			if (tout > .9) {
				tout += 1.;
			}
			if (tout < .9) {
				tout *= 10.;
			}
			if (iout == 0) {
				info[10] = 0;
				tout = 1e-8;
				nli = 0;
				nni = 0;
			}
/* L200: */
		}

L210:
		lenrw = iwork[17];
		leniw = iwork[16];
		nst = iwork[10];
		nre = iwork[11];
		npe = iwork[12];
		nni = iwork[18];
		nli = iwork[19];
		nps = iwork[20];
		if (nni > 0) {
			avlin = (real_number) nli / (real_number) nni;
		}
		ncfn = iwork[14];
		ncfl = iwork[15];
		nrte = iwork[35];
		fprintf(outFilewd, fmt_220, lenrw, leniw, nst, nre, nrte, npe, nps, nni, nli, avlin, ncfn, ncfl);

/* L300: */
    }
    exit(0);
/* ------  End of main program for DWEB example program ------------------ */
    return 0;
} /* MAIN__ */

/* Subroutine */ int setpar_(void)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static real_number pi;

/* ----------------------------------------------------------------------- */
/* This routine sets the basic problem parameters, namely */
/* AX, AY, NS, MX, MY,  problem coefficients ACOEF, BCOEF, DIFF, */
/* ALPH, BETA, using parameters NP, AA, EE, GG, BB, DPREY, DPRED. */
/* ----------------------------------------------------------------------- */

    ppar2_1.ax = 1.;
    ppar2_1.ay = 1.;
    ppar2_1.np = 1;
    ppar2_1.mx = 20;
    ppar2_1.my = 20;
    ppar1_1.aa = 1.;
    ppar1_1.ee = 1e4;
    ppar1_1.gg = 5e-7;
    ppar1_1.bb = 1.;
    ppar1_1.dprey = 1.;
    ppar1_1.dpred = .05;
    ppar2_1.alph = 50.;
    ppar2_1.beta = 100.;
    ppar2_1.ns = ppar2_1.np << 1;
    i__1 = ppar2_1.np;
    for (j = 1; j <= i__1; ++j) {
	i__2 = ppar2_1.np;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ppar2_1.acoef[ppar2_1.np + i__ + (j << 1) - 3] = ppar1_1.ee;
	    ppar2_1.acoef[i__ + (ppar2_1.np + j << 1) - 3] = -ppar1_1.gg;
/* L10: */
	}
	ppar2_1.acoef[j + (j << 1) - 3] = -ppar1_1.aa;
	ppar2_1.acoef[ppar2_1.np + j + (ppar2_1.np + j << 1) - 3] = 
		-ppar1_1.aa;
	ppar2_1.bcoef[j - 1] = ppar1_1.bb;
	ppar2_1.bcoef[ppar2_1.np + j - 1] = -ppar1_1.bb;
	ppar2_1.diff[j - 1] = ppar1_1.dprey;
	ppar2_1.diff[ppar2_1.np + j - 1] = ppar1_1.dpred;
/* L20: */
    }
    pi = 3.141592653589793;
    ppar2_1.fpi = pi * 4.;

    return 0;
/* ------------  End of Subroutine SETPAR  ------------------------------- */
} /* setpar_ */

/* Subroutine */ int setid_(integer *mx, integer *my, integer *ns, integer *
	nsd, integer *lid, integer *iwork)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, i0, i00, jx, jy, nsdp1;

/* ----------------------------------------------------------------------- */
/* This routine sets the ID array in IWORK, indicating which components */
/* are differential and which are algebraic. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --iwork;

    /* Function Body */
    nsdp1 = *nsd + 1;
    i__1 = *my;
    for (jy = 1; jy <= i__1; ++jy) {
	i00 = *mx * *ns * (jy - 1) + *lid;
	i__2 = *mx;
	for (jx = 1; jx <= i__2; ++jx) {
	    i0 = i00 + *ns * (jx - 1);
	    i__3 = *nsd;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/* L10: */
		iwork[i0 + i__] = 1;
	    }
	    i__3 = *ns;
	    for (i__ = nsdp1; i__ <= i__3; ++i__) {
/* L20: */
		iwork[i0 + i__] = -1;
	    }
/* L30: */
	}
/* L40: */
    }

    return 0;
/* ------------  End of Subroutine SETID  -------------------------------- */
} /* setid_ */

/* Subroutine */ int cinit_(real_number *cc, real_number *ccprime, real_number *
	predic, real_number *rpar)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__;
    static real_number t, x, y;
    static integer jx, jy;
    static real_number fac;
    static integer npp1, ioff;
    extern /* Subroutine */ int fweb_(real_number *, real_number *, real_number *
	    , real_number *);
    static real_number argx, argy;
    static integer iyoff;

/* ----------------------------------------------------------------------- */
/* This routine computes and loads the vectors of initial values. */
/* ----------------------------------------------------------------------- */

/* Load CC. */
    /* Parameter adjustments */
    --rpar;
    --ccprime;
    --cc;

    /* Function Body */
    npp1 = ppar2_1.np + 1;
    i__1 = ppar2_1.my;
    for (jy = 1; jy <= i__1; ++jy) {
	y = (real_number) (jy - 1) * ppar2_1.dy;
	argy = y * 16. * y * (ppar2_1.ay - y) * (ppar2_1.ay - y);
	iyoff = ppar2_1.mxns * (jy - 1);
	i__2 = ppar2_1.mx;
	for (jx = 1; jx <= i__2; ++jx) {
	    x = (real_number) (jx - 1) * ppar2_1.dx;
	    argx = x * 16. * x * (ppar2_1.ax - x) * (ppar2_1.ax - x);
	    ioff = iyoff + ppar2_1.ns * (jx - 1);
	    fac = ppar2_1.alph * x * y + 1. + ppar2_1.beta * sin(ppar2_1.fpi *
		     x) * sin(ppar2_1.fpi * y);
	    i__3 = ppar2_1.np;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/* L10: */
		cc[ioff + i__] = (real_number) i__ * argx * argy + 10.;
	    }
	    i__3 = ppar2_1.ns;
	    for (i__ = npp1; i__ <= i__3; ++i__) {
/* L15: */
		cc[ioff + i__] = *predic;
	    }
/* L20: */
	}
/* L30: */
    }

/* Load CCPRIME. */
    t = 0.;
    fweb_(&t, &cc[1], &ccprime[1], &rpar[1]);
    i__1 = ppar2_1.my;
    for (jy = 1; jy <= i__1; ++jy) {
	iyoff = ppar2_1.mxns * (jy - 1);
	i__2 = ppar2_1.mx;
	for (jx = 1; jx <= i__2; ++jx) {
	    ioff = iyoff + ppar2_1.ns * (jx - 1);
	    i__3 = ppar2_1.ns;
	    for (i__ = npp1; i__ <= i__3; ++i__) {
/* L40: */
		ccprime[ioff + i__] = 0.;
	    }
/* L50: */
	}
/* L60: */
    }

    return 0;
/* ------------  End of Subroutine CINIT  -------------------------------- */
} /* cinit_ */

/* Subroutine */ int outweb_(real_number *t, real_number *c__, integer *ns,
	integer *mx, integer *my, integer *lun, FILE* outFile)
{
    /* Format strings */
    static char fmt_10[] =
    	" -------------------------------------------------------------------------------\n"
    	"                        At time t = %16.E\n"
    	" -------------------------------------------------------------------------------\n";
    static char fmt_20[] = " the species c(%2d) values are \n";
    static char fmt_25[] = " %12.6g ";
    static char fmt_35[] =
    	" -------------------------------------------------------------------------------\n";

    /* System generated locals */
    integer c_dim1, c_dim2, c_offset, i__1, i__2;

    /* Local variables */
    static integer i__, jx, jy;


/* ----------------------------------------------------------------------- */
/* This routine prints the values of the individual species densities */
/* at the current time T, to logical unit LUN. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    c_dim1 = *ns;
    c_dim2 = *mx;
    c_offset = 1 + c_dim1 * (1 + c_dim2);
    c__ -= c_offset;

    /* Function Body */
    fprintf(outFile, fmt_10, *t);

    i__1 = *ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
    	fprintf(outFile, fmt_20, i__);
		for (jy = *my; jy >= 1; --jy) {
			i__2 = *mx;
			for (jx = 1; jx <= i__2; ++jx) {
				fprintf(outFile, fmt_25, c__[i__ + (jx + jy * c_dim2) * c_dim1]);
				if (jx % 6 == 0){
					fprintf(outFile, "\n");
				}
			}
			fprintf(outFile, "\n");
/* L30: */
		}
		fputs(fmt_35, outFile);
/* L40: */
    }

    return 0;
/* ------------  End of Subroutine OUTWEB  ------------------------------- */
} /* outweb_ */

/* Subroutine */ int resweb_(real_number *t, real_number *u, real_number *uprime,
	 real_number *cj, real_number *delta, integer *ires, real_number *rpar,
	integer *ipar)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, jx, jy, ic0, ici;
    extern /* Subroutine */ int fweb_(real_number *, real_number *, real_number *
	    , real_number *);
    static integer iyoff;

/* ----------------------------------------------------------------------- */
/* This routine computes the residual vector, using Subroutine FWEB */
/* for the right-hand sides. */
/* ----------------------------------------------------------------------- */


    /* Parameter adjustments */
    --ipar;
    --rpar;
    --delta;
    --uprime;
    --u;

    /* Function Body */
    fweb_(t, &u[1], &delta[1], &rpar[1]);

    i__1 = ppar2_1.my;
    for (jy = 1; jy <= i__1; ++jy) {
	iyoff = ppar2_1.mxns * (jy - 1);
	i__2 = ppar2_1.mx;
	for (jx = 1; jx <= i__2; ++jx) {
	    ic0 = iyoff + ppar2_1.ns * (jx - 1);
	    i__3 = ppar2_1.ns;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		ici = ic0 + i__;
		if (i__ > ppar2_1.np) {
		    delta[ici] = -delta[ici];
		} else {
		    delta[ici] = uprime[ici] - delta[ici];
		}
/* L10: */
	    }
/* L20: */
	}
/* L30: */
    }

    return 0;
/* ------------  End of Subroutine RESWEB  ------------------------------- */
} /* resweb_ */

/* Subroutine */ int fweb_(real_number *t, real_number *cc, real_number *crate,
	real_number *rpar)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, ic, jx, jy, ici;
    extern /* Subroutine */ int webr_(real_number *, integer *, integer *,
	    real_number *, real_number *);
    static integer idxl, idyl, idxu, idyu;
    static real_number dcxli, dcyli;
    static integer iyoff;
    static real_number dcyui, dcxui;

/* ----------------------------------------------------------------------- */
/* This routine computes the right-hand sides of all the equations */
/* and returns them in the array CRATE. */
/* The interaction rates are computed by calls to WEBR, and these are */
/* saved in RPAR(1),...,RPAR(NEQ) for use later in preconditioning. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --rpar;
    --crate;
    --cc;

    /* Function Body */
    i__1 = ppar2_1.my;
    for (jy = 1; jy <= i__1; ++jy) {
	iyoff = ppar2_1.mxns * (jy - 1);
	idyu = ppar2_1.mxns;
	if (jy == ppar2_1.my) {
	    idyu = -ppar2_1.mxns;
	}
	idyl = ppar2_1.mxns;
	if (jy == 1) {
	    idyl = -ppar2_1.mxns;
	}
	i__2 = ppar2_1.mx;
	for (jx = 1; jx <= i__2; ++jx) {
	    ic = iyoff + ppar2_1.ns * (jx - 1) + 1;
/* Get interaction rates at one point (X,Y). */
	    webr_(t, &jx, &jy, &cc[ic], &rpar[ic]);
	    idxu = ppar2_1.ns;
	    if (jx == ppar2_1.mx) {
		idxu = -ppar2_1.ns;
	    }
	    idxl = ppar2_1.ns;
	    if (jx == 1) {
		idxl = -ppar2_1.ns;
	    }
	    i__3 = ppar2_1.ns;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		ici = ic + i__ - 1;
/* Do differencing in Y. */
		dcyli = cc[ici] - cc[ici - idyl];
		dcyui = cc[ici + idyu] - cc[ici];
/* Do differencing in X. */
		dcxli = cc[ici] - cc[ici - idxl];
		dcxui = cc[ici + idxu] - cc[ici];
/* Collect terms and load CRATE elements. */
		crate[ici] = ppar2_1.coy[i__ - 1] * (dcyui - dcyli) + 
			ppar2_1.cox[i__ - 1] * (dcxui - dcxli) + rpar[ici];
/* L20: */
	    }
/* L40: */
	}
/* L60: */
    }
    return 0;
/* ------------  End of Subroutine FWEB  --------------------------------- */
} /* fweb_ */

/* Subroutine */ int webr_(real_number *t, integer *jx, integer *jy,
	real_number *c__, real_number *crate)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sin(real_number);

    /* Local variables */
    static integer i__, j;
    static real_number x, y, fac;
    extern /* Subroutine */ int _daskr_daxpy_(integer *, real_number *, real_number *,
	    integer *, real_number *, integer *);

/* ----------------------------------------------------------------------- */
/* This routine computes one block of the interaction term R of the */
/* system, namely block (JX,JY), for use in preconditioning. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --crate;
    --c__;

    /* Function Body */
    y = (real_number) (*jy - 1) * ppar2_1.dy;
    x = (real_number) (*jx - 1) * ppar2_1.dx;
    i__1 = ppar2_1.ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	crate[i__] = 0.;
    }
    i__1 = ppar2_1.ns;
    for (j = 1; j <= i__1; ++j) {
    	_daskr_daxpy_(&ppar2_1.ns, &c__[j], &ppar2_1.acoef[(j << 1) - 2], &c__1, &
		crate[1], &c__1);
/* L15: */
    }
    fac = ppar2_1.alph * x * y + 1. + ppar2_1.beta * sin(ppar2_1.fpi * x) * 
	    sin(ppar2_1.fpi * y);
    i__1 = ppar2_1.ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
	crate[i__] = c__[i__] * (ppar2_1.bcoef[i__ - 1] * fac + crate[i__]);
    }
    return 0;
/* ------------  End of Subroutine WEBR  --------------------------------- */
} /* webr_ */

/* Subroutine */ int jacrs_(real_number *res, integer *ires, integer *neq,
	real_number *t, real_number *cc, real_number *ccprime, real_number *rewt,
	real_number *savr, real_number *wk, real_number *h__, real_number *cj,
	real_number *wp, integer *iwp, integer *ier, real_number *rpar, integer
	*ipar)
{
    static integer jbg;
    extern /* Subroutine */ int webr_(real_number *, integer *, integer *,
	    real_number *, real_number *), _daskr_drbdja_(real_number *, real_number *,
	    real_number *, Unknown_fp, real_number *, real_number *, real_number *,
	    real_number *, integer *, integer *), _daskr_drbgja_(real_number *,
	    real_number *, real_number *, Unknown_fp, real_number *, real_number *,
	    real_number *, real_number *, integer *, integer *);

/* ----------------------------------------------------------------------- */
/* This routine interfaces to Subroutine DRBDJA or Subroutine DRBGJA, */
/* depending on the flag JBG = IPAR(2), to generate and preprocess the */
/* block-diagonal Jacobian corresponding to the reaction term R. */
/* If JBG = 0, we call DRBDJA, with no block-grouping. */
/* If JBG = 1, we call DRBGJA, and use block-grouping. */
/* Array RPAR, containing the current R vector, is passed to DRBDJA and */
/* DRBGJA as argument R0, consistent with the loading of RPAR in FWEB. */
/* The external name WEBR is passed, as the routine which computes the */
/* individual blocks of the R vector. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --ipar;
    --rpar;
    --iwp;
    --wp;
    --wk;
    --savr;
    --rewt;
    --ccprime;
    --cc;

    /* Function Body */
    jbg = ipar[2];
    if (jbg == 0) {
    	_daskr_drbdja_(t, &cc[1], &rpar[1], (Unknown_fp)webr_, &wk[1], &rewt[1], cj, &wp[1]
		, &iwp[1], ier);
    } else {
    	_daskr_drbgja_(t, &cc[1], &rpar[1], (Unknown_fp)webr_, &wk[1], &rewt[1], cj, &wp[1]
		, &iwp[1], ier);
    }

    return 0;
/* ------------  End of Subroutine JACRS  -------------------------------- */
} /* jacrs_ */

/* Subroutine */ int psolrs_(integer *neq, real_number *t, real_number *cc,
	real_number *ccprime, real_number *savr, real_number *wk, real_number *cj,
	 real_number *wt, real_number *wp, integer *iwp, real_number *b,
	real_number *eplin, integer *ier, real_number *rpar, integer *ipar)
{
    extern /* Subroutine */ int gs_(integer *, real_number *, real_number *,
	    real_number *);
    static real_number hl0;
    static integer jbg, jpre;
    extern /* Subroutine */ int _daskr_drbdps_(real_number *, real_number *, integer *)
	    , _daskr_drbgps_(real_number *, real_number *, integer *);

/* ----------------------------------------------------------------------- */
/* This routine applies the inverse of a product preconditioner matrix */
/* to the vector in the array B.  Depending on the flag JPRE, this */
/* involves a call to GS, for the inverse of the spatial factor, and/or */
/* a call to DRBDPS or DRBGPS for the inverse of the reaction-based */
/* factor (CJ*I_d - dR/dy).  The latter factor uses block-grouping */
/* (with a call to DRBGPS) if JBG = 1, and does not (with a call to */
/* DRBDPS) if JBG = 0.  JBG is communicated as IPAR(2). */
/* The array B is overwritten with the solution. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --ipar;
    --rpar;
    --b;
    --iwp;
    --wp;
    --wk;
    --savr;
    --ccprime;
    --cc;

    /* Function Body */
    jpre = ipar[1];
    *ier = 0;
    hl0 = 1. / *cj;

    jbg = ipar[2];

    if (jpre == 2 || jpre == 3) {
	gs_(neq, &hl0, &b[1], &wk[1]);
    }

    if (jpre != 2) {
	if (jbg == 0) {
		_daskr_drbdps_(&b[1], &wp[1], &iwp[1]);
	}
	if (jbg == 1) {
		_daskr_drbgps_(&b[1], &wp[1], &iwp[1]);
	}
    }

    if (jpre == 4) {
	gs_(neq, &hl0, &b[1], &wk[1]);
    }

    return 0;
/* ------------  End of Subroutine PSOLRS  ------------------------------- */
} /* psolrs_ */

/* Subroutine */ int gs_(integer *n, real_number *hl0, real_number *z__,
	real_number *x)
{
    /* Initialized data */

    static integer itmax = 5;

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, ic, ii, jx, jy, ici;
    static real_number dinv[2];
    static integer iter;
    static real_number beta1[2], beta2[2];
    static integer iyoff;
    static real_number gamma1[2], gamma2[2], elamda;

/* ----------------------------------------------------------------------- */
/* This routine provides the inverse of the spatial factor for a */
/* product preconditoner in an NS-species reaction-diffusion problem. */
/* It performs ITMAX = 5 Gauss-Seidel iterations to compute an */
/* approximation to (A_S)-inverse * Z, where A_S = I - hl0*Jd, and Jd */
/* represents the diffusion contributions to the Jacobian. */
/* The solution is vector returned in Z. */
/* ----------------------------------------------------------------------- */


    /* Parameter adjustments */
    --x;
    --z__;

    /* Function Body */
/* ----------------------------------------------------------------------- */
/* Write matrix as A = D - L - U. */
/* Load local arrays BETA, BETA2, GAMMA1, GAMMA2, and DINV. */
/* ----------------------------------------------------------------------- */
    i__1 = ppar2_1.ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	elamda = 1. / (*hl0 * 2. * (ppar2_1.cox[i__ - 1] + ppar2_1.coy[i__ - 
		1]) + 1.);
	beta1[i__ - 1] = *hl0 * ppar2_1.cox[i__ - 1] * elamda;
	beta2[i__ - 1] = beta1[i__ - 1] * 2.;
	gamma1[i__ - 1] = *hl0 * ppar2_1.coy[i__ - 1] * elamda;
	gamma2[i__ - 1] = gamma1[i__ - 1] * 2.;
	dinv[i__ - 1] = elamda;
/* L10: */
    }
/* ----------------------------------------------------------------------- */
/* Begin iteration loop. */
/* Load array X with (D-inverse)*Z for first iteration. */
/* ----------------------------------------------------------------------- */
    iter = 1;
/* Zero X in all its components, since X is added to Z at the end. */
    i__1 = *n;
    for (ii = 1; ii <= i__1; ++ii) {
/* L15: */
	x[ii] = 0.;
    }

    i__1 = ppar2_1.my;
    for (jy = 1; jy <= i__1; ++jy) {
	iyoff = ppar2_1.mxns * (jy - 1);
	i__2 = ppar2_1.mx;
	for (jx = 1; jx <= i__2; ++jx) {
	    ic = iyoff + ppar2_1.ns * (jx - 1);
	    i__3 = ppar2_1.ns;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		ici = ic + i__;
		x[ici] = dinv[i__ - 1] * z__[ici];
		z__[ici] = 0.;
/* L30: */
	    }
/* L40: */
	}
/* L50: */
    }
    goto L160;
/* ----------------------------------------------------------------------- */
/* Calculate (D-inverse)*U*X. */
/* ----------------------------------------------------------------------- */
L70:
    ++iter;
    jy = 1;
    jx = 1;
    ic = ppar2_1.ns * (jx - 1);
    i__1 = ppar2_1.ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ici = ic + i__;
/* L75: */
	x[ici] = beta2[i__ - 1] * x[ici + ppar2_1.ns] + gamma2[i__ - 1] * x[
		ici + ppar2_1.mxns];
    }
    i__1 = ppar2_1.mx - 1;
    for (jx = 2; jx <= i__1; ++jx) {
	ic = ppar2_1.ns * (jx - 1);
	i__2 = ppar2_1.ns;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ici = ic + i__;
/* L80: */
	    x[ici] = beta1[i__ - 1] * x[ici + ppar2_1.ns] + gamma2[i__ - 1] * 
		    x[ici + ppar2_1.mxns];
	}
/* L85: */
    }
    jx = ppar2_1.mx;
    ic = ppar2_1.ns * (jx - 1);
    i__1 = ppar2_1.ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ici = ic + i__;
/* L90: */
	x[ici] = gamma2[i__ - 1] * x[ici + ppar2_1.mxns];
    }
    i__1 = ppar2_1.my - 1;
    for (jy = 2; jy <= i__1; ++jy) {
	iyoff = ppar2_1.mxns * (jy - 1);
	jx = 1;
	ic = iyoff;
	i__2 = ppar2_1.ns;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ici = ic + i__;
/* L95: */
	    x[ici] = beta2[i__ - 1] * x[ici + ppar2_1.ns] + gamma1[i__ - 1] * 
		    x[ici + ppar2_1.mxns];
	}
	i__2 = ppar2_1.mx - 1;
	for (jx = 2; jx <= i__2; ++jx) {
	    ic = iyoff + ppar2_1.ns * (jx - 1);
	    i__3 = ppar2_1.ns;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		ici = ic + i__;
/* L100: */
		x[ici] = beta1[i__ - 1] * x[ici + ppar2_1.ns] + gamma1[i__ - 
			1] * x[ici + ppar2_1.mxns];
	    }
/* L105: */
	}
	jx = ppar2_1.mx;
	ic = iyoff + ppar2_1.ns * (jx - 1);
	i__2 = ppar2_1.ns;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ici = ic + i__;
/* L110: */
	    x[ici] = gamma1[i__ - 1] * x[ici + ppar2_1.mxns];
	}
/* L115: */
    }
    jy = ppar2_1.my;
    iyoff = ppar2_1.mxns * (jy - 1);
    jx = 1;
    ic = iyoff;
    i__1 = ppar2_1.ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ici = ic + i__;
/* L120: */
	x[ici] = beta2[i__ - 1] * x[ici + ppar2_1.ns];
    }
    i__1 = ppar2_1.mx - 1;
    for (jx = 2; jx <= i__1; ++jx) {
	ic = iyoff + ppar2_1.ns * (jx - 1);
	i__2 = ppar2_1.ns;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ici = ic + i__;
/* L125: */
	    x[ici] = beta1[i__ - 1] * x[ici + ppar2_1.ns];
	}
/* L130: */
    }
    jx = ppar2_1.mx;
    ic = iyoff + ppar2_1.ns * (jx - 1);
    i__1 = ppar2_1.ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ici = ic + i__;
/* L135: */
	x[ici] = 0.;
    }
/* ----------------------------------------------------------------------- */
/* Calculate [(I - (D-inverse)*L)]-inverse * X. */
/* ----------------------------------------------------------------------- */
L160:
    jy = 1;
    i__1 = ppar2_1.mx - 1;
    for (jx = 2; jx <= i__1; ++jx) {
	ic = ppar2_1.ns * (jx - 1);
	i__2 = ppar2_1.ns;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ici = ic + i__;
/* L170: */
	    x[ici] += beta1[i__ - 1] * x[ici - ppar2_1.ns];
	}
/* L175: */
    }
    jx = ppar2_1.mx;
    ic = ppar2_1.ns * (jx - 1);
    i__1 = ppar2_1.ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ici = ic + i__;
/* L180: */
	x[ici] += beta2[i__ - 1] * x[ici - ppar2_1.ns];
    }
    i__1 = ppar2_1.my - 1;
    for (jy = 2; jy <= i__1; ++jy) {
	iyoff = ppar2_1.mxns * (jy - 1);
	jx = 1;
	ic = iyoff;
	i__2 = ppar2_1.ns;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ici = ic + i__;
/* L185: */
	    x[ici] += gamma1[i__ - 1] * x[ici - ppar2_1.mxns];
	}
	i__2 = ppar2_1.mx - 1;
	for (jx = 2; jx <= i__2; ++jx) {
	    ic = iyoff + ppar2_1.ns * (jx - 1);
	    i__3 = ppar2_1.ns;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		ici = ic + i__;
		x[ici] = x[ici] + beta1[i__ - 1] * x[ici - ppar2_1.ns] + 
			gamma1[i__ - 1] * x[ici - ppar2_1.mxns];
/* L195: */
	    }
/* L200: */
	}
	jx = ppar2_1.mx;
	ic = iyoff + ppar2_1.ns * (jx - 1);
	i__2 = ppar2_1.ns;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ici = ic + i__;
	    x[ici] = x[ici] + beta2[i__ - 1] * x[ici - ppar2_1.ns] + gamma1[
		    i__ - 1] * x[ici - ppar2_1.mxns];
/* L205: */
	}
/* L210: */
    }
    jy = ppar2_1.my;
    iyoff = ppar2_1.mxns * (jy - 1);
    jx = 1;
    ic = iyoff;
    i__1 = ppar2_1.ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ici = ic + i__;
/* L215: */
	x[ici] += gamma2[i__ - 1] * x[ici - ppar2_1.mxns];
    }
    i__1 = ppar2_1.mx - 1;
    for (jx = 2; jx <= i__1; ++jx) {
	ic = iyoff + ppar2_1.ns * (jx - 1);
	i__2 = ppar2_1.ns;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ici = ic + i__;
	    x[ici] = x[ici] + beta1[i__ - 1] * x[ici - ppar2_1.ns] + gamma2[
		    i__ - 1] * x[ici - ppar2_1.mxns];
/* L220: */
	}
/* L225: */
    }
    jx = ppar2_1.mx;
    ic = iyoff + ppar2_1.ns * (jx - 1);
    i__1 = ppar2_1.ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ici = ic + i__;
	x[ici] = x[ici] + beta2[i__ - 1] * x[ici - ppar2_1.ns] + gamma2[i__ - 
		1] * x[ici - ppar2_1.mxns];
/* L230: */
    }
/* ----------------------------------------------------------------------- */
/* Add increment X to Z. */
/* ----------------------------------------------------------------------- */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L300: */
	z__[i__] += x[i__];
    }

    if (iter < itmax) {
	goto L70;
    }
    return 0;
/* ------------  End of Subroutine GS  ----------------------------------- */
} /* gs_ */

/* Subroutine */ int avc1_(real_number *cc, real_number *c1ave)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer jx, jy;
    static real_number sum;
    static integer npp1, ioff, iyoff;

/* ----------------------------------------------------------------------- */
/* This routine computes the average value of c1. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --cc;

    /* Function Body */
    sum = 0.;
    npp1 = ppar2_1.np + 1;
    i__1 = ppar2_1.my;
    for (jy = 1; jy <= i__1; ++jy) {
	iyoff = ppar2_1.mxns * (jy - 1);
	i__2 = ppar2_1.mx;
	for (jx = 1; jx <= i__2; ++jx) {
	    ioff = iyoff + ppar2_1.ns * (jx - 1);
	    sum += cc[ioff + 1];
/* L20: */
	}
/* L30: */
    }

    *c1ave = sum / (ppar2_1.mx * ppar2_1.my);

    return 0;
/* ------------  End of Subroutine AVC1  --------------------------------- */
} /* avc1_ */

/* Subroutine */ int rtweb_(integer *neq, real_number *t, real_number *cc,
	real_number *cp, integer *nrt, real_number *rval, real_number *rpar,
	integer *ipar)
{
    extern /* Subroutine */ int avc1_(real_number *, real_number *);
    static real_number c1ave;


/* This routine sets RVAL = average(c1) - 20.0. */


    /* Parameter adjustments */
    --cp;
    --cc;

    /* Function Body */
    avc1_(&cc[1], &c1ave);
    *rval = c1ave - 20.;

    return 0;
/* ------------  End of Subroutine RTWEB  -------------------------------- */
} /* rtweb_ */

