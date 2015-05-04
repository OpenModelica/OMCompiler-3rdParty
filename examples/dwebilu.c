/* dwebilu.f -- translated by f2c (version 20100827).
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
static integer c__20002 = 20002;
static integer c__33604 = 33604;
static integer c__2 = 2;
static integer c__40 = 40;
static integer c__35306 = 35306;
static integer c__35244 = 35244;

/* ***BEGIN PROLOGUE  DWEBILU */
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

/* The root function is R(Y) = average(c1) - 20. */

/* The DAE system is solved by DDASKR with a sparse approximate Jacobian */
/* matrix as the preconditioner for the preconditioned Krylov option of */
/* DDASKR.  The Jacobian is formed by using finite-differences of residual */
/* values, and is stored in sparse format (compressed sparse row). */
/* An incomplete factorization is then performed using the ILUT(P) routine */
/* of the SPARSKIT library of Yousef Saad at the University of Minnesota. */

/* The preconditioner is set up and solved in the three main */
/* subroutines --  DSPSETUP, DJACILU and DPSOLILU. */
/* These routines are provided separately for use on general DAE */
/* problems. */

/* Two output files are written.. one with the problem description and */
/* performance statistics on unit LOUT = 9, and one with solution */
/* profiles at selected output times on unit LCOUT = 10. */
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
/*   SETPAR, CINIT, DDASKR, OUTWEB */

/* ***END PROLOGUE  DWEBILU */

/* Main program */ int main(void)
{
    /* Initialized data */

	static integer lout = 9;
    static integer lcout = 10;

    /* Format strings */
    static char fmt_20[] = " DWEBILU: Example program for DDASKR package\n\n"
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
    static char fmt_45[] = " Error return from DSPSETUP: IERR = %5d\n";
    static char fmt_70[] =
    	" Tolerance parameters.. RTOL =%10.2E   ATOL =%10.2E\n"
    	" Internal I.C. calculation flag INFO(11) =%2d  (0 = off, 1 = on)\n"
    	" Predator I.C. guess =%10.2E\n"
    	" Alternate error test flag INFO(16) =%2d  (0 = off, 1 = on)\n";
    static char fmt_100[] = " Linear solver is: Krylov with ILU preconditioner\n"
    	" Preconditioner flag is IPREMETH =%3d  (1 = ILUT, 2 = ILUTP)\n";
    static char fmt_140[] = "\n   t             Ave.c1   NSTEP   NRE    NNI    NLI "
    		"   NPE    NQ     H          AVLIN\n";
    static char fmt_160[] = "%13.5E %10.5f %5d %6d  %5d  %5d  %5d  %4d %11.2E %9.4f\n";
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
    static char fmt_230[] = " Minimum lengths for work arrays WP and IWP: %7d %7d\n";

    /* System generated locals */
    integer i__1;
    real_number d__1;

    /* Local variables */
    extern /* Subroutine */ int _daskr_dpsolilu_();
    extern /* Subroutine */ int _daskr_dspsetup_(integer *, integer *, integer *,
	    real_number *, integer *, integer *, integer *, integer *);
    static integer i__;
    static real_number t, cc[800];
    static integer ml;
    static real_number hu;
    static integer mu, nli, neq, nni, nre, npe, nps, nrt, nst, nqu;
    extern /* Subroutine */ int avc1_(real_number *, real_number *);
    static integer idid, ncfl, ncfn, info[20], ipar[30];
    static real_number atol;
    static integer ierr;
    static real_number rpar[802];
    static integer lniw, nrte;
    static real_number rtol;
    static integer iout, lnrw, nout;
    static real_number tout, c1ave;
    static integer imod3;
    extern /* Subroutine */ int cinit_(real_number *, real_number *, real_number
	    *, real_number *), setid_(integer *, integer *, integer *, integer
	    *, integer *, integer *);
    static real_number avlin;
    extern /* Subroutine */ int rtweb_();
    static integer iwork[35244], jroot;
    static real_number rwork[35306];
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
	    real_number *, integer *, integer *, integer *, integer *, FILE *);
    static integer lwpmin;
    extern /* Subroutine */ int _daskr_xsetun_(integer *);
    extern /* Subroutine */ int _daskr_djacilu_();
    static real_number ccprime[800];
    static integer liwpmin;

    FILE* outFilewc, *outFilewd;



/* Load parameters for sparse preconditioner. */
/* =2 means ILUTP preconditioner used */
/* =1 means ILUT preconditioner used */

/* Dimension solution arrays and work arrays. */

/* When INFO(12) = 1, with INFO(5) = 0, INFO(6) = 1: */
/*    The length required for RWORK is */
/*       101 + 3*NRT + 19*NEQ + LENWP, */
/*    where LENWP is given by */
/*       LENWP =  length of RWORK segment WP = */
/*              2*LENPFAC*NEQ + LENPLUFAC*NEQ + ISRNORM*NEQ + 2*(NEQ+1) */
/*    The length required for IWORK is */
/*       40 + NEQ + LENIWP, */
/*    where LENIWP is given by */
/*       LENIWP = length of IWORK segment IWP = */
/*              4*(NEQ+1) + 3*LENPFAC*NEQ + 2*LENPLUFAC*NEQ */
/*                 + 2*IREORDER*NEQ + 2*(IPREMETH-1)*NEQ */

/* The dimensions for the various arrays are set below using parameters */
/*   MAXN    which must be .ge. NEQ = NS*MX*MY, */
/*   MAXS    which must be .ge. NS. */



/* The COMMON blocks /PPAR1/ and /PPAR2/ contain problem parameters. */


/* Set output unit numbers for main output and tabulated solution. */

/* Set unit number for DDASKR error message output. */
    _daskr_xsetun_(&lout);

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

/* Load values into IPAR and RPAR for sparse preconditioner. */
    ml = ppar2_1.ns * ppar2_1.mx + 1;
    mu = ml;
    ipar[0] = ml;
    ipar[1] = mu;
    ipar[2] = 10;
    ipar[3] = 2;
    ipar[4] = 2;
    ipar[5] = 10;
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
    _daskr_dspsetup_(&neq, &c__20002, &c__33604, rpar, ipar, &ierr, &lwpmin, &
	    liwpmin);
    if (ierr != 0) {
    	printf(fmt_45, ierr);
		if (lwpmin > 20002) {
			puts(" More WP work array length needed");
		}
		if (liwpmin > 33604) {
			puts(" More IWP work array length needed");
		}
		exit(0);
    }

/* Set remaining method parameters for DDASKR. */
/* These include the INFO array and tolerances. */

    for (i__ = 1; i__ <= 20; ++i__) {
/* L50: */
	info[i__ - 1] = 0;
    }

/* Here set INFO(11) = 1, indicating I.C. calculation requested. */
    info[10] = 1;

/* Here set INFO(12) = 1, indicating the Krylov linear system method. */
    info[11] = 1;

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

/* Set and print various preconditioner parameters. */
    iwork[26] = 20002;
    iwork[27] = 33604;
    fprintf(outFilewd, fmt_100, c__2);
/* Here call SETID to set the IWORK segment ID  indicating the */
/* differential and algebraic components. */
    setid_(&ppar2_1.mx, &ppar2_1.my, &ppar2_1.ns, &ppar2_1.np, &c__40, iwork);
/* Set the initial T and TOUT, and call CINIT to set initial values. */
    t = 0.;
    tout = 1e-8;
    cinit_(cc, ccprime, &predic, rpar);

    nli = 0;
    nni = 0;

    fputs(fmt_140, outFilewd);

/* Loop over output times and call DDASKR.  At each output time, */
/* print average c1 value and performance data. */
/* The first call, with IOUT = 0, is to calculate initial values only. */
/* After the first call, reset INFO(11) = 0 and the initial TOUT. */
/* If a root was found, we flag this, and return to the DDASKR call. */

    i__1 = nout;
    for (iout = 0; iout <= i__1; ++iout) {

L150:
    _daskr_ddaskr_((Unknown_fp)resweb_, &neq, &t, cc, ccprime, &tout, info, &rtol, &
		atol, &idid, rwork, &c__35306, iwork, &c__35244, rpar, ipar, (
		Unknown_fp)_daskr_djacilu_, (Unknown_fp)_daskr_dpsolilu_, (Unknown_fp)rtweb_, &nrt, &jroot);

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

	imod3 = iout - iout / 3 * 3;
	if (imod3 == 0) {
	    outweb_(&t, cc, &ppar2_1.ns, &ppar2_1.mx, &ppar2_1.my, &lcout, outFilewc);
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
    lnrw = iwork[17];
    lniw = iwork[16];
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
    fprintf(outFilewd, fmt_220, lnrw, lniw, nst, nre, nrte, npe, nps, nni, nli, avlin, ncfn, ncfl);
    fprintf(outFilewd, fmt_230, lwpmin, liwpmin);

    exit(0);
/* ------  End of main program for DWEBILU example program --------------- */
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

    /* Builtin functions */
    double sin(real_number);

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
    fprintf(outFile, fmt_10, *t);;

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
	    webr_(t, &jx, &jy, &cc[ic], &rpar[ic + 2]);
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
			ppar2_1.cox[i__ - 1] * (dcxui - dcxli) + rpar[ici + 2]
			;
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

