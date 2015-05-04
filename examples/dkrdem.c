/* dkrdem.f -- translated by f2c (version 20100827).
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../solver/ddaskr_types.h"

static integer global_neq;

/* DKRDEM program: DDASKR demonstration program */
/* ----------------------------------------------------------------------- */

/* ***BEGIN PROLOGUE  DKRDEM */
/* ***DATE WRITTEN   020813     (YYMMDD) */
/* ***REVISION DATE  021217   Added JROOT output value in Problem 2. */
/* ***AUTHORS  Linda R. Petzold and Alan C. Hindmarsh */
/*             LAWRENCE LIVERMORE NATIONAL LABORATORY */
/*             LIVERMORE, CA    94550 */
/* ***PURPOSE  Quick check for routine DDASKR. */
/* ***DESCRIPTION */
/*       Demonstration program for DDASKR. */
/*       This version is in double precision. */

/*       DDASKR is used to solve two simple problems, */
/*       one nonstiff and one intermittently stiff. */
/*       If the errors are too large, or other difficulty occurs, */
/*       a warning message is printed.  All output is on unit LUN. */

/*       To run the demonstration problems with full printing, */
/*       set KPRINT = 3. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  DDASKR,RES1,RT1,RES2,JAC2,RT2 */
/* ***END PROLOGUE DKRDEM */

/* ----------------------------------------------------------------------- */

/* Main program */ int main(void)
{
    /* Format strings */
    static char fmt_110[] = " DKRDEM: Demonstration Program for DDASKR\n\n\n"
    	" Problem 1..\n\n"
    	" Problem is  dY/dT = ((2*LOG(Y)+8)/T - 5)*Y,  Y(1) = 1\n"
    	" Solution is  Y(T) = EXP(-T**2 + 5*T - 4)\n"
    	" Root functions are..\n"
    	" R1 = dY/dT  (root at T = 2.5)\n"
    	" R2 = LOG(Y) - 2.2491  (roots at T = 2.47 and T = 2.53)\n"
        " RTOL =%10.1E   ATOL =%10.1E    JTYPE =%4d\n\n";
    static char fmt_130[] = " At t =%15.7E     y =%15.7E     error =%12.4E\n";
    static char fmt_135[] = "\n\n WARNING.. Error exceeds 1000 * tolerance\n\n";
    static char fmt_150[] = "\n      Root found at t =%15.7E      JROOT = %3d  %3d\n";
    static char fmt_160[] = "      Error in t location of root is %12.4E\n";
    static char fmt_165[] = "\n\n WARNING.. Root error exceeds 1.0D-3\n\n";
    static char fmt_190[] = "\n Final statistics for this run..\n"
	    " number of steps =%5d\n"
	    " number of Gs    =%5d\n"
    	" (excluding Js)  =%5d\n"
    	" number of Js    =%5d\n"
    	" number of Rs    =%5d\n";
    static char fmt_195[] = " error overrun   =%10.2E\n";
    static char fmt_200[] = "\n--------------------------------------------------------------------------------\n"
    	" Problem 2.. Van Der Pol oscillator\n"
    	" Problem is dY1/dT = Y2,  dY2/dT = 100*(1-Y1**2)*Y2 - Y1\n"
    	"            Y1(0) = 2,  Y2(0) = 0\n"
    	" Root function is R(T,Y,YP) = Y1\n"
        " RTOL =%10.1E   ATOL =%10.1E%10.1E\n\n";
    static char fmt_210[] = "\n--------------------------------------------------------------------------------\n"
        " Solution with JTYPE =%4d\n";
    static char fmt_230[] = " At t =%15.7E     y1 =%15.7E     y2 =%12.4E\n";
    static char fmt_240[] = "\n      Root found at t =%15.7E      JROOT = %3d\n";
    static char fmt_300[] = "\n--------------------------------------------------------------------------------\n"
    	"\n Number of errors encountered =%4d";
    static char fmt_700[] = "\n **********DDASKR passed all tests**********\n";
    static char fmt_800[] = "\n **********DDASKR failed some tests**********\n";

    /* Local variables */
    static integer i__;
    static real_number t, y[2], er, yt;
    extern /* Subroutine */ int rt1_(), rt2_();
    static integer nje, nre;
    static real_number ero;
    static integer liw, lun, nrt, lrw, nst;
    extern /* Subroutine */ int jac2_(), res1_(), res2_();
    static integer idid, nrea, info[20], ipar;
    static real_number atol[2];
    static integer jdum;
    static real_number rpar;
    static integer nerr, nrte;
    static real_number errt;
    static integer iout;
    static real_number rtol[2], tout;
    static integer leniw, ipass, lenrw;
    static real_number psdum;
    static integer iwork[100], jtype, jroot[2], kroot;
    static real_number tzero, rwork[100];
    extern /* Subroutine */ int _daskr_ddaskr_();
    static real_number yprime[2];
    static integer kprint;



    lun = 6;
    kprint = 3;
    ipass = 1;
    nerr = 0;
/* ----------------------------------------------------------------------- */
/* First problem. */
/* The initial value problem is.. */
/*   dY/dT = ((2*LOG(Y) + 8)/T - 5)*Y,  Y(1) = 1,  1 .LE. T .LE. 6 */
/* The solution is  Y(T) = EXP(-T**2 + 5*T - 4), YPRIME(1) = 3 */
/* The two root functions are.. */
/*   R1(T,Y,Y') = dY/dT  (with root at T = 2.5), */
/*   R2(T,Y,Y') = LOG(Y) - 2.2491  (with roots at T = 2.47 and 2.53) */
/* ----------------------------------------------------------------------- */
/* Set all input parameters and print heading. */
    for (i__ = 1; i__ <= 20; ++i__) {
/* L10: */
	info[i__ - 1] = 0;
    }
    global_neq = 1;
    t = 1.;
    y[0] = 1.;
    tout = 2.;
    rtol[0] = 0.;
    atol[0] = 1e-6;
    lrw = 100;
    liw = 100;
    idid = 0;

/* Set INFO(11) = 1 if DDASKR is to compute the initial YPRIME, and */
/* generate an initial guess for YPRIME.  Otherwise, set INFO(11) = 0 */
/* and supply the correct initial value for YPRIME. */

    info[10] = 0;
    yprime[0] = 3.;

/* Note: JTYPE indicates the Jacobian type: */
/*       JTYPE = 1 ==> Jacobian is dense and user-supplied */
/*       JTYPE = 2 ==> Jacobian is dense and computed internally */

    jtype = 2;
    info[4] = 2 - jtype;
    nrt = 2;
    if (kprint >= 2) {
	  printf(fmt_110, rtol[0], atol[0], jtype);
    }

/* Call DDASKR in loop over TOUT values = 2, 3, 4, 5, 6. */
    ero = 0.;
    for (iout = 1; iout <= 5; ++iout) {

L120:
	_daskr_ddaskr_((Unknown_fp)res1_, &global_neq, &t, y, yprime, &tout, info, rtol,
		atol, &idid, rwork, &lrw, iwork, &liw, &rpar, &ipar, &jdum, &
		psdum, (Unknown_fp)rt1_, &nrt, jroot);

/* Print Y and error in Y, and print warning if error too large. */
	yt = exp(-t * t + t * 5. - 4.);
	er = y[0] - yt;
	if (kprint > 2) {
	    printf(fmt_130, t, y[0], er);
	}
	if (idid < 0) {
	    goto L185;
	}
	er = fabs(er) / atol[0];
	ero = MAX(ero,er);
	if (er < 1e3) {
	    goto L140;
	}
	if (kprint >= 2) {
	    printf(fmt_135);
	}
	++nerr;
L140:
	if (idid != 5) {
	    goto L175;
	}

/* If a root was found, write results and check root location. */
/* Then return to DDASKR to continue the integration. */
	if (kprint > 2) {
	    printf(fmt_150, t, jroot[0], jroot[1]);
	}
	if (jroot[0] != 0) {
	    errt = t - 2.5;
	}
	if (jroot[1] != 0 && t <= 2.5) {
	    errt = t - 2.47;
	}
	if (jroot[1] != 0 && t > 2.5) {
	    errt = t - 2.53;
	}
	if (kprint > 2) {
	    printf(fmt_160,  errt);
	}
	if (fabs(errt) < .001) {
	    goto L170;
	}
	if (kprint >= 2) {
	    ipass = 0;
	    printf(fmt_165);
	}
	++nerr;
L170:
	goto L120;

/* If no root found, increment TOUT and loop back. */
L175:
	tout += 1.;
/* L180: */
    }

/* Problem complete.  Print final statistics. */
L185:
    if (idid < 0) {
	++nerr;
    }
    nst = iwork[10];
    nre = iwork[11];
    nje = iwork[12];
    nrte = iwork[35];
    lenrw = 0;
    leniw = 0;
    nrea = nre;
    if (jtype == 2) {
	nre += global_neq * nje;
    }

    if (kprint > 2) {
		printf(fmt_190, nst, nre, nrea, nje, nrte);
		printf(fmt_195, ero);
    }

/* ----------------------------------------------------------------------- */
/* Second problem (Van Der Pol oscillator). */
/* The initial value problem is.. */
/*   dY1/dT = Y2,  dY2/dT = 100*(1 - Y1**2)*Y2 - Y1, */
/*   Y1(0) = 2,  Y2(0) = 0,  0 .LE. T .LE. 200 */
/*   Y1PRIME(0) = 0, Y2PRIME(0) = -2 */
/* The root function is  R(t,Y,Y') = Y1. */
/* An analytic solution is not known, but the zeros of Y1 are known */
/* to 15 figures for purposes of checking the accuracy. */
/* ----------------------------------------------------------------------- */

/* Reset INFO array */

    for (i__ = 1; i__ <= 20; ++i__) {
/* L195: */
	info[i__ - 1] = 0;
    }

/* Set tolerance parameters and print heading. */
/* Note that INFO(2) is set to 1, indicating that RTOL and ATOL */
/* are arrays.  Each entry of RTOL and ATOL must then be defined. */

    info[1] = 1;
    rtol[0] = 1e-6;
    rtol[1] = 1e-6;
    atol[0] = 1e-6;
    atol[1] = 1e-4;
    if (kprint >= 2) {
		printf(fmt_200, rtol[0], atol[0], atol[1]);
    }

/* Note: JTYPE indicates the Jacobian type: */
/*       JTYPE = 1 ==> Jacobian is dense and user-supplied */
/*       JTYPE = 2 ==> Jacobian is dense and computed internally */

/* Loop over JTYPE = 1, 2.  Set remaining parameters and print JTYPE. */
    for (jtype = 1; jtype <= 2; ++jtype) {

/*     Set INFO(1) = 0 to indicate start of a new problem */
/*     Set INFO(5) = 2-JTYPE to tell DDASKR the Jacobian type. */

	info[0] = 0;
	info[4] = 2 - jtype;
	global_neq = 2;
	t = 0.;
	y[0] = 2.;
	y[1] = 0.;
	yprime[0] = 0.;
	yprime[1] = -2.;
	tout = 20.;
	nrt = 1;
	if (kprint > 2) {
	    printf(fmt_210, jtype);
	}

/* Call DDASKR in loop over TOUT values = 20, 40, ..., 200. */
	for (iout = 1; iout <= 10; ++iout) {

L220:
	    _daskr_ddaskr_((Unknown_fp)res2_, &global_neq, &t, y, yprime, &tout, info,
		    rtol, atol, &idid, rwork, &lrw, iwork, &liw, &rpar, &ipar,
		     (Unknown_fp)jac2_, &psdum, (Unknown_fp)rt2_, &nrt, jroot);

/* Print Y1 and Y2. */
	    if (kprint > 2) {
			printf(fmt_230, t, y[0], y[1]);
	    }
	    if (idid < 0) {
		goto L275;
	    }
	    if (idid != 5) {
		goto L265;
	    }

/* If a root was found, write results and check root location. */
/* Then return to DDASKR to continue the integration. */
	    if (kprint > 2) {
			printf(fmt_240, t, jroot[0]);
	    }
	    kroot = (integer) (t / 81.2 + .5);
	    tzero = (real_number) (kroot - 1) * 81.41853556212 +
		    81.17237787055;
	    errt = t - tzero;
	    if (kprint > 2) {
			printf(fmt_160,  errt);
	    }
	    if (errt < 1.) {
		goto L260;
	    }
	    if (kprint >= 2) {
	    	ipass = 0;
			printf(fmt_165);
	    }
	    ++nerr;
L260:
	    goto L220;

/* If no root found, increment TOUT and loop back. */
L265:
	    tout += 20.;
/* L270: */
	}

/* Problem complete.  Print final statistics. */
L275:
	if (idid < 0) {
	    ++nerr;
	}
	nst = iwork[10];
	nre = iwork[11];
	nje = iwork[12];
	nrte = iwork[35];
	lenrw = 0;
	leniw = 0;
	nrea = nre;
	if (jtype == 2) {
	    nre += global_neq * nje;
	}
	if (kprint >= 2) {
	    printf(fmt_190, nst, nre, nrea, nje, nrte);
	}
/* L290: */
    }


    if (kprint >= 2) {
		printf(fmt_300, nerr);
    }

    if (nerr > 0) {
    	ipass = 0;
    }
    if (ipass == 1 && kprint > 1) {
		printf(fmt_700, nerr);
    }
    if (ipass == 0 && kprint != 0) {
		printf(fmt_800, nerr);
    }
    exit(0);
    return 0;
} /* MAIN__ */


/* Subroutine */ int res1_(real_number *t, real_number *y, real_number *yprime,
	real_number *cj, real_number *delta, integer *ires, real_number *rpar,
	integer *ipar)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int f1_(integer *, real_number *, real_number *,
	    real_number *);


/*     Check Y to make sure that it is valid input. */
/*     If Y is less than or equal to zero, this is invalid input. */

    /* Parameter adjustments */
    --delta;
    --yprime;
    --y;

    /* Function Body */
    if (y[1] <= 0.) {
	*ires = -1;
	return 0;
    } else {

/*        Call F1 to obtain F(T,Y) */

	f1_(&global_neq, t, &y[1], &delta[1]);

/*        Form G = Y' - F(T,Y) */

	i__1 = global_neq;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    delta[i__] = yprime[i__] - delta[i__];
/* L10: */
	}
    }

    return 0;
} /* res1_ */


/* Subroutine */ int f1_(integer *global_neq, real_number *t, real_number *y,
	real_number *ydot)
{
    /* Builtin functions */
    double log(real_number);

    /* Parameter adjustments */
    --ydot;
    --y;

    /* Function Body */
    ydot[1] = ((log(y[1]) * 2. + 8.) / *t - 5.) * y[1];
    return 0;
} /* f1_ */


/* Subroutine */ int rt1_(integer *neq, real_number *t, real_number *y,
	real_number *yp, integer *nrt, real_number *rval, real_number *rpar,
	integer *ipar)
{
    /* Builtin functions */
    double log(real_number);

    /* Parameter adjustments */
    --rval;

    /* Function Body */
    rval[1] = *yp;
    rval[2] = log(*y) - 2.2491;
    return 0;
} /* rt1_ */


/* Subroutine */ int res2_(real_number *t, real_number *y, real_number *yprime,
	real_number *cj, real_number *delta, integer *ires, real_number *rpar,
	integer *ipar)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int f2_(integer *, real_number *, real_number *,
	    real_number *);


/*     CALL F2 TO OBTAIN F(T,Y) */

    /* Parameter adjustments */
    --delta;
    --yprime;
    --y;

    /* Function Body */
    f2_(&global_neq, t, &y[1], &delta[1]);

/*     FORM G = Y' - F(T,Y) */

    i__1 = global_neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	delta[i__] = yprime[i__] - delta[i__];
    }

    return 0;
} /* res2_ */


/* Subroutine */ int f2_(integer *global_neq, real_number *t, real_number *y,
	real_number *ydot)
{
    /* Parameter adjustments */
    --ydot;
    --y;

    /* Function Body */
    ydot[1] = y[2];
    ydot[2] = (1. - y[1] * y[1]) * 100. * y[2] - y[1];
    return 0;
} /* f2_ */


/* Subroutine */ int jac2_(real_number *t, real_number *y, real_number *yprime,
	real_number *pd, real_number *cj, real_number *rpar, integer *ipar)
{

/* First define the Jacobian matrix for the right-hand side */
/* of the ODE: Y' = F(T,Y) , i.e. dF/dY. */

    /* Parameter adjustments */
    pd -= 3;
    --y;

    /* Function Body */
    pd[3] = 0.;
    pd[5] = 1.;
    pd[4] = y[1] * -200. * y[2] - 1.;
    pd[6] = (1. - y[1] * y[1]) * 100.;

/* Next update the Jacobian with the right-hand side to form the */
/* DAE Jacobian: CJ*dR/dY' + dR/dY = CJ*I - dF/dY. */

    pd[3] = *cj - pd[3];
    pd[5] = -pd[5];
    pd[4] = -pd[4];
    pd[6] = *cj - pd[6];

    return 0;
} /* jac2_ */


/* Subroutine */ int rt2_(integer *neq, real_number *t, real_number *y,
	real_number *yp, integer *nrt, real_number *rval, real_number *rpar,
	integer *ipar)
{
    /* Parameter adjustments */
    --rval;
    --yp;
    --y;

    /* Function Body */
    rval[1] = y[1];
    return 0;
} /* rt2_ */

