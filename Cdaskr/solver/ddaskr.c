
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ddaskr_types.h"

/* Table of constant values */

static integer c__49 = 49;
static integer c__201 = 201;
static integer c__0 = 0;
static real_number c_b38 = 0.;
static integer c__47 = 47;
static integer c__202 = 202;
static integer c__1 = 1;
static integer c__41 = 41;
static integer c__203 = 203;
static integer c__4 = 4;
static real_number c_b68 = .6667;
static integer c__2 = 2;
static integer c__56 = 56;
static integer c__501 = 501;
static integer c__502 = 502;
static integer c__503 = 503;
static integer c__3 = 3;
static integer c__38 = 38;
static integer c__610 = 610;
static integer c__48 = 48;
static integer c__611 = 611;
static integer c__620 = 620;
static integer c__621 = 621;
static integer c__43 = 43;
static integer c__622 = 622;
static integer c__630 = 630;
static integer c__28 = 28;
static integer c__631 = 631;
static integer c__44 = 44;
static integer c__640 = 640;
static integer c__57 = 57;
static integer c__641 = 641;
static integer c__650 = 650;
static integer c__651 = 651;
static integer c__40 = 40;
static integer c__652 = 652;
static integer c__655 = 655;
static integer c__46 = 46;
static integer c__656 = 656;
static integer c__660 = 660;
static integer c__661 = 661;
static integer c__670 = 670;
static integer c__45 = 45;
static integer c__671 = 671;
static integer c__672 = 672;
static integer c__675 = 675;
static integer c__51 = 51;
static integer c__676 = 676;
static integer c__677 = 677;
static integer c__680 = 680;
static integer c__36 = 36;
static integer c__681 = 681;
static integer c__685 = 685;
static integer c__686 = 686;
static integer c__690 = 690;
static integer c__35 = 35;
static integer c__691 = 691;
static integer c__695 = 695;
static integer c__50 = 50;
static integer c__696 = 696;
static integer c__25 = 25;
static integer c__34 = 34;
static integer c__60 = 60;
static integer c__5 = 5;
static integer c__39 = 39;
static integer c__6 = 6;
static integer c__7 = 7;
static integer c__8 = 8;
static integer c__54 = 54;
static integer c__9 = 9;
static integer c__10 = 10;
static integer c__11 = 11;
static integer c__29 = 29;
static integer c__12 = 12;
static integer c__13 = 13;
static integer c__14 = 14;
static integer c__15 = 15;
static integer c__52 = 52;
static integer c__17 = 17;
static integer c__18 = 18;
static integer c__19 = 19;
static integer c__20 = 20;
static integer c__21 = 21;
static integer c__22 = 22;
static integer c__58 = 58;
static integer c__23 = 23;
static integer c__24 = 24;
static integer c__26 = 26;
static integer c__27 = 27;
static integer c__30 = 30;
static integer c__31 = 31;
static integer c__701 = 701;
static integer c__702 = 702;
static real_number c_b758 = 1.;
static integer c__901 = 901;
static integer c__902 = 902;
static integer c__903 = 903;
static integer c__904 = 904;
static integer c__905 = 905;
static integer c__42 = 42;
static integer c__906 = 906;
static integer c__921 = 921;
static integer c__922 = 922;
static integer c__923 = 923;
static integer c__924 = 924;
static integer c__925 = 925;
static integer c__926 = 926;

/* Subroutine */ int _daskr_ddaskr_(Unknown_fp res, integer *neq, real_number *t,
	real_number *y, real_number *yprime, real_number *tout, integer *info,
	real_number *rtol, real_number *atol, integer *idid, real_number *rwork,
	integer *lrw, integer *iwork, integer *liw, real_number *rpar, integer
	*ipar, Unknown_fp jac, Unknown_fp psol, Unknown_fp rt, integer *nrt, integer *jroot)
{
    /* System generated locals */
    integer i__1, i__2;
    real_number d__1, d__2;

    /* Builtin functions */
    /* Subroutine */ int _daskr_str_copy(char *, char *, integer, integer);
    double _daskr_real_pow(real_number *, real_number *), _daskr_real_sign(
	    real_number *, real_number *);

    /* Local variables */
    static real_number h__;
    static integer i__;
    static real_number r__, h0;
    static integer le;
    static real_number rh, tn;
    static integer lr0, lr1, ici, idi, lid, ier;
    static char msg[80];
    static integer lwm, irt, lvt, lwt, lrx, nwt, nli0, nni0;
    static integer lcfl, lcfn, done;
    static real_number rcfl;
    static integer nnid;
    static integer lavl;
    static integer maxl, iret;
    static real_number hmax;
    static integer lphi;
    static real_number hmin;
    static integer lyic, lpwk, nstd;
    static real_number rcfn;
    static integer ncfl0, ncfn0;
    extern /* Subroutine */ int _daskr_dnedd_(real_number *x, real_number *y, real_number *yprime,
	integer *neq, Unknown_fp res, Unknown_fp jacd, real_number *pdum, real_number *h__,
	real_number *wt, integer *jstart, integer *idid, real_number *rpar,
	integer *ipar, real_number *phi, real_number *gamma, real_number *dumsvr,
	 real_number *delta, real_number *e, real_number *wm, integer *iwm,
	real_number *cj, real_number *cjold, real_number *cjlast, real_number *s,
	real_number *uround, real_number *dume, real_number *dums, real_number *
	dumr, real_number *epcon, integer *jcalc, integer *jfdum, integer *kp1,
	 integer *nonneg, integer *ntype, integer *iernls);
    static integer mband;
    extern /* Subroutine */ int _daskr_dnedk_(real_number *x, real_number *y, real_number *yprime,
	integer *neq, Unknown_fp res, Unknown_fp jack, Unknown_fp psol, real_number *h__,
	real_number *wt, integer *jstart, integer *idid, real_number *rpar,
	integer *ipar, real_number *phi, real_number *gamma, real_number *savr,
	real_number *delta, real_number *e, real_number *wm, integer *iwm,
	real_number *cj, real_number *cjold, real_number *cjlast, real_number *s,
	real_number *uround, real_number *epli, real_number *sqrtn, real_number *
	rsqrtn, real_number *epcon, integer *jcalc, integer *jflg, integer *
	kp1, integer *nonneg, integer *ntype, integer *iernls);
    static integer lenic, lenid, ncphi, lenpd, lsoff, msave, index, itemp, 
	    leniw, nzflg;
    static real_number atoli;
    static integer lypic;
    static integer lwarn;
    static integer lenwp, lenrw, mxord, nwarn;
    static real_number rtoli;
    static integer lsavr;
    extern real_number _daskr_d1mach_(integer *);
    static real_number tdist, tnext, avlin, fmaxl, tstop;
    extern /* Subroutine */ int _daskr_ddstp_(real_number *, real_number *, real_number
	    *, integer *, Unknown_fp, Unknown_fp, Unknown_fp, real_number *, real_number *,
	    real_number *, integer *, integer *, real_number *, integer *,
	    real_number *, real_number *, real_number *, real_number *,
	    real_number *, integer *, real_number *, real_number *, real_number *,
	     real_number *, real_number *, real_number *, real_number *,
	    real_number *, real_number *, real_number *, real_number *,
	    real_number *, real_number *, real_number *, real_number *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, Unknown_fp), _daskr_dcnst0_(integer *, real_number *, integer *,
	    integer *), _daskr_ddasic_(real_number *, real_number *, real_number *,
	    integer *, integer *, integer *, Unknown_fp, Unknown_fp, Unknown_fp, real_number *,
	    real_number *, real_number *, integer *, integer *, real_number *,
	    integer *, real_number *, real_number *, real_number *, real_number *,
	     real_number *, real_number *, real_number *, real_number *, integer *
	    , real_number *, real_number *, real_number *, real_number *,
	    real_number *, real_number *, integer *, integer *, integer *, Unknown_fp)
	    ;
    extern /* Subroutine */ int _daskr_ddasid_(real_number *x, real_number *y, real_number *yprime,
	 integer *neq, integer *icopt, integer *id, Unknown_fp res, Unknown_fp jacd,
	real_number *pdum, real_number *h__, real_number *tscale, real_number *wt,
	 integer *jsdum, real_number *rpar, integer *ipar, real_number *dumsvr,
	real_number *delta, real_number *r__, real_number *yic, real_number *ypic,
	 real_number *dumpwk, real_number *wm, integer *iwm, real_number *cj,
	real_number *uround, real_number *dume, real_number *dums, real_number *
	dumr, real_number *epcon, real_number *ratemx, real_number *stptol,
	integer *jfdum, integer *icnflg, integer *icnstr, integer *iernls);
    extern /* Subroutine */ int _daskr_ddasik_(real_number *x, real_number *y, real_number *yprime,
	 integer *neq, integer *icopt, integer *id, Unknown_fp res, Unknown_fp jack, Unknown_fp
	psol, real_number *h__, real_number *tscale, real_number *wt, integer *
	jskip, real_number *rpar, integer *ipar, real_number *savr, real_number *
	delta, real_number *r__, real_number *yic, real_number *ypic, real_number
	*pwk, real_number *wm, integer *iwm, real_number *cj, real_number *
	uround, real_number *epli, real_number *sqrtn, real_number *rsqrtn,
	real_number *epcon, real_number *ratemx, real_number *stptol, integer *
	jflg, integer *icnflg, integer *icnstr, integer *iernls);
    extern /* Subroutine */ int _daskr_drchek_(integer *, Unknown_fp, integer *, integer *,
	     real_number *, real_number *, real_number *, real_number *,
	    real_number *, real_number *, integer *, real_number *, real_number *,
	     real_number *, integer *, integer *, real_number *, integer *,
	    real_number *, integer *, real_number *, integer *);
    static integer icnflg;
    static real_number tscale, epconi;
    extern /* Subroutine */ int _daskr_ddatrp_(real_number *, real_number *,
	    real_number *, real_number *, integer *, integer *, real_number *,
	    real_number *);
    static real_number floatn;
    static integer nonneg;
    extern /* Subroutine */ int _daskr_ddawts_(integer *, integer *, real_number *,
	    real_number *, real_number *, real_number *, real_number *, integer *)
	    ;
    extern real_number _daskr_ddwnrm_(integer *, real_number *, real_number *,
	    real_number *, integer *);
    static integer leniwp;
    extern /* Subroutine */ int _daskr_xerrwd_(char *, integer *, integer *, integer
	    *, integer *, integer *, integer *, integer *, real_number *,
	    real_number *, integer), _daskr_dinvwt_(integer *, real_number *, integer *)
	    ;
    static real_number uround, ypnorm;


/* ***BEGIN PROLOGUE  DDASKR */
/* ***REVISION HISTORY  (YYMMDD) */
/*   020815  DATE WRITTEN */
/*   021105  Changed yprime argument in DRCHEK calls to YPRIME. */
/*   021217  Modified error return for zeros found too close together. */
/*   021217  Added root direction output in JROOT. */
/*   040518  Changed adjustment to X2 in Subr. DROOTS. */
/*   050511  Revised stopping tests in statements 530 - 580; reordered */
/*           to test for tn at tstop before testing for tn past tout. */
/*   060712  In DMATD, changed minimum D.Q. increment to 1/EWT(j). */
/*   071003  In DRCHEK, fixed bug in TEMP2 (HMINR) below 110. */
/*   110608  In DRCHEK, fixed bug in setting of T1 at 300. */
/* ***CATEGORY NO.  I1A2 */
/* ***KEYWORDS  DIFFERENTIAL/ALGEBRAIC, BACKWARD DIFFERENTIATION FORMULAS, */
/*             IMPLICIT DIFFERENTIAL SYSTEMS, KRYLOV ITERATION */
/* ***AUTHORS   Linda R. Petzold, Peter N. Brown, Alan C. Hindmarsh, and */
/*                  Clement W. Ulrich */
/*             Center for Computational Sciences & Engineering, L-316 */
/*             Lawrence Livermore National Laboratory */
/*             P.O. Box 808, */
/*             Livermore, CA 94551 */
/* ***PURPOSE  This code solves a system of differential/algebraic */
/*            equations of the form */
/*               G(t,y,y') = 0 , */
/*            using a combination of Backward Differentiation Formula */
/*            (BDF) methods and a choice of two linear system solution */
/*            methods: direct (dense or band) or Krylov (iterative). */
/*            This version is in double precision. */
/* ----------------------------------------------------------------------- */
/* ***DESCRIPTION */

/* *Usage: */

/*      IMPLICIT DOUBLE PRECISION(A-H,O-Z) */
/*      INTEGER NEQ, INFO(N), IDID, LRW, LIW, IWORK(LIW), IPAR(*) */
/*      DOUBLE PRECISION T, Y(*), YPRIME(*), TOUT, RTOL(*), ATOL(*), */
/*         RWORK(LRW), RPAR(*) */
/*      EXTERNAL RES, JAC, PSOL, RT */

/*      CALL DDASKR (RES, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL, */
/*     *             IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC, PSOL, */
/*     *             RT, NRT, JROOT) */

/*  Quantities which may be altered by the code are: */
/*     T, Y(*), YPRIME(*), INFO(1), RTOL, ATOL, IDID, RWORK(*), IWORK(*) */


/* *Arguments: */

/*  RES:EXT          This is the name of a subroutine which you */
/*                   provide to define the residual function G(t,y,y') */
/*                   of the differential/algebraic system. */

/*  NEQ:IN           This is the number of equations in the system. */

/*  T:INOUT          This is the current value of the independent */
/*                   variable. */

/*  Y(*):INOUT       This array contains the solution components at T. */

/*  YPRIME(*):INOUT  This array contains the derivatives of the solution */
/*                   components at T. */

/*  TOUT:IN          This is a point at which a solution is desired. */

/*  INFO(N):IN       This is an integer array used to communicate details */
/*                   of how the solution is to be carried out, such as */
/*                   tolerance type, matrix structure, step size and */
/*                   order limits, and choice of nonlinear system method. */
/*                   N must be at least 20. */

/*  RTOL,ATOL:INOUT  These quantities represent absolute and relative */
/*                   error tolerances (on local error) which you provide */
/*                   to indicate how accurately you wish the solution to */
/*                   be computed.  You may choose them to be both scalars */
/*                   or else both arrays of length NEQ. */

/*  IDID:OUT         This integer scalar is an indicator reporting what */
/*                   the code did.  You must monitor this variable to */
/*                   decide what action to take next. */

/*  RWORK:WORK       A real work array of length LRW which provides the */
/*                   code with needed storage space. */

/*  LRW:IN           The length of RWORK. */

/*  IWORK:WORK       An integer work array of length LIW which provides */
/*                   the code with needed storage space. */

/*  LIW:IN           The length of IWORK. */

/*  RPAR,IPAR:IN     These are real and integer parameter arrays which */
/*                   you can use for communication between your calling */
/*                   program and the RES, JAC, and PSOL subroutines. */

/*  JAC:EXT          This is the name of a subroutine which you may */
/*                   provide (optionally) for calculating Jacobian */
/*                   (partial derivative) data involved in solving linear */
/*                   systems within DDASKR. */

/*  PSOL:EXT         This is the name of a subroutine which you must */
/*                   provide for solving linear systems if you selected */
/*                   a Krylov method.  The purpose of PSOL is to solve */
/*                   linear systems involving a left preconditioner P. */

/*  RT:EXT           This is the name of the subroutine for defining */
/*                   constraint functions Ri(T,Y,Y')) whose roots are */
/*                   desired during the integration.  This name must be */
/*                   declared external in the calling program. */

/*  NRT:IN           This is the number of constraint functions */
/*                   Ri(T,Y,Y').  If there are no constraints, set */
/*                   NRT = 0, and pass a dummy name for RT. */

/*  JROOT:OUT        This is an integer array of length NRT for output */
/*                   of root information. */

/* *Overview */

/*  The DDASKR solver uses the backward differentiation formulas of */
/*  orders one through five to solve a system of the form G(t,y,y') = 0 */
/*  for y = Y and y' = YPRIME.  Values for Y and YPRIME at the initial */
/*  time must be given as input.  These values should be consistent, */
/*  that is, if T, Y, YPRIME are the given initial values, they should */
/*  satisfy G(T,Y,YPRIME) = 0.  However, if consistent values are not */
/*  known, in many cases you can have DDASKR solve for them -- see */
/*  INFO(11). (This and other options are described in detail below.) */

/*  Normally, DDASKR solves the system from T to TOUT.  It is easy to */
/*  continue the solution to get results at additional TOUT.  This is */
/*  the interval mode of operation.  Intermediate results can also be */
/*  obtained easily by specifying INFO(3). */

/*  On each step taken by DDASKR, a sequence of nonlinear algebraic */
/*  systems arises.  These are solved by one of two types of */
/*  methods: */
/*    * a Newton iteration with a direct method for the linear */
/*      systems involved (INFO(12) = 0), or */
/*    * a Newton iteration with a preconditioned Krylov iterative */
/*      method for the linear systems involved (INFO(12) = 1). */

/*  The direct method choices are dense and band matrix solvers, */
/*  with either a user-supplied or an internal difference quotient */
/*  Jacobian matrix, as specified by INFO(5) and INFO(6). */
/*  In the band case, INFO(6) = 1, you must supply half-bandwidths */
/*  in IWORK(1) and IWORK(2). */

/*  The Krylov method is the Generalized Minimum Residual (GMRES) */
/*  method, in either complete or incomplete form, and with */
/*  scaling and preconditioning.  The method is implemented */
/*  in an algorithm called SPIGMR.  Certain options in the Krylov */
/*  method case are specified by INFO(13) and INFO(15). */

/*  If the Krylov method is chosen, you may supply a pair of routines, */
/*  JAC and PSOL, to apply preconditioning to the linear system. */
/*  If the system is A*x = b, the matrix is A = dG/dY + CJ*dG/dYPRIME */
/*  (of order NEQ).  This system can then be preconditioned in the form */
/*  (P-inverse)*A*x = (P-inverse)*b, with left preconditioner P. */
/*  (DDASKR does not allow right preconditioning.) */
/*  Then the Krylov method is applied to this altered, but equivalent, */
/*  linear system, hopefully with much better performance than without */
/*  preconditioning.  (In addition, a diagonal scaling matrix based on */
/*  the tolerances is also introduced into the altered system.) */

/*  The JAC routine evaluates any data needed for solving systems */
/*  with coefficient matrix P, and PSOL carries out that solution. */
/*  In any case, in order to improve convergence, you should try to */
/*  make P approximate the matrix A as much as possible, while keeping */
/*  the system P*x = b reasonably easy and inexpensive to solve for x, */
/*  given a vector b. */

/*  While integrating the given DAE system, DDASKR also searches for */
/*  roots of the given constraint functions Ri(T,Y,Y') given by RT. */
/*  If DDASKR detects a sign change in any Ri(T,Y,Y'), it will return */
/*  the intermediate value of T and Y for which Ri(T,Y,Y') = 0. */
/*  Caution: If some Ri has a root at or very near the initial time, */
/*  DDASKR may fail to find it, or may find extraneous roots there, */
/*  because it does not yet have a sufficient history of the solution. */

/* *Description */

/* ------INPUT - WHAT TO DO ON THE FIRST CALL TO DDASKR------------------- */


/*  The first call of the code is defined to be the start of each new */
/*  problem.  Read through the descriptions of all the following items, */
/*  provide sufficient storage space for designated arrays, set */
/*  appropriate variables for the initialization of the problem, and */
/*  give information about how you want the problem to be solved. */


/*  RES -- Provide a subroutine of the form */

/*             SUBROUTINE RES (T, Y, YPRIME, CJ, DELTA, IRES, RPAR, IPAR) */

/*         to define the system of differential/algebraic */
/*         equations which is to be solved. For the given values */
/*         of T, Y and YPRIME, the subroutine should return */
/*         the residual of the differential/algebraic system */
/*             DELTA = G(T,Y,YPRIME) */
/*         DELTA is a vector of length NEQ which is output from RES. */

/*         Subroutine RES must not alter T, Y, YPRIME, or CJ. */
/*         You must declare the name RES in an EXTERNAL */
/*         statement in your program that calls DDASKR. */
/*         You must dimension Y, YPRIME, and DELTA in RES. */

/*         The input argument CJ can be ignored, or used to rescale */
/*         constraint equations in the system (see Ref. 2, p. 145). */
/*         Note: In this respect, DDASKR is not downward-compatible */
/*         with DDASSL, which does not have the RES argument CJ. */

/*         IRES is an integer flag which is always equal to zero */
/*         on input.  Subroutine RES should alter IRES only if it */
/*         encounters an illegal value of Y or a stop condition. */
/*         Set IRES = -1 if an input value is illegal, and DDASKR */
/*         will try to solve the problem without getting IRES = -1. */
/*         If IRES = -2, DDASKR will return control to the calling */
/*         program with IDID = -11. */

/*         RPAR and IPAR are real and integer parameter arrays which */
/*         you can use for communication between your calling program */
/*         and subroutine RES. They are not altered by DDASKR. If you */
/*         do not need RPAR or IPAR, ignore these parameters by treat- */
/*         ing them as dummy arguments. If you do choose to use them, */
/*         dimension them in your calling program and in RES as arrays */
/*         of appropriate length. */

/*  NEQ -- Set it to the number of equations in the system (NEQ .GE. 1). */

/*  T -- Set it to the initial point of the integration. (T must be */
/*       a variable.) */

/*  Y(*) -- Set this array to the initial values of the NEQ solution */
/*          components at the initial point.  You must dimension Y of */
/*          length at least NEQ in your calling program. */

/*  YPRIME(*) -- Set this array to the initial values of the NEQ first */
/*               derivatives of the solution components at the initial */
/*               point.  You must dimension YPRIME at least NEQ in your */
/*               calling program. */

/*  TOUT - Set it to the first point at which a solution is desired. */
/*         You cannot take TOUT = T.  Integration either forward in T */
/*         (TOUT .GT. T) or backward in T (TOUT .LT. T) is permitted. */

/*         The code advances the solution from T to TOUT using step */
/*         sizes which are automatically selected so as to achieve the */
/*         desired accuracy.  If you wish, the code will return with the */
/*         solution and its derivative at intermediate steps (the */
/*         intermediate-output mode) so that you can monitor them, */
/*         but you still must provide TOUT in accord with the basic */
/*         aim of the code. */

/*         The first step taken by the code is a critical one because */
/*         it must reflect how fast the solution changes near the */
/*         initial point.  The code automatically selects an initial */
/*         step size which is practically always suitable for the */
/*         problem.  By using the fact that the code will not step past */
/*         TOUT in the first step, you could, if necessary, restrict the */
/*         length of the initial step. */

/*         For some problems it may not be permissible to integrate */
/*         past a point TSTOP, because a discontinuity occurs there */
/*         or the solution or its derivative is not defined beyond */
/*         TSTOP.  When you have declared a TSTOP point (see INFO(4) */
/*         and RWORK(1)), you have told the code not to integrate past */
/*         TSTOP.  In this case any tout beyond TSTOP is invalid input. */

/*  INFO(*) - Use the INFO array to give the code more details about */
/*            how you want your problem solved.  This array should be */
/*            dimensioned of length 20, though DDASKR uses only the */
/*            first 15 entries.  You must respond to all of the following */
/*            items, which are arranged as questions.  The simplest use */
/*            of DDASKR corresponds to setting all entries of INFO to 0. */

/*       INFO(1) - This parameter enables the code to initialize itself. */
/*              You must set it to indicate the start of every new */
/*              problem. */

/*          **** Is this the first call for this problem ... */
/*                yes - set INFO(1) = 0 */
/*                 no - not applicable here. */
/*                      See below for continuation calls.  **** */

/*       INFO(2) - How much accuracy you want of your solution */
/*              is specified by the error tolerances RTOL and ATOL. */
/*              The simplest use is to take them both to be scalars. */
/*              To obtain more flexibility, they can both be arrays. */
/*              The code must be told your choice. */

/*          **** Are both error tolerances RTOL, ATOL scalars ... */
/*                yes - set INFO(2) = 0 */
/*                      and input scalars for both RTOL and ATOL */
/*                 no - set INFO(2) = 1 */
/*                      and input arrays for both RTOL and ATOL **** */

/*       INFO(3) - The code integrates from T in the direction of TOUT */
/*              by steps.  If you wish, it will return the computed */
/*              solution and derivative at the next intermediate step */
/*              (the intermediate-output mode) or TOUT, whichever comes */
/*              first.  This is a good way to proceed if you want to */
/*              see the behavior of the solution.  If you must have */
/*              solutions at a great many specific TOUT points, this */
/*              code will compute them efficiently. */

/*          **** Do you want the solution only at */
/*               TOUT (and not at the next intermediate step) ... */
/*                yes - set INFO(3) = 0 (interval-output mode) */
/*                 no - set INFO(3) = 1 (intermediate-output mode) **** */

/*       INFO(4) - To handle solutions at a great many specific */
/*              values TOUT efficiently, this code may integrate past */
/*              TOUT and interpolate to obtain the result at TOUT. */
/*              Sometimes it is not possible to integrate beyond some */
/*              point TSTOP because the equation changes there or it is */
/*              not defined past TSTOP.  Then you must tell the code */
/*              this stop condition. */

/*           **** Can the integration be carried out without any */
/*                restrictions on the independent variable T ... */
/*                 yes - set INFO(4) = 0 */
/*                  no - set INFO(4) = 1 */
/*                       and define the stopping point TSTOP by */
/*                       setting RWORK(1) = TSTOP **** */

/*       INFO(5) - used only when INFO(12) = 0 (direct methods). */
/*              To solve differential/algebraic systems you may wish */
/*              to use a matrix of partial derivatives of the */
/*              system of differential equations.  If you do not */
/*              provide a subroutine to evaluate it analytically (see */
/*              description of the item JAC in the call list), it will */
/*              be approximated by numerical differencing in this code. */
/*              Although it is less trouble for you to have the code */
/*              compute partial derivatives by numerical differencing, */
/*              the solution will be more reliable if you provide the */
/*              derivatives via JAC.  Usually numerical differencing is */
/*              more costly than evaluating derivatives in JAC, but */
/*              sometimes it is not - this depends on your problem. */

/*           **** Do you want the code to evaluate the partial deriv- */
/*                atives automatically by numerical differences ... */
/*                 yes - set INFO(5) = 0 */
/*                  no - set INFO(5) = 1 */
/*                       and provide subroutine JAC for evaluating the */
/*                       matrix of partial derivatives **** */

/*       INFO(6) - used only when INFO(12) = 0 (direct methods). */
/*              DDASKR will perform much better if the matrix of */
/*              partial derivatives, dG/dY + CJ*dG/dYPRIME (here CJ is */
/*              a scalar determined by DDASKR), is banded and the code */
/*              is told this.  In this case, the storage needed will be */
/*              greatly reduced, numerical differencing will be performed */
/*              much cheaper, and a number of important algorithms will */
/*              execute much faster.  The differential equation is said */
/*              to have half-bandwidths ML (lower) and MU (upper) if */
/*              equation i involves only unknowns Y(j) with */
/*                             i-ML .le. j .le. i+MU . */
/*              For all i=1,2,...,NEQ.  Thus, ML and MU are the widths */
/*              of the lower and upper parts of the band, respectively, */
/*              with the main diagonal being excluded.  If you do not */
/*              indicate that the equation has a banded matrix of partial */
/*              derivatives the code works with a full matrix of NEQ**2 */
/*              elements (stored in the conventional way).  Computations */
/*              with banded matrices cost less time and storage than with */
/*              full matrices if  2*ML+MU .lt. NEQ.  If you tell the */
/*              code that the matrix of partial derivatives has a banded */
/*              structure and you want to provide subroutine JAC to */
/*              compute the partial derivatives, then you must be careful */
/*              to store the elements of the matrix in the special form */
/*              indicated in the description of JAC. */

/*          **** Do you want to solve the problem using a full (dense) */
/*               matrix (and not a special banded structure) ... */
/*                yes - set INFO(6) = 0 */
/*                 no - set INFO(6) = 1 */
/*                       and provide the lower (ML) and upper (MU) */
/*                       bandwidths by setting */
/*                       IWORK(1)=ML */
/*                       IWORK(2)=MU **** */

/*       INFO(7) - You can specify a maximum (absolute value of) */
/*              stepsize, so that the code will avoid passing over very */
/*              large regions. */

/*          ****  Do you want the code to decide on its own the maximum */
/*                stepsize ... */
/*                 yes - set INFO(7) = 0 */
/*                  no - set INFO(7) = 1 */
/*                       and define HMAX by setting */
/*                       RWORK(2) = HMAX **** */

/*       INFO(8) -  Differential/algebraic problems may occasionally */
/*              suffer from severe scaling difficulties on the first */
/*              step.  If you know a great deal about the scaling of */
/*              your problem, you can help to alleviate this problem */
/*              by specifying an initial stepsize H0. */

/*          ****  Do you want the code to define its own initial */
/*                stepsize ... */
/*                 yes - set INFO(8) = 0 */
/*                  no - set INFO(8) = 1 */
/*                       and define H0 by setting */
/*                       RWORK(3) = H0 **** */

/*       INFO(9) -  If storage is a severe problem, you can save some */
/*              storage by restricting the maximum method order MAXORD. */
/*              The default value is 5.  For each order decrease below 5, */
/*              the code requires NEQ fewer locations, but it is likely */
/*              to be slower.  In any case, you must have */
/*              1 .le. MAXORD .le. 5. */
/*          ****  Do you want the maximum order to default to 5 ... */
/*                 yes - set INFO(9) = 0 */
/*                  no - set INFO(9) = 1 */
/*                       and define MAXORD by setting */
/*                       IWORK(3) = MAXORD **** */

/*       INFO(10) - If you know that certain components of the */
/*              solutions to your equations are always nonnegative */
/*              (or nonpositive), it may help to set this */
/*              parameter.  There are three options that are */
/*              available: */
/*              1.  To have constraint checking only in the initial */
/*                  condition calculation. */
/*              2.  To enforce nonnegativity in Y during the integration. */
/*              3.  To enforce both options 1 and 2. */

/*              When selecting option 2 or 3, it is probably best to try */
/*              the code without using this option first, and only use */
/*              this option if that does not work very well. */

/*          ****  Do you want the code to solve the problem without */
/*                invoking any special inequality constraints ... */
/*                 yes - set INFO(10) = 0 */
/*                  no - set INFO(10) = 1 to have option 1 enforced */
/*                  no - set INFO(10) = 2 to have option 2 enforced */
/*                  no - set INFO(10) = 3 to have option 3 enforced **** */

/*                  If you have specified INFO(10) = 1 or 3, then you */
/*                  will also need to identify how each component of Y */
/*                  in the initial condition calculation is constrained. */
/*                  You must set: */
/*                  IWORK(40+I) = +1 if Y(I) must be .GE. 0, */
/*                  IWORK(40+I) = +2 if Y(I) must be .GT. 0, */
/*                  IWORK(40+I) = -1 if Y(I) must be .LE. 0, while */
/*                  IWORK(40+I) = -2 if Y(I) must be .LT. 0, while */
/*                  IWORK(40+I) =  0 if Y(I) is not constrained. */

/*       INFO(11) - DDASKR normally requires the initial T, Y, and */
/*              YPRIME to be consistent.  That is, you must have */
/*              G(T,Y,YPRIME) = 0 at the initial T.  If you do not know */
/*              the initial conditions precisely, in some cases */
/*              DDASKR may be able to compute it. */

/*              Denoting the differential variables in Y by Y_d */
/*              and the algebraic variables by Y_a, DDASKR can solve */
/*              one of two initialization problems: */
/*              1.  Given Y_d, calculate Y_a and Y'_d, or */
/*              2.  Given Y', calculate Y. */
/*              In either case, initial values for the given */
/*              components are input, and initial guesses for */
/*              the unknown components must also be provided as input. */

/*          ****  Are the initial T, Y, YPRIME consistent ... */

/*                 yes - set INFO(11) = 0 */
/*                  no - set INFO(11) = 1 to calculate option 1 above, */
/*                    or set INFO(11) = 2 to calculate option 2 **** */

/*                  If you have specified INFO(11) = 1, then you */
/*                  will also need to identify  which are the */
/*                  differential and which are the algebraic */
/*                  components (algebraic components are components */
/*                  whose derivatives do not appear explicitly */
/*                  in the function G(T,Y,YPRIME)).  You must set: */
/*                  IWORK(LID+I) = +1 if Y(I) is a differential variable */
/*                  IWORK(LID+I) = -1 if Y(I) is an algebraic variable, */
/*                  where LID = 40 if INFO(10) = 0 or 2 and LID = 40+NEQ */
/*                  if INFO(10) = 1 or 3. */

/*       INFO(12) - Except for the addition of the RES argument CJ, */
/*              DDASKR by default is downward-compatible with DDASSL, */
/*              which uses only direct (dense or band) methods to solve */
/*              the linear systems involved.  You must set INFO(12) to */
/*              indicate whether you want the direct methods or the */
/*              Krylov iterative method. */
/*          ****   Do you want DDASKR to use standard direct methods */
/*                 (dense or band) or the Krylov (iterative) method ... */
/*                   direct methods - set INFO(12) = 0. */
/*                   Krylov method  - set INFO(12) = 1, */
/*                       and check the settings of INFO(13) and INFO(15). */

/*       INFO(13) - used when INFO(12) = 1 (Krylov methods). */
/*              DDASKR uses scalars MAXL, KMP, NRMAX, and EPLI for the */
/*              iterative solution of linear systems.  INFO(13) allows */
/*              you to override the default values of these parameters. */
/*              These parameters and their defaults are as follows: */
/*              MAXL = maximum number of iterations in the SPIGMR */
/*                 algorithm (MAXL .le. NEQ).  The default is */
/*                 MAXL = MIN(5,NEQ). */
/*              KMP = number of vectors on which orthogonalization is */
/*                 done in the SPIGMR algorithm.  The default is */
/*                 KMP = MAXL, which corresponds to complete GMRES */
/*                 iteration, as opposed to the incomplete form. */
/*              NRMAX = maximum number of restarts of the SPIGMR */
/*                 algorithm per nonlinear iteration.  The default is */
/*                 NRMAX = 5. */
/*              EPLI = convergence test constant in SPIGMR algorithm. */
/*                 The default is EPLI = 0.05. */
/*              Note that the length of RWORK depends on both MAXL */
/*              and KMP.  See the definition of LRW below. */
/*          ****   Are MAXL, KMP, and EPLI to be given their */
/*                 default values ... */
/*                  yes - set INFO(13) = 0 */
/*                   no - set INFO(13) = 1, */
/*                        and set all of the following: */
/*                        IWORK(24) = MAXL (1 .le. MAXL .le. NEQ) */
/*                        IWORK(25) = KMP  (1 .le. KMP .le. MAXL) */
/*                        IWORK(26) = NRMAX  (NRMAX .ge. 0) */
/*                        RWORK(10) = EPLI (0 .lt. EPLI .lt. 1.0) **** */

/*        INFO(14) - used with INFO(11) > 0 (initial condition */
/*               calculation is requested).  In this case, you may */
/*               request control to be returned to the calling program */
/*               immediately after the initial condition calculation, */
/*               before proceeding to the integration of the system */
/*               (e.g. to examine the computed Y and YPRIME). */
/*               If this is done, and if the initialization succeeded */
/*               (IDID = 4), you should reset INFO(11) to 0 for the */
/*               next call, to prevent the solver from repeating the */
/*               initialization (and to avoid an infinite loop). */
/*          ****   Do you want to proceed to the integration after */
/*                 the initial condition calculation is done ... */
/*                 yes - set INFO(14) = 0 */
/*                  no - set INFO(14) = 1                        **** */

/*        INFO(15) - used when INFO(12) = 1 (Krylov methods). */
/*               When using preconditioning in the Krylov method, */
/*               you must supply a subroutine, PSOL, which solves the */
/*               associated linear systems using P. */
/*               The usage of DDASKR is simpler if PSOL can carry out */
/*               the solution without any prior calculation of data. */
/*               However, if some partial derivative data is to be */
/*               calculated in advance and used repeatedly in PSOL, */
/*               then you must supply a JAC routine to do this, */
/*               and set INFO(15) to indicate that JAC is to be called */
/*               for this purpose.  For example, P might be an */
/*               approximation to a part of the matrix A which can be */
/*               calculated and LU-factored for repeated solutions of */
/*               the preconditioner system.  The arrays WP and IWP */
/*               (described under JAC and PSOL) can be used to */
/*               communicate data between JAC and PSOL. */
/*          ****   Does PSOL operate with no prior preparation ... */
/*                 yes - set INFO(15) = 0 (no JAC routine) */
/*                  no - set INFO(15) = 1 */
/*                       and supply a JAC routine to evaluate and */
/*                       preprocess any required Jacobian data.  **** */

/*         INFO(16) - option to exclude algebraic variables from */
/*               the error test. */
/*          ****   Do you wish to control errors locally on */
/*                 all the variables... */
/*                 yes - set INFO(16) = 0 */
/*                  no - set INFO(16) = 1 */
/*                       If you have specified INFO(16) = 1, then you */
/*                       will also need to identify  which are the */
/*                       differential and which are the algebraic */
/*                       components (algebraic components are components */
/*                       whose derivatives do not appear explicitly */
/*                       in the function G(T,Y,YPRIME)).  You must set: */
/*                       IWORK(LID+I) = +1 if Y(I) is a differential */
/*                                      variable, and */
/*                       IWORK(LID+I) = -1 if Y(I) is an algebraic */
/*                                      variable, */
/*                       where LID = 40 if INFO(10) = 0 or 2 and */
/*                       LID = 40 + NEQ if INFO(10) = 1 or 3. */

/*       INFO(17) - used when INFO(11) > 0 (DDASKR is to do an */
/*              initial condition calculation). */
/*              DDASKR uses several heuristic control quantities in the */
/*              initial condition calculation.  They have default values, */
/*              but can  also be set by the user using INFO(17). */
/*              These parameters and their defaults are as follows: */
/*              MXNIT  = maximum number of Newton iterations */
/*                 per Jacobian or preconditioner evaluation. */
/*                 The default is: */
/*                 MXNIT =  5 in the direct case (INFO(12) = 0), and */
/*                 MXNIT = 15 in the Krylov case (INFO(12) = 1). */
/*              MXNJ   = maximum number of Jacobian or preconditioner */
/*                 evaluations.  The default is: */
/*                 MXNJ = 6 in the direct case (INFO(12) = 0), and */
/*                 MXNJ = 2 in the Krylov case (INFO(12) = 1). */
/*              MXNH   = maximum number of values of the artificial */
/*                 stepsize parameter H to be tried if INFO(11) = 1. */
/*                 The default is MXNH = 5. */
/*                 NOTE: the maximum number of Newton iterations */
/*                 allowed in all is MXNIT*MXNJ*MXNH if INFO(11) = 1, */
/*                 and MXNIT*MXNJ if INFO(11) = 2. */
/*              LSOFF  = flag to turn off the linesearch algorithm */
/*                 (LSOFF = 0 means linesearch is on, LSOFF = 1 means */
/*                 it is turned off).  The default is LSOFF = 0. */
/*              STPTOL = minimum scaled step in linesearch algorithm. */
/*                 The default is STPTOL = (unit roundoff)**(2/3). */
/*              EPINIT = swing factor in the Newton iteration convergence */
/*                 test.  The test is applied to the residual vector, */
/*                 premultiplied by the approximate Jacobian (in the */
/*                 direct case) or the preconditioner (in the Krylov */
/*                 case).  For convergence, the weighted RMS norm of */
/*                 this vector (scaled by the error weights) must be */
/*                 less than EPINIT*EPCON, where EPCON = .33 is the */
/*                 analogous test constant used in the time steps. */
/*                 The default is EPINIT = .01. */
/*          ****   Are the initial condition heuristic controls to be */
/*                 given their default values... */
/*                  yes - set INFO(17) = 0 */
/*                   no - set INFO(17) = 1, */
/*                        and set all of the following: */
/*                        IWORK(32) = MXNIT (.GT. 0) */
/*                        IWORK(33) = MXNJ (.GT. 0) */
/*                        IWORK(34) = MXNH (.GT. 0) */
/*                        IWORK(35) = LSOFF ( = 0 or 1) */
/*                        RWORK(14) = STPTOL (.GT. 0.0) */
/*                        RWORK(15) = EPINIT (.GT. 0.0)  **** */

/*         INFO(18) - option to get extra printing in initial condition */
/*                calculation. */
/*          ****   Do you wish to have extra printing... */
/*                 no  - set INFO(18) = 0 */
/*                 yes - set INFO(18) = 1 for minimal printing, or */
/*                       set INFO(18) = 2 for full printing. */
/*                       If you have specified INFO(18) .ge. 1, data */
/*                       will be printed with the error handler routines. */
/*                       To print to a non-default unit number L, include */
/*                       the line  CALL XSETUN(L)  in your program.  **** */

/*   RTOL, ATOL -- You must assign relative (RTOL) and absolute (ATOL) */
/*               error tolerances to tell the code how accurately you */
/*               want the solution to be computed.  They must be defined */
/*               as variables because the code may change them. */
/*               you have two choices -- */
/*                     Both RTOL and ATOL are scalars (INFO(2) = 0), or */
/*                     both RTOL and ATOL are vectors (INFO(2) = 1). */
/*               In either case all components must be non-negative. */

/*               The tolerances are used by the code in a local error */
/*               test at each step which requires roughly that */
/*                        abs(local error in Y(i)) .le. EWT(i) , */
/*               where EWT(i) = RTOL*abs(Y(i)) + ATOL is an error weight */
/*               quantity, for each vector component. */
/*               (More specifically, a root-mean-square norm is used to */
/*               measure the size of vectors, and the error test uses the */
/*               magnitude of the solution at the beginning of the step.) */

/*               The true (global) error is the difference between the */
/*               true solution of the initial value problem and the */
/*               computed approximation.  Practically all present day */
/*               codes, including this one, control the local error at */
/*               each step and do not even attempt to control the global */
/*               error directly. */

/*               Usually, but not always, the true accuracy of */
/*               the computed Y is comparable to the error tolerances. */
/*               This code will usually, but not always, deliver a more */
/*               accurate solution if you reduce the tolerances and */
/*               integrate again.  By comparing two such solutions you */
/*               can get a fairly reliable idea of the true error in the */
/*               solution at the larger tolerances. */

/*               Setting ATOL = 0. results in a pure relative error test */
/*               on that component.  Setting RTOL = 0. results in a pure */
/*               absolute error test on that component.  A mixed test */
/*               with non-zero RTOL and ATOL corresponds roughly to a */
/*               relative error test when the solution component is */
/*               much bigger than ATOL and to an absolute error test */
/*               when the solution component is smaller than the */
/*               threshold ATOL. */

/*               The code will not attempt to compute a solution at an */
/*               accuracy unreasonable for the machine being used.  It */
/*               will advise you if you ask for too much accuracy and */
/*               inform you as to the maximum accuracy it believes */
/*               possible. */

/*  RWORK(*) -- a real work array, which should be dimensioned in your */
/*               calling program with a length equal to the value of */
/*               LRW (or greater). */

/*  LRW -- Set it to the declared length of the RWORK array.  The */
/*               minimum length depends on the options you have selected, */
/*               given by a base value plus additional storage as */
/*               described below. */

/*               If INFO(12) = 0 (standard direct method), the base value */
/*               is BASE = 60 + MAX(MAXORD+4,7)*NEQ + 3*NRT. */
/*               The default value is MAXORD = 5 (see INFO(9)).  With the */
/*               default MAXORD, BASE = 60 + 9*NEQ + 3*NRT. */
/*               Additional storage must be added to the base value for */
/*               any or all of the following options: */
/*                 If INFO(6) = 0 (dense matrix), add NEQ**2. */
/*                 If INFO(6) = 1 (banded matrix), then: */
/*                    if INFO(5) = 0, add (2*ML+MU+1)*NEQ */
/*                                           + 2*[NEQ/(ML+MU+1) + 1], and */
/*                    if INFO(5) = 1, add (2*ML+MU+1)*NEQ. */
/*                 If INFO(16) = 1, add NEQ. */

/*               If INFO(12) = 1 (Krylov method), the base value is */
/*               BASE = 60 + (MAXORD+5)*NEQ + 3*NRT */
/*                         + [MAXL + 3 + MIN(1,MAXL-KMP)]*NEQ */
/*                         + (MAXL+3)*MAXL + 1 + LENWP. */
/*               See PSOL for description of LENWP.  The default values */
/*               are: MAXORD = 5 (see INFO(9)), MAXL = MIN(5,NEQ) and */
/*               KMP = MAXL  (see INFO(13)).  With these default values, */
/*               BASE = 101 + 18*NEQ + 3*NRT + LENWP. */
/*               Additional storage must be added to the base value for */
/*               the following option: */
/*                 If INFO(16) = 1, add NEQ. */


/*  IWORK(*) -- an integer work array, which should be dimensioned in */
/*              your calling program with a length equal to the value */
/*              of LIW (or greater). */

/*  LIW -- Set it to the declared length of the IWORK array.  The */
/*             minimum length depends on the options you have selected, */
/*             given by a base value plus additions as described below. */

/*             If INFO(12) = 0 (standard direct method), the base value */
/*             is BASE = 40 + NEQ. */
/*             IF INFO(10) = 1 or 3, add NEQ to the base value. */
/*             If INFO(11) = 1 or INFO(16) =1, add NEQ to the base value. */

/*             If INFO(12) = 1 (Krylov method), the base value is */
/*             BASE = 40 + LENIWP.  See PSOL for description of LENIWP. */
/*             If INFO(10) = 1 or 3, add NEQ to the base value. */
/*             If INFO(11) = 1 or INFO(16) =1, add NEQ to the base value. */


/*  RPAR, IPAR -- These are arrays of double precision and integer type, */
/*             respectively, which are available for you to use */
/*             for communication between your program that calls */
/*             DDASKR and the RES subroutine (and the JAC and PSOL */
/*             subroutines).  They are not altered by DDASKR. */
/*             If you do not need RPAR or IPAR, ignore these */
/*             parameters by treating them as dummy arguments. */
/*             If you do choose to use them, dimension them in */
/*             your calling program and in RES (and in JAC and PSOL) */
/*             as arrays of appropriate length. */

/*  JAC -- This is the name of a routine that you may supply */
/*         (optionally) that relates to the Jacobian matrix of the */
/*         nonlinear system that the code must solve at each T step. */
/*         The role of JAC (and its call sequence) depends on whether */
/*         a direct (INFO(12) = 0) or Krylov (INFO(12) = 1) method */
/*         is selected. */

/*         **** INFO(12) = 0 (direct methods): */
/*           If you are letting the code generate partial derivatives */
/*           numerically (INFO(5) = 0), then JAC can be absent */
/*           (or perhaps a dummy routine to satisfy the loader). */
/*           Otherwise you must supply a JAC routine to compute */
/*           the matrix A = dG/dY + CJ*dG/dYPRIME.  It must have */
/*           the form */

/*           SUBROUTINE JAC (T, Y, YPRIME, DELTA, PD, CJ, H, WT, RPAR, IPAR) */

/*           The JAC routine must dimension Y, YPRIME, and PD (and RPAR */
/*           and IPAR if used).  CJ is a scalar which is input to JAC. */
/*           For the given values of T, Y, and YPRIME, the JAC routine */
/*           must evaluate the nonzero elements of the matrix A, and */
/*           store these values in the array PD.  The elements of PD are */
/*           set to zero before each call to JAC, so that only nonzero */
/*           elements need to be defined. */
/*			 Also DELTA contains the last evaluation delta of RES. */
/*		     H is a scalar and contains current step size in integration. */
/*		     EWT Vector of error weights for computing norms. */
/*           The way you store the elements into the PD array depends */
/*           on the structure of the matrix indicated by INFO(6). */
/*           *** INFO(6) = 0 (full or dense matrix) *** */
/*               Give PD a first dimension of NEQ.  When you evaluate the */
/*               nonzero partial derivatives of equation i (i.e. of G(i)) */
/*               with respect to component j (of Y and YPRIME), you must */
/*               store the element in PD according to */
/*                  PD(i,j) = dG(i)/dY(j) + CJ*dG(i)/dYPRIME(j). */
/*           *** INFO(6) = 1 (banded matrix with half-bandwidths ML, MU */
/*                            as described under INFO(6)) *** */
/*               Give PD a first dimension of 2*ML+MU+1.  When you */
/*               evaluate the nonzero partial derivatives of equation i */
/*               (i.e. of G(i)) with respect to component j (of Y and */
/*               YPRIME), you must store the element in PD according to */
/*                  IROW = i - j + ML + MU + 1 */
/*                  PD(IROW,j) = dG(i)/dY(j) + CJ*dG(i)/dYPRIME(j). */

/*          **** INFO(12) = 1 (Krylov method): */
/*            If you are not calculating Jacobian data in advance for use */
/*            in PSOL (INFO(15) = 0), JAC can be absent (or perhaps a */
/*            dummy routine to satisfy the loader).  Otherwise, you may */
/*            supply a JAC routine to compute and preprocess any parts of */
/*            of the Jacobian matrix  A = dG/dY + CJ*dG/dYPRIME that are */
/*            involved in the preconditioner matrix P. */
/*            It is to have the form */

/*            SUBROUTINE JAC (RES, IRES, NEQ, T, Y, YPRIME, REWT, SAVR, */
/*                            WK, H, CJ, WP, IWP, IER, RPAR, IPAR) */

/*           The JAC routine must dimension Y, YPRIME, REWT, SAVR, WK, */
/*           and (if used) WP, IWP, RPAR, and IPAR. */
/*           The Y, YPRIME, and SAVR arrays contain the current values */
/*           of Y, YPRIME, and the residual G, respectively. */
/*           The array WK is work space of length NEQ. */
/*           H is the step size.  CJ is a scalar, input to JAC, that is */
/*           normally proportional to 1/H.  REWT is an array of */
/*           reciprocal error weights, 1/EWT(i), where EWT(i) is */
/*           RTOL*abs(Y(i)) + ATOL (unless you supplied routine DDAWTS */
/*           instead), for use in JAC if needed.  For example, if JAC */
/*           computes difference quotient approximations to partial */
/*           derivatives, the REWT array may be useful in setting the */
/*           increments used.  The JAC routine should do any */
/*           factorization operations called for, in preparation for */
/*           solving linear systems in PSOL.  The matrix P should */
/*           be an approximation to the Jacobian, */
/*           A = dG/dY + CJ*dG/dYPRIME. */

/*           WP and IWP are real and integer work arrays which you may */
/*           use for communication between your JAC routine and your */
/*           PSOL routine.  These may be used to store elements of the */
/*           preconditioner P, or related matrix data (such as factored */
/*           forms).  They are not altered by DDASKR. */
/*           If you do not need WP or IWP, ignore these parameters by */
/*           treating them as dummy arguments.  If you do use them, */
/*           dimension them appropriately in your JAC and PSOL routines. */
/*           See the PSOL description for instructions on setting */
/*           the lengths of WP and IWP. */

/*           On return, JAC should set the error flag IER as follows.. */
/*             IER = 0    if JAC was successful, */
/*             IER .ne. 0 if JAC was unsuccessful (e.g. if Y or YPRIME */
/*                        was illegal, or a singular matrix is found). */
/*           (If IER .ne. 0, a smaller stepsize will be tried.) */
/*           IER = 0 on entry to JAC, so need be reset only on a failure. */
/*           If RES is used within JAC, then a nonzero value of IRES will */
/*           override any nonzero value of IER (see the RES description). */

/*         Regardless of the method type, subroutine JAC must not */
/*         alter T, Y(*), YPRIME(*), H, CJ, or REWT(*). */
/*         You must declare the name JAC in an EXTERNAL statement in */
/*         your program that calls DDASKR. */

/* PSOL --  This is the name of a routine you must supply if you have */
/*         selected a Krylov method (INFO(12) = 1) with preconditioning. */
/*         In the direct case (INFO(12) = 0), PSOL can be absent */
/*         (a dummy routine may have to be supplied to satisfy the */
/*         loader).  Otherwise, you must provide a PSOL routine to */
/*         solve linear systems arising from preconditioning. */
/*         When supplied with INFO(12) = 1, the PSOL routine is to */
/*         have the form */

/*         SUBROUTINE PSOL (NEQ, T, Y, YPRIME, SAVR, WK, CJ, WGHT, */
/*                          WP, IWP, B, EPLIN, IER, RPAR, IPAR) */

/*         The PSOL routine must solve linear systems of the form */
/*         P*x = b where P is the left preconditioner matrix. */

/*         The right-hand side vector b is in the B array on input, and */
/*         PSOL must return the solution vector x in B. */
/*         The Y, YPRIME, and SAVR arrays contain the current values */
/*         of Y, YPRIME, and the residual G, respectively. */

/*         Work space required by JAC and/or PSOL, and space for data to */
/*         be communicated from JAC to PSOL is made available in the form */
/*         of arrays WP and IWP, which are parts of the RWORK and IWORK */
/*         arrays, respectively.  The lengths of these real and integer */
/*         work spaces WP and IWP must be supplied in LENWP and LENIWP, */
/*         respectively, as follows.. */
/*           IWORK(27) = LENWP = length of real work space WP */
/*           IWORK(28) = LENIWP = length of integer work space IWP. */

/*         WK is a work array of length NEQ for use by PSOL. */
/*         CJ is a scalar, input to PSOL, that is normally proportional */
/*         to 1/H (H = stepsize).  If the old value of CJ */
/*         (at the time of the last JAC call) is needed, it must have */
/*         been saved by JAC in WP. */

/*         WGHT is an array of weights, to be used if PSOL uses an */
/*         iterative method and performs a convergence test.  (In terms */
/*         of the argument REWT to JAC, WGHT is REWT/sqrt(NEQ).) */
/*         If PSOL uses an iterative method, it should use EPLIN */
/*         (a heuristic parameter) as the bound on the weighted norm of */
/*         the residual for the computed solution.  Specifically, the */
/*         residual vector R should satisfy */
/*              SQRT (SUM ( (R(i)*WGHT(i))**2 ) ) .le. EPLIN */

/*         PSOL must not alter NEQ, T, Y, YPRIME, SAVR, CJ, WGHT, EPLIN. */

/*         On return, PSOL should set the error flag IER as follows.. */
/*           IER = 0 if PSOL was successful, */
/*           IER .lt. 0 if an unrecoverable error occurred, meaning */
/*                 control will be passed to the calling routine, */
/*           IER .gt. 0 if a recoverable error occurred, meaning that */
/*                 the step will be retried with the same step size */
/*                 but with a call to JAC to update necessary data, */
/*                 unless the Jacobian data is current, in which case */
/*                 the step will be retried with a smaller step size. */
/*           IER = 0 on entry to PSOL so need be reset only on a failure. */

/*         You must declare the name PSOL in an EXTERNAL statement in */
/*         your program that calls DDASKR. */

/* RT --   This is the name of the subroutine for defining the vector */
/*         R(T,Y,Y') of constraint functions Ri(T,Y,Y'), whose roots */
/*         are desired during the integration.  It is to have the form */
/*             SUBROUTINE RT(NEQ, T, Y, YP, NRT, RVAL, RPAR, IPAR) */
/*             DIMENSION Y(NEQ), YP(NEQ), RVAL(NRT), */
/*         where NEQ, T, Y and NRT are INPUT, and the array RVAL is */
/*         output.  NEQ, T, Y, and YP have the same meaning as in the */
/*         RES routine, and RVAL is an array of length NRT. */
/*         For i = 1,...,NRT, this routine is to load into RVAL(i) the */
/*         value at (T,Y,Y') of the i-th constraint function Ri(T,Y,Y'). */
/*         DDASKR will find roots of the Ri of odd multiplicity */
/*         (that is, sign changes) as they occur during the integration. */
/*         RT must be declared EXTERNAL in the calling program. */

/*         CAUTION.. Because of numerical errors in the functions Ri */
/*         due to roundoff and integration error, DDASKR may return */
/*         false roots, or return the same root at two or more nearly */
/*         equal values of T.  If such false roots are suspected, */
/*         the user should consider smaller error tolerances and/or */
/*         higher precision in the evaluation of the Ri. */

/*         If a root of some Ri defines the end of the problem, */
/*         the input to DDASKR should nevertheless allow */
/*         integration to a point slightly past that root, so */
/*         that DDASKR can locate the root by interpolation. */

/* NRT --  The number of constraint functions Ri(T,Y,Y').  If there are */
/*         no constraints, set NRT = 0 and pass a dummy name for RT. */

/* JROOT -- This is an integer array of length NRT, used only for output. */
/*         On a return where one or more roots were found (IDID = 5), */
/*         JROOT(i) = 1 or -1 if Ri(T,Y,Y') has a root at T, and */
/*         JROOT(i) = 0 if not.  If nonzero, JROOT(i) shows the direction */
/*         of the sign change in Ri in the direction of integration: */
/*         JROOT(i) = 1  means Ri changed from negative to positive. */
/*         JROOT(i) = -1 means Ri changed from positive to negative. */


/*  OPTIONALLY REPLACEABLE SUBROUTINE: */

/*  DDASKR uses a weighted root-mean-square norm to measure the */
/*  size of various error vectors.  The weights used in this norm */
/*  are set in the following subroutine: */

/*    SUBROUTINE DDAWTS (NEQ, IWT, RTOL, ATOL, Y, EWT, RPAR, IPAR) */
/*    DIMENSION RTOL(*), ATOL(*), Y(*), EWT(*), RPAR(*), IPAR(*) */

/*  A DDAWTS routine has been included with DDASKR which sets the */
/*  weights according to */
/*    EWT(I) = RTOL*ABS(Y(I)) + ATOL */
/*  in the case of scalar tolerances (IWT = 0) or */
/*    EWT(I) = RTOL(I)*ABS(Y(I)) + ATOL(I) */
/*  in the case of array tolerances (IWT = 1).  (IWT is INFO(2).) */
/*  In some special cases, it may be appropriate for you to define */
/*  your own error weights by writing a subroutine DDAWTS to be */
/*  called instead of the version supplied.  However, this should */
/*  be attempted only after careful thought and consideration. */
/*  If you supply this routine, you may use the tolerances and Y */
/*  as appropriate, but do not overwrite these variables.  You */
/*  may also use RPAR and IPAR to communicate data as appropriate. */
/*  ***Note: Aside from the values of the weights, the choice of */
/*  norm used in DDASKR (weighted root-mean-square) is not subject */
/*  to replacement by the user.  In this respect, DDASKR is not */
/*  downward-compatible with the original DDASSL solver (in which */
/*  the norm routine was optionally user-replaceable). */


/* ------OUTPUT - AFTER ANY RETURN FROM DDASKR---------------------------- */

/*  The principal aim of the code is to return a computed solution at */
/*  T = TOUT, although it is also possible to obtain intermediate */
/*  results along the way.  To find out whether the code achieved its */
/*  goal or if the integration process was interrupted before the task */
/*  was completed, you must check the IDID parameter. */


/*   T -- The output value of T is the point to which the solution */
/*        was successfully advanced. */

/*   Y(*) -- contains the computed solution approximation at T. */

/*   YPRIME(*) -- contains the computed derivative approximation at T. */

/*   IDID -- reports what the code did, described as follows: */

/*                     *** TASK COMPLETED *** */
/*                Reported by positive values of IDID */

/*           IDID = 1 -- A step was successfully taken in the */
/*                   interval-output mode.  The code has not */
/*                   yet reached TOUT. */

/*           IDID = 2 -- The integration to TSTOP was successfully */
/*                   completed (T = TSTOP) by stepping exactly to TSTOP. */

/*           IDID = 3 -- The integration to TOUT was successfully */
/*                   completed (T = TOUT) by stepping past TOUT. */
/*                   Y(*) and YPRIME(*) are obtained by interpolation. */

/*           IDID = 4 -- The initial condition calculation, with */
/*                   INFO(11) > 0, was successful, and INFO(14) = 1. */
/*                   No integration steps were taken, and the solution */
/*                   is not considered to have been started. */

/*           IDID = 5 -- The integration was successfully completed */
/*                   by finding one or more roots of R(T,Y,Y') at T. */

/*                    *** TASK INTERRUPTED *** */
/*                Reported by negative values of IDID */

/*           IDID = -1 -- A large amount of work has been expended */
/*                     (about 500 steps). */

/*           IDID = -2 -- The error tolerances are too stringent. */

/*           IDID = -3 -- The local error test cannot be satisfied */
/*                     because you specified a zero component in ATOL */
/*                     and the corresponding computed solution component */
/*                     is zero.  Thus, a pure relative error test is */
/*                     impossible for this component. */

/*           IDID = -5 -- There were repeated failures in the evaluation */
/*                     or processing of the preconditioner (in JAC). */

/*           IDID = -6 -- DDASKR had repeated error test failures on the */
/*                     last attempted step. */

/*           IDID = -7 -- The nonlinear system solver in the time */
/*                     integration could not converge. */

/*           IDID = -8 -- The matrix of partial derivatives appears */
/*                     to be singular (direct method). */

/*           IDID = -9 -- The nonlinear system solver in the integration */
/*                     failed to achieve convergence, and there were */
/*                     repeated  error test failures in this step. */

/*           IDID =-10 -- The nonlinear system solver in the integration */
/*                     failed to achieve convergence because IRES was */
/*                     equal  to -1. */

/*           IDID =-11 -- IRES = -2 was encountered and control is */
/*                     being returned to the calling program. */

/*           IDID =-12 -- DDASKR failed to compute the initial Y, YPRIME. */

/*           IDID =-13 -- An unrecoverable error was encountered inside */
/*                     the user's PSOL routine, and control is being */
/*                     returned to the calling program. */

/*           IDID =-14 -- The Krylov linear system solver could not */
/*                     achieve convergence. */

/*           IDID =-15,..,-32 -- Not applicable for this code. */

/*                    *** TASK TERMINATED *** */
/*                reported by the value of IDID=-33 */

/*           IDID = -33 -- The code has encountered trouble from which */
/*                   it cannot recover.  A message is printed */
/*                   explaining the trouble and control is returned */
/*                   to the calling program.  For example, this occurs */
/*                   when invalid input is detected. */

/*   RTOL, ATOL -- these quantities remain unchanged except when */
/*               IDID = -2.  In this case, the error tolerances have been */
/*               increased by the code to values which are estimated to */
/*               be appropriate for continuing the integration.  However, */
/*               the reported solution at T was obtained using the input */
/*               values of RTOL and ATOL. */

/*   RWORK, IWORK -- contain information which is usually of no interest */
/*               to the user but necessary for subsequent calls. */
/*               However, you may be interested in the performance data */
/*               listed below.  These quantities are accessed in RWORK */
/*               and IWORK but have internal mnemonic names, as follows.. */

/*               RWORK(3)--contains H, the step size h to be attempted */
/*                        on the next step. */

/*               RWORK(4)--contains TN, the current value of the */
/*                        independent variable, i.e. the farthest point */
/*                        integration has reached.  This will differ */
/*                        from T if interpolation has been performed */
/*                        (IDID = 3). */

/*               RWORK(7)--contains HOLD, the stepsize used on the last */
/*                        successful step.  If INFO(11) = INFO(14) = 1, */
/*                        this contains the value of H used in the */
/*                        initial condition calculation. */

/*               IWORK(7)--contains K, the order of the method to be */
/*                        attempted on the next step. */

/*               IWORK(8)--contains KOLD, the order of the method used */
/*                        on the last step. */

/*               IWORK(11)--contains NST, the number of steps (in T) */
/*                        taken so far. */

/*               IWORK(12)--contains NRE, the number of calls to RES */
/*                        so far. */

/*               IWORK(13)--contains NJE, the number of calls to JAC so */
/*                        far (Jacobian or preconditioner evaluations). */

/*               IWORK(14)--contains NETF, the total number of error test */
/*                        failures so far. */

/*               IWORK(15)--contains NCFN, the total number of nonlinear */
/*                        convergence failures so far (includes counts */
/*                        of singular iteration matrix or singular */
/*                        preconditioners). */

/*               IWORK(16)--contains NCFL, the number of convergence */
/*                        failures of the linear iteration so far. */

/*               IWORK(17)--contains LENIW, the length of IWORK actually */
/*                        required.  This is defined on normal returns */
/*                        and on an illegal input return for */
/*                        insufficient storage. */

/*               IWORK(18)--contains LENRW, the length of RWORK actually */
/*                        required.  This is defined on normal returns */
/*                        and on an illegal input return for */
/*                        insufficient storage. */

/*               IWORK(19)--contains NNI, the total number of nonlinear */
/*                        iterations so far (each of which calls a */
/*                        linear solver). */

/*               IWORK(20)--contains NLI, the total number of linear */
/*                        (Krylov) iterations so far. */

/*               IWORK(21)--contains NPS, the number of PSOL calls so */
/*                        far, for preconditioning solve operations or */
/*                        for solutions with the user-supplied method. */

/*               IWORK(36)--contains the total number of calls to the */
/*                        constraint function routine RT so far. */

/*               Note: The various counters in IWORK do not include */
/*               counts during a prior call made with INFO(11) > 0 and */
/*               INFO(14) = 1. */


/* ------INPUT - WHAT TO DO TO CONTINUE THE INTEGRATION  ----------------- */
/*              (CALLS AFTER THE FIRST) */

/*     This code is organized so that subsequent calls to continue the */
/*     integration involve little (if any) additional effort on your */
/*     part.  You must monitor the IDID parameter in order to determine */
/*     what to do next. */

/*     Recalling that the principal task of the code is to integrate */
/*     from T to TOUT (the interval mode), usually all you will need */
/*     to do is specify a new TOUT upon reaching the current TOUT. */

/*     Do not alter any quantity not specifically permitted below.  In */
/*     particular do not alter NEQ, T, Y(*), YPRIME(*), RWORK(*), */
/*     IWORK(*), or the differential equation in subroutine RES.  Any */
/*     such alteration constitutes a new problem and must be treated */
/*     as such, i.e. you must start afresh. */

/*     You cannot change from array to scalar error control or vice */
/*     versa (INFO(2)), but you can change the size of the entries of */
/*     RTOL or ATOL.  Increasing a tolerance makes the equation easier */
/*     to integrate.  Decreasing a tolerance will make the equation */
/*     harder to integrate and should generally be avoided. */

/*     You can switch from the intermediate-output mode to the */
/*     interval mode (INFO(3)) or vice versa at any time. */

/*     If it has been necessary to prevent the integration from going */
/*     past a point TSTOP (INFO(4), RWORK(1)), keep in mind that the */
/*     code will not integrate to any TOUT beyond the currently */
/*     specified TSTOP.  Once TSTOP has been reached, you must change */
/*     the value of TSTOP or set INFO(4) = 0.  You may change INFO(4) */
/*     or TSTOP at any time but you must supply the value of TSTOP in */
/*     RWORK(1) whenever you set INFO(4) = 1. */

/*     Do not change INFO(5), INFO(6), INFO(12-17) or their associated */
/*     IWORK/RWORK locations unless you are going to restart the code. */

/*                    *** FOLLOWING A COMPLETED TASK *** */

/*     If.. */
/*     IDID = 1, call the code again to continue the integration */
/*                  another step in the direction of TOUT. */

/*     IDID = 2 or 3, define a new TOUT and call the code again. */
/*                  TOUT must be different from T.  You cannot change */
/*                  the direction of integration without restarting. */

/*     IDID = 4, reset INFO(11) = 0 and call the code again to begin */
/*                  the integration.  (If you leave INFO(11) > 0 and */
/*                  INFO(14) = 1, you may generate an infinite loop.) */
/*                  In this situation, the next call to DDASKR is */
/*                  considered to be the first call for the problem, */
/*                  in that all initializations are done. */

/*     IDID = 5, call the code again to continue the integration in the */
/*                  direction of TOUT.  You may change the functions */
/*                  Ri defined by RT after a return with IDID = 5, but */
/*                  the number of constraint functions NRT must remain */
/*                  the same.  If you wish to change the functions in */
/*                  RES or in RT, then you must restart the code. */

/*                    *** FOLLOWING AN INTERRUPTED TASK *** */

/*     To show the code that you realize the task was interrupted and */
/*     that you want to continue, you must take appropriate action and */
/*     set INFO(1) = 1. */

/*     If.. */
/*     IDID = -1, the code has taken about 500 steps.  If you want to */
/*                  continue, set INFO(1) = 1 and call the code again. */
/*                  An additional 500 steps will be allowed. */


/*     IDID = -2, the error tolerances RTOL, ATOL have been increased */
/*                  to values the code estimates appropriate for */
/*                  continuing.  You may want to change them yourself. */
/*                  If you are sure you want to continue with relaxed */
/*                  error tolerances, set INFO(1) = 1 and call the code */
/*                  again. */

/*     IDID = -3, a solution component is zero and you set the */
/*                  corresponding component of ATOL to zero.  If you */
/*                  are sure you want to continue, you must first alter */
/*                  the error criterion to use positive values of ATOL */
/*                  for those components corresponding to zero solution */
/*                  components, then set INFO(1) = 1 and call the code */
/*                  again. */

/*     IDID = -4  --- cannot occur with this code. */

/*     IDID = -5, your JAC routine failed with the Krylov method.  Check */
/*                  for errors in JAC and restart the integration. */

/*     IDID = -6, repeated error test failures occurred on the last */
/*                  attempted step in DDASKR.  A singularity in the */
/*                  solution may be present.  If you are absolutely */
/*                  certain you want to continue, you should restart */
/*                  the integration.  (Provide initial values of Y and */
/*                  YPRIME which are consistent.) */

/*     IDID = -7, repeated convergence test failures occurred on the last */
/*                  attempted step in DDASKR.  An inaccurate or ill- */
/*                  conditioned Jacobian or preconditioner may be the */
/*                  problem.  If you are absolutely certain you want */
/*                  to continue, you should restart the integration. */


/*     IDID = -8, the matrix of partial derivatives is singular, with */
/*                  the use of direct methods.  Some of your equations */
/*                  may be redundant.  DDASKR cannot solve the problem */
/*                  as stated.  It is possible that the redundant */
/*                  equations could be removed, and then DDASKR could */
/*                  solve the problem.  It is also possible that a */
/*                  solution to your problem either does not exist */
/*                  or is not unique. */

/*     IDID = -9, DDASKR had multiple convergence test failures, preceded */
/*                  by multiple error test failures, on the last */
/*                  attempted step.  It is possible that your problem is */
/*                  ill-posed and cannot be solved using this code.  Or, */
/*                  there may be a discontinuity or a singularity in the */
/*                  solution.  If you are absolutely certain you want to */
/*                  continue, you should restart the integration. */

/*     IDID = -10, DDASKR had multiple convergence test failures */
/*                  because IRES was equal to -1.  If you are */
/*                  absolutely certain you want to continue, you */
/*                  should restart the integration. */

/*     IDID = -11, there was an unrecoverable error (IRES = -2) from RES */
/*                  inside the nonlinear system solver.  Determine the */
/*                  cause before trying again. */

/*     IDID = -12, DDASKR failed to compute the initial Y and YPRIME */
/*                  vectors.  This could happen because the initial */
/*                  approximation to Y or YPRIME was not very good, or */
/*                  because no consistent values of these vectors exist. */
/*                  The problem could also be caused by an inaccurate or */
/*                  singular iteration matrix, or a poor preconditioner. */

/*     IDID = -13, there was an unrecoverable error encountered inside */
/*                  your PSOL routine.  Determine the cause before */
/*                  trying again. */

/*     IDID = -14, the Krylov linear system solver failed to achieve */
/*                  convergence.  This may be due to ill-conditioning */
/*                  in the iteration matrix, or a singularity in the */
/*                  preconditioner (if one is being used). */
/*                  Another possibility is that there is a better */
/*                  choice of Krylov parameters (see INFO(13)). */
/*                  Possibly the failure is caused by redundant equations */
/*                  in the system, or by inconsistent equations. */
/*                  In that case, reformulate the system to make it */
/*                  consistent and non-redundant. */

/*     IDID = -15,..,-32 --- Cannot occur with this code. */

/*                       *** FOLLOWING A TERMINATED TASK *** */

/*     If IDID = -33, you cannot continue the solution of this problem. */
/*                  An attempt to do so will result in your run being */
/*                  terminated. */

/*  --------------------------------------------------------------------- */

/* ***REFERENCES */
/*  1.  L. R. Petzold, A Description of DASSL: A Differential/Algebraic */
/*      System Solver, in Scientific Computing, R. S. Stepleman et al. */
/*      (Eds.), North-Holland, Amsterdam, 1983, pp. 65-68. */
/*  2.  K. E. Brenan, S. L. Campbell, and L. R. Petzold, Numerical */
/*      Solution of Initial-Value Problems in Differential-Algebraic */
/*      Equations, Elsevier, New York, 1989. */
/*  3.  P. N. Brown and A. C. Hindmarsh, Reduced Storage Matrix Methods */
/*      in Stiff ODE Systems, J. Applied Mathematics and Computation, */
/*      31 (1989), pp. 40-91. */
/*  4.  P. N. Brown, A. C. Hindmarsh, and L. R. Petzold, Using Krylov */
/*      Methods in the Solution of Large-Scale Differential-Algebraic */
/*      Systems, SIAM J. Sci. Comp., 15 (1994), pp. 1467-1488. */
/*  5.  P. N. Brown, A. C. Hindmarsh, and L. R. Petzold, Consistent */
/*      Initial Condition Calculation for Differential-Algebraic */
/*      Systems, SIAM J. Sci. Comp. 19 (1998), pp. 1495-1512. */

/* ***ROUTINES CALLED */

/*   The following are all the subordinate routines used by DDASKR. */

/*   DRCHEK does preliminary checking for roots, and serves as an */
/*          interface between Subroutine DDASKR and Subroutine DROOTS. */
/*   DROOTS finds the leftmost root of a set of functions. */
/*   DDASIC computes consistent initial conditions. */
/*   DYYPNW updates Y and YPRIME in linesearch for initial condition */
/*          calculation. */
/*   DDSTP  carries out one step of the integration. */
/*   DCNSTR/DCNST0 check the current solution for constraint violations. */
/*   DDAWTS sets error weight quantities. */
/*   DINVWT tests and inverts the error weights. */
/*   DDATRP performs interpolation to get an output solution. */
/*   DDWNRM computes the weighted root-mean-square norm of a vector. */
/*   D1MACH provides the unit roundoff of the computer. */
/*   XERRWD/XSETF/XSETUN/IXSAV is a package to handle error messages. */
/*   DDASID nonlinear equation driver to initialize Y and YPRIME using */
/*          direct linear system solver methods.  Interfaces to Newton */
/*          solver (direct case). */
/*   DNSID  solves the nonlinear system for unknown initial values by */
/*          modified Newton iteration and direct linear system methods. */
/*   DLINSD carries out linesearch algorithm for initial condition */
/*          calculation (direct case). */
/*   DFNRMD calculates weighted norm of preconditioned residual in */
/*          initial condition calculation (direct case). */
/*   DNEDD  nonlinear equation driver for direct linear system solver */
/*          methods.  Interfaces to Newton solver (direct case). */
/*   DMATD  assembles the iteration matrix (direct case). */
/*   DNSD   solves the associated nonlinear system by modified */
/*          Newton iteration and direct linear system methods. */
/*   DSLVD  interfaces to linear system solver (direct case). */
/*   DDASIK nonlinear equation driver to initialize Y and YPRIME using */
/*          Krylov iterative linear system methods.  Interfaces to */
/*          Newton solver (Krylov case). */
/*   DNSIK  solves the nonlinear system for unknown initial values by */
/*          Newton iteration and Krylov iterative linear system methods. */
/*   DLINSK carries out linesearch algorithm for initial condition */
/*          calculation (Krylov case). */
/*   DFNRMK calculates weighted norm of preconditioned residual in */
/*          initial condition calculation (Krylov case). */
/*   DNEDK  nonlinear equation driver for iterative linear system solver */
/*          methods.  Interfaces to Newton solver (Krylov case). */
/*   DNSK   solves the associated nonlinear system by Inexact Newton */
/*          iteration and (linear) Krylov iteration. */
/*   DSLVK  interfaces to linear system solver (Krylov case). */
/*   DSPIGM solves a linear system by SPIGMR algorithm. */
/*   DATV   computes matrix-vector product in Krylov algorithm. */
/*   DORTH  performs orthogonalization of Krylov basis vectors. */
/*   DHEQR  performs QR factorization of Hessenberg matrix. */
/*   DHELS  finds least-squares solution of Hessenberg linear system. */
/*   DGEFA, DGESL, DGBFA, DGBSL are LINPACK routines for solving */
/*          linear systems (dense or band direct methods). */
/*   DAXPY, DCOPY, DDOT, DNRM2, DSCAL are Basic Linear Algebra (BLAS) */
/*          routines. */

/* The routines called directly by DDASKR are: */
/*   DCNST0, DDAWTS, DINVWT, D1MACH, DDWNRM, DDASIC, DDATRP, DDSTP, */
/*   DRCHEK, XERRWD */

/* ***END PROLOGUE DDASKR */



/*     Set pointers into IWORK. */


/*     Set pointers into RWORK. */




/* ***FIRST EXECUTABLE STATEMENT  DDASKR */


    /* Parameter adjustments */
    --y;
    --yprime;
    --info;
    --rtol;
    --atol;
    --rwork;
    --iwork;
    --rpar;
    --ipar;

    /* Function Body */
    if (info[1] != 0) {
	goto L100;
    }

/* ----------------------------------------------------------------------- */
/*     This block is executed for the initial call only. */
/*     It contains checking of inputs and initializations. */
/* ----------------------------------------------------------------------- */

/*     First check INFO array to make sure all elements of INFO */
/*     Are within the proper range.  (INFO(1) is checked later, because */
/*     it must be tested on every call.) ITEMP holds the location */
/*     within INFO which may be out of range. */

    for (i__ = 2; i__ <= 9; ++i__) {
	itemp = i__;
	if (info[i__] != 0 && info[i__] != 1) {
	    goto L701;
	}
/* L10: */
    }
    itemp = 10;
    if (info[10] < 0 || info[10] > 3) {
	goto L701;
    }
    itemp = 11;
    if (info[11] < 0 || info[11] > 2) {
	goto L701;
    }
    for (i__ = 12; i__ <= 17; ++i__) {
	itemp = i__;
	if (info[i__] != 0 && info[i__] != 1) {
	    goto L701;
	}
/* L15: */
    }
    itemp = 18;
    if (info[18] < 0 || info[18] > 2) {
	goto L701;
    }

/*     Check NEQ to see if it is positive. */

    if (*neq <= 0) {
	goto L702;
    }

/*     Check and compute maximum order. */

    mxord = 5;
    if (info[9] != 0) {
	mxord = iwork[3];
	if (mxord < 1 || mxord > 5) {
	    goto L703;
	}
    }
    iwork[3] = mxord;

/*     Set and/or check inputs for constraint checking (INFO(10) .NE. 0). */
/*     Set values for ICNFLG, NONNEG, and pointer LID. */

    icnflg = 0;
    nonneg = 0;
    lid = 41;
    if (info[10] == 0) {
	goto L20;
    }
    if (info[10] == 1) {
	icnflg = 1;
	nonneg = 0;
	lid = *neq + 41;
    } else if (info[10] == 2) {
	icnflg = 0;
	nonneg = 1;
    } else {
	icnflg = 1;
	nonneg = 1;
	lid = *neq + 41;
    }

L20:

/*     Set and/or check inputs for Krylov solver (INFO(12) .NE. 0). */
/*     If indicated, set default values for MAXL, KMP, NRMAX, and EPLI. */
/*     Otherwise, verify inputs required for iterative solver. */

    if (info[12] == 0) {
	goto L25;
    }

    iwork[23] = info[12];
    if (info[13] == 0) {
	iwork[24] = MIN(5,*neq);
	iwork[25] = iwork[24];
	iwork[26] = 5;
	rwork[10] = .05;
    } else {
	if (iwork[24] < 1 || iwork[24] > *neq) {
	    goto L720;
	}
	if (iwork[25] < 1 || iwork[25] > iwork[24]) {
	    goto L721;
	}
	if (iwork[26] < 0) {
	    goto L722;
	}
	if (rwork[10] <= 0. || rwork[10] >= 1.) {
	    goto L723;
	}
    }

L25:

/*     Set and/or check controls for the initial condition calculation */
/*     (INFO(11) .GT. 0).  If indicated, set default values. */
/*     Otherwise, verify inputs required for iterative solver. */

    if (info[11] == 0) {
	goto L30;
    }
    if (info[17] == 0) {
	iwork[32] = 5;
	if (info[12] > 0) {
	    iwork[32] = 15;
	}
	iwork[33] = 6;
	if (info[12] > 0) {
	    iwork[33] = 2;
	}
	iwork[34] = 5;
	iwork[35] = 0;
	rwork[15] = .01;
    } else {
	if (iwork[32] <= 0) {
	    goto L725;
	}
	if (iwork[33] <= 0) {
	    goto L725;
	}
	if (iwork[34] <= 0) {
	    goto L725;
	}
	lsoff = iwork[35];
	if (lsoff < 0 || lsoff > 1) {
	    goto L725;
	}
	if (rwork[15] <= 0.) {
	    goto L725;
	}
    }

L30:

/*     Below is the computation and checking of the work array lengths */
/*     LENIW and LENRW, using direct methods (INFO(12) = 0) or */
/*     the Krylov methods (INFO(12) = 1). */

    lenic = 0;
    if (info[10] == 1 || info[10] == 3) {
	lenic = *neq;
    }
    lenid = 0;
    if (info[11] == 1 || info[16] == 1) {
	lenid = *neq;
    }
    if (info[12] == 0) {

/*        Compute MTYPE, etc.  Check ML and MU. */

/* Computing MAX */
	i__1 = mxord + 1;
	ncphi = MAX(i__1,4);
	if (info[6] == 0) {
/* Computing 2nd power */
	    i__1 = *neq;
	    lenpd = i__1 * i__1;
	    lenrw = *nrt * 3 + 60 + (ncphi + 3) * *neq + lenpd;
	    if (info[5] == 0) {
		iwork[4] = 2;
	    } else {
		iwork[4] = 1;
	    }
	} else {
	    if (iwork[1] < 0 || iwork[1] >= *neq) {
		goto L717;
	    }
	    if (iwork[2] < 0 || iwork[2] >= *neq) {
		goto L718;
	    }
	    lenpd = ((iwork[1] << 1) + iwork[2] + 1) * *neq;
	    if (info[5] == 0) {
		iwork[4] = 5;
		mband = iwork[1] + iwork[2] + 1;
		msave = *neq / mband + 1;
		lenrw = *nrt * 3 + 60 + (ncphi + 3) * *neq + lenpd + (msave <<
			 1);
	    } else {
		iwork[4] = 4;
		lenrw = *nrt * 3 + 60 + (ncphi + 3) * *neq + lenpd;
	    }
	}

/*        Compute LENIW, LENWP, LENIWP. */

	leniw = lenic + 40 + lenid + *neq;
	lenwp = 0;
	leniwp = 0;

    } else if (info[12] == 1) {
	ncphi = mxord + 1;
	maxl = iwork[24];
	lenwp = iwork[27];
	leniwp = iwork[28];
/* Computing MIN */
	i__1 = 1, i__2 = maxl - iwork[25];
	lenpd = (maxl + 3 + MIN(i__1,i__2)) * *neq + (maxl + 3) * maxl + 1 +
		lenwp;
	lenrw = *nrt * 3 + 60 + (mxord + 5) * *neq + lenpd;
	leniw = lenic + 40 + lenid + leniwp;

    }
    if (info[16] != 0) {
	lenrw += *neq;
    }

/*     Check lengths of RWORK and IWORK. */

    iwork[17] = leniw;
    iwork[18] = lenrw;
    iwork[22] = lenpd;
    iwork[29] = lenpd - lenwp + 1;
    if (*lrw < lenrw) {
	goto L704;
    }
    if (*liw < leniw) {
	goto L705;
    }

/*     Check ICNSTR for legality. */

    if (lenic > 0) {
	i__1 = *neq;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ici = iwork[i__ + 40];
	    if (ici < -2 || ici > 2) {
		goto L726;
	    }
/* L40: */
	}
    }

/*     Check Y for consistency with constraints. */

    if (lenic > 0) {
	_daskr_dcnst0_(neq, &y[1], &iwork[41], &iret);
	if (iret != 0) {
	    goto L727;
	}
    }

/*     Check ID for legality and set INDEX = 0 or 1. */

    index = 1;
    if (lenid > 0) {
	index = 0;
	i__1 = *neq;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    idi = iwork[lid - 1 + i__];
	    if (idi != 1 && idi != -1) {
		goto L724;
	    }
	    if (idi == -1) {
		index = 1;
	    }
/* L50: */
	}
    }

/*     Check to see that TOUT is different from T, and NRT .ge. 0. */

    if (*tout == *t) {
	goto L719;
    }
    if (*nrt < 0) {
	goto L730;
    }

/*     Check HMAX. */

    if (info[7] != 0) {
	hmax = rwork[2];
	if (hmax <= 0.) {
	    goto L710;
	}
    }

/*     Initialize counters and other flags. */

    iwork[11] = 0;
    iwork[12] = 0;
    iwork[13] = 0;
    iwork[14] = 0;
    iwork[15] = 0;
    iwork[19] = 0;
    iwork[20] = 0;
    iwork[21] = 0;
    iwork[16] = 0;
    iwork[36] = 0;
    iwork[31] = info[18];
    *idid = 1;
    goto L200;

/* ----------------------------------------------------------------------- */
/*     This block is for continuation calls only. */
/*     Here we check INFO(1), and if the last step was interrupted, */
/*     we check whether appropriate action was taken. */
/* ----------------------------------------------------------------------- */

L100:
    if (info[1] == 1) {
	goto L110;
    }
    itemp = 1;
    if (info[1] != -1) {
	goto L701;
    }

/*     If we are here, the last step was interrupted by an error */
/*     condition from DDSTP, and appropriate action was not taken. */
/*     This is a fatal error. */

    _daskr_str_copy(msg, "DASKR--  THE LAST STEP TERMINATED WITH A NEGATIVE", (integer)
	    80, (integer)49);
    _daskr_xerrwd_(msg, &c__49, &c__201, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    _daskr_str_copy(msg, "DASKR--  VALUE (=I1) OF IDID AND NO APPROPRIATE", (integer)80,
	     (integer)47);
    _daskr_xerrwd_(msg, &c__47, &c__202, &c__0, &c__1, idid, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    _daskr_str_copy(msg, "DASKR--  ACTION WAS TAKEN. RUN TERMINATED", (integer)80, (
	    integer)41);
    _daskr_xerrwd_(msg, &c__41, &c__203, &c__1, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    return 0;
L110:

/* ----------------------------------------------------------------------- */
/*     This block is executed on all calls. */

/*     Counters are saved for later checks of performance. */
/*     Then the error tolerance parameters are checked, and the */
/*     work array pointers are set. */
/* ----------------------------------------------------------------------- */

L200:

/*     Save counters for use later. */

    iwork[10] = iwork[11];
    nli0 = iwork[20];
    nni0 = iwork[19];
    ncfn0 = iwork[15];
    ncfl0 = iwork[16];
    nwarn = 0;

/*     Check RTOL and ATOL. */

    nzflg = 0;
    rtoli = rtol[1];
    atoli = atol[1];
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (info[2] == 1) {
	    rtoli = rtol[i__];
	}
	if (info[2] == 1) {
	    atoli = atol[i__];
	}
	if (rtoli > 0. || atoli > 0.) {
	    nzflg = 1;
	}
	if (rtoli < 0.) {
	    goto L706;
	}
	if (atoli < 0.) {
	    goto L707;
	}
/* L210: */
    }
    if (nzflg == 0) {
	goto L708;
    }

/*     Set pointers to RWORK and IWORK segments. */
/*     For direct methods, SAVR is not used. */

    iwork[30] = lid + lenid;
    lsavr = 61;
    if (info[12] != 0) {
	lsavr = *neq + 61;
    }
    le = lsavr + *neq;
    lwt = le + *neq;
    lvt = lwt;
    if (info[16] != 0) {
	lvt = lwt + *neq;
    }
    lphi = lvt + *neq;
    lr0 = lphi + ncphi * *neq;
    lr1 = lr0 + *nrt;
    lrx = lr1 + *nrt;
    lwm = lrx + *nrt;
    if (info[1] == 1) {
	goto L400;
    }

/* ----------------------------------------------------------------------- */
/*     This block is executed on the initial call only. */
/*     Set the initial step size, the error weight vector, and PHI. */
/*     Compute unknown initial components of Y and YPRIME, if requested. */
/* ----------------------------------------------------------------------- */

/* L300: */
    tn = *t;
    *idid = 1;

/*     Set error weight array WT and altered weight array VT. */

    _daskr_ddawts_(neq, &info[2], &rtol[1], &atol[1], &y[1], &rwork[lwt], &rpar[1], &
	    ipar[1]);
    _daskr_dinvwt_(neq, &rwork[lwt], &ier);
    if (ier != 0) {
	goto L713;
    }
    if (info[16] != 0) {
	i__1 = *neq;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L305: */
/* Computing MAX */
	    i__2 = iwork[lid + i__ - 1];
	    rwork[lvt + i__ - 1] = MAX(i__2,0) * rwork[lwt + i__ - 1];
	}
    }

/*     Compute unit roundoff and HMIN. */

    uround = _daskr_d1mach_(&c__4);
    rwork[9] = uround;
/* Computing MAX */
    d__1 = fabs(*t), d__2 = fabs(*tout);
    hmin = uround * 4. * MAX(d__1,d__2);

/*     Set/check STPTOL control for initial condition calculation. */

    if (info[11] != 0) {
	if (info[17] == 0) {
	    rwork[14] = _daskr_real_pow(&uround, &c_b68);
	} else {
	    if (rwork[14] <= 0.) {
		goto L725;
	    }
	}
    }

/*     Compute EPCON and square root of NEQ and its reciprocal, used */
/*     inside iterative solver. */

    rwork[13] = .33;
    floatn = (real_number) (*neq);
    rwork[11] = sqrt(floatn);
    rwork[12] = 1. / rwork[11];

/*     Check initial interval to see that it is long enough. */

    tdist = (d__1 = *tout - *t, fabs(d__1));
    if (tdist < hmin) {
	goto L714;
    }

/*     Check H0, if this was input. */

    if (info[8] == 0) {
	goto L310;
    }
    h0 = rwork[3];
    if ((*tout - *t) * h0 < 0.) {
	goto L711;
    }
    if (h0 == 0.) {
	goto L712;
    }
    goto L320;
L310:

/*     Compute initial stepsize, to be used by either */
/*     DDSTP or DDASIC, depending on INFO(11). */

    h0 = tdist * .001;
    ypnorm = _daskr_ddwnrm_(neq, &yprime[1], &rwork[lvt], &rpar[1], &ipar[1]);
    if (ypnorm > .5 / h0) {
	h0 = .5 / ypnorm;
    }
    d__1 = *tout - *t;
    h0 = _daskr_real_sign(&h0, &d__1);

/*     Adjust H0 if necessary to meet HMAX bound. */

L320:
    if (info[7] == 0) {
	goto L330;
    }
    rh = fabs(h0) / rwork[2];
    if (rh > 1.) {
	h0 /= rh;
    }

/*     Check against TSTOP, if applicable. */

L330:
    if (info[4] == 0) {
	goto L340;
    }
    tstop = rwork[1];
    if ((tstop - *t) * h0 < 0.) {
	goto L715;
    }
    if ((*t + h0 - tstop) * h0 > 0.) {
	h0 = tstop - *t;
    }
    if ((tstop - *tout) * h0 < 0.) {
	goto L709;
    }

L340:
    if (info[11] == 0) {
	goto L370;
    }

/*     Compute unknown components of initial Y and YPRIME, depending */
/*     on INFO(11) and INFO(12).  INFO(12) represents the nonlinear */
/*     solver type (direct/Krylov).  Pass the name of the specific */
/*     nonlinear solver, depending on INFO(12).  The location of the work */
/*     arrays SAVR, YIC, YPIC, PWK also differ in the two cases. */
/*     For use in stopping tests, pass TSCALE = TDIST if INDEX = 0. */

    nwt = 1;
    epconi = rwork[15] * rwork[13];
    tscale = 0.;
    if (index == 0) {
	tscale = tdist;
    }
L350:
    if (info[12] == 0) {
	lyic = lphi + (*neq << 1);
	lypic = lyic + *neq;
	lpwk = lypic;
	_daskr_ddasic_(&tn, &y[1], &yprime[1], neq, &info[11], &iwork[lid], (Unknown_fp)
		res, (Unknown_fp)jac, (Unknown_fp)psol, &h0, &tscale, &rwork[lwt], &nwt,
		idid, &rpar[1], &ipar[1], &rwork[lphi], &rwork[lsavr], &rwork[
		61], &rwork[le], &rwork[lyic], &rwork[lypic], &rwork[lpwk], &
		rwork[lwm], &iwork[1], &rwork[9], &rwork[10], &rwork[11], &
		rwork[12], &epconi, &rwork[14], &info[15], &icnflg, &iwork[41]
		, (Unknown_fp)_daskr_ddasid_);
    } else if (info[12] == 1) {
	lyic = lwm;
	lypic = lyic + *neq;
	lpwk = lypic + *neq;
	_daskr_ddasic_(&tn, &y[1], &yprime[1], neq, &info[11], &iwork[lid], (Unknown_fp)
		res, (Unknown_fp)jac, (Unknown_fp)psol, &h0, &tscale, &rwork[lwt], &nwt,
		idid, &rpar[1], &ipar[1], &rwork[lphi], &rwork[lsavr], &rwork[
		61], &rwork[le], &rwork[lyic], &rwork[lypic], &rwork[lpwk], &
		rwork[lwm], &iwork[1], &rwork[9], &rwork[10], &rwork[11], &
		rwork[12], &epconi, &rwork[14], &info[15], &icnflg, &iwork[41]
		, (Unknown_fp)_daskr_ddasik_);
    }

    if (*idid < 0) {
	goto L600;
    }

/*     DDASIC was successful.  If this was the first call to DDASIC, */
/*     update the WT array (with the current Y) and call it again. */

    if (nwt == 2) {
	goto L355;
    }
    nwt = 2;
    _daskr_ddawts_(neq, &info[2], &rtol[1], &atol[1], &y[1], &rwork[lwt], &rpar[1], &
	    ipar[1]);
    _daskr_dinvwt_(neq, &rwork[lwt], &ier);
    if (ier != 0) {
	goto L713;
    }
    goto L350;

/*     If INFO(14) = 1, return now with IDID = 4. */

L355:
    if (info[14] == 1) {
	*idid = 4;
	h__ = h0;
	if (info[11] == 1) {
	    rwork[7] = h0;
	}
	goto L590;
    }

/*     Update the WT and VT arrays one more time, with the new Y. */

    _daskr_ddawts_(neq, &info[2], &rtol[1], &atol[1], &y[1], &rwork[lwt], &rpar[1], &
	    ipar[1]);
    _daskr_dinvwt_(neq, &rwork[lwt], &ier);
    if (ier != 0) {
	goto L713;
    }
    if (info[16] != 0) {
	i__2 = *neq;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L357: */
/* Computing MAX */
	    i__1 = iwork[lid + i__ - 1];
	    rwork[lvt + i__ - 1] = MAX(i__1,0) * rwork[lwt + i__ - 1];
	}
    }

/*     Reset the initial stepsize to be used by DDSTP. */
/*     Use H0, if this was input.  Otherwise, recompute H0, */
/*     and adjust it if necessary to meet HMAX bound. */

    if (info[8] != 0) {
	h0 = rwork[3];
	goto L360;
    }

    h0 = tdist * .001;
    ypnorm = _daskr_ddwnrm_(neq, &yprime[1], &rwork[lvt], &rpar[1], &ipar[1]);
    if (ypnorm > .5 / h0) {
	h0 = .5 / ypnorm;
    }
    d__1 = *tout - *t;
    h0 = _daskr_real_sign(&h0, &d__1);

L360:
    if (info[7] != 0) {
	rh = fabs(h0) / rwork[2];
	if (rh > 1.) {
	    h0 /= rh;
	}
    }

/*     Check against TSTOP, if applicable. */

    if (info[4] != 0) {
	tstop = rwork[1];
	if ((*t + h0 - tstop) * h0 > 0.) {
	    h0 = tstop - *t;
	}
    }

/*     Load H and RWORK(LH) with H0. */

L370:
    h__ = h0;
    rwork[3] = h__;

/*     Load Y and H*YPRIME into PHI(*,1) and PHI(*,2). */

    itemp = lphi + *neq;
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rwork[lphi + i__ - 1] = y[i__];
/* L380: */
	rwork[itemp + i__ - 1] = h__ * yprime[i__];
    }

/*     Initialize T0 in RWORK; check for a zero of R near initial T. */

    rwork[51] = *t;
    iwork[37] = 0;
    rwork[39] = h__;
    rwork[40] = h__ * 2.;
    iwork[8] = 1;
    if (*nrt == 0) {
	goto L390;
    }
    _daskr_drchek_(&c__1, (Unknown_fp)rt, nrt, neq, t, tout, &y[1], &yprime[1], &rwork[
	    lphi], &rwork[39], &iwork[8], &rwork[lr0], &rwork[lr1], &rwork[
	    lrx], jroot, &irt, &rwork[9], &info[3], &rwork[1], &iwork[1], &
	    rpar[1], &ipar[1]);
    if (irt < 0) {
	goto L731;
    }

L390:
    goto L500;

/* ----------------------------------------------------------------------- */
/*     This block is for continuation calls only. */
/*     Its purpose is to check stop conditions before taking a step. */
/*     Adjust H if necessary to meet HMAX bound. */
/* ----------------------------------------------------------------------- */

L400:
    uround = rwork[9];
    done = _FALSE_;
    tn = rwork[4];
    h__ = rwork[3];
    if (*nrt == 0) {
	goto L405;
    }

/*     Check for a zero of R near TN. */

    _daskr_drchek_(&c__2, (Unknown_fp)rt, nrt, neq, &tn, tout, &y[1], &yprime[1], &rwork[
	    lphi], &rwork[39], &iwork[8], &rwork[lr0], &rwork[lr1], &rwork[
	    lrx], jroot, &irt, &rwork[9], &info[3], &rwork[1], &iwork[1], &
	    rpar[1], &ipar[1]);
    if (irt < 0) {
	goto L731;
    }
    if (irt != 1) {
	goto L405;
    }
    iwork[37] = 1;
    *idid = 5;
    *t = rwork[51];
    done = _TRUE_;
    goto L490;
L405:

    if (info[7] == 0) {
	goto L410;
    }
    rh = fabs(h__) / rwork[2];
    if (rh > 1.) {
	h__ /= rh;
    }
L410:
    if (*t == *tout) {
	goto L719;
    }
    if ((*t - *tout) * h__ > 0.) {
	goto L711;
    }
    if (info[4] == 1) {
	goto L430;
    }
    if (info[3] == 1) {
	goto L420;
    }
    if ((tn - *tout) * h__ < 0.) {
	goto L490;
    }
    _daskr_ddatrp_(&tn, tout, &y[1], &yprime[1], neq, &iwork[8], &rwork[lphi], &
	    rwork[39]);
    *t = *tout;
    *idid = 3;
    done = _TRUE_;
    goto L490;
L420:
    if ((tn - *t) * h__ <= 0.) {
	goto L490;
    }
    if ((tn - *tout) * h__ >= 0.) {
	goto L425;
    }
    _daskr_ddatrp_(&tn, &tn, &y[1], &yprime[1], neq, &iwork[8], &rwork[lphi], &rwork[
	    39]);
    *t = tn;
    *idid = 1;
    done = _TRUE_;
    goto L490;
L425:
    _daskr_ddatrp_(&tn, tout, &y[1], &yprime[1], neq, &iwork[8], &rwork[lphi], &
	    rwork[39]);
    *t = *tout;
    *idid = 3;
    done = _TRUE_;
    goto L490;
L430:
    if (info[3] == 1) {
	goto L440;
    }
    tstop = rwork[1];
    if ((tn - tstop) * h__ > 0.) {
	goto L715;
    }
    if ((tstop - *tout) * h__ < 0.) {
	goto L709;
    }
    if ((tn - *tout) * h__ < 0.) {
	goto L450;
    }
    _daskr_ddatrp_(&tn, tout, &y[1], &yprime[1], neq, &iwork[8], &rwork[lphi], &
	    rwork[39]);
    *t = *tout;
    *idid = 3;
    done = _TRUE_;
    goto L490;
L440:
    tstop = rwork[1];
    if ((tn - tstop) * h__ > 0.) {
	goto L715;
    }
    if ((tstop - *tout) * h__ < 0.) {
	goto L709;
    }
    if ((tn - *t) * h__ <= 0.) {
	goto L450;
    }
    if ((tn - *tout) * h__ >= 0.) {
	goto L445;
    }
    _daskr_ddatrp_(&tn, &tn, &y[1], &yprime[1], neq, &iwork[8], &rwork[lphi], &rwork[
	    39]);
    *t = tn;
    *idid = 1;
    done = _TRUE_;
    goto L490;
L445:
    _daskr_ddatrp_(&tn, tout, &y[1], &yprime[1], neq, &iwork[8], &rwork[lphi], &
	    rwork[39]);
    *t = *tout;
    *idid = 3;
    done = _TRUE_;
    goto L490;
L450:

/*     Check whether we are within roundoff of TSTOP. */

    if ((d__1 = tn - tstop, fabs(d__1)) > uround * 100. * (fabs(tn) + fabs(h__)))
	     {
	goto L460;
    }
    _daskr_ddatrp_(&tn, &tstop, &y[1], &yprime[1], neq, &iwork[8], &rwork[lphi], &
	    rwork[39]);
    *idid = 2;
    *t = tstop;
    done = _TRUE_;
    goto L490;
L460:
    tnext = tn + h__;
    if ((tnext - tstop) * h__ <= 0.) {
	goto L490;
    }
    h__ = tstop - tn;
    rwork[3] = h__;

L490:
    if (done) {
	goto L590;
    }

/* ----------------------------------------------------------------------- */
/*     The next block contains the call to the one-step integrator DDSTP. */
/*     This is a looping point for the integration steps. */
/*     Check for too many steps. */
/*     Check for poor Newton/Krylov performance. */
/*     Update WT.  Check for too much accuracy requested. */
/*     Compute minimum stepsize. */
/* ----------------------------------------------------------------------- */

L500:

/*     Check for too many steps. */

    if (iwork[11] - iwork[10] < 500) {
	goto L505;
    }
    *idid = -1;
    goto L527;

/* Check for poor Newton/Krylov performance. */

L505:
    if (info[12] == 0) {
	goto L510;
    }
    nstd = iwork[11] - iwork[10];
    nnid = iwork[19] - nni0;
    if (nstd < 10 || nnid == 0) {
	goto L510;
    }
    avlin = (real_number) (iwork[20] - nli0) / (real_number) nnid;
    rcfn = (real_number) (iwork[15] - ncfn0) / (real_number) nstd;
    rcfl = (real_number) (iwork[16] - ncfl0) / (real_number) nnid;
    fmaxl = (real_number) iwork[24];
    lavl = avlin > fmaxl;
    lcfn = rcfn > .9;
    lcfl = rcfl > .9;
    lwarn = lavl || lcfn || lcfl;
    if (! lwarn) {
	goto L510;
    }
    ++nwarn;
    if (nwarn > 10) {
	goto L510;
    }
    if (lavl) {
	_daskr_str_copy(msg, "DASKR-- Warning. Poor iterative algorithm performance   "
		, (integer)80, (integer)56);
	_daskr_xerrwd_(msg, &c__56, &c__501, &c__0, &c__0, &c__0, &c__0, &c__0, &
		c_b38, &c_b38, (integer)80);
	_daskr_str_copy(msg, "      at T = R1. Average no. of linear iterations = R2  "
		, (integer)80, (integer)56);
	_daskr_xerrwd_(msg, &c__56, &c__501, &c__0, &c__0, &c__0, &c__0, &c__2, &tn,
		&avlin, (integer)80);
    }
    if (lcfn) {
	_daskr_str_copy(msg, "DASKR-- Warning. Poor iterative algorithm performance   "
		, (integer)80, (integer)56);
	_daskr_xerrwd_(msg, &c__56, &c__502, &c__0, &c__0, &c__0, &c__0, &c__0, &
		c_b38, &c_b38, (integer)80);
	_daskr_str_copy(msg, "      at T = R1. Nonlinear convergence failure rate = R2"
		, (integer)80, (integer)56);
	_daskr_xerrwd_(msg, &c__56, &c__502, &c__0, &c__0, &c__0, &c__0, &c__2, &tn,
		&rcfn, (integer)80);
    }
    if (lcfl) {
	_daskr_str_copy(msg, "DASKR-- Warning. Poor iterative algorithm performance   "
		, (integer)80, (integer)56);
	_daskr_xerrwd_(msg, &c__56, &c__503, &c__0, &c__0, &c__0, &c__0, &c__0, &
		c_b38, &c_b38, (integer)80);
	_daskr_str_copy(msg, "      at T = R1. Linear convergence failure rate = R2   "
		, (integer)80, (integer)56);
	_daskr_xerrwd_(msg, &c__56, &c__503, &c__0, &c__0, &c__0, &c__0, &c__2, &tn,
		&rcfl, (integer)80);
    }

/*     Update WT and VT, if this is not the first call. */

L510:
    _daskr_ddawts_(neq, &info[2], &rtol[1], &atol[1], &rwork[lphi], &rwork[lwt], &
	    rpar[1], &ipar[1]);
    _daskr_dinvwt_(neq, &rwork[lwt], &ier);
    if (ier != 0) {
	*idid = -3;
	goto L527;
    }
    if (info[16] != 0) {
	i__1 = *neq;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L515: */
/* Computing MAX */
	    i__2 = iwork[lid + i__ - 1];
	    rwork[lvt + i__ - 1] = MAX(i__2,0) * rwork[lwt + i__ - 1];
	}
    }

/*     Test for too much accuracy requested. */

    r__ = _daskr_ddwnrm_(neq, &rwork[lphi], &rwork[lwt], &rpar[1], &ipar[1]) * 100. *
	     uround;
    if (r__ <= 1.) {
	goto L525;
    }

/*     Multiply RTOL and ATOL by R and return. */

    if (info[2] == 1) {
	goto L523;
    }
    rtol[1] = r__ * rtol[1];
    atol[1] = r__ * atol[1];
    *idid = -2;
    goto L527;
L523:
    i__2 = *neq;
    for (i__ = 1; i__ <= i__2; ++i__) {
	rtol[i__] = r__ * rtol[i__];
/* L524: */
	atol[i__] = r__ * atol[i__];
    }
    *idid = -2;
    goto L527;
L525:

/*     Compute minimum stepsize. */

/* Computing MAX */
    d__1 = fabs(tn), d__2 = fabs(*tout);
    hmin = uround * 4. * MAX(d__1,d__2);

/*     Test H vs. HMAX */
    if (info[7] != 0) {
	rh = fabs(h__) / rwork[2];
	if (rh > 1.) {
	    h__ /= rh;
	}
    }

/*     Call the one-step integrator. */
/*     Note that INFO(12) represents the nonlinear solver type. */
/*     Pass the required nonlinear solver, depending upon INFO(12). */

    if (info[12] == 0) {
	_daskr_ddstp_(&tn, &y[1], &yprime[1], neq, (Unknown_fp)res, (Unknown_fp)jac, (Unknown_fp)psol,
		&h__, &rwork[lwt], &rwork[lvt], &info[1], idid, &rpar[1], &
		ipar[1], &rwork[lphi], &rwork[lsavr], &rwork[61], &rwork[le], 
		&rwork[lwm], &iwork[1], &rwork[21], &rwork[27], &rwork[33], &
		rwork[39], &rwork[45], &rwork[5], &rwork[6], &rwork[7], &
		rwork[8], &hmin, &rwork[9], &rwork[10], &rwork[11], &rwork[12]
		, &rwork[13], &iwork[6], &iwork[5], &info[15], &iwork[7], &
		iwork[8], &iwork[9], &nonneg, &info[12], (Unknown_fp)_daskr_dnedd_);
    } else if (info[12] == 1) {
	_daskr_ddstp_(&tn, &y[1], &yprime[1], neq, (Unknown_fp)res, (Unknown_fp)jac, (Unknown_fp)psol,
		&h__, &rwork[lwt], &rwork[lvt], &info[1], idid, &rpar[1], &
		ipar[1], &rwork[lphi], &rwork[lsavr], &rwork[61], &rwork[le], 
		&rwork[lwm], &iwork[1], &rwork[21], &rwork[27], &rwork[33], &
		rwork[39], &rwork[45], &rwork[5], &rwork[6], &rwork[7], &
		rwork[8], &hmin, &rwork[9], &rwork[10], &rwork[11], &rwork[12]
		, &rwork[13], &iwork[6], &iwork[5], &info[15], &iwork[7], &
		iwork[8], &iwork[9], &nonneg, &info[12], (Unknown_fp)_daskr_dnedk_);
    }

L527:
    if (*idid < 0) {
	goto L600;
    }

/* ----------------------------------------------------------------------- */
/*     This block handles the case of a successful return from DDSTP */
/*     (IDID=1).  Test for stop conditions. */
/* ----------------------------------------------------------------------- */

    if (*nrt == 0) {
	goto L530;
    }

/*     Check for a zero of R near TN. */

    _daskr_drchek_(&c__3, (Unknown_fp)rt, nrt, neq, &tn, tout, &y[1], &yprime[1], &rwork[
	    lphi], &rwork[39], &iwork[8], &rwork[lr0], &rwork[lr1], &rwork[
	    lrx], jroot, &irt, &rwork[9], &info[3], &rwork[1], &iwork[1], &
	    rpar[1], &ipar[1]);
    if (irt != 1) {
	goto L530;
    }
    iwork[37] = 1;
    *idid = 5;
    *t = rwork[51];
    goto L580;

L530:
    if (info[4] == 0) {
/*        Stopping tests for the case of no TSTOP. ---------------------- */
	if ((tn - *tout) * h__ >= 0.) {
	    _daskr_ddatrp_(&tn, tout, &y[1], &yprime[1], neq, &iwork[8], &rwork[lphi]
		    , &rwork[39]);
	    *t = *tout;
	    *idid = 3;
	    goto L580;
	}
	if (info[3] == 0) {
	    goto L500;
	}
	*t = tn;
	*idid = 1;
	goto L580;
    }

/* L540: */
    if (info[3] != 0) {
	goto L550;
    }
/*     Stopping tests for the TSTOP case, interval-output mode. --------- */
    if ((d__1 = tn - tstop, fabs(d__1)) <= uround * 100. * (fabs(tn) + fabs(h__))
	    ) {
	_daskr_ddatrp_(&tn, &tstop, &y[1], &yprime[1], neq, &iwork[8], &rwork[lphi],
		&rwork[39]);
	*t = tstop;
	*idid = 2;
	goto L580;
    }
    if ((tn - *tout) * h__ >= 0.) {
	_daskr_ddatrp_(&tn, tout, &y[1], &yprime[1], neq, &iwork[8], &rwork[lphi], &
		rwork[39]);
	*t = *tout;
	*idid = 3;
	goto L580;
    }
    tnext = tn + h__;
    if ((tnext - tstop) * h__ <= 0.) {
	goto L500;
    }
    h__ = tstop - tn;
    goto L500;

L550:
/*     Stopping tests for the TSTOP case, intermediate-output mode. ----- */
    if ((d__1 = tn - tstop, fabs(d__1)) <= uround * 100. * (fabs(tn) + fabs(h__))
	    ) {
	_daskr_ddatrp_(&tn, &tstop, &y[1], &yprime[1], neq, &iwork[8], &rwork[lphi],
		&rwork[39]);
	*t = tstop;
	*idid = 2;
	goto L580;
    }
    if ((tn - *tout) * h__ >= 0.) {
	_daskr_ddatrp_(&tn, tout, &y[1], &yprime[1], neq, &iwork[8], &rwork[lphi], &
		rwork[39]);
	*t = *tout;
	*idid = 3;
	goto L580;
    }
    *t = tn;
    *idid = 1;

L580:

/* ----------------------------------------------------------------------- */
/*     All successful returns from DDASKR are made from this block. */
/* ----------------------------------------------------------------------- */

L590:
    rwork[4] = tn;
    rwork[52] = *t;
    rwork[3] = h__;
    return 0;

/* ----------------------------------------------------------------------- */
/*     This block handles all unsuccessful returns other than for */
/*     illegal input. */
/* ----------------------------------------------------------------------- */

L600:
    itemp = -(*idid);
    switch (itemp) {
	case 1:  goto L610;
	case 2:  goto L620;
	case 3:  goto L630;
	case 4:  goto L700;
	case 5:  goto L655;
	case 6:  goto L640;
	case 7:  goto L650;
	case 8:  goto L660;
	case 9:  goto L670;
	case 10:  goto L675;
	case 11:  goto L680;
	case 12:  goto L685;
	case 13:  goto L690;
	case 14:  goto L695;
    }

/*     The maximum number of steps was taken before */
/*     reaching tout. */

L610:
    _daskr_str_copy(msg, "DASKR--  AT CURRENT T (=R1)  500 STEPS", (integer)80, (integer)
	    38);
    _daskr_xerrwd_(msg, &c__38, &c__610, &c__0, &c__0, &c__0, &c__0, &c__1, &tn, &
	    c_b38, (integer)80);
    _daskr_str_copy(msg, "DASKR--  TAKEN ON THIS CALL BEFORE REACHING TOUT", (integer)
	    80, (integer)48);
    _daskr_xerrwd_(msg, &c__48, &c__611, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    goto L700;

/*     Too much accuracy for machine precision. */

L620:
    _daskr_str_copy(msg, "DASKR--  AT T (=R1) TOO MUCH ACCURACY REQUESTED", (integer)80,
	     (integer)47);
    _daskr_xerrwd_(msg, &c__47, &c__620, &c__0, &c__0, &c__0, &c__0, &c__1, &tn, &
	    c_b38, (integer)80);
    _daskr_str_copy(msg, "DASKR--  FOR PRECISION OF MACHINE. RTOL AND ATOL", (integer)
	    80, (integer)48);
    _daskr_xerrwd_(msg, &c__48, &c__621, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    _daskr_str_copy(msg, "DASKR--  WERE INCREASED BY A FACTOR R (=R1)", (integer)80, (
	    integer)43);
    _daskr_xerrwd_(msg, &c__43, &c__622, &c__0, &c__0, &c__0, &c__0, &c__1, &r__, &
	    c_b38, (integer)80);
    goto L700;

/*     WT(I) .LE. 0.0D0 for some I (not at start of problem). */

L630:
    _daskr_str_copy(msg, "DASKR--  AT T (=R1) SOME ELEMENT OF WT", (integer)80, (integer)
	    38);
    _daskr_xerrwd_(msg, &c__38, &c__630, &c__0, &c__0, &c__0, &c__0, &c__1, &tn, &
	    c_b38, (integer)80);
    _daskr_str_copy(msg, "DASKR--  HAS BECOME .LE. 0.0", (integer)80, (integer)28);
    _daskr_xerrwd_(msg, &c__28, &c__631, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    goto L700;

/*     Error test failed repeatedly or with H=HMIN. */

L640:
    _daskr_str_copy(msg, "DASKR--  AT T (=R1) AND STEPSIZE H (=R2) THE", (integer)80, (
	    integer)44);
    _daskr_xerrwd_(msg, &c__44, &c__640, &c__0, &c__0, &c__0, &c__0, &c__2, &tn, &
	    h__, (integer)80);
    _daskr_str_copy(msg, "DASKR--  ERROR TEST FAILED REPEATEDLY OR WITH ABS(H)=HMIN", (
	    integer)80, (integer)57);
    _daskr_xerrwd_(msg, &c__57, &c__641, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    goto L700;

/*     Nonlinear solver failed to converge repeatedly or with H=HMIN. */

L650:
    _daskr_str_copy(msg, "DASKR--  AT T (=R1) AND STEPSIZE H (=R2) THE", (integer)80, (
	    integer)44);
    _daskr_xerrwd_(msg, &c__44, &c__650, &c__0, &c__0, &c__0, &c__0, &c__2, &tn, &
	    h__, (integer)80);
    _daskr_str_copy(msg, "DASKR--  NONLINEAR SOLVER FAILED TO CONVERGE", (integer)80, (
	    integer)44);
    _daskr_xerrwd_(msg, &c__44, &c__651, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    _daskr_str_copy(msg, "DASKR--  REPEATEDLY OR WITH ABS(H)=HMIN", (integer)80, (
	    integer)39);
    _daskr_xerrwd_(msg, &c__40, &c__652, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    goto L700;

/*     The preconditioner had repeated failures. */

L655:
    _daskr_str_copy(msg, "DASKR--  AT T (=R1) AND STEPSIZE H (=R2) THE", (integer)80, (
	    integer)44);
    _daskr_xerrwd_(msg, &c__44, &c__655, &c__0, &c__0, &c__0, &c__0, &c__2, &tn, &
	    h__, (integer)80);
    _daskr_str_copy(msg, "DASKR--  PRECONDITIONER HAD REPEATED FAILURES.", (integer)80,
	    (integer)46);
    _daskr_xerrwd_(msg, &c__46, &c__656, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    goto L700;

/*     The iteration matrix is singular. */

L660:
    _daskr_str_copy(msg, "DASKR--  AT T (=R1) AND STEPSIZE H (=R2) THE", (integer)80, (
	    integer)44);
    _daskr_xerrwd_(msg, &c__44, &c__660, &c__0, &c__0, &c__0, &c__0, &c__2, &tn, &
	    h__, (integer)80);
    _daskr_str_copy(msg, "DASKR--  ITERATION MATRIX IS SINGULAR.", (integer)80, (integer)
	    38);
    _daskr_xerrwd_(msg, &c__38, &c__661, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    goto L700;

/*     Nonlinear system failure preceded by error test failures. */

L670:
    _daskr_str_copy(msg, "DASKR--  AT T (=R1) AND STEPSIZE H (=R2) THE", (integer)80, (
	    integer)44);
    _daskr_xerrwd_(msg, &c__44, &c__670, &c__0, &c__0, &c__0, &c__0, &c__2, &tn, &
	    h__, (integer)80);
    _daskr_str_copy(msg, "DASKR--  NONLINEAR SOLVER COULD NOT CONVERGE.", (integer)80, (
	    integer)45);
    _daskr_xerrwd_(msg, &c__45, &c__671, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    _daskr_str_copy(msg, "DASKR--  ALSO, THE ERROR TEST FAILED REPEATEDLY.", (integer)
	    80, (integer)48);
    _daskr_xerrwd_(msg, &c__49, &c__672, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    goto L700;

/*     Nonlinear system failure because IRES = -1. */

L675:
    _daskr_str_copy(msg, "DASKR--  AT T (=R1) AND STEPSIZE H (=R2) THE", (integer)80, (
	    integer)44);
    _daskr_xerrwd_(msg, &c__44, &c__675, &c__0, &c__0, &c__0, &c__0, &c__2, &tn, &
	    h__, (integer)80);
    _daskr_str_copy(msg, "DASKR--  NONLINEAR SYSTEM SOLVER COULD NOT CONVERGE", (
	    integer)80, (integer)51);
    _daskr_xerrwd_(msg, &c__51, &c__676, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    _daskr_str_copy(msg, "DASKR--  BECAUSE IRES WAS EQUAL TO MINUS ONE", (integer)80, (
	    integer)44);
    _daskr_xerrwd_(msg, &c__44, &c__677, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    goto L700;

/*     Failure because IRES = -2. */

L680:
    _daskr_str_copy(msg, "DASKR--  AT T (=R1) AND STEPSIZE H (=R2)", (integer)80, (
	    integer)40);
    _daskr_xerrwd_(msg, &c__40, &c__680, &c__0, &c__0, &c__0, &c__0, &c__2, &tn, &
	    h__, (integer)80);
    _daskr_str_copy(msg, "DASKR--  IRES WAS EQUAL TO MINUS TWO", (integer)80, (integer)
	    36);
    _daskr_xerrwd_(msg, &c__36, &c__681, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    goto L700;

/*     Failed to compute initial YPRIME. */

L685:
    _daskr_str_copy(msg, "DASKR--  AT T (=R1) AND STEPSIZE H (=R2) THE", (integer)80, (
	    integer)44);
    _daskr_xerrwd_(msg, &c__44, &c__685, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    _daskr_str_copy(msg, "DASKR--  INITIAL (Y,YPRIME) COULD NOT BE COMPUTED", (integer)
	    80, (integer)49);
    _daskr_xerrwd_(msg, &c__49, &c__686, &c__0, &c__0, &c__0, &c__0, &c__2, &tn, &h0,
	     (integer)80);
    goto L700;

/*     Failure because IER was negative from PSOL. */

L690:
    _daskr_str_copy(msg, "DASKR--  AT T (=R1) AND STEPSIZE H (=R2)", (integer)80, (
	    integer)40);
    _daskr_xerrwd_(msg, &c__40, &c__690, &c__0, &c__0, &c__0, &c__0, &c__2, &tn, &
	    h__, (integer)80);
    _daskr_str_copy(msg, "DASKR--  IER WAS NEGATIVE FROM PSOL", (integer)80, (integer)35)
	    ;
    _daskr_xerrwd_(msg, &c__35, &c__691, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    goto L700;

/*     Failure because the linear system solver could not converge. */

L695:
    _daskr_str_copy(msg, "DASKR--  AT T (=R1) AND STEPSIZE H (=R2) THE", (integer)80, (
	    integer)44);
    _daskr_xerrwd_(msg, &c__44, &c__695, &c__0, &c__0, &c__0, &c__0, &c__2, &tn, &
	    h__, (integer)80);
    _daskr_str_copy(msg, "DASKR--  LINEAR SYSTEM SOLVER COULD NOT CONVERGE.", (integer)
	    80, (integer)49);
    _daskr_xerrwd_(msg, &c__50, &c__696, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    goto L700;


L700:
    info[1] = -1;
    *t = tn;
    rwork[4] = tn;
    rwork[3] = h__;
    return 0;

/* ----------------------------------------------------------------------- */
/*     This block handles all error returns due to illegal input, */
/*     as detected before calling DDSTP. */
/*     First the error message routine is called.  If this happens */
/*     twice in succession, execution is terminated. */
/* ----------------------------------------------------------------------- */

L701:
    _daskr_str_copy(msg, "DASKR--  ELEMENT (=I1) OF INFO VECTOR IS NOT VALID", (integer)
	    80, (integer)50);
    _daskr_xerrwd_(msg, &c__50, &c__1, &c__0, &c__1, &itemp, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    goto L750;
L702:
    _daskr_str_copy(msg, "DASKR--  NEQ (=I1) .LE. 0", (integer)80, (integer)25);
    _daskr_xerrwd_(msg, &c__25, &c__2, &c__0, &c__1, neq, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    goto L750;
L703:
    _daskr_str_copy(msg, "DASKR--  MAXORD (=I1) NOT IN RANGE", (integer)80, (integer)34);
    _daskr_xerrwd_(msg, &c__34, &c__3, &c__0, &c__1, &mxord, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    goto L750;
L704:
    _daskr_str_copy(msg, "DASKR--  RWORK LENGTH NEEDED, LENRW (=I1), EXCEEDS LRW (=I2)"
	    , (integer)80, (integer)60);
    _daskr_xerrwd_(msg, &c__60, &c__4, &c__0, &c__2, &lenrw, lrw, &c__0, &c_b38, &
	    c_b38, (integer)80);
    goto L750;
L705:
    _daskr_str_copy(msg, "DASKR--  IWORK LENGTH NEEDED, LENIW (=I1), EXCEEDS LIW (=I2)"
	    , (integer)80, (integer)60);
    _daskr_xerrwd_(msg, &c__60, &c__5, &c__0, &c__2, &leniw, liw, &c__0, &c_b38, &
	    c_b38, (integer)80);
    goto L750;
L706:
    _daskr_str_copy(msg, "DASKR--  SOME ELEMENT OF RTOL IS .LT. 0", (integer)80, (
	    integer)39);
    _daskr_xerrwd_(msg, &c__39, &c__6, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    goto L750;
L707:
    _daskr_str_copy(msg, "DASKR--  SOME ELEMENT OF ATOL IS .LT. 0", (integer)80, (
	    integer)39);
    _daskr_xerrwd_(msg, &c__39, &c__7, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    goto L750;
L708:
    _daskr_str_copy(msg, "DASKR--  ALL ELEMENTS OF RTOL AND ATOL ARE ZERO", (integer)80,
	     (integer)47);
    _daskr_xerrwd_(msg, &c__47, &c__8, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    goto L750;
L709:
    _daskr_str_copy(msg, "DASKR--  INFO(4) = 1 AND TSTOP (=R1) BEHIND TOUT (=R2)", (
	    integer)80, (integer)54);
    _daskr_xerrwd_(msg, &c__54, &c__9, &c__0, &c__0, &c__0, &c__0, &c__2, &tstop,
	    tout, (integer)80);
    goto L750;
L710:
    _daskr_str_copy(msg, "DASKR--  HMAX (=R1) .LT. 0.0", (integer)80, (integer)28);
    _daskr_xerrwd_(msg, &c__28, &c__10, &c__0, &c__0, &c__0, &c__0, &c__1, &hmax, &
	    c_b38, (integer)80);
    goto L750;
L711:
    _daskr_str_copy(msg, "DASKR--  TOUT (=R1) BEHIND T (=R2)", (integer)80, (integer)34);
    _daskr_xerrwd_(msg, &c__34, &c__11, &c__0, &c__0, &c__0, &c__0, &c__2, tout, t, (
	    integer)80);
    goto L750;
L712:
    _daskr_str_copy(msg, "DASKR--  INFO(8)=1 AND H0=0.0", (integer)80, (integer)29);
    _daskr_xerrwd_(msg, &c__29, &c__12, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    goto L750;
L713:
    _daskr_str_copy(msg, "DASKR--  SOME ELEMENT OF WT IS .LE. 0.0", (integer)80, (
	    integer)39);
    _daskr_xerrwd_(msg, &c__39, &c__13, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    goto L750;
L714:
    _daskr_str_copy(msg, "DASKR-- TOUT (=R1) TOO CLOSE TO T (=R2) TO START INTEGRATION"
	    , (integer)80, (integer)60);
    _daskr_xerrwd_(msg, &c__60, &c__14, &c__0, &c__0, &c__0, &c__0, &c__2, tout, t, (
	    integer)80);
    goto L750;
L715:
    _daskr_str_copy(msg, "DASKR--  INFO(4)=1 AND TSTOP (=R1) BEHIND T (=R2)", (integer)
	    80, (integer)49);
    _daskr_xerrwd_(msg, &c__49, &c__15, &c__0, &c__0, &c__0, &c__0, &c__2, &tstop, t,
	     (integer)80);
    goto L750;
L717:
    _daskr_str_copy(msg, "DASKR--  ML (=I1) ILLEGAL. EITHER .LT. 0 OR .GT. NEQ", (
	    integer)80, (integer)52);
    _daskr_xerrwd_(msg, &c__52, &c__17, &c__0, &c__1, &iwork[1], &c__0, &c__0, &
	    c_b38, &c_b38, (integer)80);
    goto L750;
L718:
    _daskr_str_copy(msg, "DASKR--  MU (=I1) ILLEGAL. EITHER .LT. 0 OR .GT. NEQ", (
	    integer)80, (integer)52);
    _daskr_xerrwd_(msg, &c__52, &c__18, &c__0, &c__1, &iwork[2], &c__0, &c__0, &
	    c_b38, &c_b38, (integer)80);
    goto L750;
L719:
    _daskr_str_copy(msg, "DASKR--  TOUT (=R1) IS EQUAL TO T (=R2)", (integer)80, (
	    integer)39);
    _daskr_xerrwd_(msg, &c__39, &c__19, &c__0, &c__0, &c__0, &c__0, &c__2, tout, t, (
	    integer)80);
    goto L750;
L720:
    _daskr_str_copy(msg, "DASKR--  MAXL (=I1) ILLEGAL. EITHER .LT. 1 OR .GT. NEQ", (
	    integer)80, (integer)54);
    _daskr_xerrwd_(msg, &c__54, &c__20, &c__0, &c__1, &iwork[24], &c__0, &c__0, &
	    c_b38, &c_b38, (integer)80);
    goto L750;
L721:
    _daskr_str_copy(msg, "DASKR--  KMP (=I1) ILLEGAL. EITHER .LT. 1 OR .GT. MAXL", (
	    integer)80, (integer)54);
    _daskr_xerrwd_(msg, &c__54, &c__21, &c__0, &c__1, &iwork[25], &c__0, &c__0, &
	    c_b38, &c_b38, (integer)80);
    goto L750;
L722:
    _daskr_str_copy(msg, "DASKR--  NRMAX (=I1) ILLEGAL. .LT. 0", (integer)80, (integer)
	    36);
    _daskr_xerrwd_(msg, &c__36, &c__22, &c__0, &c__1, &iwork[26], &c__0, &c__0, &
	    c_b38, &c_b38, (integer)80);
    goto L750;
L723:
    _daskr_str_copy(msg, "DASKR--  EPLI (=R1) ILLEGAL. EITHER .LE. 0.D0 OR .GE. 1.D0",
	    (integer)80, (integer)58);
    _daskr_xerrwd_(msg, &c__58, &c__23, &c__0, &c__0, &c__0, &c__0, &c__1, &rwork[10]
	    , &c_b38, (integer)80);
    goto L750;
L724:
    _daskr_str_copy(msg, "DASKR--  ILLEGAL IWORK VALUE FOR INFO(11) .NE. 0", (integer)
	    80, (integer)48);
    _daskr_xerrwd_(msg, &c__48, &c__24, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    goto L750;
L725:
    _daskr_str_copy(msg, "DASKR--  ONE OF THE INPUTS FOR INFO(17) = 1 IS ILLEGAL", (
	    integer)80, (integer)54);
    _daskr_xerrwd_(msg, &c__54, &c__25, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    goto L750;
L726:
    _daskr_str_copy(msg, "DASKR--  ILLEGAL IWORK VALUE FOR INFO(10) .NE. 0", (integer)
	    80, (integer)48);
    _daskr_xerrwd_(msg, &c__48, &c__26, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    goto L750;
L727:
    _daskr_str_copy(msg, "DASKR--  Y(I) AND IWORK(40+I) (I=I1) INCONSISTENT", (integer)
	    80, (integer)49);
    _daskr_xerrwd_(msg, &c__49, &c__27, &c__0, &c__1, &iret, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    goto L750;
L730:
    _daskr_str_copy(msg, "DASKR--  NRT (=I1) .LT. 0", (integer)80, (integer)25);
    _daskr_xerrwd_(msg, &c__25, &c__30, &c__1, &c__1, nrt, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    goto L750;
L731:
    _daskr_str_copy(msg, "DASKR--  R IS ILL-DEFINED.  ZERO VALUES WERE FOUND AT TWO", (
	    integer)80, (integer)57);
    _daskr_xerrwd_(msg, &c__57, &c__31, &c__1, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    _daskr_str_copy(msg, "         VERY CLOSE T VALUES, AT T = R1", (integer)80, (
	    integer)39);
    _daskr_xerrwd_(msg, &c__39, &c__31, &c__1, &c__0, &c__0, &c__0, &c__1, &rwork[51]
	    , &c_b38, (integer)80);

L750:
    if (info[1] == -1) {
	goto L760;
    }
    info[1] = -1;
    *idid = -33;
    return 0;
L760:
    _daskr_str_copy(msg, "DASKR--  REPEATED OCCURRENCES OF ILLEGAL INPUT", (integer)80,
	    (integer)46);
    _daskr_xerrwd_(msg, &c__46, &c__701, &c__0, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
/* L770: */
    _daskr_str_copy(msg, "DASKR--  RUN TERMINATED. APPARENT INFINITE LOOP", (integer)80,
	     (integer)47);
    _daskr_xerrwd_(msg, &c__47, &c__702, &c__1, &c__0, &c__0, &c__0, &c__0, &c_b38, &
	    c_b38, (integer)80);
    return 0;

/* ------END OF SUBROUTINE DDASKR----------------------------------------- */
} /* ddaskr_ */

/* Subroutine */ int _daskr_drchek_(integer *job, Unknown_fp rt, integer *nrt, integer *
	neq, real_number *tn, real_number *tout, real_number *y, real_number *yp,
	real_number *phi, real_number *psi, integer *kold, real_number *r0,
	real_number *r1, real_number *rx, integer *jroot, integer *irt,
	real_number *uround, integer *info3, real_number *rwork, integer *iwork,
	 real_number *rpar, integer *ipar)
{
    /* Initialized data */

    static real_number zero = 0.;

    /* System generated locals */
    integer phi_dim1, phi_offset, i__1;
    real_number d__1;

    /* Builtin functions */
    double _daskr_real_sign(real_number *, real_number *);

    /* Local variables */
    static real_number h__;
    static integer i__;
    static real_number x, t1, temp1, temp2;
    static integer jflag;
    static real_number hminr;
    extern /* Subroutine */ int _daskr_dcopy_(integer *, real_number *, integer *,
	    real_number *, integer *);
    static integer zroot;
    extern /* Subroutine */ int _daskr_ddatrp_(real_number *, real_number *,
	    real_number *, real_number *, integer *, integer *, real_number *,
	    real_number *), _daskr_droots_(integer *, real_number *, integer *,
	    real_number *, real_number *, real_number *, real_number *,
	    real_number *, real_number *, integer *);


/* ***BEGIN PROLOGUE  DRCHEK */
/* ***REFER TO DDASKR */
/* ***ROUTINES CALLED  DDATRP, DROOTS, DCOPY, RT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   020815  DATE WRITTEN */
/*   021217  Added test for roots close when JOB = 2. */
/*   050510  Changed T increment after 110 so that TEMP1/H .ge. 0.1. */
/*   071003  Fixed bug in TEMP2 (HMINR) below 110. */
/*   110608  Fixed bug in setting of T1 at 300. */
/* ***END PROLOGUE  DRCHEK */

/* Pointers into IWORK: */
/* Pointers into RWORK: */
    /* Parameter adjustments */
    phi_dim1 = *neq;
    phi_offset = 1 + phi_dim1;
    phi -= phi_offset;
    --y;
    --yp;
    --psi;
    --r0;
    --r1;
    --rx;
    --jroot;
    --rwork;
    --iwork;

    /* Function Body */
/* ----------------------------------------------------------------------- */
/* This routine checks for the presence of a root of R(T,Y,Y') in the */
/* vicinity of the current T, in a manner depending on the */
/* input flag JOB.  It calls subroutine DROOTS to locate the root */
/* as precisely as possible. */

/* In addition to variables described previously, DRCHEK */
/* uses the following for communication.. */
/* JOB    = integer flag indicating type of call.. */
/*          JOB = 1 means the problem is being initialized, and DRCHEK */
/*                  is to look for a root at or very near the initial T. */
/*          JOB = 2 means a continuation call to the solver was just */
/*                  made, and DRCHEK is to check for a root in the */
/*                  relevant part of the step last taken. */
/*          JOB = 3 means a successful step was just taken, and DRCHEK */
/*                  is to look for a root in the interval of the step. */
/* R0     = array of length NRT, containing the value of R at T = T0. */
/*          R0 is input for JOB .ge. 2 and on output in all cases. */
/* R1,RX  = arrays of length NRT for work space. */
/* IRT    = completion flag.. */
/*          IRT = 0  means no root was found. */
/*          IRT = -1 means JOB = 1 and a zero was found both at T0 and */
/*                   and very close to T0. */
/*          IRT = -2 means JOB = 2 and some Ri was found to have a zero */
/*                   both at T0 and very close to T0. */
/*          IRT = 1  means a legitimate root was found (JOB = 2 or 3). */
/*                   On return, T0 is the root location, and Y is the */
/*                   corresponding solution vector. */
/* T0     = value of T at one endpoint of interval of interest.  Only */
/*          roots beyond T0 in the direction of integration are sought. */
/*          T0 is input if JOB .ge. 2, and output in all cases. */
/*          T0 is updated by DRCHEK, whether a root is found or not. */
/*          Stored in the global array RWORK. */
/* TLAST  = last value of T returned by the solver (input only). */
/*          Stored in the global array RWORK. */
/* TOUT   = final output time for the solver. */
/* IRFND  = input flag showing whether the last step taken had a root. */
/*          IRFND = 1 if it did, = 0 if not. */
/*          Stored in the global array IWORK. */
/* INFO3  = copy of INFO(3) (input only). */
/* ----------------------------------------------------------------------- */

    h__ = psi[1];
    *irt = 0;
    i__1 = *nrt;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	jroot[i__] = 0;
    }
    hminr = (fabs(*tn) + fabs(h__)) * *uround * 100.;

    switch (*job) {
	case 1:  goto L100;
	case 2:  goto L200;
	case 3:  goto L300;
    }

/* Evaluate R at initial T (= RWORK(LT0)); check for zero values.-------- */
L100:
    _daskr_ddatrp_(tn, &rwork[51], &y[1], &yp[1], neq, kold, &phi[phi_offset], &psi[
	    1]);
    ((int (*)(integer*, real_number *, real_number*, real_number*, integer *, real_number *, real_number *, integer *))*rt)
      (neq, &rwork[51], &y[1], &yp[1], nrt, &r0[1], rpar, ipar);
    iwork[36] = 1;
    zroot = _FALSE_;
    i__1 = *nrt;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L110: */
	if ((d__1 = r0[i__], fabs(d__1)) == zero) {
	    zroot = _TRUE_;
	}
    }
    if (! zroot) {
	goto L190;
    }
/* R has a zero at T.  Look at R at T + (small increment). -------------- */
/* Computing MAX */
    d__1 = hminr / fabs(h__);
    temp2 = MAX(d__1,.1);
    temp1 = temp2 * h__;
    rwork[51] += temp1;
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L120: */
	y[i__] += temp2 * phi[i__ + (phi_dim1 << 1)];
    }
    ((int (*)(integer*, real_number *, real_number*, real_number*, integer *, real_number *, real_number *, integer *))*rt)
      (neq, &rwork[51], &y[1], &yp[1], nrt, &r0[1], rpar, ipar);
    ++iwork[36];
    zroot = _FALSE_;
    i__1 = *nrt;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L130: */
	if ((d__1 = r0[i__], fabs(d__1)) == zero) {
	    zroot = _TRUE_;
	}
    }
    if (! zroot) {
	goto L190;
    }
/* R has a zero at T and also close to T.  Take error return. ----------- */
    *irt = -1;
    return 0;

L190:
    return 0;

L200:
    if (iwork[37] == 0) {
	goto L260;
    }
/* If a root was found on the previous step, evaluate R0 = R(T0). ------- */
    _daskr_ddatrp_(tn, &rwork[51], &y[1], &yp[1], neq, kold, &phi[phi_offset], &psi[
	    1]);
    ((int (*)(integer*, real_number *, real_number*, real_number*, integer *, real_number *, real_number *, integer *))*rt)
      (neq, &rwork[51], &y[1], &yp[1], nrt, &r0[1], rpar, ipar);
    ++iwork[36];
    zroot = _FALSE_;
    i__1 = *nrt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if ((d__1 = r0[i__], fabs(d__1)) == zero) {
	    zroot = _TRUE_;
	    jroot[i__] = 1;
	}
/* L210: */
    }
    if (! zroot) {
	goto L260;
    }
/* R has a zero at T0.  Look at R at T0+ = T0 + (small increment). ------ */
    temp1 = _daskr_real_sign(&hminr, &h__);
    rwork[51] += temp1;
    if ((rwork[51] - *tn) * h__ < zero) {
	goto L230;
    }
    temp2 = temp1 / h__;
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L220: */
	y[i__] += temp2 * phi[i__ + (phi_dim1 << 1)];
    }
    goto L240;
L230:
    _daskr_ddatrp_(tn, &rwork[51], &y[1], &yp[1], neq, kold, &phi[phi_offset], &psi[
	    1]);
L240:
    ((int (*)(integer*, real_number *, real_number*, real_number*, integer *, real_number *, real_number *, integer *))*rt)
      (neq, &rwork[51], &y[1], &yp[1], nrt, &r0[1], rpar, ipar);
    ++iwork[36];
    i__1 = *nrt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if ((d__1 = r0[i__], fabs(d__1)) > zero) {
	    goto L250;
	}
/* If Ri has a zero at both T0+ and T0, return an error flag. ----------- */
	if (jroot[i__] == 1) {
	    *irt = -2;
	    return 0;
	} else {
/* If Ri has a zero at T0+, but not at T0, return valid root. ----------- */
	    jroot[i__] = (integer) (-_daskr_real_sign(&c_b758, &r0[i__]));
	    *irt = 1;
	}
L250:
	;
    }
    if (*irt == 1) {
	return 0;
    }
/* R0 has no zero components.  Proceed to check relevant interval. ------ */
L260:
    if (*tn == rwork[52]) {
	return 0;
    }

L300:
/* Set T1 to TN or TOUT, whichever comes first, and get R at T1. -------- */
    if ((*tout - *tn) * h__ >= zero) {
	t1 = *tn;
	goto L330;
    }
    t1 = *tout;
    if ((t1 - rwork[51]) * h__ <= zero) {
	goto L390;
    }
L330:
    _daskr_ddatrp_(tn, &t1, &y[1], &yp[1], neq, kold, &phi[phi_offset], &psi[1]);
    ((int (*)(integer*, real_number *, real_number*, real_number*, integer *, real_number *, real_number *, integer *))*rt)
      (neq, &t1, &y[1], &yp[1], nrt, &r1[1], rpar, ipar);
    ++iwork[36];
/* Call DROOTS to search for root in interval from T0 to T1. ------------ */
    jflag = 0;
L350:
    _daskr_droots_(nrt, &hminr, &jflag, &rwork[51], &t1, &r0[1], &r1[1], &rx[1], &x,
	    &jroot[1]);
    if (jflag > 1) {
	goto L360;
    }
    _daskr_ddatrp_(tn, &x, &y[1], &yp[1], neq, kold, &phi[phi_offset], &psi[1]);
    ((int (*)(integer*, real_number *, real_number*, real_number*, integer *, real_number *, real_number *, integer *))*rt)
      (neq, &x, &y[1], &yp[1], nrt, &rx[1], rpar, ipar);
    ++iwork[36];
    goto L350;
L360:
    rwork[51] = x;
    _daskr_dcopy_(nrt, &rx[1], &c__1, &r0[1], &c__1);
    if (jflag == 4) {
	goto L390;
    }
/* Found a root.  Interpolate to X and return. -------------------------- */
    _daskr_ddatrp_(tn, &x, &y[1], &yp[1], neq, kold, &phi[phi_offset], &psi[1]);
    *irt = 1;
    return 0;

L390:
    return 0;
/* ---------------------- END OF SUBROUTINE DRCHEK ----------------------- */
} /* drchek_ */

/* Subroutine */ int _daskr_droots_(integer *nrt, real_number *hmin, integer *jflag,
	real_number *x0, real_number *x1, real_number *r0, real_number *r1,
	real_number *rx, real_number *x, integer *jroot)
{
    /* Initialized data */

    static real_number zero = 0.;
    static real_number tenth = .1;
    static real_number half = .5;
    static real_number five = 5.;

    /* System generated locals */
    integer i__1;
    real_number d__1;

    /* Builtin functions */
    double _daskr_real_sign(real_number *, real_number *);

    /* Local variables */
    static integer i__;
    static real_number t2, x2;
    static integer imax, last;
    static real_number tmax, alpha;
    extern /* Subroutine */ int _daskr_dcopy_(integer *, real_number *, integer *,
	    real_number *, integer *);
    static integer xroot, zroot, sgnchg;
    static integer imxold, nxlast;
    static real_number fracsub, fracint;


/* ***BEGIN PROLOGUE  DROOTS */
/* ***REFER TO DRCHEK */
/* ***ROUTINES CALLED DCOPY */
/* ***REVISION HISTORY  (YYMMDD) */
/*   020815  DATE WRITTEN */
/*   021217  Added root direction information in JROOT. */
/*   040518  Changed adjustment to X2 at 180 to avoid infinite loop. */
/* ***END PROLOGUE  DROOTS */

/* ----------------------------------------------------------------------- */
/* This subroutine finds the leftmost root of a set of arbitrary */
/* functions Ri(x) (i = 1,...,NRT) in an interval (X0,X1).  Only roots */
/* of odd multiplicity (i.e. changes of sign of the Ri) are found. */
/* Here the sign of X1 - X0 is arbitrary, but is constant for a given */
/* problem, and -leftmost- means nearest to X0. */
/* The values of the vector-valued function R(x) = (Ri, i=1...NRT) */
/* are communicated through the call sequence of DROOTS. */
/* The method used is the Illinois algorithm. */

/* Reference: */
/* Kathie L. Hiebert and Lawrence F. Shampine, Implicitly Defined */
/* Output Points for Solutions of ODEs, Sandia Report SAND80-0180, */
/* February 1980. */

/* Description of parameters. */

/* NRT    = number of functions Ri, or the number of components of */
/*          the vector valued function R(x).  Input only. */

/* HMIN   = resolution parameter in X.  Input only.  When a root is */
/*          found, it is located only to within an error of HMIN in X. */
/*          Typically, HMIN should be set to something on the order of */
/*               100 * UROUND * MAX(ABS(X0),ABS(X1)), */
/*          where UROUND is the unit roundoff of the machine. */

/* JFLAG  = integer flag for input and output communication. */

/*          On input, set JFLAG = 0 on the first call for the problem, */
/*          and leave it unchanged until the problem is completed. */
/*          (The problem is completed when JFLAG .ge. 2 on return.) */

/*          On output, JFLAG has the following values and meanings: */
/*          JFLAG = 1 means DROOTS needs a value of R(x).  Set RX = R(X) */
/*                    and call DROOTS again. */
/*          JFLAG = 2 means a root has been found.  The root is */
/*                    at X, and RX contains R(X).  (Actually, X is the */
/*                    rightmost approximation to the root on an interval */
/*                    (X0,X1) of size HMIN or less.) */
/*          JFLAG = 3 means X = X1 is a root, with one or more of the Ri */
/*                    being zero at X1 and no sign changes in (X0,X1). */
/*                    RX contains R(X) on output. */
/*          JFLAG = 4 means no roots (of odd multiplicity) were */
/*                    found in (X0,X1) (no sign changes). */

/* X0,X1  = endpoints of the interval where roots are sought. */
/*          X1 and X0 are input when JFLAG = 0 (first call), and */
/*          must be left unchanged between calls until the problem is */
/*          completed.  X0 and X1 must be distinct, but X1 - X0 may be */
/*          of either sign.  However, the notion of -left- and -right- */
/*          will be used to mean nearer to X0 or X1, respectively. */
/*          When JFLAG .ge. 2 on return, X0 and X1 are output, and */
/*          are the endpoints of the relevant interval. */

/* R0,R1  = arrays of length NRT containing the vectors R(X0) and R(X1), */
/*          respectively.  When JFLAG = 0, R0 and R1 are input and */
/*          none of the R0(i) should be zero. */
/*          When JFLAG .ge. 2 on return, R0 and R1 are output. */

/* RX     = array of length NRT containing R(X).  RX is input */
/*          when JFLAG = 1, and output when JFLAG .ge. 2. */

/* X      = independent variable value.  Output only. */
/*          When JFLAG = 1 on output, X is the point at which R(x) */
/*          is to be evaluated and loaded into RX. */
/*          When JFLAG = 2 or 3, X is the root. */
/*          When JFLAG = 4, X is the right endpoint of the interval, X1. */

/* JROOT  = integer array of length NRT.  Output only. */
/*          When JFLAG = 2 or 3, JROOT indicates which components */
/*          of R(x) have a root at X, and the direction of the sign */
/*          change across the root in the direction of integration. */
/*          JROOT(i) =  1 if Ri has a root and changes from - to +. */
/*          JROOT(i) = -1 if Ri has a root and changes from + to -. */
/*          Otherwise JROOT(i) = 0. */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --jroot;
    --rx;
    --r1;
    --r0;

    /* Function Body */

    if (*jflag == 1) {
	goto L200;
    }
/* JFLAG .ne. 1.  Check for change in sign of R or zero at X1. ---------- */
    imax = 0;
    tmax = zero;
    zroot = _FALSE_;
    i__1 = *nrt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if ((d__1 = r1[i__], fabs(d__1)) > zero) {
	    goto L110;
	}
	zroot = _TRUE_;
	goto L120;
/* At this point, R0(i) has been checked and cannot be zero. ------------ */
L110:
	if (_daskr_real_sign(&c_b758, &r0[i__]) == _daskr_real_sign(&c_b758, &r1[i__])) {
	    goto L120;
	}
	t2 = (d__1 = r1[i__] / (r1[i__] - r0[i__]), fabs(d__1));
	if (t2 <= tmax) {
	    goto L120;
	}
	tmax = t2;
	imax = i__;
L120:
	;
    }
    if (imax > 0) {
	goto L130;
    }
    sgnchg = _FALSE_;
    goto L140;
L130:
    sgnchg = _TRUE_;
L140:
    if (! sgnchg) {
	goto L400;
    }
/* There is a sign change.  Find the first root in the interval. -------- */
    xroot = _FALSE_;
    nxlast = 0;
    last = 1;

/* Repeat until the first root in the interval is found.  Loop point. --- */
L150:
    if (xroot) {
	goto L300;
    }
    if (nxlast == last) {
	goto L160;
    }
    alpha = 1.;
    goto L180;
L160:
    if (last == 0) {
	goto L170;
    }
    alpha *= .5;
    goto L180;
L170:
    alpha *= 2.;
L180:
    x2 = *x1 - (*x1 - *x0) * r1[imax] / (r1[imax] - alpha * r0[imax]);
    if ((d__1 = x2 - *x0, fabs(d__1)) < half * *hmin) {
	fracint = (d__1 = *x1 - *x0, fabs(d__1)) / *hmin;
	if (fracint > five) {
	    fracsub = tenth;
	} else {
	    fracsub = half / fracint;
	}
	x2 = *x0 + fracsub * (*x1 - *x0);
    }
    if ((d__1 = *x1 - x2, fabs(d__1)) < half * *hmin) {
	fracint = (d__1 = *x1 - *x0, fabs(d__1)) / *hmin;
	if (fracint > five) {
	    fracsub = tenth;
	} else {
	    fracsub = half / fracint;
	}
	x2 = *x1 - fracsub * (*x1 - *x0);
    }
    *jflag = 1;
    *x = x2;
/* Return to the calling routine to get a value of RX = R(X). ----------- */
    return 0;
/* Check to see in which interval R changes sign. ----------------------- */
L200:
    imxold = imax;
    imax = 0;
    tmax = zero;
    zroot = _FALSE_;
    i__1 = *nrt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if ((d__1 = rx[i__], fabs(d__1)) > zero) {
	    goto L210;
	}
	zroot = _TRUE_;
	goto L220;
/* Neither R0(i) nor RX(i) can be zero at this point. ------------------- */
L210:
	if (_daskr_real_sign(&c_b758, &r0[i__]) == _daskr_real_sign(&c_b758, &rx[i__])) {
	    goto L220;
	}
	t2 = (d__1 = rx[i__] / (rx[i__] - r0[i__]), fabs(d__1));
	if (t2 <= tmax) {
	    goto L220;
	}
	tmax = t2;
	imax = i__;
L220:
	;
    }
    if (imax > 0) {
	goto L230;
    }
    sgnchg = _FALSE_;
    imax = imxold;
    goto L240;
L230:
    sgnchg = _TRUE_;
L240:
    nxlast = last;
    if (! sgnchg) {
	goto L250;
    }
/* Sign change between X0 and X2, so replace X1 with X2. ---------------- */
    *x1 = x2;
    _daskr_dcopy_(nrt, &rx[1], &c__1, &r1[1], &c__1);
    last = 1;
    xroot = _FALSE_;
    goto L270;
L250:
    if (! zroot) {
	goto L260;
    }
/* Zero value at X2 and no sign change in (X0,X2), so X2 is a root. ----- */
    *x1 = x2;
    _daskr_dcopy_(nrt, &rx[1], &c__1, &r1[1], &c__1);
    xroot = _TRUE_;
    goto L270;
/* No sign change between X0 and X2.  Replace X0 with X2. --------------- */
L260:
    _daskr_dcopy_(nrt, &rx[1], &c__1, &r0[1], &c__1);
    *x0 = x2;
    last = 0;
    xroot = _FALSE_;
L270:
    if ((d__1 = *x1 - *x0, fabs(d__1)) <= *hmin) {
	xroot = _TRUE_;
    }
    goto L150;

/* Return with X1 as the root.  Set JROOT.  Set X = X1 and RX = R1. ----- */
L300:
    *jflag = 2;
    *x = *x1;
    _daskr_dcopy_(nrt, &r1[1], &c__1, &rx[1], &c__1);
    i__1 = *nrt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	jroot[i__] = 0;
	if ((d__1 = r1[i__], fabs(d__1)) == zero) {
	    jroot[i__] = (integer) (-_daskr_real_sign(&c_b758, &r0[i__]));
	    goto L320;
	}
	if (_daskr_real_sign(&c_b758, &r0[i__]) != _daskr_real_sign(&c_b758, &r1[i__])) {
	    d__1 = r1[i__] - r0[i__];
	    jroot[i__] = (integer) _daskr_real_sign(&c_b758, &d__1);
	}
L320:
	;
    }
    return 0;

/* No sign change in the interval.  Check for zero at right endpoint. --- */
L400:
    if (! zroot) {
	goto L420;
    }

/* Zero value at X1 and no sign change in (X0,X1).  Return JFLAG = 3. --- */
    *x = *x1;
    _daskr_dcopy_(nrt, &r1[1], &c__1, &rx[1], &c__1);
    i__1 = *nrt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	jroot[i__] = 0;
	if ((d__1 = r1[i__], fabs(d__1)) == zero) {
	    jroot[i__] = (integer) (-_daskr_real_sign(&c_b758, &r0[i__]));
	}
/* L410: */
    }
    *jflag = 3;
    return 0;

/* No sign changes in this interval.  Set X = X1, return JFLAG = 4. ----- */
L420:
    _daskr_dcopy_(nrt, &r1[1], &c__1, &rx[1], &c__1);
    *x = *x1;
    *jflag = 4;
    return 0;
/* ----------------------- END OF SUBROUTINE DROOTS ---------------------- */
} /* droots_ */

/* Subroutine */ int _daskr_ddasic_(real_number *x, real_number *y, real_number *yprime,
	 integer *neq, integer *icopt, integer *id, Unknown_fp res, Unknown_fp jac, Unknown_fp
	psol, real_number *h__, real_number *tscale, real_number *wt, integer *
	nic, integer *idid, real_number *rpar, integer *ipar, real_number *phi,
	real_number *savr, real_number *delta, real_number *e, real_number *yic,
	real_number *ypic, real_number *pwk, real_number *wm, integer *iwm,
	real_number *uround, real_number *epli, real_number *sqrtn, real_number *
	rsqrtn, real_number *epconi, real_number *stptol, integer *jflg,
	integer *icnflg, integer *icnstr, Unknown_fp nlsic)
{
    /* Initialized data */

    static real_number rhcut = .1;
    static real_number ratemx = .8;

    /* System generated locals */
    integer phi_dim1, phi_offset;

    /* Local variables */
    static real_number cj;
    static integer nh, mxnh;
    extern /* Subroutine */ int _daskr_dcopy_(integer *, real_number *, integer *,
	    real_number *, integer *);
    static integer jskip, iernls;


/* ***BEGIN PROLOGUE  DDASIC */
/* ***REFER TO  DDASPK */
/* ***DATE WRITTEN   940628   (YYMMDD) */
/* ***REVISION DATE  941206   (YYMMDD) */
/* ***REVISION DATE  950714   (YYMMDD) */
/* ***REVISION DATE  000628   TSCALE argument added. */

/* ----------------------------------------------------------------------- */
/* ***DESCRIPTION */

/*     DDASIC is a driver routine to compute consistent initial values */
/*     for Y and YPRIME.  There are two different options: */
/*     Denoting the differential variables in Y by Y_d, and */
/*     the algebraic variables by Y_a, the problem solved is either: */
/*     1.  Given Y_d, calculate Y_a and Y_d', or */
/*     2.  Given Y', calculate Y. */
/*     In either case, initial values for the given components */
/*     are input, and initial guesses for the unknown components */
/*     must also be provided as input. */

/*     The external routine NLSIC solves the resulting nonlinear system. */

/*     The parameters represent */

/*     X  --        Independent variable. */
/*     Y  --        Solution vector at X. */
/*     YPRIME --    Derivative of solution vector. */
/*     NEQ --       Number of equations to be integrated. */
/*     ICOPT     -- Flag indicating initial condition option chosen. */
/*                    ICOPT = 1 for option 1 above. */
/*                    ICOPT = 2 for option 2. */
/*     ID        -- Array of dimension NEQ, which must be initialized */
/*                  if option 1 is chosen. */
/*                    ID(i) = +1 if Y_i is a differential variable, */
/*                    ID(i) = -1 if Y_i is an algebraic variable. */
/*     RES --       External user-supplied subroutine to evaluate the */
/*                  residual.  See RES description in DDASPK prologue. */
/*     JAC --       External user-supplied routine to update Jacobian */
/*                  or preconditioner information in the nonlinear solver */
/*                  (optional).  See JAC description in DDASPK prologue. */
/*     PSOL --      External user-supplied routine to solve */
/*                  a linear system using preconditioning. */
/*                  See PSOL in DDASPK prologue. */
/*     H --         Scaling factor in iteration matrix.  DDASIC may */
/*                  reduce H to achieve convergence. */
/*     TSCALE --    Scale factor in T, used for stopping tests if nonzero. */
/*     WT --        Vector of weights for error criterion. */
/*     NIC --       Input number of initial condition calculation call */
/*                  (= 1 or 2). */
/*     IDID --      Completion code.  See IDID in DDASPK prologue. */
/*     RPAR,IPAR -- Real and integer parameter arrays that */
/*                  are used for communication between the */
/*                  calling program and external user routines. */
/*                  They are not altered by DNSK */
/*     PHI --       Work space for DDASIC of length at least 2*NEQ. */
/*     SAVR --      Work vector for DDASIC of length NEQ. */
/*     DELTA --     Work vector for DDASIC of length NEQ. */
/*     E --         Work vector for DDASIC of length NEQ. */
/*     YIC,YPIC --  Work vectors for DDASIC, each of length NEQ. */
/*     PWK --       Work vector for DDASIC of length NEQ. */
/*     WM,IWM --    Real and integer arrays storing */
/*                  information required by the linear solver. */
/*     EPCONI --    Test constant for Newton iteration convergence. */
/*     ICNFLG --    Flag showing whether constraints on Y are to apply. */
/*     ICNSTR --    Integer array of length NEQ with constraint types. */

/*     The other parameters are for use internally by DDASIC. */

/* ----------------------------------------------------------------------- */
/* ***ROUTINES CALLED */
/*   DCOPY, NLSIC */

/* ***END PROLOGUE  DDASIC */




/* The following parameters are data-loaded here: */
/*     RHCUT  = factor by which H is reduced on retry of Newton solve. */
/*     RATEMX = maximum convergence rate for which Newton iteration */
/*              is considered converging. */

    /* Parameter adjustments */
    --y;
    --yprime;
    phi_dim1 = *neq;
    phi_offset = 1 + phi_dim1;
    phi -= phi_offset;
    --id;
    --wt;
    --rpar;
    --ipar;
    --savr;
    --delta;
    --e;
    --yic;
    --ypic;
    --pwk;
    --wm;
    --iwm;
    --icnstr;

    /* Function Body */


/* ----------------------------------------------------------------------- */
/*     BLOCK 1. */
/*     Initializations. */
/*     JSKIP is a flag set to 1 when NIC = 2 and NH = 1, to signal that */
/*     the initial call to the JAC routine is to be skipped then. */
/*     Save Y and YPRIME in PHI.  Initialize IDID, NH, and CJ. */
/* ----------------------------------------------------------------------- */

    mxnh = iwm[34];
    *idid = 1;
    nh = 1;
    jskip = 0;
    if (*nic == 2) {
	jskip = 1;
    }
    _daskr_dcopy_(neq, &y[1], &c__1, &phi[phi_dim1 + 1], &c__1);
    _daskr_dcopy_(neq, &yprime[1], &c__1, &phi[(phi_dim1 << 1) + 1], &c__1);

    if (*icopt == 2) {
	cj = 0.;
    } else {
	cj = 1. / *h__;
    }

/* ----------------------------------------------------------------------- */
/*     BLOCK 2 */
/*     Call the nonlinear system solver to obtain */
/*     consistent initial values for Y and YPRIME. */
/* ----------------------------------------------------------------------- */

L200:
    ((int (*)(real_number *x, real_number *y, real_number *yprime,
	 integer *neq, integer *icopt, integer *id, Unknown_fp res, Unknown_fp jack, Unknown_fp
	psol, real_number *h__, real_number *tscale, real_number *wt, integer *
	jskip, real_number *rpar, integer *ipar, real_number *savr, real_number *
	delta, real_number *r__, real_number *yic, real_number *ypic, real_number
	*pwk, real_number *wm, integer *iwm, real_number *cj, real_number *
	uround, real_number *epli, real_number *sqrtn, real_number *rsqrtn,
	real_number *epcon, real_number *ratemx, real_number *stptol, integer *
	jflg, integer *icnflg, integer *icnstr, integer *iernls))
    *nlsic)
    (x, &y[1], &yprime[1], neq, icopt, &id[1], (Unknown_fp)res, (Unknown_fp)jac, (
	    Unknown_fp)psol, h__, tscale, &wt[1], &jskip, &rpar[1], &ipar[1], &savr[
	    1], &delta[1], &e[1], &yic[1], &ypic[1], &pwk[1], &wm[1], &iwm[1],
	     &cj, uround, epli, sqrtn, rsqrtn, epconi, &ratemx, stptol, jflg, 
	    icnflg, &icnstr[1], &iernls);

    if (iernls == 0) {
	return 0;
    }

/* ----------------------------------------------------------------------- */
/*     BLOCK 3 */
/*     The nonlinear solver was unsuccessful.  Increment NCFN. */
/*     Return with IDID = -12 if either */
/*       IERNLS = -1: error is considered unrecoverable, */
/*       ICOPT = 2: we are doing initialization problem type 2, or */
/*       NH = MXNH: the maximum number of H values has been tried. */
/*     Otherwise (problem 1 with IERNLS .GE. 1), reduce H and try again. */
/*     If IERNLS > 1, restore Y and YPRIME to their original values. */
/* ----------------------------------------------------------------------- */

    ++iwm[15];
    jskip = 0;

    if (iernls == -1) {
	goto L350;
    }
    if (*icopt == 2) {
	goto L350;
    }
    if (nh == mxnh) {
	goto L350;
    }

    ++nh;
    *h__ *= rhcut;
    cj = 1. / *h__;

    if (iernls == 1) {
	goto L200;
    }

    _daskr_dcopy_(neq, &phi[phi_dim1 + 1], &c__1, &y[1], &c__1);
    _daskr_dcopy_(neq, &phi[(phi_dim1 << 1) + 1], &c__1, &yprime[1], &c__1);
    goto L200;

L350:
    *idid = -12;
    return 0;

/* ------END OF SUBROUTINE DDASIC----------------------------------------- */
} /* ddasic_ */

/* Subroutine */ int _daskr_dyypnw_(integer *neq, real_number *y, real_number *yprime,
	real_number *cj, real_number *rl, real_number *p, integer *icopt,
	integer *id, real_number *ynew, real_number *ypnew)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;


/* ***BEGIN PROLOGUE  DYYPNW */
/* ***REFER TO  DLINSK */
/* ***DATE WRITTEN   940830   (YYMMDD) */


/* ----------------------------------------------------------------------- */
/* ***DESCRIPTION */

/*     DYYPNW calculates the new (Y,YPRIME) pair needed in the */
/*     linesearch algorithm based on the current lambda value.  It is */
/*     called by DLINSK and DLINSD.  Based on the ICOPT and ID values, */
/*     the corresponding entry in Y or YPRIME is updated. */

/*     In addition to the parameters described in the calling programs, */
/*     the parameters represent */

/*     P      -- Array of length NEQ that contains the current */
/*               approximate Newton step. */
/*     RL     -- Scalar containing the current lambda value. */
/*     YNEW   -- Array of length NEQ containing the updated Y vector. */
/*     YPNEW  -- Array of length NEQ containing the updated YPRIME */
/*               vector. */
/* ----------------------------------------------------------------------- */

/* ***ROUTINES CALLED (NONE) */

/* ***END PROLOGUE  DYYPNW */



    /* Parameter adjustments */
    --ypnew;
    --ynew;
    --id;
    --p;
    --yprime;
    --y;

    /* Function Body */
    if (*icopt == 1) {
	i__1 = *neq;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (id[i__] < 0) {
		ynew[i__] = y[i__] - *rl * p[i__];
		ypnew[i__] = yprime[i__];
	    } else {
		ynew[i__] = y[i__];
		ypnew[i__] = yprime[i__] - *rl * *cj * p[i__];
	    }
/* L10: */
	}
    } else {
	i__1 = *neq;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ynew[i__] = y[i__] - *rl * p[i__];
	    ypnew[i__] = yprime[i__];
/* L20: */
	}
    }
    return 0;
/* ----------------------- END OF SUBROUTINE DYYPNW ---------------------- */
} /* dyypnw_ */

/* Subroutine */ int _daskr_ddstp_(real_number *x, real_number *y, real_number *yprime,
	integer *neq, Unknown_fp res, Unknown_fp jac, Unknown_fp psol, real_number *h__,
	real_number *wt, real_number *vt, integer *jstart, integer *idid,
	real_number *rpar, integer *ipar, real_number *phi, real_number *savr,
	real_number *delta, real_number *e, real_number *wm, integer *iwm,
	real_number *alpha, real_number *beta, real_number *gamma, real_number *
	psi, real_number *sigma, real_number *cj, real_number *cjold, real_number
	*hold, real_number *s, real_number *hmin, real_number *uround,
	real_number *epli, real_number *sqrtn, real_number *rsqrtn, real_number *
	epcon, integer *iphase, integer *jcalc, integer *jflg, integer *k, 
	integer *kold, integer *ns, integer *nonneg, integer *ntype, Unknown_fp nls)
{
    /* System generated locals */
    integer phi_dim1, phi_offset, i__1, i__2;
    real_number d__1, d__2;

    /* Builtin functions */
    double _daskr_real_pow(real_number *, real_number *);

    /* Local variables */
    static integer i__, j;
    static real_number r__;
    static integer j1;
    static real_number ck;
    static integer km1, kp1, kp2, ncf, nef;
    static real_number erk, err, est;
    static integer nsp1;
    static real_number hnew, terk, xold;
    static integer knew;
    static real_number erkm1, erkm2, erkp1, temp1, temp2;
    static integer kdiff;
    static real_number enorm, alpha0, terkm1, terkm2, terkp1, alphas;
    extern /* Subroutine */ int _daskr_ddatrp_(real_number *, real_number *,
	    real_number *, real_number *, integer *, integer *, real_number *,
	    real_number *);
    static real_number cjlast;
    extern real_number _daskr_ddwnrm_(integer *, real_number *, real_number *,
	    real_number *, integer *);
    static integer iernls;


/* ***BEGIN PROLOGUE  DDSTP */
/* ***REFER TO  DDASPK */
/* ***DATE WRITTEN   890101   (YYMMDD) */
/* ***REVISION DATE  900926   (YYMMDD) */
/* ***REVISION DATE  940909   (YYMMDD) (Reset PSI(1), PHI(*,2) at 690) */


/* ----------------------------------------------------------------------- */
/* ***DESCRIPTION */

/*     DDSTP solves a system of differential/algebraic equations of */
/*     the form G(X,Y,YPRIME) = 0, for one step (normally from X to X+H). */

/*     The methods used are modified divided difference, fixed leading */
/*     coefficient forms of backward differentiation formulas. */
/*     The code adjusts the stepsize and order to control the local error */
/*     per step. */


/*     The parameters represent */
/*     X  --        Independent variable. */
/*     Y  --        Solution vector at X. */
/*     YPRIME --    Derivative of solution vector */
/*                  after successful step. */
/*     NEQ --       Number of equations to be integrated. */
/*     RES --       External user-supplied subroutine */
/*                  to evaluate the residual.  See RES description */
/*                  in DDASPK prologue. */
/*     JAC --       External user-supplied routine to update */
/*                  Jacobian or preconditioner information in the */
/*                  nonlinear solver.  See JAC description in DDASPK */
/*                  prologue. */
/*     PSOL --      External user-supplied routine to solve */
/*                  a linear system using preconditioning. */
/*                  (This is optional).  See PSOL in DDASPK prologue. */
/*     H --         Appropriate step size for next step. */
/*                  Normally determined by the code. */
/*     WT --        Vector of weights for error criterion used in Newton test. */
/*     VT --        Masked vector of weights used in error test. */
/*     JSTART --    Integer variable set 0 for */
/*                  first step, 1 otherwise. */
/*     IDID --      Completion code returned from the nonlinear solver. */
/*                  See IDID description in DDASPK prologue. */
/*     RPAR,IPAR -- Real and integer parameter arrays that */
/*                  are used for communication between the */
/*                  calling program and external user routines. */
/*                  They are not altered by DNSK */
/*     PHI --       Array of divided differences used by */
/*                  DDSTP. The length is NEQ*(K+1), where */
/*                  K is the maximum order. */
/*     SAVR --      Work vector for DDSTP of length NEQ. */
/*     DELTA,E --   Work vectors for DDSTP of length NEQ. */
/*     WM,IWM --    Real and integer arrays storing */
/*                  information required by the linear solver. */

/*     The other parameters are information */
/*     which is needed internally by DDSTP to */
/*     continue from step to step. */

/* ----------------------------------------------------------------------- */
/* ***ROUTINES CALLED */
/*   NLS, DDWNRM, DDATRP */

/* ***END PROLOGUE  DDSTP */





/* ----------------------------------------------------------------------- */
/*     BLOCK 1. */
/*     Initialize.  On the first call, set */
/*     the order to 1 and initialize */
/*     other variables. */
/* ----------------------------------------------------------------------- */

/*     Initializations for all calls */

    /* Parameter adjustments */
    --y;
    --yprime;
    phi_dim1 = *neq;
    phi_offset = 1 + phi_dim1;
    phi -= phi_offset;
    --wt;
    --vt;
    --rpar;
    --ipar;
    --savr;
    --delta;
    --e;
    --wm;
    --iwm;
    --alpha;
    --beta;
    --gamma;
    --psi;
    --sigma;

    /* Function Body */
    xold = *x;
    ncf = 0;
    nef = 0;
    if (*jstart != 0) {
	goto L120;
    }

/*     If this is the first step, perform */
/*     other initializations */

    *k = 1;
    *kold = 0;
    *hold = 0.;
    psi[1] = *h__;
    *cj = 1. / *h__;
    *iphase = 0;
    *ns = 0;
L120:





/* ----------------------------------------------------------------------- */
/*     BLOCK 2 */
/*     Compute coefficients of formulas for */
/*     this step. */
/* ----------------------------------------------------------------------- */
L200:
    kp1 = *k + 1;
    kp2 = *k + 2;
    km1 = *k - 1;
    if (*h__ != *hold || *k != *kold) {
	*ns = 0;
    }
/* Computing MIN */
    i__1 = *ns + 1, i__2 = *kold + 2;
    *ns = MIN(i__1,i__2);
    nsp1 = *ns + 1;
    if (kp1 < *ns) {
	goto L230;
    }

    beta[1] = 1.;
    alpha[1] = 1.;
    temp1 = *h__;
    gamma[1] = 0.;
    sigma[1] = 1.;
    i__1 = kp1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	temp2 = psi[i__ - 1];
	psi[i__ - 1] = temp1;
	beta[i__] = beta[i__ - 1] * psi[i__ - 1] / temp2;
	temp1 = temp2 + *h__;
	alpha[i__] = *h__ / temp1;
	sigma[i__] = (i__ - 1) * sigma[i__ - 1] * alpha[i__];
	gamma[i__] = gamma[i__ - 1] + alpha[i__ - 1] / *h__;
/* L210: */
    }
    psi[kp1] = temp1;
L230:

/*     Compute ALPHAS, ALPHA0 */

    alphas = 0.;
    alpha0 = 0.;
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	alphas -= 1. / i__;
	alpha0 -= alpha[i__];
/* L240: */
    }

/*     Compute leading coefficient CJ */

    cjlast = *cj;
    *cj = -alphas / *h__;

/*     Compute variable stepsize error coefficient CK */

    ck = (d__1 = alpha[kp1] + alphas - alpha0, fabs(d__1));
/* Computing MAX */
    d__1 = ck, d__2 = alpha[kp1];
    ck = MAX(d__1,d__2);

/*     Change PHI to PHI STAR */

    if (kp1 < nsp1) {
	goto L280;
    }
    i__1 = kp1;
    for (j = nsp1; j <= i__1; ++j) {
	i__2 = *neq;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L260: */
	    phi[i__ + j * phi_dim1] = beta[j] * phi[i__ + j * phi_dim1];
	}
/* L270: */
    }
L280:

/*     Update time */

    *x += *h__;

/*     Initialize IDID to 1 */

    *idid = 1;





/* ----------------------------------------------------------------------- */
/*     BLOCK 3 */
/*     Call the nonlinear system solver to obtain the solution and */
/*     derivative. */
/* ----------------------------------------------------------------------- */

    ((int (*)(real_number *x, real_number *y, real_number *yprime,
	integer *neq, Unknown_fp res, Unknown_fp jack, Unknown_fp psol, real_number *h__,
	real_number *wt, integer *jstart, integer *idid, real_number *rpar,
	integer *ipar, real_number *phi, real_number *gamma, real_number *savr,
	real_number *delta, real_number *e, real_number *wm, integer *iwm,
	real_number *cj, real_number *cjold, real_number *cjlast, real_number *s,
	real_number *uround, real_number *epli, real_number *sqrtn, real_number *
	rsqrtn, real_number *epcon, integer *jcalc, integer *jflg, integer *
	kp1, integer *nonneg, integer *ntype, integer *iernls))(*nls))

   (x, &y[1], &yprime[1], neq, (Unknown_fp)res, (Unknown_fp)jac, (Unknown_fp)psol, h__, &
	    wt[1], jstart, idid, &rpar[1], &ipar[1], &phi[phi_offset], &gamma[
	    1], &savr[1], &delta[1], &e[1], &wm[1], &iwm[1], cj, cjold, &
	    cjlast, s, uround, epli, sqrtn, rsqrtn, epcon, jcalc, jflg, &kp1, 
	    nonneg, ntype, &iernls);

    if (iernls != 0) {
	goto L600;
    }





/* ----------------------------------------------------------------------- */
/*     BLOCK 4 */
/*     Estimate the errors at orders K,K-1,K-2 */
/*     as if constant stepsize was used. Estimate */
/*     the local error at order K and test */
/*     whether the current step is successful. */
/* ----------------------------------------------------------------------- */

/*     Estimate errors at orders K,K-1,K-2 */

    enorm = _daskr_ddwnrm_(neq, &e[1], &vt[1], &rpar[1], &ipar[1]);
    erk = sigma[*k + 1] * enorm;
    terk = (*k + 1) * erk;
    est = erk;
    knew = *k;
    if (*k == 1) {
	goto L430;
    }
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L405: */
	delta[i__] = phi[i__ + kp1 * phi_dim1] + e[i__];
    }
    erkm1 = sigma[*k] * _daskr_ddwnrm_(neq, &delta[1], &vt[1], &rpar[1], &ipar[1]);
    terkm1 = *k * erkm1;
    if (*k > 2) {
	goto L410;
    }
    if (terkm1 <= terk * .5f) {
	goto L420;
    }
    goto L430;
L410:
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L415: */
	delta[i__] = phi[i__ + *k * phi_dim1] + delta[i__];
    }
    erkm2 = sigma[*k - 1] * _daskr_ddwnrm_(neq, &delta[1], &vt[1], &rpar[1], &ipar[1]
	    );
    terkm2 = (*k - 1) * erkm2;
    if (MAX(terkm1,terkm2) > terk) {
	goto L430;
    }

/*     Lower the order */

L420:
    knew = *k - 1;
    est = erkm1;


/*     Calculate the local error for the current step */
/*     to see if the step was successful */

L430:
    err = ck * enorm;
    if (err > 1.) {
	goto L600;
    }





/* ----------------------------------------------------------------------- */
/*     BLOCK 5 */
/*     The step is successful. Determine */
/*     the best order and stepsize for */
/*     the next step. Update the differences */
/*     for the next step. */
/* ----------------------------------------------------------------------- */
    *idid = 1;
    ++iwm[11];
    kdiff = *k - *kold;
    *kold = *k;
    *hold = *h__;


/*     Estimate the error at order K+1 unless */
/*        already decided to lower order, or */
/*        already using maximum order, or */
/*        stepsize not constant, or */
/*        order raised in previous step */

    if (knew == km1 || *k == iwm[3]) {
	*iphase = 1;
    }
    if (*iphase == 0) {
	goto L545;
    }
    if (knew == km1) {
	goto L540;
    }
    if (*k == iwm[3]) {
	goto L550;
    }
    if (kp1 >= *ns || kdiff == 1) {
	goto L550;
    }
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L510: */
	delta[i__] = e[i__] - phi[i__ + kp2 * phi_dim1];
    }
    erkp1 = 1. / (*k + 2) * _daskr_ddwnrm_(neq, &delta[1], &vt[1], &rpar[1], &ipar[1]
	    );
    terkp1 = (*k + 2) * erkp1;
    if (*k > 1) {
	goto L520;
    }
    if (terkp1 >= terk * .5) {
	goto L550;
    }
    goto L530;
L520:
    if (terkm1 <= MIN(terk,terkp1)) {
	goto L540;
    }
    if (terkp1 >= terk || *k == iwm[3]) {
	goto L550;
    }

/*     Raise order */

L530:
    *k = kp1;
    est = erkp1;
    goto L550;

/*     Lower order */

L540:
    *k = km1;
    est = erkm1;
    goto L550;

/*     If IPHASE = 0, increase order by one and multiply stepsize by */
/*     factor two */

L545:
    *k = kp1;
    hnew = *h__ * 2.;
    *h__ = hnew;
    goto L575;


/*     Determine the appropriate stepsize for */
/*     the next step. */

L550:
    hnew = *h__;
    temp2 = (real_number) (*k + 1);
    d__1 = est * 2. + 1e-4;
    d__2 = -1. / temp2;
    r__ = _daskr_real_pow(&d__1, &d__2);
    if (r__ < 2.) {
	goto L555;
    }
    hnew = *h__ * 2.;
    goto L560;
L555:
    if (r__ > 1.) {
	goto L560;
    }
/* Computing MAX */
    d__1 = .5, d__2 = MIN(.9,r__);
    r__ = MAX(d__1,d__2);
    hnew = *h__ * r__;
L560:
    *h__ = hnew;


/*     Update differences for next step */

L575:
    if (*kold == iwm[3]) {
	goto L585;
    }
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L580: */
	phi[i__ + kp2 * phi_dim1] = e[i__];
    }
L585:
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L590: */
	phi[i__ + kp1 * phi_dim1] += e[i__];
    }
    i__1 = kp1;
    for (j1 = 2; j1 <= i__1; ++j1) {
	j = kp1 - j1 + 1;
	i__2 = *neq;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L595: */
	    phi[i__ + j * phi_dim1] += phi[i__ + (j + 1) * phi_dim1];
	}
    }
    *jstart = 1;
    return 0;





/* ----------------------------------------------------------------------- */
/*     BLOCK 6 */
/*     The step is unsuccessful. Restore X,PSI,PHI */
/*     Determine appropriate stepsize for */
/*     continuing the integration, or exit with */
/*     an error flag if there have been many */
/*     failures. */
/* ----------------------------------------------------------------------- */
L600:
    *iphase = 1;

/*     Restore X,PHI,PSI */

    *x = xold;
    if (kp1 < nsp1) {
	goto L630;
    }
    i__2 = kp1;
    for (j = nsp1; j <= i__2; ++j) {
	temp1 = 1. / beta[j];
	i__1 = *neq;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L610: */
	    phi[i__ + j * phi_dim1] = temp1 * phi[i__ + j * phi_dim1];
	}
/* L620: */
    }
L630:
    i__2 = kp1;
    for (i__ = 2; i__ <= i__2; ++i__) {
/* L640: */
	psi[i__ - 1] = psi[i__] - *h__;
    }


/*     Test whether failure is due to nonlinear solver */
/*     or error test */

    if (iernls == 0) {
	goto L660;
    }
    ++iwm[15];


/*     The nonlinear solver failed to converge. */
/*     Determine the cause of the failure and take appropriate action. */
/*     If IERNLS .LT. 0, then return.  Otherwise, reduce the stepsize */
/*     and try again, unless too many failures have occurred. */

    if (iernls < 0) {
	goto L675;
    }
    ++ncf;
    r__ = .25;
    *h__ *= r__;
    if (ncf < 10 && fabs(*h__) >= *hmin) {
	goto L690;
    }
    if (*idid == 1) {
	*idid = -7;
    }
    if (nef >= 3) {
	*idid = -9;
    }
    goto L675;


/*     The nonlinear solver converged, and the cause */
/*     of the failure was the error estimate */
/*     exceeding the tolerance. */

L660:
    ++nef;
    ++iwm[14];
    if (nef > 1) {
	goto L665;
    }

/*     On first error test failure, keep current order or lower */
/*     order by one.  Compute new stepsize based on differences */
/*     of the solution. */

    *k = knew;
    temp2 = (real_number) (*k + 1);
    d__1 = est * 2. + 1e-4;
    d__2 = -1. / temp2;
    r__ = _daskr_real_pow(&d__1, &d__2) * .9;
/* Computing MAX */
    d__1 = .25, d__2 = MIN(.9,r__);
    r__ = MAX(d__1,d__2);
    *h__ *= r__;
    if (fabs(*h__) >= *hmin) {
	goto L690;
    }
    *idid = -6;
    goto L675;

/*     On second error test failure, use the current order or */
/*     decrease order by one.  Reduce the stepsize by a factor of */
/*     one quarter. */

L665:
    if (nef > 2) {
	goto L670;
    }
    *k = knew;
    r__ = .25;
    *h__ = r__ * *h__;
    if (fabs(*h__) >= *hmin) {
	goto L690;
    }
    *idid = -6;
    goto L675;

/*     On third and subsequent error test failures, set the order to */
/*     one, and reduce the stepsize by a factor of one quarter. */

L670:
    *k = 1;
    r__ = .25;
    *h__ = r__ * *h__;
    if (fabs(*h__) >= *hmin) {
	goto L690;
    }
    *idid = -6;
    goto L675;




/*     For all crashes, restore Y to its last value, */
/*     interpolate to find YPRIME at last X, and return. */

/*     Before returning, verify that the user has not set */
/*     IDID to a nonnegative value.  If the user has set IDID */
/*     to a nonnegative value, then reset IDID to be -7, indicating */
/*     a failure in the nonlinear system solver. */

L675:
    _daskr_ddatrp_(x, x, &y[1], &yprime[1], neq, k, &phi[phi_offset], &psi[1]);
    *jstart = 1;
    if (*idid >= 0) {
	*idid = -7;
    }
    return 0;


/*     Go back and try this step again. */
/*     If this is the first step, reset PSI(1) and rescale PHI(*,2). */

L690:
    if (*kold == 0) {
	psi[1] = *h__;
	i__2 = *neq;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L695: */
	    phi[i__ + (phi_dim1 << 1)] = r__ * phi[i__ + (phi_dim1 << 1)];
	}
    }
    goto L200;

/* ------END OF SUBROUTINE DDSTP------------------------------------------ */
} /* ddstp_ */

/* Subroutine */ int _daskr_dcnstr_(integer *neq, real_number *y, real_number *ynew,
	integer *icnstr, real_number *tau, real_number *rlx, integer *iret,
	integer *ivar)
{
    /* Initialized data */

    static real_number fac = .6;
    static real_number fac2 = .9;
    static real_number zero = 0.;

    /* System generated locals */
    integer i__1;
    real_number d__1;

    /* Local variables */
    static integer i__;
    static real_number rdy, rdymx;


/* ***BEGIN PROLOGUE  DCNSTR */
/* ***DATE WRITTEN   950808   (YYMMDD) */
/* ***REVISION DATE  950814   (YYMMDD) */


/* ----------------------------------------------------------------------- */
/* ***DESCRIPTION */

/* This subroutine checks for constraint violations in the proposed */
/* new approximate solution YNEW. */
/* If a constraint violation occurs, then a new step length, TAU, */
/* is calculated, and this value is to be given to the linesearch routine */
/* to calculate a new approximate solution YNEW. */

/* On entry: */

/*   NEQ    -- size of the nonlinear system, and the length of arrays */
/*             Y, YNEW and ICNSTR. */

/*   Y      -- real array containing the current approximate y. */

/*   YNEW   -- real array containing the new approximate y. */

/*   ICNSTR -- INTEGER array of length NEQ containing flags indicating */
/*             which entries in YNEW are to be constrained. */
/*             if ICNSTR(I) =  2, then YNEW(I) must be .GT. 0, */
/*             if ICNSTR(I) =  1, then YNEW(I) must be .GE. 0, */
/*             if ICNSTR(I) = -1, then YNEW(I) must be .LE. 0, while */
/*             if ICNSTR(I) = -2, then YNEW(I) must be .LT. 0, while */
/*             if ICNSTR(I) =  0, then YNEW(I) is not constrained. */

/*   RLX    -- real scalar restricting update, if ICNSTR(I) = 2 or -2, */
/*             to ABS( (YNEW-Y)/Y ) < FAC2*RLX in component I. */

/*   TAU    -- the current size of the step length for the linesearch. */

/* On return */

/*   TAU    -- the adjusted size of the step length if a constraint */
/*             violation occurred (otherwise, it is unchanged).  it is */
/*             the step length to give to the linesearch routine. */

/*   IRET   -- output flag. */
/*             IRET=0 means that YNEW satisfied all constraints. */
/*             IRET=1 means that YNEW failed to satisfy all the */
/*                    constraints, and a new linesearch step */
/*                    must be computed. */

/*   IVAR   -- index of variable causing constraint to be violated. */

/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --icnstr;
    --ynew;
    --y;

    /* Function Body */
/* ----------------------------------------------------------------------- */
/* Check constraints for proposed new step YNEW.  If a constraint has */
/* been violated, then calculate a new step length, TAU, to be */
/* used in the linesearch routine. */
/* ----------------------------------------------------------------------- */
    *iret = 0;
    rdymx = zero;
    *ivar = 0;
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {

	if (icnstr[i__] == 2) {
	    rdy = (d__1 = (ynew[i__] - y[i__]) / y[i__], fabs(d__1));
	    if (rdy > rdymx) {
		rdymx = rdy;
		*ivar = i__;
	    }
	    if (ynew[i__] <= zero) {
		*tau = fac * *tau;
		*ivar = i__;
		*iret = 1;
		return 0;
	    }

	} else if (icnstr[i__] == 1) {
	    if (ynew[i__] < zero) {
		*tau = fac * *tau;
		*ivar = i__;
		*iret = 1;
		return 0;
	    }

	} else if (icnstr[i__] == -1) {
	    if (ynew[i__] > zero) {
		*tau = fac * *tau;
		*ivar = i__;
		*iret = 1;
		return 0;
	    }

	} else if (icnstr[i__] == -2) {
	    rdy = (d__1 = (ynew[i__] - y[i__]) / y[i__], fabs(d__1));
	    if (rdy > rdymx) {
		rdymx = rdy;
		*ivar = i__;
	    }
	    if (ynew[i__] >= zero) {
		*tau = fac * *tau;
		*ivar = i__;
		*iret = 1;
		return 0;
	    }

	}
/* L100: */
    }
    if (rdymx >= *rlx) {
	*tau = fac2 * *tau * *rlx / rdymx;
	*iret = 1;
    }

    return 0;
/* ----------------------- END OF SUBROUTINE DCNSTR ---------------------- */
} /* dcnstr_ */

/* Subroutine */ int _daskr_dcnst0_(integer *neq, real_number *y, integer *icnstr,
	integer *iret)
{
    /* Initialized data */

    static real_number zero = 0.;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;


/* ***BEGIN PROLOGUE  DCNST0 */
/* ***DATE WRITTEN   950808   (YYMMDD) */
/* ***REVISION DATE  950808   (YYMMDD) */


/* ----------------------------------------------------------------------- */
/* ***DESCRIPTION */

/* This subroutine checks for constraint violations in the initial */
/* approximate solution u. */

/* On entry */

/*   NEQ    -- size of the nonlinear system, and the length of arrays */
/*             Y and ICNSTR. */

/*   Y      -- real array containing the initial approximate root. */

/*   ICNSTR -- INTEGER array of length NEQ containing flags indicating */
/*             which entries in Y are to be constrained. */
/*             if ICNSTR(I) =  2, then Y(I) must be .GT. 0, */
/*             if ICNSTR(I) =  1, then Y(I) must be .GE. 0, */
/*             if ICNSTR(I) = -1, then Y(I) must be .LE. 0, while */
/*             if ICNSTR(I) = -2, then Y(I) must be .LT. 0, while */
/*             if ICNSTR(I) =  0, then Y(I) is not constrained. */

/* On return */

/*   IRET   -- output flag. */
/*             IRET=0    means that u satisfied all constraints. */
/*             IRET.NE.0 means that Y(IRET) failed to satisfy its */
/*                       constraint. */

/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --icnstr;
    --y;

    /* Function Body */
/* ----------------------------------------------------------------------- */
/* Check constraints for initial Y.  If a constraint has been violated, */
/* set IRET = I to signal an error return to calling routine. */
/* ----------------------------------------------------------------------- */
    *iret = 0;
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (icnstr[i__] == 2) {
	    if (y[i__] <= zero) {
		*iret = i__;
		return 0;
	    }
	} else if (icnstr[i__] == 1) {
	    if (y[i__] < zero) {
		*iret = i__;
		return 0;
	    }
	} else if (icnstr[i__] == -1) {
	    if (y[i__] > zero) {
		*iret = i__;
		return 0;
	    }
	} else if (icnstr[i__] == -2) {
	    if (y[i__] >= zero) {
		*iret = i__;
		return 0;
	    }
	}
/* L100: */
    }
    return 0;
/* ----------------------- END OF SUBROUTINE DCNST0 ---------------------- */
} /* dcnst0_ */

/* Subroutine */ int _daskr_ddawts_(integer *neq, integer *iwt, real_number *rtol,
	real_number *atol, real_number *y, real_number *wt, real_number *rpar,
	integer *ipar)
{
    /* System generated locals */
    integer i__1;
    real_number d__1;

    /* Local variables */
    static integer i__;
    static real_number atoli, rtoli;


/* ***BEGIN PROLOGUE  DDAWTS */
/* ***REFER TO  DDASPK */
/* ***ROUTINES CALLED  (NONE) */
/* ***DATE WRITTEN   890101   (YYMMDD) */
/* ***REVISION DATE  900926   (YYMMDD) */
/* ***END PROLOGUE  DDAWTS */
/* ----------------------------------------------------------------------- */
/*     This subroutine sets the error weight vector, */
/*     WT, according to WT(I)=RTOL(I)*ABS(Y(I))+ATOL(I), */
/*     I = 1 to NEQ. */
/*     RTOL and ATOL are scalars if IWT = 0, */
/*     and vectors if IWT = 1. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --ipar;
    --rpar;
    --wt;
    --y;
    --atol;
    --rtol;

    /* Function Body */
    rtoli = rtol[1];
    atoli = atol[1];
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*iwt == 0) {
	    goto L10;
	}
	rtoli = rtol[i__];
	atoli = atol[i__];
L10:
	wt[i__] = rtoli * (d__1 = y[i__], fabs(d__1)) + atoli;
/* L20: */
    }
    return 0;

/* ------END OF SUBROUTINE DDAWTS----------------------------------------- */
} /* ddawts_ */

/* Subroutine */ int _daskr_dinvwt_(integer *neq, real_number *wt, integer *ier)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;


/* ***BEGIN PROLOGUE  DINVWT */
/* ***REFER TO  DDASPK */
/* ***ROUTINES CALLED  (NONE) */
/* ***DATE WRITTEN   950125   (YYMMDD) */
/* ***END PROLOGUE  DINVWT */
/* ----------------------------------------------------------------------- */
/*     This subroutine checks the error weight vector WT, of length NEQ, */
/*     for components that are .le. 0, and if none are found, it */
/*     inverts the WT(I) in place.  This replaces division operations */
/*     with multiplications in all norm evaluations. */
/*     IER is returned as 0 if all WT(I) were found positive, */
/*     and the first I with WT(I) .le. 0.0 otherwise. */
/* ----------------------------------------------------------------------- */


    /* Parameter adjustments */
    --wt;

    /* Function Body */
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (wt[i__] <= 0.) {
	    goto L30;
	}
/* L10: */
    }
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
	wt[i__] = 1. / wt[i__];
    }
    *ier = 0;
    return 0;

L30:
    *ier = i__;
    return 0;

/* ------END OF SUBROUTINE DINVWT----------------------------------------- */
} /* dinvwt_ */

/* Subroutine */ int _daskr_ddatrp_(real_number *x, real_number *xout, real_number *
	yout, real_number *ypout, integer *neq, integer *kold, real_number *phi,
	 real_number *psi)
{
    /* System generated locals */
    integer phi_dim1, phi_offset, i__1, i__2;

    /* Local variables */
    static real_number c__, d__;
    static integer i__, j;
    static real_number temp1, gamma;
    static integer koldp1;


/* ***BEGIN PROLOGUE  DDATRP */
/* ***REFER TO  DDASPK */
/* ***ROUTINES CALLED  (NONE) */
/* ***DATE WRITTEN   890101   (YYMMDD) */
/* ***REVISION DATE  900926   (YYMMDD) */
/* ***END PROLOGUE  DDATRP */

/* ----------------------------------------------------------------------- */
/*     The methods in subroutine DDSTP use polynomials */
/*     to approximate the solution.  DDATRP approximates the */
/*     solution and its derivative at time XOUT by evaluating */
/*     one of these polynomials, and its derivative, there. */
/*     Information defining this polynomial is passed from */
/*     DDSTP, so DDATRP cannot be used alone. */

/*     The parameters are */

/*     X     The current time in the integration. */
/*     XOUT  The time at which the solution is desired. */
/*     YOUT  The interpolated approximation to Y at XOUT. */
/*           (This is output.) */
/*     YPOUT The interpolated approximation to YPRIME at XOUT. */
/*           (This is output.) */
/*     NEQ   Number of equations. */
/*     KOLD  Order used on last successful step. */
/*     PHI   Array of scaled divided differences of Y. */
/*     PSI   Array of past stepsize history. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --yout;
    --ypout;
    phi_dim1 = *neq;
    phi_offset = 1 + phi_dim1;
    phi -= phi_offset;
    --psi;

    /* Function Body */
    koldp1 = *kold + 1;
    temp1 = *xout - *x;
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	yout[i__] = phi[i__ + phi_dim1];
/* L10: */
	ypout[i__] = 0.;
    }
    c__ = 1.;
    d__ = 0.;
    gamma = temp1 / psi[1];
    i__1 = koldp1;
    for (j = 2; j <= i__1; ++j) {
	d__ = d__ * gamma + c__ / psi[j - 1];
	c__ *= gamma;
	gamma = (temp1 + psi[j - 1]) / psi[j];
	i__2 = *neq;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    yout[i__] += c__ * phi[i__ + j * phi_dim1];
/* L20: */
	    ypout[i__] += d__ * phi[i__ + j * phi_dim1];
	}
/* L30: */
    }
    return 0;

/* ------END OF SUBROUTINE DDATRP----------------------------------------- */
} /* ddatrp_ */

real_number _daskr_ddwnrm_(integer *neq, real_number *v, real_number *rwt, real_number *
	rpar, integer *ipar)
{
    /* System generated locals */
    integer i__1;
    real_number ret_val, d__1, d__2;

    /* Local variables */
    static integer i__;
    static real_number sum, vmax;


/* ***BEGIN PROLOGUE  DDWNRM */
/* ***ROUTINES CALLED  (NONE) */
/* ***DATE WRITTEN   890101   (YYMMDD) */
/* ***REVISION DATE  900926   (YYMMDD) */
/* ***END PROLOGUE  DDWNRM */
/* ----------------------------------------------------------------------- */
/*     This function routine computes the weighted */
/*     root-mean-square norm of the vector of length */
/*     NEQ contained in the array V, with reciprocal weights */
/*     contained in the array RWT of length NEQ. */
/*        DDWNRM=SQRT((1/NEQ)*SUM(V(I)*RWT(I))**2) */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --ipar;
    --rpar;
    --rwt;
    --v;

    /* Function Body */
    ret_val = 0.;
    vmax = 0.;
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if ((d__1 = v[i__] * rwt[i__], fabs(d__1)) > vmax) {
	    vmax = (d__2 = v[i__] * rwt[i__], fabs(d__2));
	}
/* L10: */
    }
    if (vmax <= 0.) {
	goto L30;
    }
    sum = 0.;
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
/* Computing 2nd power */
	d__1 = v[i__] * rwt[i__] / vmax;
	sum += d__1 * d__1;
    }
    ret_val = vmax * sqrt(sum / *neq);
L30:
    return ret_val;

/* ------END OF FUNCTION DDWNRM------------------------------------------- */
} /* ddwnrm_ */

/* Subroutine */ int _daskr_ddasid_(real_number *x, real_number *y, real_number *yprime,
	 integer *neq, integer *icopt, integer *id, Unknown_fp res, Unknown_fp jacd,
	real_number *pdum, real_number *h__, real_number *tscale, real_number *wt,
	 integer *jsdum, real_number *rpar, integer *ipar, real_number *dumsvr,
	real_number *delta, real_number *r__, real_number *yic, real_number *ypic,
	 real_number *dumpwk, real_number *wm, integer *iwm, real_number *cj,
	real_number *uround, real_number *dume, real_number *dums, real_number *
	dumr, real_number *epcon, real_number *ratemx, real_number *stptol,
	integer *jfdum, integer *icnflg, integer *icnstr, integer *iernls)
{
    static integer nj, ierj, ires, mxnj;
    extern /* Subroutine */ int _daskr_dmatd_(integer *, real_number *, real_number *,
	    real_number *, real_number *, real_number *, real_number *, integer *,
	     real_number *, real_number *, real_number *, integer *, Unknown_fp,
	    integer *, real_number *, Unknown_fp, real_number *, integer *), _daskr_dnsid_(
	    real_number *, real_number *, real_number *, integer *, integer *,
	    integer *, Unknown_fp, real_number *, real_number *, integer *,
	    real_number *, real_number *, real_number *, real_number *,
	    real_number *, integer *, real_number *, real_number *, real_number *,
	     real_number *, integer *, real_number *, integer *, integer *,
	    integer *);
    static integer mxnit, iernew;


/* ***BEGIN PROLOGUE  DDASID */
/* ***REFER TO  DDASPK */
/* ***DATE WRITTEN   940701   (YYMMDD) */
/* ***REVISION DATE  950808   (YYMMDD) */
/* ***REVISION DATE  951110   Removed unreachable block 390. */
/* ***REVISION DATE  000628   TSCALE argument added. */


/* ----------------------------------------------------------------------- */
/* ***DESCRIPTION */


/*     DDASID solves a nonlinear system of algebraic equations of the */
/*     form G(X,Y,YPRIME) = 0 for the unknown parts of Y and YPRIME in */
/*     the initial conditions. */

/*     The method used is a modified Newton scheme. */

/*     The parameters represent */

/*     X         -- Independent variable. */
/*     Y         -- Solution vector. */
/*     YPRIME    -- Derivative of solution vector. */
/*     NEQ       -- Number of unknowns. */
/*     ICOPT     -- Initial condition option chosen (1 or 2). */
/*     ID        -- Array of dimension NEQ, which must be initialized */
/*                  if ICOPT = 1.  See DDASIC. */
/*     RES       -- External user-supplied subroutine to evaluate the */
/*                  residual.  See RES description in DDASPK prologue. */
/*     JACD      -- External user-supplied routine to evaluate the */
/*                  Jacobian.  See JAC description for the case */
/*                  INFO(12) = 0 in the DDASPK prologue. */
/*     PDUM      -- Dummy argument. */
/*     H         -- Scaling factor for this initial condition calc. */
/*     TSCALE    -- Scale factor in T, used for stopping tests if nonzero. */
/*     WT        -- Vector of weights for error criterion. */
/*     JSDUM     -- Dummy argument. */
/*     RPAR,IPAR -- Real and integer arrays used for communication */
/*                  between the calling program and external user */
/*                  routines.  They are not altered within DASPK. */
/*     DUMSVR    -- Dummy argument. */
/*     DELTA     -- Work vector for NLS of length NEQ. */
/*     R         -- Work vector for NLS of length NEQ. */
/*     YIC,YPIC  -- Work vectors for NLS, each of length NEQ. */
/*     DUMPWK    -- Dummy argument. */
/*     WM,IWM    -- Real and integer arrays storing matrix information */
/*                  such as the matrix of partial derivatives, */
/*                  permutation vector, and various other information. */
/*     CJ        -- Matrix parameter = 1/H (ICOPT = 1) or 0 (ICOPT = 2). */
/*     UROUND    -- Unit roundoff. */
/*     DUME      -- Dummy argument. */
/*     DUMS      -- Dummy argument. */
/*     DUMR      -- Dummy argument. */
/*     EPCON     -- Tolerance to test for convergence of the Newton */
/*                  iteration. */
/*     RATEMX    -- Maximum convergence rate for which Newton iteration */
/*                  is considered converging. */
/*     JFDUM     -- Dummy argument. */
/*     STPTOL    -- Tolerance used in calculating the minimum lambda */
/*                  value allowed. */
/*     ICNFLG    -- Integer scalar.  If nonzero, then constraint */
/*                  violations in the proposed new approximate solution */
/*                  will be checked for, and the maximum step length */
/*                  will be adjusted accordingly. */
/*     ICNSTR    -- Integer array of length NEQ containing flags for */
/*                  checking constraints. */
/*     IERNLS    -- Error flag for nonlinear solver. */
/*                   0   ==> nonlinear solver converged. */
/*                   1,2 ==> recoverable error inside nonlinear solver. */
/*                           1 => retry with current Y, YPRIME */
/*                           2 => retry with original Y, YPRIME */
/*                  -1   ==> unrecoverable error in nonlinear solver. */

/*     All variables with "DUM" in their names are dummy variables */
/*     which are not used in this routine. */

/* ----------------------------------------------------------------------- */

/* ***ROUTINES CALLED */
/*   RES, DMATD, DNSID */

/* ***END PROLOGUE  DDASID */





/*     Perform initializations. */

    /* Parameter adjustments */
    --icnstr;
    --iwm;
    --wm;
    --ypic;
    --yic;
    --r__;
    --delta;
    --ipar;
    --rpar;
    --wt;
    --id;
    --yprime;
    --y;

    /* Function Body */
    mxnit = iwm[32];
    mxnj = iwm[33];
    *iernls = 0;
    nj = 0;

/*     Call RES to initialize DELTA. */

    ires = 0;
    ++iwm[12];
    ((int (*)(real_number *, real_number *, real_number*, real_number*, real_number *, integer *, real_number *, integer *))*res)
      (x, &y[1], &yprime[1], cj, &delta[1], &ires, &rpar[1], &ipar[1]);
    if (ires < 0) {
	goto L370;
    }

/*     Looping point for updating the Jacobian. */

L300:

/*     Initialize all error flags to zero. */

    ierj = 0;
    ires = 0;
    iernew = 0;

/*     Reevaluate the iteration matrix, J = dG/dY + CJ*dG/dYPRIME, */
/*     where G(X,Y,YPRIME) = 0. */

    ++nj;
    ++iwm[13];
    _daskr_dmatd_(neq, x, &y[1], &yprime[1], &delta[1], cj, h__, &ierj, &wt[1], &r__[
	    1], &wm[1], &iwm[1], (Unknown_fp)res, &ires, uround, (Unknown_fp)jacd, &rpar[
	    1], &ipar[1]);
    if (ires < 0 || ierj != 0) {
	goto L370;
    }

/*     Call the nonlinear Newton solver for up to MXNIT iterations. */

    _daskr_dnsid_(x, &y[1], &yprime[1], neq, icopt, &id[1], (Unknown_fp)res, &wt[1], &rpar[
	    1], &ipar[1], &delta[1], &r__[1], &yic[1], &ypic[1], &wm[1], &iwm[
	    1], cj, tscale, epcon, ratemx, &mxnit, stptol, icnflg, &icnstr[1],
	     &iernew);

    if (iernew == 1 && nj < mxnj) {

/*        MXNIT iterations were done, the convergence rate is < 1, */
/*        and the number of Jacobian evaluations is less than MXNJ. */
/*        Call RES, reevaluate the Jacobian, and try again. */

	++iwm[12];
	((int (*)(real_number *, real_number *, real_number*, real_number*, real_number *, integer *, real_number *, integer *))*res)
    (x, &y[1], &yprime[1], cj, &delta[1], &ires, &rpar[1], &ipar[1])
		;
	if (ires < 0) {
	    goto L370;
	}
	goto L300;
    }

    if (iernew != 0) {
	goto L380;
    }
    return 0;


/*     Unsuccessful exits from nonlinear solver. */
/*     Compute IERNLS accordingly. */

L370:
    *iernls = 2;
    if (ires <= -2) {
	*iernls = -1;
    }
    return 0;

L380:
    *iernls = MIN(iernew,2);
    return 0;

/* ------END OF SUBROUTINE DDASID----------------------------------------- */
} /* ddasid_ */

/* Subroutine */ int _daskr_dnsid_(real_number *x, real_number *y, real_number *yprime,
	integer *neq, integer *icopt, integer *id, Unknown_fp res, real_number *wt,
	real_number *rpar, integer *ipar, real_number *delta, real_number *r__,
	real_number *yic, real_number *ypic, real_number *wm, integer *iwm,
	real_number *cj, real_number *tscale, real_number *epcon, real_number *
	ratemx, integer *maxit, real_number *stptol, integer *icnflg, integer *
	icnstr, integer *iernew)
{
    static integer m;
    static real_number rlx, rate, fnrm;
    static integer iret, ires, lsoff;
    extern /* Subroutine */ int _daskr_dslvd_(integer *, real_number *, real_number *,
	    integer *), _daskr_dcopy_(integer *, real_number *, integer *, real_number
	    *, integer *), _daskr_dlinsd_(integer *, real_number *, real_number *,
	    real_number *, real_number *, real_number *, real_number *,
	    real_number *, real_number *, integer *, real_number *, integer *,
	    Unknown_fp, integer *, real_number *, integer *, real_number *, integer *,
	     integer *, real_number *, real_number *, real_number *, integer *,
	    integer *, real_number *, real_number *, integer *);
    static real_number oldfnm, delnrm;
    extern real_number _daskr_ddwnrm_(integer *, real_number *, real_number *,
	    real_number *, integer *);


/* ***BEGIN PROLOGUE  DNSID */
/* ***REFER TO  DDASPK */
/* ***DATE WRITTEN   940701   (YYMMDD) */
/* ***REVISION DATE  950713   (YYMMDD) */
/* ***REVISION DATE  000628   TSCALE argument added. */


/* ----------------------------------------------------------------------- */
/* ***DESCRIPTION */

/*     DNSID solves a nonlinear system of algebraic equations of the */
/*     form G(X,Y,YPRIME) = 0 for the unknown parts of Y and YPRIME */
/*     in the initial conditions. */

/*     The method used is a modified Newton scheme. */

/*     The parameters represent */

/*     X         -- Independent variable. */
/*     Y         -- Solution vector. */
/*     YPRIME    -- Derivative of solution vector. */
/*     NEQ       -- Number of unknowns. */
/*     ICOPT     -- Initial condition option chosen (1 or 2). */
/*     ID        -- Array of dimension NEQ, which must be initialized */
/*                  if ICOPT = 1.  See DDASIC. */
/*     RES       -- External user-supplied subroutine to evaluate the */
/*                  residual.  See RES description in DDASPK prologue. */
/*     WT        -- Vector of weights for error criterion. */
/*     RPAR,IPAR -- Real and integer arrays used for communication */
/*                  between the calling program and external user */
/*                  routines.  They are not altered within DASPK. */
/*     DELTA     -- Residual vector on entry, and work vector of */
/*                  length NEQ for DNSID. */
/*     WM,IWM    -- Real and integer arrays storing matrix information */
/*                  such as the matrix of partial derivatives, */
/*                  permutation vector, and various other information. */
/*     CJ        -- Matrix parameter = 1/H (ICOPT = 1) or 0 (ICOPT = 2). */
/*     TSCALE    -- Scale factor in T, used for stopping tests if nonzero. */
/*     R         -- Array of length NEQ used as workspace by the */
/*                  linesearch routine DLINSD. */
/*     YIC,YPIC  -- Work vectors for DLINSD, each of length NEQ. */
/*     EPCON     -- Tolerance to test for convergence of the Newton */
/*                  iteration. */
/*     RATEMX    -- Maximum convergence rate for which Newton iteration */
/*                  is considered converging. */
/*     MAXIT     -- Maximum allowed number of Newton iterations. */
/*     STPTOL    -- Tolerance used in calculating the minimum lambda */
/*                  value allowed. */
/*     ICNFLG    -- Integer scalar.  If nonzero, then constraint */
/*                  violations in the proposed new approximate solution */
/*                  will be checked for, and the maximum step length */
/*                  will be adjusted accordingly. */
/*     ICNSTR    -- Integer array of length NEQ containing flags for */
/*                  checking constraints. */
/*     IERNEW    -- Error flag for Newton iteration. */
/*                   0  ==> Newton iteration converged. */
/*                   1  ==> failed to converge, but RATE .le. RATEMX. */
/*                   2  ==> failed to converge, RATE .gt. RATEMX. */
/*                   3  ==> other recoverable error (IRES = -1, or */
/*                          linesearch failed). */
/*                  -1  ==> unrecoverable error (IRES = -2). */

/* ----------------------------------------------------------------------- */

/* ***ROUTINES CALLED */
/*   DSLVD, DDWNRM, DLINSD, DCOPY */

/* ***END PROLOGUE  DNSID */





/*     Initializations.  M is the Newton iteration counter. */

    /* Parameter adjustments */
    --icnstr;
    --iwm;
    --wm;
    --ypic;
    --yic;
    --r__;
    --delta;
    --ipar;
    --rpar;
    --wt;
    --id;
    --yprime;
    --y;

    /* Function Body */
    lsoff = iwm[35];
    m = 0;
    rate = 1.;
    rlx = .4;

/*     Compute a new step vector DELTA by back-substitution. */

    _daskr_dslvd_(neq, &delta[1], &wm[1], &iwm[1]);

/*     Get norm of DELTA.  Return now if norm(DELTA) .le. EPCON. */

    delnrm = _daskr_ddwnrm_(neq, &delta[1], &wt[1], &rpar[1], &ipar[1]);
    fnrm = delnrm;
    if (*tscale > 0.) {
	fnrm = fnrm * *tscale * fabs(*cj);
    }
    if (fnrm <= *epcon) {
	return 0;
    }

/*     Newton iteration loop. */

L300:
    ++iwm[19];

/*     Call linesearch routine for global strategy and set RATE */

    oldfnm = fnrm;

    _daskr_dlinsd_(neq, &y[1], x, &yprime[1], cj, tscale, &delta[1], &delnrm, &wt[1],
	     &lsoff, stptol, &iret, (Unknown_fp)res, &ires, &wm[1], &iwm[1], &fnrm,
	    icopt, &id[1], &r__[1], &yic[1], &ypic[1], icnflg, &icnstr[1], &
	    rlx, &rpar[1], &ipar[1]);

    rate = fnrm / oldfnm;

/*     Check for error condition from linesearch. */
    if (iret != 0) {
	goto L390;
    }

/*     Test for convergence of the iteration, and return or loop. */

    if (fnrm <= *epcon) {
	return 0;
    }

/*     The iteration has not yet converged.  Update M. */
/*     Test whether the maximum number of iterations have been tried. */

    ++m;
    if (m >= *maxit) {
	goto L380;
    }

/*     Copy the residual to DELTA and its norm to DELNRM, and loop for */
/*     another iteration. */

    _daskr_dcopy_(neq, &r__[1], &c__1, &delta[1], &c__1);
    delnrm = fnrm;
    goto L300;

/*     The maximum number of iterations was done.  Set IERNEW and return. */

L380:
    if (rate <= *ratemx) {
	*iernew = 1;
    } else {
	*iernew = 2;
    }
    return 0;

L390:
    if (ires <= -2) {
	*iernew = -1;
    } else {
	*iernew = 3;
    }
    return 0;


/* ------END OF SUBROUTINE DNSID------------------------------------------ */
} /* dnsid_ */

/* Subroutine */ int _daskr_dlinsd_(integer *neq, real_number *y, real_number *t,
	real_number *yprime, real_number *cj, real_number *tscale, real_number *p,
	 real_number *pnrm, real_number *wt, integer *lsoff, real_number *stptol,
	 integer *iret, Unknown_fp res, integer *ires, real_number *wm, integer *iwm,
	 real_number *fnrm, integer *icopt, integer *id, real_number *r__,
	real_number *ynew, real_number *ypnew, integer *icnflg, integer *icnstr,
	 real_number *rlx, real_number *rpar, integer *ipar)
{
    /* Initialized data */

    static real_number alpha = 1e-4;
    static real_number one = 1.;
    static real_number two = 2.;

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int _daskr_str_copy(char *, char *, integer, integer);

    /* Local variables */
    static integer i__;
    static real_number rl;
    static char msg[80];
    static real_number tau;
    static integer ivar;
    static real_number slpi, f1nrm, ratio;
    extern /* Subroutine */ int _daskr_dcopy_(integer *, real_number *, integer *,
	    real_number *, integer *);
    static real_number rlmin, fnrmp;
    static integer kprin;
    static real_number ratio1, f1nrmp;
    extern /* Subroutine */ int _daskr_dfnrmd_(integer *, real_number *, real_number *,
	     real_number *, real_number *, real_number *, real_number *,
	    real_number *, Unknown_fp, integer *, real_number *, real_number *,
	    integer *, real_number *, integer *), _daskr_dcnstr_(integer *,
	    real_number *, real_number *, integer *, real_number *, real_number *,
	     integer *, integer *), _daskr_xerrwd_(char *, integer *, integer *,
	    integer *, integer *, integer *, integer *, integer *, real_number
	    *, real_number *, integer), _daskr_dyypnw_(integer *, real_number *,
	    real_number *, real_number *, real_number *, real_number *, integer *,
	     integer *, real_number *, real_number *);


/* ***BEGIN PROLOGUE  DLINSD */
/* ***REFER TO  DNSID */
/* ***DATE WRITTEN   941025   (YYMMDD) */
/* ***REVISION DATE  941215   (YYMMDD) */
/* ***REVISION DATE  960129   Moved line RL = ONE to top block. */
/* ***REVISION DATE  000628   TSCALE argument added. */


/* ----------------------------------------------------------------------- */
/* ***DESCRIPTION */

/*     DLINSD uses a linesearch algorithm to calculate a new (Y,YPRIME) */
/*     pair (YNEW,YPNEW) such that */

/*     f(YNEW,YPNEW) .le. (1 - 2*ALPHA*RL)*f(Y,YPRIME) , */

/*     where 0 < RL <= 1.  Here, f(y,y') is defined as */

/*      f(y,y') = (1/2)*norm( (J-inverse)*G(t,y,y') )**2 , */

/*     where norm() is the weighted RMS vector norm, G is the DAE */
/*     system residual function, and J is the system iteration matrix */
/*     (Jacobian). */

/*     In addition to the parameters defined elsewhere, we have */

/*     TSCALE  --  Scale factor in T, used for stopping tests if nonzero. */
/*     P       -- Approximate Newton step used in backtracking. */
/*     PNRM    -- Weighted RMS norm of P. */
/*     LSOFF   -- Flag showing whether the linesearch algorithm is */
/*                to be invoked.  0 means do the linesearch, and */
/*                1 means turn off linesearch. */
/*     STPTOL  -- Tolerance used in calculating the minimum lambda */
/*                value allowed. */
/*     ICNFLG  -- Integer scalar.  If nonzero, then constraint violations */
/*                in the proposed new approximate solution will be */
/*                checked for, and the maximum step length will be */
/*                adjusted accordingly. */
/*     ICNSTR  -- Integer array of length NEQ containing flags for */
/*                checking constraints. */
/*     RLX     -- Real scalar restricting update size in DCNSTR. */
/*     YNEW    -- Array of length NEQ used to hold the new Y in */
/*                performing the linesearch. */
/*     YPNEW   -- Array of length NEQ used to hold the new YPRIME in */
/*                performing the linesearch. */
/*     Y       -- Array of length NEQ containing the new Y (i.e.,=YNEW). */
/*     YPRIME  -- Array of length NEQ containing the new YPRIME */
/*                (i.e.,=YPNEW). */
/*     FNRM    -- Real scalar containing SQRT(2*f(Y,YPRIME)) for the */
/*                current (Y,YPRIME) on input and output. */
/*     R       -- Work array of length NEQ, containing the scaled */
/*                residual (J-inverse)*G(t,y,y') on return. */
/*     IRET    -- Return flag. */
/*                IRET=0 means that a satisfactory (Y,YPRIME) was found. */
/*                IRET=1 means that the routine failed to find a new */
/*                       (Y,YPRIME) that was sufficiently distinct from */
/*                       the current (Y,YPRIME) pair. */
/*                IRET=2 means IRES .ne. 0 from RES. */
/* ----------------------------------------------------------------------- */

/* ***ROUTINES CALLED */
/*   DFNRMD, DYYPNW, DCNSTR, DCOPY, XERRWD */

/* ***END PROLOGUE  DLINSD */



    /* Parameter adjustments */
    --ipar;
    --rpar;
    --icnstr;
    --ypnew;
    --ynew;
    --r__;
    --id;
    --iwm;
    --wm;
    --wt;
    --p;
    --yprime;
    --y;

    /* Function Body */

    kprin = iwm[31];

    f1nrm = *fnrm * *fnrm / two;
    ratio = one;
    if (kprin >= 2) {
	_daskr_str_copy(msg, "------ IN ROUTINE DLINSD-- PNRM = (R1)", (integer)80, (
		integer)38);
	_daskr_xerrwd_(msg, &c__38, &c__901, &c__0, &c__0, &c__0, &c__0, &c__1, pnrm,
		 &c_b38, (integer)80);
    }
    tau = *pnrm;
    rl = one;
/* ----------------------------------------------------------------------- */
/* Check for violations of the constraints, if any are imposed. */
/* If any violations are found, the step vector P is rescaled, and the */
/* constraint check is repeated, until no violations are found. */
/* ----------------------------------------------------------------------- */
    if (*icnflg != 0) {
L10:
	_daskr_dyypnw_(neq, &y[1], &yprime[1], cj, &rl, &p[1], icopt, &id[1], &ynew[
		1], &ypnew[1]);
	_daskr_dcnstr_(neq, &y[1], &ynew[1], &icnstr[1], &tau, rlx, iret, &ivar);
	if (*iret == 1) {
	    ratio1 = tau / *pnrm;
	    ratio *= ratio1;
	    i__1 = *neq;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
		p[i__] *= ratio1;
	    }
	    *pnrm = tau;
	    if (kprin >= 2) {
		_daskr_str_copy(msg, "------ CONSTRAINT VIOL., PNRM = (R1), INDEX = ("
			"I1)", (integer)80, (integer)50);
		_daskr_xerrwd_(msg, &c__50, &c__902, &c__0, &c__1, &ivar, &c__0, &
			c__1, pnrm, &c_b38, (integer)80);
	    }
	    if (*pnrm <= *stptol) {
		*iret = 1;
		return 0;
	    }
	    goto L10;
	}
    }

    slpi = -two * f1nrm * ratio;
    rlmin = *stptol / *pnrm;
    if (*lsoff == 0 && kprin >= 2) {
	_daskr_str_copy(msg, "------ MIN. LAMBDA = (R1)", (integer)80, (integer)25);
	_daskr_xerrwd_(msg, &c__25, &c__903, &c__0, &c__0, &c__0, &c__0, &c__1, &
		rlmin, &c_b38, (integer)80);
    }
/* ----------------------------------------------------------------------- */
/* Begin iteration to find RL value satisfying alpha-condition. */
/* If RL becomes less than RLMIN, then terminate with IRET = 1. */
/* ----------------------------------------------------------------------- */
L100:
    _daskr_dyypnw_(neq, &y[1], &yprime[1], cj, &rl, &p[1], icopt, &id[1], &ynew[1], &
	    ypnew[1]);
    _daskr_dfnrmd_(neq, &ynew[1], t, &ypnew[1], &r__[1], cj, tscale, &wt[1], (Unknown_fp)
	    res, ires, &fnrmp, &wm[1], &iwm[1], &rpar[1], &ipar[1]);
    ++iwm[12];
    if (*ires != 0) {
	*iret = 2;
	return 0;
    }
    if (*lsoff == 1) {
	goto L150;
    }

    f1nrmp = fnrmp * fnrmp / two;
    if (kprin >= 2) {
	_daskr_str_copy(msg, "------ LAMBDA = (R1)", (integer)80, (integer)20);
	_daskr_xerrwd_(msg, &c__20, &c__904, &c__0, &c__0, &c__0, &c__0, &c__1, &rl,
		&c_b38, (integer)80);
	_daskr_str_copy(msg, "------ NORM(F1) = (R1),  NORM(F1NEW) = (R2)", (integer)80,
		 (integer)43);
	_daskr_xerrwd_(msg, &c__43, &c__905, &c__0, &c__0, &c__0, &c__0, &c__2, &
		f1nrm, &f1nrmp, (integer)80);
    }
    if (f1nrmp > f1nrm + alpha * slpi * rl) {
	goto L200;
    }
/* ----------------------------------------------------------------------- */
/* Alpha-condition is satisfied, or linesearch is turned off. */
/* Copy YNEW,YPNEW to Y,YPRIME and return. */
/* ----------------------------------------------------------------------- */
L150:
    *iret = 0;
    _daskr_dcopy_(neq, &ynew[1], &c__1, &y[1], &c__1);
    _daskr_dcopy_(neq, &ypnew[1], &c__1, &yprime[1], &c__1);
    *fnrm = fnrmp;
    if (kprin >= 1) {
	_daskr_str_copy(msg, "------ LEAVING ROUTINE DLINSD, FNRM = (R1)", (integer)80,
		(integer)42);
	_daskr_xerrwd_(msg, &c__42, &c__906, &c__0, &c__0, &c__0, &c__0, &c__1, fnrm,
		 &c_b38, (integer)80);
    }
    return 0;
/* ----------------------------------------------------------------------- */
/* Alpha-condition not satisfied.  Perform backtrack to compute new RL */
/* value.  If no satisfactory YNEW,YPNEW can be found sufficiently */
/* distinct from Y,YPRIME, then return IRET = 1. */
/* ----------------------------------------------------------------------- */
L200:
    if (rl < rlmin) {
	*iret = 1;
	return 0;
    }

    rl /= two;
    goto L100;

/* ----------------------- END OF SUBROUTINE DLINSD ---------------------- */
} /* dlinsd_ */

/* Subroutine */ int _daskr_dfnrmd_(integer *neq, real_number *y, real_number *t,
	real_number *yprime, real_number *r__, real_number *cj, real_number *
	tscale, real_number *wt, Unknown_fp res, integer *ires, real_number *fnorm,
	real_number *wm, integer *iwm, real_number *rpar, integer *ipar)
{
    extern /* Subroutine */ int _daskr_dslvd_(integer *, real_number *, real_number *,
	    integer *);
    extern real_number _daskr_ddwnrm_(integer *, real_number *, real_number *,
	    real_number *, integer *);


/* ***BEGIN PROLOGUE  DFNRMD */
/* ***REFER TO  DLINSD */
/* ***DATE WRITTEN   941025   (YYMMDD) */
/* ***REVISION DATE  000628   TSCALE argument added. */


/* ----------------------------------------------------------------------- */
/* ***DESCRIPTION */

/*     DFNRMD calculates the scaled preconditioned norm of the nonlinear */
/*     function used in the nonlinear iteration for obtaining consistent */
/*     initial conditions.  Specifically, DFNRMD calculates the weighted */
/*     root-mean-square norm of the vector (J-inverse)*G(T,Y,YPRIME), */
/*     where J is the Jacobian matrix. */

/*     In addition to the parameters described in the calling program */
/*     DLINSD, the parameters represent */

/*     R      -- Array of length NEQ that contains */
/*               (J-inverse)*G(T,Y,YPRIME) on return. */
/*     TSCALE -- Scale factor in T, used for stopping tests if nonzero. */
/*     FNORM  -- Scalar containing the weighted norm of R on return. */
/* ----------------------------------------------------------------------- */

/* ***ROUTINES CALLED */
/*   RES, DSLVD, DDWNRM */

/* ***END PROLOGUE  DFNRMD */


/* ----------------------------------------------------------------------- */
/*     Call RES routine. */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --ipar;
    --rpar;
    --iwm;
    --wm;
    --wt;
    --r__;
    --yprime;
    --y;

    /* Function Body */
    *ires = 0;
    ((int (*)(real_number *, real_number *, real_number*, real_number*, real_number *, integer *, real_number *, integer *))*res)
      (t, &y[1], &yprime[1], cj, &r__[1], ires, &rpar[1], &ipar[1]);
    if (*ires < 0) {
	return 0;
    }
/* ----------------------------------------------------------------------- */
/*     Apply inverse of Jacobian to vector R. */
/* ----------------------------------------------------------------------- */
    _daskr_dslvd_(neq, &r__[1], &wm[1], &iwm[1]);
/* ----------------------------------------------------------------------- */
/*     Calculate norm of R. */
/* ----------------------------------------------------------------------- */
    *fnorm = _daskr_ddwnrm_(neq, &r__[1], &wt[1], &rpar[1], &ipar[1]);
    if (*tscale > 0.) {
	*fnorm = *fnorm * *tscale * fabs(*cj);
    }

    return 0;
/* ----------------------- END OF SUBROUTINE DFNRMD ---------------------- */
} /* dfnrmd_ */

/* Subroutine */ int _daskr_dnedd_(real_number *x, real_number *y, real_number *yprime,
	integer *neq, Unknown_fp res, Unknown_fp jacd, real_number *pdum, real_number *h__,
	real_number *wt, integer *jstart, integer *idid, real_number *rpar,
	integer *ipar, real_number *phi, real_number *gamma, real_number *dumsvr,
	 real_number *delta, real_number *e, real_number *wm, integer *iwm,
	real_number *cj, real_number *cjold, real_number *cjlast, real_number *s,
	real_number *uround, real_number *dume, real_number *dums, real_number *
	dumr, real_number *epcon, integer *jcalc, integer *jfdum, integer *kp1,
	 integer *nonneg, integer *ntype, integer *iernls)
{
    /* Initialized data */

    static integer muldel = 1;
    static integer maxit = 4;
    static real_number xrate = .25;

    /* System generated locals */
    integer phi_dim1, phi_offset, i__1, i__2;
    real_number d__1;

    /* Local variables */
    static integer i__, j, ierj;
    extern /* Subroutine */ int _daskr_dnsd_(real_number *, real_number *, real_number *
	    , integer *, Unknown_fp, real_number *, real_number *, real_number *,
	    integer *, real_number *, real_number *, real_number *, real_number *,
	     integer *, real_number *, real_number *, real_number *, real_number *
	    , real_number *, real_number *, real_number *, real_number *, integer
	    *, integer *, integer *, integer *, integer *);
    static integer idum, ires;
    static real_number temp1, temp2;
    extern /* Subroutine */ int _daskr_dmatd_(integer *, real_number *, real_number *,
	    real_number *, real_number *, real_number *, real_number *, integer *,
	     real_number *, real_number *, real_number *, integer *, Unknown_fp,
	    integer *, real_number *, Unknown_fp, real_number *, integer *);
    static real_number pnorm, delnrm;
    static integer iernew;
    extern real_number _daskr_ddwnrm_(integer *, real_number *, real_number *,
	    real_number *, integer *);
    static real_number tolnew;
    static integer iertyp;


/* ***BEGIN PROLOGUE  DNEDD */
/* ***REFER TO  DDASPK */
/* ***DATE WRITTEN   891219   (YYMMDD) */
/* ***REVISION DATE  900926   (YYMMDD) */


/* ----------------------------------------------------------------------- */
/* ***DESCRIPTION */

/*     DNEDD solves a nonlinear system of */
/*     algebraic equations of the form */
/*     G(X,Y,YPRIME) = 0 for the unknown Y. */

/*     The method used is a modified Newton scheme. */

/*     The parameters represent */

/*     X         -- Independent variable. */
/*     Y         -- Solution vector. */
/*     YPRIME    -- Derivative of solution vector. */
/*     NEQ       -- Number of unknowns. */
/*     RES       -- External user-supplied subroutine */
/*                  to evaluate the residual.  See RES description */
/*                  in DDASPK prologue. */
/*     JACD      -- External user-supplied routine to evaluate the */
/*                  Jacobian.  See JAC description for the case */
/*                  INFO(12) = 0 in the DDASPK prologue. */
/*     PDUM      -- Dummy argument. */
/*     H         -- Appropriate step size for next step. */
/*     WT        -- Vector of weights for error criterion. */
/*     JSTART    -- Indicates first call to this routine. */
/*                  If JSTART = 0, then this is the first call, */
/*                  otherwise it is not. */
/*     IDID      -- Completion flag, output by DNEDD. */
/*                  See IDID description in DDASPK prologue. */
/*     RPAR,IPAR -- Real and integer arrays used for communication */
/*                  between the calling program and external user */
/*                  routines.  They are not altered within DASPK. */
/*     PHI       -- Array of divided differences used by */
/*                  DNEDD.  The length is NEQ*(K+1),where */
/*                  K is the maximum order. */
/*     GAMMA     -- Array used to predict Y and YPRIME.  The length */
/*                  is MAXORD+1 where MAXORD is the maximum order. */
/*     DUMSVR    -- Dummy argument. */
/*     DELTA     -- Work vector for NLS of length NEQ. */
/*     E         -- Error accumulation vector for NLS of length NEQ. */
/*     WM,IWM    -- Real and integer arrays storing */
/*                  matrix information such as the matrix */
/*                  of partial derivatives, permutation */
/*                  vector, and various other information. */
/*     CJ        -- Parameter always proportional to 1/H. */
/*     CJOLD     -- Saves the value of CJ as of the last call to DMATD. */
/*                  Accounts for changes in CJ needed to */
/*                  decide whether to call DMATD. */
/*     CJLAST    -- Previous value of CJ. */
/*     S         -- A scalar determined by the approximate rate */
/*                  of convergence of the Newton iteration and used */
/*                  in the convergence test for the Newton iteration. */

/*                  If RATE is defined to be an estimate of the */
/*                  rate of convergence of the Newton iteration, */
/*                  then S = RATE/(1.D0-RATE). */

/*                  The closer RATE is to 0., the faster the Newton */
/*                  iteration is converging; the closer RATE is to 1., */
/*                  the slower the Newton iteration is converging. */

/*                  On the first Newton iteration with an up-dated */
/*                  preconditioner S = 100.D0, Thus the initial */
/*                  RATE of convergence is approximately 1. */

/*                  S is preserved from call to call so that the rate */
/*                  estimate from a previous step can be applied to */
/*                  the current step. */
/*     UROUND    -- Unit roundoff. */
/*     DUME      -- Dummy argument. */
/*     DUMS      -- Dummy argument. */
/*     DUMR      -- Dummy argument. */
/*     EPCON     -- Tolerance to test for convergence of the Newton */
/*                  iteration. */
/*     JCALC     -- Flag used to determine when to update */
/*                  the Jacobian matrix.  In general: */

/*                  JCALC = -1 ==> Call the DMATD routine to update */
/*                                 the Jacobian matrix. */
/*                  JCALC =  0 ==> Jacobian matrix is up-to-date. */
/*                  JCALC =  1 ==> Jacobian matrix is out-dated, */
/*                                 but DMATD will not be called unless */
/*                                 JCALC is set to -1. */
/*     JFDUM     -- Dummy argument. */
/*     KP1       -- The current order(K) + 1;  updated across calls. */
/*     NONNEG    -- Flag to determine nonnegativity constraints. */
/*     NTYPE     -- Identification code for the NLS routine. */
/*                   0  ==> modified Newton; direct solver. */
/*     IERNLS    -- Error flag for nonlinear solver. */
/*                   0  ==> nonlinear solver converged. */
/*                   1  ==> recoverable error inside nonlinear solver. */
/*                  -1  ==> unrecoverable error inside nonlinear solver. */

/*     All variables with "DUM" in their names are dummy variables */
/*     which are not used in this routine. */

/*     Following is a list and description of local variables which */
/*     may not have an obvious usage.  They are listed in roughly the */
/*     order they occur in this subroutine. */

/*     The following group of variables are passed as arguments to */
/*     the Newton iteration solver.  They are explained in greater detail */
/*     in DNSD: */
/*        TOLNEW, MULDEL, MAXIT, IERNEW */

/*     IERTYP -- Flag which tells whether this subroutine is correct. */
/*               0 ==> correct subroutine. */
/*               1 ==> incorrect subroutine. */

/* ----------------------------------------------------------------------- */
/* ***ROUTINES CALLED */
/*   DDWNRM, RES, DMATD, DNSD */

/* ***END PROLOGUE  DNEDD */




    /* Parameter adjustments */
    --y;
    --yprime;
    phi_dim1 = *neq;
    phi_offset = 1 + phi_dim1;
    phi -= phi_offset;
    --wt;
    --rpar;
    --ipar;
    --gamma;
    --delta;
    --e;
    --wm;
    --iwm;

    /* Function Body */

/*     Verify that this is the correct subroutine. */

    iertyp = 0;
    if (*ntype != 0) {
	iertyp = 1;
	goto L380;
    }

/*     If this is the first step, perform initializations. */

    if (*jstart == 0) {
	*cjold = *cj;
	*jcalc = -1;
    }

/*     Perform all other initializations. */

    *iernls = 0;

/*     Decide whether new Jacobian is needed. */

    temp1 = (1. - xrate) / (xrate + 1.);
    temp2 = 1. / temp1;
    if (*cj / *cjold < temp1 || *cj / *cjold > temp2) {
	*jcalc = -1;
    }
    if (*cj != *cjlast) {
	*s = 100.;
    }

/* ----------------------------------------------------------------------- */
/*     Entry point for updating the Jacobian with current */
/*     stepsize. */
/* ----------------------------------------------------------------------- */
L300:

/*     Initialize all error flags to zero. */

    ierj = 0;
    ires = 0;
    iernew = 0;

/*     Predict the solution and derivative and compute the tolerance */
/*     for the Newton iteration. */

    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] = phi[i__ + phi_dim1];
/* L310: */
	yprime[i__] = 0.;
    }
    i__1 = *kp1;
    for (j = 2; j <= i__1; ++j) {
	i__2 = *neq;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y[i__] += phi[i__ + j * phi_dim1];
/* L320: */
	    yprime[i__] += gamma[j] * phi[i__ + j * phi_dim1];
	}
/* L330: */
    }
    pnorm = _daskr_ddwnrm_(neq, &y[1], &wt[1], &rpar[1], &ipar[1]);
    tolnew = *uround * 100. * pnorm;

/*     Call RES to initialize DELTA. */

    ++iwm[12];
    ((int (*)(real_number *, real_number *, real_number*, real_number*, real_number *, integer *, real_number *, integer *))*res)
      (x, &y[1], &yprime[1], cj, &delta[1], &ires, &rpar[1], &ipar[1]);
    if (ires < 0) {
	goto L380;
    }

/*     If indicated, reevaluate the iteration matrix */
/*     J = dG/dY + CJ*dG/dYPRIME (where G(X,Y,YPRIME)=0). */
/*     Set JCALC to 0 as an indicator that this has been done. */

    if (*jcalc == -1) {
	++iwm[13];
	*jcalc = 0;
	_daskr_dmatd_(neq, x, &y[1], &yprime[1], &delta[1], cj, h__, &ierj, &wt[1], &
		e[1], &wm[1], &iwm[1], (Unknown_fp)res, &ires, uround, (Unknown_fp)jacd, &
		rpar[1], &ipar[1]);
	*cjold = *cj;
	*s = 100.;
	if (ires < 0) {
	    goto L380;
	}
	if (ierj != 0) {
	    goto L380;
	}
    }

/*     Call the nonlinear Newton solver. */

    temp1 = 2. / (*cj / *cjold + 1.);
    _daskr_dnsd_(x, &y[1], &yprime[1], neq, (Unknown_fp)res, pdum, &wt[1], &rpar[1], &ipar[
	    1], dumsvr, &delta[1], &e[1], &wm[1], &iwm[1], cj, dums, dumr, 
	    dume, epcon, s, &temp1, &tolnew, &muldel, &maxit, &ires, &idum, &
	    iernew);

    if (iernew > 0 && *jcalc != 0) {

/*        The Newton iteration had a recoverable failure with an old */
/*        iteration matrix.  Retry the step with a new iteration matrix. */

	*jcalc = -1;
	goto L300;
    }

    if (iernew != 0) {
	goto L380;
    }

/*     The Newton iteration has converged.  If nonnegativity of */
/*     solution is required, set the solution nonnegative, if the */
/*     perturbation to do it is small enough.  If the change is too */
/*     large, then consider the corrector iteration to have failed. */

/* L375: */
    if (*nonneg == 0) {
	goto L390;
    }
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L377: */
/* Computing MIN */
	d__1 = y[i__];
	delta[i__] = MIN(d__1,0.);
    }
    delnrm = _daskr_ddwnrm_(neq, &delta[1], &wt[1], &rpar[1], &ipar[1]);
    if (delnrm > *epcon) {
	goto L380;
    }
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L378: */
	e[i__] -= delta[i__];
    }
    goto L390;


/*     Exits from nonlinear solver. */
/*     No convergence with current iteration */
/*     matrix, or singular iteration matrix. */
/*     Compute IERNLS and IDID accordingly. */

L380:
    if (ires <= -2 || iertyp != 0) {
	*iernls = -1;
	if (ires <= -2) {
	    *idid = -11;
	}
	if (iertyp != 0) {
	    *idid = -15;
	}
    } else {
	*iernls = 1;
	if (ires < 0) {
	    *idid = -10;
	}
	if (ierj != 0) {
	    *idid = -8;
	}
    }

L390:
    *jcalc = 1;
    return 0;

/* ------END OF SUBROUTINE DNEDD------------------------------------------ */
} /* dnedd_ */

/* Subroutine */ int _daskr_dnsd_(real_number *x, real_number *y, real_number *yprime,
	integer *neq, Unknown_fp res, real_number *pdum, real_number *wt, real_number *
	rpar, integer *ipar, real_number *dumsvr, real_number *delta,
	real_number *e, real_number *wm, integer *iwm, real_number *cj,
	real_number *dums, real_number *dumr, real_number *dume, real_number *
	epcon, real_number *s, real_number *confac, real_number *tolnew, integer
	*muldel, integer *maxit, integer *ires, integer *idum, integer *
	iernew)
{
    /* System generated locals */
    integer i__1;
    real_number d__1, d__2;

    /* Builtin functions */
    double _daskr_real_pow(real_number *, real_number *);

    /* Local variables */
    static integer i__, m;
    static real_number rate;
    extern /* Subroutine */ int _daskr_dslvd_(integer *, real_number *, real_number *,
	    integer *);
    static real_number delnrm;
    extern real_number _daskr_ddwnrm_(integer *, real_number *, real_number *,
	    real_number *, integer *);
    static real_number oldnrm;


/* ***BEGIN PROLOGUE  DNSD */
/* ***REFER TO  DDASPK */
/* ***DATE WRITTEN   891219   (YYMMDD) */
/* ***REVISION DATE  900926   (YYMMDD) */
/* ***REVISION DATE  950126   (YYMMDD) */
/* ***REVISION DATE  000711   (YYMMDD) */


/* ----------------------------------------------------------------------- */
/* ***DESCRIPTION */

/*     DNSD solves a nonlinear system of */
/*     algebraic equations of the form */
/*     G(X,Y,YPRIME) = 0 for the unknown Y. */

/*     The method used is a modified Newton scheme. */

/*     The parameters represent */

/*     X         -- Independent variable. */
/*     Y         -- Solution vector. */
/*     YPRIME    -- Derivative of solution vector. */
/*     NEQ       -- Number of unknowns. */
/*     RES       -- External user-supplied subroutine */
/*                  to evaluate the residual.  See RES description */
/*                  in DDASPK prologue. */
/*     PDUM      -- Dummy argument. */
/*     WT        -- Vector of weights for error criterion. */
/*     RPAR,IPAR -- Real and integer arrays used for communication */
/*                  between the calling program and external user */
/*                  routines.  They are not altered within DASPK. */
/*     DUMSVR    -- Dummy argument. */
/*     DELTA     -- Work vector for DNSD of length NEQ. */
/*     E         -- Error accumulation vector for DNSD of length NEQ. */
/*     WM,IWM    -- Real and integer arrays storing */
/*                  matrix information such as the matrix */
/*                  of partial derivatives, permutation */
/*                  vector, and various other information. */
/*     CJ        -- Parameter always proportional to 1/H (step size). */
/*     DUMS      -- Dummy argument. */
/*     DUMR      -- Dummy argument. */
/*     DUME      -- Dummy argument. */
/*     EPCON     -- Tolerance to test for convergence of the Newton */
/*                  iteration. */
/*     S         -- Used for error convergence tests. */
/*                  In the Newton iteration: S = RATE/(1 - RATE), */
/*                  where RATE is the estimated rate of convergence */
/*                  of the Newton iteration. */
/*                  The calling routine passes the initial value */
/*                  of S to the Newton iteration. */
/*     CONFAC    -- A residual scale factor to improve convergence. */
/*     TOLNEW    -- Tolerance on the norm of Newton correction in */
/*                  alternative Newton convergence test. */
/*     MULDEL    -- A flag indicating whether or not to multiply */
/*                  DELTA by CONFAC. */
/*                  0  ==> do not scale DELTA by CONFAC. */
/*                  1  ==> scale DELTA by CONFAC. */
/*     MAXIT     -- Maximum allowed number of Newton iterations. */
/*     IRES      -- Error flag returned from RES.  See RES description */
/*                  in DDASPK prologue.  If IRES = -1, then IERNEW */
/*                  will be set to 1. */
/*                  If IRES < -1, then IERNEW will be set to -1. */
/*     IDUM      -- Dummy argument. */
/*     IERNEW    -- Error flag for Newton iteration. */
/*                   0  ==> Newton iteration converged. */
/*                   1  ==> recoverable error inside Newton iteration. */
/*                  -1  ==> unrecoverable error inside Newton iteration. */

/*     All arguments with "DUM" in their names are dummy arguments */
/*     which are not used in this routine. */
/* ----------------------------------------------------------------------- */

/* ***ROUTINES CALLED */
/*   DSLVD, DDWNRM, RES */

/* ***END PROLOGUE  DNSD */




/*     Initialize Newton counter M and accumulation vector E. */

    /* Parameter adjustments */
    --iwm;
    --wm;
    --e;
    --delta;
    --ipar;
    --rpar;
    --wt;
    --yprime;
    --y;

    /* Function Body */
    m = 0;
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L100: */
	e[i__] = 0.;
    }

/*     Corrector loop. */

L300:
    ++iwm[19];

/*     If necessary, multiply residual by convergence factor. */

    if (*muldel == 1) {
	i__1 = *neq;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L320: */
	    delta[i__] *= *confac;
	}
    }

/*     Compute a new iterate (back-substitution). */
/*     Store the correction in DELTA. */

    _daskr_dslvd_(neq, &delta[1], &wm[1], &iwm[1]);

/*     Update Y, E, and YPRIME. */

    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] -= delta[i__];
	e[i__] -= delta[i__];
/* L340: */
	yprime[i__] -= *cj * delta[i__];
    }

/*     Test for convergence of the iteration. */

    delnrm = _daskr_ddwnrm_(neq, &delta[1], &wt[1], &rpar[1], &ipar[1]);
    if (m == 0) {
	oldnrm = delnrm;
	if (delnrm <= *tolnew) {
	    goto L370;
	}
    } else {
	d__1 = delnrm / oldnrm;
	d__2 = 1. / m;
	rate = _daskr_real_pow(&d__1, &d__2);
	if (rate > .9) {
	    goto L380;
	}
	*s = rate / (1. - rate);
    }
    if (*s * delnrm <= *epcon) {
	goto L370;
    }

/*     The corrector has not yet converged. */
/*     Update M and test whether the */
/*     maximum number of iterations have */
/*     been tried. */

    ++m;
    if (m >= *maxit) {
	goto L380;
    }

/*     Evaluate the residual, */
/*     and go back to do another iteration. */

    ++iwm[12];
    ((int (*)(real_number *, real_number *, real_number*, real_number*, real_number *, integer *, real_number *, integer *))*res)
      (x, &y[1], &yprime[1], cj, &delta[1], ires, &rpar[1], &ipar[1]);
    if (*ires < 0) {
	goto L380;
    }
    goto L300;

/*     The iteration has converged. */

L370:
    return 0;

/*     The iteration has not converged.  Set IERNEW appropriately. */

L380:
    if (*ires <= -2) {
	*iernew = -1;
    } else {
	*iernew = 1;
    }
    return 0;


/* ------END OF SUBROUTINE DNSD------------------------------------------- */
} /* dnsd_ */

/* Subroutine */ int _daskr_dmatd_(integer *neq, real_number *x, real_number *y,
	real_number *yprime, real_number *delta, real_number *cj, real_number *
	h__, integer *ier, real_number *ewt, real_number *e, real_number *wm,
	integer *iwm, Unknown_fp res, integer *ires, real_number *uround, Unknown_fp jacd,
	real_number *rpar, integer *ipar)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    real_number d__1, d__2, d__3, d__4, d__5, d__6;

    /* Builtin functions */
    double sqrt(real_number), _daskr_real_sign(real_number *, real_number *);

    /* Local variables */
    static integer i__, j, k, l, n, i1, i2, ii, mba;
    static real_number del;
    static integer meb1, nrow;
    static real_number squr;
    extern /* Subroutine */ int _daskr_dgbfa_(real_number *, integer *, integer *,
	    integer *, integer *, integer *, integer *), _daskr_dgefa_(real_number *,
	    integer *, integer *, integer *, integer *);
    static integer mband, lenpd, isave, msave;
    static real_number ysave;
    static integer lipvt, mtype, meband;
    static real_number delinv;
    static integer ipsave;
    static real_number ypsave;


/* ***BEGIN PROLOGUE  DMATD */
/* ***REFER TO  DDASPK */
/* ***DATE WRITTEN   890101   (YYMMDD) */
/* ***REVISION DATE  900926   (YYMMDD) */
/* ***REVISION DATE  940701   (new LIPVT) */
/* ***REVISION DATE  060712   (Changed minimum D.Q. increment to 1/EWT(j)) */

/* ----------------------------------------------------------------------- */
/* ***DESCRIPTION */

/*     This routine computes the iteration matrix */
/*     J = dG/dY+CJ*dG/dYPRIME (where G(X,Y,YPRIME)=0). */
/*     Here J is computed by: */
/*       the user-supplied routine JACD if IWM(MTYPE) is 1 or 4, or */
/*       by numerical difference quotients if IWM(MTYPE) is 2 or 5. */

/*     The parameters have the following meanings. */
/*     X        = Independent variable. */
/*     Y        = Array containing predicted values. */
/*     YPRIME   = Array containing predicted derivatives. */
/*     DELTA    = Residual evaluated at (X,Y,YPRIME). */
/*                (Used only if IWM(MTYPE)=2 or 5). */
/*     CJ       = Scalar parameter defining iteration matrix. */
/*     H        = Current stepsize in integration. */
/*     IER      = Variable which is .NE. 0 if iteration matrix */
/*                is singular, and 0 otherwise. */
/*     EWT      = Vector of error weights for computing norms. */
/*     E        = Work space (temporary) of length NEQ. */
/*     WM       = Real work space for matrices.  On output */
/*                it contains the LU decomposition */
/*                of the iteration matrix. */
/*     IWM      = Integer work space containing */
/*                matrix information. */
/*     RES      = External user-supplied subroutine */
/*                to evaluate the residual.  See RES description */
/*                in DDASPK prologue. */
/*     IRES     = Flag which is equal to zero if no illegal values */
/*                in RES, and less than zero otherwise.  (If IRES */
/*                is less than zero, the matrix was not completed). */
/*                In this case (if IRES .LT. 0), then IER = 0. */
/*     UROUND   = The unit roundoff error of the machine being used. */
/*     JACD     = Name of the external user-supplied routine */
/*                to evaluate the iteration matrix.  (This routine */
/*                is only used if IWM(MTYPE) is 1 or 4) */
/*                See JAC description for the case INFO(12) = 0 */
/*                in DDASPK prologue. */
/*     RPAR,IPAR= Real and integer parameter arrays that */
/*                are used for communication between the */
/*                calling program and external user routines. */
/*                They are not altered by DMATD. */
/* ----------------------------------------------------------------------- */
/* ***ROUTINES CALLED */
/*   JACD, RES, DGEFA, DGBFA */

/* ***END PROLOGUE  DMATD */




    /* Parameter adjustments */
    --ipar;
    --rpar;
    --iwm;
    --wm;
    --e;
    --ewt;
    --delta;
    --yprime;
    --y;

    /* Function Body */
    lipvt = iwm[30];
    *ier = 0;
    mtype = iwm[4];
    switch (mtype) {
	case 1:  goto L100;
	case 2:  goto L200;
	case 3:  goto L300;
	case 4:  goto L400;
	case 5:  goto L500;
    }


/*     Dense user-supplied matrix. */

L100:
    lenpd = iwm[22];
    i__1 = lenpd;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L110: */
	wm[i__] = 0.;
    }
    ((int (*)(real_number *, real_number *, real_number*, real_number*, real_number *, real_number*, real_number *, real_number *, real_number *, integer *))*jacd)
      (x, &y[1], &yprime[1], &delta[1], &wm[1], cj, h__, &ewt[1], &rpar[1], &ipar[1]);
    goto L230;


/*     Dense finite-difference-generated matrix. */

L200:
    *ires = 0;
    nrow = 0;
    squr = sqrt(*uround);
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
/* Computing MAX */
	d__5 = (d__1 = y[i__], fabs(d__1)), d__6 = (d__2 = *h__ * yprime[i__],
		fabs(d__2));
	d__3 = squr * MAX(d__5,d__6), d__4 = 1. / ewt[i__];
	del = MAX(d__3,d__4);
	d__1 = *h__ * yprime[i__];
	del = _daskr_real_sign(&del, &d__1);
	del = y[i__] + del - y[i__];
	ysave = y[i__];
	ypsave = yprime[i__];
	y[i__] += del;
	yprime[i__] += *cj * del;
	++iwm[12];
	((int (*)(real_number *, real_number *, real_number*, real_number*, real_number *, integer *, real_number *, integer *))*res)
    (x, &y[1], &yprime[1], cj, &e[1], ires, &rpar[1], &ipar[1]);
	if (*ires < 0) {
	    return 0;
	}
	delinv = 1. / del;
	i__2 = *neq;
	for (l = 1; l <= i__2; ++l) {
/* L220: */
	    wm[nrow + l] = (e[l] - delta[l]) * delinv;
	}
	nrow += *neq;
	y[i__] = ysave;
	yprime[i__] = ypsave;
/* L210: */
    }


/*     Do dense-matrix LU decomposition on J. */

L230:
    _daskr_dgefa_(&wm[1], neq, neq, &iwm[lipvt], ier);
    return 0;


/*     Dummy section for IWM(MTYPE)=3. */

L300:
    return 0;


/*     Banded user-supplied matrix. */

L400:
    lenpd = iwm[22];
    i__1 = lenpd;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L410: */
	wm[i__] = 0.;
    }
    ((int (*)(real_number *, real_number *, real_number*, real_number*, real_number *, real_number*, real_number *, real_number *, real_number *, integer *))*jacd)
      (x, &y[1], &yprime[1], &delta[1], &wm[1], cj, h__, &ewt[1], &rpar[1], &ipar[1]);
    meband = (iwm[1] << 1) + iwm[2] + 1;
    goto L550;


/*     Banded finite-difference-generated matrix. */

L500:
    mband = iwm[1] + iwm[2] + 1;
    mba = MIN(mband,*neq);
    meband = mband + iwm[1];
    meb1 = meband - 1;
    msave = *neq / mband + 1;
    isave = iwm[22];
    ipsave = isave + msave;
    *ires = 0;
    squr = sqrt(*uround);
    i__1 = mba;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *neq;
	i__3 = mband;
	for (n = j; i__3 < 0 ? n >= i__2 : n <= i__2; n += i__3) {
	    k = (n - j) / mband + 1;
	    wm[isave + k] = y[n];
	    wm[ipsave + k] = yprime[n];
/* Computing MAX */
/* Computing MAX */
	    d__5 = (d__1 = y[n], fabs(d__1)), d__6 = (d__2 = *h__ * yprime[n],
		    fabs(d__2));
	    d__3 = squr * MAX(d__5,d__6), d__4 = 1. / ewt[n];
	    del = MAX(d__3,d__4);
	    d__1 = *h__ * yprime[n];
	    del = _daskr_real_sign(&del, &d__1);
	    del = y[n] + del - y[n];
	    y[n] += del;
/* L510: */
	    yprime[n] += *cj * del;
	}
	++iwm[12];
	((int (*)(real_number *, real_number *, real_number*, real_number*, real_number *, integer *, real_number *, integer *))*res)
    (x, &y[1], &yprime[1], cj, &e[1], ires, &rpar[1], &ipar[1]);
	if (*ires < 0) {
	    return 0;
	}
	i__3 = *neq;
	i__2 = mband;
	for (n = j; i__2 < 0 ? n >= i__3 : n <= i__3; n += i__2) {
	    k = (n - j) / mband + 1;
	    y[n] = wm[isave + k];
	    yprime[n] = wm[ipsave + k];
/* Computing MAX */
/* Computing MAX */
	    d__5 = (d__1 = y[n], fabs(d__1)), d__6 = (d__2 = *h__ * yprime[n],
		    fabs(d__2));
	    d__3 = squr * MAX(d__5,d__6), d__4 = 1. / ewt[n];
	    del = MAX(d__3,d__4);
	    d__1 = *h__ * yprime[n];
	    del = _daskr_real_sign(&del, &d__1);
	    del = y[n] + del - y[n];
	    delinv = 1. / del;
/* Computing MAX */
	    i__4 = 1, i__5 = n - iwm[2];
	    i1 = MAX(i__4,i__5);
/* Computing MIN */
	    i__4 = *neq, i__5 = n + iwm[1];
	    i2 = MIN(i__4,i__5);
	    ii = n * meb1 - iwm[1];
	    i__4 = i2;
	    for (i__ = i1; i__ <= i__4; ++i__) {
/* L520: */
		wm[ii + i__] = (e[i__] - delta[i__]) * delinv;
	    }
/* L530: */
	}
/* L540: */
    }


/*     Do LU decomposition of banded J. */

L550:
    _daskr_dgbfa_(&wm[1], &meband, neq, &iwm[1], &iwm[2], &iwm[lipvt], ier);
    return 0;

/* ------END OF SUBROUTINE DMATD------------------------------------------ */
} /* dmatd_ */

/* Subroutine */ int _daskr_dslvd_(integer *neq, real_number *delta, real_number *wm,
	integer *iwm)
{
    extern /* Subroutine */ int _daskr_dgbsl_(real_number *, integer *, integer *,
	    integer *, integer *, integer *, real_number *, integer *), _daskr_dgesl_(
	    real_number *, integer *, integer *, integer *, real_number *,
	    integer *);
    static integer lipvt, mtype, meband;


/* ***BEGIN PROLOGUE  DSLVD */
/* ***REFER TO  DDASPK */
/* ***DATE WRITTEN   890101   (YYMMDD) */
/* ***REVISION DATE  900926   (YYMMDD) */
/* ***REVISION DATE  940701   (YYMMDD) (new LIPVT) */

/* ----------------------------------------------------------------------- */
/* ***DESCRIPTION */

/*     This routine manages the solution of the linear */
/*     system arising in the Newton iteration. */
/*     Real matrix information and real temporary storage */
/*     is stored in the array WM. */
/*     Integer matrix information is stored in the array IWM. */
/*     For a dense matrix, the LINPACK routine DGESL is called. */
/*     For a banded matrix, the LINPACK routine DGBSL is called. */
/* ----------------------------------------------------------------------- */
/* ***ROUTINES CALLED */
/*   DGESL, DGBSL */

/* ***END PROLOGUE  DSLVD */




    /* Parameter adjustments */
    --iwm;
    --wm;
    --delta;

    /* Function Body */
    lipvt = iwm[30];
    mtype = iwm[4];
    switch (mtype) {
	case 1:  goto L100;
	case 2:  goto L100;
	case 3:  goto L300;
	case 4:  goto L400;
	case 5:  goto L400;
    }

/*     Dense matrix. */

L100:
    _daskr_dgesl_(&wm[1], neq, neq, &iwm[lipvt], &delta[1], &c__0);
    return 0;

/*     Dummy section for MTYPE=3. */

L300:
    return 0;

/*     Banded matrix. */

L400:
    meband = (iwm[1] << 1) + iwm[2] + 1;
    _daskr_dgbsl_(&wm[1], &meband, neq, &iwm[1], &iwm[2], &iwm[lipvt], &delta[1], &
	    c__0);
    return 0;

/* ------END OF SUBROUTINE DSLVD------------------------------------------ */
} /* dslvd_ */

/* Subroutine */ int _daskr_ddasik_(real_number *x, real_number *y, real_number *yprime,
	 integer *neq, integer *icopt, integer *id, Unknown_fp res, Unknown_fp jack, Unknown_fp
	psol, real_number *h__, real_number *tscale, real_number *wt, integer *
	jskip, real_number *rpar, integer *ipar, real_number *savr, real_number *
	delta, real_number *r__, real_number *yic, real_number *ypic, real_number
	*pwk, real_number *wm, integer *iwm, real_number *cj, real_number *
	uround, real_number *epli, real_number *sqrtn, real_number *rsqrtn,
	real_number *epcon, real_number *ratemx, real_number *stptol, integer *
	jflg, integer *icnflg, integer *icnstr, integer *iernls)
{
    static integer nj, lwp, ires, liwp, mxnj;
    static real_number eplin;
    extern /* Subroutine */ int _daskr_dnsik_(real_number *, real_number *, real_number
	    *, integer *, integer *, integer *, Unknown_fp, Unknown_fp, real_number *,
	    real_number *, integer *, real_number *, real_number *, real_number *,
	     real_number *, real_number *, real_number *, real_number *, integer *
	    , real_number *, real_number *, real_number *, real_number *,
	    real_number *, real_number *, real_number *, integer *, real_number *,
	     integer *, integer *, integer *);
    static integer ierpj;
    extern /* Subroutine */ int _daskr_dcopy_(integer *, real_number *, integer *,
	    real_number *, integer *);
    static integer mxnit, iernew;


/* ***BEGIN PROLOGUE  DDASIK */
/* ***REFER TO  DDASPK */
/* ***DATE WRITTEN   941026   (YYMMDD) */
/* ***REVISION DATE  950808   (YYMMDD) */
/* ***REVISION DATE  951110   Removed unreachable block 390. */
/* ***REVISION DATE  000628   TSCALE argument added. */


/* ----------------------------------------------------------------------- */
/* ***DESCRIPTION */


/*     DDASIK solves a nonlinear system of algebraic equations of the */
/*     form G(X,Y,YPRIME) = 0 for the unknown parts of Y and YPRIME in */
/*     the initial conditions. */

/*     An initial value for Y and initial guess for YPRIME are input. */

/*     The method used is a Newton scheme with Krylov iteration and a */
/*     linesearch algorithm. */

/*     The parameters represent */

/*     X         -- Independent variable. */
/*     Y         -- Solution vector at x. */
/*     YPRIME    -- Derivative of solution vector. */
/*     NEQ       -- Number of equations to be integrated. */
/*     ICOPT     -- Initial condition option chosen (1 or 2). */
/*     ID        -- Array of dimension NEQ, which must be initialized */
/*                  if ICOPT = 1.  See DDASIC. */
/*     RES       -- External user-supplied subroutine */
/*                  to evaluate the residual.  See RES description */
/*                  in DDASPK prologue. */
/*     JACK     --  External user-supplied routine to update */
/*                  the preconditioner.  (This is optional). */
/*                  See JAC description for the case */
/*                  INFO(12) = 1 in the DDASPK prologue. */
/*     PSOL      -- External user-supplied routine to solve */
/*                  a linear system using preconditioning. */
/*                  (This is optional).  See explanation inside DDASPK. */
/*     H         -- Scaling factor for this initial condition calc. */
/*     TSCALE    -- Scale factor in T, used for stopping tests if nonzero. */
/*     WT        -- Vector of weights for error criterion. */
/*     JSKIP     -- input flag to signal if initial JAC call is to be */
/*                  skipped.  1 => skip the call, 0 => do not skip call. */
/*     RPAR,IPAR -- Real and integer arrays used for communication */
/*                  between the calling program and external user */
/*                  routines.  They are not altered within DASPK. */
/*     SAVR      -- Work vector for DDASIK of length NEQ. */
/*     DELTA     -- Work vector for DDASIK of length NEQ. */
/*     R         -- Work vector for DDASIK of length NEQ. */
/*     YIC,YPIC  -- Work vectors for DDASIK, each of length NEQ. */
/*     PWK       -- Work vector for DDASIK of length NEQ. */
/*     WM,IWM    -- Real and integer arrays storing */
/*                  matrix information for linear system */
/*                  solvers, and various other information. */
/*     CJ        -- Matrix parameter = 1/H (ICOPT = 1) or 0 (ICOPT = 2). */
/*     UROUND    -- Unit roundoff.  Not used here. */
/*     EPLI      -- convergence test constant. */
/*                  See DDASPK prologue for more details. */
/*     SQRTN     -- Square root of NEQ. */
/*     RSQRTN    -- reciprical of square root of NEQ. */
/*     EPCON     -- Tolerance to test for convergence of the Newton */
/*                  iteration. */
/*     RATEMX    -- Maximum convergence rate for which Newton iteration */
/*                  is considered converging. */
/*     JFLG      -- Flag showing whether a Jacobian routine is supplied. */
/*     ICNFLG    -- Integer scalar.  If nonzero, then constraint */
/*                  violations in the proposed new approximate solution */
/*                  will be checked for, and the maximum step length */
/*                  will be adjusted accordingly. */
/*     ICNSTR    -- Integer array of length NEQ containing flags for */
/*                  checking constraints. */
/*     IERNLS    -- Error flag for nonlinear solver. */
/*                   0   ==> nonlinear solver converged. */
/*                   1,2 ==> recoverable error inside nonlinear solver. */
/*                           1 => retry with current Y, YPRIME */
/*                           2 => retry with original Y, YPRIME */
/*                  -1   ==> unrecoverable error in nonlinear solver. */

/* ----------------------------------------------------------------------- */

/* ***ROUTINES CALLED */
/*   RES, JACK, DNSIK, DCOPY */

/* ***END PROLOGUE  DDASIK */





/*     Perform initializations. */

    /* Parameter adjustments */
    --icnstr;
    --iwm;
    --wm;
    --pwk;
    --ypic;
    --yic;
    --r__;
    --delta;
    --savr;
    --ipar;
    --rpar;
    --wt;
    --id;
    --yprime;
    --y;

    /* Function Body */
    lwp = iwm[29];
    liwp = iwm[30];
    mxnit = iwm[32];
    mxnj = iwm[33];
    *iernls = 0;
    nj = 0;
    eplin = *epli * *epcon;

/*     Call RES to initialize DELTA. */

    ires = 0;
    ++iwm[12];
    ((int (*)(real_number *, real_number *, real_number*, real_number*, real_number *, integer *, real_number *, integer *))*res)
      (x, &y[1], &yprime[1], cj, &delta[1], &ires, &rpar[1], &ipar[1]);
    if (ires < 0) {
	goto L370;
    }

/*     Looping point for updating the preconditioner. */

L300:

/*     Initialize all error flags to zero. */

    ierpj = 0;
    ires = 0;
    iernew = 0;

/*     If a Jacobian routine was supplied, call it. */

    if (*jflg == 1 && *jskip == 0) {
	++nj;
	++iwm[13];
	((int (*)(Unknown_fp, integer *, integer *, real_number*, real_number *, real_number *, real_number *, real_number *,real_number *, real_number *, real_number*, real_number*, integer *, integer *, real_number *, integer *))*jack)
    ((Unknown_fp)res, &ires, neq, x, &y[1], &yprime[1], &wt[1], &delta[1]
		, &r__[1], h__, cj, &wm[lwp], &iwm[liwp], &ierpj, &rpar[1], &
		ipar[1]);
	if (ires < 0 || ierpj != 0) {
	    goto L370;
	}
    }
    *jskip = 0;

/*     Call the nonlinear Newton solver for up to MXNIT iterations. */

    _daskr_dnsik_(x, &y[1], &yprime[1], neq, icopt, &id[1], (Unknown_fp)res, (Unknown_fp)psol, &
	    wt[1], &rpar[1], &ipar[1], &savr[1], &delta[1], &r__[1], &yic[1], 
	    &ypic[1], &pwk[1], &wm[1], &iwm[1], cj, tscale, sqrtn, rsqrtn, &
	    eplin, epcon, ratemx, &mxnit, stptol, icnflg, &icnstr[1], &iernew)
	    ;

    if (iernew == 1 && nj < mxnj && *jflg == 1) {

/*       Up to MXNIT iterations were done, the convergence rate is < 1, */
/*       a Jacobian routine is supplied, and the number of JACK calls */
/*       is less than MXNJ. */
/*       Copy the residual SAVR to DELTA, call JACK, and try again. */

	_daskr_dcopy_(neq, &savr[1], &c__1, &delta[1], &c__1);
	goto L300;
    }

    if (iernew != 0) {
	goto L380;
    }
    return 0;


/*     Unsuccessful exits from nonlinear solver. */
/*     Set IERNLS accordingly. */

L370:
    *iernls = 2;
    if (ires <= -2) {
	*iernls = -1;
    }
    return 0;

L380:
    *iernls = MIN(iernew,2);
    return 0;

/* ----------------------- END OF SUBROUTINE DDASIK----------------------- */
} /* ddasik_ */

/* Subroutine */ int _daskr_dnsik_(real_number *x, real_number *y, real_number *yprime,
	integer *neq, integer *icopt, integer *id, Unknown_fp res, Unknown_fp psol,
	real_number *wt, real_number *rpar, integer *ipar, real_number *savr,
	real_number *delta, real_number *r__, real_number *yic, real_number *ypic,
	 real_number *pwk, real_number *wm, integer *iwm, real_number *cj,
	real_number *tscale, real_number *sqrtn, real_number *rsqrtn, real_number
	*eplin, real_number *epcon, real_number *ratemx, integer *maxit,
	real_number *stptol, integer *icnflg, integer *icnstr, integer *iernew)
{
    static integer m, ier, lwp;
    static real_number rlx, rate;
    static integer ires;
    static real_number fnrm, rhok;
    static integer iret, liwp;
    static real_number fnrm0;
    static integer lsoff;
    extern /* Subroutine */ int _daskr_dcopy_(integer *, real_number *, integer *,
	    real_number *, integer *);
    static integer iersl;
    extern /* Subroutine */ int _daskr_dslvk_(integer *, real_number *, real_number *,
	    real_number *, real_number *, real_number *, real_number *,
	    real_number *, integer *, Unknown_fp, integer *, Unknown_fp, integer *,
	    real_number *, real_number *, real_number *, real_number *,
	    real_number *, real_number *, integer *);
    static real_number oldfnm;
    extern /* Subroutine */ int _daskr_dfnrmk_(integer *, real_number *, real_number *,
	     real_number *, real_number *, real_number *, real_number *,
	    real_number *, real_number *, real_number *, real_number *, Unknown_fp,
	    integer *, Unknown_fp, integer *, integer *, real_number *, real_number *,
	     real_number *, integer *, real_number *, real_number *, integer *);
    static real_number delnrm;
    extern /* Subroutine */ int _daskr_dlinsk_(integer *, real_number *, real_number *,
	     real_number *, real_number *, real_number *, real_number *,
	    real_number *, real_number *, real_number *, real_number *,
	    real_number *, integer *, real_number *, integer *, Unknown_fp, integer *,
	     Unknown_fp, real_number *, integer *, real_number *, real_number *,
	    integer *, integer *, real_number *, integer *, real_number *,
	    real_number *, real_number *, real_number *, real_number *, integer *,
	     integer *, real_number *, real_number *, integer *);
    extern real_number _daskr_ddwnrm_(integer *, real_number *, real_number *,
	    real_number *, integer *);


/* ***BEGIN PROLOGUE  DNSIK */
/* ***REFER TO  DDASPK */
/* ***DATE WRITTEN   940701   (YYMMDD) */
/* ***REVISION DATE  950714   (YYMMDD) */
/* ***REVISION DATE  000628   TSCALE argument added. */
/* ***REVISION DATE  000628   Added criterion for IERNEW = 1 return. */


/* ----------------------------------------------------------------------- */
/* ***DESCRIPTION */

/*     DNSIK solves a nonlinear system of algebraic equations of the */
/*     form G(X,Y,YPRIME) = 0 for the unknown parts of Y and YPRIME in */
/*     the initial conditions. */

/*     The method used is a Newton scheme combined with a linesearch */
/*     algorithm, using Krylov iterative linear system methods. */

/*     The parameters represent */

/*     X         -- Independent variable. */
/*     Y         -- Solution vector. */
/*     YPRIME    -- Derivative of solution vector. */
/*     NEQ       -- Number of unknowns. */
/*     ICOPT     -- Initial condition option chosen (1 or 2). */
/*     ID        -- Array of dimension NEQ, which must be initialized */
/*                  if ICOPT = 1.  See DDASIC. */
/*     RES       -- External user-supplied subroutine */
/*                  to evaluate the residual.  See RES description */
/*                  in DDASPK prologue. */
/*     PSOL      -- External user-supplied routine to solve */
/*                  a linear system using preconditioning. */
/*                  See explanation inside DDASPK. */
/*     WT        -- Vector of weights for error criterion. */
/*     RPAR,IPAR -- Real and integer arrays used for communication */
/*                  between the calling program and external user */
/*                  routines.  They are not altered within DASPK. */
/*     SAVR      -- Work vector for DNSIK of length NEQ. */
/*     DELTA     -- Residual vector on entry, and work vector of */
/*                  length NEQ for DNSIK. */
/*     R         -- Work vector for DNSIK of length NEQ. */
/*     YIC,YPIC  -- Work vectors for DNSIK, each of length NEQ. */
/*     PWK       -- Work vector for DNSIK of length NEQ. */
/*     WM,IWM    -- Real and integer arrays storing */
/*                  matrix information such as the matrix */
/*                  of partial derivatives, permutation */
/*                  vector, and various other information. */
/*     CJ        -- Matrix parameter = 1/H (ICOPT = 1) or 0 (ICOPT = 2). */
/*     TSCALE    -- Scale factor in T, used for stopping tests if nonzero. */
/*     SQRTN     -- Square root of NEQ. */
/*     RSQRTN    -- reciprical of square root of NEQ. */
/*     EPLIN     -- Tolerance for linear system solver. */
/*     EPCON     -- Tolerance to test for convergence of the Newton */
/*                  iteration. */
/*     RATEMX    -- Maximum convergence rate for which Newton iteration */
/*                  is considered converging. */
/*     MAXIT     -- Maximum allowed number of Newton iterations. */
/*     STPTOL    -- Tolerance used in calculating the minimum lambda */
/*                  value allowed. */
/*     ICNFLG    -- Integer scalar.  If nonzero, then constraint */
/*                  violations in the proposed new approximate solution */
/*                  will be checked for, and the maximum step length */
/*                  will be adjusted accordingly. */
/*     ICNSTR    -- Integer array of length NEQ containing flags for */
/*                  checking constraints. */
/*     IERNEW    -- Error flag for Newton iteration. */
/*                   0  ==> Newton iteration converged. */
/*                   1  ==> failed to converge, but RATE .lt. 1, or the */
/*                          residual norm was reduced by a factor of .1. */
/*                   2  ==> failed to converge, RATE .gt. RATEMX. */
/*                   3  ==> other recoverable error. */
/*                  -1  ==> unrecoverable error inside Newton iteration. */
/* ----------------------------------------------------------------------- */

/* ***ROUTINES CALLED */
/*   DFNRMK, DSLVK, DDWNRM, DLINSK, DCOPY */

/* ***END PROLOGUE  DNSIK */





/*     Initializations.  M is the Newton iteration counter. */

    /* Parameter adjustments */
    --icnstr;
    --iwm;
    --wm;
    --pwk;
    --ypic;
    --yic;
    --r__;
    --delta;
    --savr;
    --ipar;
    --rpar;
    --wt;
    --id;
    --yprime;
    --y;

    /* Function Body */
    lsoff = iwm[35];
    m = 0;
    rate = 1.;
    lwp = iwm[29];
    liwp = iwm[30];
    rlx = .4;

/*     Save residual in SAVR. */

    _daskr_dcopy_(neq, &delta[1], &c__1, &savr[1], &c__1);

/*     Compute norm of (P-inverse)*(residual). */

    _daskr_dfnrmk_(neq, &y[1], x, &yprime[1], &savr[1], &r__[1], cj, tscale, &wt[1],
	    sqrtn, rsqrtn, (Unknown_fp)res, &ires, (Unknown_fp)psol, &c__1, &ier, &fnrm,
	    eplin, &wm[lwp], &iwm[liwp], &pwk[1], &rpar[1], &ipar[1]);
    ++iwm[21];
    if (ier != 0) {
	*iernew = 3;
	return 0;
    }

/*     Return now if residual norm is .le. EPCON. */

    if (fnrm <= *epcon) {
	return 0;
    }

/*     Newton iteration loop. */

    fnrm0 = fnrm;
L300:
    ++iwm[19];

/*     Compute a new step vector DELTA. */

    _daskr_dslvk_(neq, &y[1], x, &yprime[1], &savr[1], &delta[1], &wt[1], &wm[1], &
	    iwm[1], (Unknown_fp)res, &ires, (Unknown_fp)psol, &iersl, cj, eplin, sqrtn,
	    rsqrtn, &rhok, &rpar[1], &ipar[1]);
    if (ires != 0 || iersl != 0) {
	goto L390;
    }

/*     Get norm of DELTA.  Return now if DELTA is zero. */

    delnrm = _daskr_ddwnrm_(neq, &delta[1], &wt[1], &rpar[1], &ipar[1]);
    if (delnrm == 0.) {
	return 0;
    }

/*     Call linesearch routine for global strategy and set RATE. */

    oldfnm = fnrm;

    _daskr_dlinsk_(neq, &y[1], x, &yprime[1], &savr[1], cj, tscale, &delta[1], &
	    delnrm, &wt[1], sqrtn, rsqrtn, &lsoff, stptol, &iret, (Unknown_fp)res, &
	    ires, (Unknown_fp)psol, &wm[1], &iwm[1], &rhok, &fnrm, icopt, &id[1], &
	    wm[lwp], &iwm[liwp], &r__[1], eplin, &yic[1], &ypic[1], &pwk[1], 
	    icnflg, &icnstr[1], &rlx, &rpar[1], &ipar[1]);

    rate = fnrm / oldfnm;

/*     Check for error condition from linesearch. */
    if (iret != 0) {
	goto L390;
    }

/*     Test for convergence of the iteration, and return or loop. */

    if (fnrm <= *epcon) {
	return 0;
    }

/*     The iteration has not yet converged.  Update M. */
/*     Test whether the maximum number of iterations have been tried. */

    ++m;
    if (m >= *maxit) {
	goto L380;
    }

/*     Copy the residual SAVR to DELTA and loop for another iteration. */

    _daskr_dcopy_(neq, &savr[1], &c__1, &delta[1], &c__1);
    goto L300;

/*     The maximum number of iterations was done.  Set IERNEW and return. */

L380:
    if (rate <= *ratemx || fnrm <= fnrm0 * .1) {
	*iernew = 1;
    } else {
	*iernew = 2;
    }
    return 0;

L390:
    if (ires <= -2 || iersl < 0) {
	*iernew = -1;
    } else {
	*iernew = 3;
	if (ires == 0 && iersl == 1 && m >= 2 && rate < 1.) {
	    *iernew = 1;
	}
    }
    return 0;


/* ----------------------- END OF SUBROUTINE DNSIK------------------------ */
} /* dnsik_ */

/* Subroutine */ int _daskr_dlinsk_(integer *neq, real_number *y, real_number *t,
	real_number *yprime, real_number *savr, real_number *cj, real_number *
	tscale, real_number *p, real_number *pnrm, real_number *wt, real_number *
	sqrtn, real_number *rsqrtn, integer *lsoff, real_number *stptol,
	integer *iret, Unknown_fp res, integer *ires, Unknown_fp psol, real_number *wm,
	integer *iwm, real_number *rhok, real_number *fnrm, integer *icopt,
	integer *id, real_number *wp, integer *iwp, real_number *r__,
	real_number *eplin, real_number *ynew, real_number *ypnew, real_number *
	pwk, integer *icnflg, integer *icnstr, real_number *rlx, real_number *
	rpar, integer *ipar)
{
    /* Initialized data */

    static real_number alpha = 1e-4;
    static real_number one = 1.;
    static real_number two = 2.;

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int _daskr_str_copy(char *, char *, integer, integer);

    /* Local variables */
    static integer i__;
    static real_number rl;
    static integer ier;
    static char msg[80];
    static real_number tau;
    static integer ivar;
    static real_number slpi, f1nrm, ratio;
    extern /* Subroutine */ int _daskr_dcopy_(integer *, real_number *, integer *,
	    real_number *, integer *);
    static real_number rlmin, fnrmp;
    static integer kprin;
    static real_number ratio1, f1nrmp;
    extern /* Subroutine */ int _daskr_dfnrmk_(integer *, real_number *, real_number *,
	     real_number *, real_number *, real_number *, real_number *,
	    real_number *, real_number *, real_number *, real_number *, Unknown_fp,
	    integer *, Unknown_fp, integer *, integer *, real_number *, real_number *,
	     real_number *, integer *, real_number *, real_number *, integer *),
	    _daskr_dcnstr_(integer *, real_number *, real_number *, integer *,
	    real_number *, real_number *, integer *, integer *), _daskr_xerrwd_(char *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, real_number *, real_number *, integer), _daskr_dyypnw_(integer *
	    , real_number *, real_number *, real_number *, real_number *,
	    real_number *, integer *, integer *, real_number *, real_number *);


/* ***BEGIN PROLOGUE  DLINSK */
/* ***REFER TO  DNSIK */
/* ***DATE WRITTEN   940830   (YYMMDD) */
/* ***REVISION DATE  951006   (Arguments SQRTN, RSQRTN added.) */
/* ***REVISION DATE  960129   Moved line RL = ONE to top block. */
/* ***REVISION DATE  000628   TSCALE argument added. */
/* ***REVISION DATE  000628   RHOK*RHOK term removed in alpha test. */


/* ----------------------------------------------------------------------- */
/* ***DESCRIPTION */

/*     DLINSK uses a linesearch algorithm to calculate a new (Y,YPRIME) */
/*     pair (YNEW,YPNEW) such that */

/*     f(YNEW,YPNEW) .le. (1 - 2*ALPHA*RL)*f(Y,YPRIME) */

/*     where 0 < RL <= 1, and RHOK is the scaled preconditioned norm of */
/*     the final residual vector in the Krylov iteration. */
/*     Here, f(y,y') is defined as */

/*      f(y,y') = (1/2)*norm( (P-inverse)*G(t,y,y') )**2 , */

/*     where norm() is the weighted RMS vector norm, G is the DAE */
/*     system residual function, and P is the preconditioner used */
/*     in the Krylov iteration. */

/*     In addition to the parameters defined elsewhere, we have */

/*     SAVR    -- Work array of length NEQ, containing the residual */
/*                vector G(t,y,y') on return. */
/*     TSCALE  -- Scale factor in T, used for stopping tests if nonzero. */
/*     P       -- Approximate Newton step used in backtracking. */
/*     PNRM    -- Weighted RMS norm of P. */
/*     LSOFF   -- Flag showing whether the linesearch algorithm is */
/*                to be invoked.  0 means do the linesearch, */
/*                1 means turn off linesearch. */
/*     STPTOL  -- Tolerance used in calculating the minimum lambda */
/*                value allowed. */
/*     ICNFLG  -- Integer scalar.  If nonzero, then constraint violations */
/*                in the proposed new approximate solution will be */
/*                checked for, and the maximum step length will be */
/*                adjusted accordingly. */
/*     ICNSTR  -- Integer array of length NEQ containing flags for */
/*                checking constraints. */
/*     RHOK    -- Weighted norm of preconditioned Krylov residual. */
/*     RLX     -- Real scalar restricting update size in DCNSTR. */
/*     YNEW    -- Array of length NEQ used to hold the new Y in */
/*                performing the linesearch. */
/*     YPNEW   -- Array of length NEQ used to hold the new YPRIME in */
/*                performing the linesearch. */
/*     PWK     -- Work vector of length NEQ for use in PSOL. */
/*     Y       -- Array of length NEQ containing the new Y (i.e.,=YNEW). */
/*     YPRIME  -- Array of length NEQ containing the new YPRIME */
/*                (i.e.,=YPNEW). */
/*     FNRM    -- Real scalar containing SQRT(2*f(Y,YPRIME)) for the */
/*                current (Y,YPRIME) on input and output. */
/*     R       -- Work space length NEQ for residual vector. */
/*     IRET    -- Return flag. */
/*                IRET=0 means that a satisfactory (Y,YPRIME) was found. */
/*                IRET=1 means that the routine failed to find a new */
/*                       (Y,YPRIME) that was sufficiently distinct from */
/*                       the current (Y,YPRIME) pair. */
/*                IRET=2 means a failure in RES or PSOL. */
/* ----------------------------------------------------------------------- */

/* ***ROUTINES CALLED */
/*   DFNRMK, DYYPNW, DCNSTR, DCOPY, XERRWD */

/* ***END PROLOGUE  DLINSK */



    /* Parameter adjustments */
    --ipar;
    --rpar;
    --icnstr;
    --pwk;
    --ypnew;
    --ynew;
    --r__;
    --iwp;
    --wp;
    --id;
    --iwm;
    --wm;
    --wt;
    --p;
    --savr;
    --yprime;
    --y;

    /* Function Body */

    kprin = iwm[31];
    f1nrm = *fnrm * *fnrm / two;
    ratio = one;

    if (kprin >= 2) {
	_daskr_str_copy(msg, "------ IN ROUTINE DLINSK-- PNRM = (R1)", (integer)80, (
		integer)38);
	_daskr_xerrwd_(msg, &c__38, &c__921, &c__0, &c__0, &c__0, &c__0, &c__1, pnrm,
		 &c_b38, (integer)80);
    }
    tau = *pnrm;
    rl = one;
/* ----------------------------------------------------------------------- */
/* Check for violations of the constraints, if any are imposed. */
/* If any violations are found, the step vector P is rescaled, and the */
/* constraint check is repeated, until no violations are found. */
/* ----------------------------------------------------------------------- */
    if (*icnflg != 0) {
L10:
	_daskr_dyypnw_(neq, &y[1], &yprime[1], cj, &rl, &p[1], icopt, &id[1], &ynew[
		1], &ypnew[1]);
	_daskr_dcnstr_(neq, &y[1], &ynew[1], &icnstr[1], &tau, rlx, iret, &ivar);
	if (*iret == 1) {
	    ratio1 = tau / *pnrm;
	    ratio *= ratio1;
	    i__1 = *neq;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
		p[i__] *= ratio1;
	    }
	    *pnrm = tau;
	    if (kprin >= 2) {
		_daskr_str_copy(msg, "------ CONSTRAINT VIOL., PNRM = (R1), INDEX = ("
			"I1)", (integer)80, (integer)50);
		_daskr_xerrwd_(msg, &c__50, &c__922, &c__0, &c__1, &ivar, &c__0, &
			c__1, pnrm, &c_b38, (integer)80);
	    }
	    if (*pnrm <= *stptol) {
		*iret = 1;
		return 0;
	    }
	    goto L10;
	}
    }

    slpi = -two * f1nrm * ratio;
    rlmin = *stptol / *pnrm;
    if (*lsoff == 0 && kprin >= 2) {
	_daskr_str_copy(msg, "------ MIN. LAMBDA = (R1)", (integer)80, (integer)25);
	_daskr_xerrwd_(msg, &c__25, &c__923, &c__0, &c__0, &c__0, &c__0, &c__1, &
		rlmin, &c_b38, (integer)80);
    }
/* ----------------------------------------------------------------------- */
/* Begin iteration to find RL value satisfying alpha-condition. */
/* Update YNEW and YPNEW, then compute norm of new scaled residual and */
/* perform alpha condition test. */
/* ----------------------------------------------------------------------- */
L100:
    _daskr_dyypnw_(neq, &y[1], &yprime[1], cj, &rl, &p[1], icopt, &id[1], &ynew[1], &
	    ypnew[1]);
    _daskr_dfnrmk_(neq, &ynew[1], t, &ypnew[1], &savr[1], &r__[1], cj, tscale, &wt[1]
	    , sqrtn, rsqrtn, (Unknown_fp)res, ires, (Unknown_fp)psol, &c__0, &ier, &fnrmp,
	     eplin, &wp[1], &iwp[1], &pwk[1], &rpar[1], &ipar[1]);
    ++iwm[12];
    if (*ires >= 0) {
	++iwm[21];
    }
    if (*ires != 0 || ier != 0) {
	*iret = 2;
	return 0;
    }
    if (*lsoff == 1) {
	goto L150;
    }

    f1nrmp = fnrmp * fnrmp / two;
    if (kprin >= 2) {
	_daskr_str_copy(msg, "------ LAMBDA = (R1)", (integer)80, (integer)20);
	_daskr_xerrwd_(msg, &c__20, &c__924, &c__0, &c__0, &c__0, &c__0, &c__1, &rl,
		&c_b38, (integer)80);
	_daskr_str_copy(msg, "------ NORM(F1) = (R1),  NORM(F1NEW) = (R2)", (integer)80,
		 (integer)43);
	_daskr_xerrwd_(msg, &c__43, &c__925, &c__0, &c__0, &c__0, &c__0, &c__2, &
		f1nrm, &f1nrmp, (integer)80);
    }
    if (f1nrmp > f1nrm + alpha * slpi * rl) {
	goto L200;
    }
/* ----------------------------------------------------------------------- */
/* Alpha-condition is satisfied, or linesearch is turned off. */
/* Copy YNEW,YPNEW to Y,YPRIME and return. */
/* ----------------------------------------------------------------------- */
L150:
    *iret = 0;
    _daskr_dcopy_(neq, &ynew[1], &c__1, &y[1], &c__1);
    _daskr_dcopy_(neq, &ypnew[1], &c__1, &yprime[1], &c__1);
    *fnrm = fnrmp;
    if (kprin >= 1) {
	_daskr_str_copy(msg, "------ LEAVING ROUTINE DLINSK, FNRM = (R1)", (integer)80,
		(integer)42);
	_daskr_xerrwd_(msg, &c__42, &c__926, &c__0, &c__0, &c__0, &c__0, &c__1, fnrm,
		 &c_b38, (integer)80);
    }
    return 0;
/* ----------------------------------------------------------------------- */
/* Alpha-condition not satisfied.  Perform backtrack to compute new RL */
/* value.  If RL is less than RLMIN, i.e. no satisfactory YNEW,YPNEW can */
/* be found sufficiently distinct from Y,YPRIME, then return IRET = 1. */
/* ----------------------------------------------------------------------- */
L200:
    if (rl < rlmin) {
	*iret = 1;
	return 0;
    }

    rl /= two;
    goto L100;

/* ----------------------- END OF SUBROUTINE DLINSK ---------------------- */
} /* dlinsk_ */

/* Subroutine */ int _daskr_dfnrmk_(integer *neq, real_number *y, real_number *t,
	real_number *yprime, real_number *savr, real_number *r__, real_number *cj,
	 real_number *tscale, real_number *wt, real_number *sqrtn, real_number *
	rsqrtn, Unknown_fp res, integer *ires, Unknown_fp psol, integer *irin, integer *
	ier, real_number *fnorm, real_number *eplin, real_number *wp, integer *
	iwp, real_number *pwk, real_number *rpar, integer *ipar)
{
    extern /* Subroutine */ int _daskr_dscal_(integer *, real_number *, real_number *,
	    integer *), _daskr_dcopy_(integer *, real_number *, integer *, real_number
	    *, integer *);
    extern real_number _daskr_ddwnrm_(integer *, real_number *, real_number *,
	    real_number *, integer *);


/* ***BEGIN PROLOGUE  DFNRMK */
/* ***REFER TO  DLINSK */
/* ***DATE WRITTEN   940830   (YYMMDD) */
/* ***REVISION DATE  951006   (SQRTN, RSQRTN, and scaling of WT added.) */
/* ***REVISION DATE  000628   TSCALE argument added. */


/* ----------------------------------------------------------------------- */
/* ***DESCRIPTION */

/*     DFNRMK calculates the scaled preconditioned norm of the nonlinear */
/*     function used in the nonlinear iteration for obtaining consistent */
/*     initial conditions.  Specifically, DFNRMK calculates the weighted */
/*     root-mean-square norm of the vector (P-inverse)*G(T,Y,YPRIME), */
/*     where P is the preconditioner matrix. */

/*     In addition to the parameters described in the calling program */
/*     DLINSK, the parameters represent */

/*     TSCALE -- Scale factor in T, used for stopping tests if nonzero. */
/*     IRIN   -- Flag showing whether the current residual vector is */
/*               input in SAVR.  1 means it is, 0 means it is not. */
/*     R      -- Array of length NEQ that contains */
/*               (P-inverse)*G(T,Y,YPRIME) on return. */
/*     FNORM  -- Scalar containing the weighted norm of R on return. */
/* ----------------------------------------------------------------------- */

/* ***ROUTINES CALLED */
/*   RES, DCOPY, DSCAL, PSOL, DDWNRM */

/* ***END PROLOGUE  DFNRMK */


/* ----------------------------------------------------------------------- */
/*     Call RES routine if IRIN = 0. */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --ipar;
    --rpar;
    --pwk;
    --iwp;
    --wp;
    --wt;
    --r__;
    --savr;
    --yprime;
    --y;

    /* Function Body */
    if (*irin == 0) {
	*ires = 0;
	((int (*)(real_number *, real_number *, real_number*, real_number*, real_number *, integer *, real_number *, integer *))*res)
    (t, &y[1], &yprime[1], cj, &savr[1], ires, &rpar[1], &ipar[1]);
	if (*ires < 0) {
	    return 0;
	}
    }
/* ----------------------------------------------------------------------- */
/*     Apply inverse of left preconditioner to vector R. */
/*     First scale WT array by 1/sqrt(N), and undo scaling afterward. */
/* ----------------------------------------------------------------------- */
    _daskr_dcopy_(neq, &savr[1], &c__1, &r__[1], &c__1);
    _daskr_dscal_(neq, rsqrtn, &wt[1], &c__1);
    *ier = 0;
    ((int (*)(integer*, real_number *, real_number *, real_number*, real_number *, real_number *, real_number *, real_number *, real_number *, integer *, real_number*, real_number*, integer *, real_number *, integer *))*psol)
      (neq, t, &y[1], &yprime[1], &savr[1], &pwk[1], cj, &wt[1], &wp[1], 
	    &iwp[1], &r__[1], eplin, ier, &rpar[1], &ipar[1]);
    _daskr_dscal_(neq, sqrtn, &wt[1], &c__1);
    if (*ier != 0) {
	return 0;
    }
/* ----------------------------------------------------------------------- */
/*     Calculate norm of R. */
/* ----------------------------------------------------------------------- */
    *fnorm = _daskr_ddwnrm_(neq, &r__[1], &wt[1], &rpar[1], &ipar[1]);
    if (*tscale > 0.) {
	*fnorm = *fnorm * *tscale * fabs(*cj);
    }

    return 0;
/* ----------------------- END OF SUBROUTINE DFNRMK ---------------------- */
} /* dfnrmk_ */

/* Subroutine */ int _daskr_dnedk_(real_number *x, real_number *y, real_number *yprime,
	integer *neq, Unknown_fp res, Unknown_fp jack, Unknown_fp psol, real_number *h__,
	real_number *wt, integer *jstart, integer *idid, real_number *rpar,
	integer *ipar, real_number *phi, real_number *gamma, real_number *savr,
	real_number *delta, real_number *e, real_number *wm, integer *iwm,
	real_number *cj, real_number *cjold, real_number *cjlast, real_number *s,
	real_number *uround, real_number *epli, real_number *sqrtn, real_number *
	rsqrtn, real_number *epcon, integer *jcalc, integer *jflg, integer *
	kp1, integer *nonneg, integer *ntype, integer *iernls)
{
    /* Initialized data */

    static integer muldel = 0;
    static integer maxit = 4;
    static real_number xrate = .25;

    /* System generated locals */
    integer phi_dim1, phi_offset, i__1, i__2;
    real_number d__1;

    /* Local variables */
    static integer i__, j, lwp;
    extern /* Subroutine */ int _daskr_dnsk_(real_number *, real_number *, real_number *
	    , integer *, Unknown_fp, Unknown_fp, real_number *, real_number *, integer *,
	    real_number *, real_number *, real_number *, real_number *, integer *,
	     real_number *, real_number *, real_number *, real_number *,
	    real_number *, real_number *, real_number *, real_number *, integer *,
	     integer *, integer *, integer *, integer *);
    static integer ires, liwp;
    static real_number temp1, temp2, eplin;
    static integer ierpj, iersl;
    static real_number delnrm;
    static integer iernew;
    extern real_number _daskr_ddwnrm_(integer *, real_number *, real_number *,
	    real_number *, integer *);
    static real_number tolnew;
    static integer iertyp;


/* ***BEGIN PROLOGUE  DNEDK */
/* ***REFER TO  DDASPK */
/* ***DATE WRITTEN   891219   (YYMMDD) */
/* ***REVISION DATE  900926   (YYMMDD) */
/* ***REVISION DATE  940701   (YYMMDD) */


/* ----------------------------------------------------------------------- */
/* ***DESCRIPTION */

/*     DNEDK solves a nonlinear system of */
/*     algebraic equations of the form */
/*     G(X,Y,YPRIME) = 0 for the unknown Y. */

/*     The method used is a matrix-free Newton scheme. */

/*     The parameters represent */
/*     X         -- Independent variable. */
/*     Y         -- Solution vector at x. */
/*     YPRIME    -- Derivative of solution vector */
/*                  after successful step. */
/*     NEQ       -- Number of equations to be integrated. */
/*     RES       -- External user-supplied subroutine */
/*                  to evaluate the residual.  See RES description */
/*                  in DDASPK prologue. */
/*     JACK     --  External user-supplied routine to update */
/*                  the preconditioner.  (This is optional). */
/*                  See JAC description for the case */
/*                  INFO(12) = 1 in the DDASPK prologue. */
/*     PSOL      -- External user-supplied routine to solve */
/*                  a linear system using preconditioning. */
/*                  (This is optional).  See explanation inside DDASPK. */
/*     H         -- Appropriate step size for this step. */
/*     WT        -- Vector of weights for error criterion. */
/*     JSTART    -- Indicates first call to this routine. */
/*                  If JSTART = 0, then this is the first call, */
/*                  otherwise it is not. */
/*     IDID      -- Completion flag, output by DNEDK. */
/*                  See IDID description in DDASPK prologue. */
/*     RPAR,IPAR -- Real and integer arrays used for communication */
/*                  between the calling program and external user */
/*                  routines.  They are not altered within DASPK. */
/*     PHI       -- Array of divided differences used by */
/*                  DNEDK.  The length is NEQ*(K+1), where */
/*                  K is the maximum order. */
/*     GAMMA     -- Array used to predict Y and YPRIME.  The length */
/*                  is K+1, where K is the maximum order. */
/*     SAVR      -- Work vector for DNEDK of length NEQ. */
/*     DELTA     -- Work vector for DNEDK of length NEQ. */
/*     E         -- Error accumulation vector for DNEDK of length NEQ. */
/*     WM,IWM    -- Real and integer arrays storing */
/*                  matrix information for linear system */
/*                  solvers, and various other information. */
/*     CJ        -- Parameter always proportional to 1/H. */
/*     CJOLD     -- Saves the value of CJ as of the last call to DITMD. */
/*                  Accounts for changes in CJ needed to */
/*                  decide whether to call DITMD. */
/*     CJLAST    -- Previous value of CJ. */
/*     S         -- A scalar determined by the approximate rate */
/*                  of convergence of the Newton iteration and used */
/*                  in the convergence test for the Newton iteration. */

/*                  If RATE is defined to be an estimate of the */
/*                  rate of convergence of the Newton iteration, */
/*                  then S = RATE/(1.D0-RATE). */

/*                  The closer RATE is to 0., the faster the Newton */
/*                  iteration is converging; the closer RATE is to 1., */
/*                  the slower the Newton iteration is converging. */

/*                  On the first Newton iteration with an up-dated */
/*                  preconditioner S = 100.D0, Thus the initial */
/*                  RATE of convergence is approximately 1. */

/*                  S is preserved from call to call so that the rate */
/*                  estimate from a previous step can be applied to */
/*                  the current step. */
/*     UROUND    -- Unit roundoff.  Not used here. */
/*     EPLI      -- convergence test constant. */
/*                  See DDASPK prologue for more details. */
/*     SQRTN     -- Square root of NEQ. */
/*     RSQRTN    -- reciprical of square root of NEQ. */
/*     EPCON     -- Tolerance to test for convergence of the Newton */
/*                  iteration. */
/*     JCALC     -- Flag used to determine when to update */
/*                  the Jacobian matrix.  In general: */

/*                  JCALC = -1 ==> Call the DITMD routine to update */
/*                                 the Jacobian matrix. */
/*                  JCALC =  0 ==> Jacobian matrix is up-to-date. */
/*                  JCALC =  1 ==> Jacobian matrix is out-dated, */
/*                                 but DITMD will not be called unless */
/*                                 JCALC is set to -1. */
/*     JFLG      -- Flag showing whether a Jacobian routine is supplied. */
/*     KP1       -- The current order + 1;  updated across calls. */
/*     NONNEG    -- Flag to determine nonnegativity constraints. */
/*     NTYPE     -- Identification code for the DNEDK routine. */
/*                   1 ==> modified Newton; iterative linear solver. */
/*                   2 ==> modified Newton; user-supplied linear solver. */
/*     IERNLS    -- Error flag for nonlinear solver. */
/*                   0 ==> nonlinear solver converged. */
/*                   1 ==> recoverable error inside non-linear solver. */
/*                  -1 ==> unrecoverable error inside non-linear solver. */

/*     The following group of variables are passed as arguments to */
/*     the Newton iteration solver.  They are explained in greater detail */
/*     in DNSK: */
/*        TOLNEW, MULDEL, MAXIT, IERNEW */

/*     IERTYP -- Flag which tells whether this subroutine is correct. */
/*               0 ==> correct subroutine. */
/*               1 ==> incorrect subroutine. */

/* ----------------------------------------------------------------------- */
/* ***ROUTINES CALLED */
/*   RES, JACK, DDWNRM, DNSK */

/* ***END PROLOGUE  DNEDK */




    /* Parameter adjustments */
    --y;
    --yprime;
    phi_dim1 = *neq;
    phi_offset = 1 + phi_dim1;
    phi -= phi_offset;
    --wt;
    --rpar;
    --ipar;
    --gamma;
    --savr;
    --delta;
    --e;
    --wm;
    --iwm;

    /* Function Body */

/*     Verify that this is the correct subroutine. */

    iertyp = 0;
    if (*ntype != 1) {
	iertyp = 1;
	goto L380;
    }

/*     If this is the first step, perform initializations. */

    if (*jstart == 0) {
	*cjold = *cj;
	*jcalc = -1;
	*s = 100.;
    }

/*     Perform all other initializations. */

    *iernls = 0;
    lwp = iwm[29];
    liwp = iwm[30];

/*     Decide whether to update the preconditioner. */

    if (*jflg != 0) {
	temp1 = (1. - xrate) / (xrate + 1.);
	temp2 = 1. / temp1;
	if (*cj / *cjold < temp1 || *cj / *cjold > temp2) {
	    *jcalc = -1;
	}
	if (*cj != *cjlast) {
	    *s = 100.;
	}
    } else {
	*jcalc = 0;
    }

/*     Looping point for updating preconditioner with current stepsize. */

L300:

/*     Initialize all error flags to zero. */

    ierpj = 0;
    ires = 0;
    iersl = 0;
    iernew = 0;

/*     Predict the solution and derivative and compute the tolerance */
/*     for the Newton iteration. */

    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] = phi[i__ + phi_dim1];
/* L310: */
	yprime[i__] = 0.;
    }
    i__1 = *kp1;
    for (j = 2; j <= i__1; ++j) {
	i__2 = *neq;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y[i__] += phi[i__ + j * phi_dim1];
/* L320: */
	    yprime[i__] += gamma[j] * phi[i__ + j * phi_dim1];
	}
/* L330: */
    }
    eplin = *epli * *epcon;
    tolnew = eplin;

/*     Call RES to initialize DELTA. */

    ++iwm[12];
    ((int (*)(real_number *, real_number *, real_number*, real_number*, real_number *, integer *, real_number *, integer *))*res)
      (x, &y[1], &yprime[1], cj, &delta[1], &ires, &rpar[1], &ipar[1]);
    if (ires < 0) {
	goto L380;
    }


/*     If indicated, update the preconditioner. */
/*     Set JCALC to 0 as an indicator that this has been done. */

    if (*jcalc == -1) {
	++iwm[13];
	*jcalc = 0;
	((int (*)(Unknown_fp, integer *, integer *, real_number*, real_number *, real_number *, real_number *, real_number *,real_number *, real_number *, real_number*, real_number*, integer *, integer *, real_number *, integer *))*jack)
  ((Unknown_fp)res, &ires, neq, x, &y[1], &yprime[1], &wt[1], &delta[1]
		, &e[1], h__, cj, &wm[lwp], &iwm[liwp], &ierpj, &rpar[1], &
		ipar[1]);
	*cjold = *cj;
	*s = 100.;
	if (ires < 0) {
	    goto L380;
	}
	if (ierpj != 0) {
	    goto L380;
	}
    }

/*     Call the nonlinear Newton solver. */

    _daskr_dnsk_(x, &y[1], &yprime[1], neq, (Unknown_fp)res, (Unknown_fp)psol, &wt[1], &rpar[1],
	    &ipar[1], &savr[1], &delta[1], &e[1], &wm[1], &iwm[1], cj, sqrtn, 
	    rsqrtn, &eplin, epcon, s, &temp1, &tolnew, &muldel, &maxit, &ires,
	     &iersl, &iernew);

    if (iernew > 0 && *jcalc != 0) {

/*     The Newton iteration had a recoverable failure with an old */
/*     preconditioner.  Retry the step with a new preconditioner. */

	*jcalc = -1;
	goto L300;
    }

    if (iernew != 0) {
	goto L380;
    }

/*     The Newton iteration has converged.  If nonnegativity of */
/*     solution is required, set the solution nonnegative, if the */
/*     perturbation to do it is small enough.  If the change is too */
/*     large, then consider the corrector iteration to have failed. */

    if (*nonneg == 0) {
	goto L390;
    }
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L360: */
/* Computing MIN */
	d__1 = y[i__];
	delta[i__] = MIN(d__1,0.);
    }
    delnrm = _daskr_ddwnrm_(neq, &delta[1], &wt[1], &rpar[1], &ipar[1]);
    if (delnrm > *epcon) {
	goto L380;
    }
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L370: */
	e[i__] -= delta[i__];
    }
    goto L390;


/*     Exits from nonlinear solver. */
/*     No convergence with current preconditioner. */
/*     Compute IERNLS and IDID accordingly. */

L380:
    if (ires <= -2 || iersl < 0 || iertyp != 0) {
	*iernls = -1;
	if (ires <= -2) {
	    *idid = -11;
	}
	if (iersl < 0) {
	    *idid = -13;
	}
	if (iertyp != 0) {
	    *idid = -15;
	}
    } else {
	*iernls = 1;
	if (ires == -1) {
	    *idid = -10;
	}
	if (ierpj != 0) {
	    *idid = -5;
	}
	if (iersl > 0) {
	    *idid = -14;
	}
    }


L390:
    *jcalc = 1;
    return 0;

/* ------END OF SUBROUTINE DNEDK------------------------------------------ */
} /* dnedk_ */

/* Subroutine */ int _daskr_dnsk_(real_number *x, real_number *y, real_number *yprime,
	integer *neq, Unknown_fp res, Unknown_fp psol, real_number *wt, real_number *rpar,
	integer *ipar, real_number *savr, real_number *delta, real_number *e,
	real_number *wm, integer *iwm, real_number *cj, real_number *sqrtn,
	real_number *rsqrtn, real_number *eplin, real_number *epcon, real_number *
	s, real_number *confac, real_number *tolnew, integer *muldel, integer *
	maxit, integer *ires, integer *iersl, integer *iernew)
{
    /* System generated locals */
    integer i__1;
    real_number d__1, d__2;

    /* Builtin functions */
    double _daskr_real_pow(real_number *, real_number *);

    /* Local variables */
    static integer i__, m;
    static real_number rate, rhok;
    extern /* Subroutine */ int _daskr_dslvk_(integer *, real_number *, real_number *,
	    real_number *, real_number *, real_number *, real_number *,
	    real_number *, integer *, Unknown_fp, integer *, Unknown_fp, integer *,
	    real_number *, real_number *, real_number *, real_number *,
	    real_number *, real_number *, integer *);
    static real_number delnrm;
    extern real_number _daskr_ddwnrm_(integer *, real_number *, real_number *,
	    real_number *, integer *);
    static real_number oldnrm;


/* ***BEGIN PROLOGUE  DNSK */
/* ***REFER TO  DDASPK */
/* ***DATE WRITTEN   891219   (YYMMDD) */
/* ***REVISION DATE  900926   (YYMMDD) */
/* ***REVISION DATE  950126   (YYMMDD) */
/* ***REVISION DATE  000711   (YYMMDD) */


/* ----------------------------------------------------------------------- */
/* ***DESCRIPTION */

/*     DNSK solves a nonlinear system of */
/*     algebraic equations of the form */
/*     G(X,Y,YPRIME) = 0 for the unknown Y. */

/*     The method used is a modified Newton scheme. */

/*     The parameters represent */

/*     X         -- Independent variable. */
/*     Y         -- Solution vector. */
/*     YPRIME    -- Derivative of solution vector. */
/*     NEQ       -- Number of unknowns. */
/*     RES       -- External user-supplied subroutine */
/*                  to evaluate the residual.  See RES description */
/*                  in DDASPK prologue. */
/*     PSOL      -- External user-supplied routine to solve */
/*                  a linear system using preconditioning. */
/*                  See explanation inside DDASPK. */
/*     WT        -- Vector of weights for error criterion. */
/*     RPAR,IPAR -- Real and integer arrays used for communication */
/*                  between the calling program and external user */
/*                  routines.  They are not altered within DASPK. */
/*     SAVR      -- Work vector for DNSK of length NEQ. */
/*     DELTA     -- Work vector for DNSK of length NEQ. */
/*     E         -- Error accumulation vector for DNSK of length NEQ. */
/*     WM,IWM    -- Real and integer arrays storing */
/*                  matrix information such as the matrix */
/*                  of partial derivatives, permutation */
/*                  vector, and various other information. */
/*     CJ        -- Parameter always proportional to 1/H (step size). */
/*     SQRTN     -- Square root of NEQ. */
/*     RSQRTN    -- reciprical of square root of NEQ. */
/*     EPLIN     -- Tolerance for linear system solver. */
/*     EPCON     -- Tolerance to test for convergence of the Newton */
/*                  iteration. */
/*     S         -- Used for error convergence tests. */
/*                  In the Newton iteration: S = RATE/(1.D0-RATE), */
/*                  where RATE is the estimated rate of convergence */
/*                  of the Newton iteration. */

/*                  The closer RATE is to 0., the faster the Newton */
/*                  iteration is converging; the closer RATE is to 1., */
/*                  the slower the Newton iteration is converging. */

/*                  The calling routine sends the initial value */
/*                  of S to the Newton iteration. */
/*     CONFAC    -- A residual scale factor to improve convergence. */
/*     TOLNEW    -- Tolerance on the norm of Newton correction in */
/*                  alternative Newton convergence test. */
/*     MULDEL    -- A flag indicating whether or not to multiply */
/*                  DELTA by CONFAC. */
/*                  0  ==> do not scale DELTA by CONFAC. */
/*                  1  ==> scale DELTA by CONFAC. */
/*     MAXIT     -- Maximum allowed number of Newton iterations. */
/*     IRES      -- Error flag returned from RES.  See RES description */
/*                  in DDASPK prologue.  If IRES = -1, then IERNEW */
/*                  will be set to 1. */
/*                  If IRES < -1, then IERNEW will be set to -1. */
/*     IERSL     -- Error flag for linear system solver. */
/*                  See IERSL description in subroutine DSLVK. */
/*                  If IERSL = 1, then IERNEW will be set to 1. */
/*                  If IERSL < 0, then IERNEW will be set to -1. */
/*     IERNEW    -- Error flag for Newton iteration. */
/*                   0  ==> Newton iteration converged. */
/*                   1  ==> recoverable error inside Newton iteration. */
/*                  -1  ==> unrecoverable error inside Newton iteration. */
/* ----------------------------------------------------------------------- */

/* ***ROUTINES CALLED */
/*   RES, DSLVK, DDWNRM */

/* ***END PROLOGUE  DNSK */




/*     Initialize Newton counter M and accumulation vector E. */

    /* Parameter adjustments */
    --iwm;
    --wm;
    --e;
    --delta;
    --savr;
    --ipar;
    --rpar;
    --wt;
    --yprime;
    --y;

    /* Function Body */
    m = 0;
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L100: */
	e[i__] = 0.;
    }

/*     Corrector loop. */

L300:
    ++iwm[19];

/*     If necessary, multiply residual by convergence factor. */

    if (*muldel == 1) {
	i__1 = *neq;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L320: */
	    delta[i__] *= *confac;
	}
    }

/*     Save residual in SAVR. */

    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L340: */
	savr[i__] = delta[i__];
    }

/*     Compute a new iterate.  Store the correction in DELTA. */

    _daskr_dslvk_(neq, &y[1], x, &yprime[1], &savr[1], &delta[1], &wt[1], &wm[1], &
	    iwm[1], (Unknown_fp)res, ires, (Unknown_fp)psol, iersl, cj, eplin, sqrtn,
	    rsqrtn, &rhok, &rpar[1], &ipar[1]);
    if (*ires != 0 || *iersl != 0) {
	goto L380;
    }

/*     Update Y, E, and YPRIME. */

    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] -= delta[i__];
	e[i__] -= delta[i__];
/* L360: */
	yprime[i__] -= *cj * delta[i__];
    }

/*     Test for convergence of the iteration. */

    delnrm = _daskr_ddwnrm_(neq, &delta[1], &wt[1], &rpar[1], &ipar[1]);
    if (m == 0) {
	oldnrm = delnrm;
	if (delnrm <= *tolnew) {
	    goto L370;
	}
    } else {
	d__1 = delnrm / oldnrm;
	d__2 = 1. / m;
	rate = _daskr_real_pow(&d__1, &d__2);
	if (rate > .9) {
	    goto L380;
	}
	*s = rate / (1. - rate);
    }
    if (*s * delnrm <= *epcon) {
	goto L370;
    }

/*     The corrector has not yet converged.  Update M and test whether */
/*     the maximum number of iterations have been tried. */

    ++m;
    if (m >= *maxit) {
	goto L380;
    }

/*     Evaluate the residual, and go back to do another iteration. */

    ++iwm[12];
    ((int (*)(real_number *, real_number *, real_number*, real_number*, real_number *, integer *, real_number *, integer *))*res)
      (x, &y[1], &yprime[1], cj, &delta[1], ires, &rpar[1], &ipar[1]);
    if (*ires < 0) {
	goto L380;
    }
    goto L300;

/*     The iteration has converged. */

L370:
    return 0;

/*     The iteration has not converged.  Set IERNEW appropriately. */

L380:
    if (*ires <= -2 || *iersl < 0) {
	*iernew = -1;
    } else {
	*iernew = 1;
    }
    return 0;


/* ------END OF SUBROUTINE DNSK------------------------------------------- */
} /* dnsk_ */

/* Subroutine */ int _daskr_dslvk_(integer *neq, real_number *y, real_number *tn,
	real_number *yprime, real_number *savr, real_number *x, real_number *ewt,
	real_number *wm, integer *iwm, Unknown_fp res, integer *ires, Unknown_fp psol,
	integer *iersl, real_number *cj, real_number *eplin, real_number *sqrtn,
	real_number *rsqrtn, real_number *rhok, real_number *rpar, integer *ipar)
{
    /* Initialized data */

    static integer irst = 1;

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, lq, lr, lv, lz, ldl, nli, nre, kmp, lwk, nps, lwp, 
	    ncfl, lhes, lgmr, maxl, nres, npsl, liwp, iflag;
    extern /* Subroutine */ int _daskr_dscal_(integer *, real_number *, real_number *,
	    integer *), _daskr_dcopy_(integer *, real_number *, integer *, real_number
	    *, integer *);
    static integer miter, nrmax, nrsts, maxlp1;
    extern /* Subroutine */ int _daskr_dspigm_(integer *, real_number *, real_number *,
	     real_number *, real_number *, real_number *, real_number *, integer *
	    , integer *, integer *, real_number *, real_number *, Unknown_fp, integer
	    *, integer *, Unknown_fp, integer *, real_number *, real_number *,
	    real_number *, real_number *, integer *, real_number *, integer *,
	    real_number *, real_number *, real_number *, integer *, integer *,
	    integer *, real_number *, integer *);


/* ***BEGIN PROLOGUE  DSLVK */
/* ***REFER TO  DDASPK */
/* ***DATE WRITTEN   890101   (YYMMDD) */
/* ***REVISION DATE  900926   (YYMMDD) */
/* ***REVISION DATE  940928   Removed MNEWT and added RHOK in call list. */


/* ----------------------------------------------------------------------- */
/* ***DESCRIPTION */

/* DSLVK uses a restart algorithm and interfaces to DSPIGM for */
/* the solution of the linear system arising from a Newton iteration. */

/* In addition to variables described elsewhere, */
/* communication with DSLVK uses the following variables.. */
/* WM    = Real work space containing data for the algorithm */
/*         (Krylov basis vectors, Hessenberg matrix, etc.). */
/* IWM   = Integer work space containing data for the algorithm. */
/* X     = The right-hand side vector on input, and the solution vector */
/*         on output, of length NEQ. */
/* IRES  = Error flag from RES. */
/* IERSL = Output flag .. */
/*         IERSL =  0 means no trouble occurred (or user RES routine */
/*                    returned IRES < 0) */
/*         IERSL =  1 means the iterative method failed to converge */
/*                    (DSPIGM returned IFLAG > 0.) */
/*         IERSL = -1 means there was a nonrecoverable error in the */
/*                    iterative solver, and an error exit will occur. */
/* ----------------------------------------------------------------------- */
/* ***ROUTINES CALLED */
/*   DSCAL, DCOPY, DSPIGM */

/* ***END PROLOGUE  DSLVK */




/* ----------------------------------------------------------------------- */
/* IRST is set to 1, to indicate restarting is in effect. */
/* NRMAX is the maximum number of restarts. */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --ipar;
    --rpar;
    --iwm;
    --wm;
    --ewt;
    --x;
    --savr;
    --yprime;
    --y;

    /* Function Body */

    liwp = iwm[30];
    nli = iwm[20];
    nps = iwm[21];
    ncfl = iwm[16];
    nre = iwm[12];
    lwp = iwm[29];
    maxl = iwm[24];
    kmp = iwm[25];
    nrmax = iwm[26];
    miter = iwm[23];
    *iersl = 0;
    *ires = 0;
/* ----------------------------------------------------------------------- */
/* Use a restarting strategy to solve the linear system */
/* P*X = -F.  Parse the work vector, and perform initializations. */
/* Note that zero is the initial guess for X. */
/* ----------------------------------------------------------------------- */
    maxlp1 = maxl + 1;
    lv = 1;
    lr = lv + *neq * maxl;
    lhes = lr + *neq + 1;
    lq = lhes + maxl * maxlp1;
    lwk = lq + (maxl << 1);
/* Computing MIN */
    i__1 = 1, i__2 = maxl - kmp;
    ldl = lwk + MIN(i__1,i__2) * *neq;
    lz = ldl + *neq;
    _daskr_dscal_(neq, rsqrtn, &ewt[1], &c__1);
    _daskr_dcopy_(neq, &x[1], &c__1, &wm[lr], &c__1);
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L110: */
	x[i__] = 0.;
    }
/* ----------------------------------------------------------------------- */
/* Top of loop for the restart algorithm.  Initial pass approximates */
/* X and sets up a transformed system to perform subsequent restarts */
/* to update X.  NRSTS is initialized to -1, because restarting */
/* does not occur until after the first pass. */
/* Update NRSTS; conditionally copy DL to R; call the DSPIGM */
/* algorithm to solve A*Z = R;  updated counters;  update X with */
/* the residual solution. */
/* Note:  if convergence is not achieved after NRMAX restarts, */
/* then the linear solver is considered to have failed. */
/* ----------------------------------------------------------------------- */
    nrsts = -1;
L115:
    ++nrsts;
    if (nrsts > 0) {
	_daskr_dcopy_(neq, &wm[ldl], &c__1, &wm[lr], &c__1);
    }
    _daskr_dspigm_(neq, tn, &y[1], &yprime[1], &savr[1], &wm[lr], &ewt[1], &maxl, &
	    maxlp1, &kmp, eplin, cj, (Unknown_fp)res, ires, &nres, (Unknown_fp)psol, &
	    npsl, &wm[lz], &wm[lv], &wm[lhes], &wm[lq], &lgmr, &wm[lwp], &iwm[
	    liwp], &wm[lwk], &wm[ldl], rhok, &iflag, &irst, &nrsts, &rpar[1], 
	    &ipar[1]);
    nli += lgmr;
    nps += npsl;
    nre += nres;
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L120: */
	x[i__] += wm[lz + i__ - 1];
    }
    if (iflag == 1 && nrsts < nrmax && *ires == 0) {
	goto L115;
    }
/* ----------------------------------------------------------------------- */
/* The restart scheme is finished.  Test IRES and IFLAG to see if */
/* convergence was not achieved, and set flags accordingly. */
/* ----------------------------------------------------------------------- */
    if (*ires < 0) {
	++ncfl;
    } else if (iflag != 0) {
	++ncfl;
	if (iflag > 0) {
	    *iersl = 1;
	}
	if (iflag < 0) {
	    *iersl = -1;
	}
    }
/* ----------------------------------------------------------------------- */
/* Update IWM with counters, rescale EWT, and return. */
/* ----------------------------------------------------------------------- */
    iwm[20] = nli;
    iwm[21] = nps;
    iwm[16] = ncfl;
    iwm[12] = nre;
    _daskr_dscal_(neq, sqrtn, &ewt[1], &c__1);
    return 0;

/* ------END OF SUBROUTINE DSLVK------------------------------------------ */
} /* dslvk_ */

/* Subroutine */ int _daskr_dspigm_(integer *neq, real_number *tn, real_number *y,
	real_number *yprime, real_number *savr, real_number *r__, real_number *
	wght, integer *maxl, integer *maxlp1, integer *kmp, real_number *eplin,
	 real_number *cj, Unknown_fp res, integer *ires, integer *nre, Unknown_fp psol,
	integer *npsl, real_number *z__, real_number *v, real_number *hes,
	real_number *q, integer *lgmr, real_number *wp, integer *iwp,
	real_number *wk, real_number *dl, real_number *rhok, integer *iflag,
	integer *irst, integer *nrsts, real_number *rpar, integer *ipar)
{
    /* System generated locals */
    integer v_dim1, v_offset, hes_dim1, hes_offset, i__1, i__2, i__3;
    real_number d__1;

    /* Local variables */
    static real_number c__;
    static integer i__, j, k;
    static real_number s;
    static integer i2, ll, ip1, ier;
    static real_number tem, rho;
    static integer llp1, info;
    extern /* Subroutine */ int _daskr_datv_(integer *, real_number *, real_number *,
	    real_number *, real_number *, real_number *, real_number *,
	    real_number *, Unknown_fp, integer *, Unknown_fp, real_number *, real_number *,
	    real_number *, integer *, real_number *, real_number *, integer *,
	    integer *, integer *, real_number *, integer *);
    static real_number prod, rnrm;
    extern real_number _daskr_dnrm2_(integer *, real_number *, integer *);
    extern /* Subroutine */ int _daskr_dscal_(integer *, real_number *, real_number *,
	    integer *), _daskr_dhels_(real_number *, integer *, integer *, real_number
	    *, real_number *), _daskr_dheqr_(real_number *, integer *, integer *,
	    real_number *, integer *, integer *);
    static real_number dlnrm;
    extern /* Subroutine */ int _daskr_dcopy_(integer *, real_number *, integer *,
	    real_number *, integer *), _daskr_dorth_(real_number *, real_number *,
	    real_number *, integer *, integer *, integer *, integer *,
	    real_number *), _daskr_daxpy_(integer *, real_number *, real_number *,
	    integer *, real_number *, integer *);
    static integer maxlm1;
    static real_number snormw;


/* ***BEGIN PROLOGUE  DSPIGM */
/* ***DATE WRITTEN   890101   (YYMMDD) */
/* ***REVISION DATE  900926   (YYMMDD) */
/* ***REVISION DATE  940927   Removed MNEWT and added RHOK in call list. */


/* ----------------------------------------------------------------------- */
/* ***DESCRIPTION */

/* This routine solves the linear system A * Z = R using a scaled */
/* preconditioned version of the generalized minimum residual method. */
/* An initial guess of Z = 0 is assumed. */

/*      On entry */

/*          NEQ = Problem size, passed to PSOL. */

/*           TN = Current Value of T. */

/*            Y = Array Containing current dependent variable vector. */

/*       YPRIME = Array Containing current first derivative of Y. */

/*         SAVR = Array containing current value of G(T,Y,YPRIME). */

/*            R = The right hand side of the system A*Z = R. */
/*                R is also used as work space when computing */
/*                the final approximation and will therefore be */
/*                destroyed. */
/*                (R is the same as V(*,MAXL+1) in the call to DSPIGM.) */

/*         WGHT = The vector of length NEQ containing the nonzero */
/*                elements of the diagonal scaling matrix. */

/*         MAXL = The maximum allowable order of the matrix H. */

/*       MAXLP1 = MAXL + 1, used for dynamic dimensioning of HES. */

/*          KMP = The number of previous vectors the new vector, VNEW, */
/*                must be made orthogonal to.  (KMP .LE. MAXL.) */

/*        EPLIN = Tolerance on residuals R-A*Z in weighted rms norm. */

/*           CJ = Scalar proportional to current value of */
/*                1/(step size H). */

/*           WK = Real work array used by routine DATV and PSOL. */

/*           DL = Real work array used for calculation of the residual */
/*                norm RHO when the method is incomplete (KMP.LT.MAXL) */
/*                and/or when using restarting. */

/*           WP = Real work array used by preconditioner PSOL. */

/*          IWP = Integer work array used by preconditioner PSOL. */

/*         IRST = Method flag indicating if restarting is being */
/*                performed.  IRST .GT. 0 means restarting is active, */
/*                while IRST = 0 means restarting is not being used. */

/*        NRSTS = Counter for the number of restarts on the current */
/*                call to DSPIGM.  If NRSTS .GT. 0, then the residual */
/*                R is already scaled, and so scaling of R is not */
/*                necessary. */


/*      On Return */

/*         Z    = The final computed approximation to the solution */
/*                of the system A*Z = R. */

/*         LGMR = The number of iterations performed and */
/*                the current order of the upper Hessenberg */
/*                matrix HES. */

/*         NRE  = The number of calls to RES (i.e. DATV) */

/*         NPSL = The number of calls to PSOL. */

/*         V    = The neq by (LGMR+1) array containing the LGMR */
/*                orthogonal vectors V(*,1) to V(*,LGMR). */

/*         HES  = The upper triangular factor of the QR decomposition */
/*                of the (LGMR+1) by LGMR upper Hessenberg matrix whose */
/*                entries are the scaled inner-products of A*V(*,I) */
/*                and V(*,K). */

/*         Q    = Real array of length 2*MAXL containing the components */
/*                of the givens rotations used in the QR decomposition */
/*                of HES.  It is loaded in DHEQR and used in DHELS. */

/*         IRES = Error flag from RES. */

/*           DL = Scaled preconditioned residual, */
/*                (D-inverse)*(P-inverse)*(R-A*Z). Only loaded when */
/*                performing restarts of the Krylov iteration. */

/*         RHOK = Weighted norm of final preconditioned residual. */

/*        IFLAG = Integer error flag.. */
/*                0 Means convergence in LGMR iterations, LGMR.LE.MAXL. */
/*                1 Means the convergence test did not pass in MAXL */
/*                  iterations, but the new residual norm (RHO) is */
/*                  .LT. the old residual norm (RNRM), and so Z is */
/*                  computed. */
/*                2 Means the convergence test did not pass in MAXL */
/*                  iterations, new residual norm (RHO) .GE. old residual */
/*                  norm (RNRM), and the initial guess, Z = 0, is */
/*                  returned. */
/*                3 Means there was a recoverable error in PSOL */
/*                  caused by the preconditioner being out of date. */
/*               -1 Means there was an unrecoverable error in PSOL. */

/* ----------------------------------------------------------------------- */
/* ***ROUTINES CALLED */
/*   PSOL, DNRM2, DSCAL, DATV, DORTH, DHEQR, DCOPY, DHELS, DAXPY */

/* ***END PROLOGUE  DSPIGM */


    /* Parameter adjustments */
    v_dim1 = *neq;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --y;
    --yprime;
    --savr;
    --r__;
    --wght;
    hes_dim1 = *maxlp1;
    hes_offset = 1 + hes_dim1;
    hes -= hes_offset;
    --z__;
    --q;
    --wp;
    --iwp;
    --wk;
    --dl;
    --rpar;
    --ipar;

    /* Function Body */
    ier = 0;
    *iflag = 0;
    *lgmr = 0;
    *npsl = 0;
    *nre = 0;
/* ----------------------------------------------------------------------- */
/* The initial guess for Z is 0.  The initial residual is therefore */
/* the vector R.  Initialize Z to 0. */
/* ----------------------------------------------------------------------- */
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	z__[i__] = 0.;
    }
/* ----------------------------------------------------------------------- */
/* Apply inverse of left preconditioner to vector R if NRSTS .EQ. 0. */
/* Form V(*,1), the scaled preconditioned right hand side. */
/* ----------------------------------------------------------------------- */
    if (*nrsts == 0) {
   ((int (*)(integer*, real_number *, real_number *, real_number*, real_number *, real_number *, real_number *, real_number *, real_number *, integer *, real_number*, real_number*, integer *, real_number *, integer *))*psol)
    (neq, tn, &y[1], &yprime[1], &savr[1], &wk[1], cj, &wght[1], &
		wp[1], &iwp[1], &r__[1], eplin, &ier, &rpar[1], &ipar[1]);
	*npsl = 1;
	if (ier != 0) {
	    goto L300;
	}
	i__1 = *neq;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L30: */
	    v[i__ + v_dim1] = r__[i__] * wght[i__];
	}
    } else {
	i__1 = *neq;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L35: */
	    v[i__ + v_dim1] = r__[i__];
	}
    }
/* ----------------------------------------------------------------------- */
/* Calculate norm of scaled vector V(*,1) and normalize it */
/* If, however, the norm of V(*,1) (i.e. the norm of the preconditioned */
/* residual) is .le. EPLIN, then return with Z=0. */
/* ----------------------------------------------------------------------- */
    rnrm = _daskr_dnrm2_(neq, &v[v_offset], &c__1);
    if (rnrm <= *eplin) {
	*rhok = rnrm;
	return 0;
    }
    tem = 1. / rnrm;
    _daskr_dscal_(neq, &tem, &v[v_dim1 + 1], &c__1);
/* ----------------------------------------------------------------------- */
/* Zero out the HES array. */
/* ----------------------------------------------------------------------- */
    i__1 = *maxl;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *maxlp1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L60: */
	    hes[i__ + j * hes_dim1] = 0.;
	}
/* L65: */
    }
/* ----------------------------------------------------------------------- */
/* Main loop to compute the vectors V(*,2) to V(*,MAXL). */
/* The running product PROD is needed for the convergence test. */
/* ----------------------------------------------------------------------- */
    prod = 1.;
    i__1 = *maxl;
    for (ll = 1; ll <= i__1; ++ll) {
	*lgmr = ll;
/* ----------------------------------------------------------------------- */
/* Call routine DATV to compute VNEW = ABAR*V(LL), where ABAR is */
/* the matrix A with scaling and inverse preconditioner factors applied. */
/* Call routine DORTH to orthogonalize the new vector VNEW = V(*,LL+1). */
/* call routine DHEQR to update the factors of HES. */
/* ----------------------------------------------------------------------- */
	_daskr_datv_(neq, &y[1], tn, &yprime[1], &savr[1], &v[ll * v_dim1 + 1], &
		wght[1], &z__[1], (Unknown_fp)res, ires, (Unknown_fp)psol, &v[(ll + 1) *
		v_dim1 + 1], &wk[1], &wp[1], &iwp[1], cj, eplin, &ier, nre, 
		npsl, &rpar[1], &ipar[1]);
	if (*ires < 0) {
	    return 0;
	}
	if (ier != 0) {
	    goto L300;
	}
	_daskr_dorth_(&v[(ll + 1) * v_dim1 + 1], &v[v_offset], &hes[hes_offset], neq,
		 &ll, maxlp1, kmp, &snormw);
	hes[ll + 1 + ll * hes_dim1] = snormw;
	_daskr_dheqr_(&hes[hes_offset], maxlp1, &ll, &q[1], &info, &ll);
	if (info == ll) {
	    goto L120;
	}
/* ----------------------------------------------------------------------- */
/* Update RHO, the estimate of the norm of the residual R - A*ZL. */
/* If KMP .LT. MAXL, then the vectors V(*,1),...,V(*,LL+1) are not */
/* necessarily orthogonal for LL .GT. KMP.  The vector DL must then */
/* be computed, and its norm used in the calculation of RHO. */
/* ----------------------------------------------------------------------- */
	prod *= q[ll * 2];
	rho = (d__1 = prod * rnrm, fabs(d__1));
	if (ll > *kmp && *kmp < *maxl) {
	    if (ll == *kmp + 1) {
		_daskr_dcopy_(neq, &v[v_dim1 + 1], &c__1, &dl[1], &c__1);
		i__2 = *kmp;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    ip1 = i__ + 1;
		    i2 = i__ << 1;
		    s = q[i2];
		    c__ = q[i2 - 1];
		    i__3 = *neq;
		    for (k = 1; k <= i__3; ++k) {
/* L70: */
			dl[k] = s * dl[k] + c__ * v[k + ip1 * v_dim1];
		    }
/* L75: */
		}
	    }
	    s = q[ll * 2];
	    c__ = q[(ll << 1) - 1] / snormw;
	    llp1 = ll + 1;
	    i__2 = *neq;
	    for (k = 1; k <= i__2; ++k) {
/* L80: */
		dl[k] = s * dl[k] + c__ * v[k + llp1 * v_dim1];
	    }
	    dlnrm = _daskr_dnrm2_(neq, &dl[1], &c__1);
	    rho *= dlnrm;
	}
/* ----------------------------------------------------------------------- */
/* Test for convergence.  If passed, compute approximation ZL. */
/* If failed and LL .LT. MAXL, then continue iterating. */
/* ----------------------------------------------------------------------- */
	if (rho <= *eplin) {
	    goto L200;
	}
	if (ll == *maxl) {
	    goto L100;
	}
/* ----------------------------------------------------------------------- */
/* Rescale so that the norm of V(1,LL+1) is one. */
/* ----------------------------------------------------------------------- */
	tem = 1. / snormw;
	_daskr_dscal_(neq, &tem, &v[(ll + 1) * v_dim1 + 1], &c__1);
/* L90: */
    }
L100:
    if (rho < rnrm) {
	goto L150;
    }
L120:
    *iflag = 2;
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L130: */
	z__[i__] = 0.;
    }
    return 0;
L150:
    *iflag = 1;
/* ----------------------------------------------------------------------- */
/* The tolerance was not met, but the residual norm was reduced. */
/* If performing restarting (IRST .gt. 0) calculate the residual vector */
/* RL and store it in the DL array.  If the incomplete version is */
/* being used (KMP .lt. MAXL) then DL has already been calculated. */
/* ----------------------------------------------------------------------- */
    if (*irst > 0) {
	if (*kmp == *maxl) {

/*           Calculate DL from the V(I)'s. */

	    _daskr_dcopy_(neq, &v[v_dim1 + 1], &c__1, &dl[1], &c__1);
	    maxlm1 = *maxl - 1;
	    i__1 = maxlm1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		ip1 = i__ + 1;
		i2 = i__ << 1;
		s = q[i2];
		c__ = q[i2 - 1];
		i__2 = *neq;
		for (k = 1; k <= i__2; ++k) {
/* L170: */
		    dl[k] = s * dl[k] + c__ * v[k + ip1 * v_dim1];
		}
/* L175: */
	    }
	    s = q[*maxl * 2];
	    c__ = q[(*maxl << 1) - 1] / snormw;
	    i__1 = *neq;
	    for (k = 1; k <= i__1; ++k) {
/* L180: */
		dl[k] = s * dl[k] + c__ * v[k + *maxlp1 * v_dim1];
	    }
	}

/*        Scale DL by RNRM*PROD to obtain the residual RL. */

	tem = rnrm * prod;
	_daskr_dscal_(neq, &tem, &dl[1], &c__1);
    }
/* ----------------------------------------------------------------------- */
/* Compute the approximation ZL to the solution. */
/* Since the vector Z was used as work space, and the initial guess */
/* of the Newton correction is zero, Z must be reset to zero. */
/* ----------------------------------------------------------------------- */
L200:
    ll = *lgmr;
    llp1 = ll + 1;
    i__1 = llp1;
    for (k = 1; k <= i__1; ++k) {
/* L210: */
	r__[k] = 0.;
    }
    r__[1] = rnrm;
    _daskr_dhels_(&hes[hes_offset], maxlp1, &ll, &q[1], &r__[1]);
    i__1 = *neq;
    for (k = 1; k <= i__1; ++k) {
/* L220: */
	z__[k] = 0.;
    }
    i__1 = ll;
    for (i__ = 1; i__ <= i__1; ++i__) {
	_daskr_daxpy_(neq, &r__[i__], &v[i__ * v_dim1 + 1], &c__1, &z__[1], &c__1);
/* L230: */
    }
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L240: */
	z__[i__] /= wght[i__];
    }
/* Load RHO into RHOK. */
    *rhok = rho;
    return 0;
/* ----------------------------------------------------------------------- */
/* This block handles error returns forced by routine PSOL. */
/* ----------------------------------------------------------------------- */
L300:
    if (ier < 0) {
	*iflag = -1;
    }
    if (ier > 0) {
	*iflag = 3;
    }

    return 0;

/* ------END OF SUBROUTINE DSPIGM----------------------------------------- */
} /* dspigm_ */

/* Subroutine */ int _daskr_datv_(integer *neq, real_number *y, real_number *tn,
	real_number *yprime, real_number *savr, real_number *v, real_number *wght,
	 real_number *yptem, Unknown_fp res, integer *ires, Unknown_fp psol, real_number *
	z__, real_number *vtem, real_number *wp, integer *iwp, real_number *cj,
	real_number *eplin, integer *ier, integer *nre, integer *npsl,
	real_number *rpar, integer *ipar)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;


/* ***BEGIN PROLOGUE  DATV */
/* ***DATE WRITTEN   890101   (YYMMDD) */
/* ***REVISION DATE  900926   (YYMMDD) */


/* ----------------------------------------------------------------------- */
/* ***DESCRIPTION */

/* This routine computes the product */

/*   Z = (D-inverse)*(P-inverse)*(dF/dY)*(D*V), */

/* where F(Y) = G(T, Y, CJ*(Y-A)), CJ is a scalar proportional to 1/H, */
/* and A involves the past history of Y.  The quantity CJ*(Y-A) is */
/* an approximation to the first derivative of Y and is stored */
/* in the array YPRIME.  Note that dF/dY = dG/dY + CJ*dG/dYPRIME. */

/* D is a diagonal scaling matrix, and P is the left preconditioning */
/* matrix.  V is assumed to have L2 norm equal to 1. */
/* The product is stored in Z and is computed by means of a */
/* difference quotient, a call to RES, and one call to PSOL. */

/*      On entry */

/*          NEQ = Problem size, passed to RES and PSOL. */

/*            Y = Array containing current dependent variable vector. */

/*       YPRIME = Array containing current first derivative of y. */

/*         SAVR = Array containing current value of G(T,Y,YPRIME). */

/*            V = Real array of length NEQ (can be the same array as Z). */

/*         WGHT = Array of length NEQ containing scale factors. */
/*                1/WGHT(I) are the diagonal elements of the matrix D. */

/*        YPTEM = Work array of length NEQ. */

/*         VTEM = Work array of length NEQ used to store the */
/*                unscaled version of V. */

/*         WP = Real work array used by preconditioner PSOL. */

/*         IWP = Integer work array used by preconditioner PSOL. */

/*           CJ = Scalar proportional to current value of */
/*                1/(step size H). */


/*      On return */

/*            Z = Array of length NEQ containing desired scaled */
/*                matrix-vector product. */

/*         IRES = Error flag from RES. */

/*          IER = Error flag from PSOL. */

/*         NRE  = The number of calls to RES. */

/*         NPSL = The number of calls to PSOL. */

/* ----------------------------------------------------------------------- */
/* ***ROUTINES CALLED */
/*   RES, PSOL */

/* ***END PROLOGUE  DATV */


    /* Parameter adjustments */
    --ipar;
    --rpar;
    --iwp;
    --wp;
    --vtem;
    --z__;
    --yptem;
    --wght;
    --v;
    --savr;
    --yprime;
    --y;

    /* Function Body */
    *ires = 0;
/* ----------------------------------------------------------------------- */
/* Set VTEM = D * V. */
/* ----------------------------------------------------------------------- */
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	vtem[i__] = v[i__] / wght[i__];
    }
    *ier = 0;
/* ----------------------------------------------------------------------- */
/* Store Y in Z and increment Z by VTEM. */
/* Store YPRIME in YPTEM and increment YPTEM by VTEM*CJ. */
/* ----------------------------------------------------------------------- */
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	yptem[i__] = yprime[i__] + vtem[i__] * *cj;
/* L20: */
	z__[i__] = y[i__] + vtem[i__];
    }
/* ----------------------------------------------------------------------- */
/* Call RES with incremented Y, YPRIME arguments */
/* stored in Z, YPTEM.  VTEM is overwritten with new residual. */
/* ----------------------------------------------------------------------- */
    ((int (*)(real_number *, real_number *, real_number*, real_number*, real_number *, integer *, real_number *, integer *))*res)
      (tn, &z__[1], &yptem[1], cj, &vtem[1], ires, &rpar[1], &ipar[1]);
    ++(*nre);
    if (*ires < 0) {
	return 0;
    }
/* ----------------------------------------------------------------------- */
/* Set Z = (dF/dY) * VBAR using difference quotient. */
/* (VBAR is old value of VTEM before calling RES) */
/* ----------------------------------------------------------------------- */
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L70: */
	z__[i__] = vtem[i__] - savr[i__];
    }
/* ----------------------------------------------------------------------- */
/* Apply inverse of left preconditioner to Z. */
/* ----------------------------------------------------------------------- */
    ((int (*)(integer*, real_number *, real_number *, real_number*, real_number *, real_number *, real_number *, real_number *, real_number *, integer *, real_number*, real_number*, integer *, real_number *, integer *))*psol)
      (neq, tn, &y[1], &yprime[1], &savr[1], &yptem[1], cj, &wght[1], &
	    wp[1], &iwp[1], &z__[1], eplin, ier, &rpar[1], &ipar[1]);
    ++(*npsl);
    if (*ier != 0) {
	return 0;
    }
/* ----------------------------------------------------------------------- */
/* Apply D-inverse to Z and return. */
/* ----------------------------------------------------------------------- */
    i__1 = *neq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L90: */
	z__[i__] *= wght[i__];
    }
    return 0;

/* ------END OF SUBROUTINE DATV------------------------------------------- */
} /* datv_ */

/* Subroutine */ int _daskr_dorth_(real_number *vnew, real_number *v, real_number *hes,
	integer *n, integer *ll, integer *ldhes, integer *kmp, real_number *
	snormw)
{
    /* System generated locals */
    integer v_dim1, v_offset, hes_dim1, hes_offset, i__1, i__2;
    real_number d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(real_number);

    /* Local variables */
    static integer i__, i0;
    static real_number arg, tem;
    extern real_number _daskr_ddot_(integer *, real_number *, integer *, real_number *,
	    integer *);
    static real_number vnrm;
    extern real_number _daskr_dnrm2_(integer *, real_number *, integer *);
    extern /* Subroutine */ int _daskr_daxpy_(integer *, real_number *, real_number *,
	    integer *, real_number *, integer *);
    static real_number sumdsq;


/* ***BEGIN PROLOGUE  DORTH */
/* ***DATE WRITTEN   890101   (YYMMDD) */
/* ***REVISION DATE  900926   (YYMMDD) */


/* ----------------------------------------------------------------------- */
/* ***DESCRIPTION */

/* This routine orthogonalizes the vector VNEW against the previous */
/* KMP vectors in the V array.  It uses a modified Gram-Schmidt */
/* orthogonalization procedure with conditional reorthogonalization. */

/*      On entry */

/*         VNEW = The vector of length N containing a scaled product */
/*                OF The Jacobian and the vector V(*,LL). */

/*         V    = The N x LL array containing the previous LL */
/*                orthogonal vectors V(*,1) to V(*,LL). */

/*         HES  = An LL x LL upper Hessenberg matrix containing, */
/*                in HES(I,K), K.LT.LL, scaled inner products of */
/*                A*V(*,K) and V(*,I). */

/*        LDHES = The leading dimension of the HES array. */

/*         N    = The order of the matrix A, and the length of VNEW. */

/*         LL   = The current order of the matrix HES. */

/*          KMP = The number of previous vectors the new vector VNEW */
/*                must be made orthogonal to (KMP .LE. MAXL). */


/*      On return */

/*         VNEW = The new vector orthogonal to V(*,I0), */
/*                where I0 = MAX(1, LL-KMP+1). */

/*         HES  = Upper Hessenberg matrix with column LL filled in with */
/*                scaled inner products of A*V(*,LL) and V(*,I). */

/*       SNORMW = L-2 norm of VNEW. */

/* ----------------------------------------------------------------------- */
/* ***ROUTINES CALLED */
/*   DDOT, DNRM2, DAXPY */

/* ***END PROLOGUE  DORTH */


/* ----------------------------------------------------------------------- */
/* Get norm of unaltered VNEW for later use. */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --vnew;
    v_dim1 = *n;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    hes_dim1 = *ldhes;
    hes_offset = 1 + hes_dim1;
    hes -= hes_offset;

    /* Function Body */
    vnrm = _daskr_dnrm2_(n, &vnew[1], &c__1);
/* ----------------------------------------------------------------------- */
/* Do Modified Gram-Schmidt on VNEW = A*V(LL). */
/* Scaled inner products give new column of HES. */
/* Projections of earlier vectors are subtracted from VNEW. */
/* ----------------------------------------------------------------------- */
/* Computing MAX */
    i__1 = 1, i__2 = *ll - *kmp + 1;
    i0 = MAX(i__1,i__2);
    i__1 = *ll;
    for (i__ = i0; i__ <= i__1; ++i__) {
	hes[i__ + *ll * hes_dim1] = _daskr_ddot_(n, &v[i__ * v_dim1 + 1], &c__1, &
		vnew[1], &c__1);
	tem = -hes[i__ + *ll * hes_dim1];
	_daskr_daxpy_(n, &tem, &v[i__ * v_dim1 + 1], &c__1, &vnew[1], &c__1);
/* L10: */
    }
/* ----------------------------------------------------------------------- */
/* Compute SNORMW = norm of VNEW. */
/* If VNEW is small compared to its input value (in norm), then */
/* Reorthogonalize VNEW to V(*,1) through V(*,LL). */
/* Correct if relative correction exceeds 1000*(unit roundoff). */
/* Finally, correct SNORMW using the dot products involved. */
/* ----------------------------------------------------------------------- */
    *snormw = _daskr_dnrm2_(n, &vnew[1], &c__1);
    if (vnrm + *snormw * .001 != vnrm) {
	return 0;
    }
    sumdsq = 0.;
    i__1 = *ll;
    for (i__ = i0; i__ <= i__1; ++i__) {
	tem = -_daskr_ddot_(n, &v[i__ * v_dim1 + 1], &c__1, &vnew[1], &c__1);
	if (hes[i__ + *ll * hes_dim1] + tem * .001 == hes[i__ + *ll * 
		hes_dim1]) {
	    goto L30;
	}
	hes[i__ + *ll * hes_dim1] -= tem;
	_daskr_daxpy_(n, &tem, &v[i__ * v_dim1 + 1], &c__1, &vnew[1], &c__1);
/* Computing 2nd power */
	d__1 = tem;
	sumdsq += d__1 * d__1;
L30:
	;
    }
    if (sumdsq == 0.) {
	return 0;
    }
/* Computing MAX */
/* Computing 2nd power */
    d__3 = *snormw;
    d__1 = 0., d__2 = d__3 * d__3 - sumdsq;
    arg = MAX(d__1,d__2);
    *snormw = sqrt(arg);
    return 0;

/* ------END OF SUBROUTINE DORTH------------------------------------------ */
} /* dorth_ */

/* Subroutine */ int _daskr_dheqr_(real_number *a, integer *lda, integer *n,
	real_number *q, integer *info, integer *ijob)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(real_number);

    /* Local variables */
    static real_number c__;
    static integer i__, j, k;
    static real_number s, t, t1, t2;
    static integer iq, km1, kp1, nm1;


/* ***BEGIN PROLOGUE  DHEQR */
/* ***DATE WRITTEN   890101   (YYMMDD) */
/* ***REVISION DATE  900926   (YYMMDD) */

/* ----------------------------------------------------------------------- */
/* ***DESCRIPTION */

/*     This routine performs a QR decomposition of an upper */
/*     Hessenberg matrix A.  There are two options available: */

/*          (1)  performing a fresh decomposition */
/*          (2)  updating the QR factors by adding a row and A */
/*               column to the matrix A. */

/*     DHEQR decomposes an upper Hessenberg matrix by using Givens */
/*     rotations. */

/*     On entry */

/*        A       DOUBLE PRECISION(LDA, N) */
/*                The matrix to be decomposed. */

/*        LDA     INTEGER */
/*                The leading dimension of the array A. */

/*        N       INTEGER */
/*                A is an (N+1) by N Hessenberg matrix. */

/*        IJOB    INTEGER */
/*                = 1     Means that a fresh decomposition of the */
/*                        matrix A is desired. */
/*                .GE. 2  Means that the current decomposition of A */
/*                        will be updated by the addition of a row */
/*                        and a column. */
/*     On return */

/*        A       The upper triangular matrix R. */
/*                The factorization can be written Q*A = R, where */
/*                Q is a product of Givens rotations and R is upper */
/*                triangular. */

/*        Q       DOUBLE PRECISION(2*N) */
/*                The factors C and S of each Givens rotation used */
/*                in decomposing A. */

/*        INFO    INTEGER */
/*                = 0  normal value. */
/*                = K  If  A(K,K) .EQ. 0.0.  This is not an error */
/*                     condition for this subroutine, but it does */
/*                     indicate that DHELS will divide by zero */
/*                     if called. */

/*     Modification of LINPACK. */
/*     Peter Brown, Lawrence Livermore Natl. Lab. */

/* ----------------------------------------------------------------------- */
/* ***ROUTINES CALLED (NONE) */

/* ***END PROLOGUE  DHEQR */


    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --q;

    /* Function Body */
    if (*ijob > 1) {
	goto L70;
    }
/* ----------------------------------------------------------------------- */
/* A new factorization is desired. */
/* ----------------------------------------------------------------------- */

/*     QR decomposition without pivoting. */

    *info = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	km1 = k - 1;
	kp1 = k + 1;

/*           Compute Kth column of R. */
/*           First, multiply the Kth column of A by the previous */
/*           K-1 Givens rotations. */

	if (km1 < 1) {
	    goto L20;
	}
	i__2 = km1;
	for (j = 1; j <= i__2; ++j) {
	    i__ = (j - 1 << 1) + 1;
	    t1 = a[j + k * a_dim1];
	    t2 = a[j + 1 + k * a_dim1];
	    c__ = q[i__];
	    s = q[i__ + 1];
	    a[j + k * a_dim1] = c__ * t1 - s * t2;
	    a[j + 1 + k * a_dim1] = s * t1 + c__ * t2;
/* L10: */
	}

/*           Compute Givens components C and S. */

L20:
	iq = (km1 << 1) + 1;
	t1 = a[k + k * a_dim1];
	t2 = a[kp1 + k * a_dim1];
	if (t2 != 0.) {
	    goto L30;
	}
	c__ = 1.;
	s = 0.;
	goto L50;
L30:
	if (fabs(t2) < fabs(t1)) {
	    goto L40;
	}
	t = t1 / t2;
	s = -1. / sqrt(t * t + 1.);
	c__ = -s * t;
	goto L50;
L40:
	t = t2 / t1;
	c__ = 1. / sqrt(t * t + 1.);
	s = -c__ * t;
L50:
	q[iq] = c__;
	q[iq + 1] = s;
	a[k + k * a_dim1] = c__ * t1 - s * t2;
	if (a[k + k * a_dim1] == 0.) {
	    *info = k;
	}
/* L60: */
    }
    return 0;
/* ----------------------------------------------------------------------- */
/* The old factorization of A will be updated.  A row and a column */
/* has been added to the matrix A. */
/* N by N-1 is now the old size of the matrix. */
/* ----------------------------------------------------------------------- */
L70:
    nm1 = *n - 1;
/* ----------------------------------------------------------------------- */
/* Multiply the new column by the N previous Givens rotations. */
/* ----------------------------------------------------------------------- */
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	i__ = (k - 1 << 1) + 1;
	t1 = a[k + *n * a_dim1];
	t2 = a[k + 1 + *n * a_dim1];
	c__ = q[i__];
	s = q[i__ + 1];
	a[k + *n * a_dim1] = c__ * t1 - s * t2;
	a[k + 1 + *n * a_dim1] = s * t1 + c__ * t2;
/* L100: */
    }
/* ----------------------------------------------------------------------- */
/* Complete update of decomposition by forming last Givens rotation, */
/* and multiplying it times the column vector (A(N,N),A(NP1,N)). */
/* ----------------------------------------------------------------------- */
    *info = 0;
    t1 = a[*n + *n * a_dim1];
    t2 = a[*n + 1 + *n * a_dim1];
    if (t2 != 0.) {
	goto L110;
    }
    c__ = 1.;
    s = 0.;
    goto L130;
L110:
    if (fabs(t2) < fabs(t1)) {
	goto L120;
    }
    t = t1 / t2;
    s = -1. / sqrt(t * t + 1.);
    c__ = -s * t;
    goto L130;
L120:
    t = t2 / t1;
    c__ = 1. / sqrt(t * t + 1.);
    s = -c__ * t;
L130:
    iq = (*n << 1) - 1;
    q[iq] = c__;
    q[iq + 1] = s;
    a[*n + *n * a_dim1] = c__ * t1 - s * t2;
    if (a[*n + *n * a_dim1] == 0.) {
	*info = *n;
    }
    return 0;

/* ------END OF SUBROUTINE DHEQR------------------------------------------ */
} /* dheqr_ */

/* Subroutine */ int _daskr_dhels_(real_number *a, integer *lda, integer *n,
	real_number *q, real_number *b)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static real_number c__;
    static integer k;
    static real_number s, t, t1, t2;
    static integer kb, iq, kp1;
    extern /* Subroutine */ int _daskr_daxpy_(integer *, real_number *, real_number *,
	    integer *, real_number *, integer *);


/* ***BEGIN PROLOGUE  DHELS */
/* ***DATE WRITTEN   890101   (YYMMDD) */
/* ***REVISION DATE  900926   (YYMMDD) */


/* ----------------------------------------------------------------------- */
/* ***DESCRIPTION */

/* This is similar to the LINPACK routine DGESL except that */
/* A is an upper Hessenberg matrix. */

/*     DHELS solves the least squares problem */

/*           MIN (B-A*X,B-A*X) */

/*     using the factors computed by DHEQR. */

/*     On entry */

/*        A       DOUBLE PRECISION (LDA, N) */
/*                The output from DHEQR which contains the upper */
/*                triangular factor R in the QR decomposition of A. */

/*        LDA     INTEGER */
/*                The leading dimension of the array  A . */

/*        N       INTEGER */
/*                A is originally an (N+1) by N matrix. */

/*        Q       DOUBLE PRECISION(2*N) */
/*                The coefficients of the N givens rotations */
/*                used in the QR factorization of A. */

/*        B       DOUBLE PRECISION(N+1) */
/*                The right hand side vector. */


/*     On return */

/*        B       The solution vector X. */


/*     Modification of LINPACK. */
/*     Peter Brown, Lawrence Livermore Natl. Lab. */

/* ----------------------------------------------------------------------- */
/* ***ROUTINES CALLED */
/*   DAXPY */

/* ***END PROLOGUE  DHELS */


/*        Minimize (B-A*X,B-A*X). */
/*        First form Q*B. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --q;
    --b;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;
	iq = (k - 1 << 1) + 1;
	c__ = q[iq];
	s = q[iq + 1];
	t1 = b[k];
	t2 = b[kp1];
	b[k] = c__ * t1 - s * t2;
	b[kp1] = s * t1 + c__ * t2;
/* L20: */
    }

/*        Now solve R*X = Q*B. */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	b[k] /= a[k + k * a_dim1];
	t = -b[k];
	i__2 = k - 1;
	_daskr_daxpy_(&i__2, &t, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
/* L40: */
    }
    return 0;

/* ------END OF SUBROUTINE DHELS------------------------------------------ */
} /* dhels_ */

