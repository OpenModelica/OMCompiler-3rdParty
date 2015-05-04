/* dsparsk.f -- translated by f2c (version 20100827).
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../solver/ddaskr_types.h"


/* ----------------------------------------------------------------------c */
/*                          S P A R S K I T                             c */
/* ----------------------------------------------------------------------c */
/*        BASIC LINEAR ALGEBRA FOR SPARSE MATRICES. BLASSM MODULE       c */
/* ----------------------------------------------------------------------c */
/* aplb   :   computes     C = A+B                                      c */
/* aplb1  :   computes     C = A+B  [Sorted version: A, B, C sorted]    c */
/* aplsb  :   computes     C = A + s B                                  c */
/* diamua :   Computes     C = Diag * A                                 c */
/* amudia :   Computes     C = A* Diag                                  c */
/* aplsca :   Computes     A:= A + s I    (s = scalar)                  c */
/* ----------------------------------------------------------------------c */
/* Subroutine */ int _daskr_aplb_(integer *nrow, integer *ncol, integer *job,
	real_number *a, integer *ja, integer *ia, real_number *b, integer *jb,
	integer *ib, real_number *c__, integer *jc, integer *ic, integer *
	nzmax, integer *iw, integer *ierr)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j, k, ka, kb, ii, len, jcol, jpos;
    static integer values;

/* ----------------------------------------------------------------------- */
/* performs the matrix sum  C = A+B. */
/* ----------------------------------------------------------------------- */
/* on entry: */
/* --------- */
/* nrow	= integer. The row dimension of A and B */
/* ncol  = integer. The column dimension of A and B. */
/* job   = integer. Job indicator. When job = 0, only the structure */
/*                  (i.e. the arrays jc, ic) is computed and the */
/*                  real values are ignored. */

/* a, */
/* ja, */
/* ia   = Matrix A in compressed sparse row format. */

/* b, */
/* jb, */
/* ib	=  Matrix B in compressed sparse row format. */

/* nzmax	= integer. The  length of the arrays c and jc. */
/*         amub will stop if the result matrix C  has a number */
/*         of elements that exceeds exceeds nzmax. See ierr. */

/* on return: */
/* ---------- */
/* c, */
/* jc, */
/* ic	= resulting matrix C in compressed sparse row sparse format. */

/* ierr	= integer. serving as error message. */
/*         ierr = 0 means normal return, */
/*         ierr .gt. 0 means that amub stopped while computing the */
/*         i-th row  of C with i=ierr, because the number */
/*         of elements in C exceeds nzmax. */

/* work arrays: */
/* ------------ */
/* iw	= integer work array of length equal to the number of */
/*         columns in A. */

/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --ic;
    --ib;
    --ia;
    --iw;
    --a;
    --ja;
    --b;
    --jb;
    --c__;
    --jc;

    /* Function Body */
    values = *job != 0;
    *ierr = 0;
    len = 0;
    ic[1] = 1;
    i__1 = *ncol;
    for (j = 1; j <= i__1; ++j) {
	iw[j] = 0;
/* L1: */
    }

    i__1 = *nrow;
    for (ii = 1; ii <= i__1; ++ii) {
/*     row i */
	i__2 = ia[ii + 1] - 1;
	for (ka = ia[ii]; ka <= i__2; ++ka) {
	    ++len;
	    jcol = ja[ka];
	    if (len > *nzmax) {
		goto L999;
	    }
	    jc[len] = jcol;
	    if (values) {
		c__[len] = a[ka];
	    }
	    iw[jcol] = len;
/* L200: */
	}

	i__2 = ib[ii + 1] - 1;
	for (kb = ib[ii]; kb <= i__2; ++kb) {
	    jcol = jb[kb];
	    jpos = iw[jcol];
	    if (jpos == 0) {
		++len;
		if (len > *nzmax) {
		    goto L999;
		}
		jc[len] = jcol;
		if (values) {
		    c__[len] = b[kb];
		}
		iw[jcol] = len;
	    } else {
		if (values) {
		    c__[jpos] += b[kb];
		}
	    }
/* L300: */
	}
	i__2 = len;
	for (k = ic[ii]; k <= i__2; ++k) {
	    iw[jc[k]] = 0;
/* L301: */
	}
	ic[ii + 1] = len + 1;
/* L500: */
    }
    return 0;
L999:
    *ierr = ii;
    return 0;
/* ------------end of aplb ----------------------------------------------- */
/* ----------------------------------------------------------------------- */
} /* aplb_ */

/* Subroutine */ int _daskr_aplb1_(integer *nrow, integer *ncol, integer *job,
	real_number *a, integer *ja, integer *ia, real_number *b, integer *jb,
	integer *ib, real_number *c__, integer *jc, integer *ic, integer *
	nzmax, integer *ierr)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j1, j2, kc, ka, kb, kamax, kbmax;
    static integer values;

/* ----------------------------------------------------------------------- */
/* performs the matrix sum  C = A+B for matrices in sorted CSR format. */
/* the difference with aplb  is that the resulting matrix is such that */
/* the elements of each row are sorted with increasing column indices in */
/* each row, provided the original matrices are sorted in the same way. */
/* ----------------------------------------------------------------------- */
/* on entry: */
/* --------- */
/* nrow	= integer. The row dimension of A and B */
/* ncol  = integer. The column dimension of A and B. */
/* job   = integer. Job indicator. When job = 0, only the structure */
/*                  (i.e. the arrays jc, ic) is computed and the */
/*                  real values are ignored. */

/* a, */
/* ja, */
/* ia   = Matrix A in compressed sparse row format with entries sorted */

/* b, */
/* jb, */
/* ib	=  Matrix B in compressed sparse row format with entries sorted */
/*        ascendly in each row */

/* nzmax	= integer. The  length of the arrays c and jc. */
/*         amub will stop if the result matrix C  has a number */
/*         of elements that exceeds exceeds nzmax. See ierr. */

/* on return: */
/* ---------- */
/* c, */
/* jc, */
/* ic	= resulting matrix C in compressed sparse row sparse format */
/*         with entries sorted ascendly in each row. */

/* ierr	= integer. serving as error message. */
/*         ierr = 0 means normal return, */
/*         ierr .gt. 0 means that amub stopped while computing the */
/*         i-th row  of C with i=ierr, because the number */
/*         of elements in C exceeds nzmax. */

/* Notes: */
/* ------- */
/*     this will not work if any of the two input matrices is not sorted */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --ic;
    --ib;
    --ia;
    --a;
    --ja;
    --b;
    --jb;
    --c__;
    --jc;

    /* Function Body */
    values = *job != 0;
    *ierr = 0;
    kc = 1;
    ic[1] = kc;

    i__1 = *nrow;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ka = ia[i__];
	kb = ib[i__];
	kamax = ia[i__ + 1] - 1;
	kbmax = ib[i__ + 1] - 1;
L5:
	if (ka <= kamax) {
	    j1 = ja[ka];
	} else {
	    j1 = *ncol + 1;
	}
	if (kb <= kbmax) {
	    j2 = jb[kb];
	} else {
	    j2 = *ncol + 1;
	}

/*     three cases */

	if (j1 == j2) {
	    if (values) {
		c__[kc] = a[ka] + b[kb];
	    }
	    jc[kc] = j1;
	    ++ka;
	    ++kb;
	    ++kc;
	} else if (j1 < j2) {
	    jc[kc] = j1;
	    if (values) {
		c__[kc] = a[ka];
	    }
	    ++ka;
	    ++kc;
	} else if (j1 > j2) {
	    jc[kc] = j2;
	    if (values) {
		c__[kc] = b[kb];
	    }
	    ++kb;
	    ++kc;
	}
	if (kc > *nzmax) {
	    goto L999;
	}
	if (ka <= kamax || kb <= kbmax) {
	    goto L5;
	}
	ic[i__ + 1] = kc;
/* L6: */
    }
    return 0;
L999:
    *ierr = i__;
    return 0;
/* ------------end-of-aplb1----------------------------------------------- */
/* ----------------------------------------------------------------------- */
} /* aplb1_ */

/* Subroutine */ int _daskr_aplsb_(integer *nrow, integer *ncol, real_number *a,
	integer *ja, integer *ia, real_number *s, real_number *b, integer *jb,
	integer *ib, real_number *c__, integer *jc, integer *ic, integer *
	nzmax, integer *ierr)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j1, j2, kc, ka, kb, kamax, kbmax;

/* ----------------------------------------------------------------------- */
/* performs the operation C = A+s B for matrices in sorted CSR format. */
/* the difference with aplsb is that the resulting matrix is such that */
/* the elements of each row are sorted with increasing column indices in */
/* each row, provided the original matrices are sorted in the same way. */
/* ----------------------------------------------------------------------- */
/* on entry: */
/* --------- */
/* nrow	= integer. The row dimension of A and B */
/* ncol  = integer. The column dimension of A and B. */

/* a, */
/* ja, */
/* ia   = Matrix A in compressed sparse row format with entries sorted */

/* s	= real. scalar factor for B. */

/* b, */
/* jb, */
/* ib	=  Matrix B in compressed sparse row format with entries sorted */
/*        ascendly in each row */

/* nzmax	= integer. The  length of the arrays c and jc. */
/*         amub will stop if the result matrix C  has a number */
/*         of elements that exceeds exceeds nzmax. See ierr. */

/* on return: */
/* ---------- */
/* c, */
/* jc, */
/* ic	= resulting matrix C in compressed sparse row sparse format */
/*         with entries sorted ascendly in each row. */

/* ierr	= integer. serving as error message. */
/*         ierr = 0 means normal return, */
/*         ierr .gt. 0 means that amub stopped while computing the */
/*         i-th row  of C with i=ierr, because the number */
/*         of elements in C exceeds nzmax. */

/* Notes: */
/* ------- */
/*     this will not work if any of the two input matrices is not sorted */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --ic;
    --ib;
    --ia;
    --a;
    --ja;
    --b;
    --jb;
    --c__;
    --jc;

    /* Function Body */
    *ierr = 0;
    kc = 1;
    ic[1] = kc;

/*     the following loop does a merge of two sparse rows + adds  them. */

    i__1 = *nrow;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ka = ia[i__];
	kb = ib[i__];
	kamax = ia[i__ + 1] - 1;
	kbmax = ib[i__ + 1] - 1;
L5:

/*     this is a while  -- do loop -- */

	if (ka <= kamax || kb <= kbmax) {

	    if (ka <= kamax) {
		j1 = ja[ka];
	    } else {
/*     take j1 large enough  that always j2 .lt. j1 */
		j1 = *ncol + 1;
	    }
	    if (kb <= kbmax) {
		j2 = jb[kb];
	    } else {
/*     similarly take j2 large enough  that always j1 .lt. j2 */
		j2 = *ncol + 1;
	    }

/*     three cases */

	    if (j1 == j2) {
		c__[kc] = a[ka] + *s * b[kb];
		jc[kc] = j1;
		++ka;
		++kb;
		++kc;
	    } else if (j1 < j2) {
		jc[kc] = j1;
		c__[kc] = a[ka];
		++ka;
		++kc;
	    } else if (j1 > j2) {
		jc[kc] = j2;
		c__[kc] = *s * b[kb];
		++kb;
		++kc;
	    }
	    if (kc > *nzmax) {
		goto L999;
	    }
	    goto L5;

/*     end while loop */

	}
	ic[i__ + 1] = kc;
/* L6: */
    }
    return 0;
L999:
    *ierr = i__;
    return 0;
/* ------------end-of-aplsb --------------------------------------------- */
/* ----------------------------------------------------------------------- */
} /* aplsb_ */

/* Subroutine */ int _daskr_diamua_(integer *nrow, integer *job, real_number *a,
	integer *ja, integer *ia, real_number *diag, real_number *b, integer *
	jb, integer *ib)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer k, k1, k2, ii;
    static real_number scal;

/* ----------------------------------------------------------------------- */
/* performs the matrix by matrix product B = Diag * A  (in place) */
/* ----------------------------------------------------------------------- */
/* on entry: */
/* --------- */
/* nrow	= integer. The row dimension of A */

/* job   = integer. job indicator. Job=0 means get array b only */
/*         job = 1 means get b, and the integer arrays ib, jb. */

/* a, */
/* ja, */
/* ia   = Matrix A in compressed sparse row format. */

/* diag = diagonal matrix stored as a vector dig(1:n) */

/* on return: */
/* ---------- */

/* b, */
/* jb, */
/* ib	= resulting matrix B in compressed sparse row sparse format. */

/* Notes: */
/* ------- */
/* 1)        The column dimension of A is not needed. */
/* 2)        algorithm in place (B can take the place of A). */
/*           in this case use job=0. */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --ib;
    --diag;
    --ia;
    --a;
    --ja;
    --b;
    --jb;

    /* Function Body */
    i__1 = *nrow;
    for (ii = 1; ii <= i__1; ++ii) {

/*     normalize each row */

	k1 = ia[ii];
	k2 = ia[ii + 1] - 1;
	scal = diag[ii];
	i__2 = k2;
	for (k = k1; k <= i__2; ++k) {
	    b[k] = a[k] * scal;
/* L2: */
	}
/* L1: */
    }

    if (*job == 0) {
	return 0;
    }

    i__1 = *nrow + 1;
    for (ii = 1; ii <= i__1; ++ii) {
	ib[ii] = ia[ii];
/* L3: */
    }
    i__1 = ia[*nrow + 1] - 1;
    for (k = ia[1]; k <= i__1; ++k) {
	jb[k] = ja[k];
/* L31: */
    }
    return 0;
/* ----------end-of-diamua------------------------------------------------ */
/* ----------------------------------------------------------------------- */
} /* diamua_ */

/* Subroutine */ int _daskr_amudia_(integer *nrow, integer *job, real_number *a,
	integer *ja, integer *ia, real_number *diag, real_number *b, integer *
	jb, integer *ib)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer k, k1, k2, ii;

/* ----------------------------------------------------------------------- */
/* performs the matrix by matrix product B = A * Diag  (in place) */
/* ----------------------------------------------------------------------- */
/* on entry: */
/* --------- */
/* nrow	= integer. The row dimension of A */

/* job   = integer. job indicator. Job=0 means get array b only */
/*         job = 1 means get b, and the integer arrays ib, jb. */

/* a, */
/* ja, */
/* ia   = Matrix A in compressed sparse row format. */

/* diag = diagonal matrix stored as a vector dig(1:n) */

/* on return: */
/* ---------- */

/* b, */
/* jb, */
/* ib	= resulting matrix B in compressed sparse row sparse format. */

/* Notes: */
/* ------- */
/* 1)        The column dimension of A is not needed. */
/* 2)        algorithm in place (B can take the place of A). */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --ib;
    --diag;
    --ia;
    --a;
    --ja;
    --b;
    --jb;

    /* Function Body */
    i__1 = *nrow;
    for (ii = 1; ii <= i__1; ++ii) {

/*     scale each element */

	k1 = ia[ii];
	k2 = ia[ii + 1] - 1;
	i__2 = k2;
	for (k = k1; k <= i__2; ++k) {
	    b[k] = a[k] * diag[ja[k]];
/* L2: */
	}
/* L1: */
    }

    if (*job == 0) {
	return 0;
    }

    i__1 = *nrow + 1;
    for (ii = 1; ii <= i__1; ++ii) {
	ib[ii] = ia[ii];
/* L3: */
    }
    i__1 = ia[*nrow + 1] - 1;
    for (k = ia[1]; k <= i__1; ++k) {
	jb[k] = ja[k];
/* L31: */
    }
    return 0;
/* ----------------------------------------------------------------------- */
/* -----------end-of-amudiag---------------------------------------------- */
} /* amudia_ */

/* Subroutine */ int _daskr_aplsca_(integer *nrow, real_number *a, integer *ja,
	integer *ia, real_number *scal, integer *iw)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, k, k1, k2, ii, ko;
    static integer test;
    extern /* Subroutine */ int _daskr_diapos_(integer *, integer *, integer *,
	    integer *);
    static integer icount;

/* ----------------------------------------------------------------------- */
/* Adds a scalar to the diagonal entries of a sparse matrix A :=A + s I */
/* ----------------------------------------------------------------------- */
/* on entry: */
/* --------- */
/* nrow	= integer. The row dimension of A */

/* a, */
/* ja, */
/* ia    = Matrix A in compressed sparse row format. */

/* scal  = real. scalar to add to the diagonal entries. */

/* on return: */
/* ---------- */

/* a, */
/* ja, */
/* ia	= matrix A with diagonal elements shifted (or created). */

/* iw    = integer work array of length n. On return iw will */
/*         contain  the positions of the diagonal entries in the */
/*         output matrix. (i.e., a(iw(k)), ja(iw(k)), k=1,...n, */
/*         are the values/column indices of the diagonal elements */
/*         of the output matrix. ). */

/* Notes: */
/* ------- */
/*     The column dimension of A is not needed. */
/*     important: the matrix a may be expanded slightly to allow for */
/*     additions of nonzero elements to previously nonexisting diagonals. */
/*     The is no checking as to whether there is enough space appended */
/*     to the arrays a and ja. if not sure allow for n additional */
/*     elemnts. */
/*     coded by Y. Saad. Latest version July, 19, 1990 */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --ia;
    --a;
    --ja;
    --iw;

    /* Function Body */
    _daskr_diapos_(nrow, &ja[1], &ia[1], &iw[1]);
    icount = 0;
    i__1 = *nrow;
    for (j = 1; j <= i__1; ++j) {
	if (iw[j] == 0) {
	    ++icount;
	} else {
	    a[iw[j]] += *scal;
	}
/* L1: */
    }

/*     if no diagonal elements to insert in data structure return. */

    if (icount == 0) {
	return 0;
    }

/* shift the nonzero elements if needed, to allow for created */
/* diagonal elements. */

    ko = ia[*nrow + 1] + icount;

/*     copy rows backward */

    for (ii = *nrow; ii >= 1; --ii) {

/*     go through  row ii */

	k1 = ia[ii];
	k2 = ia[ii + 1] - 1;
	ia[ii + 1] = ko;
	test = iw[ii] == 0;
	i__1 = k1;
	for (k = k2; k >= i__1; --k) {
	    j = ja[k];
	    if (test && j < ii) {
		test = _FALSE_;
		--ko;
		a[ko] = *scal;
		ja[ko] = ii;
		iw[ii] = ko;
	    }
	    --ko;
	    a[ko] = a[k];
	    ja[ko] = j;
/* L4: */
	}
/*     diagonal element has not been added yet. */
	if (test) {
	    --ko;
	    a[ko] = *scal;
	    ja[ko] = ii;
	    iw[ii] = ko;
	}
/* L5: */
    }
    ia[1] = ko;
    return 0;
/* ----------------------------------------------------------------------- */
/* ----------end-of-aplsca------------------------------------------------ */
} /* aplsca_ */

/* ----------------------------------------------------------------------c */
/*                          S P A R S K I T                             c */
/* ----------------------------------------------------------------------c */
/*          BASIC MATRIX-VECTOR OPERATIONS - MATVEC MODULE              c */
/* ----------------------------------------------------------------------c */
/* amux  : A times a vector. Compressed Sparse Row (CSR) format.        c */
/* ----------------------------------------------------------------------c */
/* Subroutine */ int _daskr_amux_(integer *n, real_number *x, real_number *y,
	real_number *a, integer *ja, integer *ia)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, k;
    static real_number t;

/* ----------------------------------------------------------------------- */
/*         A times a vector */
/* ----------------------------------------------------------------------- */
/* multiplies a matrix by a vector using the dot product form */
/* Matrix A is stored in compressed sparse row storage. */

/* on entry: */
/* ---------- */
/* n     = row dimension of A */
/* x     = real array of length equal to the column dimension of */
/*         the A matrix. */
/* a, ja, */
/*    ia = input matrix in compressed sparse row format. */

/* on return: */
/* ----------- */
/* y     = real array of length n, containing the product y=Ax */

/* ----------------------------------------------------------------------- */
/* local variables */

/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --ia;
    --ja;
    --a;
    --y;
    --x;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {

/*     compute the inner product of row i with vector x */

	t = 0.;
	i__2 = ia[i__ + 1] - 1;
	for (k = ia[i__]; k <= i__2; ++k) {
	    t += a[k] * x[ja[k]];
/* L99: */
	}

/*     store result in y(i) */

	y[i__] = t;
/* L100: */
    }

    return 0;
/* ---------end-of-amux--------------------------------------------------- */
/* ----------------------------------------------------------------------- */
} /* amux_ */

/* ----------------------------------------------------------------------c */
/*                          S P A R S K I T                             c */
/* ----------------------------------------------------------------------c */
/*                        INPUT-OUTPUT MODULE                           c */
/* ----------------------------------------------------------------------c */
/*  prtmt  : prints matrices in the Boeing/Harwell format.              c */
/* ----------------------------------------------------------------------c */

/* Feature writing Jacobian Matrix to file is deactivated 				  */

/* Subroutine */ int _daskr_prtmt_(integer *nrow, integer *ncol, real_number *a,
	integer *ja, integer *ia, real_number *rhs, char *guesol, char *title,
	char *key, char *type__, integer *ifmt, integer *job, integer *iounit,
	 integer guesol_len, integer title_len, integer key_len, integer type_len)
{
    /* Format strings */
    static char fmt_101[] = "(\002(\002,i2,\002I\002,i2,\002)\002)";
    static char fmt_100[] = "(\002(\002,i2,\002I\002,i1,\002)\002)";
    static char fmt_100_new[] = "  (%2d I %1d)";
    static char fmt_101_new[] = "  (%2d I %2d)";
    static char fmt_102[] = "(\002(\002,i2,\002F\002,i1,\002.\002,i1,\002"
	    ")\002)";
    static char fmt_103[] = "(\002(\002,i2,\002F\002,i2,\002.\002,i1,\002"
	    ")\002)";
    static char fmt_104[] = "(\002(\002,i2,\002F\002,i2,\002.\002,i2,\002"
	    ")\002)";
    static char fmt_105[] = "(\002(\002,i1,\002D\002,i1,\002.\002,i1,\002"
	    ")\002)";
    static char fmt_106[] = "(\002(\002,i1,\002D\002,i2,\002.\002,i1,\002"
	    ")\002)";
    static char fmt_107[] = "(\002(\002,i1,\002D\002,i2,\002.\002,i2,\002"
	    ")\002)";
    static char fmt_108[] = "(\002(\002,i2,\002D\002,i1,\002.\002,i1,\002"
	    ")\002)";
    static char fmt_109[] = "(\002(\002,i2,\002D\002,i2,\002.\002,i1,\002"
	    ")\002)";
    static char fmt_110[] = "(\002(\002,i2,\002D\002,i2,\002.\002,i2,\002"
	    ")\002)";
    static char fmt_10[] = "(a72,a8/5i14/a3,11x,4i14/2a16,2a20)";
    static char fmt_11[] = "(a3,11x,i14,i14)";

    /* System generated locals */
    char a__1[2];
    integer i__1, i__2[2];
    real_number r__1;
    /* deactivate printing */
    /* icilist ici__1; */

    /* Builtin functions
    double r_lg10(real *);
    integer s_wsfi(icilist *), do_fio(integer *, char *, integer), e_wsfi(void)
	    ;
	*/
    /* Subroutine  int s_cat(char *, char **, integer *, integer *, integer);
    integer s_wsfe(cilist *), e_wsfe(void); */

    /* Local variables */
    static integer i__, ix, len, nnz, iend, nrhs, next, ihead, indcrd, valcrd;
    static char indfmt[16];
    static integer rhscrd;
    static char valfmt[20];
    static integer nperli, totcrd, ptrcrd;
    static char ptrfmt[16], rhstyp[3];
    static integer nrwindx;


    /* Fortran I/O blocks */
    /*
    static icilist io___55 = { 0, indfmt, 0, fmt_100, 16, 1 };
    static cilist io___64 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___66 = { 0, 0, 0, fmt_11, 0 };
    static cilist io___67 = { 0, 0, 0, ptrfmt, 0 };
    static cilist io___68 = { 0, 0, 0, indfmt, 0 };
    static cilist io___69 = { 0, 0, 0, valfmt, 0 };
    static cilist io___72 = { 0, 0, 0, valfmt, 0 };
    static cilist io___73 = { 0, 0, 0, valfmt, 0 };
    static cilist io___74 = { 0, 0, 0, valfmt, 0 };
    */


    /* Assigned format variables */
    static char *ix_fmt;

/* ----------------------------------------------------------------------- */
/* writes a matrix in Harwell-Boeing format into a file. */
/* assumes that the matrix is stored in COMPRESSED SPARSE COLUMN FORMAT. */
/* some limited functionality for right hand sides. */
/* Author: Youcef Saad - Date: Sept., 1989 - updated Oct. 31, 1989 to */
/* cope with new format. */
/* ----------------------------------------------------------------------- */
/* on entry: */
/* --------- */
/* nrow   = number of rows in matrix */
/* ncol	 = number of columns in matrix */
/* a	 = real*8 array containing the values of the matrix stored */
/*          columnwise */
/* ja 	 = integer array of the same length as a containing the column */
/*          indices of the corresponding matrix elements of array a. */
/* ia     = integer array of containing the pointers to the beginning of */
/* 	   the row in arrays a and ja. */
/* rhs    = real array  containing the right-hand-side (s) and optionally */
/*          the associated initial guesses and/or exact solutions */
/*          in this order. See also guesol for details. the vector rhs will */
/*          be used only if job .gt. 2 (see below). Only full storage for */
/*          the right hand sides is supported. */

/* guesol = a 2-character string indicating whether an initial guess */
/*          (1-st character) and / or the exact solution (2-nd) */
/*          character) is provided with the right hand side. */
/* 	   if the first character of guesol is 'G' it means that an */
/*          an intial guess is provided for each right-hand sides. */
/*          These are assumed to be appended to the right hand-sides in */
/*          the array rhs. */
/* 	   if the second character of guesol is 'X' it means that an */
/*          exact solution is provided for each right-hand side. */
/*          These are assumed to be appended to the right hand-sides */
/*          and the initial guesses (if any) in the array rhs. */

/* title  = character*72 = title of matrix test ( character a*72 ). */
/* key    = character*8  = key of matrix */
/* type   = charatcer*3  = type of matrix. */

/* ifmt	 = integer specifying the format chosen for the real values */
/* 	   to be output (i.e., for a, and for rhs-guess-sol if */
/*          applicable). The meaning of ifmt is as follows. */
/* 	  * if (ifmt .lt. 100) then the D descriptor is used, */
/*           format Dd.m, in which the length (m) of the mantissa is */
/*           precisely the integer ifmt (and d = ifmt+6) */
/* 	  * if (ifmt .gt. 100) then prtmt will use the */
/*           F- descriptor (format Fd.m) in which the length of the */
/*           mantissa (m) is the integer mod(ifmt,100) and the length */
/*           of the integer part is k=ifmt/100 (and d = k+m+2) */
/* 	    Thus  ifmt= 4   means  D10.4  +.xxxxD+ee    while */
/* 	          ifmt=104  means  F7.4   +x.xxxx */
/* 	          ifmt=205  means  F9.5   +xx.xxxxx */
/* 	    Note: formats for ja, and ia are internally computed. */

/* job	 = integer to indicate whether matrix values and */
/* 	   a right-hand-side is available to be written */
/*          job = 1   write srtucture only, i.e., the arrays ja and ia. */
/*          job = 2   write matrix including values, i.e., a, ja, ia */
/*          job = 3   write matrix and one right hand side: a,ja,ia,rhs. */
/* 	   job = nrhs+2 write matrix and nrhs successive right hand sides */
/* 	   Note that there cannot be any right-hand-side if the matrix */
/* 	   has no values. Also the initial guess and exact solutions when */
/*          provided are for each right hand side. For example if nrhs=2 */
/*          and guesol='GX' there are 6 vectors to write. */


/* iounit = logical unit number where to write the matrix into. */

/* on return: */
/* ---------- */
/* the matrix a, ja, ia will be written in output unit iounit */
/* in the Harwell-Boeing format. None of the inputs is modofied. */

/* Notes: 1) This code attempts to pack as many elements as possible per */
/*        80-character line. */
/*        2) this code attempts to avoid as much as possible to put */
/*        blanks in the formats that are written in the 4-line header */
/* 	 (This is done for purely esthetical reasons since blanks */
/*        are ignored in format descriptors.) */
/*        3) sparse formats for right hand sides and guesses are not */
/*        supported. */
/* ----------------------------------------------------------------------- */
/* -------------- */
/*     compute pointer format */
/* -------------- */
    /* Parameter adjustments */
    --rhs;
    --ia;
    --ja;
    --a;

    /* Function Body */
    nnz = ia[*ncol + 1] - 1;
    r__1 = (real_number) (nnz + 1) + .1f;
    len = (integer) log10(r__1) + 1;
    nperli = 80 / len;
    ptrcrd = *ncol / nperli + 1;
    if (len > 9) {
	ix = 0;
	ix_fmt = fmt_101;
    } else {
	ix = 1;
	ix_fmt = fmt_100;
    }
    /*
    ici__1.icierr = 0;
    ici__1.icirnum = 1;
    ici__1.icirlen = 16;
    ici__1.iciunit = ptrfmt;
    ici__1.icifmt = ix_fmt;
    s_wsfi(&ici__1);
    do_fio(&c__1, (char *)&nperli, (integer)sizeof(integer));
    do_fio(&c__1, (char *)&len, (integer)sizeof(integer));
    e_wsfi();
    */
L100:
L101:
/* ---------------------------- */
/* compute ROW index format */
/* ---------------------------- */
    r__1 = (real_number) (*nrow) + .1f;
    len = (integer) log10(r__1) + 1;
/* Computing MIN */
    i__1 = 80 / len;
    nperli = MIN(i__1,nnz);
    indcrd = (nnz - 1) / nperli + 1;
    /*
    s_wsfi(&io___55);
    do_fio(&c__1, (char *)&nperli, (integer)sizeof(integer));
    do_fio(&c__1, (char *)&len, (integer)sizeof(integer));
    e_wsfi();
    */
/* --------------- */
/* compute values and rhs format (using the same for both) */
/* --------------- */
    valcrd = 0;
    rhscrd = 0;
/* quit this part if no values provided. */
    if (*job <= 1) {
	goto L20;
    }

    if (*ifmt >= 100) {
	ihead = *ifmt / 100;
	*ifmt -= ihead * 100;
	len = ihead + *ifmt + 2;
	nperli = 80 / len;

	if (len <= 9) {
	    ix = 2;
	    ix_fmt = fmt_102;
	} else if (*ifmt <= 9) {
	    ix = 3;
	    ix_fmt = fmt_103;
	} else {
	    ix = 4;
	    ix_fmt = fmt_104;
	}
	/*
	ici__1.icierr = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 20;
	ici__1.iciunit = valfmt;
	ici__1.icifmt = ix_fmt;
	s_wsfi(&ici__1);
	do_fio(&c__1, (char *)&nperli, (integer)sizeof(integer));
	do_fio(&c__1, (char *)&len, (integer)sizeof(integer));
	do_fio(&c__1, (char *)&(*ifmt), (integer)sizeof(integer));
	e_wsfi();
	*/
L102:
L103:
L104:

	;
    } else {
	len = *ifmt + 6;
	nperli = 80 / len;
/*     try to minimize the blanks in the format strings. */
	if (nperli <= 9) {
	    if (len <= 9) {
		ix = 5;
		ix_fmt = fmt_105;
	    } else if (*ifmt <= 9) {
		ix = 6;
		ix_fmt = fmt_106;
	    } else {
		ix = 7;
		ix_fmt = fmt_107;
	    }
	} else {
	    if (len <= 9) {
		ix = 8;
		ix_fmt = fmt_108;
	    } else if (*ifmt <= 9) {
		ix = 9;
		ix_fmt = fmt_109;
	    } else {
		ix = 10;
		ix_fmt = fmt_110;
	    }
	}
/* ----------- */
	/*
	ici__1.icierr = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 20;
	ici__1.iciunit = valfmt;
	ici__1.icifmt = ix_fmt;
	s_wsfi(&ici__1);
	do_fio(&c__1, (char *)&nperli, (integer)sizeof(integer));
	do_fio(&c__1, (char *)&len, (integer)sizeof(integer));
	do_fio(&c__1, (char *)&(*ifmt), (integer)sizeof(integer));
	e_wsfi();
	*/
L105:
L106:
L107:
L108:
L109:
L110:

	;
    }
    valcrd = (nnz - 1) / nperli + 1;
    nrhs = *job - 2;
    if (nrhs >= 1) {
	i__ = (nrhs * *nrow - 1) / nperli + 1;
	rhscrd = i__;
	if (*(unsigned char *)guesol == 'G' || *(unsigned char *)guesol == 
		'g') {
	    rhscrd += i__;
	}
	if (*(unsigned char *)&guesol[1] == 'X' || *(unsigned char *)&guesol[
		1] == 'x') {
	    rhscrd += i__;
	}
/* Writing concatenation */
	i__2[0] = 1, a__1[0] = "F";
	i__2[1] = 2, a__1[1] = guesol;
	/* s_cat(rhstyp, a__1, i__2, &c__2, (integer)3); */
    }
L20:

    totcrd = ptrcrd + indcrd + valcrd + rhscrd;
/*     write 4-line or five line header */
    /*
    io___64.ciunit = *iounit;
    s_wsfe(&io___64);
    do_fio(&c__1, title, (integer)72);
    do_fio(&c__1, key, (integer)8);
    do_fio(&c__1, (char *)&totcrd, (integer)sizeof(integer));
    do_fio(&c__1, (char *)&ptrcrd, (integer)sizeof(integer));
    do_fio(&c__1, (char *)&indcrd, (integer)sizeof(integer));
    do_fio(&c__1, (char *)&valcrd, (integer)sizeof(integer));
    do_fio(&c__1, (char *)&rhscrd, (integer)sizeof(integer));
    do_fio(&c__1, type__, (integer)3);
    do_fio(&c__1, (char *)&(*nrow), (integer)sizeof(integer));
    do_fio(&c__1, (char *)&(*ncol), (integer)sizeof(integer));
    do_fio(&c__1, (char *)&nnz, (integer)sizeof(integer));
    do_fio(&c__1, (char *)&nrhs, (integer)sizeof(integer));
    do_fio(&c__1, ptrfmt, (integer)16);
    do_fio(&c__1, indfmt, (integer)16);
    do_fio(&c__1, valfmt, (integer)20);
    do_fio(&c__1, valfmt, (integer)20);
    e_wsfe();
    */
/* ----------------------------------------------------------------------- */
    /*
    nrwindx = 0;
    if (nrhs >= 1) {
	io___66.ciunit = *iounit;
	s_wsfe(&io___66);
	do_fio(&c__1, rhstyp, (integer)3);
	do_fio(&c__1, (char *)&nrhs, (integer)sizeof(integer));
	do_fio(&c__1, (char *)&nrwindx, (integer)sizeof(integer));
	e_wsfe();
    }

    io___67.ciunit = *iounit;
    s_wsfe(&io___67);
    i__1 = *ncol + 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&ia[i__], (integer)sizeof(integer));
    }
    e_wsfe();
    io___68.ciunit = *iounit;
    s_wsfe(&io___68);
    i__1 = nnz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&ja[i__], (integer)sizeof(integer));
    }
    e_wsfe();
    if (*job <= 1) {
	return 0;
    }
    io___69.ciunit = *iounit;
    s_wsfe(&io___69);
    i__1 = nnz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&a[i__], (integer)sizeof(real_number));
    }
    e_wsfe();
    if (*job <= 2) {
	return 0;
    }
    len = *nrow * nrhs;
    next = 1;
    iend = len;
    io___72.ciunit = *iounit;
    s_wsfe(&io___72);
    i__1 = iend;
    for (i__ = next; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&rhs[i__], (integer)sizeof(real_number));
    }
    e_wsfe();
    */

/*     write initial guesses if available */
    /*
    if (*(unsigned char *)guesol == 'G' || *(unsigned char *)guesol == 'g') {
	next += len;
	iend += len;
	io___73.ciunit = *iounit;
	s_wsfe(&io___73);
	i__1 = iend;
	for (i__ = next; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&rhs[i__], (integer)sizeof(real_number));
	}
	e_wsfe();
    }
    */

/*     write exact solutions if available */
    /*
    if (*(unsigned char *)&guesol[1] == 'X' || *(unsigned char *)&guesol[1] ==
	     'x') {
	next += len;
	iend += len;
	io___74.ciunit = *iounit;
	s_wsfe(&io___74);
	i__1 = iend;
	for (i__ = next; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&rhs[i__], (integer)sizeof(real_number));
	}
	e_wsfe();
    }
    */

    return 0;
/* ----------end of prtmt ------------------------------------------------ */
/* ----------------------------------------------------------------------- */
} /* prtmt_ */

/* ----------------------------------------------------------------------c */
/*                          S P A R S K I T                             c */
/* ----------------------------------------------------------------------c */
/*                    FORMAT CONVERSION MODULE                          c */
/* ----------------------------------------------------------------------c */
/* csrdns  : converts a row-stored sparse matrix into the dense format. c */
/* coocsr  : converts coordinate to  to csr format                      c */
/* coicsr  : in-place conversion of coordinate to csr format            c */
/* csrcoo  : converts compressed sparse row to coordinate.              c */
/* csrcsc  : converts compressed sparse row format to compressed sparse c */
/*           column format (transposition)                              c */
/* csrcsc2 : rectangular version of csrcsc                              c */
/* csrdia  : converts a compressed sparse row format into a diagonal    c */
/*           format.                                                    c */
/* csrbnd  : converts a compressed sparse row format into a banded      c */
/*           format (linpack style).                                    c */
/* ----------------------------------------------------------------------c */
/* Subroutine */ int _daskr_csrdns_(integer *nrow, integer *ncol, real_number *a,
	integer *ja, integer *ia, real_number *dns, integer *ndns, integer *
	ierr)
{
    /* System generated locals */
    integer dns_dim1, dns_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, k;

/* ----------------------------------------------------------------------- */
/* Compressed Sparse Row    to    Dense */
/* ----------------------------------------------------------------------- */

/* converts a row-stored sparse matrix into a densely stored one */

/* On entry: */
/* ---------- */

/* nrow	= row-dimension of a */
/* ncol	= column dimension of a */
/* a, */
/* ja, */
/* ia    = input matrix in compressed sparse row format. */
/*         (a=value array, ja=column array, ia=pointer array) */
/* dns   = array where to store dense matrix */
/* ndns	= first dimension of array dns */

/* on return: */
/* ----------- */
/* dns   = the sparse matrix a, ja, ia has been stored in dns(ndns,*) */

/* ierr  = integer error indicator. */
/*         ierr .eq. 0  means normal return */
/*         ierr .eq. i  means that the code has stopped when processing */
/*         row number i, because it found a column number .gt. ncol. */

/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --a;
    --ja;
    --ia;
    dns_dim1 = *ndns;
    dns_offset = 1 + dns_dim1;
    dns -= dns_offset;

    /* Function Body */
    *ierr = 0;
    i__1 = *nrow;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *ncol;
	for (j = 1; j <= i__2; ++j) {
	    dns[i__ + j * dns_dim1] = 0.;
/* L2: */
	}
/* L1: */
    }

    i__1 = *nrow;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ia[i__ + 1] - 1;
	for (k = ia[i__]; k <= i__2; ++k) {
	    j = ja[k];
	    if (j > *ncol) {
		*ierr = i__;
		return 0;
	    }
	    dns[i__ + j * dns_dim1] = a[k];
/* L3: */
	}
/* L4: */
    }
    return 0;
/* ---- end of csrdns ---------------------------------------------------- */
/* ----------------------------------------------------------------------- */
} /* csrdns_ */

/* Subroutine */ int _daskr_coocsr_(integer *nrow, integer *nnz, real_number *a,
	integer *ir, integer *jc, real_number *ao, integer *jao, integer *iao)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k;
    static real_number x;
    static integer k0, iad;

/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*  Coordinate     to   Compressed Sparse Row */
/* ----------------------------------------------------------------------- */
/* converts a matrix that is stored in coordinate format */
/*  a, ir, jc into a row general sparse ao, jao, iao format. */

/* on entry: */
/* --------- */
/* nrow	= dimension of the matrix */
/* nnz	= number of nonzero elements in matrix */
/* a, */
/* ir, */
/* jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz */
/*         nonzero elements of the matrix with a(k) = actual real value of */
/* 	  the elements, ir(k) = its row number and jc(k) = its column */
/* 	  number. The order of the elements is arbitrary. */

/* on return: */
/* ----------- */
/* ir 	is destroyed */

/* ao, jao, iao = matrix in general sparse matrix format with ao */
/* 	continung the real values, jao containing the column indices, */
/* 	and iao being the pointer to the beginning of the row, */
/* 	in arrays ao, jao. */

/* Notes: */
/* ------ This routine is NOT in place.  See coicsr */

/* ------------------------------------------------------------------------ */
    /* Parameter adjustments */
    --iao;
    --jao;
    --ao;
    --jc;
    --ir;
    --a;

    /* Function Body */
    i__1 = *nrow + 1;
    for (k = 1; k <= i__1; ++k) {
	iao[k] = 0;
/* L1: */
    }
/* determine row-lengths. */
    i__1 = *nnz;
    for (k = 1; k <= i__1; ++k) {
	++iao[ir[k]];
/* L2: */
    }
/* starting position of each row.. */
    k = 1;
    i__1 = *nrow + 1;
    for (j = 1; j <= i__1; ++j) {
	k0 = iao[j];
	iao[j] = k;
	k += k0;
/* L3: */
    }
/* go through the structure  once more. Fill in output matrix. */
    i__1 = *nnz;
    for (k = 1; k <= i__1; ++k) {
	i__ = ir[k];
	j = jc[k];
	x = a[k];
	iad = iao[i__];
	ao[iad] = x;
	jao[iad] = j;
	iao[i__] = iad + 1;
/* L4: */
    }
/* shift back iao */
    for (j = *nrow; j >= 1; --j) {
	iao[j + 1] = iao[j];
/* L5: */
    }
    iao[1] = 1;
    return 0;
/* ------------- end of coocsr ------------------------------------------- */
/* ----------------------------------------------------------------------- */
} /* coocsr_ */

/* Subroutine */ int _daskr_coicsr_(integer *n, integer *nnz, integer *job,
	real_number *a, integer *ja, integer *ia, integer *iwk)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k;
    static real_number t;
    static integer init, ipos, inext, jnext;
    static real_number tnext;
    static integer values;

/* ------------------------------------------------------------------------ */
/* IN-PLACE coo-csr conversion routine. */
/* ------------------------------------------------------------------------ */
/* this subroutine converts a matrix stored in coordinate format into */
/* the csr format. The conversion is done in place in that the arrays */
/* a,ja,ia of the result are overwritten onto the original arrays. */
/* ------------------------------------------------------------------------ */
/* on entry: */
/* --------- */
/* n	= integer. row dimension of A. */
/* nnz	= integer. number of nonzero elements in A. */
/* job   = integer. Job indicator. when job=1, the real values in a are */
/*         filled. Otherwise a is not touched and the structure of the */
/*         array only (i.e. ja, ia)  is obtained. */
/* a	= real array of size nnz (number of nonzero elements in A) */
/*         containing the nonzero elements */
/* ja	= integer array of length nnz containing the column positions */
/* 	  of the corresponding elements in a. */
/* ia	= integer array of length nnz containing the row positions */
/* 	  of the corresponding elements in a. */
/* iwk	= integer work array of length n+1 */
/* on return: */
/* ---------- */
/* a */
/* ja */
/* ia	= contains the compressed sparse row data structure for the */
/*         resulting matrix. */
/* Note: */
/* ------- */
/*         the entries of the output matrix are not sorted (the column */
/*         indices in each are not in increasing order) use coocsr */
/*         if you want them sorted. */
/* ----------------------------------------------------------------------c */
/*  Coded by Y. Saad, Sep. 26 1989                                      c */
/* ----------------------------------------------------------------------c */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --iwk;
    --ia;
    --ja;
    --a;

    /* Function Body */
    values = *job == 1;
/* find pointer array for resulting matrix. */
    i__1 = *n + 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iwk[i__] = 0;
/* L35: */
    }
    i__1 = *nnz;
    for (k = 1; k <= i__1; ++k) {
	i__ = ia[k];
	++iwk[i__ + 1];
/* L4: */
    }
/* ------------------------------------------------------------------------ */
    iwk[1] = 1;
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	iwk[i__] = iwk[i__ - 1] + iwk[i__];
/* L44: */
    }

/*     loop for a cycle in chasing process. */

    init = 1;
    k = 0;
L5:
    if (values) {
	t = a[init];
    }
    i__ = ia[init];
    j = ja[init];
    ia[init] = -1;
/* ------------------------------------------------------------------------ */
L6:
    ++k;
/*     current row number is i.  determine  where to go. */
    ipos = iwk[i__];
/*     save the chased element. */
    if (values) {
	tnext = a[ipos];
    }
    inext = ia[ipos];
    jnext = ja[ipos];
/*     then occupy its location. */
    if (values) {
	a[ipos] = t;
    }
    ja[ipos] = j;
/*     update pointer information for next element to come in row i. */
    iwk[i__] = ipos + 1;
/*     determine  next element to be chased, */
    if (ia[ipos] < 0) {
	goto L65;
    }
    t = tnext;
    i__ = inext;
    j = jnext;
    ia[ipos] = -1;
    if (k < *nnz) {
	goto L6;
    }
    goto L70;
L65:
    ++init;
    if (init > *nnz) {
	goto L70;
    }
    if (ia[init] < 0) {
	goto L65;
    }
/*     restart chasing -- */
    goto L5;
L70:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ia[i__ + 1] = iwk[i__];
/* L80: */
    }
    ia[1] = 1;
    return 0;
/* ----------------- end of coicsr ---------------------------------------- */
/* ------------------------------------------------------------------------ */
} /* coicsr_ */

/* Subroutine */ int _daskr_csrcoo_(integer *nrow, integer *job, integer *nzmax,
	real_number *a, integer *ja, integer *ia, integer *nnz, real_number *ao,
	 integer *ir, integer *jc, integer *ierr)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k, k1, k2;

/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/*  Compressed Sparse Row      to      Coordinate */
/* ----------------------------------------------------------------------- */
/* converts a matrix that is stored in coordinate format */
/*  a, ir, jc into a row general sparse ao, jao, iao format. */

/* on entry: */
/* --------- */
/* nrow	= dimension of the matrix. */
/* job   = integer serving as a job indicator. */
/*         if job = 1 fill in only the array ir, ignore jc, and ao. */
/*         if job = 2 fill in ir, and jc but not ao */
/*         if job = 3 fill in everything. */
/*         The reason why these options are provided is that on return */
/*         ao and jc are the same as a, ja. So when job = 3, a and ja are */
/*         simply copied into ao, jc.  When job=2, only jc and ir are */
/*         returned. With job=1 only the array ir is returned. Moreover, */
/*         the algorithm is in place: */
/* 	     call csrcoo (nrow,1,nzmax,a,ja,ia,nnz,a,ia,ja,ierr) */
/*         will write the output matrix in coordinate format on a, ja,ia. */

/* a, */
/* ja, */
/* ia    = matrix in compressed sparse row format. */
/* nzmax = length of space available in ao, ir, jc. */
/*         the code will stop immediatly if the number of */
/*         nonzero elements found in input matrix exceeds nzmax. */

/* on return: */
/* ----------- */
/* ao, ir, jc = matrix in coordinate format. */

/* nnz        = number of nonzero elements in matrix. */
/* ierr       = integer error indicator. */
/*         ierr .eq. 0 means normal retur */
/*         ierr .eq. 1 means that the the code stopped */
/*         because there was no space in ao, ir, jc */
/*         (according to the value of  nzmax). */

/* NOTES: 1)This routine is PARTIALLY in place: csrcoo can be called with */
/*         ao being the same array as as a, and jc the same array as ja. */
/*         but ir CANNOT be the same as ia. */
/*         2) note the order in the output arrays, */
/* ------------------------------------------------------------------------ */
    /* Parameter adjustments */
    --ia;
    --a;
    --ja;
    --ao;
    --ir;
    --jc;

    /* Function Body */
    *ierr = 0;
    *nnz = ia[*nrow + 1] - 1;
    if (*nnz > *nzmax) {
	*ierr = 1;
	return 0;
    }
/* ------------------------------------------------------------------------ */
    switch (*job) {
	case 1:  goto L3;
	case 2:  goto L2;
	case 3:  goto L1;
    }
L1:
    i__1 = *nnz;
    for (k = 1; k <= i__1; ++k) {
	ao[k] = a[k];
/* L10: */
    }
L2:
    i__1 = *nnz;
    for (k = 1; k <= i__1; ++k) {
	jc[k] = ja[k];
/* L11: */
    }

/*     copy backward to allow for in-place processing. */

L3:
    for (i__ = *nrow; i__ >= 1; --i__) {
	k1 = ia[i__ + 1] - 1;
	k2 = ia[i__];
	i__1 = k2;
	for (k = k1; k >= i__1; --k) {
	    ir[k] = i__;
/* L12: */
	}
/* L13: */
    }
    return 0;
/* ------------- end-of-csrcoo ------------------------------------------- */
/* ----------------------------------------------------------------------- */
} /* csrcoo_ */

/* Subroutine */ int _daskr_csrcsc_(integer *n, integer *job, integer *ipos,
	real_number *a, integer *ja, integer *ia, real_number *ao, integer *jao,
	 integer *iao)
{
    extern /* Subroutine */ int _daskr_csrcsc2_(integer *, integer *, integer *,
	    integer *, real_number *, integer *, integer *, real_number *,
	    integer *, integer *);

/* ----------------------------------------------------------------------- */
/* Compressed Sparse Row     to      Compressed Sparse Column */

/* (transposition operation)   Not in place. */
/* ----------------------------------------------------------------------- */
/* -- not in place -- */
/* this subroutine transposes a matrix stored in a, ja, ia format. */
/* --------------- */
/* on entry: */
/* ---------- */
/* n	= dimension of A. */
/* job	= integer to indicate whether to fill the values (job.eq.1) of the */
/*         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1) */

/* ipos  = starting position in ao, jao of the transposed matrix. */
/*         the iao array takes this into account (thus iao(1) is set to ipos.) */
/*         Note: this may be useful if one needs to append the data structure */
/*         of the transpose to that of A. In this case use for example */
/*                call csrcsc (n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2)) */
/* 	  for any other normal usage, enter ipos=1. */
/* a	= real array of length nnz (nnz=number of nonzero elements in input */
/*         matrix) containing the nonzero elements. */
/* ja	= integer array of length nnz containing the column positions */
/* 	  of the corresponding elements in a. */
/* ia	= integer of size n+1. ia(k) contains the position in a, ja of */
/* 	  the beginning of the k-th row. */

/* on return: */
/* ---------- */
/* output arguments: */
/* ao	= real array of size nzz containing the "a" part of the transpose */
/* jao	= integer array of size nnz containing the column indices. */
/* iao	= integer array of size n+1 containing the "ia" index array of */
/* 	  the transpose. */

/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --iao;
    --ia;
    --a;
    --ja;
    --ao;
    --jao;

    /* Function Body */
    _daskr_csrcsc2_(n, n, job, ipos, &a[1], &ja[1], &ia[1], &ao[1], &jao[1], &iao[1])
	    ;
    return 0;
} /* csrcsc_ */

/* Subroutine */ int _daskr_csrcsc2_(integer *n, integer *n2, integer *job, integer *
	ipos, real_number *a, integer *ja, integer *ia, real_number *ao,
	integer *jao, integer *iao)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, next;

/* ----------------------------------------------------------------------- */
/* Compressed Sparse Row     to      Compressed Sparse Column */

/* (transposition operation)   Not in place. */
/* ----------------------------------------------------------------------- */
/* Rectangular version.  n is number of rows of CSR matrix, */
/*                       n2 (input) is number of columns of CSC matrix. */
/* ----------------------------------------------------------------------- */
/* -- not in place -- */
/* this subroutine transposes a matrix stored in a, ja, ia format. */
/* --------------- */
/* on entry: */
/* ---------- */
/* n	= number of rows of CSR matrix. */
/* n2    = number of columns of CSC matrix. */
/* job	= integer to indicate whether to fill the values (job.eq.1) of the */
/*         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1) */

/* ipos  = starting position in ao, jao of the transposed matrix. */
/*         the iao array takes this into account (thus iao(1) is set to ipos.) */
/*         Note: this may be useful if one needs to append the data structure */
/*         of the transpose to that of A. In this case use for example */
/*                call csrcsc2 (n,n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2)) */
/* 	  for any other normal usage, enter ipos=1. */
/* a	= real array of length nnz (nnz=number of nonzero elements in input */
/*         matrix) containing the nonzero elements. */
/* ja	= integer array of length nnz containing the column positions */
/* 	  of the corresponding elements in a. */
/* ia	= integer of size n+1. ia(k) contains the position in a, ja of */
/* 	  the beginning of the k-th row. */

/* on return: */
/* ---------- */
/* output arguments: */
/* ao	= real array of size nzz containing the "a" part of the transpose */
/* jao	= integer array of size nnz containing the column indices. */
/* iao	= integer array of size n+1 containing the "ia" index array of */
/* 	  the transpose. */

/* ----------------------------------------------------------------------- */
/* ----------------- compute lengths of rows of transp(A) ---------------- */
    /* Parameter adjustments */
    --ia;
    --iao;
    --a;
    --ja;
    --ao;
    --jao;

    /* Function Body */
    i__1 = *n2 + 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iao[i__] = 0;
/* L1: */
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ia[i__ + 1] - 1;
	for (k = ia[i__]; k <= i__2; ++k) {
	    j = ja[k] + 1;
	    ++iao[j];
/* L2: */
	}
/* L3: */
    }
/* ---------- compute pointers from lengths ------------------------------ */
    iao[1] = *ipos;
    i__1 = *n2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iao[i__ + 1] = iao[i__] + iao[i__ + 1];
/* L4: */
    }
/* --------------- now do the actual copying ----------------------------- */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ia[i__ + 1] - 1;
	for (k = ia[i__]; k <= i__2; ++k) {
	    j = ja[k];
	    next = iao[j];
	    if (*job == 1) {
		ao[next] = a[k];
	    }
	    jao[next] = i__;
	    iao[j] = next + 1;
/* L62: */
	}
/* L6: */
    }
/* -------------------------- reshift iao and leave ---------------------- */
    for (i__ = *n2; i__ >= 1; --i__) {
	iao[i__ + 1] = iao[i__];
/* L7: */
    }
    iao[1] = *ipos;
/* --------------- end of csrcsc2 ---------------------------------------- */
/* ----------------------------------------------------------------------- */
    return 0;
} /* csrcsc2_ */

/* Subroutine */ int _daskr_csrdia_(integer *n, integer *idiag, integer *job,
	real_number *a, integer *ja, integer *ia, integer *ndiag, real_number *
	diag, integer *ioff, real_number *ao, integer *jao, integer *iao,
	integer *ind)
{
    /* System generated locals */
    integer diag_dim1, diag_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, l, n2, ii, ko, job1, job2, idum, jmax;
    extern /* Subroutine */ int _daskr_infdia_(integer *, integer *, integer *,
	    integer *, integer *);

/* ----------------------------------------------------------------------- */
/* Compressed sparse row     to    diagonal format */
/* ----------------------------------------------------------------------- */
/* this subroutine extracts  idiag diagonals  from the  input matrix a, */
/* a, ia, and puts the rest of  the matrix  in the  output matrix ao, */
/* jao, iao.  The diagonals to be extracted depend  on the  value of job */
/* (see below for details.)  In  the first  case, the  diagonals to be */
/* extracted are simply identified by  their offsets  provided in ioff */
/* by the caller.  In the second case, the  code internally determines */
/* the idiag most significant diagonals, i.e., those  diagonals of the */
/* matrix which  have  the  largest  number  of  nonzero elements, and */
/* extracts them. */
/* ----------------------------------------------------------------------- */
/* on entry: */
/* ---------- */
/* n	= dimension of the matrix a. */
/* idiag = integer equal to the number of diagonals to be extracted. */
/*         Note: on return idiag may be modified. */
/* a, ja, */
/*    ia = matrix stored in a, ja, ia, format */
/* job	= integer. serves as a job indicator.  Job is better thought */
/*         of as a two-digit number job=xy. If the first (x) digit */
/*         is one on entry then the diagonals to be extracted are */
/*         internally determined. In this case csrdia exctracts the */
/*         idiag most important diagonals, i.e. those having the largest */
/*         number on nonzero elements. If the first digit is zero */
/*         then csrdia assumes that ioff(*) contains the offsets */
/*         of the diagonals to be extracted. there is no verification */
/*         that ioff(*) contains valid entries. */
/*         The second (y) digit of job determines whether or not */
/*         the remainder of the matrix is to be written on ao,jao,iao. */
/*         If it is zero  then ao, jao, iao is not filled, i.e., */
/*         the diagonals are found  and put in array diag and the rest is */
/*         is discarded. if it is one, ao, jao, iao contains matrix */
/*         of the remaining elements. */
/*         Thus: */
/*         job= 0 means do not select diagonals internally (pick those */
/*                defined by ioff) and do not fill ao,jao,iao */
/*         job= 1 means do not select diagonals internally */
/*                      and fill ao,jao,iao */
/*         job=10 means  select diagonals internally */
/*                      and do not fill ao,jao,iao */
/*         job=11 means select diagonals internally */
/*                      and fill ao,jao,iao */

/* ndiag = integer equal to the first dimension of array diag. */

/* on return: */
/* ----------- */

/* idiag = number of diagonals found. This may be smaller than its value */
/*         on entry. */
/* diag  = real array of size (ndiag x idiag) containing the diagonals */
/*         of A on return */

/* ioff  = integer array of length idiag, containing the offsets of the */
/*   	  diagonals to be extracted. */
/* ao, jao */
/*  iao  = remainder of the matrix in a, ja, ia format. */
/* work arrays: */
/* ------------ */
/* ind   = integer array of length 2*n-1 used as integer work space. */
/*         needed only when job.ge.10 i.e., in case the diagonals are to */
/*         be selected internally. */

/* Notes: */
/* ------- */
/*    1) The algorithm is in place: ao, jao, iao can be overwritten on */
/*       a, ja, ia if desired */
/*    2) When the code is required to select the diagonals (job .ge. 10) */
/*       the selection of the diagonals is done from left to right */
/*       as a result if several diagonals have the same weight (number */
/*       of nonzero elemnts) the leftmost one is selected first. */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --a;
    --ja;
    --ia;
    diag_dim1 = *ndiag;
    diag_offset = 1 + diag_dim1;
    diag -= diag_offset;
    --ioff;
    --ao;
    --jao;
    --iao;
    --ind;

    /* Function Body */
    job1 = *job / 10;
    job2 = *job - job1 * 10;
    if (job1 == 0) {
	goto L50;
    }
    n2 = *n + *n - 1;
    _daskr_infdia_(n, &ja[1], &ia[1], &ind[1], &idum);
/* ----------- determine diagonals to  accept.---------------------------- */
/* ----------------------------------------------------------------------- */
    ii = 0;
L4:
    ++ii;
    jmax = 0;
    i__1 = n2;
    for (k = 1; k <= i__1; ++k) {
	j = ind[k];
	if (j <= jmax) {
	    goto L41;
	}
	i__ = k;
	jmax = j;
L41:
	;
    }
    if (jmax <= 0) {
	--ii;
	goto L42;
    }
    ioff[ii] = i__ - *n;
    ind[i__] = -jmax;
    if (ii < *idiag) {
	goto L4;
    }
L42:
    *idiag = ii;
/* ---------------- initialize diago to zero ----------------------------- */
L50:
    i__1 = *idiag;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    diag[i__ + j * diag_dim1] = 0.;
/* L54: */
	}
/* L55: */
    }
/* ----------------------------------------------------------------------- */
    ko = 1;
/* ----------------------------------------------------------------------- */
/* extract diagonals and accumulate remaining matrix. */
/* ----------------------------------------------------------------------- */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ia[i__ + 1] - 1;
	for (k = ia[i__]; k <= i__2; ++k) {
	    j = ja[k];
	    i__3 = *idiag;
	    for (l = 1; l <= i__3; ++l) {
		if (j - i__ != ioff[l]) {
		    goto L52;
		}
		diag[i__ + l * diag_dim1] = a[k];
		goto L51;
L52:
		;
	    }
/* --------------- append element not in any diagonal to ao,jao,iao ----- */
	    if (job2 == 0) {
		goto L51;
	    }
	    ao[ko] = a[k];
	    jao[ko] = j;
	    ++ko;
L51:
	    ;
	}
	if (job2 != 0) {
	    ind[i__ + 1] = ko;
	}
/* L6: */
    }
    if (job2 == 0) {
	return 0;
    }
/*     finish with iao */
    iao[1] = 1;
    i__1 = *n + 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	iao[i__] = ind[i__];
/* L7: */
    }
    return 0;
/* ----------- end of csrdia --------------------------------------------- */
/* ----------------------------------------------------------------------- */
} /* csrdia_ */

/* Subroutine */ int _daskr_csrbnd_(integer *n, real_number *a, integer *ja, integer *
	ia, integer *job, real_number *abd, integer *nabd, integer *lowd,
	integer *ml, integer *mu, integer *ierr)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, k, m, ii, mdiag;
    extern /* Subroutine */ int _daskr_getbwd_(integer *, real_number *, integer *,
	    integer *, integer *, integer *);

/* ----------------------------------------------------------------------- */
/*   Compressed Sparse Row  to  Banded (Linpack ) format. */
/* ----------------------------------------------------------------------- */
/* this subroutine converts a general sparse matrix stored in */
/* compressed sparse row format into the banded format. for the */
/* banded format,the Linpack conventions are assumed (see below). */
/* ----------------------------------------------------------------------- */
/* on entry: */
/* ---------- */
/* n	= integer,the actual row dimension of the matrix. */

/* a, */
/* ja, */
/* ia    = input matrix stored in compressed sparse row format. */

/* job	= integer. if job=1 then the values of the lower bandwith ml */
/*         and the upper bandwidth mu are determined internally. */
/*         otherwise it is assumed that the values of ml and mu */
/*         are the correct bandwidths on input. See ml and mu below. */

/* nabd  = integer. first dimension of array abd. */

/* lowd  = integer. this should be set to the row number in abd where */
/*         the lowest diagonal (leftmost) of A is located. */
/*         lowd should be  ( 1  .le.  lowd  .le. nabd). */
/*         if it is not known in advance what lowd should be */
/*         enter lowd = 0 and the default value lowd = ml+mu+1 */
/*         will be chosen. Alternative: call routine getbwd from unary */
/*         first to detrermione ml and mu then define lowd accordingly. */
/*         (Note: the banded solvers in linpack use lowd=2*ml+mu+1. ) */

/* ml	= integer. equal to the bandwidth of the strict lower part of A */
/* mu	= integer. equal to the bandwidth of the strict upper part of A */
/*         thus the total bandwidth of A is ml+mu+1. */
/*         if ml+mu+1 is found to be larger than lowd then an error */
/*         flag is raised (unless lowd = 0). see ierr. */

/* note:   ml and mu are assumed to have	 the correct bandwidth values */
/*         as defined above if job is set to zero on entry. */

/* on return: */
/* ----------- */

/* abd   = real array of dimension abd(nabd,n). */
/*         on return contains the values of the matrix stored in */
/*         banded form. The j-th column of abd contains the elements */
/*         of the j-th column of  the original matrix comprised in the */
/*         band ( i in (j-ml,j+mu) ) with the lowest diagonal at */
/*         the bottom row (row lowd). See details below for this format. */

/* ml	= integer. equal to the bandwidth of the strict lower part of A */
/* mu	= integer. equal to the bandwidth of the strict upper part of A */
/*         if job=1 on entry then these two values are internally computed. */

/* lowd  = integer. row number in abd where the lowest diagonal */
/*         (leftmost) of A is located on return. In case lowd = 0 */
/*         on return, then it is defined to ml+mu+1 on return and the */
/*         lowd will contain this value on return. ` */

/* ierr  = integer. used for error messages. On return: */
/*         ierr .eq. 0  :means normal return */
/*         ierr .eq. -1 : means invalid value for lowd. (either .lt. 0 */
/*         or larger than nabd). */
/*         ierr .eq. -2 : means that lowd is not large enough and as */
/*         result the matrix cannot be stored in array abd. */
/*         lowd should be at least ml+mu+1, where ml and mu are as */
/*         provided on output. */

/* ----------------------------------------------------------------------* */
/* Additional details on banded format.  (this closely follows the      * */
/* format used in linpack. may be useful for converting a matrix into   * */
/* this storage format in order to use the linpack  banded solvers).    * */
/* ----------------------------------------------------------------------* */
/*             ---  band storage format  for matrix abd ---             * */
/* uses ml+mu+1 rows of abd(nabd,*) to store the diagonals of           * */
/* a in rows of abd starting from the lowest (sub)-diagonal  which  is  * */
/* stored in row number lowd of abd. the minimum number of rows needed  * */
/* in abd is ml+mu+1, i.e., the minimum value for lowd is ml+mu+1. the  * */
/* j-th  column  of  abd contains the elements of the j-th column of a, * */
/* from bottom to top: the element a(j+ml,j) is stored in  position     * */
/* abd(lowd,j), then a(j+ml-1,j) in position abd(lowd-1,j) and so on.   * */
/* Generally, the element a(j+k,j) of original matrix a is stored in    * */
/* position abd(lowd+k-ml,j), for k=ml,ml-1,..,0,-1, -mu.               * */
/* The first dimension nabd of abd must be .ge. lowd                    * */
/*                                                                      * */
/*     example [from linpack ]:   if the original matrix is             * */
/*                                                                      * */
/*              11 12 13  0  0  0                                       * */
/*              21 22 23 24  0  0                                       * */
/*               0 32 33 34 35  0     original banded matrix            * */
/*               0  0 43 44 45 46                                       * */
/*               0  0  0 54 55 56                                       * */
/*               0  0  0  0 65 66                                       * */
/*                                                                      * */
/* then  n = 6, ml = 1, mu = 2. lowd should be .ge. 4 (=ml+mu+1)  and   * */
/* if lowd = 5 for example, abd  should be:                             * */
/*                                                                      * */
/* untouched --> x  x  x  x  x  x                                       * */
/*               *  * 13 24 35 46                                       * */
/*               * 12 23 34 45 56    resulting abd matrix in banded     * */
/*              11 22 33 44 55 66    format                             * */
/*  row lowd--> 21 32 43 54 65  *                                       * */
/*                                                                      * */
/* * = not used                                                         * */


/* ----------------------------------------------------------------------* */
/* first determine ml and mu. */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --ia;
    --a;
    --ja;
    abd_dim1 = *nabd;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;

    /* Function Body */
    *ierr = 0;
/* ----------- */
    if (*job == 1) {
    	_daskr_getbwd_(n, &a[1], &ja[1], &ia[1], ml, mu);
    }
    m = *ml + *mu + 1;
    if (*lowd == 0) {
	*lowd = m;
    }
    if (m > *lowd) {
	*ierr = -2;
    }
    if (*lowd > *nabd || *lowd < 0) {
	*ierr = -1;
    }
    if (*ierr < 0) {
	return 0;
    }
/* ------------ */
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = *lowd - i__ + 1;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    abd[ii + j * abd_dim1] = 0.;
/* L10: */
	}
/* L15: */
    }
/* --------------------------------------------------------------------- */
    mdiag = *lowd - *ml;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ia[i__ + 1] - 1;
	for (k = ia[i__]; k <= i__2; ++k) {
	    j = ja[k];
	    abd[i__ - j + mdiag + j * abd_dim1] = a[k];
/* L20: */
	}
/* L30: */
    }
    return 0;
/* ------------- end of csrbnd ------------------------------------------- */
/* ----------------------------------------------------------------------- */
} /* csrbnd_ */

/* ----------------------------------------------------------------------c */
/*                          S P A R S K I T                             c */
/* ----------------------------------------------------------------------c */
/*                     UNARY SUBROUTINES MODULE                         c */
/* ----------------------------------------------------------------------c */
/* rperm  : permutes the rows of a matrix (B = P A)                     c */
/* cperm  : permutes the columns of a matrix (B = A Q)                  c */
/* dperm  : permutes both the rows and columns of a matrix (B = P A Q ) c */
/* dvperm : permutes a real vector (in-place)                           c */
/* ivperm : permutes an integer vector (in-place)                       c */
/* diapos : returns the positions of the diagonal elements in A.        c */
/* getbwd : returns the bandwidth information on a matrix.              c */
/* infdia : obtains information on the diagonals of A.                  c */
/* rnrms  : computes the norms of the rows of A                         c */
/* roscal : scales the rows of a matrix by their norms.                 c */
/* ----------------------------------------------------------------------c */
/* Subroutine */ int _daskr_rperm_(integer *nrow, real_number *a, integer *ja,
	integer *ia, real_number *ao, integer *jao, integer *iao, integer *
	perm, integer *job)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, ii, ko;
    static integer values;

/* ----------------------------------------------------------------------- */
/* this subroutine permutes the rows of a matrix in CSR format. */
/* rperm  computes B = P A  where P is a permutation matrix. */
/* the permutation P is defined through the array perm: for each j, */
/* perm(j) represents the destination row number of row number j. */
/* Youcef Saad -- recoded Jan 28, 1991. */
/* ----------------------------------------------------------------------- */
/* on entry: */
/* ---------- */
/* n 	= dimension of the matrix */
/* a, ja, ia = input matrix in csr format */
/* perm 	= integer array of length nrow containing the permutation arrays */
/* 	  for the rows: perm(i) is the destination of row i in the */
/*         permuted matrix. */
/*         ---> a(i,j) in the original matrix becomes a(perm(i),j) */
/*         in the output  matrix. */

/* job	= integer indicating the work to be done: */
/* 		job = 1	permute a, ja, ia into ao, jao, iao */
/*                       (including the copying of real values ao and */
/*                       the array iao). */
/* 		job .ne. 1 :  ignore real values. */
/*                     (in which case arrays a and ao are not needed nor */
/*                      used). */

/* ------------ */
/* on return: */
/* ------------ */
/* ao, jao, iao = input matrix in a, ja, ia format */
/* note : */
/*        if (job.ne.1)  then the arrays a and ao are not used. */
/* ----------------------------------------------------------------------c */
/*           Y. Saad, May  2, 1990                                      c */
/* ----------------------------------------------------------------------c */
    /* Parameter adjustments */
    --perm;
    --iao;
    --ia;
    --a;
    --ja;
    --ao;
    --jao;

    /* Function Body */
    values = *job == 1;

/*     determine pointers for output matix. */

    i__1 = *nrow;
    for (j = 1; j <= i__1; ++j) {
	i__ = perm[j];
	iao[i__ + 1] = ia[j + 1] - ia[j];
/* L50: */
    }

/* get pointers from lengths */

    iao[1] = 1;
    i__1 = *nrow;
    for (j = 1; j <= i__1; ++j) {
	iao[j + 1] += iao[j];
/* L51: */
    }

/* copying */

    i__1 = *nrow;
    for (ii = 1; ii <= i__1; ++ii) {

/* old row = ii  -- new row = iperm(ii) -- ko = new pointer */

	ko = iao[perm[ii]];
	i__2 = ia[ii + 1] - 1;
	for (k = ia[ii]; k <= i__2; ++k) {
	    jao[ko] = ja[k];
	    if (values) {
		ao[ko] = a[k];
	    }
	    ++ko;
/* L60: */
	}
/* L100: */
    }

    return 0;
/* ---------end-of-rperm ------------------------------------------------- */
/* ----------------------------------------------------------------------- */
} /* rperm_ */

/* Subroutine */ int _daskr_cperm_(integer *nrow, real_number *a, integer *ja,
	integer *ia, real_number *ao, integer *jao, integer *iao, integer *
	perm, integer *job)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k, nnz;

/* ----------------------------------------------------------------------- */
/* this subroutine permutes the columns of a matrix a, ja, ia. */
/* the result is written in the output matrix  ao, jao, iao. */
/* cperm computes B = A P, where  P is a permutation matrix */
/* that maps column j into column perm(j), i.e., on return */
/*      a(i,j) becomes a(i,perm(j)) in new matrix */
/* Y. Saad, May 2, 1990 / modified Jan. 28, 1991. */
/* ----------------------------------------------------------------------- */
/* on entry: */
/* ---------- */
/* nrow 	= row dimension of the matrix */

/* a, ja, ia = input matrix in csr format. */

/* perm	= integer array of length ncol (number of columns of A */
/*         containing the permutation array  the columns: */
/*         a(i,j) in the original matrix becomes a(i,perm(j)) */
/*         in the output matrix. */

/* job	= integer indicating the work to be done: */
/* 		job = 1	permute a, ja, ia into ao, jao, iao */
/*                       (including the copying of real values ao and */
/*                       the array iao). */
/* 		job .ne. 1 :  ignore real values ao and ignore iao. */

/* ------------ */
/* on return: */
/* ------------ */
/* ao, jao, iao = input matrix in a, ja, ia format (array ao not needed) */

/* Notes: */
/* ------- */
/* 1. if job=1 then ao, iao are not used. */
/* 2. This routine is in place: ja, jao can be the same. */
/* 3. If the matrix is initially sorted (by increasing column number) */
/*    then ao,jao,iao  may not be on return. */

/* ----------------------------------------------------------------------c */
/* local parameters: */

    /* Parameter adjustments */
    --iao;
    --ia;
    --a;
    --ja;
    --ao;
    --jao;
    --perm;

    /* Function Body */
    nnz = ia[*nrow + 1] - 1;
    i__1 = nnz;
    for (k = 1; k <= i__1; ++k) {
	jao[k] = perm[ja[k]];
/* L100: */
    }

/*     done with ja array. return if no need to touch values. */

    if (*job != 1) {
	return 0;
    }

/* else get new pointers -- and copy values too. */

    i__1 = *nrow + 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iao[i__] = ia[i__];
/* L1: */
    }

    i__1 = nnz;
    for (k = 1; k <= i__1; ++k) {
	ao[k] = a[k];
/* L2: */
    }

    return 0;
/* ---------end-of-cperm-------------------------------------------------- */
/* ----------------------------------------------------------------------- */
} /* cperm_ */

/* Subroutine */ int _daskr_dperm_(integer *nrow, real_number *a, integer *ja,
	integer *ia, real_number *ao, integer *jao, integer *iao, integer *
	perm, integer *qperm, integer *job)
{
    extern /* Subroutine */ int _daskr_cperm_(integer *, real_number *, integer *,
	    integer *, real_number *, integer *, integer *, integer *, integer
	    *), _daskr_rperm_(integer *, real_number *, integer *, integer *,
	    real_number *, integer *, integer *, integer *, integer *);
    static integer locjob;

/* ----------------------------------------------------------------------- */
/* This routine permutes the rows and columns of a matrix stored in CSR */
/* format. i.e., it computes P A Q, where P, Q are permutation matrices. */
/* P maps row i into row perm(i) and Q maps column j into column qperm(j): */
/*      a(i,j)    becomes   a(perm(i),qperm(j)) in new matrix */
/* In the particular case where Q is the transpose of P (symmetric */
/* permutation of A) then qperm is not needed. */
/* note that qperm should be of length ncol (number of columns) but this */
/* is not checked. */
/* ----------------------------------------------------------------------- */
/* Y. Saad, Sep. 21 1989 / recoded Jan. 28 1991. */
/* ----------------------------------------------------------------------- */
/* on entry: */
/* ---------- */
/* n 	= dimension of the matrix */
/* a, ja, */
/*    ia = input matrix in a, ja, ia format */
/* perm 	= integer array of length n containing the permutation arrays */
/* 	  for the rows: perm(i) is the destination of row i in the */
/*         permuted matrix -- also the destination of column i in case */
/*         permutation is symmetric (job .le. 2) */

/* qperm	= same thing for the columns. This should be provided only */
/*         if job=3 or job=4, i.e., only in the case of a nonsymmetric */
/* 	  permutation of rows and columns. Otherwise qperm is a dummy */

/* job	= integer indicating the work to be done: */
/* * job = 1,2 permutation is symmetric  Ao :== P * A * transp(P) */
/* 		job = 1	permute a, ja, ia into ao, jao, iao */
/* 		job = 2 permute matrix ignoring real values. */
/* * job = 3,4 permutation is non-symmetric  Ao :== P * A * Q */
/* 		job = 3	permute a, ja, ia into ao, jao, iao */
/* 		job = 4 permute matrix ignoring real values. */

/* on return: */
/* ----------- */
/* ao, jao, iao = input matrix in a, ja, ia format */

/* in case job .eq. 2 or job .eq. 4, a and ao are never referred to */
/* and can be dummy arguments. */
/* Notes: */
/* ------- */
/*  1) algorithm is in place */
/*  2) column indices may not be sorted on return even  though they may be */
/*     on entry. */
/* ----------------------------------------------------------------------c */
/* local variables */

/*     locjob indicates whether or not real values must be copied. */

    /* Parameter adjustments */
    --perm;
    --iao;
    --ia;
    --a;
    --ja;
    --ao;
    --jao;
    --qperm;

    /* Function Body */
    locjob = *job % 2;

/* permute rows first */

    _daskr_rperm_(nrow, &a[1], &ja[1], &ia[1], &ao[1], &jao[1], &iao[1], &perm[1], &
	    locjob);

/* then permute columns */

    locjob = 0;

    if (*job <= 2) {
    	_daskr_cperm_(nrow, &ao[1], &jao[1], &iao[1], &ao[1], &jao[1], &iao[1], &
		perm[1], &locjob);
    } else {
    	_daskr_cperm_(nrow, &ao[1], &jao[1], &iao[1], &ao[1], &jao[1], &iao[1], &
		qperm[1], &locjob);
    }

    return 0;
/* -------end-of-dperm---------------------------------------------------- */
/* ----------------------------------------------------------------------- */
} /* dperm_ */

/* Subroutine */ int _daskr_dvperm_(integer *n, real_number *x, integer *perm)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, k, ii;
    static real_number tmp, tmp1;
    static integer init, next;

/* ----------------------------------------------------------------------- */
/* this subroutine performs an in-place permutation of a real vector x */
/* according to the permutation array perm(*), i.e., on return, */
/* the vector x satisfies, */

/* 	x(perm(j)) :== x(j), j=1,2,.., n */

/* ----------------------------------------------------------------------- */
/* on entry: */
/* --------- */
/* n 	= length of vector x. */
/* perm 	= integer array of length n containing the permutation  array. */
/* x	= input vector */

/* on return: */
/* ---------- */
/* x	= vector x permuted according to x(perm(*)) :=  x(*) */

/* ----------------------------------------------------------------------c */
/*           Y. Saad, Sep. 21 1989                                      c */
/* ----------------------------------------------------------------------c */
/* local variables */

    /* Parameter adjustments */
    --perm;
    --x;

    /* Function Body */
    init = 1;
    tmp = x[init];
    ii = perm[init];
    perm[init] = -perm[init];
    k = 0;

/* loop */

L6:
    ++k;

/* save the chased element -- */

    tmp1 = x[ii];
    x[ii] = tmp;
    next = perm[ii];
    if (next < 0) {
	goto L65;
    }

/* test for end */

    if (k > *n) {
	goto L101;
    }
    tmp = tmp1;
    perm[ii] = -perm[ii];
    ii = next;

/* end loop */

    goto L6;

/* reinitilaize cycle -- */

L65:
    ++init;
    if (init > *n) {
	goto L101;
    }
    if (perm[init] < 0) {
	goto L65;
    }
    tmp = x[init];
    ii = perm[init];
    perm[init] = -perm[init];
    goto L6;

L101:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	perm[j] = -perm[j];
/* L200: */
    }

    return 0;
/* -------------------end-of-dvperm--------------------------------------- */
/* ----------------------------------------------------------------------- */
} /* dvperm_ */

/* Subroutine */ int _daskr_ivperm_(integer *n, integer *ix, integer *perm)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, k, ii, tmp, tmp1, init, next;

/* ----------------------------------------------------------------------- */
/* this subroutine performs an in-place permutation of an integer vector */
/* ix according to the permutation array perm(*), i.e., on return, */
/* the vector x satisfies, */

/* 	ix(perm(j)) :== ix(j), j=1,2,.., n */

/* ----------------------------------------------------------------------- */
/* on entry: */
/* --------- */
/* n 	= length of vector x. */
/* perm 	= integer array of length n containing the permutation  array. */
/* ix	= input vector */

/* on return: */
/* ---------- */
/* ix	= vector x permuted according to ix(perm(*)) :=  ix(*) */

/* ----------------------------------------------------------------------c */
/*           Y. Saad, Sep. 21 1989                                      c */
/* ----------------------------------------------------------------------c */
/* local variables */

    /* Parameter adjustments */
    --perm;
    --ix;

    /* Function Body */
    init = 1;
    tmp = ix[init];
    ii = perm[init];
    perm[init] = -perm[init];
    k = 0;

/* loop */

L6:
    ++k;

/* save the chased element -- */

    tmp1 = ix[ii];
    ix[ii] = tmp;
    next = perm[ii];
    if (next < 0) {
	goto L65;
    }

/* test for end */

    if (k > *n) {
	goto L101;
    }
    tmp = tmp1;
    perm[ii] = -perm[ii];
    ii = next;

/* end loop */

    goto L6;

/* reinitilaize cycle -- */

L65:
    ++init;
    if (init > *n) {
	goto L101;
    }
    if (perm[init] < 0) {
	goto L65;
    }
    tmp = ix[init];
    ii = perm[init];
    perm[init] = -perm[init];
    goto L6;

L101:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	perm[j] = -perm[j];
/* L200: */
    }

    return 0;
/* -------------------end-of-ivperm--------------------------------------- */
/* ----------------------------------------------------------------------- */
} /* ivperm_ */

/* Subroutine */ int _daskr_diapos_(integer *n, integer *ja, integer *ia, integer *
	idiag)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, k;

/* ----------------------------------------------------------------------- */
/* this subroutine returns the positions of the diagonal elements of a */
/* sparse matrix a, ja, ia, in the array idiag. */
/* ----------------------------------------------------------------------- */
/* on entry: */
/* ---------- */

/* n	= integer. row dimension of the matrix a. */
/* a,ja, */
/*    ia = matrix stored compressed sparse row format. a array skipped. */

/* on return: */
/* ----------- */
/* idiag  = integer array of length n. The i-th entry of idiag */
/*          points to the diagonal element a(i,i) in the arrays */
/*          a, ja. (i.e., a(idiag(i)) = element A(i,i) of matrix A) */
/*          if no diagonal element is found the entry is set to 0. */
/* ----------------------------------------------------------------------c */
/*           Y. Saad, March, 1990 */
/* ----------------------------------------------------------------------c */
    /* Parameter adjustments */
    --idiag;
    --ia;
    --ja;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	idiag[i__] = 0;
/* L1: */
    }

/*     sweep through data structure. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ia[i__ + 1] - 1;
	for (k = ia[i__]; k <= i__2; ++k) {
	    if (ja[k] == i__) {
		idiag[i__] = k;
	    }
/* L51: */
	}
/* L6: */
    }
/* ----------- -end-of-diapos--------------------------------------------- */
/* ----------------------------------------------------------------------- */
    return 0;
} /* diapos_ */

/* Subroutine */ int _daskr_getbwd_(integer *n, real_number *a, integer *ja, integer *
	ia, integer *ml, integer *mu)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, k, ldist;

/* ----------------------------------------------------------------------- */
/* gets the bandwidth of lower part and upper part of A. */
/* does not assume that A is sorted. */
/* ----------------------------------------------------------------------- */
/* on entry: */
/* ---------- */
/* n	= integer = the row dimension of the matrix */
/* a, ja, */
/*    ia = matrix in compressed sparse row format. */

/* on return: */
/* ----------- */
/* ml	= integer. The bandwidth of the strict lower part of A */
/* mu	= integer. The bandwidth of the strict upper part of A */

/* Notes: */
/* ===== ml and mu are allowed to be negative or return. This may be */
/*       useful since it will tell us whether a band is confined */
/*       in the strict  upper/lower triangular part. */
/*       indeed the definitions of ml and mu are */

/*       ml = max ( (i-j)  s.t. a(i,j) .ne. 0  ) */
/*       mu = max ( (j-i)  s.t. a(i,j) .ne. 0  ) */
/* ----------------------------------------------------------------------c */
/* Y. Saad, Sep. 21 1989                                                c */
/* ----------------------------------------------------------------------c */
    /* Parameter adjustments */
    --ia;
    --a;
    --ja;

    /* Function Body */
    *ml = -(*n);
    *mu = -(*n);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ia[i__ + 1] - 1;
	for (k = ia[i__]; k <= i__2; ++k) {
	    ldist = i__ - ja[k];
	    *ml = MAX(*ml,ldist);
/* Computing MAX */
	    i__3 = *mu, i__4 = -ldist;
	    *mu = MAX(i__3,i__4);
/* L31: */
	}
/* L3: */
    }
    return 0;
/* ---------------end-of-getbwd ------------------------------------------ */
/* ----------------------------------------------------------------------- */
} /* getbwd_ */

/* Subroutine */ int _daskr_infdia_(integer *n, integer *ja, integer *ia, integer *
	ind, integer *idiag)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, n2;

/* ----------------------------------------------------------------------- */
/*     obtains information on the diagonals of A. */
/* ----------------------------------------------------------------------- */
/* this subroutine finds the lengths of each of the 2*n-1 diagonals of A */
/* it also outputs the number of nonzero diagonals found. */
/* ----------------------------------------------------------------------- */
/* on entry: */
/* ---------- */
/* n	= dimension of the matrix a. */

/* a,    ..... not needed here. */
/* ja, */
/* ia    = matrix stored in csr format */

/* on return: */
/* ----------- */

/* idiag = integer. number of nonzero diagonals found. */

/* ind   = integer array of length at least 2*n-1. The k-th entry in */
/*         ind contains the number of nonzero elements in the diagonal */
/*         number k, the numbering beeing from the lowermost diagonal */
/*         (bottom-left). In other words ind(k) = length of diagonal */
/*         whose offset wrt the main diagonal is = - n + k. */
/* ----------------------------------------------------------------------c */
/*           Y. Saad, Sep. 21 1989                                      c */
/* ----------------------------------------------------------------------c */
    /* Parameter adjustments */
    --ind;
    --ia;
    --ja;

    /* Function Body */
    n2 = *n + *n - 1;
    i__1 = n2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ind[i__] = 0;
/* L1: */
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ia[i__ + 1] - 1;
	for (k = ia[i__]; k <= i__2; ++k) {
	    j = ja[k];
	    ++ind[*n + j - i__];
/* L2: */
	}
/* L3: */
    }
/*     count the nonzero ones. */
    *idiag = 0;
    i__1 = n2;
    for (k = 1; k <= i__1; ++k) {
	if (ind[k] != 0) {
	    ++(*idiag);
	}
/* L41: */
    }
    return 0;
/* done */
/* ------end-of-infdia --------------------------------------------------- */
/* ----------------------------------------------------------------------- */
} /* infdia_ */

/* Subroutine */ int _daskr_rnrms_(integer *nrow, integer *nrm, real_number *a,
	integer *ja, integer *ia, real_number *diag)
{
    /* System generated locals */
    integer i__1, i__2;
    real_number d__1, d__2, d__3;

    /* Local variables */
    static integer k, k1, k2, ii;
    static real_number scal;

/* ----------------------------------------------------------------------- */
/* gets the norms of each row of A. (choice of three norms) */
/* ----------------------------------------------------------------------- */
/* on entry: */
/* --------- */
/* nrow	= integer. The row dimension of A */

/* nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2 */
/*                  means the 2-nrm, nrm = 0 means max norm */

/* a, */
/* ja, */
/* ia   = Matrix A in compressed sparse row format. */

/* on return: */
/* ---------- */

/* diag = real vector of length nrow containing the norms */

/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --diag;
    --ia;
    --a;
    --ja;

    /* Function Body */
    i__1 = *nrow;
    for (ii = 1; ii <= i__1; ++ii) {

/*     compute the norm if each element. */

	scal = 0.;
	k1 = ia[ii];
	k2 = ia[ii + 1] - 1;
	if (*nrm == 0) {
	    i__2 = k2;
	    for (k = k1; k <= i__2; ++k) {
/* Computing MAX */
		d__2 = scal, d__3 = (d__1 = a[k], fabs(d__1));
		scal = MAX(d__2,d__3);
/* L2: */
	    }
	} else if (*nrm == 1) {
	    i__2 = k2;
	    for (k = k1; k <= i__2; ++k) {
		scal += (d__1 = a[k], fabs(d__1));
/* L3: */
	    }
	} else {
	    i__2 = k2;
	    for (k = k1; k <= i__2; ++k) {
/* Computing 2nd power */
		d__1 = a[k];
		scal += d__1 * d__1;
/* L4: */
	    }
	}
	if (*nrm == 2) {
	    scal = sqrt(scal);
	}
	diag[ii] = scal;
/* L1: */
    }
    return 0;
/* ----------------------------------------------------------------------- */
/* -------------end-of-rnrms---------------------------------------------- */
} /* rnrms_ */

/* Subroutine */ int _daskr_roscal_(integer *nrow, integer *job, integer *nrm,
	real_number *a, integer *ja, integer *ia, real_number *diag, real_number
	*b, integer *jb, integer *ib, integer *ierr)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;
    extern /* Subroutine */ int _daskr_rnrms_(integer *, integer *, real_number *,
	    integer *, integer *, real_number *), _daskr_diamua_(integer *, integer *,
	     real_number *, integer *, integer *, real_number *, real_number *,
	    integer *, integer *);

/* ----------------------------------------------------------------------- */
/* scales the rows of A such that their norms are one on return */
/* 3 choices of norms: 1-norm, 2-norm, max-norm. */
/* ----------------------------------------------------------------------- */
/* on entry: */
/* --------- */
/* nrow	= integer. The row dimension of A */

/* job   = integer. job indicator. Job=0 means get array b only */
/*         job = 1 means get b, and the integer arrays ib, jb. */

/* nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2 */
/*                  means the 2-nrm, nrm = 0 means max norm */

/* a, */
/* ja, */
/* ia   = Matrix A in compressed sparse row format. */

/* on return: */
/* ---------- */

/* diag = diagonal matrix stored as a vector containing the matrix */
/*        by which the rows have been scaled, i.e., on return */
/*        we have B = Diag*A. */

/* b, */
/* jb, */
/* ib	= resulting matrix B in compressed sparse row sparse format. */

/* ierr  = error message. ierr=0     : Normal return */
/*                        ierr=i > 0 : Row number i is a zero row. */
/* Notes: */
/* ------- */
/* 1)        The column dimension of A is not needed. */
/* 2)        algorithm in place (B can take the place of A). */
/* ----------------------------------------------------------------- */
    /* Parameter adjustments */
    --ib;
    --diag;
    --ia;
    --a;
    --ja;
    --b;
    --jb;

    /* Function Body */
    _daskr_rnrms_(nrow, nrm, &a[1], &ja[1], &ia[1], &diag[1]);
    *ierr = 0;
    i__1 = *nrow;
    for (j = 1; j <= i__1; ++j) {
	if (diag[j] == 0.) {
	    *ierr = j;
	    return 0;
	} else {
	    diag[j] = 1. / diag[j];
	}
/* L1: */
    }
    _daskr_diamua_(nrow, job, &a[1], &ja[1], &ia[1], &diag[1], &b[1], &jb[1], &ib[1])
	    ;
    return 0;
/* -------end-of-roscal--------------------------------------------------- */
/* ----------------------------------------------------------------------- */
} /* roscal_ */

/* ----------------------------------------------------------------------c */
/*                          S P A R S K I T                             c */
/* ----------------------------------------------------------------------c */
/*                   ITERATIVE SOLVERS MODULE                           c */
/* ----------------------------------------------------------------------c */
/* ILUT    : Incomplete LU factorization with dual truncation strategy  c */
/* ILUTP   : ILUT with column  pivoting                                 c */
/* LUSOL   : forward followed by backward triangular solve (Precond.)   c */
/* QSPLIT  : quick split routine used by ilut to sort out the k largest c */
/*           elements in absolute value                                 c */
/* ----------------------------------------------------------------------c */
/* Subroutine */ int _daskr_ilut_(integer *n, real_number *a, integer *ja, integer *
	ia, integer *lfil, real_number *droptol, real_number *alu, integer *jlu,
	 integer *ju, integer *iwk, real_number *w, integer *jw, integer *ierr)
{
    /* System generated locals */
    integer i__1, i__2;
    real_number d__1;

    /* Local variables */
    static integer i__, j, k;
    static real_number s, t;
    static integer j1, j2, ii, jj, ju0, len;
    static real_number fact;
    static integer lenl, lenu, jpos, jrow;
    static real_number tnorm;
    extern /* Subroutine */ int _daskr_qsplit_(real_number *, integer *, integer *,
	    integer *);

/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------* */
/*                      *** ILUT preconditioner ***                     * */
/*      incomplete LU factorization with dual truncation mechanism      * */
/* ----------------------------------------------------------------------* */
/*     Author: Yousef Saad *May, 5, 1990, Latest revision, August 1996  * */
/* ----------------------------------------------------------------------* */
/* PARAMETERS */
/* ----------- */

/* on entry: */
/* ========== */
/* n       = integer. The row dimension of the matrix A. The matrix */

/* a,ja,ia = matrix stored in Compressed Sparse Row format. */

/* lfil    = integer. The fill-in parameter. Each row of L and each row */
/*           of U will have a maximum of lfil elements (excluding the */
/*           diagonal element). lfil must be .ge. 0. */
/*           ** WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO */
/*           EARLIER VERSIONS. */

/* droptol = real*8. Sets the threshold for dropping small terms in the */
/*           factorization. See below for details on dropping strategy. */


/* iwk     = integer. The lengths of arrays alu and jlu. If the arrays */
/*           are not big enough to store the ILU factorizations, ilut */
/*           will stop with an error message. */

/* On return: */
/* =========== */

/* alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing */
/*           the L and U factors together. The diagonal (stored in */
/*           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix */
/*           contains the i-th row of L (excluding the diagonal entry=1) */
/*           followed by the i-th row of U. */

/* ju      = integer array of length n containing the pointers to */
/*           the beginning of each row of U in the matrix alu,jlu. */

/* ierr    = integer. Error message with the following meaning. */
/*           ierr  = 0    --> successful return. */
/*           ierr .gt. 0  --> zero pivot encountered at step number ierr. */
/*           ierr  = -1   --> Error. input matrix may be wrong. */
/*                            (The elimination process has generated a */
/*                            row in L or U whose length is .gt.  n.) */
/*           ierr  = -2   --> The matrix L overflows the array al. */
/*           ierr  = -3   --> The matrix U overflows the array alu. */
/*           ierr  = -4   --> Illegal value for lfil. */
/*           ierr  = -5   --> zero row encountered. */

/* work arrays: */
/* ============= */
/* jw      = integer work array of length 2*n. */
/* w       = real work array of length n */

/* ---------------------------------------------------------------------- */
/* w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u] */
/* jw(n+1:2n)  stores nonzero indicators */

/* Notes: */
/* ------ */
/* The diagonal elements of the input matrix must be  nonzero (at least */
/* 'structurally'). */

/* ----------------------------------------------------------------------* */
/* ---- Dual drop strategy works as follows.                             * */
/*                                                                      * */
/*     1) Theresholding in L and U as set by droptol. Any element whose * */
/*        magnitude is less than some tolerance (relative to the abs    * */
/*        value of diagonal element in u) is dropped.                   * */
/*                                                                      * */
/*     2) Keeping only the largest lfil elements in the i-th row of L   * */
/*        and the largest lfil elements in the i-th row of U (excluding * */
/*        diagonal elements).                                           * */
/*                                                                      * */
/* Flexibility: one  can use  droptol=0  to get  a strategy  based on   * */
/* keeping  the largest  elements in  each row  of L  and U.   Taking   * */
/* droptol .ne.  0 but lfil=n will give  the usual threshold strategy   * */
/* (however, fill-in is then mpredictible).                             * */
/* ----------------------------------------------------------------------* */
/*     locals */
    /* Parameter adjustments */
    --jw;
    --w;
    --ju;
    --ia;
    --a;
    --ja;
    --alu;
    --jlu;

    /* Function Body */
    if (*lfil < 0) {
	goto L998;
    }
/* ----------------------------------------------------------------------- */
/*     initialize ju0 (points to next element to be added to alu,jlu) */
/*     and pointer array. */
/* ----------------------------------------------------------------------- */
    ju0 = *n + 2;
    jlu[1] = ju0;

/*     initialize nonzero indicator array. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	jw[*n + j] = 0;
/* L1: */
    }
/* ----------------------------------------------------------------------- */
/*     beginning of main loop. */
/* ----------------------------------------------------------------------- */
    i__1 = *n;
    for (ii = 1; ii <= i__1; ++ii) {
	j1 = ia[ii];
	j2 = ia[ii + 1] - 1;
	tnorm = 0.;
	i__2 = j2;
	for (k = j1; k <= i__2; ++k) {
	    tnorm += (d__1 = a[k], fabs(d__1));
/* L501: */
	}
	if (tnorm == 0.) {
	    goto L999;
	}
	tnorm /= (real_number) (j2 - j1 + 1);

/*     unpack L-part and U-part of row of A in arrays w */

	lenu = 1;
	lenl = 0;
	jw[ii] = ii;
	w[ii] = 0.f;
	jw[*n + ii] = ii;

	i__2 = j2;
	for (j = j1; j <= i__2; ++j) {
	    k = ja[j];
	    t = a[j];
	    if (k < ii) {
		++lenl;
		jw[lenl] = k;
		w[lenl] = t;
		jw[*n + k] = lenl;
	    } else if (k == ii) {
		w[ii] = t;
	    } else {
		++lenu;
		jpos = ii + lenu - 1;
		jw[jpos] = k;
		w[jpos] = t;
		jw[*n + k] = jpos;
	    }
/* L170: */
	}
	jj = 0;
	len = 0;

/*     eliminate previous rows */

L150:
	++jj;
	if (jj > lenl) {
	    goto L160;
	}
/* ----------------------------------------------------------------------- */
/*     in order to do the elimination in the correct order we must select */
/*     the smallest column index among jw(k), k=jj+1, ..., lenl. */
/* ----------------------------------------------------------------------- */
	jrow = jw[jj];
	k = jj;

/*     determine smallest column index */

	i__2 = lenl;
	for (j = jj + 1; j <= i__2; ++j) {
	    if (jw[j] < jrow) {
		jrow = jw[j];
		k = j;
	    }
/* L151: */
	}

	if (k != jj) {
/*     exchange in jw */
	    j = jw[jj];
	    jw[jj] = jw[k];
	    jw[k] = j;
/*     exchange in jr */
	    jw[*n + jrow] = jj;
	    jw[*n + j] = k;
/*     exchange in w */
	    s = w[jj];
	    w[jj] = w[k];
	    w[k] = s;
	}

/*     zero out element in row by setting jw(n+jrow) to zero. */

	jw[*n + jrow] = 0;

/*     get the multiplier for row to be eliminated (jrow). */

	fact = w[jj] * alu[jrow];
	if (fabs(fact) <= *droptol) {
	    goto L150;
	}

/*     combine current row and row jrow */

	i__2 = jlu[jrow + 1] - 1;
	for (k = ju[jrow]; k <= i__2; ++k) {
	    s = fact * alu[k];
	    j = jlu[k];
	    jpos = jw[*n + j];
	    if (j >= ii) {

/*     dealing with upper part. */

		if (jpos == 0) {

/*     this is a fill-in element */

		    ++lenu;
		    if (lenu > *n) {
			goto L995;
		    }
		    i__ = ii + lenu - 1;
		    jw[i__] = j;
		    jw[*n + j] = i__;
		    w[i__] = -s;
		} else {

/*     this is not a fill-in element */

		    w[jpos] -= s;
		}
	    } else {

/*     dealing  with lower part. */

		if (jpos == 0) {

/*     this is a fill-in element */

		    ++lenl;
		    if (lenl > *n) {
			goto L995;
		    }
		    jw[lenl] = j;
		    jw[*n + j] = lenl;
		    w[lenl] = -s;
		} else {

/*     this is not a fill-in element */

		    w[jpos] -= s;
		}
	    }
/* L203: */
	}

/*     store this pivot element -- (from left to right -- no danger of */
/*     overlap with the working elements in L (pivots). */

	++len;
	w[len] = fact;
	jw[len] = jrow;
	goto L150;
L160:

/*     reset double-pointer to zero (U-part) */

	i__2 = lenu;
	for (k = 1; k <= i__2; ++k) {
	    jw[*n + jw[ii + k - 1]] = 0;
/* L308: */
	}

/*     update L-matrix */

	lenl = len;
	len = MIN(lenl,*lfil);

/*     sort by quick-split */

	_daskr_qsplit_(&w[1], &jw[1], &lenl, &len);

/*     store L-part */

	i__2 = len;
	for (k = 1; k <= i__2; ++k) {
	    if (ju0 > *iwk) {
		goto L996;
	    }
	    alu[ju0] = w[k];
	    jlu[ju0] = jw[k];
	    ++ju0;
/* L204: */
	}

/*     save pointer to beginning of row ii of U */

	ju[ii] = ju0;

/*     update U-matrix -- first apply dropping strategy */

	len = 0;
	i__2 = lenu - 1;
	for (k = 1; k <= i__2; ++k) {
	    if ((d__1 = w[ii + k], fabs(d__1)) > *droptol * tnorm) {
		++len;
		w[ii + len] = w[ii + k];
		jw[ii + len] = jw[ii + k];
	    }
	}
	lenu = len + 1;
	len = MIN(lenu,*lfil);

	if (lenu > 1) {
	    i__2 = lenu - 1;
	    _daskr_qsplit_(&w[ii + 1], &jw[ii + 1], &i__2, &len);
	}

/*     copy */

	t = (d__1 = w[ii], fabs(d__1));
	if (len + ju0 > *iwk) {
	    goto L997;
	}
	i__2 = ii + len - 1;
	for (k = ii + 1; k <= i__2; ++k) {
	    jlu[ju0] = jw[k];
	    alu[ju0] = w[k];
	    t += (d__1 = w[k], fabs(d__1));
	    ++ju0;
/* L302: */
	}

/*     store inverse of diagonal element of u */

	if (w[ii] == 0.) {
	    w[ii] = (*droptol + 1e-4f) * tnorm;
	}

	alu[ii] = 1. / w[ii];

/*     update pointer to beginning of next row of U. */

	jlu[ii + 1] = ju0;
/* ----------------------------------------------------------------------- */
/*     end main loop */
/* ----------------------------------------------------------------------- */
/* L500: */
    }
    *ierr = 0;
    return 0;

/*     incomprehensible error. Matrix must be wrong. */

L995:
    *ierr = -1;
    return 0;

/*     insufficient storage in L. */

L996:
    *ierr = -2;
    return 0;

/*     insufficient storage in U. */

L997:
    *ierr = -3;
    return 0;

/*     illegal lfil entered. */

L998:
    *ierr = -4;
    return 0;

/*     zero row encountered */

L999:
    *ierr = -5;
    return 0;
/* ----------------end-of-ilut-------------------------------------------- */
/* ----------------------------------------------------------------------- */
} /* ilut_ */

/* Subroutine */ int _daskr_ilutp_(integer *n, real_number *a, integer *ja, integer *
	ia, integer *lfil, real_number *droptol, real_number *permtol, integer *
	mbloc, real_number *alu, integer *jlu, integer *ju, integer *iwk,
	real_number *w, integer *jw, integer *iperm, integer *ierr)
{
    /* System generated locals */
    integer i__1, i__2;
    real_number d__1;

    /* Local variables */
    static integer i__, j, k;
    static real_number s, t;
    static integer j1, j2, ii, jj, ju0, len;
    static real_number tmp, fact;
    static integer lenl, imax, lenu, icut, jpos;
    static real_number xmax;
    static integer jrow;
    static real_number xmax0, tnorm;
    extern /* Subroutine */ int _daskr_qsplit_(real_number *, integer *, integer *,
	    integer *);

/* ----------------------------------------------------------------------- */
/*     implicit none */
/* ----------------------------------------------------------------------* */
/*       *** ILUTP preconditioner -- ILUT with pivoting  ***            * */
/*      incomplete LU factorization with dual truncation mechanism      * */
/* ----------------------------------------------------------------------* */
/* author Yousef Saad *Sep 8, 1993 -- Latest revision, August 1996.     * */
/* ----------------------------------------------------------------------* */
/* on entry: */
/* ========== */
/* n       = integer. The dimension of the matrix A. */

/* a,ja,ia = matrix stored in Compressed Sparse Row format. */
/*           ON RETURN THE COLUMNS OF A ARE PERMUTED. SEE BELOW FOR */
/*           DETAILS. */

/* lfil    = integer. The fill-in parameter. Each row of L and each row */
/*           of U will have a maximum of lfil elements (excluding the */
/*           diagonal element). lfil must be .ge. 0. */
/*           ** WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO */
/*           EARLIER VERSIONS. */

/* droptol = real*8. Sets the threshold for dropping small terms in the */
/*           factorization. See below for details on dropping strategy. */

/* lfil    = integer. The fill-in parameter. Each row of L and */
/*           each row of U will have a maximum of lfil elements. */
/*           WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO */
/*           EARLIER VERSIONS. */
/*           lfil must be .ge. 0. */

/* permtol = tolerance ratio used to  determne whether or not to permute */
/*           two columns.  At step i columns i and j are permuted when */

/*                     abs(a(i,j))*permtol .gt. abs(a(i,i)) */

/*           [0 --> never permute; good values 0.1 to 0.01] */

/* mbloc   = if desired, permuting can be done only within the diagonal */
/*           blocks of size mbloc. Useful for PDE problems with several */
/*           degrees of freedom.. If feature not wanted take mbloc=n. */


/* iwk     = integer. The lengths of arrays alu and jlu. If the arrays */
/*           are not big enough to store the ILU factorizations, ilut */
/*           will stop with an error message. */

/* On return: */
/* =========== */

/* alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing */
/*           the L and U factors together. The diagonal (stored in */
/*           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix */
/*           contains the i-th row of L (excluding the diagonal entry=1) */
/*           followed by the i-th row of U. */

/* ju      = integer array of length n containing the pointers to */
/*           the beginning of each row of U in the matrix alu,jlu. */

/* iperm   = contains the permutation arrays. */
/*           iperm(1:n) = old numbers of unknowns */
/*           iperm(n+1:2*n) = reverse permutation = new unknowns. */

/* ierr    = integer. Error message with the following meaning. */
/*           ierr  = 0    --> successful return. */
/*           ierr .gt. 0  --> zero pivot encountered at step number ierr. */
/*           ierr  = -1   --> Error. input matrix may be wrong. */
/*                            (The elimination process has generated a */
/*                            row in L or U whose length is .gt.  n.) */
/*           ierr  = -2   --> The matrix L overflows the array al. */
/*           ierr  = -3   --> The matrix U overflows the array alu. */
/*           ierr  = -4   --> Illegal value for lfil. */
/*           ierr  = -5   --> zero row encountered. */

/* work arrays: */
/* ============= */
/* jw      = integer work array of length 2*n. */
/* w       = real work array of length n */

/* IMPORTANR NOTE: */
/* -------------- */
/* TO AVOID PERMUTING THE SOLUTION VECTORS ARRAYS FOR EACH LU-SOLVE, */
/* THE MATRIX A IS PERMUTED ON RETURN. [all column indices are */
/* changed]. SIMILARLY FOR THE U MATRIX. */
/* To permute the matrix back to its original state use the loop: */

/*      do k=ia(1), ia(n+1)-1 */
/*         ja(k) = iperm(ja(k)) */
/*      enddo */

/* ----------------------------------------------------------------------- */
/*     local variables */


    /* Parameter adjustments */
    --iperm;
    --jw;
    --w;
    --ju;
    --ia;
    --a;
    --ja;
    --alu;
    --jlu;

    /* Function Body */
    if (*lfil < 0) {
	goto L998;
    }
/* ----------------------------------------------------------------------- */
/*     initialize ju0 (points to next element to be added to alu,jlu) */
/*     and pointer array. */
/* ----------------------------------------------------------------------- */
    ju0 = *n + 2;
    jlu[1] = ju0;

/*  integer double pointer array. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	jw[*n + j] = 0;
	iperm[j] = j;
	iperm[*n + j] = j;
/* L1: */
    }
/* ----------------------------------------------------------------------- */
/*     beginning of main loop. */
/* ----------------------------------------------------------------------- */
    i__1 = *n;
    for (ii = 1; ii <= i__1; ++ii) {
	j1 = ia[ii];
	j2 = ia[ii + 1] - 1;
	tnorm = 0.;
	i__2 = j2;
	for (k = j1; k <= i__2; ++k) {
	    tnorm += (d__1 = a[k], fabs(d__1));
/* L501: */
	}
	if (tnorm == 0.) {
	    goto L999;
	}
	tnorm /= j2 - j1 + 1;

/*     unpack L-part and U-part of row of A in arrays  w  -- */

	lenu = 1;
	lenl = 0;
	jw[ii] = ii;
	w[ii] = 0.f;
	jw[*n + ii] = ii;

	i__2 = j2;
	for (j = j1; j <= i__2; ++j) {
	    k = iperm[*n + ja[j]];
	    t = a[j];
	    if (k < ii) {
		++lenl;
		jw[lenl] = k;
		w[lenl] = t;
		jw[*n + k] = lenl;
	    } else if (k == ii) {
		w[ii] = t;
	    } else {
		++lenu;
		jpos = ii + lenu - 1;
		jw[jpos] = k;
		w[jpos] = t;
		jw[*n + k] = jpos;
	    }
/* L170: */
	}
	jj = 0;
	len = 0;

/*     eliminate previous rows */

L150:
	++jj;
	if (jj > lenl) {
	    goto L160;
	}
/* ----------------------------------------------------------------------- */
/*     in order to do the elimination in the correct order we must select */
/*     the smallest column index among jw(k), k=jj+1, ..., lenl. */
/* ----------------------------------------------------------------------- */
	jrow = jw[jj];
	k = jj;

/*     determine smallest column index */

	i__2 = lenl;
	for (j = jj + 1; j <= i__2; ++j) {
	    if (jw[j] < jrow) {
		jrow = jw[j];
		k = j;
	    }
/* L151: */
	}

	if (k != jj) {
/*     exchange in jw */
	    j = jw[jj];
	    jw[jj] = jw[k];
	    jw[k] = j;
/*     exchange in jr */
	    jw[*n + jrow] = jj;
	    jw[*n + j] = k;
/*     exchange in w */
	    s = w[jj];
	    w[jj] = w[k];
	    w[k] = s;
	}

/*     zero out element in row by resetting jw(n+jrow) to zero. */

	jw[*n + jrow] = 0;

/*     get the multiplier for row to be eliminated: jrow */

	fact = w[jj] * alu[jrow];

/*     drop term if small */

	if (fabs(fact) <= *droptol) {
	    goto L150;
	}

/*     combine current row and row jrow */

	i__2 = jlu[jrow + 1] - 1;
	for (k = ju[jrow]; k <= i__2; ++k) {
	    s = fact * alu[k];
/*     new column number */
	    j = iperm[*n + jlu[k]];
	    jpos = jw[*n + j];
	    if (j >= ii) {

/*     dealing with upper part. */

		if (jpos == 0) {

/*     this is a fill-in element */

		    ++lenu;
		    i__ = ii + lenu - 1;
		    if (lenu > *n) {
			goto L995;
		    }
		    jw[i__] = j;
		    jw[*n + j] = i__;
		    w[i__] = -s;
		} else {
/*     no fill-in element -- */
		    w[jpos] -= s;
		}
	    } else {

/*     dealing with lower part. */

		if (jpos == 0) {

/*     this is a fill-in element */

		    ++lenl;
		    if (lenl > *n) {
			goto L995;
		    }
		    jw[lenl] = j;
		    jw[*n + j] = lenl;
		    w[lenl] = -s;
		} else {

/*     this is not a fill-in element */

		    w[jpos] -= s;
		}
	    }
/* L203: */
	}

/*     store this pivot element -- (from left to right -- no danger of */
/*     overlap with the working elements in L (pivots). */

	++len;
	w[len] = fact;
	jw[len] = jrow;
	goto L150;
L160:

/*     reset double-pointer to zero (U-part) */

	i__2 = lenu;
	for (k = 1; k <= i__2; ++k) {
	    jw[*n + jw[ii + k - 1]] = 0;
/* L308: */
	}

/*     update L-matrix */

	lenl = len;
	len = MIN(lenl,*lfil);

/*     sort by quick-split */

	_daskr_qsplit_(&w[1], &jw[1], &lenl, &len);

/*     store L-part -- in original coordinates .. */

	i__2 = len;
	for (k = 1; k <= i__2; ++k) {
	    if (ju0 > *iwk) {
		goto L996;
	    }
	    alu[ju0] = w[k];
	    jlu[ju0] = iperm[jw[k]];
	    ++ju0;
/* L204: */
	}

/*     save pointer to beginning of row ii of U */

	ju[ii] = ju0;

/*     update U-matrix -- first apply dropping strategy */

	len = 0;
	i__2 = lenu - 1;
	for (k = 1; k <= i__2; ++k) {
	    if ((d__1 = w[ii + k], fabs(d__1)) > *droptol * tnorm) {
		++len;
		w[ii + len] = w[ii + k];
		jw[ii + len] = jw[ii + k];
	    }
	}
	lenu = len + 1;
	len = MIN(lenu,*lfil);
	if (lenu > 1) {
	    i__2 = lenu - 1;
	    _daskr_qsplit_(&w[ii + 1], &jw[ii + 1], &i__2, &len);
	}

/*     determine next pivot -- */

	imax = ii;
	xmax = (d__1 = w[imax], fabs(d__1));
	xmax0 = xmax;
	icut = ii - 1 + *mbloc - (ii - 1) % *mbloc;
	i__2 = ii + len - 1;
	for (k = ii + 1; k <= i__2; ++k) {
	    t = (d__1 = w[k], fabs(d__1));
	    if (t > xmax && t * *permtol > xmax0 && jw[k] <= icut) {
		imax = k;
		xmax = t;
	    }
	}

/*     exchange w's */

	tmp = w[ii];
	w[ii] = w[imax];
	w[imax] = tmp;

/*     update iperm and reverse iperm */

	j = jw[imax];
	i__ = iperm[ii];
	iperm[ii] = iperm[j];
	iperm[j] = i__;

/*     reverse iperm */

	iperm[*n + iperm[ii]] = ii;
	iperm[*n + iperm[j]] = j;
/* ----------------------------------------------------------------------- */

	if (len + ju0 > *iwk) {
	    goto L997;
	}

/*     copy U-part in original coordinates */

	i__2 = ii + len - 1;
	for (k = ii + 1; k <= i__2; ++k) {
	    jlu[ju0] = iperm[jw[k]];
	    alu[ju0] = w[k];
	    ++ju0;
/* L302: */
	}

/*     store inverse of diagonal element of u */

	if (w[ii] == 0.) {
	    w[ii] = (*droptol + 1e-4) * tnorm;
	}
	alu[ii] = 1. / w[ii];

/*     update pointer to beginning of next row of U. */

	jlu[ii + 1] = ju0;
/* ----------------------------------------------------------------------- */
/*     end main loop */
/* ----------------------------------------------------------------------- */
/* L500: */
    }

/*     permute all column indices of LU ... */

    i__1 = jlu[*n + 1] - 1;
    for (k = jlu[1]; k <= i__1; ++k) {
	jlu[k] = iperm[*n + jlu[k]];
    }

/*     ...and of A */

    i__1 = ia[*n + 1] - 1;
    for (k = ia[1]; k <= i__1; ++k) {
	ja[k] = iperm[*n + ja[k]];
    }

    *ierr = 0;
    return 0;

/*     incomprehensible error. Matrix must be wrong. */

L995:
    *ierr = -1;
    return 0;

/*     insufficient storage in L. */

L996:
    *ierr = -2;
    return 0;

/*     insufficient storage in U. */

L997:
    *ierr = -3;
    return 0;

/*     illegal lfil entered. */

L998:
    *ierr = -4;
    return 0;

/*     zero row encountered */

L999:
    *ierr = -5;
    return 0;
/* ----------------end-of-ilutp------------------------------------------- */
/* ----------------------------------------------------------------------- */
} /* ilutp_ */

/* Subroutine */ int _daskr_qsplit_(real_number *a, integer *ind, integer *n, integer
	*ncut)
{
    /* System generated locals */
    integer i__1;
    real_number d__1;

    /* Local variables */
    static integer j, mid;
    static real_number tmp;
    static integer last, itmp, first;
    static real_number abskey;

/* ----------------------------------------------------------------------- */
/*     does a quick-sort split of a real array. */
/*     on input a(1:n). is a real array */
/*     on output a(1:n) is permuted such that its elements satisfy: */

/*     abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and */
/*     abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut */

/*     ind(1:n) is an integer array which permuted in the same way as a(*). */
/* ----------------------------------------------------------------------- */
/* ----- */
    /* Parameter adjustments */
    --ind;
    --a;

    /* Function Body */
    first = 1;
    last = *n;
    if (*ncut < first || *ncut > last) {
	return 0;
    }

/*     outer loop -- while mid .ne. ncut do */

L1:
    mid = first;
    abskey = (d__1 = a[mid], fabs(d__1));
    i__1 = last;
    for (j = first + 1; j <= i__1; ++j) {
	if ((d__1 = a[j], fabs(d__1)) > abskey) {
	    ++mid;
/*     interchange */
	    tmp = a[mid];
	    itmp = ind[mid];
	    a[mid] = a[j];
	    ind[mid] = ind[j];
	    a[j] = tmp;
	    ind[j] = itmp;
	}
/* L2: */
    }

/*     interchange */

    tmp = a[mid];
    a[mid] = a[first];
    a[first] = tmp;

    itmp = ind[mid];
    ind[mid] = ind[first];
    ind[first] = itmp;

/*     test for while loop */

    if (mid == *ncut) {
	return 0;
    }
    if (mid > *ncut) {
	last = mid - 1;
    } else {
	first = mid + 1;
    }
    goto L1;
/* ----------------end-of-qsplit------------------------------------------ */
    return 0;
/* ----------------------------------------------------------------------- */
} /* qsplit_ */

/* Subroutine */ int _daskr_lusol_(integer *n, real_number *y, real_number *x,
	real_number *alu, integer *jlu, integer *ju)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, k;

/* ----------------------------------------------------------------------- */

/* This routine solves the system (LU) x = y, */
/* given an LU decomposition of a matrix stored in (alu, jlu, ju) */
/* modified sparse row format */

/* ----------------------------------------------------------------------- */
/* on entry: */
/* n   = dimension of system */
/* y   = the right-hand-side vector */
/* alu, jlu, ju */
/*     = the LU matrix as provided from the ILU routines. */

/* on return */
/* x   = solution of LU x = y. */
/* ----------------------------------------------------------------------- */

/* Note: routine is in place: call lusol (n, x, x, alu, jlu, ju) */
/*       will solve the system with rhs x and overwrite the result on x . */

/* ----------------------------------------------------------------------- */
/* local variables */


/* forward solve */

    /* Parameter adjustments */
    --x;
    --y;
    --alu;
    --jlu;
    --ju;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = y[i__];
	i__2 = ju[i__] - 1;
	for (k = jlu[i__]; k <= i__2; ++k) {
	    x[i__] -= alu[k] * x[jlu[k]];
/* L41: */
	}
/* L40: */
    }

/*     backward solve. */

    for (i__ = *n; i__ >= 1; --i__) {
	i__1 = jlu[i__ + 1] - 1;
	for (k = ju[i__]; k <= i__1; ++k) {
	    x[i__] -= alu[k] * x[jlu[k]];
/* L91: */
	}
	x[i__] = alu[i__] * x[i__];
/* L90: */
    }

    return 0;
/* ----------------end of lusol ------------------------------------------ */
/* ----------------------------------------------------------------------- */
} /* lusol_ */

/* ----------------------------------------------------------------------c */
/*                          S P A R S K I T                             c */
/* ----------------------------------------------------------------------c */
/*               REORDERING ROUTINES -- LEVEL SET BASED ROUTINES        c */
/* ----------------------------------------------------------------------c */
/* dblstr   : doubled stripe partitioner */
/* BFS      : Breadth-First search traversal algorithm */
/* add_lvst : routine to add a level -- used by BFS */
/* stripes  : finds the level set structure */
/* perphn   : finds a pseudo-peripheral node and performs a BFS from it. */
/* rversp   : routine to reverse a given permutation (e.g., for RCMK) */
/* maskdeg  : integer function to compute the `masked' of a node */
/* ----------------------------------------------------------------------- */
/* Subroutine */ int _daskr_dblstr_(integer *n, integer *ja, integer *ia, integer *
	ip1, integer *ip2, integer *nfirst, integer *riord, integer *ndom, 
	integer *map, integer *mapptr, integer *mask, integer *levels, 
	integer *iwk)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j, k;
    extern /* Subroutine */ int _daskr_bfs_(integer *, integer *, integer *, integer
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *);
    static integer ndp1, idom, jdom, kdom, init, nlev;
    extern /* Subroutine */ int _daskr_perphn_(integer *, integer *, integer *,
	    integer *, integer *, integer *, integer *, integer *, integer *);
    static integer numnod, maskval, nextdom;
    extern /* Subroutine */ int _daskr_stripes_(integer *, integer *, integer *,
	    integer *, integer *, integer *, integer *);

/* ----------------------------------------------------------------------- */
/*     this routine does a two-way partitioning of a graph using */
/*     level sets recursively. First a coarse set is found by a */
/*     simple cuthill-mc Kee type algorithm. Them each of the large */
/*     domains is further partitioned into subsets using the same */
/*     technique. The ip1 and ip2 parameters indicate the desired number */
/*     number of partitions 'in each direction'. So the total number of */
/*     partitions on return ought to be equal (or close) to ip1*ip2 */
/* ----------------------parameters---------------------------------------- */
/* on entry: */
/* --------- */
/* n      = row dimension of matrix == number of vertices in graph */
/* ja, ia = pattern of matrix in CSR format (the ja,ia arrays of csr data */
/*          structure) */
/* ip1    = integer indicating the number of large partitions ('number of */
/*          paritions in first direction') */
/* ip2    = integer indicating the number of smaller partitions, per */
/*          large partition, ('number of partitions in second direction') */
/* nfirst = number of nodes in the first level that is input in riord */
/* riord  = (also an ouput argument). on entry riord contains the labels */
/*          of the nfirst nodes that constitute the first level. */
/* on return: */
/* ----------- */
/* ndom   = total number of partitions found */
/* map    = list of nodes listed partition by partition from partition 1 */
/*          to paritition ndom. */
/* mapptr = pointer array for map. All nodes from position */
/*          k1=mapptr(idom),to position k2=mapptr(idom+1)-1 in map belong */
/*          to partition idom. */
/* work arrays: */
/* ------------- */
/* mask   = array of length n, used to hold the partition number of each */
/*          node for the first (large) partitioning. */
/*          mask is also used as a marker of  visited nodes. */
/* levels = integer array of length .le. n used to hold the pointer */
/*          arrays for the various level structures obtained from BFS. */

/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --iwk;
    --levels;
    --mask;
    --mapptr;
    --map;
    --riord;
    --ia;
    --ja;

    /* Function Body */
    maskval = 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	mask[j] = maskval;
    }
    iwk[1] = 0;
    _daskr_bfs_(n, &ja[1], &ia[1], nfirst, &iwk[1], &mask[1], &maskval, &riord[1], &
	    levels[1], &nlev);

/*     init = riord(1) */
/*     call perphn (ja,ia,mask,maskval,init,nlev,riord,levels) */
    _daskr_stripes_(&nlev, &riord[1], &levels[1], ip1, &map[1], &mapptr[1], ndom);
/* ----------------------------------------------------------------------- */
    if (*ip2 == 1) {
	return 0;
    }
    ndp1 = *ndom + 1;

/*     pack info into array iwk */

    i__1 = *ndom + 1;
    for (j = 1; j <= i__1; ++j) {
	iwk[j] = ndp1 + mapptr[j];
    }
    i__1 = mapptr[*ndom + 1] - 1;
    for (j = 1; j <= i__1; ++j) {
	iwk[ndp1 + j] = map[j];
    }
    i__1 = *ndom;
    for (idom = 1; idom <= i__1; ++idom) {
	j = iwk[idom];
	numnod = iwk[idom + 1] - iwk[idom];
	init = iwk[j];
	i__2 = iwk[idom + 1] - 1;
	for (k = j; k <= i__2; ++k) {
	}
    }
    i__1 = *ndom;
    for (idom = 1; idom <= i__1; ++idom) {
	i__2 = mapptr[idom + 1] - 1;
	for (k = mapptr[idom]; k <= i__2; ++k) {
	    mask[map[k]] = idom;
	}
    }
    nextdom = 1;

/*     jdom = counter for total number of (small) subdomains */

    jdom = 1;
    mapptr[jdom] = 1;
/* ----------------------------------------------------------------------- */
    i__1 = *ndom;
    for (idom = 1; idom <= i__1; ++idom) {
	maskval = idom;
	*nfirst = 1;
	numnod = iwk[idom + 1] - iwk[idom];
	j = iwk[idom];
	init = iwk[j];
	nextdom = mapptr[jdom];
/*  note:    old version uses iperm array */
	_daskr_perphn_(&numnod, &ja[1], &ia[1], &init, &mask[1], &maskval, &nlev, &
		riord[1], &levels[1]);

	_daskr_stripes_(&nlev, &riord[1], &levels[1], ip2, &map[nextdom], &mapptr[
		jdom], &kdom);

	mapptr[jdom] = nextdom;
	i__2 = jdom + kdom - 1;
	for (j = jdom; j <= i__2; ++j) {
	    mapptr[j + 1] = nextdom + mapptr[j + 1] - 1;
	}
	jdom += kdom;
    }

    *ndom = jdom - 1;
    return 0;
} /* dblstr_ */

/* Subroutine */ int _daskr_bfs_(integer *n, integer *ja, integer *ia, integer *
	nfirst, integer *iperm, integer *mask, integer *maskval, integer *
	riord, integer *levels, integer *nlev)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int _daskr_add_lvst__(integer *, integer *, integer *,
	    integer *, integer *, integer *, integer *, integer *);
    static integer j, ii, nod, iend, istart;
    static integer permut;

/* ----------------------------------------------------------------------- */
/* finds the level-structure (breadth-first-search or CMK) ordering for a */
/* given sparse matrix. Uses add_lvst. Allows an set of nodes to be */
/* the initial level (instead of just one node). */
/* -------------------------parameters------------------------------------ */
/* on entry: */
/* --------- */
/*     n      = number of nodes in the graph */
/*     ja, ia = pattern of matrix in CSR format (the ja,ia arrays of csr data */
/*     structure) */
/*     nfirst = number of nodes in the first level that is input in riord */
/*     iperm  = integer array indicating in which order to  traverse the graph */
/*     in order to generate all connected components. */
/*     if iperm(1) .eq. 0 on entry then BFS will traverse the nodes */
/*     in the  order 1,2,...,n. */

/*     riord  = (also an ouput argument). On entry riord contains the labels */
/*     of the nfirst nodes that constitute the first level. */

/*     mask   = array used to indicate whether or not a node should be */
/*     condidered in the graph. see maskval. */
/*     mask is also used as a marker of  visited nodes. */

/*     maskval= consider node i only when:  mask(i) .eq. maskval */
/*     maskval must be .gt. 0. */
/*     thus, to consider all nodes, take mask(1:n) = 1. */
/*     maskval=1 (for example) */

/*     on return */
/*     --------- */
/*     mask   = on return mask is restored to its initial state. */
/*     riord  = `reverse permutation array'. Contains the labels of the nodes */
/*     constituting all the levels found, from the first level to */
/*     the last. */
/*     levels = pointer array for the level structure. If lev is a level */
/*     number, and k1=levels(lev),k2=levels(lev+1)-1, then */
/*     all the nodes of level number lev are: */
/*     riord(k1),riord(k1+1),...,riord(k2) */
/*     nlev   = number of levels found */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --mask;
    --iperm;
    --ja;
    --ia;
    --riord;
    --levels;

    /* Function Body */
    permut = iperm[1] != 0;

/*     start pointer structure to levels */

    *nlev = 0;

/*     previous end */

    istart = 0;
    ii = 0;

/*     current end */

    iend = *nfirst;

/*     intialize masks to zero -- except nodes of first level -- */

    i__1 = *nfirst;
    for (j = 1; j <= i__1; ++j) {
	mask[riord[j]] = 0;
/* L12: */
    }
/* ----------------------------------------------------------------------- */
/* L13: */

L1:
    ++(*nlev);
    levels[*nlev] = istart + 1;
    _daskr_add_lvst__(&istart, &iend, nlev, &riord[1], &ja[1], &ia[1], &mask[1],
	    maskval);
    if (istart < iend) {
	goto L1;
    }
L2:
    ++ii;
    if (ii <= *n) {
	nod = ii;
	if (permut) {
	    nod = iperm[nod];
	}
	if (mask[nod] == *maskval) {

/*     start a new level */

	    istart = iend;
	    ++iend;
	    riord[iend] = nod;
	    mask[nod] = 0;
	    goto L1;
	} else {
	    goto L2;
	}
    }
/* ----------------------------------------------------------------------- */
/* L3: */
    levels[*nlev + 1] = iend + 1;
    i__1 = iend;
    for (j = 1; j <= i__1; ++j) {
	mask[riord[j]] = *maskval;
    }
/* ----------------------------------------------------------------------- */
    return 0;
} /* bfs_ */

/* Subroutine */ int _daskr_add_lvst__(integer *istart, integer *iend, integer *nlev,
	 integer *riord, integer *ja, integer *ia, integer *mask, integer *
	maskval)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, ir, nod;

/* ------------------------------------------------------------- */
/*     adds one level set to the previous sets.. */
/*     span all nodes of previous mask */
/* ------------------------------------------------------------- */
    /* Parameter adjustments */
    --mask;
    --ia;
    --ja;
    --riord;

    /* Function Body */
    nod = *iend;
    i__1 = *iend;
    for (ir = *istart + 1; ir <= i__1; ++ir) {
	i__ = riord[ir];
	i__2 = ia[i__ + 1] - 1;
	for (k = ia[i__]; k <= i__2; ++k) {
	    j = ja[k];
	    if (mask[j] == *maskval) {
		++nod;
		mask[j] = 0;
		riord[nod] = j;
	    }
/* L24: */
	}
/* L25: */
    }
    *istart = *iend;
    *iend = nod;
    return 0;
} /* add_lvst__ */

/* Subroutine */ int _daskr_stripes_(integer *nlev, integer *riord, integer *levels,
	integer *ip, integer *map, integer *mapptr, integer *ndom)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer k, ib, ktr, ilev, nsiz, psiz;

/* ----------------------------------------------------------------------- */
/*    this is a post processor to BFS. stripes uses the output of BFS to */
/*    find a decomposition of the adjacency graph by stripes. It fills */
/*    the stripes level by level until a number of nodes .gt. ip is */
/*    is reached. */
/* ---------------------------parameters----------------------------------- */
/* on entry: */
/* -------- */
/* nlev   = number of levels as found by BFS */
/* riord  = reverse permutation array produced by BFS -- */
/* levels = pointer array for the level structure as computed by BFS. If */
/*          lev is a level number, and k1=levels(lev),k2=levels(lev+1)-1, */
/*          then all the nodes of level number lev are: */
/*                      riord(k1),riord(k1+1),...,riord(k2) */
/*  ip    = number of desired partitions (subdomains) of about equal size. */

/* on return */
/* --------- */
/* ndom     = number of subgraphs (subdomains) found */
/* map      = node per processor list. The nodes are listed contiguously */
/*            from proc 1 to nproc = mpx*mpy. */
/* mapptr   = pointer array for array map. list for proc. i starts at */
/*            mapptr(i) and ends at mapptr(i+1)-1 in array map. */
/* ----------------------------------------------------------------------- */
/* local variables. */

    /* Parameter adjustments */
    --levels;
    --riord;
    --map;
    --mapptr;

    /* Function Body */
    *ndom = 1;
    ib = 1;
/* to add: if (ip .le. 1) then ... */
    nsiz = levels[*nlev + 1] - levels[1];
/* Computing MAX */
    i__1 = 1, i__2 = *ip - *ndom + 1;
    psiz = (nsiz - ib) / MAX(i__1,i__2) + 1;
    mapptr[*ndom] = ib;
    ktr = 0;
    i__1 = *nlev;
    for (ilev = 1; ilev <= i__1; ++ilev) {

/*     add all nodes of this level to domain */

	i__2 = levels[ilev + 1] - 1;
	for (k = levels[ilev]; k <= i__2; ++k) {
	    map[ib] = riord[k];
	    ++ib;
	    ++ktr;
	    if (ktr >= psiz || k >= nsiz) {
		++(*ndom);
		mapptr[*ndom] = ib;
/* Computing MAX */
		i__3 = 1, i__4 = *ip - *ndom + 1;
		psiz = (nsiz - ib) / MAX(i__3,i__4) + 1;
		ktr = 0;
	    }

/* L3: */
	}
/* L10: */
    }
    --(*ndom);
    return 0;
} /* stripes_ */

integer _daskr_maskdeg_(integer *ja, integer *ia, integer *nod, integer *mask,
	integer *maskval)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static integer k, deg;

/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --mask;
    --ia;
    --ja;

    /* Function Body */
    deg = 0;
    i__1 = ia[*nod + 1] - 1;
    for (k = ia[*nod]; k <= i__1; ++k) {
	if (mask[ja[k]] == *maskval) {
	    ++deg;
	}
    }
    ret_val = deg;
    return ret_val;
} /* maskdeg_ */

/* Subroutine */ int _daskr_perphn_(integer *n, integer *ja, integer *ia, integer *
	init, integer *mask, integer *maskval, integer *nlev, integer *riord, 
	integer *levels)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, deg;
    extern /* Subroutine */ int _daskr_bfs_(integer *, integer *, integer *, integer
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *);
    static integer nod, iperm[1], nlevp, mindeg, nfirst;
    extern integer _daskr_maskdeg_(integer *, integer *, integer *, integer *,
	    integer *);

/* ----------------------------------------------------------------------- */
/*     finds a peripheral node and does a BFS search from it. */
/* ----------------------------------------------------------------------- */
/*     see routine  dblstr for description of parameters */
/* input: */
/* ------- */
/* ja, ia  = list pointer array for the adjacency graph */
/* mask    = array used for masking nodes -- see maskval */
/* maskval = value to be checked against for determing whether or */
/*           not a node is masked. If mask(k) .ne. maskval then */
/*           node k is not considered. */
/* init    = init node in the pseudo-peripheral node algorithm. */

/* output: */
/* ------- */
/* init    = actual pseudo-peripherial node found. */
/* nlev    = number of levels in the final BFS traversal. */
/* riord   = */
/* levels  = */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --levels;
    --riord;
    --mask;
    --ia;
    --ja;

    /* Function Body */
    nlevp = 0;
L1:
    riord[1] = *init;
    nfirst = 1;
    iperm[0] = 0;

    _daskr_bfs_(n, &ja[1], &ia[1], &nfirst, iperm, &mask[1], maskval, &riord[1], &
	    levels[1], nlev);
    if (*nlev > nlevp) {
	mindeg = *n + 1;
	i__1 = levels[*nlev + 1] - 1;
	for (j = levels[*nlev]; j <= i__1; ++j) {
	    nod = riord[j];
	    deg = _daskr_maskdeg_(&ja[1], &ia[1], &nod, &mask[1], maskval);
	    if (deg < mindeg) {
		*init = nod;
		mindeg = deg;
	    }
	}
	nlevp = *nlev;
	goto L1;
    }
    return 0;
} /* perphn_ */

/* Subroutine */ int _daskr_rversp_(integer *n, integer *riord)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, k;

/* ----------------------------------------------------------------------- */
/*     this routine does an in-place reversing of the permutation array */
/*     riord -- */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --riord;

    /* Function Body */
    i__1 = *n / 2;
    for (j = 1; j <= i__1; ++j) {
	k = riord[j];
	riord[j] = riord[*n - j + 1];
	riord[*n - j + 1] = k;
/* L26: */
    }
    return 0;
} /* rversp_ */

/* ----------------------------------------------------------------------c */
/*     Non-SPARSKIT utility routine */
/* ----------------------------------------------------------------------c */
/* Subroutine */ int _daskr_atob_(integer *n, real_number *a, integer *ja, integer *
	ia, real_number *b, integer *jb, integer *ib)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/* ... Copy matrix a,ja,ia to b,jb,ib.  Both matrices are in */
/*     compressed sparse row format. */
/* ... Input arguments: */
/* ... Output arguments: */
/* ... Local variable: */
    /* Parameter adjustments */
    --ib;
    --ia;
    --a;
    --ja;
    --b;
    --jb;

    /* Function Body */
    i__1 = ia[*n + 1] - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b[i__] = a[i__];
	jb[i__] = ja[i__];
    }
    i__1 = *n + 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ib[i__] = ia[i__];
    }
    return 0;
/*  end of atob */
} /* atob_ */

