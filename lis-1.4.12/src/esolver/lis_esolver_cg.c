/* Copyright (C) The Scalable Software Infrastructure Project. All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:
   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
   2. Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
   3. Neither the name of the project nor the names of its contributors 
      may be used to endorse or promote products derived from this software 
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE SCALABLE SOFTWARE INFRASTRUCTURE PROJECT
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE SCALABLE SOFTWARE INFRASTRUCTURE
   PROJECT BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
   OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

#ifdef HAVE_CONFIG_H
  #include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
  #include "lis_config_win.h"
#endif
#endif

#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_MALLOC_H
        #include <malloc.h>
#endif
#include <math.h>
#include <string.h>
#include <stdarg.h>
#ifdef USE_SSE2
  #include <emmintrin.h>
#endif
#ifdef _OPENMP
  #include <omp.h>
#endif
#ifdef USE_MPI
  #include <mpi.h>
#endif
#include "lislib.h"

/***************************************
 * Conjugate Gradient Method           *
 ***************************************
 x(0)    = (1,...,1)^T
 ***************************************
 for k=0,1,...
   mu(k)   = <x(k),x(k)> / <x(k),A * x(k)>
   r         = x(k) - mu(k)Ax(k)
   w         = M^-1 * r
   use the Rayleigh-Ritz method for I - mu * A on 
   span {w, x(k), x(k-1)}
   x(k+1)    = Ritz vector correspoinding to the maximal Ritz value
 ***************************************/

#define NWORK 9
#undef __FUNC__
#define __FUNC__ "lis_ecg_check_params"
LIS_INT lis_ecg_check_params(LIS_ESOLVER esolver)
{
  LIS_DEBUG_FUNC_IN;
  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_ecg_malloc_work"
LIS_INT lis_ecg_malloc_work(LIS_ESOLVER esolver)
{
  LIS_VECTOR  *work;
  LIS_INT      i,j,worklen,err;

  LIS_DEBUG_FUNC_IN;

  worklen = NWORK;

  work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_ecg_malloc_work::work" );
  if( work==NULL )
  {
    LIS_SETERR_MEM(worklen*sizeof(LIS_VECTOR));
    return LIS_ERR_OUT_OF_MEMORY;
  }
  if( esolver->eprecision==LIS_PRECISION_DEFAULT )
  {
    for(i=0;i<worklen;i++)
    {
      err = lis_vector_duplicate(esolver->A,&work[i]);
      if( err ) break;
    }
  }
  else
  {
    for(i=0;i<worklen;i++)
    {
      err = lis_vector_duplicateex(LIS_PRECISION_QUAD,esolver->A,&work[i]);
      if( err ) break;
    }
  }
  if( i<worklen )
  {
    for(j=0;j<i;j++) lis_vector_destroy(work[j]);
    lis_free(work);
    return err;
  }
  esolver->worklen = worklen;
  esolver->work    = work;

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_ecg"
LIS_INT lis_ecg(LIS_ESOLVER esolver)
{
  LIS_MATRIX        A;
  LIS_VECTOR        x;
  LIS_SCALAR        evalue;
  LIS_INT               emaxiter;
  LIS_REAL          tol;
  LIS_INT               iter,iter3,nsolver,i,j,output;
  LIS_INT               nprocs,my_rank;
  LIS_REAL          nrm2,resid,resid3;
  LIS_SCALAR        lshift;
  LIS_VECTOR        b,D,r,w,p,Aw,Ax,Ap,ones,Ds;
  LIS_SCALAR        *SA, *SB, *SW, *v3, *SAv3, *SBv3, *z3, *q3, *SBz3, evalue3, ievalue3;
  LIS_SOLVER        solver;
  LIS_PRECON        precon;
  LIS_MATRIX        A0;
  LIS_VECTOR        x0,z,q;
  double      times,itimes,ptimes,p_c_times,p_i_times;
  LIS_INT           nsol, precon_type;
  char              solvername[128], preconname[128];

  A = esolver->A;
  x = esolver->x;
  
  emaxiter = esolver->options[LIS_EOPTIONS_MAXITER];
  tol = esolver->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN]; 
  output  = esolver->options[LIS_EOPTIONS_OUTPUT];
  lshift = esolver->lshift;

#ifdef _LONG__DOUBLE
  if( A->my_rank==0 ) printf("local shift = %Le\n", lshift);
#else
  if( A->my_rank==0 ) printf("local shift = %e\n", lshift);
#endif
  if (lshift != 0) lis_matrix_shift_diagonal(A, lshift);

  SA = (LIS_SCALAR *)lis_malloc(3*3*sizeof(LIS_SCALAR), "lis_ecg::SA");
  SB = (LIS_SCALAR *)lis_malloc(3*3*sizeof(LIS_SCALAR), "lis_ecg::SB");
  SW = (LIS_SCALAR *)lis_malloc(3*3*sizeof(LIS_SCALAR), "lis_ecg::SW");
  v3 = (LIS_SCALAR *)lis_malloc(3*sizeof(LIS_SCALAR), "lis_ecg::v3");
  SAv3 = (LIS_SCALAR *)lis_malloc(3*sizeof(LIS_SCALAR), "lis_ecg::SAv3");
  SBv3 = (LIS_SCALAR *)lis_malloc(3*sizeof(LIS_SCALAR), "lis_ecg::SBv3");
  SBz3 = (LIS_SCALAR *)lis_malloc(3*sizeof(LIS_SCALAR), "lis_ecg::SBz3");
  z3 = (LIS_SCALAR *)lis_malloc(3*sizeof(LIS_SCALAR), "lis_ecg::z3");
  q3 = (LIS_SCALAR *)lis_malloc(3*sizeof(LIS_SCALAR), "lis_ecg::q3");

  b = esolver->work[0];
  D = esolver->work[1];
  Ds = esolver->work[2];
  r = esolver->work[3];
  w = esolver->work[4];
  p = esolver->work[5];
  Aw = esolver->work[6];
  Ax = esolver->work[7];
  Ap = esolver->work[8];

  lis_vector_set_all(1.0,b);
  lis_vector_nrm2(b, &nrm2);
  lis_vector_scale(1/nrm2, b);
  lis_solver_create(&solver);
  lis_solver_set_option("-i cg -p ilu",solver);
  lis_solver_set_optionC(solver);
  lis_solver_get_solver(solver, &nsol);
  lis_solver_get_precon(solver, &precon_type);
  lis_solver_get_solvername(nsol, solvername);
  lis_solver_get_preconname(precon_type, preconname);
#ifdef _LONGLONG
  if( A->my_rank==0 ) printf("solver     : %s %lld\n", solvername, nsol);
  if( A->my_rank==0 ) printf("precon     : %s %lld\n", preconname, precon_type);
#else
  if( A->my_rank==0 ) printf("solver     : %s %d\n", solvername, nsol);
  if( A->my_rank==0 ) printf("precon     : %s %d\n", preconname, precon_type);
#endif
  lis_solve(A, b, x, solver);
  lis_vector_copy(b,Ax);

  lis_vector_nrm2(x, &nrm2);
  lis_vector_set_all(0.0,p);
  lis_vector_set_all(0.0,Ap);

  lis_precon_create(solver, &precon);
  solver->precon = precon;

  iter=0;

  while (iter<emaxiter)
    {
      iter = iter + 1;

      /* mu = <x, x> / <x, A * x> */
      lis_vector_dot(x,Ax,&evalue);

      /* r = x - mu * A * x */
      lis_vector_axpyz(-(evalue),x,Ax,r); 
      lis_vector_nrm2(r, &nrm2);

      /* convergence check */
      resid = fabs(nrm2/(evalue));

      if( output )
  {
    if( output & LIS_EPRINT_MEM ) esolver->residual[iter] = resid;
    if( output & LIS_EPRINT_OUT && A->my_rank==0 ) lis_print_rhistory(iter,resid);
  }

      if (resid<tol) break;  

      /* w = M^-1 * r */
      lis_psolve(solver, x, w);
      lis_vector_copy(x,Aw);
      lis_vector_nrm2(w, &nrm2);

      /* Rayleigh-Ritz method for I - mu * A on span {w, x(k), x(k-1)} */
      lis_vector_dot(w,Aw,&SA[0]);
      lis_vector_dot(x,Aw,&SA[3]);
      lis_vector_dot(p,Aw,&SA[6]);
      SA[1] = SA[3];
      lis_vector_dot(x,Ax,&SA[4]);
      lis_vector_dot(p,Ax,&SA[7]);
      SA[2] = SA[6];
      SA[5] = SA[7];
      lis_vector_dot(p,Ap,&SA[8]);

      lis_vector_dot(w,w,&SB[0]);
      lis_vector_dot(x,w,&SB[3]);
      lis_vector_dot(p,w,&SB[6]);
      SB[1] = SB[3];
      lis_vector_dot(x,x,&SB[4]);
      lis_vector_dot(p,x,&SB[7]);
      SB[2] = SB[6];
      SB[5] = SB[7];
      lis_vector_dot(p,p,&SB[8]);
      
      lis_array_set_all(3, 1.0, v3);

      iter3=0;
      while (iter3<emaxiter)
  {
    iter3 = iter3 + 1;
    lis_array_nrm2(3, v3, &nrm2); 
    lis_array_scale(3, 1/nrm2, v3);
    lis_array_matvec(3, SB, v3, SBv3, LIS_INS_VALUE);
    lis_array_invvec(3, SA, SBv3, z3);
    lis_array_dot2(3, SBv3, z3, &ievalue3);
    if (ievalue3==0) 
      {
        if( A->my_rank==0 ) printf("ievalue3 is zero\n");
        lis_precon_destroy(precon);
        lis_solver_destroy(solver);
        esolver->iter       = iter;
        esolver->resid      = resid;
        esolver->evalue[0] = evalue;

        if (lshift != 0) lis_matrix_shift_diagonal(A, -lshift);
        lis_free(SA);
        lis_free(SB);
        lis_free(SW);
        lis_free(v3);
        lis_free(SAv3);
        lis_free(SBv3);
        lis_free(SBz3);
        lis_free(z3);
        lis_free(q3);
        return LIS_BREAKDOWN;
      }
    lis_array_axpyz(3, -ievalue3, SBv3, z3, q3);
    lis_array_nrm2(3, q3, &resid3); 
    resid3 = fabs(resid3 / ievalue3);
    if (resid3<1e-12) break;   
    lis_array_copy(3,z3,v3);
  }

      evalue3 = 1 / ievalue3;
      
      lis_vector_scale(v3[0],w);  
      lis_vector_axpy(v3[2],p,w);
      lis_vector_xpay(w,v3[1],x);
      lis_vector_copy(w,p);
      
      lis_vector_scale(v3[0],Aw);  
      lis_vector_axpy(v3[2],Ap,Aw);
      lis_vector_xpay(Aw,v3[1],Ax);
      lis_vector_copy(Aw,Ap);
      
      /* x = Ritz vector correspoinding to the maximal Ritz value */
      lis_vector_nrm2(x,&nrm2);
      lis_vector_scale(1/nrm2,x);
      lis_vector_scale(1/nrm2,Ax);
      
      lis_vector_nrm2(p,&nrm2);
      lis_vector_scale(1/nrm2,p);
      lis_vector_scale(1/nrm2,Ap);
      
      lis_solver_get_timeex(solver,&times,&itimes,&ptimes,&p_c_times,&p_i_times);
      esolver->ptimes += solver->ptimes;
      esolver->itimes += solver->itimes;
      esolver->p_c_times += solver->p_c_times;
      esolver->p_i_times += solver->p_i_times;

    }

  lis_precon_destroy(precon);
  lis_solver_destroy(solver);

  esolver->iter       = iter;
  esolver->resid      = resid;
  esolver->evalue[0] = evalue;

  if (lshift != 0) lis_matrix_shift_diagonal(A, -lshift);

  lis_free(SA);
  lis_free(SB);
  lis_free(SW);
  lis_free(v3);
  lis_free(SAv3);
  lis_free(SBv3);
  lis_free(SBz3);
  lis_free(z3);
  lis_free(q3);

  if (resid<tol) 
    {
      esolver->retcode = LIS_SUCCESS;
      return LIS_SUCCESS;
    }
  else
    {
      esolver->retcode = LIS_MAXITER;
      return LIS_MAXITER;
    }
}

/***************************************
 * Conjugate Residual Method           *
 ***************************************
 x    = (1,...,1)^T
 ***************************************
 lambda = <Ax, x> / <x, x> 
 r = lambda * x - A * x
 while iter < maxiter
   alpha = (<r, Ap> - lambda * <r, p>) / (<Ap, Ap> - 2 * lambda * <Ap, p>)
   x = x + alpha * p 
   lambda = <Ax, x> / <x, x>
   r = lambda * x - A * x 
   beta = - [(<Ar, Ap> - lambda * {<Ar, p> + <Ap. r>} + lambda^2 * <r, p>] / [<Ap, Ap> - 2 * lambda * <Ap, p>] 
   p = r + beta * p
   convergence check
 ***************************************/

#undef NWORK
#define NWORK 5
#undef __FUNC__
#define __FUNC__ "lis_ecr_check_params"
LIS_INT lis_ecr_check_params(LIS_ESOLVER esolver)
{
  LIS_DEBUG_FUNC_IN;
  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_ecr_malloc_work"
LIS_INT lis_ecr_malloc_work(LIS_ESOLVER esolver)
{
  LIS_VECTOR  *work;
  LIS_INT      i,j,worklen,err;

  LIS_DEBUG_FUNC_IN;

  worklen = NWORK;

  work    = (LIS_VECTOR *)lis_malloc( worklen*sizeof(LIS_VECTOR),"lis_ecr_malloc_work::work" );
  if( work==NULL )
  {
    LIS_SETERR_MEM(worklen*sizeof(LIS_VECTOR));
    return LIS_ERR_OUT_OF_MEMORY;
  }
  if( esolver->eprecision==LIS_PRECISION_DEFAULT )
  {
    for(i=0;i<worklen;i++)
    {
      err = lis_vector_duplicate(esolver->A,&work[i]);
      if( err ) break;
    }
  }
  else
  {
    for(i=0;i<worklen;i++)
    {
      err = lis_vector_duplicateex(LIS_PRECISION_QUAD,esolver->A,&work[i]);
      if( err ) break;
    }
  }
  if( i<worklen )
  {
    for(j=0;j<i;j++) lis_vector_destroy(work[j]);
    lis_free(work);
    return err;
  }
  esolver->worklen = worklen;
  esolver->work    = work;

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_ecr"
LIS_INT lis_ecr(LIS_ESOLVER esolver)
{
  LIS_MATRIX        A;
  LIS_VECTOR        x;
  LIS_SCALAR        evalue;
  LIS_INT               emaxiter;
  LIS_REAL          tol;
  LIS_INT               iter,i,j,output;
  LIS_INT               nprocs,my_rank;
  LIS_REAL          nrm2,resid;
  LIS_SCALAR        lshift;
  LIS_VECTOR        r,p,Ax,Ar,Ap;
  LIS_SCALAR        alpha, beta;
  LIS_SCALAR        rAp, rp, ApAp, pAp, pp, ArAp, pAr;
  double      times,itimes,ptimes,p_c_times,p_i_times;

  A = esolver->A;
  x = esolver->x;

  emaxiter = esolver->options[LIS_EOPTIONS_MAXITER];
  tol = esolver->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN]; 
  output  = esolver->options[LIS_EOPTIONS_OUTPUT];
  lshift = esolver->lshift;

#ifdef _LONG__DOUBLE
  if( A->my_rank==0 ) printf("local shift = %Le\n", lshift);
#else
  if( A->my_rank==0 ) printf("local shift = %e\n", lshift);
#endif
  if (lshift != 0) lis_matrix_shift_diagonal(A, lshift);

  r = esolver->work[0];
  p = esolver->work[1];
  Ax = esolver->work[2];
  Ar = esolver->work[3];
  Ap = esolver->work[4];

  /* lambda = <Ax, x> / <x, x> */
  lis_vector_set_all(1.0,x);
  lis_vector_nrm2(x,&nrm2);
  lis_vector_scale(1/nrm2,x);
  lis_matvec(A,x,Ax);
  lis_vector_dot(x,Ax,&evalue);

  /* r = lambda * x - A * x */
  lis_vector_axpyz(-evalue,x,Ax,r); 
  lis_vector_scale(-1.0,r);

  lis_vector_copy(r,p);

  lis_matvec(A,p,Ap);

  iter=0;

  while (iter<emaxiter)
    {
      iter = iter + 1;

      /* alpha = (<r, Ap> - lambda * <r, p>) / (<Ap, Ap> - 2 * lambda * <Ap, p>) */
      lis_vector_dot(r,Ap,&rAp);
      lis_vector_dot(r,p,&rp);
      lis_vector_dot(Ap,Ap,&ApAp);
      lis_vector_dot(p,Ap,&pAp); 
      lis_vector_dot(p,p,&pp);
      alpha = (rAp - evalue * rp) / (ApAp - evalue * (2.0 * pAp - evalue * pp));

      /* x = x + alpha * p */
      lis_vector_axpy(alpha,p,x);

      /* lambda = <Ax, x> / <x, x> */
      lis_matvec(A,x,Ax);
      lis_vector_dot(x,Ax,&evalue);
      lis_vector_nrm2(x, &nrm2);
      evalue = evalue / (nrm2 * nrm2);

      /* r = lambda * x - A * x */
      lis_vector_axpyz(-evalue,x,Ax,r); 
      lis_vector_scale(-1.0,r);
      lis_matvec(A,r,Ar);

      /* beta = - [(<Ar, Ap> - lambda * {<Ar, p> + <Ap. r>} + lambda^2 * <r, p>] / [<Ap, Ap> - 2 * lambda * <Ap, p>] */
      lis_vector_dot(Ar,Ap,&ArAp);
      lis_vector_dot(p,Ar,&pAr);      
      lis_vector_dot(r,Ap,&rAp);      
      lis_vector_dot(r,p,&rp);
      beta = - (ArAp - evalue * ((pAr + rAp) - evalue * rp))/ (ApAp - evalue * (2.0 * pAp - evalue * pp));

      /* p = r + beta * p */
      lis_vector_xpay(r,beta,p);

      /* convergence check */
      lis_vector_nrm2(r,&nrm2);
      resid = fabs(nrm2 / (evalue));

      if( output )
  {
    if( output & LIS_EPRINT_MEM ) esolver->residual[iter] = resid;
    if( output & LIS_EPRINT_OUT && A->my_rank==0 ) lis_print_rhistory(iter,resid);
  }

      if (resid<tol) break;  
      
    }

  esolver->iter       = iter;
  esolver->resid      = resid;
  esolver->evalue[0] = evalue;

  if (lshift != 0) lis_matrix_shift_diagonal(A, -lshift);

  if (resid<tol) 
    {
      esolver->retcode = LIS_SUCCESS;
      return LIS_SUCCESS;
    }
  else
    {
      esolver->retcode = LIS_MAXITER;
      return LIS_MAXITER;
    }
}





