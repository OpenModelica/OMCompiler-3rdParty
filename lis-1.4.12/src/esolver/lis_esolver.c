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
#include <string.h>
#include <stdarg.h>
#include <ctype.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
#ifdef USE_MPI
  #include <mpi.h>
#endif
#include "lislib.h"

/************************************************
 * lis_esolver_init
 * lis_esolver_create
 * lis_esolver_destroy
 * lis_iesolver_destroy
 * lis_esolver_work_destroy
 * lis_esolver_set_option
 * lis_esolver_get_option
 * lis_esolve
 ************************************************/

LIS_ESOLVER_EXECUTE lis_esolver_execute[] = {
  NULL,
  lis_epi, lis_eii, lis_eaii, lis_erqi, lis_esi, lis_eli, lis_ecg, lis_ecr
};

#ifdef USE_QUAD_PRECISION
LIS_ESOLVER_EXECUTE lis_esolver_execute_quad[] = {
  NULL,
  lis_epi_quad, lis_eii_quad, NULL, lis_erqi_quad, lis_esi_quad, lis_eli_quad, NULL, NULL
};
LIS_ESOLVER_EXECUTE lis_esolver_execute_switch[] = {
  NULL,
  NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL
};
#endif

LIS_ESOLVER_CHECK_PARAMS lis_esolver_check_params[] = {
  NULL,
  lis_epi_check_params,  lis_eii_check_params, lis_eaii_check_params, 
  lis_erqi_check_params, lis_esi_check_params, lis_eli_check_params,
  lis_ecg_check_params,  lis_ecr_check_params
};

LIS_ESOLVER_MALLOC_WORK lis_esolver_malloc_work[] = {
  NULL,
  lis_epi_malloc_work,  lis_eii_malloc_work, lis_eaii_malloc_work, 
  lis_erqi_malloc_work, lis_esi_malloc_work, lis_eli_malloc_work, 
  lis_ecg_malloc_work,  lis_ecr_malloc_work
};

#define LIS_ESOLVER_OPTION_LEN  12
#define LIS_ESOLVERS_LEN 8
#define LIS_EPRINT_LEN 4
#define LIS_TRUEFALSE_LEN 2
#define LIS_ESTORAGE_LEN  11
#define LIS_PRECISION_LEN 3

char *LIS_ESOLVER_OPTNAME[] = {
  "-emaxiter", "-etol", "-e", "-ss", "-m", "-shift", "-eprint", "-initx_ones", "-ie", "-estorage", "-estorage_block", "-ef"
};

LIS_INT LIS_ESOLVER_OPTACT[] = {
LIS_EOPTIONS_MAXITER, LIS_EPARAMS_RESID, LIS_EOPTIONS_ESOLVER, 
LIS_EOPTIONS_SUBSPACE, LIS_EOPTIONS_MODE, LIS_EPARAMS_SHIFT, 
LIS_EOPTIONS_OUTPUT, LIS_EOPTIONS_INITGUESS_ONES, 
LIS_EOPTIONS_INNER_ESOLVER, LIS_EOPTIONS_STORAGE, 
LIS_EOPTIONS_STORAGE_BLOCK, LIS_EOPTIONS_PRECISION
};

char *lis_esolver_atoi[] = {"pi", "ii", "aii", "rqi", "si", "li", "cg", "cr"};
char *lis_eprint_atoi[] = {"none", "mem", "out", "all"};
char *lis_etruefalse_atoi[] = {"false", "true"};
char *lis_estorage_atoi[] = {"csr", "csc", "msr", "dia", "ell", "jad", "bsr", "bsc", "vbr", "coo", "dns"};
char *lis_eprecision_atoi[] = {"double", "quad", "switch"};

char *lis_esolvername[] = {"", "Power", "Inverse", "Approximate Inverse", "Rayleigh Quotient", "Subspace", "Lanczos", "CG", "CR"};

char *lis_estoragename[]   = {"CSR", "CSC", "MSR", "DIA", "ELL", "JAD", "BSR", "BSC", "VBR", "COO", "DNS"};

char *lis_ereturncode[] = {"LIS_SUCCESS", "LIS_ILL_OPTION", "LIS_BREAKDOWN", "LIS_OUT_OF_MEMORY", "LIS_MAXITER", "LIS_NOT_IMPLEMENTED", "LIS_ERR_FILE_IO"};

char *lis_eprecisionname[] = {"double", "quad", "switch"};

LIS_VECTOR lis_esolver_evalue = NULL;
LIS_INT lis_esolver_evalue_count = 0;

LIS_VECTOR lis_esolver_evector = NULL;
LIS_INT lis_esolver_evector_count = 0;

LIS_VECTOR lis_esolver_residual_history = NULL;
LIS_INT lis_esolver_residual_history_count = 0;

#undef __FUNC__
#define __FUNC__ "lis_esolver_init"
LIS_INT lis_esolver_init(LIS_ESOLVER esolver)
{

  LIS_DEBUG_FUNC_IN;

  esolver->A        = NULL;
  esolver->x        = NULL;
  esolver->residual = NULL;
  esolver->evalue   = NULL;
  esolver->evector  = NULL;

  esolver->worklen   = 0;
  esolver->iter      = 0;
  esolver->iter2     = 0;
  esolver->resid     = 0;
  esolver->times     = 0;
  esolver->itimes    = 0;
  esolver->ptimes    = 0;
  esolver->p_i_times    = 0;
  esolver->p_c_times    = 0;

  esolver->lshift     = 0;
  esolver->tol       = 0;
  esolver->eprecision = LIS_PRECISION_DOUBLE;

  esolver->options[LIS_EOPTIONS_ESOLVER]               = LIS_ESOLVER_PI;
  esolver->options[LIS_EOPTIONS_MAXITER]               = 1000;
  esolver->options[LIS_EOPTIONS_SUBSPACE]              = 2;
  esolver->options[LIS_EOPTIONS_MODE]                  = 0;
  esolver->options[LIS_EOPTIONS_OUTPUT]                = LIS_FALSE;
  esolver->options[LIS_EOPTIONS_INITGUESS_ONES]        = LIS_TRUE;
  esolver->options[LIS_EOPTIONS_INNER_ESOLVER]         = 2;
  esolver->options[LIS_EOPTIONS_STORAGE]               = 0;
  esolver->options[LIS_EOPTIONS_STORAGE_BLOCK]         = 2;
  esolver->options[LIS_EOPTIONS_PRECISION]             = LIS_PRECISION_DOUBLE;
  esolver->params[LIS_EPARAMS_RESID - LIS_EOPTIONS_LEN] = 1.0e-12;
  esolver->params[LIS_EPARAMS_SHIFT - LIS_EOPTIONS_LEN] = 0.0;

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_create"
LIS_INT lis_esolver_create(LIS_ESOLVER *esolver)
{
  LIS_DEBUG_FUNC_IN;

  *esolver = NULL;

  *esolver = (LIS_ESOLVER)lis_malloc( sizeof(struct LIS_ESOLVER_STRUCT),"lis_esolver_create::esolver" );
  if( NULL==*esolver )
  {
    LIS_SETERR_MEM(sizeof(struct LIS_ESOLVER_STRUCT));
    return LIS_OUT_OF_MEMORY;
  }
  lis_esolver_init(*esolver);

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_work_destroy"
LIS_INT lis_esolver_work_destroy(LIS_ESOLVER esolver)
{
  LIS_INT i;

  LIS_DEBUG_FUNC_IN;

  if( esolver && esolver->work )
  {
    for(i=0;i<esolver->worklen;i++) lis_vector_destroy(esolver->work[i]);
    lis_free(esolver->work);
    esolver->work    = NULL;
    esolver->worklen = 0;
  }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_destroy"
LIS_INT lis_esolver_destroy(LIS_ESOLVER esolver)
{
  LIS_INT i,ss;
  
  LIS_DEBUG_FUNC_IN;

  if( esolver )
  {
    lis_esolver_work_destroy(esolver);
             if( esolver->residual ) lis_free(esolver->residual);
          if( esolver->evalue ) lis_free(esolver->evalue);
             if( esolver->evector ) 
      {
        if ( esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_LI || esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_SI )
          {
      ss=esolver->options[LIS_EOPTIONS_SUBSPACE];
      for(i=0;i<ss+2;i++) lis_vector_destroy(esolver->evector[i]);
          }
        lis_free(esolver->evector);
      }
    lis_free(esolver);
  }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_iesolver_destroy"
LIS_INT lis_iesolver_destroy(LIS_ESOLVER esolver)
{
  LIS_INT i,ss;
  
  LIS_DEBUG_FUNC_IN;

  if( esolver )
  {
    lis_esolver_work_destroy(esolver);
          if( esolver->evalue ) lis_free(esolver->evalue);
             if( esolver->evector ) lis_free(esolver->evector);
    lis_free(esolver);
  }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolve"
LIS_INT lis_esolve(LIS_MATRIX A, LIS_VECTOR x, LIS_SCALAR *evalue0, LIS_ESOLVER esolver)
{
        LIS_INT      nesolver,niesolver,emaxiter; 
  LIS_SCALAR              *evalue;
  LIS_VECTOR          *evector;
  LIS_SCALAR          *residual;
  LIS_INT      err;
  LIS_INT                 output;
  LIS_INT                 ss, mode;
  double                  times;
  double                  gshift;
  LIS_INT      estorage,eblock;
  LIS_MATRIX          B;
  LIS_INT                 eprecision;
  LIS_VECTOR              xx;
  char                    buf[64];
  LIS_INT                 i;
  LIS_SOLVER              solver;

  LIS_DEBUG_FUNC_IN;

  /* begin parameter check */
  err = lis_matrix_check(A,LIS_MATRIX_CHECK_ALL);

  if( err ) return err;
  if( x==NULL )
  {
    LIS_SETERR(LIS_ERR_ILL_ARG,"vector x is undefined\n");
    return LIS_ERR_ILL_ARG;
  }
  if( A->n!=x->n )
  {
    return LIS_ERR_ILL_ARG;
  }
  if( A->gn<=0 )
  {
    LIS_SETERR1(LIS_ERR_ILL_ARG,"Size n(=%d) of matrix A is less than 0\n",A->gn);
    return LIS_ERR_ILL_ARG;
  }

  nesolver = esolver->options[LIS_EOPTIONS_ESOLVER];
  niesolver = esolver->options[LIS_EOPTIONS_INNER_ESOLVER];
  ss = esolver->options[LIS_EOPTIONS_SUBSPACE];
  mode = esolver->options[LIS_EOPTIONS_MODE];
  emaxiter = esolver->options[LIS_EOPTIONS_MAXITER];
  gshift = esolver->params[LIS_EPARAMS_SHIFT - LIS_EOPTIONS_LEN];
  output = esolver->options[LIS_EOPTIONS_OUTPUT];
  estorage = esolver->options[LIS_EOPTIONS_STORAGE];
  eblock = esolver->options[LIS_EOPTIONS_STORAGE_BLOCK];
  eprecision = esolver->options[LIS_EOPTIONS_PRECISION];
  esolver->eprecision = eprecision;

  if ( output ) if( A->my_rank==0 ) printf("shift = %e\n", gshift);

  if( nesolver < 1 || nesolver > LIS_ESOLVERS_LEN )
  {
    LIS_SETERR2(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_ESOLVER is %d (Set between 1 to %d)\n",nesolver, LIS_ESOLVERS_LEN);
    return LIS_ERR_ILL_ARG;
  }

  if( niesolver < 2 || niesolver > 4) 
  {
    LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_INNER_ESOLVER is %d (Set between 2 to 4)\n", niesolver);
    return LIS_ERR_ILL_ARG;
  }

  if ( esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_LI && niesolver == LIS_ESOLVER_PI )
  {
    LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_INNER_ESOLVER is %d (Set between 2 to 4 for Lanczos)\n", niesolver);
    return LIS_ERR_ILL_ARG;
  }

  if ( esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_SI && ss > A->gn )
  {
    LIS_SETERR2(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_SUBSPACE is %d (Set less than or equal to matrix size %d for Subspace)\n", ss, A->gn);
    return LIS_ERR_ILL_ARG;
  }

  if ( esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_LI && ss > A->gn )
  {
    LIS_SETERR2(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_SUBSPACE is %d (Set less than or equal to matrix size %d for Lanczos)\n", ss, A->gn);
    return LIS_ERR_ILL_ARG;
  }

  if ( esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_SI && mode >= ss )
  {
    LIS_SETERR2(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_MODE is %d (Set less than subspace size %d for Subspace)\n", mode, ss);
    return LIS_ERR_ILL_ARG;
  }

  if ( esolver->options[LIS_EOPTIONS_ESOLVER] == LIS_ESOLVER_LI && mode >= ss )
  {
    LIS_SETERR2(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_MODE is %d (Set less than subspace size %d for Lanczos)\n", mode, ss);
    return LIS_ERR_ILL_ARG;
  }

  #ifdef USE_QUAD_PRECISION
    if( eprecision==LIS_PRECISION_QUAD && lis_esolver_execute_quad[nesolver]==NULL )
    {
      LIS_SETERR1(LIS_ERR_NOT_IMPLEMENTED,"Quad precision eigensolver %s is not implemented\n",lis_esolvername[nesolver]);
      return LIS_ERR_NOT_IMPLEMENTED;
    }
    else if( eprecision==LIS_PRECISION_SWITCH && lis_esolver_execute_switch[nesolver]==NULL )
    {
      LIS_SETERR1(LIS_ERR_NOT_IMPLEMENTED,"Switch esolver %s is not implemented\n",lis_esolvername[nesolver]);
      return LIS_ERR_NOT_IMPLEMENTED;
    }
    if( esolver->options[LIS_EOPTIONS_SWITCH_MAXITER]==-1 )
    {
      esolver->options[LIS_EOPTIONS_SWITCH_MAXITER] = emaxiter;
    }
  #endif

  /* create eigenvalue array */
  if( esolver->evalue ) lis_free(esolver->evalue);
  evalue = (LIS_SCALAR *)lis_malloc((ss+2)*sizeof(LIS_SCALAR),"lis_esolve::evalue");
  if( evalue==NULL )
  {
    LIS_SETERR_MEM((ss+2)*sizeof(LIS_SCALAR));
    esolver->retcode = err;
    return err;
  }
  evalue[0] = 1.0;
  evalue[ss-1] = 1.0;


  /* create initial vector */
  #ifndef USE_QUAD_PRECISION
    err = lis_vector_duplicate(A,&xx);
  #else
    if( eprecision==LIS_PRECISION_DOUBLE )
    {
      err = lis_vector_duplicate(A,&xx);
    }
    else
    {
      err = lis_vector_duplicateex(LIS_PRECISION_QUAD,A,&xx);
    }
  #endif
  if( err )
  {
    esolver->retcode = err;
    return err;
  }
  if( esolver->options[LIS_EOPTIONS_INITGUESS_ONES] )
  {
    if( output ) lis_printf(A->comm,"initial vector x = 1\n");
    #ifndef USE_QUAD_PRECISION
      lis_vector_set_all(1.0,xx);
    #else
      if( eprecision==LIS_PRECISION_DOUBLE )
      {
        lis_vector_set_all(1.0,xx);
      }
      else
      {
        lis_vector_set_allex_nm(1.0,xx);
      }
    #endif
  }
  else
  {
    if( output ) lis_printf(A->comm,"initial vector x = user defined\n"); 
    #ifndef USE_QUAD_PRECISION
      lis_vector_copy(x,xx);
    #else
      if( eprecision==LIS_PRECISION_DOUBLE )
      {
        lis_vector_copy(x,xx);
      }
      else
      {
        lis_vector_copyex_nm(x,xx);
      }
    #endif
  }


  /* create eigenvector array */
  if( esolver->evector ) lis_free(esolver->evector);
  evector = (LIS_VECTOR *)lis_malloc((ss+2)*sizeof(LIS_VECTOR),"lis_esolve::evector");
  if( evector==NULL )
  {
    LIS_SETERR_MEM((ss+2)*sizeof(LIS_VECTOR));
    esolver->retcode = err;
    return err;
  }

  /* create residual history vector */
  if( esolver->residual ) lis_free(esolver->residual);
  residual = (LIS_SCALAR *)lis_malloc((emaxiter+2)*sizeof(LIS_SCALAR),"lis_esolve::residual");
  if( residual==NULL )
  {
    LIS_SETERR_MEM((emaxiter+2)*sizeof(LIS_SCALAR));
    lis_vector_destroy(xx);
    esolver->retcode = err;
    return err;
  }

  /* Matrix Convert */
  if( estorage>0 && A->matrix_type!=estorage )
  {
    err = lis_matrix_duplicate(A,&B);
    if( err ) return err;
    lis_matrix_set_blocksize(B,eblock,eblock,NULL,NULL);
    lis_matrix_set_type(B,estorage);
    err = lis_matrix_convert(A,B);
    if( err ) return err;
    lis_matrix_storage_destroy(A);
    lis_matrix_DLU_destroy(A);
    lis_matrix_diag_destroy(A->WD);
    if( A->l2g_map ) lis_free( A->l2g_map );
    if( A->commtable ) lis_commtable_destroy( A->commtable );
    if( A->ranges ) lis_free( A->ranges );
    err = lis_matrix_copy_struct(B,A);
    if( err ) return err;
    lis_free(B);
  }

  esolver->A        = A;
  esolver->evalue   = evalue;
  esolver->x        = x;
  esolver->evector  = evector;
  residual[0] = 1.0;
  esolver->residual = residual;

        if( A->my_rank==0 )
    {
#ifdef _LONG__DOUBLE
        if( output ) printf("precision  : long double\n"); 
#else
      if( output ) printf("precision  : %s\n", lis_eprecisionname[eprecision]);
#endif
#ifdef _LONGLONG
      if ( output ) printf("esolver    : %s %lld\n", lis_esolvername[nesolver],nesolver);
#else
      if ( output ) printf("esolver    : %s %d\n", lis_esolvername[nesolver],nesolver);
#endif
    }

  if( A->my_rank==0 )
    {
      if( A->matrix_type==LIS_MATRIX_BSR || A->matrix_type==LIS_MATRIX_BSC )
        {
#ifdef _LONGLONG
    if( output ) printf("storage    : %s(%lld x %lld)\n", lis_estoragename[A->matrix_type-1],eblock,eblock); 
#else
    if( output ) printf("storage    : %s(%d x %d)\n", lis_estoragename[A->matrix_type-1],eblock,eblock); 
#endif
        }
      else
        {
    if( output ) printf("storage    : %s\n", lis_estoragename[A->matrix_type-1]); 
        }
    }
  
  times = lis_wtime();

  esolver->ptimes = 0;
  esolver->itimes = 0;
  esolver->p_c_times = 0;
  esolver->p_i_times = 0;


  if (gshift != 0.0) lis_matrix_shift_diagonal(A, gshift);

  /* create work vector */
  err = lis_esolver_malloc_work[nesolver](esolver);
  if( err )
  {
    lis_vector_destroy(xx);
    esolver->retcode = err;
    return err;
  }

  esolver->x        = xx;
  esolver->xx       = x;

  /* execute esolver */

  #ifndef USE_QUAD_PRECISION
    err = lis_esolver_execute[nesolver](esolver);
  #else
    if( eprecision==LIS_PRECISION_DOUBLE )
    {
      err = lis_esolver_execute[nesolver](esolver);
    }
    else if( eprecision==LIS_PRECISION_QUAD )
    {
      err = lis_esolver_execute_quad[nesolver](esolver);
    }
    else if( eprecision==LIS_PRECISION_SWITCH )
    {
      err = lis_esolver_execute_switch[nesolver](esolver);
    }
  #endif
  esolver->retcode = err;

  *evalue0 = esolver->evalue[mode];
  lis_vector_copy(esolver->x, x);

  esolver->times = lis_wtime() - times; 

  lis_matrix_shift_diagonal(A, -gshift);

        if( A->my_rank==0 )
        {
                if( err )
                {
#ifdef _LONGLONG
                  if ( output ) printf("lis_esolve : %s(code=%lld)\n\n",lis_ereturncode[err],err);
#else
                  if ( output ) printf("lis_esolve : %s(code=%d)\n\n",lis_ereturncode[err],err);
#endif

                }
                else
                {
                  if ( output ) printf("lis_esolve : normal end\n\n");
                }
        }

  if( eprecision==LIS_PRECISION_DOUBLE )
  {
    esolver->iter2 = esolver->iter;
  }
  else if( eprecision==LIS_PRECISION_QUAD )
  {
    esolver->iter2 = 0;
  }

  lis_vector_destroy(xx);

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_optionC"
LIS_INT lis_esolver_set_optionC(LIS_ESOLVER esolver)
{
  LIS_ARGS  p;

  LIS_DEBUG_FUNC_IN;

  p = cmd_args->next;
  while( p!=cmd_args )
  {
    lis_esolver_set_option2(p->arg1,p->arg2,esolver);
    p = p->next;
  }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_option"
LIS_INT lis_esolver_set_option(char *text, LIS_ESOLVER esolver)
{
  LIS_ARGS  args,p;

  LIS_DEBUG_FUNC_IN;

  lis_text2args(text,&args);
  p = args->next;
  while( p!=args )
  {
    lis_esolver_set_option2(p->arg1,p->arg2,esolver);
    p = p->next;
  }
  lis_args_free(args);

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_option2"
LIS_INT lis_esolver_set_option2(char* arg1, char *arg2, LIS_ESOLVER esolver)
{
  LIS_INT i;

  LIS_DEBUG_FUNC_IN;
  
  for(i=0;i<LIS_ESOLVER_OPTION_LEN;i++)
  {
    if( strcmp(arg1, LIS_ESOLVER_OPTNAME[i])==0 )
    {
      switch( LIS_ESOLVER_OPTACT[i] )
      {
      case LIS_EOPTIONS_ESOLVER:
        lis_esolver_set_option_esolver(arg2,esolver);
        break;
                        case LIS_EOPTIONS_OUTPUT:
        lis_esolver_set_option_print(arg2,esolver);
        break;
                        case LIS_EOPTIONS_INITGUESS_ONES:
        lis_esolver_set_option_truefalse(arg2,LIS_EOPTIONS_INITGUESS_ONES,esolver);
        break;
      case LIS_EOPTIONS_INNER_ESOLVER:
        lis_esolver_set_option_iesolver(arg2,esolver);
        break;
      case LIS_EOPTIONS_STORAGE:
        lis_esolver_set_option_storage(arg2,esolver);
        break;
                        case LIS_EOPTIONS_PRECISION:
        lis_esolver_set_option_eprecision(arg2,LIS_EOPTIONS_PRECISION,esolver);
        break;
      default:
        if( LIS_ESOLVER_OPTACT[i] < LIS_EOPTIONS_LEN )
          {
#ifdef _LONGLONG
            sscanf(arg2, "%lld", &esolver->options[LIS_ESOLVER_OPTACT[i]]);
#else
            sscanf(arg2, "%d", &esolver->options[LIS_ESOLVER_OPTACT[i]]);
#endif
          }
        else
          {
#ifdef _LONG__DOUBLE
            sscanf(arg2, "%Lg", &esolver->params[LIS_ESOLVER_OPTACT[i]-LIS_EOPTIONS_LEN]);
#else
            sscanf(arg2, "%lg", &esolver->params[LIS_ESOLVER_OPTACT[i]-LIS_EOPTIONS_LEN]);
#endif
          }
        break;
      }
    }
  }
  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_option_esolver"
LIS_INT lis_esolver_set_option_esolver(char *argv, LIS_ESOLVER esolver)
{
  LIS_INT  i;

  LIS_DEBUG_FUNC_IN;

  if( argv[0]>='0' && argv[0]<='9' )
  {
#ifdef _LONGLONG
    sscanf(argv, "%lld", &esolver->options[LIS_EOPTIONS_ESOLVER]);
#else
    sscanf(argv, "%d", &esolver->options[LIS_EOPTIONS_ESOLVER]);
#endif
  }
  else
  {
    for(i=0;i<LIS_ESOLVER_LEN;i++)
    {
      if( strcmp(argv,lis_esolver_atoi[i])==0 )
      {
        esolver->options[LIS_EOPTIONS_ESOLVER] = i+1;
        break;
      }
    }
  }
  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_option_iesolver"
LIS_INT lis_esolver_set_option_iesolver(char *argv, LIS_ESOLVER esolver)
{
  LIS_INT  i;

  LIS_DEBUG_FUNC_IN;

  if( argv[0]>='0' && argv[0]<='9' )
  {
#ifdef _LONGLONG
    sscanf(argv, "%lld", &esolver->options[LIS_EOPTIONS_INNER_ESOLVER]);
#else
    sscanf(argv, "%d", &esolver->options[LIS_EOPTIONS_INNER_ESOLVER]);
#endif
  }
  else
  {
    for(i=0;i<LIS_ESOLVER_LEN;i++)
    {
      if( strcmp(argv,lis_esolver_atoi[i])==0 )
      {
        esolver->options[LIS_EOPTIONS_INNER_ESOLVER] = i+1;
        break;
      }
    }
  }
  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_option_print"
LIS_INT lis_esolver_set_option_print(char *argv, LIS_ESOLVER esolver)
{
  LIS_INT  i;

  LIS_DEBUG_FUNC_IN;

  if( argv[0]>='0' && argv[0]<='3' )
  {
#ifdef _LONGLONG
    sscanf(argv, "%lld", &esolver->options[LIS_EOPTIONS_OUTPUT]);
#else
    sscanf(argv, "%d", &esolver->options[LIS_EOPTIONS_OUTPUT]);
#endif
  }
  else
  {
    for(i=0;i<LIS_EPRINT_LEN;i++)
    {
      if( strcmp(argv,lis_eprint_atoi[i])==0 )
      {
        esolver->options[LIS_EOPTIONS_OUTPUT] = i;
        break;
      }
    }
  }
  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_solver_set_option_truefalse"
LIS_INT lis_esolver_set_option_truefalse(char *argv, LIS_INT opt, LIS_ESOLVER esolver)
{
        LIS_INT  i;

        LIS_DEBUG_FUNC_IN;

        if( argv[0]>='0' && argv[0]<='1' )
        {
#ifdef _LONGLONG
                sscanf(argv, "%lld", &esolver->options[opt]);
#else
                sscanf(argv, "%d", &esolver->options[opt]);
#endif
        }
        else
        {
                for(i=0;i<LIS_TRUEFALSE_LEN;i++)
                {
                        if( strcmp(argv,lis_etruefalse_atoi[i])==0 )
                        {
                                esolver->options[opt] = i;
                                break;
                        }
                }
        }
        LIS_DEBUG_FUNC_OUT;
        return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_option_eprecision"
LIS_INT lis_esolver_set_option_eprecision(char *argv, LIS_INT opt, LIS_ESOLVER esolver)
{
  LIS_INT  i;

  LIS_DEBUG_FUNC_IN;

  if( argv[0]>='0' && argv[0]<='1' )
  {
#ifdef _LONGLONG
    sscanf(argv, "%lld", &esolver->options[opt]);
#else
    sscanf(argv, "%d", &esolver->options[opt]);
#endif
  }
  else
  {
    for(i=0;i<LIS_PRECISION_LEN;i++)
    {
      if( strcmp(argv,lis_eprecision_atoi[i])==0 )
      {
        esolver->options[opt] = i;
        break;
      }
    }
  }
  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_set_option_storage"
LIS_INT lis_esolver_set_option_storage(char *argv, LIS_ESOLVER esolver)
{
  LIS_INT  i;

  LIS_DEBUG_FUNC_IN;

  if( argv[0]>='0' && argv[0]<='9' )
  {
#ifdef _LONGLONG
    sscanf(argv, "%lld", &esolver->options[LIS_EOPTIONS_STORAGE]);
#else
    sscanf(argv, "%d", &esolver->options[LIS_EOPTIONS_STORAGE]);
#endif
  }
  else
  {
    for(i=0;i<LIS_ESTORAGE_LEN;i++)
    {
      if( strcmp(argv,lis_estorage_atoi[i])==0 )
      {
        esolver->options[LIS_EOPTIONS_STORAGE] = i+1;
        break;
      }
    }
  }
  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_iters"
LIS_INT lis_esolver_get_iters(LIS_ESOLVER esolver, LIS_INT *iters)
{
  LIS_DEBUG_FUNC_IN;

  *iters = esolver->iter;

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_itersex"
LIS_INT lis_esolver_get_itersex(LIS_ESOLVER esolver, LIS_INT *iters, LIS_INT *iters_double, LIS_INT *iters_quad)
{
  LIS_DEBUG_FUNC_IN;

  *iters = esolver->iter;
  *iters_double = esolver->iter2;
  *iters_quad = esolver->iter - esolver->iter2;

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_time"
LIS_INT lis_esolver_get_time(LIS_ESOLVER esolver, double *times)
{
  LIS_DEBUG_FUNC_IN;

  *times  = esolver->times;

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_timeex"
LIS_INT lis_esolver_get_timeex(LIS_ESOLVER esolver, double *times, double *itimes, double *ptimes, double *p_c_times, double *p_i_times)
{
  LIS_DEBUG_FUNC_IN;

  *times  = esolver->times;
  if( itimes ) *itimes = esolver->itimes;
  if( ptimes ) *ptimes = esolver->ptimes;
  if( p_c_times ) *p_c_times = esolver->p_c_times;
  if( p_i_times ) *p_i_times = esolver->p_i_times;

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_residualnorm"
LIS_INT lis_esolver_get_residualnorm(LIS_ESOLVER esolver, LIS_REAL *residual)
{
  LIS_DEBUG_FUNC_IN;

  *residual  = esolver->resid;

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_esolver"
LIS_INT lis_esolver_get_esolver(LIS_ESOLVER esolver, LIS_INT *nesol)
{
  LIS_DEBUG_FUNC_IN;

  *nesol = esolver->options[LIS_EOPTIONS_ESOLVER];

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_evalues"
LIS_INT lis_esolver_get_evalues(LIS_ESOLVER esolver, LIS_VECTOR v)
{
  LIS_INT             i,ii,ss,n,gn,is,ie;

  LIS_DEBUG_FUNC_IN;

  if ( esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_SI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_LI )
  {
    LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_ESOLVER is %d (Set Subspace or Lanczos)\n", esolver->options[LIS_EOPTIONS_ESOLVER]);
    return LIS_ERR_ILL_ARG;
  }

  ss = esolver->options[LIS_EOPTIONS_SUBSPACE];
  lis_vector_set_size(v,0,ss);
  lis_vector_get_size(v,&n,&gn);
  lis_vector_get_range(v,&is,&ie);
  
  for(i=0;i<n;i++)
  {
    ii=i;
    if( v->origin ) ii++;
    lis_vector_set_value(LIS_INS_VALUE,ii+is,esolver->evalue[i+is],v);
  }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_evectors"
LIS_INT lis_esolver_get_evectors(LIS_ESOLVER esolver, LIS_MATRIX A)
{
  LIS_INT             i,ii,j,jj,n,gn,is,ie,js;
  LIS_INT             ss,lis_esolver_evector_size;

  LIS_DEBUG_FUNC_IN;

  if ( esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_SI && esolver->options[LIS_EOPTIONS_ESOLVER] != LIS_ESOLVER_LI )
  {
    LIS_SETERR1(LIS_ERR_ILL_ARG,"Parameter LIS_EOPTIONS_ESOLVER is %d (Set Subspace or Lanczos)\n", esolver->options[LIS_EOPTIONS_ESOLVER]);
    return LIS_ERR_ILL_ARG;
  }

  ss = esolver->options[LIS_EOPTIONS_SUBSPACE];
  lis_esolver_evector_size = esolver->evector[0]->gn; 
  lis_matrix_set_size(A,0,lis_esolver_evector_size);
  lis_matrix_get_size(A,&n,&gn);
  lis_matrix_get_range(A,&is,&ie);
  js=0;
  if( esolver->evector[0]->origin ) 
    {
      is++;
      js++;
    }
  for(j=0;j<ss;j++)
    {
      for(i=0;i<n;i++)
        {
    ii=i+is;
    jj=j+js;
    lis_matrix_set_value(LIS_INS_VALUE,ii,jj,esolver->evector[j]->value[i],A);
        }
    }
  lis_matrix_set_type(A,LIS_MATRIX_CSR);
  lis_matrix_assemble(A);

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_status"
LIS_INT lis_esolver_get_status(LIS_ESOLVER esolver, LIS_INT *status)
{
  LIS_DEBUG_FUNC_IN;

  *status = esolver->retcode;

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_rhistory"
LIS_INT lis_esolver_get_rhistory(LIS_ESOLVER esolver, LIS_VECTOR v)
{
  LIS_INT    i,n,maxiter;

  maxiter = esolver->iter+1;
        if( esolver->retcode!=LIS_SUCCESS )
    {
      maxiter--;
    }
  n = _min(v->n,maxiter);
  for(i=0;i<n;i++)
  {
    v->value[i] = esolver->residual[i];
  }
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_esolver_get_esolvername"
LIS_INT lis_esolver_get_esolvername(LIS_INT esolver, char *esolvername)
{
  LIS_DEBUG_FUNC_IN;

  if( esolver < 1 || esolver > LIS_ESOLVERS_LEN )
  {
    esolvername = NULL;
    return LIS_FAILS;
  }
  strcpy(esolvername,lis_esolvername[esolver]);

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

