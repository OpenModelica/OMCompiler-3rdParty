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
#include <math.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
#ifdef USE_MPI
  #include <mpi.h>
#endif
#include "lislib.h"

/************************************************
 * lis_matrix_set_blocksize
 * lis_matrix_convert
 * lis_matrix_convert_self
 * lis_matrix_get_diagonal
 * lis_matrix_scaling
 * lis_matrix_split
 * lis_matrix_merge
 * lis_matrix_copy
 * lis_matrix_copyDLU
 * lis_matrix_solve
 * lis_matrix_solvet
 * lis_array_LUdecomp
 * lis_array_invGauss
 * lis_array_matvec
 * lis_array_matvect
 * lis_array_matmat
 * lis_array_matmat2
 * lis_array_nrm2
 * lis_array_nrm1
 * lis_array_dot
 * lis_array_matinv
 * lis_array_invvec
 * lis_array_invvect
 * lis_array_matvec2
 * lis_array_solve
 * lis_array_cgs
 * lis_array_mgs
 * lis_array_qr
 * lis_array_power
 * lis_array_set_all
 * lis_array_scale
 * lis_array_dot2
 * lis_array_axpyz
 * lis_array_copy
 * lis_matrix_shift_diagonal
 ************************************************/

#undef __FUNC__
#define __FUNC__ "lis_matrix_set_blocksize"
LIS_INT lis_matrix_set_blocksize(LIS_MATRIX A, LIS_INT bnr, LIS_INT bnc, LIS_INT row[], LIS_INT col[])
{
  LIS_INT i,n;
  LIS_INT err;
  LIS_INT *conv_row,*conv_col;

  LIS_DEBUG_FUNC_IN;

  err = lis_matrix_check(A,LIS_MATRIX_CHECK_NULL);
  if( err ) return err;

  if( bnr<=0 || bnc<=0 )
  {
    LIS_SETERR2(LIS_ERR_ILL_ARG,"bnr=%d <= 0 or bnc=%d <= 0\n",bnr,bnc);
    return LIS_ERR_ILL_ARG;
  }
  if( (row==NULL && col!=NULL) || (row!=NULL && col==NULL) )
  {
    LIS_SETERR2(LIS_ERR_ILL_ARG,"either row[]=%x or col[]=%x is NULL\n",row,col);
    return LIS_ERR_ILL_ARG;
  }
  if( row==NULL )
  {
    A->conv_bnr = bnr;
    A->conv_bnc = bnc;
  }
  else
  {
    n = A->n;
    conv_row = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_matrix_set_blocksize::conv_row");
    if( conv_row==NULL )
    {
      LIS_SETERR_MEM(n*sizeof(LIS_INT));
      return LIS_OUT_OF_MEMORY;
    }
    conv_col = (LIS_INT *)lis_malloc(n*sizeof(LIS_INT),"lis_matrix_set_blocksize::conv_col");
    if( conv_col==NULL )
    {
      LIS_SETERR_MEM(n*sizeof(LIS_INT));
      lis_free(conv_row);
      return LIS_OUT_OF_MEMORY;
    }
    for(i=0;i<n;i++)
    {
      conv_row[i] = row[i];
      conv_col[i] = col[i];
    }
    A->conv_row = conv_row;
    A->conv_col = conv_col;
  }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_convert"
LIS_INT lis_matrix_convert(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
  LIS_INT err;
  LIS_INT istmp;
  LIS_INT convert_matrix_type;
  LIS_MATRIX Atmp,Atmp2;

  LIS_DEBUG_FUNC_IN;

  err = lis_matrix_check(Ain,LIS_MATRIX_CHECK_ALL);
  if( err ) return err;
  err = lis_matrix_check(Aout,LIS_MATRIX_CHECK_NULL);
  if( err ) return err;

  err = lis_matrix_merge(Ain);
  if( err ) return err;

  convert_matrix_type = Aout->matrix_type;

  if( Ain->matrix_type==convert_matrix_type && !Ain->is_block )
  {
    err = lis_matrix_copy(Ain,Aout);
    return err;
  }
  if( Ain->matrix_type!=LIS_MATRIX_CSR )
  {
    istmp = LIS_TRUE;
    switch( Ain->matrix_type )
    {
    case LIS_MATRIX_RCO:
      switch( convert_matrix_type )
      {
      case LIS_MATRIX_CSR:
        err = lis_matrix_convert_rco2csr(Ain,Aout);
        LIS_DEBUG_FUNC_OUT;
        return err;
        break;
      case LIS_MATRIX_BSR:
        err = lis_matrix_convert_rco2bsr(Ain,Aout);
        LIS_DEBUG_FUNC_OUT;
        return err;
      case LIS_MATRIX_CSC:
        err = lis_matrix_convert_rco2csc(Ain,Aout);
        LIS_DEBUG_FUNC_OUT;
        return err;
        break;
      default:
        err = lis_matrix_duplicate(Ain,&Atmp);
        if( err ) return err;
        err = lis_matrix_convert_rco2csr(Ain,Atmp);
        break;
      }
      break;
    case LIS_MATRIX_CSC:
      switch( convert_matrix_type )
      {
      case LIS_MATRIX_BSC:
        err = lis_matrix_convert_csc2bsc(Ain,Aout);
        LIS_DEBUG_FUNC_OUT;
        return err;
      default:
        err = lis_matrix_duplicate(Ain,&Atmp);
        if( err ) return err;
        err = lis_matrix_convert_csc2csr(Ain,Atmp);
        break;
      }
      break;
    case LIS_MATRIX_MSR:
      err = lis_matrix_duplicate(Ain,&Atmp);
      if( err ) return err;
      err = lis_matrix_convert_msr2csr(Ain,Atmp);
      break;
    case LIS_MATRIX_DIA:
      err = lis_matrix_duplicate(Ain,&Atmp);
      if( err ) return err;
      err = lis_matrix_convert_dia2csr(Ain,Atmp);
      break;
    case LIS_MATRIX_ELL:
      err = lis_matrix_duplicate(Ain,&Atmp);
      if( err ) return err;
      err = lis_matrix_convert_ell2csr(Ain,Atmp);
      break;
    case LIS_MATRIX_JAD:
      err = lis_matrix_duplicate(Ain,&Atmp);
      if( err ) return err;
      err = lis_matrix_convert_jad2csr(Ain,Atmp);
      break;
    case LIS_MATRIX_BSR:
      err = lis_matrix_duplicate(Ain,&Atmp);
      if( err ) return err;
      err = lis_matrix_convert_bsr2csr(Ain,Atmp);
      break;
    case LIS_MATRIX_BSC:
      err = lis_matrix_duplicate(Ain,&Atmp);
      if( err ) return err;
      err = lis_matrix_convert_bsc2csr(Ain,Atmp);
      break;
    case LIS_MATRIX_VBR:
      err = lis_matrix_duplicate(Ain,&Atmp);
      if( err ) return err;
      err = lis_matrix_convert_vbr2csr(Ain,Atmp);
      break;
    case LIS_MATRIX_DNS:
      err = lis_matrix_duplicate(Ain,&Atmp);
      if( err ) return err;
      err = lis_matrix_convert_dns2csr(Ain,Atmp);
      break;
    case LIS_MATRIX_COO:
      err = lis_matrix_duplicate(Ain,&Atmp);
      if( err ) return err;
      err = lis_matrix_convert_coo2csr(Ain,Atmp);
      break;
    default:
      LIS_SETERR_IMP;
      err = LIS_ERR_NOT_IMPLEMENTED;
    }
    if( err )
    {
      return err;
    }
    if( convert_matrix_type== LIS_MATRIX_CSR )
    {
      lis_matrix_storage_destroy(Aout);
      lis_matrix_DLU_destroy(Aout);
      lis_matrix_diag_destroy(Aout->WD);
      if( Aout->l2g_map ) lis_free( Aout->l2g_map );
      if( Aout->commtable ) lis_commtable_destroy( Aout->commtable );
      if( Aout->ranges ) lis_free( Aout->ranges );
      lis_matrix_copy_struct(Atmp,Aout);
      lis_free(Atmp);
      return LIS_SUCCESS;
    }
  }
  else
  {
    istmp = LIS_FALSE;
    Atmp = Ain;
  }

  switch( convert_matrix_type )
  {
  case LIS_MATRIX_BSR:
    err = lis_matrix_convert_csr2bsr(Atmp,Aout);
    break;
  case LIS_MATRIX_CSC:
    err = lis_matrix_convert_csr2csc(Atmp,Aout); 
    break;
  case LIS_MATRIX_MSR:
    err = lis_matrix_convert_csr2msr(Atmp,Aout);
    break;
  case LIS_MATRIX_ELL:
    err = lis_matrix_convert_csr2ell(Atmp,Aout);
    break;
  case LIS_MATRIX_DIA:
    err = lis_matrix_convert_csr2dia(Atmp,Aout);
    break;
  case LIS_MATRIX_JAD:
    err = lis_matrix_convert_csr2jad(Atmp,Aout);
    break;
  case LIS_MATRIX_BSC:
    err = lis_matrix_duplicate(Atmp,&Atmp2);
    if( err ) return err;
    err = lis_matrix_convert_csr2csc(Atmp,Atmp2);
    if( err ) return err;
    if( Atmp!=Ain )
    {
      lis_matrix_destroy(Atmp);
    }
    Atmp = Atmp2;
    istmp = LIS_TRUE;
    err = lis_matrix_convert_csc2bsc(Atmp,Aout);
    break;
  case LIS_MATRIX_VBR:
    err = lis_matrix_convert_csr2vbr(Atmp,Aout);
    break;
  case LIS_MATRIX_DNS:
    err = lis_matrix_convert_csr2dns(Atmp,Aout);
    break;
  case LIS_MATRIX_COO:
    err = lis_matrix_convert_csr2coo(Atmp,Aout);
    break;
  default:
    LIS_SETERR_IMP;
    err = LIS_ERR_NOT_IMPLEMENTED;
  }

  if( istmp ) lis_matrix_destroy(Atmp);
  if( err )
  {
    return err;
  }

  LIS_DEBUG_FUNC_OUT;
  return err;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_convert_self"
LIS_INT lis_matrix_convert_self(LIS_SOLVER solver)
{
  LIS_INT      err;
  LIS_INT      storage,block;
  LIS_MATRIX  A,B;

  LIS_DEBUG_FUNC_IN;

  A = solver->A;
  storage     = solver->options[LIS_OPTIONS_STORAGE];
  block       = solver->options[LIS_OPTIONS_STORAGE_BLOCK];

  if( storage>0 && A->matrix_type!=storage )
  {
    err = lis_matrix_duplicate(A,&B);
    if( err ) return err;
    lis_matrix_set_blocksize(B,block,block,NULL,NULL);
    lis_matrix_set_type(B,storage);
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
    if( A->matrix_type==LIS_MATRIX_JAD )
    {
      A->work = (LIS_SCALAR *)lis_malloc(A->n*sizeof(LIS_SCALAR),"lis_precon_create_bjacobi::A->work");
      if( A->work==NULL )
      {
        LIS_SETERR_MEM(A->n*sizeof(LIS_SCALAR));
        return LIS_OUT_OF_MEMORY;
      }
    }
  }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_get_diagonal"
LIS_INT lis_matrix_get_diagonal(LIS_MATRIX A, LIS_VECTOR D)
{
  LIS_SCALAR *d;

  LIS_DEBUG_FUNC_IN;

  d = D->value;
  switch( A->matrix_type )
  {
  case LIS_MATRIX_CSR:
    lis_matrix_get_diagonal_csr(A, d);
    break;
  case LIS_MATRIX_CSC:
    lis_matrix_get_diagonal_csc(A, d);
    break;
  case LIS_MATRIX_MSR:
    lis_matrix_get_diagonal_msr(A, d);
    break;
  case LIS_MATRIX_DIA:
    lis_matrix_get_diagonal_dia(A, d);
    break;
  case LIS_MATRIX_ELL:
    lis_matrix_get_diagonal_ell(A, d);
    break;
  case LIS_MATRIX_JAD:
    lis_matrix_get_diagonal_jad(A, d);
    break;
  case LIS_MATRIX_BSR:
    lis_matrix_get_diagonal_bsr(A, d);
    break;
  case LIS_MATRIX_BSC:
    lis_matrix_get_diagonal_bsc(A, d);
    break;
  case LIS_MATRIX_DNS:
    lis_matrix_get_diagonal_dns(A, d);
    break;
  case LIS_MATRIX_COO:
    lis_matrix_get_diagonal_coo(A, d);
    break;
  case LIS_MATRIX_VBR:
    lis_matrix_get_diagonal_vbr(A, d);
    break;
  default:
    LIS_SETERR_IMP;
    return LIS_ERR_NOT_IMPLEMENTED;
    break;
  }
  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_scaling"
LIS_INT lis_matrix_scaling(LIS_MATRIX A, LIS_VECTOR B, LIS_VECTOR D, LIS_INT action)
{
  LIS_INT      i,n,np;
  LIS_SCALAR  *b,*d;

  n  = A->n;
  np = A->np;
  b  = B->value;
  d  = D->value;

  lis_matrix_get_diagonal(A,D);
  if( action==LIS_SCALE_SYMM_DIAG )
  {
#ifdef USE_MPI
    if( A->np>D->np )
    {
      D->value = (LIS_SCALAR *)lis_realloc(D->value,A->np*sizeof(LIS_SCALAR));
      if( D->value==NULL )
      {
        LIS_SETERR_MEM(A->np*sizeof(LIS_SCALAR));
        return LIS_OUT_OF_MEMORY;
      }
      d = D->value;
    }
    lis_send_recv(A->commtable,d);
#endif
    #ifdef _OPENMP
    #pragma omp parallel for private(i)
    #endif
    for(i=0; i<np; i++)
    {
      d[i] = 1.0 / sqrt(fabs(d[i]));
    }

    switch( A->matrix_type )
    {
    case LIS_MATRIX_CSR:
      lis_matrix_scaling_symm_csr(A, d);
      break;
    case LIS_MATRIX_CSC:
      lis_matrix_scaling_symm_csc(A, d);
      break;
    case LIS_MATRIX_MSR:
      lis_matrix_scaling_symm_msr(A, d);
      break;
    case LIS_MATRIX_DIA:
      lis_matrix_scaling_symm_dia(A, d);
      break;
    case LIS_MATRIX_ELL:
      lis_matrix_scaling_symm_ell(A, d);
      break;
    case LIS_MATRIX_JAD:
      lis_matrix_scaling_symm_jad(A, d);
      break;
    case LIS_MATRIX_BSR:
      lis_matrix_scaling_symm_bsr(A, d);
      break;
    case LIS_MATRIX_BSC:
      lis_matrix_scaling_symm_bsc(A, d);
      break;
    case LIS_MATRIX_DNS:
      lis_matrix_scaling_symm_dns(A, d);
      break;
    case LIS_MATRIX_COO:
      lis_matrix_scaling_symm_coo(A, d);
      break;
    case LIS_MATRIX_VBR:
      lis_matrix_scaling_symm_vbr(A, d);
      break;
    default:
      LIS_SETERR_IMP;
      return LIS_ERR_NOT_IMPLEMENTED;
      break;
    }
  }
  else
  {
    #ifdef _OPENMP
    #pragma omp parallel for private(i)
    #endif
    for(i=0; i<n; i++)
    {
      d[i] = 1.0 / d[i];
    }
    switch( A->matrix_type )
    {
    case LIS_MATRIX_CSR:
      lis_matrix_scaling_csr(A, d);
      break;
    case LIS_MATRIX_CSC:
      lis_matrix_scaling_csc(A, d);
      break;
    case LIS_MATRIX_MSR:
      lis_matrix_scaling_msr(A, d);
      break;
    case LIS_MATRIX_DIA:
      lis_matrix_scaling_dia(A, d);
      break;
    case LIS_MATRIX_ELL:
      lis_matrix_scaling_ell(A, d);
      break;
    case LIS_MATRIX_JAD:
      lis_matrix_scaling_jad(A, d);
      break;
    case LIS_MATRIX_BSR:
      lis_matrix_scaling_bsr(A, d);
      break;
    case LIS_MATRIX_BSC:
      lis_matrix_scaling_bsc(A, d);
      break;
    case LIS_MATRIX_DNS:
      lis_matrix_scaling_dns(A, d);
      break;
    case LIS_MATRIX_COO:
      lis_matrix_scaling_coo(A, d);
      break;
    case LIS_MATRIX_VBR:
      lis_matrix_scaling_vbr(A, d);
      break;
    default:
      LIS_SETERR_IMP;
      return LIS_ERR_NOT_IMPLEMENTED;
      break;
    }
  }

  #ifdef _OPENMP
  #pragma omp parallel for private(i)
  #endif
  for(i=0; i<n; i++)
  {
    b[i] = b[i]*d[i];
  }
  A->is_scaled = LIS_TRUE;
  B->is_scaled = LIS_TRUE;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_split"
LIS_INT lis_matrix_split(LIS_MATRIX A)
{
  LIS_INT err;

  LIS_DEBUG_FUNC_IN;

  if( A->is_splited )
  {
    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
  }
  switch( A->matrix_type )
  {
  case LIS_MATRIX_CSR:
    err = lis_matrix_split_csr(A);
    break;
  case LIS_MATRIX_CSC:
    err = lis_matrix_split_csc(A);
    break;
  case LIS_MATRIX_BSR:
    err = lis_matrix_split_bsr(A);
    break;
  case LIS_MATRIX_MSR:
    err = lis_matrix_split_msr(A);
    break;
  case LIS_MATRIX_ELL:
    err = lis_matrix_split_ell(A);
    break;
  case LIS_MATRIX_DIA:
    err = lis_matrix_split_dia(A);
    break;
  case LIS_MATRIX_JAD:
    err = lis_matrix_split_jad(A);
    break;
  case LIS_MATRIX_BSC:
    err = lis_matrix_split_bsc(A);
    break;
  case LIS_MATRIX_DNS:
    err = lis_matrix_split_dns(A);
    break;
  case LIS_MATRIX_COO:
    err = lis_matrix_split_coo(A);
    break;
  case LIS_MATRIX_VBR:
    err = lis_matrix_split_vbr(A);
    break;
  default:
    LIS_SETERR_IMP;
    return LIS_ERR_NOT_IMPLEMENTED;
    break;
  }

  if( err ) return err;
  /*
         if( !A->is_save )
         {
           lis_matrix_storage_destroy(A);
         }
  */
  A->is_splited = LIS_TRUE;

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_merge"
LIS_INT lis_matrix_merge(LIS_MATRIX A)
{
  LIS_INT err;

  LIS_DEBUG_FUNC_IN;

  if( !A->is_splited || (A->is_save && A->is_splited) )
  {
    LIS_DEBUG_FUNC_OUT;
    return LIS_SUCCESS;
  }
  switch( A->matrix_type )
  {
  case LIS_MATRIX_CSR:
    err = lis_matrix_merge_csr(A);
    break;
  case LIS_MATRIX_CSC:
    err = lis_matrix_merge_csc(A);
    break;
  case LIS_MATRIX_MSR:
    err = lis_matrix_merge_msr(A);
    break;
  case LIS_MATRIX_BSR:
    err = lis_matrix_merge_bsr(A);
    break;
  case LIS_MATRIX_ELL:
    err = lis_matrix_merge_ell(A);
    break;
  case LIS_MATRIX_JAD:
    err = lis_matrix_merge_jad(A);
    break;
  case LIS_MATRIX_DIA:
    err = lis_matrix_merge_dia(A);
    break;
  case LIS_MATRIX_BSC:
    err = lis_matrix_merge_bsc(A);
    break;
  case LIS_MATRIX_DNS:
    err = lis_matrix_merge_dns(A);
    break;
  case LIS_MATRIX_COO:
    err = lis_matrix_merge_coo(A);
    break;
  case LIS_MATRIX_VBR:
    err = lis_matrix_merge_vbr(A);
    break;
  default:
    LIS_SETERR_IMP;
    return LIS_ERR_NOT_IMPLEMENTED;
    break;
  }

  if( err ) return err;
  if( !A->is_save )
  {
    lis_matrix_DLU_destroy(A);
    A->is_splited = LIS_FALSE;
  }


  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_copy"
LIS_INT lis_matrix_copy(LIS_MATRIX Ain, LIS_MATRIX Aout)
{
  LIS_INT err;

  LIS_DEBUG_FUNC_IN;

  err = lis_matrix_check(Ain,LIS_MATRIX_CHECK_ALL);
  if( err ) return err;
  err = lis_matrix_check(Ain,LIS_MATRIX_CHECK_NULL);
  if( err ) return err;

  switch( Ain->matrix_type )
  {
  case LIS_MATRIX_CSR:
    err = lis_matrix_copy_csr(Ain,Aout);
    break;
  case LIS_MATRIX_CSC:
    err = lis_matrix_copy_csc(Ain,Aout);
    break;
  case LIS_MATRIX_MSR:
    err = lis_matrix_copy_msr(Ain,Aout);
    break;
  case LIS_MATRIX_DIA:
    err = lis_matrix_copy_dia(Ain,Aout);
    break;
  case LIS_MATRIX_ELL:
    err = lis_matrix_copy_ell(Ain,Aout);
    break;
  case LIS_MATRIX_JAD:
    err = lis_matrix_copy_jad(Ain,Aout);
    break;
  case LIS_MATRIX_BSR:
    err = lis_matrix_copy_bsr(Ain,Aout);
    break;
  case LIS_MATRIX_VBR:
    err = lis_matrix_copy_vbr(Ain,Aout);
    break;
  case LIS_MATRIX_DNS:
    err = lis_matrix_copy_dns(Ain,Aout);
    break;
  case LIS_MATRIX_COO:
    err = lis_matrix_copy_coo(Ain,Aout);
    break;
  case LIS_MATRIX_BSC:
    err = lis_matrix_copy_bsc(Ain,Aout);
    break;
    default:
      LIS_SETERR_IMP;
      return LIS_ERR_NOT_IMPLEMENTED;
      break;
  }
  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_copyDLU"
LIS_INT lis_matrix_copyDLU(LIS_MATRIX Ain, LIS_MATRIX_DIAG *D, LIS_MATRIX *L, LIS_MATRIX *U)
{
  LIS_INT err;

  LIS_DEBUG_FUNC_IN;

  err = lis_matrix_check(Ain,LIS_MATRIX_CHECK_ALL);
  if( err ) return err;

  switch( Ain->matrix_type )
  {
  case LIS_MATRIX_CSR:
    err = lis_matrix_copyDLU_csr(Ain,D,L,U);
    break;
/*
  case LIS_MATRIX_CSC:
    err = lis_matrix_copy_csc(Ain,Aout);
    break;
  case LIS_MATRIX_MSR:
    err = lis_matrix_copy_msr(Ain,Aout);
    break;
  case LIS_MATRIX_DIA:
    err = lis_matrix_copy_dia(Ain,Aout);
    break;
  case LIS_MATRIX_ELL:
    err = lis_matrix_copy_ell(Ain,Aout);
    break;
  case LIS_MATRIX_JAD:
    err = lis_matrix_copy_jad(Ain,Aout);
    break;
  case LIS_MATRIX_BSR:
    err = lis_matrix_copy_bsr(Ain,Aout);
    break;
  case LIS_MATRIX_BSC:
    err = lis_matrix_copy_bsc(Ain,Aout);
    break;
  case LIS_MATRIX_VBR:
    err = lis_matrix_copy_vbr(Ain,Aout);
    break;
  case LIS_MATRIX_DNS:
    err = lis_matrix_copy_dns(Ain,Aout);
    break;
  case LIS_MATRIX_COO:
    err = lis_matrix_copy_coo(Ain,Aout);
    break;
*/
    default:
      LIS_SETERR_IMP;
      *D = NULL;
      *L = NULL;
      *U = NULL;
      return LIS_ERR_NOT_IMPLEMENTED;
      break;
  }
  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_solve"
LIS_INT lis_matrix_solve(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, LIS_INT flag)
{
  LIS_DEBUG_FUNC_IN;

  if( !A->is_splited ) lis_matrix_split(A);

  switch( A->matrix_type )
  {
  case LIS_MATRIX_CSR:
    lis_matrix_solve_csr(A,b,x,flag);
    break;
  case LIS_MATRIX_BSR:
    lis_matrix_solve_bsr(A,b,x,flag);
    break;
  case LIS_MATRIX_CSC:
    lis_matrix_solve_csc(A,b,x,flag);
    break;
  case LIS_MATRIX_MSR:
    lis_matrix_solve_msr(A,b,x,flag);
    break;
  case LIS_MATRIX_ELL:
    lis_matrix_solve_ell(A,b,x,flag);
    break;
  case LIS_MATRIX_JAD:
    lis_matrix_solve_jad(A,b,x,flag);
    break;
  case LIS_MATRIX_DIA:
    lis_matrix_solve_dia(A,b,x,flag);
    break;
  case LIS_MATRIX_DNS:
    lis_matrix_solve_dns(A,b,x,flag);
    break;
  case LIS_MATRIX_BSC:
    lis_matrix_solve_bsc(A,b,x,flag);
    break;
  case LIS_MATRIX_VBR:
    lis_matrix_solve_vbr(A,b,x,flag);
    break;
  default:
    LIS_SETERR_IMP;
    return LIS_ERR_NOT_IMPLEMENTED;
    break;
  }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_solvet"
LIS_INT lis_matrix_solvet(LIS_MATRIX A, LIS_VECTOR b, LIS_VECTOR x, LIS_INT flag)
{
  LIS_DEBUG_FUNC_IN;

  if( !A->is_splited ) lis_matrix_split(A);

  switch( A->matrix_type )
  {
  case LIS_MATRIX_CSR:
    lis_matrix_solvet_csr(A,b,x,flag);
    break;
  case LIS_MATRIX_BSR:
    lis_matrix_solvet_bsr(A,b,x,flag);
    break;
  case LIS_MATRIX_CSC:
    lis_matrix_solvet_csc(A,b,x,flag);
    break;
  case LIS_MATRIX_MSR:
    lis_matrix_solvet_msr(A,b,x,flag);
    break;
  case LIS_MATRIX_ELL:
    lis_matrix_solvet_ell(A,b,x,flag);
    break;
  case LIS_MATRIX_JAD:
    lis_matrix_solvet_jad(A,b,x,flag);
    break;
  case LIS_MATRIX_DIA:
    lis_matrix_solvet_dia(A,b,x,flag);
    break;
  case LIS_MATRIX_DNS:
    lis_matrix_solvet_dns(A,b,x,flag);
    break;
  case LIS_MATRIX_BSC:
    lis_matrix_solvet_bsc(A,b,x,flag);
    break;
  case LIS_MATRIX_VBR:
    lis_matrix_solvet_vbr(A,b,x,flag);
    break;
  default:
    LIS_SETERR_IMP;
    return LIS_ERR_NOT_IMPLEMENTED;
    break;
  }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

void lis_array_LUdecomp(LIS_INT n, LIS_SCALAR *a)
{
  LIS_INT i,j,k;
  LIS_SCALAR t;

  for(k=0;k<n;k++)
  {
    a[k*n+k] = 1.0 / a[k*n+k];
    for(i=k+1;i<n;i++)
    {
      t = a[k*n+i] * a[k*n+k];
      for(j=k+1;j<n;j++)
      {
        a[j*n+i] -= t * a[j*n+k];
      }
      a[k*n+i] = t;
    }
  }
}

void lis_array_invGauss(LIS_INT n, LIS_SCALAR *a)
{
  LIS_INT i,j,k;
  LIS_SCALAR t,*lu;

  lu = (LIS_SCALAR *)lis_malloc(n*n*sizeof(LIS_SCALAR), "lis_array_invGauss::lu");
  memcpy(lu,a,n*n*sizeof(LIS_SCALAR));
  for(k=0;k<n;k++)
  {
    lu[k*n+k] = 1.0 / lu[k*n+k];
    for(i=k+1;i<n;i++)
    {
      t = lu[k*n+i] * lu[k*n+k];
      for(j=k+1;j<n;j++)
      {
        lu[j*n+i] -= t * lu[j*n+k];
      }
      lu[k*n+i] = t;
    }
  }
  for(k=0;k<n;k++)
  {
    for(i=0;i<n;i++)
    {
       t = (i==k);
       for(j=0;j<i;j++)
       {
         t -= lu[j*n+i] * a[k*n+j];
       }
       a[k*n+i] = t;
    }
    for(i=n-1;i>=0;i--)
    {
      t = a[k*n+i];
      for(j=i+1;j<n;j++)
      {
        t -= lu[j*n+i] * a[k*n+j];
      }
      a[k*n+i] = t * lu[i*n+i];
    }
  }
  free(lu);
}

void lis_array_matvec(LIS_INT n, LIS_SCALAR *a, LIS_SCALAR *b, LIS_SCALAR *c, LIS_INT op)
{
  LIS_INT      i,j;
  LIS_SCALAR  t;
  /* c = A*b */

  if( op==LIS_INS_VALUE )
  {
    switch(n)
    {
    case 1:
      c[0] = a[0]*b[0];
      break;
    case 2:
      c[0] = a[0]*b[0] + a[2]*b[1];
      c[1] = a[1]*b[0] + a[3]*b[1];
      break;
    case 3:
      c[0] = a[0]*b[0] + a[3]*b[1] + a[6]*b[2];
      c[1] = a[1]*b[0] + a[4]*b[1] + a[7]*b[2];
      c[2] = a[2]*b[0] + a[5]*b[1] + a[8]*b[2];
      break;
    default:
      for(i=0;i<n;i++)
      {
        t = 0.0;
        for(j=0;j<n;j++)
        {
          t += a[j*n+i] * b[j];
        }
        c[i] = t;
      }
      break;
    }
  }
  else if( op==LIS_SUB_VALUE )
  {
    switch(n)
    {
    case 1:
      c[0] -= a[0]*b[0];
      break;
    case 2:
      c[0] -= a[0]*b[0] + a[2]*b[1];
      c[1] -= a[1]*b[0] + a[3]*b[1];
      break;
    case 3:
      c[0] -= a[0]*b[0] + a[3]*b[1] + a[6]*b[2];
      c[1] -= a[1]*b[0] + a[4]*b[1] + a[7]*b[2];
      c[2] -= a[2]*b[0] + a[5]*b[1] + a[8]*b[2];
      break;
    default:
      for(i=0;i<n;i++)
      {
        t = 0.0;
        for(j=0;j<n;j++)
        {
          t += a[j*n+i] * b[j];
        }
        c[i] -= t;
      }
      break;
    }
  }
  else
  {
    switch(n)
    {
    case 1:
      c[0] += a[0]*b[0];
      break;
    case 2:
      c[0] += a[0]*b[0] + a[2]*b[1];
      c[1] += a[1]*b[0] + a[3]*b[1];
      break;
    case 3:
      c[0] += a[0]*b[0] + a[3]*b[1] + a[6]*b[2];
      c[1] += a[1]*b[0] + a[4]*b[1] + a[7]*b[2];
      c[2] += a[2]*b[0] + a[5]*b[1] + a[8]*b[2];
      break;
    default:
      for(i=0;i<n;i++)
      {
        t = 0.0;
        for(j=0;j<n;j++)
        {
          t += a[j*n+i] * b[j];
        }
        c[i] += t;
      }
      break;
    }
  }
}

void lis_array_matvect(LIS_INT n, LIS_SCALAR *a, LIS_SCALAR *b, LIS_SCALAR *c, LIS_INT op)
{
  LIS_INT      i,j;
  LIS_SCALAR  t;
  /* c = A*b */

  if( op==LIS_INS_VALUE )
  {
    switch(n)
    {
    case 1:
      c[0] = a[0]*b[0];
      break;
    case 2:
      c[0] = a[0]*b[0] + a[1]*b[1];
      c[1] = a[2]*b[0] + a[3]*b[1];
      break;
    case 3:
      c[0] = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
      c[1] = a[3]*b[0] + a[4]*b[1] + a[5]*b[2];
      c[2] = a[6]*b[0] + a[7]*b[1] + a[8]*b[2];
      break;
    default:
      for(i=0;i<n;i++)
      {
        t = 0.0;
        for(j=0;j<n;j++)
        {
          t += a[i*n+j] * b[j];
        }
        c[i] = t;
      }
      break;
    }
  }
  else if( op==LIS_SUB_VALUE )
  {
    switch(n)
    {
    case 1:
      c[0] -= a[0]*b[0];
      break;
    case 2:
      c[0] -= a[0]*b[0] + a[1]*b[1];
      c[1] -= a[2]*b[0] + a[3]*b[1];
      break;
    case 3:
      c[0] -= a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
      c[1] -= a[3]*b[0] + a[4]*b[1] + a[5]*b[2];
      c[2] -= a[6]*b[0] + a[7]*b[1] + a[8]*b[2];
      break;
    default:
      for(i=0;i<n;i++)
      {
        t = 0.0;
        for(j=0;j<n;j++)
        {
          t += a[i*n+j] * b[j];
        }
        c[i] -= t;
      }
      break;
    }
  }
  else
  {
    switch(n)
    {
    case 1:
      c[0] += a[0]*b[0];
      break;
    case 2:
      c[0] += a[0]*b[0] + a[1]*b[1];
      c[1] += a[2]*b[0] + a[3]*b[1];
      break;
    case 3:
      c[0] += a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
      c[1] += a[3]*b[0] + a[4]*b[1] + a[5]*b[2];
      c[2] += a[6]*b[0] + a[7]*b[1] + a[8]*b[2];
      break;
    default:
      for(i=0;i<n;i++)
      {
        t = 0.0;
        for(j=0;j<n;j++)
        {
          t += a[i*n+j] * b[j];
        }
        c[i] += t;
      }
      break;
    }
  }
}

void lis_array_matmat(LIS_INT n, LIS_SCALAR *a, LIS_SCALAR *b, LIS_SCALAR *c, LIS_INT op)
{
  /* C = A*B */

  if( op==LIS_INS_VALUE )
  {
    switch(n)
    {
    case 1:
      c[0] = a[0]*b[0];
      break;
    case 2:
      c[0] = a[0]*b[0] + a[2]*b[1];
      c[1] = a[1]*b[0] + a[3]*b[1];
      c[2] = a[0]*b[2] + a[2]*b[3];
      c[3] = a[1]*b[2] + a[3]*b[3];
      break;
    case 3:
      c[0] = a[0]*b[0] + a[3]*b[1] + a[6]*b[2];
      c[1] = a[1]*b[0] + a[4]*b[1] + a[7]*b[2];
      c[2] = a[2]*b[0] + a[5]*b[1] + a[8]*b[2];
      c[3] = a[0]*b[3] + a[3]*b[4] + a[6]*b[5];
      c[4] = a[1]*b[3] + a[4]*b[4] + a[7]*b[5];
      c[5] = a[2]*b[3] + a[5]*b[4] + a[8]*b[5];
      c[6] = a[0]*b[6] + a[3]*b[7] + a[6]*b[8];
      c[7] = a[1]*b[6] + a[4]*b[7] + a[7]*b[8];
      c[8] = a[2]*b[6] + a[5]*b[7] + a[8]*b[8];
      break;
    }
  }
  else if( op==LIS_SUB_VALUE )
  {
    switch(n)
    {
    case 1:
      c[0] -= a[0]*b[0];
      break;
    case 2:
      c[0] -= a[0]*b[0] + a[2]*b[1];
      c[1] -= a[1]*b[0] + a[3]*b[1];
      c[2] -= a[0]*b[2] + a[2]*b[3];
      c[3] -= a[1]*b[2] + a[3]*b[3];
      break;
    case 3:
      c[0] -= a[0]*b[0] + a[3]*b[1] + a[6]*b[2];
      c[1] -= a[1]*b[0] + a[4]*b[1] + a[7]*b[2];
      c[2] -= a[2]*b[0] + a[5]*b[1] + a[8]*b[2];
      c[3] -= a[0]*b[3] + a[3]*b[4] + a[6]*b[5];
      c[4] -= a[1]*b[3] + a[4]*b[4] + a[7]*b[5];
      c[5] -= a[2]*b[3] + a[5]*b[4] + a[8]*b[5];
      c[6] -= a[0]*b[6] + a[3]*b[7] + a[6]*b[8];
      c[7] -= a[1]*b[6] + a[4]*b[7] + a[7]*b[8];
      c[8] -= a[2]*b[6] + a[5]*b[7] + a[8]*b[8];
      break;
    }
  }
  else
  {
    switch(n)
    {
    case 1:
      c[0] += a[0]*b[0];
      break;
    case 2:
      c[0] += a[0]*b[0] + a[2]*b[1];
      c[1] += a[1]*b[0] + a[3]*b[1];
      c[2] += a[0]*b[2] + a[2]*b[3];
      c[3] += a[1]*b[2] + a[3]*b[3];
      break;
    case 3:
      c[0] += a[0]*b[0] + a[3]*b[1] + a[6]*b[2];
      c[1] += a[1]*b[0] + a[4]*b[1] + a[7]*b[2];
      c[2] += a[2]*b[0] + a[5]*b[1] + a[8]*b[2];
      c[3] += a[0]*b[3] + a[3]*b[4] + a[6]*b[5];
      c[4] += a[1]*b[3] + a[4]*b[4] + a[7]*b[5];
      c[5] += a[2]*b[3] + a[5]*b[4] + a[8]*b[5];
      c[6] += a[0]*b[6] + a[3]*b[7] + a[6]*b[8];
      c[7] += a[1]*b[6] + a[4]*b[7] + a[7]*b[8];
      c[8] += a[2]*b[6] + a[5]*b[7] + a[8]*b[8];
      break;
    }
  }
}

void lis_array_matmat2(LIS_INT m, LIS_INT n, LIS_INT k, LIS_SCALAR *a, LIS_INT lda, LIS_SCALAR *b, LIS_INT ldb, LIS_SCALAR *c, LIS_INT ldc, LIS_INT op)
{
  LIS_INT      i,j,l;
  /* C = A*B */

  if( op==LIS_INS_VALUE )
  {
    for(j=0;j<n;j++)
    {
      for(i=0;i<m;i++)
      {
        c[j*ldc+i] = 0.0;
      }
      for(l=0;l<k;l++)
      {
        for(i=0;i<m;i++)
        {
          c[j*ldc+i] += b[j*ldb+l] * a[l*lda+i];
        }
      }
    }
  }
  else if( op==LIS_SUB_VALUE )
  {
    for(j=0;j<n;j++)
    {
      for(l=0;l<k;l++)
      {
        for(i=0;i<m;i++)
        {
          c[j*ldc+i] -= b[j*ldb+l] * a[l*lda+i];
        }
      }
    }
  }
  else
  {
    switch(n)
    {
    case 1:
      c[0] += a[0]*b[0];
      break;
    case 2:
      c[0] += a[0]*b[0] + a[2]*b[1];
      c[1] += a[1]*b[0] + a[3]*b[1];
      c[2] += a[0]*b[2] + a[2]*b[3];
      c[3] += a[1]*b[2] + a[3]*b[3];
      break;
    case 3:
      c[0] += a[0]*b[0] + a[3]*b[1] + a[6]*b[2];
      c[1] += a[1]*b[0] + a[4]*b[1] + a[7]*b[2];
      c[2] += a[2]*b[0] + a[5]*b[1] + a[8]*b[2];
      c[3] += a[0]*b[3] + a[3]*b[4] + a[6]*b[5];
      c[4] += a[1]*b[3] + a[4]*b[4] + a[7]*b[5];
      c[5] += a[2]*b[3] + a[5]*b[4] + a[8]*b[5];
      c[6] += a[0]*b[6] + a[3]*b[7] + a[6]*b[8];
      c[7] += a[1]*b[6] + a[4]*b[7] + a[7]*b[8];
      c[8] += a[2]*b[6] + a[5]*b[7] + a[8]*b[8];
      break;
    }
  }
}

void lis_array_nrm2(LIS_INT n, LIS_SCALAR *v, LIS_SCALAR *nrm2)
{
  LIS_INT i;
  LIS_SCALAR t;

  t = 0.0;
  for(i=0;i<n;i++)
  {
    t += v[i]*v[i];
  }
  *nrm2 = sqrt(t);
}

void lis_array_nrm1(LIS_INT n, LIS_SCALAR *v, LIS_SCALAR *nrm1)
{
  LIS_INT i;
  LIS_SCALAR t;

  t = 0.0;
  for(i=0;i<n;i++)
  {
    t += fabs(v[i]);
  }
  *nrm1 = t;
}

void lis_array_dot(LIS_INT n, LIS_SCALAR *v, LIS_SCALAR *dot)
{
  LIS_INT i;
  LIS_SCALAR t;

  t = 0.0;
  for(i=0;i<n;i++)
  {
    t += v[i]*v[i];
  }
  *dot = t;
}

void lis_array_matinv(LIS_INT n, LIS_SCALAR *a, LIS_SCALAR *b, LIS_SCALAR *c)
{
  LIS_INT i,j,k;
  LIS_SCALAR t;

  for(i=0;i<n;i++)
  {
    c[i] = -b[i] * a[0];
    for(j=1;j<n;j++)
    {
      t = -b[j*n+i];
      for(k=0;k<j-1;k++)
      {
        t -= c[k*n+i] * a[j*n+k];
      }
      c[j*n+i] = t * a[j*n+j];
    }
  }
  for(i=0;i<n;i++)
  {
    for(j=n-1;j>=0;j--)
    {
      t = c[j*n+i];
      for(k=j+1;k<n;k++)
      {
        t -= c[k*n+i] * a[j*n+k];
      }
      c[j*n+i] = t;
    }
  }
}

void lis_array_invvec(LIS_INT n, LIS_SCALAR *a, LIS_SCALAR *x, LIS_SCALAR *y)
{
  /* y = inv(a) * x */
  LIS_INT i,j;
  LIS_SCALAR t;

  for(i=0;i<n;i++)
  {
    t = x[i];
    for(j=0;j<i;j++)
    {
      t -= a[j*n+i] * y[j];
    }
    y[i] = t;
  }
  for(i=n-1;i>=0;i--)
  {
    t = y[i];
    for(j=i+1;j<n;j++)
    {
      t -= a[j*n+i] * y[j];
    }
    y[i] = t * a[i*n+i];
  }
}

void lis_array_invvect(LIS_INT n, LIS_SCALAR *a, LIS_SCALAR *x, LIS_SCALAR *y)
{
  /* y = inv(a) * x */
  LIS_INT i,j;

  for(i=0;i<n;i++) y[i] = x[i];
  for(i=0;i<n;i++)
  {
    y[i] = a[i*n+i] * y[i];
    for(j=i+1;j<n;j++)
    {
      y[j] -= a[j*n+i] * y[i];
    }
  }
  for(i=n-1;i>=0;i--)
  {
    for(j=0;j<i;j++)
    {
      y[j] -= a[j*n+i] * y[i];
    }
  }
}

void lis_array_matvec2(LIS_INT m, LIS_INT n, LIS_SCALAR *a, LIS_INT lda, LIS_SCALAR *b, LIS_SCALAR *c, LIS_INT op)
{
  LIS_INT      i,j;
  LIS_SCALAR  t;
  /* c = A*b */

  if( op==LIS_INS_VALUE )
  {
    for(i=0;i<m;i++)
    {
      t = 0.0;
      for(j=0;j<n;j++)
      {
        t += a[j*lda+i] * b[j];
      }
      c[i] = t;
    }
  }
  else if( op==LIS_SUB_VALUE )
  {
    for(i=0;i<m;i++)
    {
      t = 0.0;
      for(j=0;j<n;j++)
      {
        t += a[j*lda+i] * b[j];
      }
      c[i] -= t;
    }
  }
  else if( op==LIS_ADD_VALUE )
  {
    for(i=0;i<m;i++)
    {
      t = 0.0;
      for(j=0;j<n;j++)
      {
        t += a[j*lda+i] * b[j];
      }
      c[i] += t;
    }
  }
  else
  {
    switch(n)
    {
    case 1:
      c[0] += a[0]*b[0];
      break;
    case 2:
      c[0] += a[0]*b[0] + a[2]*b[1];
      c[1] += a[1]*b[0] + a[3]*b[1];
      break;
    case 3:
      c[0] += a[0]*b[0] + a[3]*b[1] + a[6]*b[2];
      c[1] += a[1]*b[0] + a[4]*b[1] + a[7]*b[2];
      c[2] += a[2]*b[0] + a[5]*b[1] + a[8]*b[2];
      break;
    default:
      for(i=0;i<n;i++)
      {
        t = 0.0;
        for(j=0;j<n;j++)
        {
          t += a[j*n+i] * b[j];
        }
        c[i] += t;
      }
      break;
    }
  }
}

void lis_array_solve(LIS_INT n, LIS_SCALAR *aa, LIS_SCALAR *b, LIS_SCALAR *x,
LIS_SCALAR *a)
{
  LIS_INT i,j,k,imax,swap;
  LIS_SCALAR t,al;

  for(i=0;i<n*n;i++) a[i] = aa[i];

  switch( n )
  {
  case 1:
    x[0] = b[0] / a[0];
    break;
  case 2:
    a[0]  = 1.0 / a[0];
    a[1] *= a[0];
    a[3] -= a[1] * a[2];
    a[3]  = 1.0 / a[3];
    /* forward sub */
    x[0] = b[0];
    x[1] = b[1] - a[1] * x[0];
    /* backward sub */
    x[1] *= a[3];
    x[0] -= a[2] * x[1];
    x[0] *= a[0];
    break;
  default:
    for(k=0;k<n;k++)
    {
      a[k*n+k] = 1.0 / a[k*n+k];
      for(i=k+1;i<n;i++)
      {
        t = a[k*n+i] * a[k*n+k];
        for(j=k+1;j<n;j++)
        {
          a[j*n+i] -= t * a[j*n+k];
        }
        a[k*n+i] = t;
      }
    }

    /* forward sub */
    for(i=0;i<n;i++)
    {
      x[i] = b[i];
      for(j=0;j<i;j++)
      {
        x[i] -= a[j*n+i] * x[j];
      }
    }
    /* backward sub */
    for(i=n-1;i>=0;i--)
    {
      for(j=i+1;j<n;j++)
      {
        x[i] -= a[j*n+i] * x[j];
      }
      x[i] *= a[i*n+i];
    }
    break;
  }
}

#undef __FUNC__
#define __FUNC__ "lis_array_cgs"
LIS_INT lis_array_cgs(LIS_INT n, LIS_SCALAR *x, LIS_SCALAR *q, LIS_SCALAR *r)
{
  LIS_INT i, j, k; 
  LIS_SCALAR *x_k, nrm2;
  LIS_REAL tol;

  tol = 1e-12;

  x_k = (LIS_SCALAR *)lis_malloc(n*sizeof(LIS_SCALAR), "lis_array_cgs::x_k");

  for (i=0;i<n*n;i++)
    {
    q[i] = 0.0;
    r[i] = 0.0;
    }

  for (k=0;k<n;k++)
    {
      for (i=0;i<n;i++)
  {
    x_k[i] = x[i*n+k];
  }
      for (j=0;j<k;j++)
  {
    r[j*n+k] = 0; 
    for (i=0;i<n;i++)
      {
        r[j*n+k] += q[i*n+j] * x[i*n+k];
      }
    for (i=0;i<n;i++)
      {
        x_k[i] -= r[j*n+k] * q[i*n+j];
      }
  }
      lis_array_nrm2(n, &x_k[0], &nrm2);
      r[k*n+k] = nrm2;
      if (nrm2<tol) break; 
      for (i=0;i<n;i++)
  {
    q[i*n+k] = x_k[i] / nrm2;
  }
      
    }

  lis_free(x_k);

  return 0;
} 

#undef __FUNC__
#define __FUNC__ "lis_array_mgs"
LIS_INT lis_array_mgs(LIS_INT n, LIS_SCALAR *x, LIS_SCALAR *q, LIS_SCALAR *r)
{
  LIS_INT i, j, k; 
  LIS_SCALAR *x_j, nrm2;
  LIS_REAL tol;

  tol = 1e-12;

  x_j = (LIS_SCALAR *)lis_malloc(n*sizeof(LIS_SCALAR), "lis_array_mgs::x_j");

  for (j=0;j<n;j++)
    {
      for (i=0;i<n;i++)
  {
    x_j[i] = x[i*n+j];
  }
      lis_array_nrm2(n, &x_j[0], &nrm2);
      r[j*n+j] = nrm2;
      for (i=0;i<n;i++)
  {
    if (nrm2<tol) break; 
    q[i*n+j] = x_j[i] / nrm2;
  }
      for (k=j+1;k<n;k++)
  {
    r[j*n+k] = 0; 
    for (i=0;i<n;i++)
      {
        r[j*n+k] += q[i*n+j] * x[i*n+k];
      }
    for (i=0;i<n;i++)
      {
        x[i*n+k] -= r[j*n+k] * q[i*n+j];
      }
  }
    }
  lis_free(x_j);
  return 0;
}

#undef __FUNC__
#define __FUNC__ "lis_array_qr"
LIS_INT lis_array_qr(LIS_INT n, LIS_SCALAR *x, LIS_SCALAR *q, LIS_SCALAR *r)
{
  LIS_INT i, j, k, iter, maxiter; 
  LIS_SCALAR *x0, x_tmp;
  LIS_REAL err, err0, tol;

  maxiter = 100000;
  tol = 1e-12;

  x0 = (LIS_SCALAR *)lis_malloc(n*n*sizeof(LIS_SCALAR), "lis_array_qr::x0");
  iter = 0;
  while (iter < maxiter)
    {
      iter = iter + 1;
      lis_array_cgs(n, x, q, r);
      for (i=0;i<n;i++)
  {
    for (j=0;j<n;j++)
      {
        x[i*n+j] = 0;
        for (k=0;k<n;k++)
    {
      x[i*n+j] += r[i*n+k] * q[k*n+j];
    }
      }
  }
      err = sqrt(x[n*n-2] * x[n*n-2]); 
      if (err<tol) break;
    }

  lis_free(x0);
  return 0; 
}

void lis_array_power(LIS_INT n, LIS_SCALAR *a, LIS_SCALAR *x, LIS_SCALAR *mu, LIS_INT maxiter, LIS_REAL tol, LIS_REAL *err)
{
  LIS_INT               i, iter;
  LIS_REAL          nrm2;
  LIS_SCALAR        *z,*q;
  double            times, ptimes;

  iter=0;
  lis_array_set_all(n, 1.0, x);
  z = (LIS_SCALAR *)lis_malloc(n*sizeof(LIS_SCALAR), "lis_array_power::z");
  q = (LIS_SCALAR *)lis_malloc(n*sizeof(LIS_SCALAR), "lis_array_power::q");
  times = lis_wtime();
  while (iter<maxiter)
    {
      iter = iter+1;
      lis_array_nrm2(n, x, &nrm2);
      lis_array_scale(n, 1/nrm2, x);
      lis_array_matvec(n, a, x, z, LIS_INS_VALUE);
      lis_array_dot2(n, x, z, mu); 
      lis_array_axpyz(n, -(*mu),x,z,q); 
      lis_array_nrm2(n, q, err); 
      *err = fabs(*err / (*mu));
      if (*err<tol) break;  
      lis_array_copy(n, z, x);
    }
  lis_free(z);
  lis_free(q);
}

void lis_array_set_all(LIS_INT n, LIS_SCALAR alpha, LIS_SCALAR *v)
{
  LIS_INT i;

  for(i=0;i<n;i++)
    {
      v[i] = alpha; 
    }
}

void lis_array_scale(LIS_INT n, LIS_SCALAR alpha, LIS_SCALAR *v)
{
  LIS_INT i;

  for(i=0;i<n;i++)
    {
      v[i] = alpha * v[i];
    }
}

void lis_array_dot2(LIS_INT n, LIS_SCALAR *x, LIS_SCALAR *y, LIS_SCALAR *alpha)
{
  LIS_INT i;

  *alpha = 0;
  for(i=0;i<n;i++)
    {
      *alpha = *alpha + x[i] * y[i];
    }
}

void lis_array_axpyz(LIS_INT n, LIS_SCALAR alpha, LIS_SCALAR *x, LIS_SCALAR *y, LIS_SCALAR *z)
{
  LIS_INT i;

  for(i=0;i<n;i++)
    {
      z[i] = alpha * x[i] + y[i];
    }
}

void lis_array_copy(LIS_INT n, LIS_SCALAR *x, LIS_SCALAR *y)
{
  LIS_INT i;

  for(i=0;i<n;i++)
    {
      y[i] = x[i];
    }
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_shift_diagonal"
LIS_INT lis_matrix_shift_diagonal(LIS_MATRIX A, LIS_SCALAR shift)
{

  LIS_DEBUG_FUNC_IN;

  switch( A->matrix_type )
  {
  case LIS_MATRIX_CSR:
    lis_matrix_shift_diagonal_csr(A, shift);
    break;
  case LIS_MATRIX_CSC:
    lis_matrix_shift_diagonal_csc(A, shift);
    break;
  case LIS_MATRIX_MSR:
    lis_matrix_shift_diagonal_msr(A, shift);
    break;
  case LIS_MATRIX_DIA:
    lis_matrix_shift_diagonal_dia(A, shift);
    break;
  case LIS_MATRIX_ELL:
    lis_matrix_shift_diagonal_ell(A, shift);
    break;
  case LIS_MATRIX_JAD:
    lis_matrix_shift_diagonal_jad(A, shift);
    break;
  case LIS_MATRIX_BSR:
    lis_matrix_shift_diagonal_bsr(A, shift);
    break;
  case LIS_MATRIX_BSC:
    lis_matrix_shift_diagonal_bsc(A, shift);
    break;
  case LIS_MATRIX_DNS:
    lis_matrix_shift_diagonal_dns(A, shift);
    break;
  case LIS_MATRIX_COO:
    lis_matrix_shift_diagonal_coo(A, shift);
    break;
  case LIS_MATRIX_VBR:
    lis_matrix_shift_diagonal_vbr(A, shift);
    break;

  default:
    LIS_SETERR_IMP;
    return LIS_ERR_NOT_IMPLEMENTED;
    break;
  }
  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_shift_diagonal_csr"
LIS_INT lis_matrix_shift_diagonal_csr(LIS_MATRIX A, LIS_SCALAR shift)
{
  LIS_INT i,j;
  LIS_INT n;

  LIS_DEBUG_FUNC_IN;

  n    = A->n;
  if( A->is_splited )
  {
    #ifdef _OPENMP
    #pragma omp parallel for private(i,j)
    #endif
    for(i=0; i<n; i++)
    {
      A->D->value[i] += shift;
    }
  }
  else
  {
    #ifdef _OPENMP
    #pragma omp parallel for private(i,j)
    #endif
    for(i=0; i<n; i++)
    {
      for(j=A->ptr[i];j<A->ptr[i+1];j++)
      {
        if( i==A->index[j] )
        {
          A->value[j] += shift;
          break;
        }
      }
    }
  }
  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_shift_diagonal_csc"
LIS_INT lis_matrix_shift_diagonal_csc(LIS_MATRIX A, LIS_SCALAR shift)
{
  LIS_INT i,j;
  LIS_INT n,np;

  LIS_DEBUG_FUNC_IN;

  n    = A->n;
  np   = A->np;
  if( A->is_splited )
  {
    #ifdef _OPENMP
    #pragma omp parallel for private(i,j)
    #endif
    for(i=0; i<n; i++)
    {
      A->D->value[i] += shift;
    }
  }
  else
  {
    #ifdef _OPENMP
    #pragma omp parallel for private(i,j)
    #endif
    for(i=0; i<np; i++)
    {
      for(j=A->ptr[i];j<A->ptr[i+1];j++)
      {
        if( i==A->index[j] )
        {
          A->value[j] += shift;
          break;
        }
      }
    }
  }
  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_shift_diagonal_msr"
LIS_INT lis_matrix_shift_diagonal_msr(LIS_MATRIX A, LIS_SCALAR shift)
{
  LIS_INT i;
  LIS_INT n;

  LIS_DEBUG_FUNC_IN;

  n    = A->n;
  if( A->is_splited )
  {
    #ifdef _OPENMP
    #pragma omp parallel for private(i)
    #endif
    for(i=0; i<n; i++)
    {
      A->D->value[i] += shift;
    }
  }
  else
  {
    #ifdef _OPENMP
    #pragma omp parallel for private(i)
    #endif
    for(i=0; i<n; i++)
    {
      A->value[i] += shift;
    }
  }
  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_shift_diagonal_dia"
LIS_INT lis_matrix_shift_diagonal_dia(LIS_MATRIX A, LIS_SCALAR shift)
{
  LIS_INT i,j,k;
  LIS_INT n,nnd;
  #ifdef _OPENMP
    LIS_INT is,ie,my_rank,nprocs;
  #endif

  LIS_DEBUG_FUNC_IN;

  n    = A->n;
  nnd  = A->nnd;
  if( A->is_splited )
  {
    #ifdef _OPENMP
    #pragma omp parallel for private(i)
    #endif
    for(i=0; i<n; i++)
    {
      A->D->value[i] += shift;
    }
  }
  else
  {
    #ifdef _OPENMP
      nprocs  = omp_get_max_threads();
      for(j=0;j<nnd;j++)
      {
        if( A->index[j]==0 ) break;
      }
      #pragma omp parallel private(is,ie,my_rank)
      {
        my_rank = omp_get_thread_num();
        LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
        /*
        memcpy(&d[is],&A->value[is*nnd+j*(ie-is)],(ie-is)*sizeof(LIS_SCALAR));
        */
        for (k=is;k<ie;k++)
          {
            A->value[is*nnd+j*k] += shift;
          }
      }
    #else
      for(j=0;j<nnd;j++)
      {
        if( A->index[j]==0 ) break;
      }
      for(i=0;i<n;i++)
      {
        A->value[j*n+i] += shift;
      }
    #endif
  }
  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_shift_diagonal_ell"
LIS_INT lis_matrix_shift_diagonal_ell(LIS_MATRIX A, LIS_SCALAR shift)
{
  LIS_INT i,j;
  LIS_INT n,maxnzr;

  LIS_DEBUG_FUNC_IN;

  n    = A->n;
  if( A->is_splited )
  {
    #ifdef _OPENMP
    #pragma omp parallel for private(i,j)
    #endif
    for(i=0; i<n; i++)
    {
      A->D->value[i] += shift;
    }
  }
  else
  {
    maxnzr = A->maxnzr;
    #ifdef _OPENMP
    #pragma omp parallel for private(i,j)
    #endif
    for(i=0; i<n; i++)
    {
      for(j=0;j<maxnzr;j++)
      {
        if( i==A->index[j*n+i] )
        {
          A->value[j*n+i] += shift;
          break;
        }
      }
    }
  }
  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_shift_diagonal_jad"
LIS_INT lis_matrix_shift_diagonal_jad(LIS_MATRIX A, LIS_SCALAR shift)
{
  LIS_INT i,j,k,l;
  LIS_INT n,maxnzr;
  #ifdef _OPENMP
    LIS_INT is,ie,js,je;
    LIS_INT nprocs,my_rank;
  #endif

  LIS_DEBUG_FUNC_IN;

  n      = A->n;
  maxnzr = A->maxnzr;
  k      = n;
  if( A->is_splited )
  {
    #ifdef _OPENMP
    #pragma omp parallel for private(i)
    #endif
    for(i=0; i<n; i++)
    {
      A->D->value[i] += shift;
    }
  }
  else
  {
    #ifdef _OPENMP
      nprocs = omp_get_max_threads();
      #pragma omp parallel private(i,j,k,l,is,ie,js,je,my_rank)
      {
        my_rank = omp_get_thread_num();
        LIS_GET_ISIE(my_rank,nprocs,n,is,ie);
        k = ie-is;
        for(j=0;j<maxnzr;j++)
        {
          l  = is;
          js = A->ptr[my_rank*(maxnzr+1) + j];
          je = A->ptr[my_rank*(maxnzr+1) + j+1];
          for(i=js;i<je;i++)
          {
            if( A->row[l]==A->index[i] )
            {
              A->value[i] += shift;
              k--;
              if( k==0 ) goto get_diag_end;
            }
            l++;
          }
        }
        get_diag_end:
        ;
      }
    #else
      for(j=0;j<maxnzr;j++)
      {
        l = 0;
        for(i=A->ptr[j];i<A->ptr[j+1];i++)
        {
          if( A->row[l]==A->index[i] )
          {
            A->value[i] += shift;
            k--;
            if( k==0 ) return LIS_SUCCESS;
          }
          l++;
        }
      }
    #endif
  }
  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_shift_diagonal_bsr"
LIS_INT lis_matrix_shift_diagonal_bsr(LIS_MATRIX A, LIS_SCALAR shift)
{
  LIS_INT i,j,k,bi,bj,bjj,nr,nc;
  LIS_INT bnr,bnc,bs;
  LIS_INT n;

  LIS_DEBUG_FUNC_IN;

  n   = A->n;
  nr  = A->nr;
  nc  = A->nc;
  bnr = A->bnr;
  bnc = A->bnc;
  bs  = bnr*bnc;
  if( A->is_splited )
  {
    #ifdef _OPENMP
    #pragma omp parallel for private(i,j)
    #endif
    for(i=0;i<nr;i++)
    {
      for(j=0;j<bnr;j++)
      {
        A->D->value[i*bs+j*bnr+j] += shift;
      }
    }
  }
  else
  {
    #ifdef _OPENMP
    #pragma omp parallel for private(bi,bj,bjj,i,j,k)
    #endif
    for(bi=0;bi<nr;bi++)
    {
      k = 0;
      i = bi*bnr;
      for(bj=A->bptr[bi];bj<A->bptr[bi+1];bj++)
      {
        bjj = A->bindex[bj];
        if( i>=bjj*bnc && i<(bjj+1)*bnc )
        {
          for(j=i%bnc;j<bnc&&k<bnr&&i<n;j++)
          {
            A->value[bj*bs + j*bnr + k] += shift;
            i++;
            k++;
          }
        }
        if( k==bnr ) break;
      }
    }
  }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_shift_diagonal_bsc"
LIS_INT lis_matrix_shift_diagonal_bsc(LIS_MATRIX A, LIS_SCALAR shift)
{
  LIS_INT i,j,k,bi,bj,bjj,nr,nc;
  LIS_INT bnr,bnc,bs;
  LIS_INT n;

  LIS_DEBUG_FUNC_IN;

  n   = A->n;
  nr  = A->nr;
  nc  = A->nc;
  bnr = A->bnr;
  bnc = A->bnc;
  bs  = bnr*bnc;
  if( A->is_splited )
  {
    #ifdef _OPENMP
    #pragma omp parallel for private(i,j)
    #endif
    for(i=0;i<nr;i++)
    {
      for(j=0;j<bnr;j++)
      {
        A->D->value[i*bs+j*bnr+j] += shift;
      }
    }
  }
  else
  {
    #ifdef _OPENMP
    #pragma omp parallel for private(bi,bj,bjj,i,j,k)
    #endif
    for(bi=0;bi<nr;bi++)
    {
      k = 0;
      i = bi*bnr;
      for(bj=A->bptr[bi];bj<A->bptr[bi+1];bj++)
      {
        bjj = A->bindex[bj];
        if( i>=bjj*bnc && i<(bjj+1)*bnc )
        {
          for(j=i%bnc;j<bnc&&k<bnr&&i<n;j++)
          {
            A->value[bj*bs + j*bnr + k] += shift;
            i++;
            k++;
          }
        }
        if( k==bnr ) break;
      }
    }
  }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_shift_diagonal_dns"
LIS_INT lis_matrix_shift_diagonal_dns(LIS_MATRIX A, LIS_SCALAR shift)
{
  LIS_INT i;
  LIS_INT n;

  LIS_DEBUG_FUNC_IN;

  n    = A->n;

  if( A->is_splited )
  {
    #ifdef _OPENMP
    #pragma omp parallel for private(i)
    #endif
    for(i=0; i<n; i++)
    {
      A->D->value[i] += shift;
    }
  }
  else
  {
    #ifdef _OPENMP
    #pragma omp parallel for private(i)
    #endif
    for(i=0; i<n; i++)
    {
      A->value[i*n + i] += shift;
    }
  }
  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_shift_diagonal_coo"
LIS_INT lis_matrix_shift_diagonal_coo(LIS_MATRIX A, LIS_SCALAR shift)
{
  LIS_INT i;
  LIS_INT n,nnz;

  LIS_DEBUG_FUNC_IN;

  n    = A->n;
  nnz  = A->nnz;
  if( A->is_splited )
  {
    #ifdef _OPENMP
    #pragma omp parallel for private(i)
    #endif
    for(i=0; i<n; i++)
    {
      A->D->value[i] += shift;
    }
  }
  else
  {
    for(i=0; i<nnz;i++)
    {
      if( A->row[i]==A->col[i] )
      {
        A->value[i] += shift;
      }
    }
  }
  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

#undef __FUNC__
#define __FUNC__ "lis_matrix_shift_diagonal_vbr"
LIS_INT lis_matrix_shift_diagonal_vbr(LIS_MATRIX A, LIS_SCALAR shift)
{
  LIS_INT i,j,k,bi,bj,bjj,nr,nc;
  LIS_INT bnr,bnc;
  LIS_INT n;

  LIS_DEBUG_FUNC_IN;


  n   = A->n;
  nr  = A->nr;
  nc  = A->nc;

  if( A->is_splited )
  {
    #ifdef _OPENMP
    #pragma omp parallel for private(i,j,bnr)
    #endif
    for(i=0;i<nr;i++)
    {
      bnr = A->D->bns[i];
      for(j=0;j<bnr;j++)
      {
        A->D->v_value[i][j*bnr+j] += shift;
      }

    }
  }
  else
  {
    #ifdef _OPENMP
    #pragma omp parallel for private(bi,bj,bjj,bnr,bnc,i,j,k)
    #endif
    for(bi=0;bi<nr;bi++)
    {
      k = 0;
      i = A->row[bi];
      bnr = A->row[bi+1] - A->row[bi];
      for(bj=A->bptr[bi];bj<A->bptr[bi+1];bj++)
      {
        bjj = A->bindex[bj];
        bnc = A->col[bjj+1] - A->col[bjj];
        if( i>=bjj*bnc && i<(bjj+1)*bnc )
        {
          for(j=i%bnc;j<bnc&&k<bnr&&i<n;j++)
          {
            A->value[A->ptr[bj] + j*bnr + k] += shift;
            i++;
            k++;
          }
        }
        if( k==bnr ) break;
      }
    }
  }

  LIS_DEBUG_FUNC_OUT;
  return LIS_SUCCESS;
}

