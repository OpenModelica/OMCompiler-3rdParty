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
#include <string.h>
#include <math.h>
#include "lis.h"


#undef __FUNC__
#define __FUNC__ "main"
LIS_INT main(LIS_INT argc, char* argv[])
{
  LIS_MATRIX    A;
  LIS_VECTOR    b,x;
  LIS_INT      nprocs,my_rank;
  int                     int_nprocs,int_my_rank;
  LIS_INT      format;
  LIS_INT      err;

  LIS_DEBUG_FUNC_IN;


  lis_initialize(&argc, &argv);

  #ifdef USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD,&int_nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&int_my_rank);
    nprocs = int_nprocs;
    my_rank = int_my_rank;
  #else
    nprocs  = 1;
    my_rank = 0;
  #endif

  if( argc < 4 )
  {
err_print:
    if( my_rank==0 ) 
      {
        printf("%s out_format out_matrix_filename in_matrix_filename [in_rhs_file in_init_file]\n", argv[0]);
      }
    CHKERR(1);
  }

  if( strncmp(argv[1], "mmb", 3)==0 )
  {
    format = LIS_FMT_MMB;
  }
  else if( strncmp(argv[1], "mm", 2)==0 )
  {
    format = LIS_FMT_MM;
  }
  else
  {
    goto err_print;
  }

  if( my_rank==0 )
    {
      printf("\n");
#ifdef _LONGLONG
      printf("number of processes = %lld\n",nprocs);
#else
      printf("number of processes = %d\n",nprocs);
#endif
    }

#ifdef _OPENMP
  if( my_rank==0 )
    {
#ifdef _LONGLONG
      printf("max number of threads = %lld\n",omp_get_num_procs());
      printf("number of threads = %lld\n",omp_get_max_threads());
#else
      printf("max number of threads = %d\n",omp_get_num_procs());
      printf("number of threads = %d\n",omp_get_max_threads());
#endif
    }
#endif
    
  /* read matrix and vectors from file */
  err = lis_matrix_create(LIS_COMM_WORLD,&A); CHKERR(err);
  err = lis_vector_create(LIS_COMM_WORLD,&b); CHKERR(err);
  err = lis_vector_create(LIS_COMM_WORLD,&x); CHKERR(err);
  err = lis_input(A,b,x,argv[3]);
  CHKERR(err);
  if( argc>4 )
  {
    if( !lis_vector_is_null(b) )
    {
      lis_vector_destroy(b);
      err = lis_vector_create(LIS_COMM_WORLD,&b); CHKERR(err);
    }
    err = lis_input_vector(b,argv[4]);
    CHKERR(err);
    if( A->n!=b->n )
    {
#ifdef _LONGLONG
      printf("different dimension: matrix A=%lld, vector b=%lld\n",A->n,b->n);
#else
      printf("different dimension: matrix A=%d, vector b=%d\n",A->n,b->n);
#endif
      goto err_to;
    }
    if( argc==6 )
    {
      if( !lis_vector_is_null(x) )
      {
        lis_vector_destroy(x);
        err = lis_vector_create(LIS_COMM_WORLD,&x); CHKERR(err);
      }
      err = lis_input_vector(x,argv[5]);
      CHKERR(err);
      if( A->n!=x->n )
      {
#ifdef _LONGLONG
        printf("different dimension: matrix A=%lld, vector x=%lld\n",A->n,x->n);
#else
        printf("different dimension: matrix A=%d, vector x=%d\n",A->n,x->n);
#endif
        goto err_to;
      }
    }
  }

  err = lis_output(A,b,x,format,argv[2]);
  CHKERR(err);

err_to:
  lis_matrix_destroy(A);
  lis_vector_destroy(b);
  lis_vector_destroy(x);

  lis_finalize();

  LIS_DEBUG_FUNC_OUT;
  return 0;
}

