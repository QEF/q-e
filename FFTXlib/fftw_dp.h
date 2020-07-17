#ifndef FFTW_DP_H
#define FFTW_DP_H

/*
 * Copyright (c) 1997 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to use, copy, modify, and distribute the Software without
 * restriction, provided the Software, including any modified copies made
 * under this license, is not distributed for a fee, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE MASSACHUSETTS INSTITUTE OF TECHNOLOGY BE LIABLE
 * FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * 
 * Except as contained in this notice, the name of the Massachusetts
 * Institute of Technology shall not be used in advertising or otherwise
 * to promote the sale, use or other dealings in this Software without
 * prior written authorization from the Massachusetts Institute of
 * Technology.
 *  
 */

/* fftw.h -- system-wide definitions */

#include <stdlib.h>
#include <stdio.h>
#include "fftw.h"

/* ----------- */

#ifdef __cplusplus
extern "C" {
#endif				/* __cplusplus */

/* our real numbers */
typedef double FFTW_REAL;

/*********************************************
 * Complex numbers and operations 
 *********************************************/
typedef struct {
     FFTW_REAL re, im;
} FFTW_COMPLEX;


/*********************************************
 *              Codelets
 *********************************************/
/*
 * There are two kinds of codelets:
 *
 * NO_TWIDDLE    computes the FFT of a certain size, operating
 *               out-of-place (i.e., take an input and produce a
 *               separate output)
 *
 * TWIDDLE       like no_twiddle, but operating in place.  Moreover,
 *               multiplies the input by twiddle factors.
 */

typedef void (notw_codelet) (const FFTW_COMPLEX *, FFTW_COMPLEX *, int, int);
typedef void (twiddle_codelet) (FFTW_COMPLEX *, const FFTW_COMPLEX *, int, int, int);
typedef void (generic_codelet) (FFTW_COMPLEX *, const FFTW_COMPLEX *, int, int, int, int);

/*********************************************
 *     Configurations
 *********************************************/
/*
 * A configuration is a database of all known codelets
 */

typedef struct {
     int size;			/* size of the problem */
     int signature;		/* unique codelet id */
     notw_codelet *codelet;	/*
				 * pointer to the codelet that solves
				 * the problem
				 */
} config_notw;

typedef struct {
     int size;			/* size of the problem */
     int signature;		/* unique codelet id */
     twiddle_codelet *codelet;
} config_twiddle;

/*****************************
 *        Plans
 *****************************/
/*
 * A plan is a sequence of reductions to compute a FFT of
 * a given size.  At each step, the FFT algorithm can:
 *
 * 1) apply a notw codelet, or
 * 2) recurse and apply a twiddle codelet, or
 * 3) apply the generic codelet.
 */

/* structure that contains twiddle factors */
typedef struct fftw_twiddle_struct {
     int n;
     int r;
     int m;
     FFTW_COMPLEX *twarray;
     struct fftw_twiddle_struct *next;
     int refcnt;
} fftw_twiddle;

/* structure that holds all the data needed for a given step */
typedef struct fftw_plan_node_struct {
     enum fftw_node_type type;

     union {
	  /* nodes of type FFTW_NOTW */
	  struct {
	       int size;
	       notw_codelet *codelet;
	  } notw;

	  /* nodes of type FFTW_TWIDDLE */
	  struct {
	       int size;
	       twiddle_codelet *codelet;
	       fftw_twiddle *tw;
	       struct fftw_plan_node_struct *recurse;
	  } twiddle;

	  /* nodes of type FFTW_GENERIC */
	  struct {
	       int size;
	       generic_codelet *codelet;
	       fftw_twiddle *tw;
	       struct fftw_plan_node_struct *recurse;
	  } generic;

     } nodeu;

     int refcnt;
} fftw_plan_node;

struct fftw_plan_struct {
     int n;
     fftw_direction dir;
     fftw_plan_node *root;

     double cost;
     int flags;

     enum fftw_node_type wisdom_type;
     int wisdom_signature;

     struct fftw_plan_struct *next;
     int refcnt;
};

/* a plan is just an array of instructions */
typedef struct fftw_plan_struct *fftw_plan;

fftw_plan qe_fftw_create_plan(int n, fftw_direction dir, int flags);
void qe_fftw_destroy_plan(fftw_plan plan);
void fftw(fftw_plan plan, int howmany, FFTW_COMPLEX *in, int istride, int idist, FFTW_COMPLEX *out, int ostride, int odist);


/*****************************
 *    N-dimensional code
 *****************************/
typedef struct {
     int is_in_place;		/* 1 if for in-place FFT's, 0 otherwise */
     int rank;			/* 
				 * the rank (number of dimensions) of the
				 * array to be FFT'ed
				 */
     int *n;			/* 
				 * the dimensions of the array to the
				 * FFT'ed 
				 */
     int *n_before;		/* 
				 * n_before[i] = product of n[j] for j < i 
				 */
     int *n_after;		/* n_after[i] = product of n[j] for j > i */
     fftw_plan *plans;		/* fftw plans for each dimension */
     FFTW_COMPLEX *work;	/* 
				 * work array for FFT when doing
				 * "in-place" FFT 
				 */
} fftwnd_aux_data;

typedef fftwnd_aux_data *fftwnd_plan;

/* Initializing the FFTWND Auxiliary Data */
fftwnd_plan qe_fftw2d_create_plan(int nx, int ny, fftw_direction dir, int flags);
fftwnd_plan qe_fftw3d_create_plan(int nx, int ny, int nz,fftw_direction dir, int flags);
fftwnd_plan qe_fftwnd_create_plan(int rank, const int *n, fftw_direction dir, int flags);
/* Freeing the FFTWND Auxiliary Data */
void qe_fftwnd_destroy_plan(fftwnd_plan plan);
/* Computing the N-Dimensional FFT */
void fftwnd(fftwnd_plan plan, int howmany, FFTW_COMPLEX *in, int istride, int idist, FFTW_COMPLEX *out, int ostride, int odist);

void fftw_naive(int n, FFTW_COMPLEX* in, FFTW_COMPLEX* out);
void fftwi_naive(int n, FFTW_COMPLEX* in, FFTW_COMPLEX* out);
void fftw_fprint_plan(FILE* f, fftw_plan plan);

/****************************************************************************/

#ifdef __cplusplus
}                               /* extern "C" */
#endif				/* __cplusplus */

#endif				/* FFTW_H */

