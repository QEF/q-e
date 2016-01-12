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

#include <stdio.h>
#include <stdlib.h>

#if defined(__QK_USER__)
#include <catamount/dclock.h>
#endif

#include "fftw.h"

/**************** import/export using file ***************/

static void file_emitter(char c, void *data)
{
     putc(c,(FILE *) data);
}

void fftw_export_wisdom_to_file(FILE *output_file)
{
     if (output_file)
	  fftw_export_wisdom(file_emitter,(void *) output_file);
}

static int file_get_input(void *data)
{
     return getc((FILE *) data);
}

fftw_status fftw_import_wisdom_from_file(FILE *input_file)
{
     if (!input_file)
	  return FFTW_FAILURE;
     return fftw_import_wisdom(file_get_input, (void *) input_file);
}

/*************** import/export using string **************/

static void emission_counter(char c, void *data)
{
     int *counter = (int *) data;

     ++*counter;
}

static void string_emitter(char c, void *data)
{
     char **output_string = (char **) data;

     *((*output_string)++) = c;
     **output_string = 0;
}

char *fftw_export_wisdom_to_string(void)
{
     int string_length = 0;
     char *s, *s2;

     fftw_export_wisdom(emission_counter, (void *) &string_length);

     s = fftw_malloc(sizeof(char) * (string_length + 1));
     if (!s)
	  return 0;
     s2 = s;

     fftw_export_wisdom(string_emitter, (void *) &s2);

     if (s + string_length != s2)
	  fftw_die("Unexpected output string length!");

     return s;
}

static int string_get_input(void *data)
{
     char **input_string = (char **) data;

     if (**input_string)
	  return *((*input_string)++);
     else
	  return 0;
}

fftw_status fftw_import_wisdom_from_string(const char *input_string)
{
     const char *s = input_string;

     if (!input_string)
	  return FFTW_FAILURE;
     return fftw_import_wisdom(string_get_input, (void *) &s);
}
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

/* config.c -- this file contains all the codelets the system knows about */

/* $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

#if defined FFTW_USING_CILK
#include <cilk.h>
#include <cilk-compat.h>
#endif

#include "fftw.h"

/* the signature is the same as the size, for now */
#define NOTW_CODELET(x)  { x, x, fftw_no_twiddle_##x }
#define NOTWI_CODELET(x)  { x, x, fftwi_no_twiddle_##x }

extern notw_codelet fftw_no_twiddle_1;
extern notw_codelet fftw_no_twiddle_2;
extern notw_codelet fftw_no_twiddle_3;
extern notw_codelet fftw_no_twiddle_4;
extern notw_codelet fftw_no_twiddle_5;
extern notw_codelet fftw_no_twiddle_6;
extern notw_codelet fftw_no_twiddle_7;
extern notw_codelet fftw_no_twiddle_8;
extern notw_codelet fftw_no_twiddle_9;
extern notw_codelet fftw_no_twiddle_10;
extern notw_codelet fftw_no_twiddle_11;
extern notw_codelet fftw_no_twiddle_12;
extern notw_codelet fftw_no_twiddle_13;
extern notw_codelet fftw_no_twiddle_14;
extern notw_codelet fftw_no_twiddle_15;
extern notw_codelet fftw_no_twiddle_16;
extern notw_codelet fftw_no_twiddle_32;
extern notw_codelet fftw_no_twiddle_64;

extern notw_codelet fftwi_no_twiddle_1;
extern notw_codelet fftwi_no_twiddle_2;
extern notw_codelet fftwi_no_twiddle_3;
extern notw_codelet fftwi_no_twiddle_4;
extern notw_codelet fftwi_no_twiddle_5;
extern notw_codelet fftwi_no_twiddle_6;
extern notw_codelet fftwi_no_twiddle_7;
extern notw_codelet fftwi_no_twiddle_8;
extern notw_codelet fftwi_no_twiddle_9;
extern notw_codelet fftwi_no_twiddle_10;
extern notw_codelet fftwi_no_twiddle_11;
extern notw_codelet fftwi_no_twiddle_12;
extern notw_codelet fftwi_no_twiddle_13;
extern notw_codelet fftwi_no_twiddle_14;
extern notw_codelet fftwi_no_twiddle_15;
extern notw_codelet fftwi_no_twiddle_16;
extern notw_codelet fftwi_no_twiddle_32;
extern notw_codelet fftwi_no_twiddle_64;

config_notw fftw_config_notw[] =
{
     NOTW_CODELET(1),
     NOTW_CODELET(2),
     NOTW_CODELET(3),
     NOTW_CODELET(4),
     NOTW_CODELET(5),
     NOTW_CODELET(6),
     NOTW_CODELET(7),
     NOTW_CODELET(8),
     NOTW_CODELET(9),
     NOTW_CODELET(10),
     NOTW_CODELET(11),
     NOTW_CODELET(12),
     NOTW_CODELET(13),
     NOTW_CODELET(14),
     NOTW_CODELET(15),
     NOTW_CODELET(16),
     NOTW_CODELET(32),
     NOTW_CODELET(64),
     {0, 0, (notw_codelet *) 0}
};

config_notw fftwi_config_notw[] =
{
     NOTWI_CODELET(1),
     NOTWI_CODELET(2),
     NOTWI_CODELET(3),
     NOTWI_CODELET(4),
     NOTWI_CODELET(5),
     NOTWI_CODELET(6),
     NOTWI_CODELET(7),
     NOTWI_CODELET(8),
     NOTWI_CODELET(9),
     NOTWI_CODELET(10),
     NOTWI_CODELET(11),
     NOTWI_CODELET(12),
     NOTWI_CODELET(13),
     NOTWI_CODELET(14),
     NOTWI_CODELET(15),
     NOTWI_CODELET(16),
     NOTWI_CODELET(32),
     NOTWI_CODELET(64),
     {0, 0, (notw_codelet *) 0}
};

/* the signature is the same as the size, for now */
#define TWIDDLE_CODELET(x)  { x, x, fftw_twiddle_##x }
#define TWIDDLEI_CODELET(x)  { x, x, fftwi_twiddle_##x }

extern twiddle_codelet fftw_twiddle_2;
extern twiddle_codelet fftw_twiddle_3;
extern twiddle_codelet fftw_twiddle_4;
extern twiddle_codelet fftw_twiddle_5;
extern twiddle_codelet fftw_twiddle_6;
extern twiddle_codelet fftw_twiddle_7;
extern twiddle_codelet fftw_twiddle_8;
extern twiddle_codelet fftw_twiddle_9;
extern twiddle_codelet fftw_twiddle_10;
extern twiddle_codelet fftw_twiddle_16;
extern twiddle_codelet fftw_twiddle_32;
extern twiddle_codelet fftw_twiddle_64;

extern twiddle_codelet fftwi_twiddle_2;
extern twiddle_codelet fftwi_twiddle_3;
extern twiddle_codelet fftwi_twiddle_4;
extern twiddle_codelet fftwi_twiddle_5;
extern twiddle_codelet fftwi_twiddle_6;
extern twiddle_codelet fftwi_twiddle_7;
extern twiddle_codelet fftwi_twiddle_8;
extern twiddle_codelet fftwi_twiddle_9;
extern twiddle_codelet fftwi_twiddle_10;
extern twiddle_codelet fftwi_twiddle_16;
extern twiddle_codelet fftwi_twiddle_32;
extern twiddle_codelet fftwi_twiddle_64;

config_twiddle fftw_config_twiddle[] =
{
     TWIDDLE_CODELET(2),
     TWIDDLE_CODELET(3),
     TWIDDLE_CODELET(4),
     TWIDDLE_CODELET(5),
     TWIDDLE_CODELET(6),
     TWIDDLE_CODELET(7),
     TWIDDLE_CODELET(8),
     TWIDDLE_CODELET(9),
     TWIDDLE_CODELET(10),
     TWIDDLE_CODELET(16),
     TWIDDLE_CODELET(32),
     TWIDDLE_CODELET(64),
     {0, 0, (twiddle_codelet *) 0}
};

config_twiddle fftwi_config_twiddle[] =
{
     TWIDDLEI_CODELET(2),
     TWIDDLEI_CODELET(3),
     TWIDDLEI_CODELET(4),
     TWIDDLEI_CODELET(5),
     TWIDDLEI_CODELET(6),
     TWIDDLEI_CODELET(7),
     TWIDDLEI_CODELET(8),
     TWIDDLEI_CODELET(9),
     TWIDDLEI_CODELET(10),
     TWIDDLEI_CODELET(16),
     TWIDDLEI_CODELET(32),
     TWIDDLEI_CODELET(64),
     {0, 0, (twiddle_codelet *) 0}
};
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

/*
 * executor.c -- execute the fft
 */

/* $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */
#include "fftw.h"
#include <stdio.h>
#include <stdlib.h>

/*
 * This function is called in other files, so we cannot declare
 * it as static. 
 */

void fftw_strided_copy(int n, FFTW_COMPLEX *in, int ostride,
		       FFTW_COMPLEX *out)
{
     int i;
     FFTW_REAL r0, r1, i0, i1;
     FFTW_REAL r2, r3, i2, i3;

     i = 0;
     if (n & 3)
	  for (; i < (n & 3); ++i) {
	       out[i * ostride] = in[i];
	  }
     for (; i < n; i += 4) {
	  r0 = c_re(in[i]);
	  i0 = c_im(in[i]);
	  r1 = c_re(in[i + 1]);
	  i1 = c_im(in[i + 1]);
	  r2 = c_re(in[i + 2]);
	  i2 = c_im(in[i + 2]);
	  r3 = c_re(in[i + 3]);
	  i3 = c_im(in[i + 3]);
	  c_re(out[i * ostride]) = r0;
	  c_im(out[i * ostride]) = i0;
	  c_re(out[(i + 1) * ostride]) = r1;
	  c_im(out[(i + 1) * ostride]) = i1;
	  c_re(out[(i + 2) * ostride]) = r2;
	  c_im(out[(i + 2) * ostride]) = i2;
	  c_re(out[(i + 3) * ostride]) = r3;
	  c_im(out[(i + 3) * ostride]) = i3;
     }
}

/*
 * Do *not* declare simple executor as static--we need to call it
 * from executor_cilk.cilk...also, preface its name with "fftw_"
 * to avoid any possible name collisions. 
 */
void fftw_executor_simple(int n, const FFTW_COMPLEX *in,
			  FFTW_COMPLEX *out,
			  fftw_plan_node *p,
			  int istride,
			  int ostride)
{
     switch (p->type) {
	 case FFTW_NOTW:
	      (p->nodeu.notw.codelet) (in, out, istride, ostride);
	      break;

	 case FFTW_TWIDDLE:
	      {
		   int r = p->nodeu.twiddle.size;
		   int m = n / r;
		   int i;
		   twiddle_codelet *codelet;
		   FFTW_COMPLEX *W;

		   for (i = 0; i < r; ++i) {
			fftw_executor_simple(m, in + i * istride,
					     out + i * (m * ostride),
					     p->nodeu.twiddle.recurse,
					     istride * r, ostride);
		   }

		   codelet = p->nodeu.twiddle.codelet;
		   W = p->nodeu.twiddle.tw->twarray;
		   codelet(out, W, m * ostride, m, ostride);

		   break;
	      }

	 case FFTW_GENERIC:
	      {
		   int r = p->nodeu.generic.size;
		   int m = n / r;
		   int i;
		   generic_codelet *codelet;
		   FFTW_COMPLEX *W;

		   for (i = 0; i < r; ++i) {
			fftw_executor_simple(m, in + i * istride,
					     out + i * (m * ostride),
					     p->nodeu.generic.recurse,
					     istride * r, ostride);
		   }

		   codelet = p->nodeu.generic.codelet;
		   W = p->nodeu.generic.tw->twarray;
		   codelet(out, W, m, r, n, ostride);

		   break;
	      }

	 default:
	      fftw_die("BUG in executor: illegal plan\n");
	      break;
     }
}

static void executor_simple_inplace(int n, FFTW_COMPLEX *in,
				    FFTW_COMPLEX *out,
				    fftw_plan_node *p,
				    int istride)
{
     switch (p->type) {
	 case FFTW_NOTW:
	      (p->nodeu.notw.codelet) (in, in, istride, istride);
	      break;

	 default:
	      {
		   FFTW_COMPLEX *tmp;

		   if (out)
			tmp = out;
		   else
			tmp = (FFTW_COMPLEX *)
			    fftw_malloc(n * sizeof(FFTW_COMPLEX));

		   fftw_executor_simple(n, in, tmp, p, istride, 1);
		   fftw_strided_copy(n, tmp, istride, in);

		   if (!out)
			fftw_free(tmp);
	      }
     }
}

static void executor_many(int n, const FFTW_COMPLEX *in,
			  FFTW_COMPLEX *out,
			  fftw_plan_node *p,
			  int istride,
			  int ostride,
			  int howmany, int idist, int odist)
{
     switch (p->type) {
	 case FFTW_NOTW:
	      {
		   int s;
		   notw_codelet *codelet = p->nodeu.notw.codelet;
		   for (s = 0; s < howmany; ++s)
			codelet(in + s * idist,
				out + s * odist,
				istride, ostride);
		   break;
	      }

	 default:
	      {
		   int s;
		   for (s = 0; s < howmany; ++s) {
			fftw_executor_simple(n, in + s * idist,
					     out + s * odist,
					     p, istride, ostride);
		   }
	      }
     }
}

static void executor_many_inplace(int n, FFTW_COMPLEX *in,
				  FFTW_COMPLEX *out,
				  fftw_plan_node *p,
				  int istride,
				  int howmany, int idist)
{
     switch (p->type) {
	 case FFTW_NOTW:
	      {
		   int s;
		   notw_codelet *codelet = p->nodeu.notw.codelet;
		   for (s = 0; s < howmany; ++s)
			codelet(in + s * idist,
				in + s * idist,
				istride, istride);
		   break;
	      }

	 default:
	      {
		   int s;
		   FFTW_COMPLEX *tmp;
		   if (out)
			tmp = out;
		   else
			tmp = (FFTW_COMPLEX *)
			    fftw_malloc(n * sizeof(FFTW_COMPLEX));

		   for (s = 0; s < howmany; ++s) {
			fftw_executor_simple(n,
					     in + s * idist,
					     tmp,
					     p, istride, 1);
			fftw_strided_copy(n, tmp, istride, in + s * idist);
		   }

		   if (!out)
			fftw_free(tmp);
	      }
     }
}

/* user interface */
void fftw(fftw_plan plan, int howmany, FFTW_COMPLEX *in, int istride,
	  int idist, FFTW_COMPLEX *out, int ostride, int odist)
{
     int n = plan->n;

     if (plan->flags & FFTW_IN_PLACE) {
	  if (howmany == 1) {
	       executor_simple_inplace(n, in, out, plan->root, istride);
	  } else {
	       executor_many_inplace(n, in, out, plan->root, istride, howmany,
				     idist);
	  }
     } else {
	  if (howmany == 1) {
	       fftw_executor_simple(n, in, out, plan->root, istride, ostride);
	  } else {
	       executor_many(n, in, out, plan->root, istride, ostride,
			     howmany, idist, odist);
	  }
     }
}

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

/* $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

#include <stdlib.h>

#include "fftw.h"

/* Prototypes for functions used internally in this file: */

static void fftw2d_out_of_place_aux(fftwnd_plan p, int howmany,
				FFTW_COMPLEX *in, int istride, int idist,
			      FFTW_COMPLEX *out, int ostride, int odist);
static void fftw3d_out_of_place_aux(fftwnd_plan p, int howmany,
				FFTW_COMPLEX *in, int istride, int idist,
			      FFTW_COMPLEX *out, int ostride, int odist);
static void fftwnd_out_of_place_aux(fftwnd_plan p, int howmany,
				FFTW_COMPLEX *in, int istride, int idist,
			      FFTW_COMPLEX *out, int ostride, int odist);

static void fftw2d_in_place_aux(fftwnd_plan p, int howmany,
			   FFTW_COMPLEX *in_out, int istride, int idist);
static void fftw3d_in_place_aux(fftwnd_plan p, int howmany,
			   FFTW_COMPLEX *in_out, int istride, int idist);
static void fftwnd_in_place_aux(fftwnd_plan p, int howmany,
			   FFTW_COMPLEX *in_out, int istride, int idist);

/*********** Initializing the FFTWND Auxiliary Data **********/

fftwnd_plan fftw2d_create_plan(int nx, int ny, fftw_direction dir, int flags)
{
     int n[2];

     n[0] = nx;
     n[1] = ny;

     return fftwnd_create_plan(2, n, dir, flags);
}

fftwnd_plan fftw3d_create_plan(int nx, int ny, int nz, fftw_direction dir,
			       int flags)
{
     int n[3];

     n[0] = nx;
     n[1] = ny;
     n[2] = nz;

     return fftwnd_create_plan(3, n, dir, flags);
}

fftwnd_plan fftwnd_create_plan(int rank, const int *n, 
			       fftw_direction dir, int flags)
{
     int i, j, max_dim = 0;
     fftwnd_plan p;
     int cur_flags;

     if (rank < 0)
	  return 0;

     for (i = 0; i < rank; ++i)
	  if (n[i] <= 0)
	       return 0;

     p = (fftwnd_plan) fftw_malloc(sizeof(fftwnd_aux_data));
     p->n = 0;
     p->n_before = 0;
     p->n_after = 0;
     p->plans = 0;
     p->work = 0;

     p->rank = rank;
     p->is_in_place = flags & FFTW_IN_PLACE;

     if (rank == 0)
	  return 0;

     p->n = (int *) fftw_malloc(sizeof(int) * rank);
     p->n_before = (int *) fftw_malloc(sizeof(int) * rank);
     p->n_after = (int *) fftw_malloc(sizeof(int) * rank);
     p->plans = (fftw_plan *) fftw_malloc(rank * sizeof(fftw_plan));
     p->n_before[0] = 1;
     p->n_after[rank - 1] = 1;

     for (i = 0; i < rank; ++i) {
	  p->n[i] = n[i];

	  if (i) {
	       p->n_before[i] = p->n_before[i - 1] * n[i - 1];
	       p->n_after[rank - 1 - i] = p->n_after[rank - i] * n[rank - i];
	  }
	  if (i < rank - 1 || (flags & FFTW_IN_PLACE)) {
	       /* fft's except the last dimension are always in-place */
	       cur_flags = flags | FFTW_IN_PLACE;
	       for (j = i - 1; j >= 0 && n[i] != n[j]; --j);

	       if (n[i] > max_dim)
		    max_dim = n[i];
	  } else {
	       cur_flags = flags;
	       /* we must create a separate plan for the last dimension */
	       j = -1;
	  }

	  if (j >= 0) {
	       /* 
	        * If a plan already exists for this size
	        * array, reuse it: 
	        */
	       p->plans[i] = p->plans[j];
	  } else {
	       /* generate a new plan: */
	       p->plans[i] = fftw_create_plan(n[i], dir, cur_flags);
	       if (!p->plans[i]) {
		    fftwnd_destroy_plan(p);
		    return 0;
	       }
	  }
     }

     /* Create work array for in-place FFTs: */
     if (max_dim > 0)
	  p->work = (FFTW_COMPLEX *)
	      fftw_malloc(sizeof(FFTW_COMPLEX) * max_dim);

     return p;
}

/************* Freeing the FFTWND Auxiliary Data *************/

void fftwnd_destroy_plan(fftwnd_plan plan)
{
     if (plan) {
	  if (plan->plans) {
	       int i, j;

	       for (i = 0; i < plan->rank; ++i) {
		    for (j = i - 1;
			 j >= 0 && plan->plans[i] != plan->plans[j];
			 --j);
		    if (j < 0 && plan->plans[i])
			 fftw_destroy_plan(plan->plans[i]);
	       }
	       fftw_free(plan->plans);
	  }
	  if (plan->n)
	       fftw_free(plan->n);

	  if (plan->n_before)
	       fftw_free(plan->n_before);

	  if (plan->n_after)
	       fftw_free(plan->n_after);

	  if (plan->work)
	       fftw_free(plan->work);

	  fftw_free(plan);
     }
}

/************** Computing the N-Dimensional FFT **************/

void fftwnd(fftwnd_plan plan, int howmany,
	    FFTW_COMPLEX *in, int istride, int idist,
	    FFTW_COMPLEX *out, int ostride, int odist)
{
     if (plan->is_in_place)	/* fft is in-place */
	  switch (plan->rank) {
	      case 0:
		   break;
	      case 1:
		   fftw(plan->plans[0], howmany, in, istride, idist,
			plan->work, 1, 0);
		   break;
	      case 2:
		   fftw2d_in_place_aux(plan, howmany, in, istride, idist);
		   break;
	      case 3:
		   fftw3d_in_place_aux(plan, howmany, in, istride, idist);
		   break;
	      default:
		   fftwnd_in_place_aux(plan, howmany, in, istride, idist);
     } else {
	  if (in == out || out == 0)
	       fftw_die("Illegal attempt to perform in-place FFT!\n");
	  switch (plan->rank) {
	      case 0:
		   break;
	      case 1:
		   fftw(plan->plans[0], howmany, in, istride, idist,
			out, ostride, odist);
		   break;
	      case 2:
		   fftw2d_out_of_place_aux(plan, howmany, in, istride,
					   idist, out, ostride, odist);
		   break;
	      case 3:
		   fftw3d_out_of_place_aux(plan, howmany, in, istride,
					   idist, out, ostride, odist);
		   break;
	      default:
		   fftwnd_out_of_place_aux(plan, howmany, in, istride,
					   idist, out, ostride, odist);
	  }
     }
}

static void fftw2d_out_of_place_aux(fftwnd_plan p, int howmany,
				FFTW_COMPLEX *in, int istride, int idist,
			       FFTW_COMPLEX *out, int ostride, int odist)
{
     int fft_iter;
     fftw_plan p0, p1;
     int n0, n1;

     p0 = p->plans[0];
     p1 = p->plans[1];
     n0 = p->n[0];
     n1 = p->n[1];

     for (fft_iter = 0; fft_iter < howmany; ++fft_iter) {
	  /* FFT y dimension (out-of-place): */
	  fftw(p1, n0,
	       in + fft_iter * idist, istride, n1 * istride,
	       out + fft_iter * odist, ostride, n1 * ostride);
	  /* FFT x dimension (in-place): */
	  fftw(p0, n1,
	       out + fft_iter * odist, n1 * ostride, ostride,
	       p->work, 1, 1);
     }
}

static void fftw3d_out_of_place_aux(fftwnd_plan p, int howmany,
				FFTW_COMPLEX *in, int istride, int idist,
			       FFTW_COMPLEX *out, int ostride, int odist)
{
     int fft_iter;
     int i;
     fftw_plan p0, p1, p2;
     int n0, n1, n2;

     p0 = p->plans[0];
     p1 = p->plans[1];
     p2 = p->plans[2];
     n0 = p->n[0];
     n1 = p->n[1];
     n2 = p->n[2];

     for (fft_iter = 0; fft_iter < howmany; ++fft_iter) {
	  /* FFT z dimension (out-of-place): */
	  fftw(p2, n0 * n1,
	       in + fft_iter * idist, istride, n2 * istride,
	       out + fft_iter * odist, ostride, n2 * ostride);
	  /* FFT y dimension (in-place): */
	  for (i = 0; i < n0; ++i)
	       fftw(p1, n2,
		    out + fft_iter * odist + i * n1 * n2 * ostride,
		    n2 * ostride, ostride, p->work, 1, 0);
	  /* FFT x dimension (in-place): */
	  fftw(p0, n1 * n2,
	       out + fft_iter * odist, n1 * n2 * ostride, ostride,
	       p->work, 1, 0);
     }
}

static void fftwnd_out_of_place_aux(fftwnd_plan p, int howmany,
				FFTW_COMPLEX *in, int istride, int idist,
			       FFTW_COMPLEX *out, int ostride, int odist)
{
     int fft_iter;
     int j, i;

     /* Do FFT for rank > 3: */

     for (fft_iter = 0; fft_iter < howmany; ++fft_iter) {
	  /* do last dimension (out-of-place): */
	  fftw(p->plans[p->rank - 1], p->n_before[p->rank - 1],
	     in + fft_iter * idist, istride, p->n[p->rank - 1] * istride,
	   out + fft_iter * odist, ostride, p->n[p->rank - 1] * ostride);

	  /* do first dimension (in-place): */
	  fftw(p->plans[0], p->n_after[0],
	       out + fft_iter * odist, p->n_after[0] * ostride, ostride,
	       p->work, 1, 0);

	  /* do other dimensions (in-place): */
	  for (j = 1; j < p->rank - 1; ++j)
	       for (i = 0; i < p->n_before[j]; ++i)
		    fftw(p->plans[j], p->n_after[j],
			 out + fft_iter * odist + i * ostride * p->n[j] *
			 p->n_after[j], p->n_after[j] * ostride,
			 ostride, p->work, 1, 0);
     }
}

static void fftw2d_in_place_aux(fftwnd_plan p, int howmany,
			    FFTW_COMPLEX *in_out, int istride, int idist)
{
     int fft_iter;
     fftw_plan p0, p1;
     int n0, n1;

     p0 = p->plans[0];
     p1 = p->plans[1];
     n0 = p->n[0];
     n1 = p->n[1];

     for (fft_iter = 0; fft_iter < howmany; ++fft_iter) {
	  /* FFT y dimension: */
	  fftw(p1, n0,
	       in_out + fft_iter * idist, istride, istride * n1,
	       p->work, 1, 0);
	  /* FFT x dimension: */
	  fftw(p0, n1,
	       in_out + fft_iter * idist, istride * n1, istride,
	       p->work, 1, 0);
     }
}

static void fftw3d_in_place_aux(fftwnd_plan p, int howmany,
			    FFTW_COMPLEX *in_out, int istride, int idist)
{
     int i;
     int fft_iter;
     fftw_plan p0, p1, p2;
     int n0, n1, n2;

     p0 = p->plans[0];
     p1 = p->plans[1];
     p2 = p->plans[2];
     n0 = p->n[0];
     n1 = p->n[1];
     n2 = p->n[2];

     for (fft_iter = 0; fft_iter < howmany; ++fft_iter) {
	  /* FFT z dimension: */
	  fftw(p2, n0 * n1,
	       in_out + fft_iter * idist, istride, n2 * istride,
	       p->work, 1, 0);
	  /* FFT y dimension: */
	  for (i = 0; i < n0; ++i)
	       fftw(p1, n2,
		    in_out + fft_iter * idist + i * n1 *
		    n2 * istride, n2 * istride, istride, p->work, 1, 0);
	  /* FFT x dimension: */
	  fftw(p0, n1 * n2,
	       in_out + fft_iter * idist, n1 * n2 * istride, istride,
	       p->work, 1, 0);
     }
}

static void fftwnd_in_place_aux(fftwnd_plan p, int howmany,
			    FFTW_COMPLEX *in_out, int istride, int idist)
/* Do FFT for rank > 3: */
{
     int fft_iter;
     int j, i;

     for (fft_iter = 0; fft_iter < howmany; ++fft_iter) {
	  /* do last dimension: */
	  fftw(p->plans[p->rank - 1], p->n_before[p->rank - 1],
	  in_out + fft_iter * idist, istride, p->n[p->rank - 1] * istride,
	       p->work, 1, 0);

	  /* do first dimension: */
	  fftw(p->plans[0], p->n_after[0],
	     in_out + fft_iter * idist, p->n_after[0] * istride, istride,
	       p->work, 1, 0);

	  /* do other dimensions: */
	  for (j = 1; j < p->rank - 1; ++j)
	       for (i = 0; i < p->n_before[j]; ++i)
		    fftw(p->plans[j], p->n_after[j],
		      in_out + fft_iter * idist + i * istride * p->n[j] *
			 p->n_after[j], p->n_after[j] * istride, istride,
			 p->work, 1, 0);
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 0 FP additions and 0 FP multiplications */

void fftw_no_twiddle_1(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     tre0_0_0 = c_re(in[0]);
     tim0_0_0 = c_im(in[0]);
     c_re(out[0]) = tre0_0_0;
     c_im(out[0]) = tim0_0_0;
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 108 FP additions and 32 FP multiplications */

void fftw_no_twiddle_10(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_0_1;
     FFTW_REAL tim0_0_1;
     FFTW_REAL tre0_0_2;
     FFTW_REAL tim0_0_2;
     FFTW_REAL tre0_0_3;
     FFTW_REAL tim0_0_3;
     FFTW_REAL tre0_0_4;
     FFTW_REAL tim0_0_4;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_1_1;
     FFTW_REAL tim0_1_1;
     FFTW_REAL tre0_1_2;
     FFTW_REAL tim0_1_2;
     FFTW_REAL tre0_1_3;
     FFTW_REAL tim0_1_3;
     FFTW_REAL tre0_1_4;
     FFTW_REAL tim0_1_4;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[0]);
	  tim1_0_0 = c_im(in[0]);
	  tre1_1_0 = c_re(in[5 * istride]);
	  tim1_1_0 = c_im(in[5 * istride]);
	  tre0_0_0 = tre1_0_0 + tre1_1_0;
	  tim0_0_0 = tim1_0_0 + tim1_1_0;
	  tre0_1_0 = tre1_0_0 - tre1_1_0;
	  tim0_1_0 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[2 * istride]);
	  tim1_0_0 = c_im(in[2 * istride]);
	  tre1_1_0 = c_re(in[7 * istride]);
	  tim1_1_0 = c_im(in[7 * istride]);
	  tre0_0_1 = tre1_0_0 + tre1_1_0;
	  tim0_0_1 = tim1_0_0 + tim1_1_0;
	  tre0_1_1 = tre1_0_0 - tre1_1_0;
	  tim0_1_1 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[4 * istride]);
	  tim1_0_0 = c_im(in[4 * istride]);
	  tre1_1_0 = c_re(in[9 * istride]);
	  tim1_1_0 = c_im(in[9 * istride]);
	  tre0_0_2 = tre1_0_0 + tre1_1_0;
	  tim0_0_2 = tim1_0_0 + tim1_1_0;
	  tre0_1_2 = tre1_0_0 - tre1_1_0;
	  tim0_1_2 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[6 * istride]);
	  tim1_0_0 = c_im(in[6 * istride]);
	  tre1_1_0 = c_re(in[istride]);
	  tim1_1_0 = c_im(in[istride]);
	  tre0_0_3 = tre1_0_0 + tre1_1_0;
	  tim0_0_3 = tim1_0_0 + tim1_1_0;
	  tre0_1_3 = tre1_0_0 - tre1_1_0;
	  tim0_1_3 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[8 * istride]);
	  tim1_0_0 = c_im(in[8 * istride]);
	  tre1_1_0 = c_re(in[3 * istride]);
	  tim1_1_0 = c_im(in[3 * istride]);
	  tre0_0_4 = tre1_0_0 + tre1_1_0;
	  tim0_0_4 = tim1_0_0 + tim1_1_0;
	  tre0_1_4 = tre1_0_0 - tre1_1_0;
	  tim0_1_4 = tim1_0_0 - tim1_1_0;
     }
     c_re(out[0]) = tre0_0_0 + tre0_0_1 + tre0_0_2 + tre0_0_3 + tre0_0_4;
     c_im(out[0]) = tim0_0_0 + tim0_0_1 + tim0_0_2 + tim0_0_3 + tim0_0_4;
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_0_1 + tre0_0_4)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_0_2 + tre0_0_3));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tim0_0_1 - tim0_0_4)) + (((FFTW_REAL) FFTW_K587785252) * (tim0_0_2 - tim0_0_3));
	  c_re(out[6 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[4 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_0_1 + tim0_0_4)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_0_2 + tim0_0_3));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tre0_0_4 - tre0_0_1)) + (((FFTW_REAL) FFTW_K587785252) * (tre0_0_3 - tre0_0_2));
	  c_im(out[6 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[4 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_0_2 + tre0_0_3)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_0_1 + tre0_0_4));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tim0_0_1 - tim0_0_4)) + (((FFTW_REAL) FFTW_K951056516) * (tim0_0_3 - tim0_0_2));
	  c_re(out[2 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[8 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_0_2 + tim0_0_3)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_0_1 + tim0_0_4));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tre0_0_4 - tre0_0_1)) + (((FFTW_REAL) FFTW_K951056516) * (tre0_0_2 - tre0_0_3));
	  c_im(out[2 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[8 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     c_re(out[5 * ostride]) = tre0_1_0 + tre0_1_1 + tre0_1_2 + tre0_1_3 + tre0_1_4;
     c_im(out[5 * ostride]) = tim0_1_0 + tim0_1_1 + tim0_1_2 + tim0_1_3 + tim0_1_4;
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_1_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_1_1 + tre0_1_4)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_1_2 + tre0_1_3));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tim0_1_1 - tim0_1_4)) + (((FFTW_REAL) FFTW_K587785252) * (tim0_1_2 - tim0_1_3));
	  c_re(out[ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[9 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_1_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_1_1 + tim0_1_4)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_1_2 + tim0_1_3));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tre0_1_4 - tre0_1_1)) + (((FFTW_REAL) FFTW_K587785252) * (tre0_1_3 - tre0_1_2));
	  c_im(out[ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[9 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_1_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_1_2 + tre0_1_3)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_1_1 + tre0_1_4));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tim0_1_1 - tim0_1_4)) + (((FFTW_REAL) FFTW_K951056516) * (tim0_1_3 - tim0_1_2));
	  c_re(out[7 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[3 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_1_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_1_2 + tim0_1_3)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_1_1 + tim0_1_4));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tre0_1_4 - tre0_1_1)) + (((FFTW_REAL) FFTW_K951056516) * (tre0_1_2 - tre0_1_3));
	  c_im(out[7 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[3 * ostride]) = tim2_0_0 - tim2_1_0;
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 230 FP additions and 100 FP multiplications */

void fftw_no_twiddle_11(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_2_0;
     FFTW_REAL tim0_2_0;
     FFTW_REAL tre0_3_0;
     FFTW_REAL tim0_3_0;
     FFTW_REAL tre0_4_0;
     FFTW_REAL tim0_4_0;
     FFTW_REAL tre0_5_0;
     FFTW_REAL tim0_5_0;
     FFTW_REAL tre0_6_0;
     FFTW_REAL tim0_6_0;
     FFTW_REAL tre0_7_0;
     FFTW_REAL tim0_7_0;
     FFTW_REAL tre0_8_0;
     FFTW_REAL tim0_8_0;
     FFTW_REAL tre0_9_0;
     FFTW_REAL tim0_9_0;
     FFTW_REAL tre0_10_0;
     FFTW_REAL tim0_10_0;
     tre0_0_0 = c_re(in[0]);
     tim0_0_0 = c_im(in[0]);
     tre0_1_0 = c_re(in[istride]);
     tim0_1_0 = c_im(in[istride]);
     tre0_2_0 = c_re(in[2 * istride]);
     tim0_2_0 = c_im(in[2 * istride]);
     tre0_3_0 = c_re(in[3 * istride]);
     tim0_3_0 = c_im(in[3 * istride]);
     tre0_4_0 = c_re(in[4 * istride]);
     tim0_4_0 = c_im(in[4 * istride]);
     tre0_5_0 = c_re(in[5 * istride]);
     tim0_5_0 = c_im(in[5 * istride]);
     tre0_6_0 = c_re(in[6 * istride]);
     tim0_6_0 = c_im(in[6 * istride]);
     tre0_7_0 = c_re(in[7 * istride]);
     tim0_7_0 = c_im(in[7 * istride]);
     tre0_8_0 = c_re(in[8 * istride]);
     tim0_8_0 = c_im(in[8 * istride]);
     tre0_9_0 = c_re(in[9 * istride]);
     tim0_9_0 = c_im(in[9 * istride]);
     tre0_10_0 = c_re(in[10 * istride]);
     tim0_10_0 = c_im(in[10 * istride]);
     c_re(out[0]) = tre0_0_0 + tre0_1_0 + tre0_2_0 + tre0_3_0 + tre0_4_0 + tre0_5_0 + tre0_6_0 + tre0_7_0 + tre0_8_0 + tre0_9_0 + tre0_10_0;
     c_im(out[0]) = tim0_0_0 + tim0_1_0 + tim0_2_0 + tim0_3_0 + tim0_4_0 + tim0_5_0 + tim0_6_0 + tim0_7_0 + tim0_8_0 + tim0_9_0 + tim0_10_0;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K841253532) * (tre0_1_0 + tre0_10_0)) + (((FFTW_REAL) FFTW_K415415013) * (tre0_2_0 + tre0_9_0)) - (((FFTW_REAL) FFTW_K959492973) * (tre0_5_0 + tre0_6_0)) - (((FFTW_REAL) FFTW_K654860733) * (tre0_4_0 + tre0_7_0)) - (((FFTW_REAL) FFTW_K142314838) * (tre0_3_0 + tre0_8_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K540640817) * (tim0_1_0 - tim0_10_0)) + (((FFTW_REAL) FFTW_K909631995) * (tim0_2_0 - tim0_9_0)) + (((FFTW_REAL) FFTW_K989821441) * (tim0_3_0 - tim0_8_0)) + (((FFTW_REAL) FFTW_K755749574) * (tim0_4_0 - tim0_7_0)) + (((FFTW_REAL) FFTW_K281732556) * (tim0_5_0 - tim0_6_0));
	  c_re(out[ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[10 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K841253532) * (tim0_1_0 + tim0_10_0)) + (((FFTW_REAL) FFTW_K415415013) * (tim0_2_0 + tim0_9_0)) - (((FFTW_REAL) FFTW_K959492973) * (tim0_5_0 + tim0_6_0)) - (((FFTW_REAL) FFTW_K654860733) * (tim0_4_0 + tim0_7_0)) - (((FFTW_REAL) FFTW_K142314838) * (tim0_3_0 + tim0_8_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K540640817) * (tre0_10_0 - tre0_1_0)) + (((FFTW_REAL) FFTW_K909631995) * (tre0_9_0 - tre0_2_0)) + (((FFTW_REAL) FFTW_K989821441) * (tre0_8_0 - tre0_3_0)) + (((FFTW_REAL) FFTW_K755749574) * (tre0_7_0 - tre0_4_0)) + (((FFTW_REAL) FFTW_K281732556) * (tre0_6_0 - tre0_5_0));
	  c_im(out[ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[10 * ostride]) = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K415415013) * (tre0_1_0 + tre0_10_0)) + (((FFTW_REAL) FFTW_K841253532) * (tre0_5_0 + tre0_6_0)) - (((FFTW_REAL) FFTW_K142314838) * (tre0_4_0 + tre0_7_0)) - (((FFTW_REAL) FFTW_K959492973) * (tre0_3_0 + tre0_8_0)) - (((FFTW_REAL) FFTW_K654860733) * (tre0_2_0 + tre0_9_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K909631995) * (tim0_1_0 - tim0_10_0)) + (((FFTW_REAL) FFTW_K755749574) * (tim0_2_0 - tim0_9_0)) + (((FFTW_REAL) FFTW_K281732556) * (tim0_8_0 - tim0_3_0)) + (((FFTW_REAL) FFTW_K989821441) * (tim0_7_0 - tim0_4_0)) + (((FFTW_REAL) FFTW_K540640817) * (tim0_6_0 - tim0_5_0));
	  c_re(out[2 * ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[9 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K415415013) * (tim0_1_0 + tim0_10_0)) + (((FFTW_REAL) FFTW_K841253532) * (tim0_5_0 + tim0_6_0)) - (((FFTW_REAL) FFTW_K142314838) * (tim0_4_0 + tim0_7_0)) - (((FFTW_REAL) FFTW_K959492973) * (tim0_3_0 + tim0_8_0)) - (((FFTW_REAL) FFTW_K654860733) * (tim0_2_0 + tim0_9_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K909631995) * (tre0_10_0 - tre0_1_0)) + (((FFTW_REAL) FFTW_K755749574) * (tre0_9_0 - tre0_2_0)) + (((FFTW_REAL) FFTW_K281732556) * (tre0_3_0 - tre0_8_0)) + (((FFTW_REAL) FFTW_K989821441) * (tre0_4_0 - tre0_7_0)) + (((FFTW_REAL) FFTW_K540640817) * (tre0_5_0 - tre0_6_0));
	  c_im(out[2 * ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[9 * ostride]) = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K415415013) * (tre0_3_0 + tre0_8_0)) + (((FFTW_REAL) FFTW_K841253532) * (tre0_4_0 + tre0_7_0)) - (((FFTW_REAL) FFTW_K654860733) * (tre0_5_0 + tre0_6_0)) - (((FFTW_REAL) FFTW_K959492973) * (tre0_2_0 + tre0_9_0)) - (((FFTW_REAL) FFTW_K142314838) * (tre0_1_0 + tre0_10_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K989821441) * (tim0_1_0 - tim0_10_0)) + (((FFTW_REAL) FFTW_K281732556) * (tim0_9_0 - tim0_2_0)) + (((FFTW_REAL) FFTW_K909631995) * (tim0_8_0 - tim0_3_0)) + (((FFTW_REAL) FFTW_K540640817) * (tim0_4_0 - tim0_7_0)) + (((FFTW_REAL) FFTW_K755749574) * (tim0_5_0 - tim0_6_0));
	  c_re(out[3 * ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[8 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K415415013) * (tim0_3_0 + tim0_8_0)) + (((FFTW_REAL) FFTW_K841253532) * (tim0_4_0 + tim0_7_0)) - (((FFTW_REAL) FFTW_K654860733) * (tim0_5_0 + tim0_6_0)) - (((FFTW_REAL) FFTW_K959492973) * (tim0_2_0 + tim0_9_0)) - (((FFTW_REAL) FFTW_K142314838) * (tim0_1_0 + tim0_10_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K989821441) * (tre0_10_0 - tre0_1_0)) + (((FFTW_REAL) FFTW_K281732556) * (tre0_2_0 - tre0_9_0)) + (((FFTW_REAL) FFTW_K909631995) * (tre0_3_0 - tre0_8_0)) + (((FFTW_REAL) FFTW_K540640817) * (tre0_7_0 - tre0_4_0)) + (((FFTW_REAL) FFTW_K755749574) * (tre0_6_0 - tre0_5_0));
	  c_im(out[3 * ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[8 * ostride]) = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K841253532) * (tre0_3_0 + tre0_8_0)) + (((FFTW_REAL) FFTW_K415415013) * (tre0_5_0 + tre0_6_0)) - (((FFTW_REAL) FFTW_K959492973) * (tre0_4_0 + tre0_7_0)) - (((FFTW_REAL) FFTW_K142314838) * (tre0_2_0 + tre0_9_0)) - (((FFTW_REAL) FFTW_K654860733) * (tre0_1_0 + tre0_10_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K755749574) * (tim0_1_0 - tim0_10_0)) + (((FFTW_REAL) FFTW_K989821441) * (tim0_9_0 - tim0_2_0)) + (((FFTW_REAL) FFTW_K540640817) * (tim0_3_0 - tim0_8_0)) + (((FFTW_REAL) FFTW_K281732556) * (tim0_4_0 - tim0_7_0)) + (((FFTW_REAL) FFTW_K909631995) * (tim0_6_0 - tim0_5_0));
	  c_re(out[4 * ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[7 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K841253532) * (tim0_3_0 + tim0_8_0)) + (((FFTW_REAL) FFTW_K415415013) * (tim0_5_0 + tim0_6_0)) - (((FFTW_REAL) FFTW_K959492973) * (tim0_4_0 + tim0_7_0)) - (((FFTW_REAL) FFTW_K142314838) * (tim0_2_0 + tim0_9_0)) - (((FFTW_REAL) FFTW_K654860733) * (tim0_1_0 + tim0_10_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K755749574) * (tre0_10_0 - tre0_1_0)) + (((FFTW_REAL) FFTW_K989821441) * (tre0_2_0 - tre0_9_0)) + (((FFTW_REAL) FFTW_K540640817) * (tre0_8_0 - tre0_3_0)) + (((FFTW_REAL) FFTW_K281732556) * (tre0_7_0 - tre0_4_0)) + (((FFTW_REAL) FFTW_K909631995) * (tre0_5_0 - tre0_6_0));
	  c_im(out[4 * ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[7 * ostride]) = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K841253532) * (tre0_2_0 + tre0_9_0)) + (((FFTW_REAL) FFTW_K415415013) * (tre0_4_0 + tre0_7_0)) - (((FFTW_REAL) FFTW_K142314838) * (tre0_5_0 + tre0_6_0)) - (((FFTW_REAL) FFTW_K654860733) * (tre0_3_0 + tre0_8_0)) - (((FFTW_REAL) FFTW_K959492973) * (tre0_1_0 + tre0_10_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K281732556) * (tim0_1_0 - tim0_10_0)) + (((FFTW_REAL) FFTW_K540640817) * (tim0_9_0 - tim0_2_0)) + (((FFTW_REAL) FFTW_K755749574) * (tim0_3_0 - tim0_8_0)) + (((FFTW_REAL) FFTW_K909631995) * (tim0_7_0 - tim0_4_0)) + (((FFTW_REAL) FFTW_K989821441) * (tim0_5_0 - tim0_6_0));
	  c_re(out[5 * ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[6 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K841253532) * (tim0_2_0 + tim0_9_0)) + (((FFTW_REAL) FFTW_K415415013) * (tim0_4_0 + tim0_7_0)) - (((FFTW_REAL) FFTW_K142314838) * (tim0_5_0 + tim0_6_0)) - (((FFTW_REAL) FFTW_K654860733) * (tim0_3_0 + tim0_8_0)) - (((FFTW_REAL) FFTW_K959492973) * (tim0_1_0 + tim0_10_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K281732556) * (tre0_10_0 - tre0_1_0)) + (((FFTW_REAL) FFTW_K540640817) * (tre0_2_0 - tre0_9_0)) + (((FFTW_REAL) FFTW_K755749574) * (tre0_8_0 - tre0_3_0)) + (((FFTW_REAL) FFTW_K909631995) * (tre0_4_0 - tre0_7_0)) + (((FFTW_REAL) FFTW_K989821441) * (tre0_6_0 - tre0_5_0));
	  c_im(out[5 * ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[6 * ostride]) = tim1_0_0 - tim1_1_0;
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 104 FP additions and 16 FP multiplications */

void fftw_no_twiddle_12(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_0_1;
     FFTW_REAL tim0_0_1;
     FFTW_REAL tre0_0_2;
     FFTW_REAL tim0_0_2;
     FFTW_REAL tre0_0_3;
     FFTW_REAL tim0_0_3;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_1_1;
     FFTW_REAL tim0_1_1;
     FFTW_REAL tre0_1_2;
     FFTW_REAL tim0_1_2;
     FFTW_REAL tre0_1_3;
     FFTW_REAL tim0_1_3;
     FFTW_REAL tre0_2_0;
     FFTW_REAL tim0_2_0;
     FFTW_REAL tre0_2_1;
     FFTW_REAL tim0_2_1;
     FFTW_REAL tre0_2_2;
     FFTW_REAL tim0_2_2;
     FFTW_REAL tre0_2_3;
     FFTW_REAL tim0_2_3;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_2_0;
	  FFTW_REAL tim1_2_0;
	  tre1_0_0 = c_re(in[0]);
	  tim1_0_0 = c_im(in[0]);
	  tre1_1_0 = c_re(in[4 * istride]);
	  tim1_1_0 = c_im(in[4 * istride]);
	  tre1_2_0 = c_re(in[8 * istride]);
	  tim1_2_0 = c_im(in[8 * istride]);
	  tre0_0_0 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	  tim0_0_0 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_1_0 - tim1_2_0);
	       tre0_1_0 = tre2_0_0 + tre2_1_0;
	       tre0_2_0 = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_2_0 - tre1_1_0);
	       tim0_1_0 = tim2_0_0 + tim2_1_0;
	       tim0_2_0 = tim2_0_0 - tim2_1_0;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_2_0;
	  FFTW_REAL tim1_2_0;
	  tre1_0_0 = c_re(in[3 * istride]);
	  tim1_0_0 = c_im(in[3 * istride]);
	  tre1_1_0 = c_re(in[7 * istride]);
	  tim1_1_0 = c_im(in[7 * istride]);
	  tre1_2_0 = c_re(in[11 * istride]);
	  tim1_2_0 = c_im(in[11 * istride]);
	  tre0_0_1 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	  tim0_0_1 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_1_0 - tim1_2_0);
	       tre0_1_1 = tre2_0_0 + tre2_1_0;
	       tre0_2_1 = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_2_0 - tre1_1_0);
	       tim0_1_1 = tim2_0_0 + tim2_1_0;
	       tim0_2_1 = tim2_0_0 - tim2_1_0;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_2_0;
	  FFTW_REAL tim1_2_0;
	  tre1_0_0 = c_re(in[6 * istride]);
	  tim1_0_0 = c_im(in[6 * istride]);
	  tre1_1_0 = c_re(in[10 * istride]);
	  tim1_1_0 = c_im(in[10 * istride]);
	  tre1_2_0 = c_re(in[2 * istride]);
	  tim1_2_0 = c_im(in[2 * istride]);
	  tre0_0_2 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	  tim0_0_2 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_1_0 - tim1_2_0);
	       tre0_1_2 = tre2_0_0 + tre2_1_0;
	       tre0_2_2 = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_2_0 - tre1_1_0);
	       tim0_1_2 = tim2_0_0 + tim2_1_0;
	       tim0_2_2 = tim2_0_0 - tim2_1_0;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_2_0;
	  FFTW_REAL tim1_2_0;
	  tre1_0_0 = c_re(in[9 * istride]);
	  tim1_0_0 = c_im(in[9 * istride]);
	  tre1_1_0 = c_re(in[istride]);
	  tim1_1_0 = c_im(in[istride]);
	  tre1_2_0 = c_re(in[5 * istride]);
	  tim1_2_0 = c_im(in[5 * istride]);
	  tre0_0_3 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	  tim0_0_3 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_1_0 - tim1_2_0);
	       tre0_1_3 = tre2_0_0 + tre2_1_0;
	       tre0_2_3 = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_2_0 - tre1_1_0);
	       tim0_1_3 = tim2_0_0 + tim2_1_0;
	       tim0_2_3 = tim2_0_0 - tim2_1_0;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  tre1_0_0 = tre0_0_0 + tre0_0_2;
	  tim1_0_0 = tim0_0_0 + tim0_0_2;
	  tre1_1_0 = tre0_0_0 - tre0_0_2;
	  tim1_1_0 = tim0_0_0 - tim0_0_2;
	  tre1_0_1 = tre0_0_1 + tre0_0_3;
	  tim1_0_1 = tim0_0_1 + tim0_0_3;
	  tre1_1_1 = tre0_0_1 - tre0_0_3;
	  tim1_1_1 = tim0_0_1 - tim0_0_3;
	  c_re(out[0]) = tre1_0_0 + tre1_0_1;
	  c_im(out[0]) = tim1_0_0 + tim1_0_1;
	  c_re(out[6 * ostride]) = tre1_0_0 - tre1_0_1;
	  c_im(out[6 * ostride]) = tim1_0_0 - tim1_0_1;
	  c_re(out[9 * ostride]) = tre1_1_0 + tim1_1_1;
	  c_im(out[9 * ostride]) = tim1_1_0 - tre1_1_1;
	  c_re(out[3 * ostride]) = tre1_1_0 - tim1_1_1;
	  c_im(out[3 * ostride]) = tim1_1_0 + tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  tre1_0_0 = tre0_1_0 + tre0_1_2;
	  tim1_0_0 = tim0_1_0 + tim0_1_2;
	  tre1_1_0 = tre0_1_0 - tre0_1_2;
	  tim1_1_0 = tim0_1_0 - tim0_1_2;
	  tre1_0_1 = tre0_1_1 + tre0_1_3;
	  tim1_0_1 = tim0_1_1 + tim0_1_3;
	  tre1_1_1 = tre0_1_1 - tre0_1_3;
	  tim1_1_1 = tim0_1_1 - tim0_1_3;
	  c_re(out[4 * ostride]) = tre1_0_0 + tre1_0_1;
	  c_im(out[4 * ostride]) = tim1_0_0 + tim1_0_1;
	  c_re(out[10 * ostride]) = tre1_0_0 - tre1_0_1;
	  c_im(out[10 * ostride]) = tim1_0_0 - tim1_0_1;
	  c_re(out[ostride]) = tre1_1_0 + tim1_1_1;
	  c_im(out[ostride]) = tim1_1_0 - tre1_1_1;
	  c_re(out[7 * ostride]) = tre1_1_0 - tim1_1_1;
	  c_im(out[7 * ostride]) = tim1_1_0 + tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  tre1_0_0 = tre0_2_0 + tre0_2_2;
	  tim1_0_0 = tim0_2_0 + tim0_2_2;
	  tre1_1_0 = tre0_2_0 - tre0_2_2;
	  tim1_1_0 = tim0_2_0 - tim0_2_2;
	  tre1_0_1 = tre0_2_1 + tre0_2_3;
	  tim1_0_1 = tim0_2_1 + tim0_2_3;
	  tre1_1_1 = tre0_2_1 - tre0_2_3;
	  tim1_1_1 = tim0_2_1 - tim0_2_3;
	  c_re(out[8 * ostride]) = tre1_0_0 + tre1_0_1;
	  c_im(out[8 * ostride]) = tim1_0_0 + tim1_0_1;
	  c_re(out[2 * ostride]) = tre1_0_0 - tre1_0_1;
	  c_im(out[2 * ostride]) = tim1_0_0 - tim1_0_1;
	  c_re(out[5 * ostride]) = tre1_1_0 + tim1_1_1;
	  c_im(out[5 * ostride]) = tim1_1_0 - tre1_1_1;
	  c_re(out[11 * ostride]) = tre1_1_0 - tim1_1_1;
	  c_im(out[11 * ostride]) = tim1_1_0 + tre1_1_1;
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 324 FP additions and 144 FP multiplications */

void fftw_no_twiddle_13(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_2_0;
     FFTW_REAL tim0_2_0;
     FFTW_REAL tre0_3_0;
     FFTW_REAL tim0_3_0;
     FFTW_REAL tre0_4_0;
     FFTW_REAL tim0_4_0;
     FFTW_REAL tre0_5_0;
     FFTW_REAL tim0_5_0;
     FFTW_REAL tre0_6_0;
     FFTW_REAL tim0_6_0;
     FFTW_REAL tre0_7_0;
     FFTW_REAL tim0_7_0;
     FFTW_REAL tre0_8_0;
     FFTW_REAL tim0_8_0;
     FFTW_REAL tre0_9_0;
     FFTW_REAL tim0_9_0;
     FFTW_REAL tre0_10_0;
     FFTW_REAL tim0_10_0;
     FFTW_REAL tre0_11_0;
     FFTW_REAL tim0_11_0;
     FFTW_REAL tre0_12_0;
     FFTW_REAL tim0_12_0;
     tre0_0_0 = c_re(in[0]);
     tim0_0_0 = c_im(in[0]);
     tre0_1_0 = c_re(in[istride]);
     tim0_1_0 = c_im(in[istride]);
     tre0_2_0 = c_re(in[2 * istride]);
     tim0_2_0 = c_im(in[2 * istride]);
     tre0_3_0 = c_re(in[3 * istride]);
     tim0_3_0 = c_im(in[3 * istride]);
     tre0_4_0 = c_re(in[4 * istride]);
     tim0_4_0 = c_im(in[4 * istride]);
     tre0_5_0 = c_re(in[5 * istride]);
     tim0_5_0 = c_im(in[5 * istride]);
     tre0_6_0 = c_re(in[6 * istride]);
     tim0_6_0 = c_im(in[6 * istride]);
     tre0_7_0 = c_re(in[7 * istride]);
     tim0_7_0 = c_im(in[7 * istride]);
     tre0_8_0 = c_re(in[8 * istride]);
     tim0_8_0 = c_im(in[8 * istride]);
     tre0_9_0 = c_re(in[9 * istride]);
     tim0_9_0 = c_im(in[9 * istride]);
     tre0_10_0 = c_re(in[10 * istride]);
     tim0_10_0 = c_im(in[10 * istride]);
     tre0_11_0 = c_re(in[11 * istride]);
     tim0_11_0 = c_im(in[11 * istride]);
     tre0_12_0 = c_re(in[12 * istride]);
     tim0_12_0 = c_im(in[12 * istride]);
     c_re(out[0]) = tre0_0_0 + tre0_1_0 + tre0_2_0 + tre0_3_0 + tre0_4_0 + tre0_5_0 + tre0_6_0 + tre0_7_0 + tre0_8_0 + tre0_9_0 + tre0_10_0 + tre0_11_0 + tre0_12_0;
     c_im(out[0]) = tim0_0_0 + tim0_1_0 + tim0_2_0 + tim0_3_0 + tim0_4_0 + tim0_5_0 + tim0_6_0 + tim0_7_0 + tim0_8_0 + tim0_9_0 + tim0_10_0 + tim0_11_0 + tim0_12_0;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K885456025) * (tre0_1_0 + tre0_12_0)) + (((FFTW_REAL) FFTW_K568064746) * (tre0_2_0 + tre0_11_0)) + (((FFTW_REAL) FFTW_K120536680) * (tre0_3_0 + tre0_10_0)) - (((FFTW_REAL) FFTW_K970941817) * (tre0_6_0 + tre0_7_0)) - (((FFTW_REAL) FFTW_K748510748) * (tre0_5_0 + tre0_8_0)) - (((FFTW_REAL) FFTW_K354604887) * (tre0_4_0 + tre0_9_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K464723172) * (tim0_1_0 - tim0_12_0)) + (((FFTW_REAL) FFTW_K822983865) * (tim0_2_0 - tim0_11_0)) + (((FFTW_REAL) FFTW_K992708874) * (tim0_3_0 - tim0_10_0)) + (((FFTW_REAL) FFTW_K935016242) * (tim0_4_0 - tim0_9_0)) + (((FFTW_REAL) FFTW_K663122658) * (tim0_5_0 - tim0_8_0)) + (((FFTW_REAL) FFTW_K239315664) * (tim0_6_0 - tim0_7_0));
	  c_re(out[ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[12 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K885456025) * (tim0_1_0 + tim0_12_0)) + (((FFTW_REAL) FFTW_K568064746) * (tim0_2_0 + tim0_11_0)) + (((FFTW_REAL) FFTW_K120536680) * (tim0_3_0 + tim0_10_0)) - (((FFTW_REAL) FFTW_K970941817) * (tim0_6_0 + tim0_7_0)) - (((FFTW_REAL) FFTW_K748510748) * (tim0_5_0 + tim0_8_0)) - (((FFTW_REAL) FFTW_K354604887) * (tim0_4_0 + tim0_9_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K464723172) * (tre0_12_0 - tre0_1_0)) + (((FFTW_REAL) FFTW_K822983865) * (tre0_11_0 - tre0_2_0)) + (((FFTW_REAL) FFTW_K992708874) * (tre0_10_0 - tre0_3_0)) + (((FFTW_REAL) FFTW_K935016242) * (tre0_9_0 - tre0_4_0)) + (((FFTW_REAL) FFTW_K663122658) * (tre0_8_0 - tre0_5_0)) + (((FFTW_REAL) FFTW_K239315664) * (tre0_7_0 - tre0_6_0));
	  c_im(out[ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[12 * ostride]) = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K568064746) * (tre0_1_0 + tre0_12_0)) + (((FFTW_REAL) FFTW_K120536680) * (tre0_5_0 + tre0_8_0)) + (((FFTW_REAL) FFTW_K885456025) * (tre0_6_0 + tre0_7_0)) - (((FFTW_REAL) FFTW_K748510748) * (tre0_4_0 + tre0_9_0)) - (((FFTW_REAL) FFTW_K970941817) * (tre0_3_0 + tre0_10_0)) - (((FFTW_REAL) FFTW_K354604887) * (tre0_2_0 + tre0_11_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K822983865) * (tim0_1_0 - tim0_12_0)) + (((FFTW_REAL) FFTW_K935016242) * (tim0_2_0 - tim0_11_0)) + (((FFTW_REAL) FFTW_K239315664) * (tim0_3_0 - tim0_10_0)) + (((FFTW_REAL) FFTW_K663122658) * (tim0_9_0 - tim0_4_0)) + (((FFTW_REAL) FFTW_K992708874) * (tim0_8_0 - tim0_5_0)) + (((FFTW_REAL) FFTW_K464723172) * (tim0_7_0 - tim0_6_0));
	  c_re(out[2 * ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[11 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K568064746) * (tim0_1_0 + tim0_12_0)) + (((FFTW_REAL) FFTW_K120536680) * (tim0_5_0 + tim0_8_0)) + (((FFTW_REAL) FFTW_K885456025) * (tim0_6_0 + tim0_7_0)) - (((FFTW_REAL) FFTW_K748510748) * (tim0_4_0 + tim0_9_0)) - (((FFTW_REAL) FFTW_K970941817) * (tim0_3_0 + tim0_10_0)) - (((FFTW_REAL) FFTW_K354604887) * (tim0_2_0 + tim0_11_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K822983865) * (tre0_12_0 - tre0_1_0)) + (((FFTW_REAL) FFTW_K935016242) * (tre0_11_0 - tre0_2_0)) + (((FFTW_REAL) FFTW_K239315664) * (tre0_10_0 - tre0_3_0)) + (((FFTW_REAL) FFTW_K663122658) * (tre0_4_0 - tre0_9_0)) + (((FFTW_REAL) FFTW_K992708874) * (tre0_5_0 - tre0_8_0)) + (((FFTW_REAL) FFTW_K464723172) * (tre0_6_0 - tre0_7_0));
	  c_im(out[2 * ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[11 * ostride]) = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K120536680) * (tre0_1_0 + tre0_12_0)) + (((FFTW_REAL) FFTW_K885456025) * (tre0_4_0 + tre0_9_0)) + (((FFTW_REAL) FFTW_K568064746) * (tre0_5_0 + tre0_8_0)) - (((FFTW_REAL) FFTW_K748510748) * (tre0_6_0 + tre0_7_0)) - (((FFTW_REAL) FFTW_K354604887) * (tre0_3_0 + tre0_10_0)) - (((FFTW_REAL) FFTW_K970941817) * (tre0_2_0 + tre0_11_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K992708874) * (tim0_1_0 - tim0_12_0)) + (((FFTW_REAL) FFTW_K239315664) * (tim0_2_0 - tim0_11_0)) + (((FFTW_REAL) FFTW_K935016242) * (tim0_10_0 - tim0_3_0)) + (((FFTW_REAL) FFTW_K464723172) * (tim0_9_0 - tim0_4_0)) + (((FFTW_REAL) FFTW_K822983865) * (tim0_5_0 - tim0_8_0)) + (((FFTW_REAL) FFTW_K663122658) * (tim0_6_0 - tim0_7_0));
	  c_re(out[3 * ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[10 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K120536680) * (tim0_1_0 + tim0_12_0)) + (((FFTW_REAL) FFTW_K885456025) * (tim0_4_0 + tim0_9_0)) + (((FFTW_REAL) FFTW_K568064746) * (tim0_5_0 + tim0_8_0)) - (((FFTW_REAL) FFTW_K748510748) * (tim0_6_0 + tim0_7_0)) - (((FFTW_REAL) FFTW_K354604887) * (tim0_3_0 + tim0_10_0)) - (((FFTW_REAL) FFTW_K970941817) * (tim0_2_0 + tim0_11_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K992708874) * (tre0_12_0 - tre0_1_0)) + (((FFTW_REAL) FFTW_K239315664) * (tre0_11_0 - tre0_2_0)) + (((FFTW_REAL) FFTW_K935016242) * (tre0_3_0 - tre0_10_0)) + (((FFTW_REAL) FFTW_K464723172) * (tre0_4_0 - tre0_9_0)) + (((FFTW_REAL) FFTW_K822983865) * (tre0_8_0 - tre0_5_0)) + (((FFTW_REAL) FFTW_K663122658) * (tre0_7_0 - tre0_6_0));
	  c_im(out[3 * ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[10 * ostride]) = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K885456025) * (tre0_3_0 + tre0_10_0)) + (((FFTW_REAL) FFTW_K120536680) * (tre0_4_0 + tre0_9_0)) + (((FFTW_REAL) FFTW_K568064746) * (tre0_6_0 + tre0_7_0)) - (((FFTW_REAL) FFTW_K970941817) * (tre0_5_0 + tre0_8_0)) - (((FFTW_REAL) FFTW_K748510748) * (tre0_2_0 + tre0_11_0)) - (((FFTW_REAL) FFTW_K354604887) * (tre0_1_0 + tre0_12_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K935016242) * (tim0_1_0 - tim0_12_0)) + (((FFTW_REAL) FFTW_K663122658) * (tim0_11_0 - tim0_2_0)) + (((FFTW_REAL) FFTW_K464723172) * (tim0_10_0 - tim0_3_0)) + (((FFTW_REAL) FFTW_K992708874) * (tim0_4_0 - tim0_9_0)) + (((FFTW_REAL) FFTW_K239315664) * (tim0_8_0 - tim0_5_0)) + (((FFTW_REAL) FFTW_K822983865) * (tim0_7_0 - tim0_6_0));
	  c_re(out[4 * ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[9 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K885456025) * (tim0_3_0 + tim0_10_0)) + (((FFTW_REAL) FFTW_K120536680) * (tim0_4_0 + tim0_9_0)) + (((FFTW_REAL) FFTW_K568064746) * (tim0_6_0 + tim0_7_0)) - (((FFTW_REAL) FFTW_K970941817) * (tim0_5_0 + tim0_8_0)) - (((FFTW_REAL) FFTW_K748510748) * (tim0_2_0 + tim0_11_0)) - (((FFTW_REAL) FFTW_K354604887) * (tim0_1_0 + tim0_12_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K935016242) * (tre0_12_0 - tre0_1_0)) + (((FFTW_REAL) FFTW_K663122658) * (tre0_2_0 - tre0_11_0)) + (((FFTW_REAL) FFTW_K464723172) * (tre0_3_0 - tre0_10_0)) + (((FFTW_REAL) FFTW_K992708874) * (tre0_9_0 - tre0_4_0)) + (((FFTW_REAL) FFTW_K239315664) * (tre0_5_0 - tre0_8_0)) + (((FFTW_REAL) FFTW_K822983865) * (tre0_6_0 - tre0_7_0));
	  c_im(out[4 * ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[9 * ostride]) = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K120536680) * (tre0_2_0 + tre0_11_0)) + (((FFTW_REAL) FFTW_K568064746) * (tre0_3_0 + tre0_10_0)) + (((FFTW_REAL) FFTW_K885456025) * (tre0_5_0 + tre0_8_0)) - (((FFTW_REAL) FFTW_K354604887) * (tre0_6_0 + tre0_7_0)) - (((FFTW_REAL) FFTW_K970941817) * (tre0_4_0 + tre0_9_0)) - (((FFTW_REAL) FFTW_K748510748) * (tre0_1_0 + tre0_12_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K663122658) * (tim0_1_0 - tim0_12_0)) + (((FFTW_REAL) FFTW_K992708874) * (tim0_11_0 - tim0_2_0)) + (((FFTW_REAL) FFTW_K822983865) * (tim0_3_0 - tim0_10_0)) + (((FFTW_REAL) FFTW_K239315664) * (tim0_9_0 - tim0_4_0)) + (((FFTW_REAL) FFTW_K464723172) * (tim0_8_0 - tim0_5_0)) + (((FFTW_REAL) FFTW_K935016242) * (tim0_6_0 - tim0_7_0));
	  c_re(out[5 * ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[8 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K120536680) * (tim0_2_0 + tim0_11_0)) + (((FFTW_REAL) FFTW_K568064746) * (tim0_3_0 + tim0_10_0)) + (((FFTW_REAL) FFTW_K885456025) * (tim0_5_0 + tim0_8_0)) - (((FFTW_REAL) FFTW_K354604887) * (tim0_6_0 + tim0_7_0)) - (((FFTW_REAL) FFTW_K970941817) * (tim0_4_0 + tim0_9_0)) - (((FFTW_REAL) FFTW_K748510748) * (tim0_1_0 + tim0_12_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K663122658) * (tre0_12_0 - tre0_1_0)) + (((FFTW_REAL) FFTW_K992708874) * (tre0_2_0 - tre0_11_0)) + (((FFTW_REAL) FFTW_K822983865) * (tre0_10_0 - tre0_3_0)) + (((FFTW_REAL) FFTW_K239315664) * (tre0_4_0 - tre0_9_0)) + (((FFTW_REAL) FFTW_K464723172) * (tre0_5_0 - tre0_8_0)) + (((FFTW_REAL) FFTW_K935016242) * (tre0_7_0 - tre0_6_0));
	  c_im(out[5 * ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[8 * ostride]) = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K885456025) * (tre0_2_0 + tre0_11_0)) + (((FFTW_REAL) FFTW_K568064746) * (tre0_4_0 + tre0_9_0)) + (((FFTW_REAL) FFTW_K120536680) * (tre0_6_0 + tre0_7_0)) - (((FFTW_REAL) FFTW_K354604887) * (tre0_5_0 + tre0_8_0)) - (((FFTW_REAL) FFTW_K748510748) * (tre0_3_0 + tre0_10_0)) - (((FFTW_REAL) FFTW_K970941817) * (tre0_1_0 + tre0_12_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K239315664) * (tim0_1_0 - tim0_12_0)) + (((FFTW_REAL) FFTW_K464723172) * (tim0_11_0 - tim0_2_0)) + (((FFTW_REAL) FFTW_K663122658) * (tim0_3_0 - tim0_10_0)) + (((FFTW_REAL) FFTW_K822983865) * (tim0_9_0 - tim0_4_0)) + (((FFTW_REAL) FFTW_K935016242) * (tim0_5_0 - tim0_8_0)) + (((FFTW_REAL) FFTW_K992708874) * (tim0_7_0 - tim0_6_0));
	  c_re(out[6 * ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[7 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K885456025) * (tim0_2_0 + tim0_11_0)) + (((FFTW_REAL) FFTW_K568064746) * (tim0_4_0 + tim0_9_0)) + (((FFTW_REAL) FFTW_K120536680) * (tim0_6_0 + tim0_7_0)) - (((FFTW_REAL) FFTW_K354604887) * (tim0_5_0 + tim0_8_0)) - (((FFTW_REAL) FFTW_K748510748) * (tim0_3_0 + tim0_10_0)) - (((FFTW_REAL) FFTW_K970941817) * (tim0_1_0 + tim0_12_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K239315664) * (tre0_12_0 - tre0_1_0)) + (((FFTW_REAL) FFTW_K464723172) * (tre0_2_0 - tre0_11_0)) + (((FFTW_REAL) FFTW_K663122658) * (tre0_10_0 - tre0_3_0)) + (((FFTW_REAL) FFTW_K822983865) * (tre0_4_0 - tre0_9_0)) + (((FFTW_REAL) FFTW_K935016242) * (tre0_8_0 - tre0_5_0)) + (((FFTW_REAL) FFTW_K992708874) * (tre0_6_0 - tre0_7_0));
	  c_im(out[6 * ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[7 * ostride]) = tim1_0_0 - tim1_1_0;
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 208 FP additions and 72 FP multiplications */

void fftw_no_twiddle_14(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_0_1;
     FFTW_REAL tim0_0_1;
     FFTW_REAL tre0_0_2;
     FFTW_REAL tim0_0_2;
     FFTW_REAL tre0_0_3;
     FFTW_REAL tim0_0_3;
     FFTW_REAL tre0_0_4;
     FFTW_REAL tim0_0_4;
     FFTW_REAL tre0_0_5;
     FFTW_REAL tim0_0_5;
     FFTW_REAL tre0_0_6;
     FFTW_REAL tim0_0_6;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_1_1;
     FFTW_REAL tim0_1_1;
     FFTW_REAL tre0_1_2;
     FFTW_REAL tim0_1_2;
     FFTW_REAL tre0_1_3;
     FFTW_REAL tim0_1_3;
     FFTW_REAL tre0_1_4;
     FFTW_REAL tim0_1_4;
     FFTW_REAL tre0_1_5;
     FFTW_REAL tim0_1_5;
     FFTW_REAL tre0_1_6;
     FFTW_REAL tim0_1_6;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[0]);
	  tim1_0_0 = c_im(in[0]);
	  tre1_1_0 = c_re(in[7 * istride]);
	  tim1_1_0 = c_im(in[7 * istride]);
	  tre0_0_0 = tre1_0_0 + tre1_1_0;
	  tim0_0_0 = tim1_0_0 + tim1_1_0;
	  tre0_1_0 = tre1_0_0 - tre1_1_0;
	  tim0_1_0 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[2 * istride]);
	  tim1_0_0 = c_im(in[2 * istride]);
	  tre1_1_0 = c_re(in[9 * istride]);
	  tim1_1_0 = c_im(in[9 * istride]);
	  tre0_0_1 = tre1_0_0 + tre1_1_0;
	  tim0_0_1 = tim1_0_0 + tim1_1_0;
	  tre0_1_1 = tre1_0_0 - tre1_1_0;
	  tim0_1_1 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[4 * istride]);
	  tim1_0_0 = c_im(in[4 * istride]);
	  tre1_1_0 = c_re(in[11 * istride]);
	  tim1_1_0 = c_im(in[11 * istride]);
	  tre0_0_2 = tre1_0_0 + tre1_1_0;
	  tim0_0_2 = tim1_0_0 + tim1_1_0;
	  tre0_1_2 = tre1_0_0 - tre1_1_0;
	  tim0_1_2 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[6 * istride]);
	  tim1_0_0 = c_im(in[6 * istride]);
	  tre1_1_0 = c_re(in[13 * istride]);
	  tim1_1_0 = c_im(in[13 * istride]);
	  tre0_0_3 = tre1_0_0 + tre1_1_0;
	  tim0_0_3 = tim1_0_0 + tim1_1_0;
	  tre0_1_3 = tre1_0_0 - tre1_1_0;
	  tim0_1_3 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[8 * istride]);
	  tim1_0_0 = c_im(in[8 * istride]);
	  tre1_1_0 = c_re(in[istride]);
	  tim1_1_0 = c_im(in[istride]);
	  tre0_0_4 = tre1_0_0 + tre1_1_0;
	  tim0_0_4 = tim1_0_0 + tim1_1_0;
	  tre0_1_4 = tre1_0_0 - tre1_1_0;
	  tim0_1_4 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[10 * istride]);
	  tim1_0_0 = c_im(in[10 * istride]);
	  tre1_1_0 = c_re(in[3 * istride]);
	  tim1_1_0 = c_im(in[3 * istride]);
	  tre0_0_5 = tre1_0_0 + tre1_1_0;
	  tim0_0_5 = tim1_0_0 + tim1_1_0;
	  tre0_1_5 = tre1_0_0 - tre1_1_0;
	  tim0_1_5 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[12 * istride]);
	  tim1_0_0 = c_im(in[12 * istride]);
	  tre1_1_0 = c_re(in[5 * istride]);
	  tim1_1_0 = c_im(in[5 * istride]);
	  tre0_0_6 = tre1_0_0 + tre1_1_0;
	  tim0_0_6 = tim1_0_0 + tim1_1_0;
	  tre0_1_6 = tre1_0_0 - tre1_1_0;
	  tim0_1_6 = tim1_0_0 - tim1_1_0;
     }
     c_re(out[0]) = tre0_0_0 + tre0_0_1 + tre0_0_2 + tre0_0_3 + tre0_0_4 + tre0_0_5 + tre0_0_6;
     c_im(out[0]) = tim0_0_0 + tim0_0_1 + tim0_0_2 + tim0_0_3 + tim0_0_4 + tim0_0_5 + tim0_0_6;
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tre0_0_1 + tre0_0_6)) - (((FFTW_REAL) FFTW_K900968867) * (tre0_0_3 + tre0_0_4)) - (((FFTW_REAL) FFTW_K222520933) * (tre0_0_2 + tre0_0_5));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K781831482) * (tim0_0_1 - tim0_0_6)) + (((FFTW_REAL) FFTW_K974927912) * (tim0_0_2 - tim0_0_5)) + (((FFTW_REAL) FFTW_K433883739) * (tim0_0_3 - tim0_0_4));
	  c_re(out[8 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[6 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tim0_0_1 + tim0_0_6)) - (((FFTW_REAL) FFTW_K900968867) * (tim0_0_3 + tim0_0_4)) - (((FFTW_REAL) FFTW_K222520933) * (tim0_0_2 + tim0_0_5));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K781831482) * (tre0_0_6 - tre0_0_1)) + (((FFTW_REAL) FFTW_K974927912) * (tre0_0_5 - tre0_0_2)) + (((FFTW_REAL) FFTW_K433883739) * (tre0_0_4 - tre0_0_3));
	  c_im(out[8 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[6 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tre0_0_3 + tre0_0_4)) - (((FFTW_REAL) FFTW_K900968867) * (tre0_0_2 + tre0_0_5)) - (((FFTW_REAL) FFTW_K222520933) * (tre0_0_1 + tre0_0_6));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K974927912) * (tim0_0_1 - tim0_0_6)) + (((FFTW_REAL) FFTW_K433883739) * (tim0_0_5 - tim0_0_2)) + (((FFTW_REAL) FFTW_K781831482) * (tim0_0_4 - tim0_0_3));
	  c_re(out[2 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[12 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tim0_0_3 + tim0_0_4)) - (((FFTW_REAL) FFTW_K900968867) * (tim0_0_2 + tim0_0_5)) - (((FFTW_REAL) FFTW_K222520933) * (tim0_0_1 + tim0_0_6));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K974927912) * (tre0_0_6 - tre0_0_1)) + (((FFTW_REAL) FFTW_K433883739) * (tre0_0_2 - tre0_0_5)) + (((FFTW_REAL) FFTW_K781831482) * (tre0_0_3 - tre0_0_4));
	  c_im(out[2 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[12 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tre0_0_2 + tre0_0_5)) - (((FFTW_REAL) FFTW_K222520933) * (tre0_0_3 + tre0_0_4)) - (((FFTW_REAL) FFTW_K900968867) * (tre0_0_1 + tre0_0_6));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K433883739) * (tim0_0_1 - tim0_0_6)) + (((FFTW_REAL) FFTW_K781831482) * (tim0_0_5 - tim0_0_2)) + (((FFTW_REAL) FFTW_K974927912) * (tim0_0_3 - tim0_0_4));
	  c_re(out[10 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[4 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tim0_0_2 + tim0_0_5)) - (((FFTW_REAL) FFTW_K222520933) * (tim0_0_3 + tim0_0_4)) - (((FFTW_REAL) FFTW_K900968867) * (tim0_0_1 + tim0_0_6));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K433883739) * (tre0_0_6 - tre0_0_1)) + (((FFTW_REAL) FFTW_K781831482) * (tre0_0_2 - tre0_0_5)) + (((FFTW_REAL) FFTW_K974927912) * (tre0_0_4 - tre0_0_3));
	  c_im(out[10 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[4 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     c_re(out[7 * ostride]) = tre0_1_0 + tre0_1_1 + tre0_1_2 + tre0_1_3 + tre0_1_4 + tre0_1_5 + tre0_1_6;
     c_im(out[7 * ostride]) = tim0_1_0 + tim0_1_1 + tim0_1_2 + tim0_1_3 + tim0_1_4 + tim0_1_5 + tim0_1_6;
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_1_0 + (((FFTW_REAL) FFTW_K623489801) * (tre0_1_1 + tre0_1_6)) - (((FFTW_REAL) FFTW_K900968867) * (tre0_1_3 + tre0_1_4)) - (((FFTW_REAL) FFTW_K222520933) * (tre0_1_2 + tre0_1_5));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K781831482) * (tim0_1_1 - tim0_1_6)) + (((FFTW_REAL) FFTW_K974927912) * (tim0_1_2 - tim0_1_5)) + (((FFTW_REAL) FFTW_K433883739) * (tim0_1_3 - tim0_1_4));
	  c_re(out[ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[13 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_1_0 + (((FFTW_REAL) FFTW_K623489801) * (tim0_1_1 + tim0_1_6)) - (((FFTW_REAL) FFTW_K900968867) * (tim0_1_3 + tim0_1_4)) - (((FFTW_REAL) FFTW_K222520933) * (tim0_1_2 + tim0_1_5));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K781831482) * (tre0_1_6 - tre0_1_1)) + (((FFTW_REAL) FFTW_K974927912) * (tre0_1_5 - tre0_1_2)) + (((FFTW_REAL) FFTW_K433883739) * (tre0_1_4 - tre0_1_3));
	  c_im(out[ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[13 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_1_0 + (((FFTW_REAL) FFTW_K623489801) * (tre0_1_3 + tre0_1_4)) - (((FFTW_REAL) FFTW_K900968867) * (tre0_1_2 + tre0_1_5)) - (((FFTW_REAL) FFTW_K222520933) * (tre0_1_1 + tre0_1_6));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K974927912) * (tim0_1_1 - tim0_1_6)) + (((FFTW_REAL) FFTW_K433883739) * (tim0_1_5 - tim0_1_2)) + (((FFTW_REAL) FFTW_K781831482) * (tim0_1_4 - tim0_1_3));
	  c_re(out[9 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[5 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_1_0 + (((FFTW_REAL) FFTW_K623489801) * (tim0_1_3 + tim0_1_4)) - (((FFTW_REAL) FFTW_K900968867) * (tim0_1_2 + tim0_1_5)) - (((FFTW_REAL) FFTW_K222520933) * (tim0_1_1 + tim0_1_6));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K974927912) * (tre0_1_6 - tre0_1_1)) + (((FFTW_REAL) FFTW_K433883739) * (tre0_1_2 - tre0_1_5)) + (((FFTW_REAL) FFTW_K781831482) * (tre0_1_3 - tre0_1_4));
	  c_im(out[9 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[5 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_1_0 + (((FFTW_REAL) FFTW_K623489801) * (tre0_1_2 + tre0_1_5)) - (((FFTW_REAL) FFTW_K222520933) * (tre0_1_3 + tre0_1_4)) - (((FFTW_REAL) FFTW_K900968867) * (tre0_1_1 + tre0_1_6));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K433883739) * (tim0_1_1 - tim0_1_6)) + (((FFTW_REAL) FFTW_K781831482) * (tim0_1_5 - tim0_1_2)) + (((FFTW_REAL) FFTW_K974927912) * (tim0_1_3 - tim0_1_4));
	  c_re(out[3 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[11 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_1_0 + (((FFTW_REAL) FFTW_K623489801) * (tim0_1_2 + tim0_1_5)) - (((FFTW_REAL) FFTW_K222520933) * (tim0_1_3 + tim0_1_4)) - (((FFTW_REAL) FFTW_K900968867) * (tim0_1_1 + tim0_1_6));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K433883739) * (tre0_1_6 - tre0_1_1)) + (((FFTW_REAL) FFTW_K781831482) * (tre0_1_2 - tre0_1_5)) + (((FFTW_REAL) FFTW_K974927912) * (tre0_1_4 - tre0_1_3));
	  c_im(out[3 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[11 * ostride]) = tim2_0_0 - tim2_1_0;
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 202 FP additions and 68 FP multiplications */

void fftw_no_twiddle_15(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_0_1;
     FFTW_REAL tim0_0_1;
     FFTW_REAL tre0_0_2;
     FFTW_REAL tim0_0_2;
     FFTW_REAL tre0_0_3;
     FFTW_REAL tim0_0_3;
     FFTW_REAL tre0_0_4;
     FFTW_REAL tim0_0_4;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_1_1;
     FFTW_REAL tim0_1_1;
     FFTW_REAL tre0_1_2;
     FFTW_REAL tim0_1_2;
     FFTW_REAL tre0_1_3;
     FFTW_REAL tim0_1_3;
     FFTW_REAL tre0_1_4;
     FFTW_REAL tim0_1_4;
     FFTW_REAL tre0_2_0;
     FFTW_REAL tim0_2_0;
     FFTW_REAL tre0_2_1;
     FFTW_REAL tim0_2_1;
     FFTW_REAL tre0_2_2;
     FFTW_REAL tim0_2_2;
     FFTW_REAL tre0_2_3;
     FFTW_REAL tim0_2_3;
     FFTW_REAL tre0_2_4;
     FFTW_REAL tim0_2_4;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_2_0;
	  FFTW_REAL tim1_2_0;
	  tre1_0_0 = c_re(in[0]);
	  tim1_0_0 = c_im(in[0]);
	  tre1_1_0 = c_re(in[5 * istride]);
	  tim1_1_0 = c_im(in[5 * istride]);
	  tre1_2_0 = c_re(in[10 * istride]);
	  tim1_2_0 = c_im(in[10 * istride]);
	  tre0_0_0 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	  tim0_0_0 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_1_0 - tim1_2_0);
	       tre0_1_0 = tre2_0_0 + tre2_1_0;
	       tre0_2_0 = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_2_0 - tre1_1_0);
	       tim0_1_0 = tim2_0_0 + tim2_1_0;
	       tim0_2_0 = tim2_0_0 - tim2_1_0;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_2_0;
	  FFTW_REAL tim1_2_0;
	  tre1_0_0 = c_re(in[3 * istride]);
	  tim1_0_0 = c_im(in[3 * istride]);
	  tre1_1_0 = c_re(in[8 * istride]);
	  tim1_1_0 = c_im(in[8 * istride]);
	  tre1_2_0 = c_re(in[13 * istride]);
	  tim1_2_0 = c_im(in[13 * istride]);
	  tre0_0_1 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	  tim0_0_1 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_1_0 - tim1_2_0);
	       tre0_1_1 = tre2_0_0 + tre2_1_0;
	       tre0_2_1 = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_2_0 - tre1_1_0);
	       tim0_1_1 = tim2_0_0 + tim2_1_0;
	       tim0_2_1 = tim2_0_0 - tim2_1_0;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_2_0;
	  FFTW_REAL tim1_2_0;
	  tre1_0_0 = c_re(in[6 * istride]);
	  tim1_0_0 = c_im(in[6 * istride]);
	  tre1_1_0 = c_re(in[11 * istride]);
	  tim1_1_0 = c_im(in[11 * istride]);
	  tre1_2_0 = c_re(in[istride]);
	  tim1_2_0 = c_im(in[istride]);
	  tre0_0_2 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	  tim0_0_2 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_1_0 - tim1_2_0);
	       tre0_1_2 = tre2_0_0 + tre2_1_0;
	       tre0_2_2 = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_2_0 - tre1_1_0);
	       tim0_1_2 = tim2_0_0 + tim2_1_0;
	       tim0_2_2 = tim2_0_0 - tim2_1_0;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_2_0;
	  FFTW_REAL tim1_2_0;
	  tre1_0_0 = c_re(in[9 * istride]);
	  tim1_0_0 = c_im(in[9 * istride]);
	  tre1_1_0 = c_re(in[14 * istride]);
	  tim1_1_0 = c_im(in[14 * istride]);
	  tre1_2_0 = c_re(in[4 * istride]);
	  tim1_2_0 = c_im(in[4 * istride]);
	  tre0_0_3 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	  tim0_0_3 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_1_0 - tim1_2_0);
	       tre0_1_3 = tre2_0_0 + tre2_1_0;
	       tre0_2_3 = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_2_0 - tre1_1_0);
	       tim0_1_3 = tim2_0_0 + tim2_1_0;
	       tim0_2_3 = tim2_0_0 - tim2_1_0;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_2_0;
	  FFTW_REAL tim1_2_0;
	  tre1_0_0 = c_re(in[12 * istride]);
	  tim1_0_0 = c_im(in[12 * istride]);
	  tre1_1_0 = c_re(in[2 * istride]);
	  tim1_1_0 = c_im(in[2 * istride]);
	  tre1_2_0 = c_re(in[7 * istride]);
	  tim1_2_0 = c_im(in[7 * istride]);
	  tre0_0_4 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	  tim0_0_4 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_1_0 - tim1_2_0);
	       tre0_1_4 = tre2_0_0 + tre2_1_0;
	       tre0_2_4 = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_2_0 - tre1_1_0);
	       tim0_1_4 = tim2_0_0 + tim2_1_0;
	       tim0_2_4 = tim2_0_0 - tim2_1_0;
	  }
     }
     c_re(out[0]) = tre0_0_0 + tre0_0_1 + tre0_0_2 + tre0_0_3 + tre0_0_4;
     c_im(out[0]) = tim0_0_0 + tim0_0_1 + tim0_0_2 + tim0_0_3 + tim0_0_4;
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_0_1 + tre0_0_4)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_0_2 + tre0_0_3));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tim0_0_1 - tim0_0_4)) + (((FFTW_REAL) FFTW_K587785252) * (tim0_0_2 - tim0_0_3));
	  c_re(out[6 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[9 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_0_1 + tim0_0_4)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_0_2 + tim0_0_3));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tre0_0_4 - tre0_0_1)) + (((FFTW_REAL) FFTW_K587785252) * (tre0_0_3 - tre0_0_2));
	  c_im(out[6 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[9 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_0_2 + tre0_0_3)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_0_1 + tre0_0_4));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tim0_0_1 - tim0_0_4)) + (((FFTW_REAL) FFTW_K951056516) * (tim0_0_3 - tim0_0_2));
	  c_re(out[12 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[3 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_0_2 + tim0_0_3)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_0_1 + tim0_0_4));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tre0_0_4 - tre0_0_1)) + (((FFTW_REAL) FFTW_K951056516) * (tre0_0_2 - tre0_0_3));
	  c_im(out[12 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[3 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     c_re(out[10 * ostride]) = tre0_1_0 + tre0_1_1 + tre0_1_2 + tre0_1_3 + tre0_1_4;
     c_im(out[10 * ostride]) = tim0_1_0 + tim0_1_1 + tim0_1_2 + tim0_1_3 + tim0_1_4;
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_1_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_1_1 + tre0_1_4)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_1_2 + tre0_1_3));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tim0_1_1 - tim0_1_4)) + (((FFTW_REAL) FFTW_K587785252) * (tim0_1_2 - tim0_1_3));
	  c_re(out[ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[4 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_1_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_1_1 + tim0_1_4)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_1_2 + tim0_1_3));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tre0_1_4 - tre0_1_1)) + (((FFTW_REAL) FFTW_K587785252) * (tre0_1_3 - tre0_1_2));
	  c_im(out[ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[4 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_1_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_1_2 + tre0_1_3)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_1_1 + tre0_1_4));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tim0_1_1 - tim0_1_4)) + (((FFTW_REAL) FFTW_K951056516) * (tim0_1_3 - tim0_1_2));
	  c_re(out[7 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[13 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_1_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_1_2 + tim0_1_3)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_1_1 + tim0_1_4));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tre0_1_4 - tre0_1_1)) + (((FFTW_REAL) FFTW_K951056516) * (tre0_1_2 - tre0_1_3));
	  c_im(out[7 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[13 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     c_re(out[5 * ostride]) = tre0_2_0 + tre0_2_1 + tre0_2_2 + tre0_2_3 + tre0_2_4;
     c_im(out[5 * ostride]) = tim0_2_0 + tim0_2_1 + tim0_2_2 + tim0_2_3 + tim0_2_4;
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_2_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_2_1 + tre0_2_4)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_2_2 + tre0_2_3));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tim0_2_1 - tim0_2_4)) + (((FFTW_REAL) FFTW_K587785252) * (tim0_2_2 - tim0_2_3));
	  c_re(out[11 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[14 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_2_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_2_1 + tim0_2_4)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_2_2 + tim0_2_3));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tre0_2_4 - tre0_2_1)) + (((FFTW_REAL) FFTW_K587785252) * (tre0_2_3 - tre0_2_2));
	  c_im(out[11 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[14 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_2_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_2_2 + tre0_2_3)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_2_1 + tre0_2_4));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tim0_2_1 - tim0_2_4)) + (((FFTW_REAL) FFTW_K951056516) * (tim0_2_3 - tim0_2_2));
	  c_re(out[2 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[8 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_2_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_2_2 + tim0_2_3)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_2_1 + tim0_2_4));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tre0_2_4 - tre0_2_1)) + (((FFTW_REAL) FFTW_K951056516) * (tre0_2_2 - tre0_2_3));
	  c_im(out[2 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[8 * ostride]) = tim2_0_0 - tim2_1_0;
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 144 FP additions and 24 FP multiplications */

void fftw_no_twiddle_16(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_0_1;
     FFTW_REAL tim0_0_1;
     FFTW_REAL tre0_0_2;
     FFTW_REAL tim0_0_2;
     FFTW_REAL tre0_0_3;
     FFTW_REAL tim0_0_3;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_1_1;
     FFTW_REAL tim0_1_1;
     FFTW_REAL tre0_1_2;
     FFTW_REAL tim0_1_2;
     FFTW_REAL tre0_1_3;
     FFTW_REAL tim0_1_3;
     FFTW_REAL tre0_2_0;
     FFTW_REAL tim0_2_0;
     FFTW_REAL tre0_2_1;
     FFTW_REAL tim0_2_1;
     FFTW_REAL tre0_2_2;
     FFTW_REAL tim0_2_2;
     FFTW_REAL tre0_2_3;
     FFTW_REAL tim0_2_3;
     FFTW_REAL tre0_3_0;
     FFTW_REAL tim0_3_0;
     FFTW_REAL tre0_3_1;
     FFTW_REAL tim0_3_1;
     FFTW_REAL tre0_3_2;
     FFTW_REAL tim0_3_2;
     FFTW_REAL tre0_3_3;
     FFTW_REAL tim0_3_3;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[0]);
	       tim2_0_0 = c_im(in[0]);
	       tre2_1_0 = c_re(in[8 * istride]);
	       tim2_1_0 = c_im(in[8 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[4 * istride]);
	       tim2_0_0 = c_im(in[4 * istride]);
	       tre2_1_0 = c_re(in[12 * istride]);
	       tim2_1_0 = c_im(in[12 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  tre0_0_0 = tre1_0_0 + tre1_0_1;
	  tim0_0_0 = tim1_0_0 + tim1_0_1;
	  tre0_2_0 = tre1_0_0 - tre1_0_1;
	  tim0_2_0 = tim1_0_0 - tim1_0_1;
	  tre0_1_0 = tre1_1_0 + tim1_1_1;
	  tim0_1_0 = tim1_1_0 - tre1_1_1;
	  tre0_3_0 = tre1_1_0 - tim1_1_1;
	  tim0_3_0 = tim1_1_0 + tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[istride]);
	       tim2_0_0 = c_im(in[istride]);
	       tre2_1_0 = c_re(in[9 * istride]);
	       tim2_1_0 = c_im(in[9 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[5 * istride]);
	       tim2_0_0 = c_im(in[5 * istride]);
	       tre2_1_0 = c_re(in[13 * istride]);
	       tim2_1_0 = c_im(in[13 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  tre0_0_1 = tre1_0_0 + tre1_0_1;
	  tim0_0_1 = tim1_0_0 + tim1_0_1;
	  tre0_2_1 = tre1_0_0 - tre1_0_1;
	  tim0_2_1 = tim1_0_0 - tim1_0_1;
	  tre0_1_1 = tre1_1_0 + tim1_1_1;
	  tim0_1_1 = tim1_1_0 - tre1_1_1;
	  tre0_3_1 = tre1_1_0 - tim1_1_1;
	  tim0_3_1 = tim1_1_0 + tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[2 * istride]);
	       tim2_0_0 = c_im(in[2 * istride]);
	       tre2_1_0 = c_re(in[10 * istride]);
	       tim2_1_0 = c_im(in[10 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[6 * istride]);
	       tim2_0_0 = c_im(in[6 * istride]);
	       tre2_1_0 = c_re(in[14 * istride]);
	       tim2_1_0 = c_im(in[14 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  tre0_0_2 = tre1_0_0 + tre1_0_1;
	  tim0_0_2 = tim1_0_0 + tim1_0_1;
	  tre0_2_2 = tre1_0_0 - tre1_0_1;
	  tim0_2_2 = tim1_0_0 - tim1_0_1;
	  tre0_1_2 = tre1_1_0 + tim1_1_1;
	  tim0_1_2 = tim1_1_0 - tre1_1_1;
	  tre0_3_2 = tre1_1_0 - tim1_1_1;
	  tim0_3_2 = tim1_1_0 + tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[3 * istride]);
	       tim2_0_0 = c_im(in[3 * istride]);
	       tre2_1_0 = c_re(in[11 * istride]);
	       tim2_1_0 = c_im(in[11 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[7 * istride]);
	       tim2_0_0 = c_im(in[7 * istride]);
	       tre2_1_0 = c_re(in[15 * istride]);
	       tim2_1_0 = c_im(in[15 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  tre0_0_3 = tre1_0_0 + tre1_0_1;
	  tim0_0_3 = tim1_0_0 + tim1_0_1;
	  tre0_2_3 = tre1_0_0 - tre1_0_1;
	  tim0_2_3 = tim1_0_0 - tim1_0_1;
	  tre0_1_3 = tre1_1_0 + tim1_1_1;
	  tim0_1_3 = tim1_1_0 - tre1_1_1;
	  tre0_3_3 = tre1_1_0 - tim1_1_1;
	  tim0_3_3 = tim1_1_0 + tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  tre1_0_0 = tre0_0_0 + tre0_0_2;
	  tim1_0_0 = tim0_0_0 + tim0_0_2;
	  tre1_1_0 = tre0_0_0 - tre0_0_2;
	  tim1_1_0 = tim0_0_0 - tim0_0_2;
	  tre1_0_1 = tre0_0_1 + tre0_0_3;
	  tim1_0_1 = tim0_0_1 + tim0_0_3;
	  tre1_1_1 = tre0_0_1 - tre0_0_3;
	  tim1_1_1 = tim0_0_1 - tim0_0_3;
	  c_re(out[0]) = tre1_0_0 + tre1_0_1;
	  c_im(out[0]) = tim1_0_0 + tim1_0_1;
	  c_re(out[8 * ostride]) = tre1_0_0 - tre1_0_1;
	  c_im(out[8 * ostride]) = tim1_0_0 - tim1_0_1;
	  c_re(out[4 * ostride]) = tre1_1_0 + tim1_1_1;
	  c_im(out[4 * ostride]) = tim1_1_0 - tre1_1_1;
	  c_re(out[12 * ostride]) = tre1_1_0 - tim1_1_1;
	  c_im(out[12 * ostride]) = tim1_1_0 + tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  {
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_1_2 + tim0_1_2);
	       tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_1_2 - tre0_1_2);
	       tre1_0_0 = tre0_1_0 + tre2_1_0;
	       tim1_0_0 = tim0_1_0 + tim2_1_0;
	       tre1_1_0 = tre0_1_0 - tre2_1_0;
	       tim1_1_0 = tim0_1_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_1_1) + (((FFTW_REAL) FFTW_K382683432) * tim0_1_1);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_1_1) - (((FFTW_REAL) FFTW_K382683432) * tre0_1_1);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_1_3) + (((FFTW_REAL) FFTW_K923879532) * tim0_1_3);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_1_3) - (((FFTW_REAL) FFTW_K923879532) * tre0_1_3);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  c_re(out[ostride]) = tre1_0_0 + tre1_0_1;
	  c_im(out[ostride]) = tim1_0_0 + tim1_0_1;
	  c_re(out[9 * ostride]) = tre1_0_0 - tre1_0_1;
	  c_im(out[9 * ostride]) = tim1_0_0 - tim1_0_1;
	  c_re(out[5 * ostride]) = tre1_1_0 + tim1_1_1;
	  c_im(out[5 * ostride]) = tim1_1_0 - tre1_1_1;
	  c_re(out[13 * ostride]) = tre1_1_0 - tim1_1_1;
	  c_im(out[13 * ostride]) = tim1_1_0 + tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  tre1_0_0 = tre0_2_0 + tim0_2_2;
	  tim1_0_0 = tim0_2_0 - tre0_2_2;
	  tre1_1_0 = tre0_2_0 - tim0_2_2;
	  tim1_1_0 = tim0_2_0 + tre0_2_2;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_2_1 + tim0_2_1);
	       tim2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_2_1 - tre0_2_1);
	       tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_2_3 - tre0_2_3);
	       tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_2_3 + tre0_2_3);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 - tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 + tim2_1_0;
	  }
	  c_re(out[2 * ostride]) = tre1_0_0 + tre1_0_1;
	  c_im(out[2 * ostride]) = tim1_0_0 + tim1_0_1;
	  c_re(out[10 * ostride]) = tre1_0_0 - tre1_0_1;
	  c_im(out[10 * ostride]) = tim1_0_0 - tim1_0_1;
	  c_re(out[6 * ostride]) = tre1_1_0 + tim1_1_1;
	  c_im(out[6 * ostride]) = tim1_1_0 - tre1_1_1;
	  c_re(out[14 * ostride]) = tre1_1_0 - tim1_1_1;
	  c_im(out[14 * ostride]) = tim1_1_0 + tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  {
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_3_2 - tre0_3_2);
	       tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_3_2 + tre0_3_2);
	       tre1_0_0 = tre0_3_0 + tre2_1_0;
	       tim1_0_0 = tim0_3_0 - tim2_1_0;
	       tre1_1_0 = tre0_3_0 - tre2_1_0;
	       tim1_1_0 = tim0_3_0 + tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_3_1) + (((FFTW_REAL) FFTW_K923879532) * tim0_3_1);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_3_1) - (((FFTW_REAL) FFTW_K923879532) * tre0_3_1);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_3_3) + (((FFTW_REAL) FFTW_K382683432) * tim0_3_3);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_3_3) - (((FFTW_REAL) FFTW_K923879532) * tim0_3_3);
	       tre1_0_1 = tre2_0_0 - tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 + tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  c_re(out[3 * ostride]) = tre1_0_0 + tre1_0_1;
	  c_im(out[3 * ostride]) = tim1_0_0 + tim1_0_1;
	  c_re(out[11 * ostride]) = tre1_0_0 - tre1_0_1;
	  c_im(out[11 * ostride]) = tim1_0_0 - tim1_0_1;
	  c_re(out[7 * ostride]) = tre1_1_0 + tim1_1_1;
	  c_im(out[7 * ostride]) = tim1_1_0 - tre1_1_1;
	  c_re(out[15 * ostride]) = tre1_1_0 - tim1_1_1;
	  c_im(out[15 * ostride]) = tim1_1_0 + tre1_1_1;
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 4 FP additions and 0 FP multiplications */

void fftw_no_twiddle_2(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     tre0_0_0 = c_re(in[0]);
     tim0_0_0 = c_im(in[0]);
     tre0_1_0 = c_re(in[istride]);
     tim0_1_0 = c_im(in[istride]);
     c_re(out[0]) = tre0_0_0 + tre0_1_0;
     c_im(out[0]) = tim0_0_0 + tim0_1_0;
     c_re(out[ostride]) = tre0_0_0 - tre0_1_0;
     c_im(out[ostride]) = tim0_0_0 - tim0_1_0;
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 14 FP additions and 4 FP multiplications */

void fftw_no_twiddle_3(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_2_0;
     FFTW_REAL tim0_2_0;
     tre0_0_0 = c_re(in[0]);
     tim0_0_0 = c_im(in[0]);
     tre0_1_0 = c_re(in[istride]);
     tim0_1_0 = c_im(in[istride]);
     tre0_2_0 = c_re(in[2 * istride]);
     tim0_2_0 = c_im(in[2 * istride]);
     c_re(out[0]) = tre0_0_0 + tre0_1_0 + tre0_2_0;
     c_im(out[0]) = tim0_0_0 + tim0_1_0 + tim0_2_0;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre0_1_0 + tre0_2_0));
	  tre1_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim0_1_0 - tim0_2_0);
	  c_re(out[ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[2 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim0_1_0 + tim0_2_0));
	  tim1_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre0_2_0 - tre0_1_0);
	  c_im(out[ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[2 * ostride]) = tim1_0_0 - tim1_1_0;
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 376 FP additions and 88 FP multiplications */

void fftw_no_twiddle_32(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_0_1;
     FFTW_REAL tim0_0_1;
     FFTW_REAL tre0_0_2;
     FFTW_REAL tim0_0_2;
     FFTW_REAL tre0_0_3;
     FFTW_REAL tim0_0_3;
     FFTW_REAL tre0_0_4;
     FFTW_REAL tim0_0_4;
     FFTW_REAL tre0_0_5;
     FFTW_REAL tim0_0_5;
     FFTW_REAL tre0_0_6;
     FFTW_REAL tim0_0_6;
     FFTW_REAL tre0_0_7;
     FFTW_REAL tim0_0_7;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_1_1;
     FFTW_REAL tim0_1_1;
     FFTW_REAL tre0_1_2;
     FFTW_REAL tim0_1_2;
     FFTW_REAL tre0_1_3;
     FFTW_REAL tim0_1_3;
     FFTW_REAL tre0_1_4;
     FFTW_REAL tim0_1_4;
     FFTW_REAL tre0_1_5;
     FFTW_REAL tim0_1_5;
     FFTW_REAL tre0_1_6;
     FFTW_REAL tim0_1_6;
     FFTW_REAL tre0_1_7;
     FFTW_REAL tim0_1_7;
     FFTW_REAL tre0_2_0;
     FFTW_REAL tim0_2_0;
     FFTW_REAL tre0_2_1;
     FFTW_REAL tim0_2_1;
     FFTW_REAL tre0_2_2;
     FFTW_REAL tim0_2_2;
     FFTW_REAL tre0_2_3;
     FFTW_REAL tim0_2_3;
     FFTW_REAL tre0_2_4;
     FFTW_REAL tim0_2_4;
     FFTW_REAL tre0_2_5;
     FFTW_REAL tim0_2_5;
     FFTW_REAL tre0_2_6;
     FFTW_REAL tim0_2_6;
     FFTW_REAL tre0_2_7;
     FFTW_REAL tim0_2_7;
     FFTW_REAL tre0_3_0;
     FFTW_REAL tim0_3_0;
     FFTW_REAL tre0_3_1;
     FFTW_REAL tim0_3_1;
     FFTW_REAL tre0_3_2;
     FFTW_REAL tim0_3_2;
     FFTW_REAL tre0_3_3;
     FFTW_REAL tim0_3_3;
     FFTW_REAL tre0_3_4;
     FFTW_REAL tim0_3_4;
     FFTW_REAL tre0_3_5;
     FFTW_REAL tim0_3_5;
     FFTW_REAL tre0_3_6;
     FFTW_REAL tim0_3_6;
     FFTW_REAL tre0_3_7;
     FFTW_REAL tim0_3_7;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[0]);
	       tim2_0_0 = c_im(in[0]);
	       tre2_1_0 = c_re(in[16 * istride]);
	       tim2_1_0 = c_im(in[16 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[8 * istride]);
	       tim2_0_0 = c_im(in[8 * istride]);
	       tre2_1_0 = c_re(in[24 * istride]);
	       tim2_1_0 = c_im(in[24 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  tre0_0_0 = tre1_0_0 + tre1_0_1;
	  tim0_0_0 = tim1_0_0 + tim1_0_1;
	  tre0_2_0 = tre1_0_0 - tre1_0_1;
	  tim0_2_0 = tim1_0_0 - tim1_0_1;
	  tre0_1_0 = tre1_1_0 + tim1_1_1;
	  tim0_1_0 = tim1_1_0 - tre1_1_1;
	  tre0_3_0 = tre1_1_0 - tim1_1_1;
	  tim0_3_0 = tim1_1_0 + tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[istride]);
	       tim2_0_0 = c_im(in[istride]);
	       tre2_1_0 = c_re(in[17 * istride]);
	       tim2_1_0 = c_im(in[17 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[9 * istride]);
	       tim2_0_0 = c_im(in[9 * istride]);
	       tre2_1_0 = c_re(in[25 * istride]);
	       tim2_1_0 = c_im(in[25 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  tre0_0_1 = tre1_0_0 + tre1_0_1;
	  tim0_0_1 = tim1_0_0 + tim1_0_1;
	  tre0_2_1 = tre1_0_0 - tre1_0_1;
	  tim0_2_1 = tim1_0_0 - tim1_0_1;
	  tre0_1_1 = tre1_1_0 + tim1_1_1;
	  tim0_1_1 = tim1_1_0 - tre1_1_1;
	  tre0_3_1 = tre1_1_0 - tim1_1_1;
	  tim0_3_1 = tim1_1_0 + tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[2 * istride]);
	       tim2_0_0 = c_im(in[2 * istride]);
	       tre2_1_0 = c_re(in[18 * istride]);
	       tim2_1_0 = c_im(in[18 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[10 * istride]);
	       tim2_0_0 = c_im(in[10 * istride]);
	       tre2_1_0 = c_re(in[26 * istride]);
	       tim2_1_0 = c_im(in[26 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  tre0_0_2 = tre1_0_0 + tre1_0_1;
	  tim0_0_2 = tim1_0_0 + tim1_0_1;
	  tre0_2_2 = tre1_0_0 - tre1_0_1;
	  tim0_2_2 = tim1_0_0 - tim1_0_1;
	  tre0_1_2 = tre1_1_0 + tim1_1_1;
	  tim0_1_2 = tim1_1_0 - tre1_1_1;
	  tre0_3_2 = tre1_1_0 - tim1_1_1;
	  tim0_3_2 = tim1_1_0 + tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[3 * istride]);
	       tim2_0_0 = c_im(in[3 * istride]);
	       tre2_1_0 = c_re(in[19 * istride]);
	       tim2_1_0 = c_im(in[19 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[11 * istride]);
	       tim2_0_0 = c_im(in[11 * istride]);
	       tre2_1_0 = c_re(in[27 * istride]);
	       tim2_1_0 = c_im(in[27 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  tre0_0_3 = tre1_0_0 + tre1_0_1;
	  tim0_0_3 = tim1_0_0 + tim1_0_1;
	  tre0_2_3 = tre1_0_0 - tre1_0_1;
	  tim0_2_3 = tim1_0_0 - tim1_0_1;
	  tre0_1_3 = tre1_1_0 + tim1_1_1;
	  tim0_1_3 = tim1_1_0 - tre1_1_1;
	  tre0_3_3 = tre1_1_0 - tim1_1_1;
	  tim0_3_3 = tim1_1_0 + tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[4 * istride]);
	       tim2_0_0 = c_im(in[4 * istride]);
	       tre2_1_0 = c_re(in[20 * istride]);
	       tim2_1_0 = c_im(in[20 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[12 * istride]);
	       tim2_0_0 = c_im(in[12 * istride]);
	       tre2_1_0 = c_re(in[28 * istride]);
	       tim2_1_0 = c_im(in[28 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  tre0_0_4 = tre1_0_0 + tre1_0_1;
	  tim0_0_4 = tim1_0_0 + tim1_0_1;
	  tre0_2_4 = tre1_0_0 - tre1_0_1;
	  tim0_2_4 = tim1_0_0 - tim1_0_1;
	  tre0_1_4 = tre1_1_0 + tim1_1_1;
	  tim0_1_4 = tim1_1_0 - tre1_1_1;
	  tre0_3_4 = tre1_1_0 - tim1_1_1;
	  tim0_3_4 = tim1_1_0 + tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[5 * istride]);
	       tim2_0_0 = c_im(in[5 * istride]);
	       tre2_1_0 = c_re(in[21 * istride]);
	       tim2_1_0 = c_im(in[21 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[13 * istride]);
	       tim2_0_0 = c_im(in[13 * istride]);
	       tre2_1_0 = c_re(in[29 * istride]);
	       tim2_1_0 = c_im(in[29 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  tre0_0_5 = tre1_0_0 + tre1_0_1;
	  tim0_0_5 = tim1_0_0 + tim1_0_1;
	  tre0_2_5 = tre1_0_0 - tre1_0_1;
	  tim0_2_5 = tim1_0_0 - tim1_0_1;
	  tre0_1_5 = tre1_1_0 + tim1_1_1;
	  tim0_1_5 = tim1_1_0 - tre1_1_1;
	  tre0_3_5 = tre1_1_0 - tim1_1_1;
	  tim0_3_5 = tim1_1_0 + tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[6 * istride]);
	       tim2_0_0 = c_im(in[6 * istride]);
	       tre2_1_0 = c_re(in[22 * istride]);
	       tim2_1_0 = c_im(in[22 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[14 * istride]);
	       tim2_0_0 = c_im(in[14 * istride]);
	       tre2_1_0 = c_re(in[30 * istride]);
	       tim2_1_0 = c_im(in[30 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  tre0_0_6 = tre1_0_0 + tre1_0_1;
	  tim0_0_6 = tim1_0_0 + tim1_0_1;
	  tre0_2_6 = tre1_0_0 - tre1_0_1;
	  tim0_2_6 = tim1_0_0 - tim1_0_1;
	  tre0_1_6 = tre1_1_0 + tim1_1_1;
	  tim0_1_6 = tim1_1_0 - tre1_1_1;
	  tre0_3_6 = tre1_1_0 - tim1_1_1;
	  tim0_3_6 = tim1_1_0 + tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[7 * istride]);
	       tim2_0_0 = c_im(in[7 * istride]);
	       tre2_1_0 = c_re(in[23 * istride]);
	       tim2_1_0 = c_im(in[23 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[15 * istride]);
	       tim2_0_0 = c_im(in[15 * istride]);
	       tre2_1_0 = c_re(in[31 * istride]);
	       tim2_1_0 = c_im(in[31 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  tre0_0_7 = tre1_0_0 + tre1_0_1;
	  tim0_0_7 = tim1_0_0 + tim1_0_1;
	  tre0_2_7 = tre1_0_0 - tre1_0_1;
	  tim0_2_7 = tim1_0_0 - tim1_0_1;
	  tre0_1_7 = tre1_1_0 + tim1_1_1;
	  tim0_1_7 = tim1_1_0 - tre1_1_1;
	  tre0_3_7 = tre1_1_0 - tim1_1_1;
	  tim0_3_7 = tim1_1_0 + tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  tre1_0_0 = tre0_0_0 + tre0_0_4;
	  tim1_0_0 = tim0_0_0 + tim0_0_4;
	  tre1_1_0 = tre0_0_0 - tre0_0_4;
	  tim1_1_0 = tim0_0_0 - tim0_0_4;
	  tre1_0_1 = tre0_0_1 + tre0_0_5;
	  tim1_0_1 = tim0_0_1 + tim0_0_5;
	  tre1_1_1 = tre0_0_1 - tre0_0_5;
	  tim1_1_1 = tim0_0_1 - tim0_0_5;
	  tre1_0_2 = tre0_0_2 + tre0_0_6;
	  tim1_0_2 = tim0_0_2 + tim0_0_6;
	  tre1_1_2 = tre0_0_2 - tre0_0_6;
	  tim1_1_2 = tim0_0_2 - tim0_0_6;
	  tre1_0_3 = tre0_0_3 + tre0_0_7;
	  tim1_0_3 = tim0_0_3 + tim0_0_7;
	  tre1_1_3 = tre0_0_3 - tre0_0_7;
	  tim1_1_3 = tim0_0_3 - tim0_0_7;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       c_re(out[0]) = tre2_0_0 + tre2_0_1;
	       c_im(out[0]) = tim2_0_0 + tim2_0_1;
	       c_re(out[16 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[16 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[8 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[8 * ostride]) = tim2_1_0 - tre2_1_1;
	       c_re(out[24 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[24 * ostride]) = tim2_1_0 + tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 + tim1_1_2;
	       tim2_0_0 = tim1_1_0 - tre1_1_2;
	       tre2_1_0 = tre1_1_0 - tim1_1_2;
	       tim2_1_0 = tim1_1_0 + tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
		    tre2_0_1 = tre3_0_0 + tre3_1_0;
		    tim2_0_1 = tim3_0_0 - tim3_1_0;
		    tre2_1_1 = tre3_0_0 - tre3_1_0;
		    tim2_1_1 = tim3_0_0 + tim3_1_0;
	       }
	       c_re(out[4 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[4 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[20 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[20 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[12 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[12 * ostride]) = tim2_1_0 - tre2_1_1;
	       c_re(out[28 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[28 * ostride]) = tim2_1_0 + tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_1_4 + tim0_1_4);
	       tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_1_4 - tre0_1_4);
	       tre1_0_0 = tre0_1_0 + tre2_1_0;
	       tim1_0_0 = tim0_1_0 + tim2_1_0;
	       tre1_1_0 = tre0_1_0 - tre2_1_0;
	       tim1_1_0 = tim0_1_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tre0_1_1) + (((FFTW_REAL) FFTW_K195090322) * tim0_1_1);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tim0_1_1) - (((FFTW_REAL) FFTW_K195090322) * tre0_1_1);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K555570233) * tre0_1_5) + (((FFTW_REAL) FFTW_K831469612) * tim0_1_5);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K555570233) * tim0_1_5) - (((FFTW_REAL) FFTW_K831469612) * tre0_1_5);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_1_2) + (((FFTW_REAL) FFTW_K382683432) * tim0_1_2);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_1_2) - (((FFTW_REAL) FFTW_K382683432) * tre0_1_2);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_1_6) + (((FFTW_REAL) FFTW_K923879532) * tim0_1_6);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_1_6) - (((FFTW_REAL) FFTW_K923879532) * tre0_1_6);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_1_3) + (((FFTW_REAL) FFTW_K555570233) * tim0_1_3);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_1_3) - (((FFTW_REAL) FFTW_K555570233) * tre0_1_3);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tre0_1_7) + (((FFTW_REAL) FFTW_K980785280) * tim0_1_7);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tim0_1_7) - (((FFTW_REAL) FFTW_K980785280) * tre0_1_7);
	       tre1_0_3 = tre2_0_0 + tre2_1_0;
	       tim1_0_3 = tim2_0_0 + tim2_1_0;
	       tre1_1_3 = tre2_0_0 - tre2_1_0;
	       tim1_1_3 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       c_re(out[ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[17 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[17 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[9 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[9 * ostride]) = tim2_1_0 - tre2_1_1;
	       c_re(out[25 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[25 * ostride]) = tim2_1_0 + tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 + tim1_1_2;
	       tim2_0_0 = tim1_1_0 - tre1_1_2;
	       tre2_1_0 = tre1_1_0 - tim1_1_2;
	       tim2_1_0 = tim1_1_0 + tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
		    tre2_0_1 = tre3_0_0 + tre3_1_0;
		    tim2_0_1 = tim3_0_0 - tim3_1_0;
		    tre2_1_1 = tre3_0_0 - tre3_1_0;
		    tim2_1_1 = tim3_0_0 + tim3_1_0;
	       }
	       c_re(out[5 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[5 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[21 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[21 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[13 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[13 * ostride]) = tim2_1_0 - tre2_1_1;
	       c_re(out[29 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[29 * ostride]) = tim2_1_0 + tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  tre1_0_0 = tre0_2_0 + tim0_2_4;
	  tim1_0_0 = tim0_2_0 - tre0_2_4;
	  tre1_1_0 = tre0_2_0 - tim0_2_4;
	  tim1_1_0 = tim0_2_0 + tre0_2_4;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_2_1) + (((FFTW_REAL) FFTW_K382683432) * tim0_2_1);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_2_1) - (((FFTW_REAL) FFTW_K382683432) * tre0_2_1);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_2_5) - (((FFTW_REAL) FFTW_K382683432) * tre0_2_5);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_2_5) + (((FFTW_REAL) FFTW_K923879532) * tre0_2_5);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 - tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 + tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_2_2 + tim0_2_2);
	       tim2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_2_2 - tre0_2_2);
	       tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_2_6 - tre0_2_6);
	       tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_2_6 + tre0_2_6);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 - tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 + tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_2_3) + (((FFTW_REAL) FFTW_K923879532) * tim0_2_3);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_2_3) - (((FFTW_REAL) FFTW_K923879532) * tre0_2_3);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_2_7) - (((FFTW_REAL) FFTW_K923879532) * tre0_2_7);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_2_7) + (((FFTW_REAL) FFTW_K382683432) * tre0_2_7);
	       tre1_0_3 = tre2_0_0 + tre2_1_0;
	       tim1_0_3 = tim2_0_0 - tim2_1_0;
	       tre1_1_3 = tre2_0_0 - tre2_1_0;
	       tim1_1_3 = tim2_0_0 + tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       c_re(out[2 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[2 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[18 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[18 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[10 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[10 * ostride]) = tim2_1_0 - tre2_1_1;
	       c_re(out[26 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[26 * ostride]) = tim2_1_0 + tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 + tim1_1_2;
	       tim2_0_0 = tim1_1_0 - tre1_1_2;
	       tre2_1_0 = tre1_1_0 - tim1_1_2;
	       tim2_1_0 = tim1_1_0 + tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
		    tre2_0_1 = tre3_0_0 + tre3_1_0;
		    tim2_0_1 = tim3_0_0 - tim3_1_0;
		    tre2_1_1 = tre3_0_0 - tre3_1_0;
		    tim2_1_1 = tim3_0_0 + tim3_1_0;
	       }
	       c_re(out[6 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[6 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[22 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[22 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[14 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[14 * ostride]) = tim2_1_0 - tre2_1_1;
	       c_re(out[30 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[30 * ostride]) = tim2_1_0 + tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_3_4 - tre0_3_4);
	       tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_3_4 + tre0_3_4);
	       tre1_0_0 = tre0_3_0 + tre2_1_0;
	       tim1_0_0 = tim0_3_0 - tim2_1_0;
	       tre1_1_0 = tre0_3_0 - tre2_1_0;
	       tim1_1_0 = tim0_3_0 + tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_3_1) + (((FFTW_REAL) FFTW_K555570233) * tim0_3_1);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_3_1) - (((FFTW_REAL) FFTW_K555570233) * tre0_3_1);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tim0_3_5) - (((FFTW_REAL) FFTW_K980785280) * tre0_3_5);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K980785280) * tim0_3_5) + (((FFTW_REAL) FFTW_K195090322) * tre0_3_5);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 - tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 + tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_3_2) + (((FFTW_REAL) FFTW_K923879532) * tim0_3_2);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_3_2) - (((FFTW_REAL) FFTW_K923879532) * tre0_3_2);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_3_6) + (((FFTW_REAL) FFTW_K382683432) * tim0_3_6);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_3_6) - (((FFTW_REAL) FFTW_K923879532) * tim0_3_6);
	       tre1_0_2 = tre2_0_0 - tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 + tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tim0_3_3) - (((FFTW_REAL) FFTW_K195090322) * tre0_3_3);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K195090322) * tim0_3_3) + (((FFTW_REAL) FFTW_K980785280) * tre0_3_3);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K555570233) * tre0_3_7) + (((FFTW_REAL) FFTW_K831469612) * tim0_3_7);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_3_7) - (((FFTW_REAL) FFTW_K555570233) * tim0_3_7);
	       tre1_0_3 = tre2_0_0 - tre2_1_0;
	       tim1_0_3 = tim2_1_0 - tim2_0_0;
	       tre1_1_3 = tre2_0_0 + tre2_1_0;
	       tim1_1_3 = (-(tim2_0_0 + tim2_1_0));
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       c_re(out[3 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[3 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[19 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[19 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[11 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[11 * ostride]) = tim2_1_0 - tre2_1_1;
	       c_re(out[27 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[27 * ostride]) = tim2_1_0 + tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 + tim1_1_2;
	       tim2_0_0 = tim1_1_0 - tre1_1_2;
	       tre2_1_0 = tre1_1_0 - tim1_1_2;
	       tim2_1_0 = tim1_1_0 + tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
		    tre2_0_1 = tre3_0_0 + tre3_1_0;
		    tim2_0_1 = tim3_0_0 - tim3_1_0;
		    tre2_1_1 = tre3_0_0 - tre3_1_0;
		    tim2_1_1 = tim3_0_0 + tim3_1_0;
	       }
	       c_re(out[7 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[7 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[23 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[23 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[15 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[15 * ostride]) = tim2_1_0 - tre2_1_1;
	       c_re(out[31 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[31 * ostride]) = tim2_1_0 + tre2_1_1;
	  }
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 16 FP additions and 0 FP multiplications */

void fftw_no_twiddle_4(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_0_1;
     FFTW_REAL tim0_0_1;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_1_1;
     FFTW_REAL tim0_1_1;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[0]);
	  tim1_0_0 = c_im(in[0]);
	  tre1_1_0 = c_re(in[2 * istride]);
	  tim1_1_0 = c_im(in[2 * istride]);
	  tre0_0_0 = tre1_0_0 + tre1_1_0;
	  tim0_0_0 = tim1_0_0 + tim1_1_0;
	  tre0_1_0 = tre1_0_0 - tre1_1_0;
	  tim0_1_0 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[istride]);
	  tim1_0_0 = c_im(in[istride]);
	  tre1_1_0 = c_re(in[3 * istride]);
	  tim1_1_0 = c_im(in[3 * istride]);
	  tre0_0_1 = tre1_0_0 + tre1_1_0;
	  tim0_0_1 = tim1_0_0 + tim1_1_0;
	  tre0_1_1 = tre1_0_0 - tre1_1_0;
	  tim0_1_1 = tim1_0_0 - tim1_1_0;
     }
     c_re(out[0]) = tre0_0_0 + tre0_0_1;
     c_im(out[0]) = tim0_0_0 + tim0_0_1;
     c_re(out[2 * ostride]) = tre0_0_0 - tre0_0_1;
     c_im(out[2 * ostride]) = tim0_0_0 - tim0_0_1;
     c_re(out[ostride]) = tre0_1_0 + tim0_1_1;
     c_im(out[ostride]) = tim0_1_0 - tre0_1_1;
     c_re(out[3 * ostride]) = tre0_1_0 - tim0_1_1;
     c_im(out[3 * ostride]) = tim0_1_0 + tre0_1_1;
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 44 FP additions and 16 FP multiplications */

void fftw_no_twiddle_5(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_2_0;
     FFTW_REAL tim0_2_0;
     FFTW_REAL tre0_3_0;
     FFTW_REAL tim0_3_0;
     FFTW_REAL tre0_4_0;
     FFTW_REAL tim0_4_0;
     tre0_0_0 = c_re(in[0]);
     tim0_0_0 = c_im(in[0]);
     tre0_1_0 = c_re(in[istride]);
     tim0_1_0 = c_im(in[istride]);
     tre0_2_0 = c_re(in[2 * istride]);
     tim0_2_0 = c_im(in[2 * istride]);
     tre0_3_0 = c_re(in[3 * istride]);
     tim0_3_0 = c_im(in[3 * istride]);
     tre0_4_0 = c_re(in[4 * istride]);
     tim0_4_0 = c_im(in[4 * istride]);
     c_re(out[0]) = tre0_0_0 + tre0_1_0 + tre0_2_0 + tre0_3_0 + tre0_4_0;
     c_im(out[0]) = tim0_0_0 + tim0_1_0 + tim0_2_0 + tim0_3_0 + tim0_4_0;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_1_0 + tre0_4_0)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_2_0 + tre0_3_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tim0_1_0 - tim0_4_0)) + (((FFTW_REAL) FFTW_K587785252) * (tim0_2_0 - tim0_3_0));
	  c_re(out[ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[4 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_1_0 + tim0_4_0)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_2_0 + tim0_3_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tre0_4_0 - tre0_1_0)) + (((FFTW_REAL) FFTW_K587785252) * (tre0_3_0 - tre0_2_0));
	  c_im(out[ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[4 * ostride]) = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_2_0 + tre0_3_0)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_1_0 + tre0_4_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tim0_1_0 - tim0_4_0)) + (((FFTW_REAL) FFTW_K951056516) * (tim0_3_0 - tim0_2_0));
	  c_re(out[2 * ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[3 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_2_0 + tim0_3_0)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_1_0 + tim0_4_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tre0_4_0 - tre0_1_0)) + (((FFTW_REAL) FFTW_K951056516) * (tre0_2_0 - tre0_3_0));
	  c_im(out[2 * ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[3 * ostride]) = tim1_0_0 - tim1_1_0;
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 40 FP additions and 8 FP multiplications */

void fftw_no_twiddle_6(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_0_1;
     FFTW_REAL tim0_0_1;
     FFTW_REAL tre0_0_2;
     FFTW_REAL tim0_0_2;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_1_1;
     FFTW_REAL tim0_1_1;
     FFTW_REAL tre0_1_2;
     FFTW_REAL tim0_1_2;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[0]);
	  tim1_0_0 = c_im(in[0]);
	  tre1_1_0 = c_re(in[3 * istride]);
	  tim1_1_0 = c_im(in[3 * istride]);
	  tre0_0_0 = tre1_0_0 + tre1_1_0;
	  tim0_0_0 = tim1_0_0 + tim1_1_0;
	  tre0_1_0 = tre1_0_0 - tre1_1_0;
	  tim0_1_0 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[2 * istride]);
	  tim1_0_0 = c_im(in[2 * istride]);
	  tre1_1_0 = c_re(in[5 * istride]);
	  tim1_1_0 = c_im(in[5 * istride]);
	  tre0_0_1 = tre1_0_0 + tre1_1_0;
	  tim0_0_1 = tim1_0_0 + tim1_1_0;
	  tre0_1_1 = tre1_0_0 - tre1_1_0;
	  tim0_1_1 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[4 * istride]);
	  tim1_0_0 = c_im(in[4 * istride]);
	  tre1_1_0 = c_re(in[istride]);
	  tim1_1_0 = c_im(in[istride]);
	  tre0_0_2 = tre1_0_0 + tre1_1_0;
	  tim0_0_2 = tim1_0_0 + tim1_1_0;
	  tre0_1_2 = tre1_0_0 - tre1_1_0;
	  tim0_1_2 = tim1_0_0 - tim1_1_0;
     }
     c_re(out[0]) = tre0_0_0 + tre0_0_1 + tre0_0_2;
     c_im(out[0]) = tim0_0_0 + tim0_0_1 + tim0_0_2;
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre0_0_1 + tre0_0_2));
	  tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim0_0_1 - tim0_0_2);
	  c_re(out[4 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[2 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim0_0_1 + tim0_0_2));
	  tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre0_0_2 - tre0_0_1);
	  c_im(out[4 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[2 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     c_re(out[3 * ostride]) = tre0_1_0 + tre0_1_1 + tre0_1_2;
     c_im(out[3 * ostride]) = tim0_1_0 + tim0_1_1 + tim0_1_2;
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_1_0 - (((FFTW_REAL) FFTW_K499999999) * (tre0_1_1 + tre0_1_2));
	  tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim0_1_1 - tim0_1_2);
	  c_re(out[ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[5 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_1_0 - (((FFTW_REAL) FFTW_K499999999) * (tim0_1_1 + tim0_1_2));
	  tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre0_1_2 - tre0_1_1);
	  c_im(out[ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[5 * ostride]) = tim2_0_0 - tim2_1_0;
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 928 FP additions and 248 FP multiplications */

void fftw_no_twiddle_64(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_0_1;
     FFTW_REAL tim0_0_1;
     FFTW_REAL tre0_0_2;
     FFTW_REAL tim0_0_2;
     FFTW_REAL tre0_0_3;
     FFTW_REAL tim0_0_3;
     FFTW_REAL tre0_0_4;
     FFTW_REAL tim0_0_4;
     FFTW_REAL tre0_0_5;
     FFTW_REAL tim0_0_5;
     FFTW_REAL tre0_0_6;
     FFTW_REAL tim0_0_6;
     FFTW_REAL tre0_0_7;
     FFTW_REAL tim0_0_7;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_1_1;
     FFTW_REAL tim0_1_1;
     FFTW_REAL tre0_1_2;
     FFTW_REAL tim0_1_2;
     FFTW_REAL tre0_1_3;
     FFTW_REAL tim0_1_3;
     FFTW_REAL tre0_1_4;
     FFTW_REAL tim0_1_4;
     FFTW_REAL tre0_1_5;
     FFTW_REAL tim0_1_5;
     FFTW_REAL tre0_1_6;
     FFTW_REAL tim0_1_6;
     FFTW_REAL tre0_1_7;
     FFTW_REAL tim0_1_7;
     FFTW_REAL tre0_2_0;
     FFTW_REAL tim0_2_0;
     FFTW_REAL tre0_2_1;
     FFTW_REAL tim0_2_1;
     FFTW_REAL tre0_2_2;
     FFTW_REAL tim0_2_2;
     FFTW_REAL tre0_2_3;
     FFTW_REAL tim0_2_3;
     FFTW_REAL tre0_2_4;
     FFTW_REAL tim0_2_4;
     FFTW_REAL tre0_2_5;
     FFTW_REAL tim0_2_5;
     FFTW_REAL tre0_2_6;
     FFTW_REAL tim0_2_6;
     FFTW_REAL tre0_2_7;
     FFTW_REAL tim0_2_7;
     FFTW_REAL tre0_3_0;
     FFTW_REAL tim0_3_0;
     FFTW_REAL tre0_3_1;
     FFTW_REAL tim0_3_1;
     FFTW_REAL tre0_3_2;
     FFTW_REAL tim0_3_2;
     FFTW_REAL tre0_3_3;
     FFTW_REAL tim0_3_3;
     FFTW_REAL tre0_3_4;
     FFTW_REAL tim0_3_4;
     FFTW_REAL tre0_3_5;
     FFTW_REAL tim0_3_5;
     FFTW_REAL tre0_3_6;
     FFTW_REAL tim0_3_6;
     FFTW_REAL tre0_3_7;
     FFTW_REAL tim0_3_7;
     FFTW_REAL tre0_4_0;
     FFTW_REAL tim0_4_0;
     FFTW_REAL tre0_4_1;
     FFTW_REAL tim0_4_1;
     FFTW_REAL tre0_4_2;
     FFTW_REAL tim0_4_2;
     FFTW_REAL tre0_4_3;
     FFTW_REAL tim0_4_3;
     FFTW_REAL tre0_4_4;
     FFTW_REAL tim0_4_4;
     FFTW_REAL tre0_4_5;
     FFTW_REAL tim0_4_5;
     FFTW_REAL tre0_4_6;
     FFTW_REAL tim0_4_6;
     FFTW_REAL tre0_4_7;
     FFTW_REAL tim0_4_7;
     FFTW_REAL tre0_5_0;
     FFTW_REAL tim0_5_0;
     FFTW_REAL tre0_5_1;
     FFTW_REAL tim0_5_1;
     FFTW_REAL tre0_5_2;
     FFTW_REAL tim0_5_2;
     FFTW_REAL tre0_5_3;
     FFTW_REAL tim0_5_3;
     FFTW_REAL tre0_5_4;
     FFTW_REAL tim0_5_4;
     FFTW_REAL tre0_5_5;
     FFTW_REAL tim0_5_5;
     FFTW_REAL tre0_5_6;
     FFTW_REAL tim0_5_6;
     FFTW_REAL tre0_5_7;
     FFTW_REAL tim0_5_7;
     FFTW_REAL tre0_6_0;
     FFTW_REAL tim0_6_0;
     FFTW_REAL tre0_6_1;
     FFTW_REAL tim0_6_1;
     FFTW_REAL tre0_6_2;
     FFTW_REAL tim0_6_2;
     FFTW_REAL tre0_6_3;
     FFTW_REAL tim0_6_3;
     FFTW_REAL tre0_6_4;
     FFTW_REAL tim0_6_4;
     FFTW_REAL tre0_6_5;
     FFTW_REAL tim0_6_5;
     FFTW_REAL tre0_6_6;
     FFTW_REAL tim0_6_6;
     FFTW_REAL tre0_6_7;
     FFTW_REAL tim0_6_7;
     FFTW_REAL tre0_7_0;
     FFTW_REAL tim0_7_0;
     FFTW_REAL tre0_7_1;
     FFTW_REAL tim0_7_1;
     FFTW_REAL tre0_7_2;
     FFTW_REAL tim0_7_2;
     FFTW_REAL tre0_7_3;
     FFTW_REAL tim0_7_3;
     FFTW_REAL tre0_7_4;
     FFTW_REAL tim0_7_4;
     FFTW_REAL tre0_7_5;
     FFTW_REAL tim0_7_5;
     FFTW_REAL tre0_7_6;
     FFTW_REAL tim0_7_6;
     FFTW_REAL tre0_7_7;
     FFTW_REAL tim0_7_7;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[0]);
	       tim2_0_0 = c_im(in[0]);
	       tre2_1_0 = c_re(in[32 * istride]);
	       tim2_1_0 = c_im(in[32 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[8 * istride]);
	       tim2_0_0 = c_im(in[8 * istride]);
	       tre2_1_0 = c_re(in[40 * istride]);
	       tim2_1_0 = c_im(in[40 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[16 * istride]);
	       tim2_0_0 = c_im(in[16 * istride]);
	       tre2_1_0 = c_re(in[48 * istride]);
	       tim2_1_0 = c_im(in[48 * istride]);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[24 * istride]);
	       tim2_0_0 = c_im(in[24 * istride]);
	       tre2_1_0 = c_re(in[56 * istride]);
	       tim2_1_0 = c_im(in[56 * istride]);
	       tre1_0_3 = tre2_0_0 + tre2_1_0;
	       tim1_0_3 = tim2_0_0 + tim2_1_0;
	       tre1_1_3 = tre2_0_0 - tre2_1_0;
	       tim1_1_3 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       tre0_0_0 = tre2_0_0 + tre2_0_1;
	       tim0_0_0 = tim2_0_0 + tim2_0_1;
	       tre0_4_0 = tre2_0_0 - tre2_0_1;
	       tim0_4_0 = tim2_0_0 - tim2_0_1;
	       tre0_2_0 = tre2_1_0 + tim2_1_1;
	       tim0_2_0 = tim2_1_0 - tre2_1_1;
	       tre0_6_0 = tre2_1_0 - tim2_1_1;
	       tim0_6_0 = tim2_1_0 + tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 + tim1_1_2;
	       tim2_0_0 = tim1_1_0 - tre1_1_2;
	       tre2_1_0 = tre1_1_0 - tim1_1_2;
	       tim2_1_0 = tim1_1_0 + tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
		    tre2_0_1 = tre3_0_0 + tre3_1_0;
		    tim2_0_1 = tim3_0_0 - tim3_1_0;
		    tre2_1_1 = tre3_0_0 - tre3_1_0;
		    tim2_1_1 = tim3_0_0 + tim3_1_0;
	       }
	       tre0_1_0 = tre2_0_0 + tre2_0_1;
	       tim0_1_0 = tim2_0_0 + tim2_0_1;
	       tre0_5_0 = tre2_0_0 - tre2_0_1;
	       tim0_5_0 = tim2_0_0 - tim2_0_1;
	       tre0_3_0 = tre2_1_0 + tim2_1_1;
	       tim0_3_0 = tim2_1_0 - tre2_1_1;
	       tre0_7_0 = tre2_1_0 - tim2_1_1;
	       tim0_7_0 = tim2_1_0 + tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[istride]);
	       tim2_0_0 = c_im(in[istride]);
	       tre2_1_0 = c_re(in[33 * istride]);
	       tim2_1_0 = c_im(in[33 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[9 * istride]);
	       tim2_0_0 = c_im(in[9 * istride]);
	       tre2_1_0 = c_re(in[41 * istride]);
	       tim2_1_0 = c_im(in[41 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[17 * istride]);
	       tim2_0_0 = c_im(in[17 * istride]);
	       tre2_1_0 = c_re(in[49 * istride]);
	       tim2_1_0 = c_im(in[49 * istride]);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[25 * istride]);
	       tim2_0_0 = c_im(in[25 * istride]);
	       tre2_1_0 = c_re(in[57 * istride]);
	       tim2_1_0 = c_im(in[57 * istride]);
	       tre1_0_3 = tre2_0_0 + tre2_1_0;
	       tim1_0_3 = tim2_0_0 + tim2_1_0;
	       tre1_1_3 = tre2_0_0 - tre2_1_0;
	       tim1_1_3 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       tre0_0_1 = tre2_0_0 + tre2_0_1;
	       tim0_0_1 = tim2_0_0 + tim2_0_1;
	       tre0_4_1 = tre2_0_0 - tre2_0_1;
	       tim0_4_1 = tim2_0_0 - tim2_0_1;
	       tre0_2_1 = tre2_1_0 + tim2_1_1;
	       tim0_2_1 = tim2_1_0 - tre2_1_1;
	       tre0_6_1 = tre2_1_0 - tim2_1_1;
	       tim0_6_1 = tim2_1_0 + tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 + tim1_1_2;
	       tim2_0_0 = tim1_1_0 - tre1_1_2;
	       tre2_1_0 = tre1_1_0 - tim1_1_2;
	       tim2_1_0 = tim1_1_0 + tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
		    tre2_0_1 = tre3_0_0 + tre3_1_0;
		    tim2_0_1 = tim3_0_0 - tim3_1_0;
		    tre2_1_1 = tre3_0_0 - tre3_1_0;
		    tim2_1_1 = tim3_0_0 + tim3_1_0;
	       }
	       tre0_1_1 = tre2_0_0 + tre2_0_1;
	       tim0_1_1 = tim2_0_0 + tim2_0_1;
	       tre0_5_1 = tre2_0_0 - tre2_0_1;
	       tim0_5_1 = tim2_0_0 - tim2_0_1;
	       tre0_3_1 = tre2_1_0 + tim2_1_1;
	       tim0_3_1 = tim2_1_0 - tre2_1_1;
	       tre0_7_1 = tre2_1_0 - tim2_1_1;
	       tim0_7_1 = tim2_1_0 + tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[2 * istride]);
	       tim2_0_0 = c_im(in[2 * istride]);
	       tre2_1_0 = c_re(in[34 * istride]);
	       tim2_1_0 = c_im(in[34 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[10 * istride]);
	       tim2_0_0 = c_im(in[10 * istride]);
	       tre2_1_0 = c_re(in[42 * istride]);
	       tim2_1_0 = c_im(in[42 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[18 * istride]);
	       tim2_0_0 = c_im(in[18 * istride]);
	       tre2_1_0 = c_re(in[50 * istride]);
	       tim2_1_0 = c_im(in[50 * istride]);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[26 * istride]);
	       tim2_0_0 = c_im(in[26 * istride]);
	       tre2_1_0 = c_re(in[58 * istride]);
	       tim2_1_0 = c_im(in[58 * istride]);
	       tre1_0_3 = tre2_0_0 + tre2_1_0;
	       tim1_0_3 = tim2_0_0 + tim2_1_0;
	       tre1_1_3 = tre2_0_0 - tre2_1_0;
	       tim1_1_3 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       tre0_0_2 = tre2_0_0 + tre2_0_1;
	       tim0_0_2 = tim2_0_0 + tim2_0_1;
	       tre0_4_2 = tre2_0_0 - tre2_0_1;
	       tim0_4_2 = tim2_0_0 - tim2_0_1;
	       tre0_2_2 = tre2_1_0 + tim2_1_1;
	       tim0_2_2 = tim2_1_0 - tre2_1_1;
	       tre0_6_2 = tre2_1_0 - tim2_1_1;
	       tim0_6_2 = tim2_1_0 + tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 + tim1_1_2;
	       tim2_0_0 = tim1_1_0 - tre1_1_2;
	       tre2_1_0 = tre1_1_0 - tim1_1_2;
	       tim2_1_0 = tim1_1_0 + tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
		    tre2_0_1 = tre3_0_0 + tre3_1_0;
		    tim2_0_1 = tim3_0_0 - tim3_1_0;
		    tre2_1_1 = tre3_0_0 - tre3_1_0;
		    tim2_1_1 = tim3_0_0 + tim3_1_0;
	       }
	       tre0_1_2 = tre2_0_0 + tre2_0_1;
	       tim0_1_2 = tim2_0_0 + tim2_0_1;
	       tre0_5_2 = tre2_0_0 - tre2_0_1;
	       tim0_5_2 = tim2_0_0 - tim2_0_1;
	       tre0_3_2 = tre2_1_0 + tim2_1_1;
	       tim0_3_2 = tim2_1_0 - tre2_1_1;
	       tre0_7_2 = tre2_1_0 - tim2_1_1;
	       tim0_7_2 = tim2_1_0 + tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[3 * istride]);
	       tim2_0_0 = c_im(in[3 * istride]);
	       tre2_1_0 = c_re(in[35 * istride]);
	       tim2_1_0 = c_im(in[35 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[11 * istride]);
	       tim2_0_0 = c_im(in[11 * istride]);
	       tre2_1_0 = c_re(in[43 * istride]);
	       tim2_1_0 = c_im(in[43 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[19 * istride]);
	       tim2_0_0 = c_im(in[19 * istride]);
	       tre2_1_0 = c_re(in[51 * istride]);
	       tim2_1_0 = c_im(in[51 * istride]);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[27 * istride]);
	       tim2_0_0 = c_im(in[27 * istride]);
	       tre2_1_0 = c_re(in[59 * istride]);
	       tim2_1_0 = c_im(in[59 * istride]);
	       tre1_0_3 = tre2_0_0 + tre2_1_0;
	       tim1_0_3 = tim2_0_0 + tim2_1_0;
	       tre1_1_3 = tre2_0_0 - tre2_1_0;
	       tim1_1_3 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       tre0_0_3 = tre2_0_0 + tre2_0_1;
	       tim0_0_3 = tim2_0_0 + tim2_0_1;
	       tre0_4_3 = tre2_0_0 - tre2_0_1;
	       tim0_4_3 = tim2_0_0 - tim2_0_1;
	       tre0_2_3 = tre2_1_0 + tim2_1_1;
	       tim0_2_3 = tim2_1_0 - tre2_1_1;
	       tre0_6_3 = tre2_1_0 - tim2_1_1;
	       tim0_6_3 = tim2_1_0 + tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 + tim1_1_2;
	       tim2_0_0 = tim1_1_0 - tre1_1_2;
	       tre2_1_0 = tre1_1_0 - tim1_1_2;
	       tim2_1_0 = tim1_1_0 + tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
		    tre2_0_1 = tre3_0_0 + tre3_1_0;
		    tim2_0_1 = tim3_0_0 - tim3_1_0;
		    tre2_1_1 = tre3_0_0 - tre3_1_0;
		    tim2_1_1 = tim3_0_0 + tim3_1_0;
	       }
	       tre0_1_3 = tre2_0_0 + tre2_0_1;
	       tim0_1_3 = tim2_0_0 + tim2_0_1;
	       tre0_5_3 = tre2_0_0 - tre2_0_1;
	       tim0_5_3 = tim2_0_0 - tim2_0_1;
	       tre0_3_3 = tre2_1_0 + tim2_1_1;
	       tim0_3_3 = tim2_1_0 - tre2_1_1;
	       tre0_7_3 = tre2_1_0 - tim2_1_1;
	       tim0_7_3 = tim2_1_0 + tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[4 * istride]);
	       tim2_0_0 = c_im(in[4 * istride]);
	       tre2_1_0 = c_re(in[36 * istride]);
	       tim2_1_0 = c_im(in[36 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[12 * istride]);
	       tim2_0_0 = c_im(in[12 * istride]);
	       tre2_1_0 = c_re(in[44 * istride]);
	       tim2_1_0 = c_im(in[44 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[20 * istride]);
	       tim2_0_0 = c_im(in[20 * istride]);
	       tre2_1_0 = c_re(in[52 * istride]);
	       tim2_1_0 = c_im(in[52 * istride]);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[28 * istride]);
	       tim2_0_0 = c_im(in[28 * istride]);
	       tre2_1_0 = c_re(in[60 * istride]);
	       tim2_1_0 = c_im(in[60 * istride]);
	       tre1_0_3 = tre2_0_0 + tre2_1_0;
	       tim1_0_3 = tim2_0_0 + tim2_1_0;
	       tre1_1_3 = tre2_0_0 - tre2_1_0;
	       tim1_1_3 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       tre0_0_4 = tre2_0_0 + tre2_0_1;
	       tim0_0_4 = tim2_0_0 + tim2_0_1;
	       tre0_4_4 = tre2_0_0 - tre2_0_1;
	       tim0_4_4 = tim2_0_0 - tim2_0_1;
	       tre0_2_4 = tre2_1_0 + tim2_1_1;
	       tim0_2_4 = tim2_1_0 - tre2_1_1;
	       tre0_6_4 = tre2_1_0 - tim2_1_1;
	       tim0_6_4 = tim2_1_0 + tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 + tim1_1_2;
	       tim2_0_0 = tim1_1_0 - tre1_1_2;
	       tre2_1_0 = tre1_1_0 - tim1_1_2;
	       tim2_1_0 = tim1_1_0 + tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
		    tre2_0_1 = tre3_0_0 + tre3_1_0;
		    tim2_0_1 = tim3_0_0 - tim3_1_0;
		    tre2_1_1 = tre3_0_0 - tre3_1_0;
		    tim2_1_1 = tim3_0_0 + tim3_1_0;
	       }
	       tre0_1_4 = tre2_0_0 + tre2_0_1;
	       tim0_1_4 = tim2_0_0 + tim2_0_1;
	       tre0_5_4 = tre2_0_0 - tre2_0_1;
	       tim0_5_4 = tim2_0_0 - tim2_0_1;
	       tre0_3_4 = tre2_1_0 + tim2_1_1;
	       tim0_3_4 = tim2_1_0 - tre2_1_1;
	       tre0_7_4 = tre2_1_0 - tim2_1_1;
	       tim0_7_4 = tim2_1_0 + tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[5 * istride]);
	       tim2_0_0 = c_im(in[5 * istride]);
	       tre2_1_0 = c_re(in[37 * istride]);
	       tim2_1_0 = c_im(in[37 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[13 * istride]);
	       tim2_0_0 = c_im(in[13 * istride]);
	       tre2_1_0 = c_re(in[45 * istride]);
	       tim2_1_0 = c_im(in[45 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[21 * istride]);
	       tim2_0_0 = c_im(in[21 * istride]);
	       tre2_1_0 = c_re(in[53 * istride]);
	       tim2_1_0 = c_im(in[53 * istride]);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[29 * istride]);
	       tim2_0_0 = c_im(in[29 * istride]);
	       tre2_1_0 = c_re(in[61 * istride]);
	       tim2_1_0 = c_im(in[61 * istride]);
	       tre1_0_3 = tre2_0_0 + tre2_1_0;
	       tim1_0_3 = tim2_0_0 + tim2_1_0;
	       tre1_1_3 = tre2_0_0 - tre2_1_0;
	       tim1_1_3 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       tre0_0_5 = tre2_0_0 + tre2_0_1;
	       tim0_0_5 = tim2_0_0 + tim2_0_1;
	       tre0_4_5 = tre2_0_0 - tre2_0_1;
	       tim0_4_5 = tim2_0_0 - tim2_0_1;
	       tre0_2_5 = tre2_1_0 + tim2_1_1;
	       tim0_2_5 = tim2_1_0 - tre2_1_1;
	       tre0_6_5 = tre2_1_0 - tim2_1_1;
	       tim0_6_5 = tim2_1_0 + tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 + tim1_1_2;
	       tim2_0_0 = tim1_1_0 - tre1_1_2;
	       tre2_1_0 = tre1_1_0 - tim1_1_2;
	       tim2_1_0 = tim1_1_0 + tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
		    tre2_0_1 = tre3_0_0 + tre3_1_0;
		    tim2_0_1 = tim3_0_0 - tim3_1_0;
		    tre2_1_1 = tre3_0_0 - tre3_1_0;
		    tim2_1_1 = tim3_0_0 + tim3_1_0;
	       }
	       tre0_1_5 = tre2_0_0 + tre2_0_1;
	       tim0_1_5 = tim2_0_0 + tim2_0_1;
	       tre0_5_5 = tre2_0_0 - tre2_0_1;
	       tim0_5_5 = tim2_0_0 - tim2_0_1;
	       tre0_3_5 = tre2_1_0 + tim2_1_1;
	       tim0_3_5 = tim2_1_0 - tre2_1_1;
	       tre0_7_5 = tre2_1_0 - tim2_1_1;
	       tim0_7_5 = tim2_1_0 + tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[6 * istride]);
	       tim2_0_0 = c_im(in[6 * istride]);
	       tre2_1_0 = c_re(in[38 * istride]);
	       tim2_1_0 = c_im(in[38 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[14 * istride]);
	       tim2_0_0 = c_im(in[14 * istride]);
	       tre2_1_0 = c_re(in[46 * istride]);
	       tim2_1_0 = c_im(in[46 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[22 * istride]);
	       tim2_0_0 = c_im(in[22 * istride]);
	       tre2_1_0 = c_re(in[54 * istride]);
	       tim2_1_0 = c_im(in[54 * istride]);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[30 * istride]);
	       tim2_0_0 = c_im(in[30 * istride]);
	       tre2_1_0 = c_re(in[62 * istride]);
	       tim2_1_0 = c_im(in[62 * istride]);
	       tre1_0_3 = tre2_0_0 + tre2_1_0;
	       tim1_0_3 = tim2_0_0 + tim2_1_0;
	       tre1_1_3 = tre2_0_0 - tre2_1_0;
	       tim1_1_3 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       tre0_0_6 = tre2_0_0 + tre2_0_1;
	       tim0_0_6 = tim2_0_0 + tim2_0_1;
	       tre0_4_6 = tre2_0_0 - tre2_0_1;
	       tim0_4_6 = tim2_0_0 - tim2_0_1;
	       tre0_2_6 = tre2_1_0 + tim2_1_1;
	       tim0_2_6 = tim2_1_0 - tre2_1_1;
	       tre0_6_6 = tre2_1_0 - tim2_1_1;
	       tim0_6_6 = tim2_1_0 + tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 + tim1_1_2;
	       tim2_0_0 = tim1_1_0 - tre1_1_2;
	       tre2_1_0 = tre1_1_0 - tim1_1_2;
	       tim2_1_0 = tim1_1_0 + tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
		    tre2_0_1 = tre3_0_0 + tre3_1_0;
		    tim2_0_1 = tim3_0_0 - tim3_1_0;
		    tre2_1_1 = tre3_0_0 - tre3_1_0;
		    tim2_1_1 = tim3_0_0 + tim3_1_0;
	       }
	       tre0_1_6 = tre2_0_0 + tre2_0_1;
	       tim0_1_6 = tim2_0_0 + tim2_0_1;
	       tre0_5_6 = tre2_0_0 - tre2_0_1;
	       tim0_5_6 = tim2_0_0 - tim2_0_1;
	       tre0_3_6 = tre2_1_0 + tim2_1_1;
	       tim0_3_6 = tim2_1_0 - tre2_1_1;
	       tre0_7_6 = tre2_1_0 - tim2_1_1;
	       tim0_7_6 = tim2_1_0 + tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[7 * istride]);
	       tim2_0_0 = c_im(in[7 * istride]);
	       tre2_1_0 = c_re(in[39 * istride]);
	       tim2_1_0 = c_im(in[39 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[15 * istride]);
	       tim2_0_0 = c_im(in[15 * istride]);
	       tre2_1_0 = c_re(in[47 * istride]);
	       tim2_1_0 = c_im(in[47 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[23 * istride]);
	       tim2_0_0 = c_im(in[23 * istride]);
	       tre2_1_0 = c_re(in[55 * istride]);
	       tim2_1_0 = c_im(in[55 * istride]);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[31 * istride]);
	       tim2_0_0 = c_im(in[31 * istride]);
	       tre2_1_0 = c_re(in[63 * istride]);
	       tim2_1_0 = c_im(in[63 * istride]);
	       tre1_0_3 = tre2_0_0 + tre2_1_0;
	       tim1_0_3 = tim2_0_0 + tim2_1_0;
	       tre1_1_3 = tre2_0_0 - tre2_1_0;
	       tim1_1_3 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       tre0_0_7 = tre2_0_0 + tre2_0_1;
	       tim0_0_7 = tim2_0_0 + tim2_0_1;
	       tre0_4_7 = tre2_0_0 - tre2_0_1;
	       tim0_4_7 = tim2_0_0 - tim2_0_1;
	       tre0_2_7 = tre2_1_0 + tim2_1_1;
	       tim0_2_7 = tim2_1_0 - tre2_1_1;
	       tre0_6_7 = tre2_1_0 - tim2_1_1;
	       tim0_6_7 = tim2_1_0 + tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 + tim1_1_2;
	       tim2_0_0 = tim1_1_0 - tre1_1_2;
	       tre2_1_0 = tre1_1_0 - tim1_1_2;
	       tim2_1_0 = tim1_1_0 + tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
		    tre2_0_1 = tre3_0_0 + tre3_1_0;
		    tim2_0_1 = tim3_0_0 - tim3_1_0;
		    tre2_1_1 = tre3_0_0 - tre3_1_0;
		    tim2_1_1 = tim3_0_0 + tim3_1_0;
	       }
	       tre0_1_7 = tre2_0_0 + tre2_0_1;
	       tim0_1_7 = tim2_0_0 + tim2_0_1;
	       tre0_5_7 = tre2_0_0 - tre2_0_1;
	       tim0_5_7 = tim2_0_0 - tim2_0_1;
	       tre0_3_7 = tre2_1_0 + tim2_1_1;
	       tim0_3_7 = tim2_1_0 - tre2_1_1;
	       tre0_7_7 = tre2_1_0 - tim2_1_1;
	       tim0_7_7 = tim2_1_0 + tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  tre1_0_0 = tre0_0_0 + tre0_0_4;
	  tim1_0_0 = tim0_0_0 + tim0_0_4;
	  tre1_1_0 = tre0_0_0 - tre0_0_4;
	  tim1_1_0 = tim0_0_0 - tim0_0_4;
	  tre1_0_1 = tre0_0_1 + tre0_0_5;
	  tim1_0_1 = tim0_0_1 + tim0_0_5;
	  tre1_1_1 = tre0_0_1 - tre0_0_5;
	  tim1_1_1 = tim0_0_1 - tim0_0_5;
	  tre1_0_2 = tre0_0_2 + tre0_0_6;
	  tim1_0_2 = tim0_0_2 + tim0_0_6;
	  tre1_1_2 = tre0_0_2 - tre0_0_6;
	  tim1_1_2 = tim0_0_2 - tim0_0_6;
	  tre1_0_3 = tre0_0_3 + tre0_0_7;
	  tim1_0_3 = tim0_0_3 + tim0_0_7;
	  tre1_1_3 = tre0_0_3 - tre0_0_7;
	  tim1_1_3 = tim0_0_3 - tim0_0_7;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       c_re(out[0]) = tre2_0_0 + tre2_0_1;
	       c_im(out[0]) = tim2_0_0 + tim2_0_1;
	       c_re(out[32 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[32 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[16 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[16 * ostride]) = tim2_1_0 - tre2_1_1;
	       c_re(out[48 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[48 * ostride]) = tim2_1_0 + tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 + tim1_1_2;
	       tim2_0_0 = tim1_1_0 - tre1_1_2;
	       tre2_1_0 = tre1_1_0 - tim1_1_2;
	       tim2_1_0 = tim1_1_0 + tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
		    tre2_0_1 = tre3_0_0 + tre3_1_0;
		    tim2_0_1 = tim3_0_0 - tim3_1_0;
		    tre2_1_1 = tre3_0_0 - tre3_1_0;
		    tim2_1_1 = tim3_0_0 + tim3_1_0;
	       }
	       c_re(out[8 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[8 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[40 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[40 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[24 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[24 * ostride]) = tim2_1_0 - tre2_1_1;
	       c_re(out[56 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[56 * ostride]) = tim2_1_0 + tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_1_4) + (((FFTW_REAL) FFTW_K382683432) * tim0_1_4);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_1_4) - (((FFTW_REAL) FFTW_K382683432) * tre0_1_4);
	       tre1_0_0 = tre0_1_0 + tre2_1_0;
	       tim1_0_0 = tim0_1_0 + tim2_1_0;
	       tre1_1_0 = tre0_1_0 - tre2_1_0;
	       tim1_1_0 = tim0_1_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K995184726) * tre0_1_1) + (((FFTW_REAL) FFTW_K098017140) * tim0_1_1);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K995184726) * tim0_1_1) - (((FFTW_REAL) FFTW_K098017140) * tre0_1_1);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K881921264) * tre0_1_5) + (((FFTW_REAL) FFTW_K471396736) * tim0_1_5);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K881921264) * tim0_1_5) - (((FFTW_REAL) FFTW_K471396736) * tre0_1_5);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tre0_1_2) + (((FFTW_REAL) FFTW_K195090322) * tim0_1_2);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tim0_1_2) - (((FFTW_REAL) FFTW_K195090322) * tre0_1_2);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_1_6) + (((FFTW_REAL) FFTW_K555570233) * tim0_1_6);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_1_6) - (((FFTW_REAL) FFTW_K555570233) * tre0_1_6);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K956940335) * tre0_1_3) + (((FFTW_REAL) FFTW_K290284677) * tim0_1_3);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K956940335) * tim0_1_3) - (((FFTW_REAL) FFTW_K290284677) * tre0_1_3);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K773010453) * tre0_1_7) + (((FFTW_REAL) FFTW_K634393284) * tim0_1_7);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K773010453) * tim0_1_7) - (((FFTW_REAL) FFTW_K634393284) * tre0_1_7);
	       tre1_0_3 = tre2_0_0 + tre2_1_0;
	       tim1_0_3 = tim2_0_0 + tim2_1_0;
	       tre1_1_3 = tre2_0_0 - tre2_1_0;
	       tim1_1_3 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       c_re(out[ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[33 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[33 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[17 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[17 * ostride]) = tim2_1_0 - tre2_1_1;
	       c_re(out[49 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[49 * ostride]) = tim2_1_0 + tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 + tim1_1_2;
	       tim2_0_0 = tim1_1_0 - tre1_1_2;
	       tre2_1_0 = tre1_1_0 - tim1_1_2;
	       tim2_1_0 = tim1_1_0 + tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
		    tre2_0_1 = tre3_0_0 + tre3_1_0;
		    tim2_0_1 = tim3_0_0 - tim3_1_0;
		    tre2_1_1 = tre3_0_0 - tre3_1_0;
		    tim2_1_1 = tim3_0_0 + tim3_1_0;
	       }
	       c_re(out[9 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[9 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[41 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[41 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[25 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[25 * ostride]) = tim2_1_0 - tre2_1_1;
	       c_re(out[57 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[57 * ostride]) = tim2_1_0 + tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_2_4 + tim0_2_4);
	       tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_2_4 - tre0_2_4);
	       tre1_0_0 = tre0_2_0 + tre2_1_0;
	       tim1_0_0 = tim0_2_0 + tim2_1_0;
	       tre1_1_0 = tre0_2_0 - tre2_1_0;
	       tim1_1_0 = tim0_2_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tre0_2_1) + (((FFTW_REAL) FFTW_K195090322) * tim0_2_1);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tim0_2_1) - (((FFTW_REAL) FFTW_K195090322) * tre0_2_1);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K555570233) * tre0_2_5) + (((FFTW_REAL) FFTW_K831469612) * tim0_2_5);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K555570233) * tim0_2_5) - (((FFTW_REAL) FFTW_K831469612) * tre0_2_5);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_2_2) + (((FFTW_REAL) FFTW_K382683432) * tim0_2_2);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_2_2) - (((FFTW_REAL) FFTW_K382683432) * tre0_2_2);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_2_6) + (((FFTW_REAL) FFTW_K923879532) * tim0_2_6);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_2_6) - (((FFTW_REAL) FFTW_K923879532) * tre0_2_6);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_2_3) + (((FFTW_REAL) FFTW_K555570233) * tim0_2_3);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_2_3) - (((FFTW_REAL) FFTW_K555570233) * tre0_2_3);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tre0_2_7) + (((FFTW_REAL) FFTW_K980785280) * tim0_2_7);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tim0_2_7) - (((FFTW_REAL) FFTW_K980785280) * tre0_2_7);
	       tre1_0_3 = tre2_0_0 + tre2_1_0;
	       tim1_0_3 = tim2_0_0 + tim2_1_0;
	       tre1_1_3 = tre2_0_0 - tre2_1_0;
	       tim1_1_3 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       c_re(out[2 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[2 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[34 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[34 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[18 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[18 * ostride]) = tim2_1_0 - tre2_1_1;
	       c_re(out[50 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[50 * ostride]) = tim2_1_0 + tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 + tim1_1_2;
	       tim2_0_0 = tim1_1_0 - tre1_1_2;
	       tre2_1_0 = tre1_1_0 - tim1_1_2;
	       tim2_1_0 = tim1_1_0 + tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
		    tre2_0_1 = tre3_0_0 + tre3_1_0;
		    tim2_0_1 = tim3_0_0 - tim3_1_0;
		    tre2_1_1 = tre3_0_0 - tre3_1_0;
		    tim2_1_1 = tim3_0_0 + tim3_1_0;
	       }
	       c_re(out[10 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[10 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[42 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[42 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[26 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[26 * ostride]) = tim2_1_0 - tre2_1_1;
	       c_re(out[58 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[58 * ostride]) = tim2_1_0 + tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_3_4) + (((FFTW_REAL) FFTW_K923879532) * tim0_3_4);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_3_4) - (((FFTW_REAL) FFTW_K923879532) * tre0_3_4);
	       tre1_0_0 = tre0_3_0 + tre2_1_0;
	       tim1_0_0 = tim0_3_0 + tim2_1_0;
	       tre1_1_0 = tre0_3_0 - tre2_1_0;
	       tim1_1_0 = tim0_3_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K956940335) * tre0_3_1) + (((FFTW_REAL) FFTW_K290284677) * tim0_3_1);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K956940335) * tim0_3_1) - (((FFTW_REAL) FFTW_K290284677) * tre0_3_1);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K098017140) * tre0_3_5) + (((FFTW_REAL) FFTW_K995184726) * tim0_3_5);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K098017140) * tim0_3_5) - (((FFTW_REAL) FFTW_K995184726) * tre0_3_5);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_3_2) + (((FFTW_REAL) FFTW_K555570233) * tim0_3_2);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_3_2) - (((FFTW_REAL) FFTW_K555570233) * tre0_3_2);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K980785280) * tim0_3_6) - (((FFTW_REAL) FFTW_K195090322) * tre0_3_6);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tim0_3_6) + (((FFTW_REAL) FFTW_K980785280) * tre0_3_6);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 - tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 + tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K634393284) * tre0_3_3) + (((FFTW_REAL) FFTW_K773010453) * tim0_3_3);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K634393284) * tim0_3_3) - (((FFTW_REAL) FFTW_K773010453) * tre0_3_3);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K881921264) * tim0_3_7) - (((FFTW_REAL) FFTW_K471396736) * tre0_3_7);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K471396736) * tim0_3_7) + (((FFTW_REAL) FFTW_K881921264) * tre0_3_7);
	       tre1_0_3 = tre2_0_0 + tre2_1_0;
	       tim1_0_3 = tim2_0_0 - tim2_1_0;
	       tre1_1_3 = tre2_0_0 - tre2_1_0;
	       tim1_1_3 = tim2_0_0 + tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       c_re(out[3 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[3 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[35 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[35 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[19 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[19 * ostride]) = tim2_1_0 - tre2_1_1;
	       c_re(out[51 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[51 * ostride]) = tim2_1_0 + tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 + tim1_1_2;
	       tim2_0_0 = tim1_1_0 - tre1_1_2;
	       tre2_1_0 = tre1_1_0 - tim1_1_2;
	       tim2_1_0 = tim1_1_0 + tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
		    tre2_0_1 = tre3_0_0 + tre3_1_0;
		    tim2_0_1 = tim3_0_0 - tim3_1_0;
		    tre2_1_1 = tre3_0_0 - tre3_1_0;
		    tim2_1_1 = tim3_0_0 + tim3_1_0;
	       }
	       c_re(out[11 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[11 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[43 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[43 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[27 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[27 * ostride]) = tim2_1_0 - tre2_1_1;
	       c_re(out[59 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[59 * ostride]) = tim2_1_0 + tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  tre1_0_0 = tre0_4_0 + tim0_4_4;
	  tim1_0_0 = tim0_4_0 - tre0_4_4;
	  tre1_1_0 = tre0_4_0 - tim0_4_4;
	  tim1_1_0 = tim0_4_0 + tre0_4_4;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_4_1) + (((FFTW_REAL) FFTW_K382683432) * tim0_4_1);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_4_1) - (((FFTW_REAL) FFTW_K382683432) * tre0_4_1);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_4_5) - (((FFTW_REAL) FFTW_K382683432) * tre0_4_5);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_4_5) + (((FFTW_REAL) FFTW_K923879532) * tre0_4_5);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 - tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 + tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_4_2 + tim0_4_2);
	       tim2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_4_2 - tre0_4_2);
	       tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_4_6 - tre0_4_6);
	       tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_4_6 + tre0_4_6);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 - tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 + tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_4_3) + (((FFTW_REAL) FFTW_K923879532) * tim0_4_3);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_4_3) - (((FFTW_REAL) FFTW_K923879532) * tre0_4_3);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_4_7) - (((FFTW_REAL) FFTW_K923879532) * tre0_4_7);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_4_7) + (((FFTW_REAL) FFTW_K382683432) * tre0_4_7);
	       tre1_0_3 = tre2_0_0 + tre2_1_0;
	       tim1_0_3 = tim2_0_0 - tim2_1_0;
	       tre1_1_3 = tre2_0_0 - tre2_1_0;
	       tim1_1_3 = tim2_0_0 + tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       c_re(out[4 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[4 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[36 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[36 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[20 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[20 * ostride]) = tim2_1_0 - tre2_1_1;
	       c_re(out[52 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[52 * ostride]) = tim2_1_0 + tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 + tim1_1_2;
	       tim2_0_0 = tim1_1_0 - tre1_1_2;
	       tre2_1_0 = tre1_1_0 - tim1_1_2;
	       tim2_1_0 = tim1_1_0 + tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
		    tre2_0_1 = tre3_0_0 + tre3_1_0;
		    tim2_0_1 = tim3_0_0 - tim3_1_0;
		    tre2_1_1 = tre3_0_0 - tre3_1_0;
		    tim2_1_1 = tim3_0_0 + tim3_1_0;
	       }
	       c_re(out[12 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[12 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[44 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[44 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[28 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[28 * ostride]) = tim2_1_0 - tre2_1_1;
	       c_re(out[60 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[60 * ostride]) = tim2_1_0 + tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_5_4) - (((FFTW_REAL) FFTW_K382683432) * tre0_5_4);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_5_4) + (((FFTW_REAL) FFTW_K923879532) * tre0_5_4);
	       tre1_0_0 = tre0_5_0 + tre2_1_0;
	       tim1_0_0 = tim0_5_0 - tim2_1_0;
	       tre1_1_0 = tre0_5_0 - tre2_1_0;
	       tim1_1_0 = tim0_5_0 + tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K881921264) * tre0_5_1) + (((FFTW_REAL) FFTW_K471396736) * tim0_5_1);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K881921264) * tim0_5_1) - (((FFTW_REAL) FFTW_K471396736) * tre0_5_1);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K634393284) * tim0_5_5) - (((FFTW_REAL) FFTW_K773010453) * tre0_5_5);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K773010453) * tim0_5_5) + (((FFTW_REAL) FFTW_K634393284) * tre0_5_5);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 - tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 + tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K555570233) * tre0_5_2) + (((FFTW_REAL) FFTW_K831469612) * tim0_5_2);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K555570233) * tim0_5_2) - (((FFTW_REAL) FFTW_K831469612) * tre0_5_2);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tim0_5_6) - (((FFTW_REAL) FFTW_K980785280) * tre0_5_6);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K980785280) * tim0_5_6) + (((FFTW_REAL) FFTW_K195090322) * tre0_5_6);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 - tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 + tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K098017140) * tre0_5_3) + (((FFTW_REAL) FFTW_K995184726) * tim0_5_3);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K098017140) * tim0_5_3) - (((FFTW_REAL) FFTW_K995184726) * tre0_5_3);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K956940335) * tre0_5_7) + (((FFTW_REAL) FFTW_K290284677) * tim0_5_7);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K290284677) * tre0_5_7) - (((FFTW_REAL) FFTW_K956940335) * tim0_5_7);
	       tre1_0_3 = tre2_0_0 - tre2_1_0;
	       tim1_0_3 = tim2_0_0 + tim2_1_0;
	       tre1_1_3 = tre2_0_0 + tre2_1_0;
	       tim1_1_3 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       c_re(out[5 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[5 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[37 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[37 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[21 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[21 * ostride]) = tim2_1_0 - tre2_1_1;
	       c_re(out[53 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[53 * ostride]) = tim2_1_0 + tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 + tim1_1_2;
	       tim2_0_0 = tim1_1_0 - tre1_1_2;
	       tre2_1_0 = tre1_1_0 - tim1_1_2;
	       tim2_1_0 = tim1_1_0 + tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
		    tre2_0_1 = tre3_0_0 + tre3_1_0;
		    tim2_0_1 = tim3_0_0 - tim3_1_0;
		    tre2_1_1 = tre3_0_0 - tre3_1_0;
		    tim2_1_1 = tim3_0_0 + tim3_1_0;
	       }
	       c_re(out[13 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[13 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[45 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[45 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[29 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[29 * ostride]) = tim2_1_0 - tre2_1_1;
	       c_re(out[61 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[61 * ostride]) = tim2_1_0 + tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_6_4 - tre0_6_4);
	       tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_6_4 + tre0_6_4);
	       tre1_0_0 = tre0_6_0 + tre2_1_0;
	       tim1_0_0 = tim0_6_0 - tim2_1_0;
	       tre1_1_0 = tre0_6_0 - tre2_1_0;
	       tim1_1_0 = tim0_6_0 + tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_6_1) + (((FFTW_REAL) FFTW_K555570233) * tim0_6_1);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_6_1) - (((FFTW_REAL) FFTW_K555570233) * tre0_6_1);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tim0_6_5) - (((FFTW_REAL) FFTW_K980785280) * tre0_6_5);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K980785280) * tim0_6_5) + (((FFTW_REAL) FFTW_K195090322) * tre0_6_5);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 - tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 + tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_6_2) + (((FFTW_REAL) FFTW_K923879532) * tim0_6_2);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_6_2) - (((FFTW_REAL) FFTW_K923879532) * tre0_6_2);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_6_6) + (((FFTW_REAL) FFTW_K382683432) * tim0_6_6);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_6_6) - (((FFTW_REAL) FFTW_K923879532) * tim0_6_6);
	       tre1_0_2 = tre2_0_0 - tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 + tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tim0_6_3) - (((FFTW_REAL) FFTW_K195090322) * tre0_6_3);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K195090322) * tim0_6_3) + (((FFTW_REAL) FFTW_K980785280) * tre0_6_3);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K555570233) * tre0_6_7) + (((FFTW_REAL) FFTW_K831469612) * tim0_6_7);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_6_7) - (((FFTW_REAL) FFTW_K555570233) * tim0_6_7);
	       tre1_0_3 = tre2_0_0 - tre2_1_0;
	       tim1_0_3 = tim2_1_0 - tim2_0_0;
	       tre1_1_3 = tre2_0_0 + tre2_1_0;
	       tim1_1_3 = (-(tim2_0_0 + tim2_1_0));
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       c_re(out[6 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[6 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[38 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[38 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[22 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[22 * ostride]) = tim2_1_0 - tre2_1_1;
	       c_re(out[54 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[54 * ostride]) = tim2_1_0 + tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 + tim1_1_2;
	       tim2_0_0 = tim1_1_0 - tre1_1_2;
	       tre2_1_0 = tre1_1_0 - tim1_1_2;
	       tim2_1_0 = tim1_1_0 + tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
		    tre2_0_1 = tre3_0_0 + tre3_1_0;
		    tim2_0_1 = tim3_0_0 - tim3_1_0;
		    tre2_1_1 = tre3_0_0 - tre3_1_0;
		    tim2_1_1 = tim3_0_0 + tim3_1_0;
	       }
	       c_re(out[14 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[14 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[46 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[46 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[30 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[30 * ostride]) = tim2_1_0 - tre2_1_1;
	       c_re(out[62 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[62 * ostride]) = tim2_1_0 + tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_7_4) - (((FFTW_REAL) FFTW_K923879532) * tre0_7_4);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_7_4) + (((FFTW_REAL) FFTW_K382683432) * tre0_7_4);
	       tre1_0_0 = tre0_7_0 + tre2_1_0;
	       tim1_0_0 = tim0_7_0 - tim2_1_0;
	       tre1_1_0 = tre0_7_0 - tre2_1_0;
	       tim1_1_0 = tim0_7_0 + tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K773010453) * tre0_7_1) + (((FFTW_REAL) FFTW_K634393284) * tim0_7_1);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K773010453) * tim0_7_1) - (((FFTW_REAL) FFTW_K634393284) * tre0_7_1);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K956940335) * tre0_7_5) + (((FFTW_REAL) FFTW_K290284677) * tim0_7_5);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K290284677) * tre0_7_5) - (((FFTW_REAL) FFTW_K956940335) * tim0_7_5);
	       tre1_0_1 = tre2_0_0 - tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 + tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K195090322) * tre0_7_2) + (((FFTW_REAL) FFTW_K980785280) * tim0_7_2);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K195090322) * tim0_7_2) - (((FFTW_REAL) FFTW_K980785280) * tre0_7_2);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K555570233) * tre0_7_6) + (((FFTW_REAL) FFTW_K831469612) * tim0_7_6);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_7_6) - (((FFTW_REAL) FFTW_K555570233) * tim0_7_6);
	       tre1_0_2 = tre2_0_0 - tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 + tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K881921264) * tim0_7_3) - (((FFTW_REAL) FFTW_K471396736) * tre0_7_3);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K471396736) * tim0_7_3) + (((FFTW_REAL) FFTW_K881921264) * tre0_7_3);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K098017140) * tre0_7_7) - (((FFTW_REAL) FFTW_K995184726) * tim0_7_7);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K098017140) * tim0_7_7) + (((FFTW_REAL) FFTW_K995184726) * tre0_7_7);
	       tre1_0_3 = tre2_0_0 + tre2_1_0;
	       tim1_0_3 = tim2_1_0 - tim2_0_0;
	       tre1_1_3 = tre2_0_0 - tre2_1_0;
	       tim1_1_3 = (-(tim2_0_0 + tim2_1_0));
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       c_re(out[7 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[7 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[39 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[39 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[23 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[23 * ostride]) = tim2_1_0 - tre2_1_1;
	       c_re(out[55 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[55 * ostride]) = tim2_1_0 + tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 + tim1_1_2;
	       tim2_0_0 = tim1_1_0 - tre1_1_2;
	       tre2_1_0 = tre1_1_0 - tim1_1_2;
	       tim2_1_0 = tim1_1_0 + tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
		    tre2_0_1 = tre3_0_0 + tre3_1_0;
		    tim2_0_1 = tim3_0_0 - tim3_1_0;
		    tre2_1_1 = tre3_0_0 - tre3_1_0;
		    tim2_1_1 = tim3_0_0 + tim3_1_0;
	       }
	       c_re(out[15 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[15 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[47 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[47 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[31 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[31 * ostride]) = tim2_1_0 - tre2_1_1;
	       c_re(out[63 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[63 * ostride]) = tim2_1_0 + tre2_1_1;
	  }
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 90 FP additions and 36 FP multiplications */

void fftw_no_twiddle_7(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_2_0;
     FFTW_REAL tim0_2_0;
     FFTW_REAL tre0_3_0;
     FFTW_REAL tim0_3_0;
     FFTW_REAL tre0_4_0;
     FFTW_REAL tim0_4_0;
     FFTW_REAL tre0_5_0;
     FFTW_REAL tim0_5_0;
     FFTW_REAL tre0_6_0;
     FFTW_REAL tim0_6_0;
     tre0_0_0 = c_re(in[0]);
     tim0_0_0 = c_im(in[0]);
     tre0_1_0 = c_re(in[istride]);
     tim0_1_0 = c_im(in[istride]);
     tre0_2_0 = c_re(in[2 * istride]);
     tim0_2_0 = c_im(in[2 * istride]);
     tre0_3_0 = c_re(in[3 * istride]);
     tim0_3_0 = c_im(in[3 * istride]);
     tre0_4_0 = c_re(in[4 * istride]);
     tim0_4_0 = c_im(in[4 * istride]);
     tre0_5_0 = c_re(in[5 * istride]);
     tim0_5_0 = c_im(in[5 * istride]);
     tre0_6_0 = c_re(in[6 * istride]);
     tim0_6_0 = c_im(in[6 * istride]);
     c_re(out[0]) = tre0_0_0 + tre0_1_0 + tre0_2_0 + tre0_3_0 + tre0_4_0 + tre0_5_0 + tre0_6_0;
     c_im(out[0]) = tim0_0_0 + tim0_1_0 + tim0_2_0 + tim0_3_0 + tim0_4_0 + tim0_5_0 + tim0_6_0;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tre0_1_0 + tre0_6_0)) - (((FFTW_REAL) FFTW_K900968867) * (tre0_3_0 + tre0_4_0)) - (((FFTW_REAL) FFTW_K222520933) * (tre0_2_0 + tre0_5_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K781831482) * (tim0_1_0 - tim0_6_0)) + (((FFTW_REAL) FFTW_K974927912) * (tim0_2_0 - tim0_5_0)) + (((FFTW_REAL) FFTW_K433883739) * (tim0_3_0 - tim0_4_0));
	  c_re(out[ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[6 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tim0_1_0 + tim0_6_0)) - (((FFTW_REAL) FFTW_K900968867) * (tim0_3_0 + tim0_4_0)) - (((FFTW_REAL) FFTW_K222520933) * (tim0_2_0 + tim0_5_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K781831482) * (tre0_6_0 - tre0_1_0)) + (((FFTW_REAL) FFTW_K974927912) * (tre0_5_0 - tre0_2_0)) + (((FFTW_REAL) FFTW_K433883739) * (tre0_4_0 - tre0_3_0));
	  c_im(out[ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[6 * ostride]) = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tre0_3_0 + tre0_4_0)) - (((FFTW_REAL) FFTW_K900968867) * (tre0_2_0 + tre0_5_0)) - (((FFTW_REAL) FFTW_K222520933) * (tre0_1_0 + tre0_6_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K974927912) * (tim0_1_0 - tim0_6_0)) + (((FFTW_REAL) FFTW_K433883739) * (tim0_5_0 - tim0_2_0)) + (((FFTW_REAL) FFTW_K781831482) * (tim0_4_0 - tim0_3_0));
	  c_re(out[2 * ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[5 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tim0_3_0 + tim0_4_0)) - (((FFTW_REAL) FFTW_K900968867) * (tim0_2_0 + tim0_5_0)) - (((FFTW_REAL) FFTW_K222520933) * (tim0_1_0 + tim0_6_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K974927912) * (tre0_6_0 - tre0_1_0)) + (((FFTW_REAL) FFTW_K433883739) * (tre0_2_0 - tre0_5_0)) + (((FFTW_REAL) FFTW_K781831482) * (tre0_3_0 - tre0_4_0));
	  c_im(out[2 * ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[5 * ostride]) = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tre0_2_0 + tre0_5_0)) - (((FFTW_REAL) FFTW_K222520933) * (tre0_3_0 + tre0_4_0)) - (((FFTW_REAL) FFTW_K900968867) * (tre0_1_0 + tre0_6_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K433883739) * (tim0_1_0 - tim0_6_0)) + (((FFTW_REAL) FFTW_K781831482) * (tim0_5_0 - tim0_2_0)) + (((FFTW_REAL) FFTW_K974927912) * (tim0_3_0 - tim0_4_0));
	  c_re(out[3 * ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[4 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tim0_2_0 + tim0_5_0)) - (((FFTW_REAL) FFTW_K222520933) * (tim0_3_0 + tim0_4_0)) - (((FFTW_REAL) FFTW_K900968867) * (tim0_1_0 + tim0_6_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K433883739) * (tre0_6_0 - tre0_1_0)) + (((FFTW_REAL) FFTW_K781831482) * (tre0_2_0 - tre0_5_0)) + (((FFTW_REAL) FFTW_K974927912) * (tre0_4_0 - tre0_3_0));
	  c_im(out[3 * ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[4 * ostride]) = tim1_0_0 - tim1_1_0;
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 52 FP additions and 4 FP multiplications */

void fftw_no_twiddle_8(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_0_1;
     FFTW_REAL tim0_0_1;
     FFTW_REAL tre0_0_2;
     FFTW_REAL tim0_0_2;
     FFTW_REAL tre0_0_3;
     FFTW_REAL tim0_0_3;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_1_1;
     FFTW_REAL tim0_1_1;
     FFTW_REAL tre0_1_2;
     FFTW_REAL tim0_1_2;
     FFTW_REAL tre0_1_3;
     FFTW_REAL tim0_1_3;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[0]);
	  tim1_0_0 = c_im(in[0]);
	  tre1_1_0 = c_re(in[4 * istride]);
	  tim1_1_0 = c_im(in[4 * istride]);
	  tre0_0_0 = tre1_0_0 + tre1_1_0;
	  tim0_0_0 = tim1_0_0 + tim1_1_0;
	  tre0_1_0 = tre1_0_0 - tre1_1_0;
	  tim0_1_0 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[istride]);
	  tim1_0_0 = c_im(in[istride]);
	  tre1_1_0 = c_re(in[5 * istride]);
	  tim1_1_0 = c_im(in[5 * istride]);
	  tre0_0_1 = tre1_0_0 + tre1_1_0;
	  tim0_0_1 = tim1_0_0 + tim1_1_0;
	  tre0_1_1 = tre1_0_0 - tre1_1_0;
	  tim0_1_1 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[2 * istride]);
	  tim1_0_0 = c_im(in[2 * istride]);
	  tre1_1_0 = c_re(in[6 * istride]);
	  tim1_1_0 = c_im(in[6 * istride]);
	  tre0_0_2 = tre1_0_0 + tre1_1_0;
	  tim0_0_2 = tim1_0_0 + tim1_1_0;
	  tre0_1_2 = tre1_0_0 - tre1_1_0;
	  tim0_1_2 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[3 * istride]);
	  tim1_0_0 = c_im(in[3 * istride]);
	  tre1_1_0 = c_re(in[7 * istride]);
	  tim1_1_0 = c_im(in[7 * istride]);
	  tre0_0_3 = tre1_0_0 + tre1_1_0;
	  tim0_0_3 = tim1_0_0 + tim1_1_0;
	  tre0_1_3 = tre1_0_0 - tre1_1_0;
	  tim0_1_3 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  tre1_0_0 = tre0_0_0 + tre0_0_2;
	  tim1_0_0 = tim0_0_0 + tim0_0_2;
	  tre1_1_0 = tre0_0_0 - tre0_0_2;
	  tim1_1_0 = tim0_0_0 - tim0_0_2;
	  tre1_0_1 = tre0_0_1 + tre0_0_3;
	  tim1_0_1 = tim0_0_1 + tim0_0_3;
	  tre1_1_1 = tre0_0_1 - tre0_0_3;
	  tim1_1_1 = tim0_0_1 - tim0_0_3;
	  c_re(out[0]) = tre1_0_0 + tre1_0_1;
	  c_im(out[0]) = tim1_0_0 + tim1_0_1;
	  c_re(out[4 * ostride]) = tre1_0_0 - tre1_0_1;
	  c_im(out[4 * ostride]) = tim1_0_0 - tim1_0_1;
	  c_re(out[2 * ostride]) = tre1_1_0 + tim1_1_1;
	  c_im(out[2 * ostride]) = tim1_1_0 - tre1_1_1;
	  c_re(out[6 * ostride]) = tre1_1_0 - tim1_1_1;
	  c_im(out[6 * ostride]) = tim1_1_0 + tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  tre1_0_0 = tre0_1_0 + tim0_1_2;
	  tim1_0_0 = tim0_1_0 - tre0_1_2;
	  tre1_1_0 = tre0_1_0 - tim0_1_2;
	  tim1_1_0 = tim0_1_0 + tre0_1_2;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_1_1 + tim0_1_1);
	       tim2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_1_1 - tre0_1_1);
	       tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_1_3 - tre0_1_3);
	       tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_1_3 + tre0_1_3);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 - tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 + tim2_1_0;
	  }
	  c_re(out[ostride]) = tre1_0_0 + tre1_0_1;
	  c_im(out[ostride]) = tim1_0_0 + tim1_0_1;
	  c_re(out[5 * ostride]) = tre1_0_0 - tre1_0_1;
	  c_im(out[5 * ostride]) = tim1_0_0 - tim1_0_1;
	  c_re(out[3 * ostride]) = tre1_1_0 + tim1_1_1;
	  c_im(out[3 * ostride]) = tim1_1_0 - tre1_1_1;
	  c_re(out[7 * ostride]) = tre1_1_0 - tim1_1_1;
	  c_im(out[7 * ostride]) = tim1_1_0 + tre1_1_1;
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 92 FP additions and 40 FP multiplications */

void fftw_no_twiddle_9(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_0_1;
     FFTW_REAL tim0_0_1;
     FFTW_REAL tre0_0_2;
     FFTW_REAL tim0_0_2;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_1_1;
     FFTW_REAL tim0_1_1;
     FFTW_REAL tre0_1_2;
     FFTW_REAL tim0_1_2;
     FFTW_REAL tre0_2_0;
     FFTW_REAL tim0_2_0;
     FFTW_REAL tre0_2_1;
     FFTW_REAL tim0_2_1;
     FFTW_REAL tre0_2_2;
     FFTW_REAL tim0_2_2;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_2_0;
	  FFTW_REAL tim1_2_0;
	  tre1_0_0 = c_re(in[0]);
	  tim1_0_0 = c_im(in[0]);
	  tre1_1_0 = c_re(in[3 * istride]);
	  tim1_1_0 = c_im(in[3 * istride]);
	  tre1_2_0 = c_re(in[6 * istride]);
	  tim1_2_0 = c_im(in[6 * istride]);
	  tre0_0_0 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	  tim0_0_0 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_1_0 - tim1_2_0);
	       tre0_1_0 = tre2_0_0 + tre2_1_0;
	       tre0_2_0 = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_2_0 - tre1_1_0);
	       tim0_1_0 = tim2_0_0 + tim2_1_0;
	       tim0_2_0 = tim2_0_0 - tim2_1_0;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_2_0;
	  FFTW_REAL tim1_2_0;
	  tre1_0_0 = c_re(in[istride]);
	  tim1_0_0 = c_im(in[istride]);
	  tre1_1_0 = c_re(in[4 * istride]);
	  tim1_1_0 = c_im(in[4 * istride]);
	  tre1_2_0 = c_re(in[7 * istride]);
	  tim1_2_0 = c_im(in[7 * istride]);
	  tre0_0_1 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	  tim0_0_1 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_1_0 - tim1_2_0);
	       tre0_1_1 = tre2_0_0 + tre2_1_0;
	       tre0_2_1 = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_2_0 - tre1_1_0);
	       tim0_1_1 = tim2_0_0 + tim2_1_0;
	       tim0_2_1 = tim2_0_0 - tim2_1_0;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_2_0;
	  FFTW_REAL tim1_2_0;
	  tre1_0_0 = c_re(in[2 * istride]);
	  tim1_0_0 = c_im(in[2 * istride]);
	  tre1_1_0 = c_re(in[5 * istride]);
	  tim1_1_0 = c_im(in[5 * istride]);
	  tre1_2_0 = c_re(in[8 * istride]);
	  tim1_2_0 = c_im(in[8 * istride]);
	  tre0_0_2 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	  tim0_0_2 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_1_0 - tim1_2_0);
	       tre0_1_2 = tre2_0_0 + tre2_1_0;
	       tre0_2_2 = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_2_0 - tre1_1_0);
	       tim0_1_2 = tim2_0_0 + tim2_1_0;
	       tim0_2_2 = tim2_0_0 - tim2_1_0;
	  }
     }
     c_re(out[0]) = tre0_0_0 + tre0_0_1 + tre0_0_2;
     c_im(out[0]) = tim0_0_0 + tim0_0_1 + tim0_0_2;
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre0_0_1 + tre0_0_2));
	  tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim0_0_1 - tim0_0_2);
	  c_re(out[3 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[6 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim0_0_1 + tim0_0_2));
	  tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre0_0_2 - tre0_0_1);
	  c_im(out[3 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[6 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     {
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_2_0;
	  FFTW_REAL tim1_2_0;
	  tre1_1_0 = (((FFTW_REAL) FFTW_K766044443) * tre0_1_1) + (((FFTW_REAL) FFTW_K642787609) * tim0_1_1);
	  tim1_1_0 = (((FFTW_REAL) FFTW_K766044443) * tim0_1_1) - (((FFTW_REAL) FFTW_K642787609) * tre0_1_1);
	  tre1_2_0 = (((FFTW_REAL) FFTW_K173648177) * tre0_1_2) + (((FFTW_REAL) FFTW_K984807753) * tim0_1_2);
	  tim1_2_0 = (((FFTW_REAL) FFTW_K173648177) * tim0_1_2) - (((FFTW_REAL) FFTW_K984807753) * tre0_1_2);
	  c_re(out[ostride]) = tre0_1_0 + tre1_1_0 + tre1_2_0;
	  c_im(out[ostride]) = tim0_1_0 + tim1_1_0 + tim1_2_0;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre0_1_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_1_0 - tim1_2_0);
	       c_re(out[4 * ostride]) = tre2_0_0 + tre2_1_0;
	       c_re(out[7 * ostride]) = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim0_1_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_2_0 - tre1_1_0);
	       c_im(out[4 * ostride]) = tim2_0_0 + tim2_1_0;
	       c_im(out[7 * ostride]) = tim2_0_0 - tim2_1_0;
	  }
     }
     {
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_2_0;
	  FFTW_REAL tim1_2_0;
	  tre1_1_0 = (((FFTW_REAL) FFTW_K173648177) * tre0_2_1) + (((FFTW_REAL) FFTW_K984807753) * tim0_2_1);
	  tim1_1_0 = (((FFTW_REAL) FFTW_K173648177) * tim0_2_1) - (((FFTW_REAL) FFTW_K984807753) * tre0_2_1);
	  tre1_2_0 = (((FFTW_REAL) FFTW_K342020143) * tim0_2_2) - (((FFTW_REAL) FFTW_K939692620) * tre0_2_2);
	  tim1_2_0 = (((FFTW_REAL) FFTW_K939692620) * tim0_2_2) + (((FFTW_REAL) FFTW_K342020143) * tre0_2_2);
	  c_re(out[2 * ostride]) = tre0_2_0 + tre1_1_0 + tre1_2_0;
	  c_im(out[2 * ostride]) = tim0_2_0 + tim1_1_0 - tim1_2_0;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre0_2_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_1_0 + tim1_2_0);
	       c_re(out[5 * ostride]) = tre2_0_0 + tre2_1_0;
	       c_re(out[8 * ostride]) = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim0_2_0 + (((FFTW_REAL) FFTW_K499999999) * (tim1_2_0 - tim1_1_0));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_2_0 - tre1_1_0);
	       c_im(out[5 * ostride]) = tim2_0_0 + tim2_1_0;
	       c_im(out[8 * ostride]) = tim2_0_0 - tim2_1_0;
	  }
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 0 FP additions and 0 FP multiplications */

void fftwi_no_twiddle_1(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     tre0_0_0 = c_re(in[0]);
     tim0_0_0 = c_im(in[0]);
     c_re(out[0]) = tre0_0_0;
     c_im(out[0]) = tim0_0_0;
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 108 FP additions and 32 FP multiplications */

void fftwi_no_twiddle_10(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_0_1;
     FFTW_REAL tim0_0_1;
     FFTW_REAL tre0_0_2;
     FFTW_REAL tim0_0_2;
     FFTW_REAL tre0_0_3;
     FFTW_REAL tim0_0_3;
     FFTW_REAL tre0_0_4;
     FFTW_REAL tim0_0_4;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_1_1;
     FFTW_REAL tim0_1_1;
     FFTW_REAL tre0_1_2;
     FFTW_REAL tim0_1_2;
     FFTW_REAL tre0_1_3;
     FFTW_REAL tim0_1_3;
     FFTW_REAL tre0_1_4;
     FFTW_REAL tim0_1_4;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[0]);
	  tim1_0_0 = c_im(in[0]);
	  tre1_1_0 = c_re(in[5 * istride]);
	  tim1_1_0 = c_im(in[5 * istride]);
	  tre0_0_0 = tre1_0_0 + tre1_1_0;
	  tim0_0_0 = tim1_0_0 + tim1_1_0;
	  tre0_1_0 = tre1_0_0 - tre1_1_0;
	  tim0_1_0 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[2 * istride]);
	  tim1_0_0 = c_im(in[2 * istride]);
	  tre1_1_0 = c_re(in[7 * istride]);
	  tim1_1_0 = c_im(in[7 * istride]);
	  tre0_0_1 = tre1_0_0 + tre1_1_0;
	  tim0_0_1 = tim1_0_0 + tim1_1_0;
	  tre0_1_1 = tre1_0_0 - tre1_1_0;
	  tim0_1_1 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[4 * istride]);
	  tim1_0_0 = c_im(in[4 * istride]);
	  tre1_1_0 = c_re(in[9 * istride]);
	  tim1_1_0 = c_im(in[9 * istride]);
	  tre0_0_2 = tre1_0_0 + tre1_1_0;
	  tim0_0_2 = tim1_0_0 + tim1_1_0;
	  tre0_1_2 = tre1_0_0 - tre1_1_0;
	  tim0_1_2 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[6 * istride]);
	  tim1_0_0 = c_im(in[6 * istride]);
	  tre1_1_0 = c_re(in[istride]);
	  tim1_1_0 = c_im(in[istride]);
	  tre0_0_3 = tre1_0_0 + tre1_1_0;
	  tim0_0_3 = tim1_0_0 + tim1_1_0;
	  tre0_1_3 = tre1_0_0 - tre1_1_0;
	  tim0_1_3 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[8 * istride]);
	  tim1_0_0 = c_im(in[8 * istride]);
	  tre1_1_0 = c_re(in[3 * istride]);
	  tim1_1_0 = c_im(in[3 * istride]);
	  tre0_0_4 = tre1_0_0 + tre1_1_0;
	  tim0_0_4 = tim1_0_0 + tim1_1_0;
	  tre0_1_4 = tre1_0_0 - tre1_1_0;
	  tim0_1_4 = tim1_0_0 - tim1_1_0;
     }
     c_re(out[0]) = tre0_0_0 + tre0_0_1 + tre0_0_2 + tre0_0_3 + tre0_0_4;
     c_im(out[0]) = tim0_0_0 + tim0_0_1 + tim0_0_2 + tim0_0_3 + tim0_0_4;
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_0_1 + tre0_0_4)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_0_2 + tre0_0_3));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tim0_0_4 - tim0_0_1)) + (((FFTW_REAL) FFTW_K587785252) * (tim0_0_3 - tim0_0_2));
	  c_re(out[6 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[4 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_0_1 + tim0_0_4)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_0_2 + tim0_0_3));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tre0_0_1 - tre0_0_4)) + (((FFTW_REAL) FFTW_K587785252) * (tre0_0_2 - tre0_0_3));
	  c_im(out[6 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[4 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_0_2 + tre0_0_3)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_0_1 + tre0_0_4));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tim0_0_4 - tim0_0_1)) + (((FFTW_REAL) FFTW_K951056516) * (tim0_0_2 - tim0_0_3));
	  c_re(out[2 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[8 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_0_2 + tim0_0_3)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_0_1 + tim0_0_4));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tre0_0_1 - tre0_0_4)) + (((FFTW_REAL) FFTW_K951056516) * (tre0_0_3 - tre0_0_2));
	  c_im(out[2 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[8 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     c_re(out[5 * ostride]) = tre0_1_0 + tre0_1_1 + tre0_1_2 + tre0_1_3 + tre0_1_4;
     c_im(out[5 * ostride]) = tim0_1_0 + tim0_1_1 + tim0_1_2 + tim0_1_3 + tim0_1_4;
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_1_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_1_1 + tre0_1_4)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_1_2 + tre0_1_3));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tim0_1_4 - tim0_1_1)) + (((FFTW_REAL) FFTW_K587785252) * (tim0_1_3 - tim0_1_2));
	  c_re(out[ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[9 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_1_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_1_1 + tim0_1_4)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_1_2 + tim0_1_3));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tre0_1_1 - tre0_1_4)) + (((FFTW_REAL) FFTW_K587785252) * (tre0_1_2 - tre0_1_3));
	  c_im(out[ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[9 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_1_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_1_2 + tre0_1_3)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_1_1 + tre0_1_4));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tim0_1_4 - tim0_1_1)) + (((FFTW_REAL) FFTW_K951056516) * (tim0_1_2 - tim0_1_3));
	  c_re(out[7 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[3 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_1_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_1_2 + tim0_1_3)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_1_1 + tim0_1_4));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tre0_1_1 - tre0_1_4)) + (((FFTW_REAL) FFTW_K951056516) * (tre0_1_3 - tre0_1_2));
	  c_im(out[7 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[3 * ostride]) = tim2_0_0 - tim2_1_0;
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 230 FP additions and 100 FP multiplications */

void fftwi_no_twiddle_11(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_2_0;
     FFTW_REAL tim0_2_0;
     FFTW_REAL tre0_3_0;
     FFTW_REAL tim0_3_0;
     FFTW_REAL tre0_4_0;
     FFTW_REAL tim0_4_0;
     FFTW_REAL tre0_5_0;
     FFTW_REAL tim0_5_0;
     FFTW_REAL tre0_6_0;
     FFTW_REAL tim0_6_0;
     FFTW_REAL tre0_7_0;
     FFTW_REAL tim0_7_0;
     FFTW_REAL tre0_8_0;
     FFTW_REAL tim0_8_0;
     FFTW_REAL tre0_9_0;
     FFTW_REAL tim0_9_0;
     FFTW_REAL tre0_10_0;
     FFTW_REAL tim0_10_0;
     tre0_0_0 = c_re(in[0]);
     tim0_0_0 = c_im(in[0]);
     tre0_1_0 = c_re(in[istride]);
     tim0_1_0 = c_im(in[istride]);
     tre0_2_0 = c_re(in[2 * istride]);
     tim0_2_0 = c_im(in[2 * istride]);
     tre0_3_0 = c_re(in[3 * istride]);
     tim0_3_0 = c_im(in[3 * istride]);
     tre0_4_0 = c_re(in[4 * istride]);
     tim0_4_0 = c_im(in[4 * istride]);
     tre0_5_0 = c_re(in[5 * istride]);
     tim0_5_0 = c_im(in[5 * istride]);
     tre0_6_0 = c_re(in[6 * istride]);
     tim0_6_0 = c_im(in[6 * istride]);
     tre0_7_0 = c_re(in[7 * istride]);
     tim0_7_0 = c_im(in[7 * istride]);
     tre0_8_0 = c_re(in[8 * istride]);
     tim0_8_0 = c_im(in[8 * istride]);
     tre0_9_0 = c_re(in[9 * istride]);
     tim0_9_0 = c_im(in[9 * istride]);
     tre0_10_0 = c_re(in[10 * istride]);
     tim0_10_0 = c_im(in[10 * istride]);
     c_re(out[0]) = tre0_0_0 + tre0_1_0 + tre0_2_0 + tre0_3_0 + tre0_4_0 + tre0_5_0 + tre0_6_0 + tre0_7_0 + tre0_8_0 + tre0_9_0 + tre0_10_0;
     c_im(out[0]) = tim0_0_0 + tim0_1_0 + tim0_2_0 + tim0_3_0 + tim0_4_0 + tim0_5_0 + tim0_6_0 + tim0_7_0 + tim0_8_0 + tim0_9_0 + tim0_10_0;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K841253532) * (tre0_1_0 + tre0_10_0)) + (((FFTW_REAL) FFTW_K415415013) * (tre0_2_0 + tre0_9_0)) - (((FFTW_REAL) FFTW_K959492973) * (tre0_5_0 + tre0_6_0)) - (((FFTW_REAL) FFTW_K654860733) * (tre0_4_0 + tre0_7_0)) - (((FFTW_REAL) FFTW_K142314838) * (tre0_3_0 + tre0_8_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K540640817) * (tim0_10_0 - tim0_1_0)) + (((FFTW_REAL) FFTW_K909631995) * (tim0_9_0 - tim0_2_0)) + (((FFTW_REAL) FFTW_K989821441) * (tim0_8_0 - tim0_3_0)) + (((FFTW_REAL) FFTW_K755749574) * (tim0_7_0 - tim0_4_0)) + (((FFTW_REAL) FFTW_K281732556) * (tim0_6_0 - tim0_5_0));
	  c_re(out[ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[10 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K841253532) * (tim0_1_0 + tim0_10_0)) + (((FFTW_REAL) FFTW_K415415013) * (tim0_2_0 + tim0_9_0)) - (((FFTW_REAL) FFTW_K959492973) * (tim0_5_0 + tim0_6_0)) - (((FFTW_REAL) FFTW_K654860733) * (tim0_4_0 + tim0_7_0)) - (((FFTW_REAL) FFTW_K142314838) * (tim0_3_0 + tim0_8_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K540640817) * (tre0_1_0 - tre0_10_0)) + (((FFTW_REAL) FFTW_K909631995) * (tre0_2_0 - tre0_9_0)) + (((FFTW_REAL) FFTW_K989821441) * (tre0_3_0 - tre0_8_0)) + (((FFTW_REAL) FFTW_K755749574) * (tre0_4_0 - tre0_7_0)) + (((FFTW_REAL) FFTW_K281732556) * (tre0_5_0 - tre0_6_0));
	  c_im(out[ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[10 * ostride]) = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K415415013) * (tre0_1_0 + tre0_10_0)) + (((FFTW_REAL) FFTW_K841253532) * (tre0_5_0 + tre0_6_0)) - (((FFTW_REAL) FFTW_K142314838) * (tre0_4_0 + tre0_7_0)) - (((FFTW_REAL) FFTW_K959492973) * (tre0_3_0 + tre0_8_0)) - (((FFTW_REAL) FFTW_K654860733) * (tre0_2_0 + tre0_9_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K909631995) * (tim0_10_0 - tim0_1_0)) + (((FFTW_REAL) FFTW_K755749574) * (tim0_9_0 - tim0_2_0)) + (((FFTW_REAL) FFTW_K281732556) * (tim0_3_0 - tim0_8_0)) + (((FFTW_REAL) FFTW_K989821441) * (tim0_4_0 - tim0_7_0)) + (((FFTW_REAL) FFTW_K540640817) * (tim0_5_0 - tim0_6_0));
	  c_re(out[2 * ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[9 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K415415013) * (tim0_1_0 + tim0_10_0)) + (((FFTW_REAL) FFTW_K841253532) * (tim0_5_0 + tim0_6_0)) - (((FFTW_REAL) FFTW_K142314838) * (tim0_4_0 + tim0_7_0)) - (((FFTW_REAL) FFTW_K959492973) * (tim0_3_0 + tim0_8_0)) - (((FFTW_REAL) FFTW_K654860733) * (tim0_2_0 + tim0_9_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K909631995) * (tre0_1_0 - tre0_10_0)) + (((FFTW_REAL) FFTW_K755749574) * (tre0_2_0 - tre0_9_0)) + (((FFTW_REAL) FFTW_K281732556) * (tre0_8_0 - tre0_3_0)) + (((FFTW_REAL) FFTW_K989821441) * (tre0_7_0 - tre0_4_0)) + (((FFTW_REAL) FFTW_K540640817) * (tre0_6_0 - tre0_5_0));
	  c_im(out[2 * ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[9 * ostride]) = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K415415013) * (tre0_3_0 + tre0_8_0)) + (((FFTW_REAL) FFTW_K841253532) * (tre0_4_0 + tre0_7_0)) - (((FFTW_REAL) FFTW_K654860733) * (tre0_5_0 + tre0_6_0)) - (((FFTW_REAL) FFTW_K959492973) * (tre0_2_0 + tre0_9_0)) - (((FFTW_REAL) FFTW_K142314838) * (tre0_1_0 + tre0_10_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K989821441) * (tim0_10_0 - tim0_1_0)) + (((FFTW_REAL) FFTW_K281732556) * (tim0_2_0 - tim0_9_0)) + (((FFTW_REAL) FFTW_K909631995) * (tim0_3_0 - tim0_8_0)) + (((FFTW_REAL) FFTW_K540640817) * (tim0_7_0 - tim0_4_0)) + (((FFTW_REAL) FFTW_K755749574) * (tim0_6_0 - tim0_5_0));
	  c_re(out[3 * ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[8 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K415415013) * (tim0_3_0 + tim0_8_0)) + (((FFTW_REAL) FFTW_K841253532) * (tim0_4_0 + tim0_7_0)) - (((FFTW_REAL) FFTW_K654860733) * (tim0_5_0 + tim0_6_0)) - (((FFTW_REAL) FFTW_K959492973) * (tim0_2_0 + tim0_9_0)) - (((FFTW_REAL) FFTW_K142314838) * (tim0_1_0 + tim0_10_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K989821441) * (tre0_1_0 - tre0_10_0)) + (((FFTW_REAL) FFTW_K281732556) * (tre0_9_0 - tre0_2_0)) + (((FFTW_REAL) FFTW_K909631995) * (tre0_8_0 - tre0_3_0)) + (((FFTW_REAL) FFTW_K540640817) * (tre0_4_0 - tre0_7_0)) + (((FFTW_REAL) FFTW_K755749574) * (tre0_5_0 - tre0_6_0));
	  c_im(out[3 * ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[8 * ostride]) = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K841253532) * (tre0_3_0 + tre0_8_0)) + (((FFTW_REAL) FFTW_K415415013) * (tre0_5_0 + tre0_6_0)) - (((FFTW_REAL) FFTW_K959492973) * (tre0_4_0 + tre0_7_0)) - (((FFTW_REAL) FFTW_K142314838) * (tre0_2_0 + tre0_9_0)) - (((FFTW_REAL) FFTW_K654860733) * (tre0_1_0 + tre0_10_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K755749574) * (tim0_10_0 - tim0_1_0)) + (((FFTW_REAL) FFTW_K989821441) * (tim0_2_0 - tim0_9_0)) + (((FFTW_REAL) FFTW_K540640817) * (tim0_8_0 - tim0_3_0)) + (((FFTW_REAL) FFTW_K281732556) * (tim0_7_0 - tim0_4_0)) + (((FFTW_REAL) FFTW_K909631995) * (tim0_5_0 - tim0_6_0));
	  c_re(out[4 * ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[7 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K841253532) * (tim0_3_0 + tim0_8_0)) + (((FFTW_REAL) FFTW_K415415013) * (tim0_5_0 + tim0_6_0)) - (((FFTW_REAL) FFTW_K959492973) * (tim0_4_0 + tim0_7_0)) - (((FFTW_REAL) FFTW_K142314838) * (tim0_2_0 + tim0_9_0)) - (((FFTW_REAL) FFTW_K654860733) * (tim0_1_0 + tim0_10_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K755749574) * (tre0_1_0 - tre0_10_0)) + (((FFTW_REAL) FFTW_K989821441) * (tre0_9_0 - tre0_2_0)) + (((FFTW_REAL) FFTW_K540640817) * (tre0_3_0 - tre0_8_0)) + (((FFTW_REAL) FFTW_K281732556) * (tre0_4_0 - tre0_7_0)) + (((FFTW_REAL) FFTW_K909631995) * (tre0_6_0 - tre0_5_0));
	  c_im(out[4 * ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[7 * ostride]) = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K841253532) * (tre0_2_0 + tre0_9_0)) + (((FFTW_REAL) FFTW_K415415013) * (tre0_4_0 + tre0_7_0)) - (((FFTW_REAL) FFTW_K142314838) * (tre0_5_0 + tre0_6_0)) - (((FFTW_REAL) FFTW_K654860733) * (tre0_3_0 + tre0_8_0)) - (((FFTW_REAL) FFTW_K959492973) * (tre0_1_0 + tre0_10_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K281732556) * (tim0_10_0 - tim0_1_0)) + (((FFTW_REAL) FFTW_K540640817) * (tim0_2_0 - tim0_9_0)) + (((FFTW_REAL) FFTW_K755749574) * (tim0_8_0 - tim0_3_0)) + (((FFTW_REAL) FFTW_K909631995) * (tim0_4_0 - tim0_7_0)) + (((FFTW_REAL) FFTW_K989821441) * (tim0_6_0 - tim0_5_0));
	  c_re(out[5 * ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[6 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K841253532) * (tim0_2_0 + tim0_9_0)) + (((FFTW_REAL) FFTW_K415415013) * (tim0_4_0 + tim0_7_0)) - (((FFTW_REAL) FFTW_K142314838) * (tim0_5_0 + tim0_6_0)) - (((FFTW_REAL) FFTW_K654860733) * (tim0_3_0 + tim0_8_0)) - (((FFTW_REAL) FFTW_K959492973) * (tim0_1_0 + tim0_10_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K281732556) * (tre0_1_0 - tre0_10_0)) + (((FFTW_REAL) FFTW_K540640817) * (tre0_9_0 - tre0_2_0)) + (((FFTW_REAL) FFTW_K755749574) * (tre0_3_0 - tre0_8_0)) + (((FFTW_REAL) FFTW_K909631995) * (tre0_7_0 - tre0_4_0)) + (((FFTW_REAL) FFTW_K989821441) * (tre0_5_0 - tre0_6_0));
	  c_im(out[5 * ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[6 * ostride]) = tim1_0_0 - tim1_1_0;
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 104 FP additions and 16 FP multiplications */

void fftwi_no_twiddle_12(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_0_1;
     FFTW_REAL tim0_0_1;
     FFTW_REAL tre0_0_2;
     FFTW_REAL tim0_0_2;
     FFTW_REAL tre0_0_3;
     FFTW_REAL tim0_0_3;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_1_1;
     FFTW_REAL tim0_1_1;
     FFTW_REAL tre0_1_2;
     FFTW_REAL tim0_1_2;
     FFTW_REAL tre0_1_3;
     FFTW_REAL tim0_1_3;
     FFTW_REAL tre0_2_0;
     FFTW_REAL tim0_2_0;
     FFTW_REAL tre0_2_1;
     FFTW_REAL tim0_2_1;
     FFTW_REAL tre0_2_2;
     FFTW_REAL tim0_2_2;
     FFTW_REAL tre0_2_3;
     FFTW_REAL tim0_2_3;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_2_0;
	  FFTW_REAL tim1_2_0;
	  tre1_0_0 = c_re(in[0]);
	  tim1_0_0 = c_im(in[0]);
	  tre1_1_0 = c_re(in[4 * istride]);
	  tim1_1_0 = c_im(in[4 * istride]);
	  tre1_2_0 = c_re(in[8 * istride]);
	  tim1_2_0 = c_im(in[8 * istride]);
	  tre0_0_0 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	  tim0_0_0 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_2_0 - tim1_1_0);
	       tre0_1_0 = tre2_0_0 + tre2_1_0;
	       tre0_2_0 = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_1_0 - tre1_2_0);
	       tim0_1_0 = tim2_0_0 + tim2_1_0;
	       tim0_2_0 = tim2_0_0 - tim2_1_0;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_2_0;
	  FFTW_REAL tim1_2_0;
	  tre1_0_0 = c_re(in[3 * istride]);
	  tim1_0_0 = c_im(in[3 * istride]);
	  tre1_1_0 = c_re(in[7 * istride]);
	  tim1_1_0 = c_im(in[7 * istride]);
	  tre1_2_0 = c_re(in[11 * istride]);
	  tim1_2_0 = c_im(in[11 * istride]);
	  tre0_0_1 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	  tim0_0_1 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_2_0 - tim1_1_0);
	       tre0_1_1 = tre2_0_0 + tre2_1_0;
	       tre0_2_1 = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_1_0 - tre1_2_0);
	       tim0_1_1 = tim2_0_0 + tim2_1_0;
	       tim0_2_1 = tim2_0_0 - tim2_1_0;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_2_0;
	  FFTW_REAL tim1_2_0;
	  tre1_0_0 = c_re(in[6 * istride]);
	  tim1_0_0 = c_im(in[6 * istride]);
	  tre1_1_0 = c_re(in[10 * istride]);
	  tim1_1_0 = c_im(in[10 * istride]);
	  tre1_2_0 = c_re(in[2 * istride]);
	  tim1_2_0 = c_im(in[2 * istride]);
	  tre0_0_2 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	  tim0_0_2 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_2_0 - tim1_1_0);
	       tre0_1_2 = tre2_0_0 + tre2_1_0;
	       tre0_2_2 = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_1_0 - tre1_2_0);
	       tim0_1_2 = tim2_0_0 + tim2_1_0;
	       tim0_2_2 = tim2_0_0 - tim2_1_0;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_2_0;
	  FFTW_REAL tim1_2_0;
	  tre1_0_0 = c_re(in[9 * istride]);
	  tim1_0_0 = c_im(in[9 * istride]);
	  tre1_1_0 = c_re(in[istride]);
	  tim1_1_0 = c_im(in[istride]);
	  tre1_2_0 = c_re(in[5 * istride]);
	  tim1_2_0 = c_im(in[5 * istride]);
	  tre0_0_3 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	  tim0_0_3 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_2_0 - tim1_1_0);
	       tre0_1_3 = tre2_0_0 + tre2_1_0;
	       tre0_2_3 = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_1_0 - tre1_2_0);
	       tim0_1_3 = tim2_0_0 + tim2_1_0;
	       tim0_2_3 = tim2_0_0 - tim2_1_0;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  tre1_0_0 = tre0_0_0 + tre0_0_2;
	  tim1_0_0 = tim0_0_0 + tim0_0_2;
	  tre1_1_0 = tre0_0_0 - tre0_0_2;
	  tim1_1_0 = tim0_0_0 - tim0_0_2;
	  tre1_0_1 = tre0_0_1 + tre0_0_3;
	  tim1_0_1 = tim0_0_1 + tim0_0_3;
	  tre1_1_1 = tre0_0_1 - tre0_0_3;
	  tim1_1_1 = tim0_0_1 - tim0_0_3;
	  c_re(out[0]) = tre1_0_0 + tre1_0_1;
	  c_im(out[0]) = tim1_0_0 + tim1_0_1;
	  c_re(out[6 * ostride]) = tre1_0_0 - tre1_0_1;
	  c_im(out[6 * ostride]) = tim1_0_0 - tim1_0_1;
	  c_re(out[9 * ostride]) = tre1_1_0 - tim1_1_1;
	  c_im(out[9 * ostride]) = tim1_1_0 + tre1_1_1;
	  c_re(out[3 * ostride]) = tre1_1_0 + tim1_1_1;
	  c_im(out[3 * ostride]) = tim1_1_0 - tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  tre1_0_0 = tre0_1_0 + tre0_1_2;
	  tim1_0_0 = tim0_1_0 + tim0_1_2;
	  tre1_1_0 = tre0_1_0 - tre0_1_2;
	  tim1_1_0 = tim0_1_0 - tim0_1_2;
	  tre1_0_1 = tre0_1_1 + tre0_1_3;
	  tim1_0_1 = tim0_1_1 + tim0_1_3;
	  tre1_1_1 = tre0_1_1 - tre0_1_3;
	  tim1_1_1 = tim0_1_1 - tim0_1_3;
	  c_re(out[4 * ostride]) = tre1_0_0 + tre1_0_1;
	  c_im(out[4 * ostride]) = tim1_0_0 + tim1_0_1;
	  c_re(out[10 * ostride]) = tre1_0_0 - tre1_0_1;
	  c_im(out[10 * ostride]) = tim1_0_0 - tim1_0_1;
	  c_re(out[ostride]) = tre1_1_0 - tim1_1_1;
	  c_im(out[ostride]) = tim1_1_0 + tre1_1_1;
	  c_re(out[7 * ostride]) = tre1_1_0 + tim1_1_1;
	  c_im(out[7 * ostride]) = tim1_1_0 - tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  tre1_0_0 = tre0_2_0 + tre0_2_2;
	  tim1_0_0 = tim0_2_0 + tim0_2_2;
	  tre1_1_0 = tre0_2_0 - tre0_2_2;
	  tim1_1_0 = tim0_2_0 - tim0_2_2;
	  tre1_0_1 = tre0_2_1 + tre0_2_3;
	  tim1_0_1 = tim0_2_1 + tim0_2_3;
	  tre1_1_1 = tre0_2_1 - tre0_2_3;
	  tim1_1_1 = tim0_2_1 - tim0_2_3;
	  c_re(out[8 * ostride]) = tre1_0_0 + tre1_0_1;
	  c_im(out[8 * ostride]) = tim1_0_0 + tim1_0_1;
	  c_re(out[2 * ostride]) = tre1_0_0 - tre1_0_1;
	  c_im(out[2 * ostride]) = tim1_0_0 - tim1_0_1;
	  c_re(out[5 * ostride]) = tre1_1_0 - tim1_1_1;
	  c_im(out[5 * ostride]) = tim1_1_0 + tre1_1_1;
	  c_re(out[11 * ostride]) = tre1_1_0 + tim1_1_1;
	  c_im(out[11 * ostride]) = tim1_1_0 - tre1_1_1;
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 324 FP additions and 144 FP multiplications */

void fftwi_no_twiddle_13(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_2_0;
     FFTW_REAL tim0_2_0;
     FFTW_REAL tre0_3_0;
     FFTW_REAL tim0_3_0;
     FFTW_REAL tre0_4_0;
     FFTW_REAL tim0_4_0;
     FFTW_REAL tre0_5_0;
     FFTW_REAL tim0_5_0;
     FFTW_REAL tre0_6_0;
     FFTW_REAL tim0_6_0;
     FFTW_REAL tre0_7_0;
     FFTW_REAL tim0_7_0;
     FFTW_REAL tre0_8_0;
     FFTW_REAL tim0_8_0;
     FFTW_REAL tre0_9_0;
     FFTW_REAL tim0_9_0;
     FFTW_REAL tre0_10_0;
     FFTW_REAL tim0_10_0;
     FFTW_REAL tre0_11_0;
     FFTW_REAL tim0_11_0;
     FFTW_REAL tre0_12_0;
     FFTW_REAL tim0_12_0;
     tre0_0_0 = c_re(in[0]);
     tim0_0_0 = c_im(in[0]);
     tre0_1_0 = c_re(in[istride]);
     tim0_1_0 = c_im(in[istride]);
     tre0_2_0 = c_re(in[2 * istride]);
     tim0_2_0 = c_im(in[2 * istride]);
     tre0_3_0 = c_re(in[3 * istride]);
     tim0_3_0 = c_im(in[3 * istride]);
     tre0_4_0 = c_re(in[4 * istride]);
     tim0_4_0 = c_im(in[4 * istride]);
     tre0_5_0 = c_re(in[5 * istride]);
     tim0_5_0 = c_im(in[5 * istride]);
     tre0_6_0 = c_re(in[6 * istride]);
     tim0_6_0 = c_im(in[6 * istride]);
     tre0_7_0 = c_re(in[7 * istride]);
     tim0_7_0 = c_im(in[7 * istride]);
     tre0_8_0 = c_re(in[8 * istride]);
     tim0_8_0 = c_im(in[8 * istride]);
     tre0_9_0 = c_re(in[9 * istride]);
     tim0_9_0 = c_im(in[9 * istride]);
     tre0_10_0 = c_re(in[10 * istride]);
     tim0_10_0 = c_im(in[10 * istride]);
     tre0_11_0 = c_re(in[11 * istride]);
     tim0_11_0 = c_im(in[11 * istride]);
     tre0_12_0 = c_re(in[12 * istride]);
     tim0_12_0 = c_im(in[12 * istride]);
     c_re(out[0]) = tre0_0_0 + tre0_1_0 + tre0_2_0 + tre0_3_0 + tre0_4_0 + tre0_5_0 + tre0_6_0 + tre0_7_0 + tre0_8_0 + tre0_9_0 + tre0_10_0 + tre0_11_0 + tre0_12_0;
     c_im(out[0]) = tim0_0_0 + tim0_1_0 + tim0_2_0 + tim0_3_0 + tim0_4_0 + tim0_5_0 + tim0_6_0 + tim0_7_0 + tim0_8_0 + tim0_9_0 + tim0_10_0 + tim0_11_0 + tim0_12_0;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K885456025) * (tre0_1_0 + tre0_12_0)) + (((FFTW_REAL) FFTW_K568064746) * (tre0_2_0 + tre0_11_0)) + (((FFTW_REAL) FFTW_K120536680) * (tre0_3_0 + tre0_10_0)) - (((FFTW_REAL) FFTW_K970941817) * (tre0_6_0 + tre0_7_0)) - (((FFTW_REAL) FFTW_K748510748) * (tre0_5_0 + tre0_8_0)) - (((FFTW_REAL) FFTW_K354604887) * (tre0_4_0 + tre0_9_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K464723172) * (tim0_12_0 - tim0_1_0)) + (((FFTW_REAL) FFTW_K822983865) * (tim0_11_0 - tim0_2_0)) + (((FFTW_REAL) FFTW_K992708874) * (tim0_10_0 - tim0_3_0)) + (((FFTW_REAL) FFTW_K935016242) * (tim0_9_0 - tim0_4_0)) + (((FFTW_REAL) FFTW_K663122658) * (tim0_8_0 - tim0_5_0)) + (((FFTW_REAL) FFTW_K239315664) * (tim0_7_0 - tim0_6_0));
	  c_re(out[ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[12 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K885456025) * (tim0_1_0 + tim0_12_0)) + (((FFTW_REAL) FFTW_K568064746) * (tim0_2_0 + tim0_11_0)) + (((FFTW_REAL) FFTW_K120536680) * (tim0_3_0 + tim0_10_0)) - (((FFTW_REAL) FFTW_K970941817) * (tim0_6_0 + tim0_7_0)) - (((FFTW_REAL) FFTW_K748510748) * (tim0_5_0 + tim0_8_0)) - (((FFTW_REAL) FFTW_K354604887) * (tim0_4_0 + tim0_9_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K464723172) * (tre0_1_0 - tre0_12_0)) + (((FFTW_REAL) FFTW_K822983865) * (tre0_2_0 - tre0_11_0)) + (((FFTW_REAL) FFTW_K992708874) * (tre0_3_0 - tre0_10_0)) + (((FFTW_REAL) FFTW_K935016242) * (tre0_4_0 - tre0_9_0)) + (((FFTW_REAL) FFTW_K663122658) * (tre0_5_0 - tre0_8_0)) + (((FFTW_REAL) FFTW_K239315664) * (tre0_6_0 - tre0_7_0));
	  c_im(out[ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[12 * ostride]) = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K568064746) * (tre0_1_0 + tre0_12_0)) + (((FFTW_REAL) FFTW_K120536680) * (tre0_5_0 + tre0_8_0)) + (((FFTW_REAL) FFTW_K885456025) * (tre0_6_0 + tre0_7_0)) - (((FFTW_REAL) FFTW_K748510748) * (tre0_4_0 + tre0_9_0)) - (((FFTW_REAL) FFTW_K970941817) * (tre0_3_0 + tre0_10_0)) - (((FFTW_REAL) FFTW_K354604887) * (tre0_2_0 + tre0_11_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K822983865) * (tim0_12_0 - tim0_1_0)) + (((FFTW_REAL) FFTW_K935016242) * (tim0_11_0 - tim0_2_0)) + (((FFTW_REAL) FFTW_K239315664) * (tim0_10_0 - tim0_3_0)) + (((FFTW_REAL) FFTW_K663122658) * (tim0_4_0 - tim0_9_0)) + (((FFTW_REAL) FFTW_K992708874) * (tim0_5_0 - tim0_8_0)) + (((FFTW_REAL) FFTW_K464723172) * (tim0_6_0 - tim0_7_0));
	  c_re(out[2 * ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[11 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K568064746) * (tim0_1_0 + tim0_12_0)) + (((FFTW_REAL) FFTW_K120536680) * (tim0_5_0 + tim0_8_0)) + (((FFTW_REAL) FFTW_K885456025) * (tim0_6_0 + tim0_7_0)) - (((FFTW_REAL) FFTW_K748510748) * (tim0_4_0 + tim0_9_0)) - (((FFTW_REAL) FFTW_K970941817) * (tim0_3_0 + tim0_10_0)) - (((FFTW_REAL) FFTW_K354604887) * (tim0_2_0 + tim0_11_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K822983865) * (tre0_1_0 - tre0_12_0)) + (((FFTW_REAL) FFTW_K935016242) * (tre0_2_0 - tre0_11_0)) + (((FFTW_REAL) FFTW_K239315664) * (tre0_3_0 - tre0_10_0)) + (((FFTW_REAL) FFTW_K663122658) * (tre0_9_0 - tre0_4_0)) + (((FFTW_REAL) FFTW_K992708874) * (tre0_8_0 - tre0_5_0)) + (((FFTW_REAL) FFTW_K464723172) * (tre0_7_0 - tre0_6_0));
	  c_im(out[2 * ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[11 * ostride]) = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K120536680) * (tre0_1_0 + tre0_12_0)) + (((FFTW_REAL) FFTW_K885456025) * (tre0_4_0 + tre0_9_0)) + (((FFTW_REAL) FFTW_K568064746) * (tre0_5_0 + tre0_8_0)) - (((FFTW_REAL) FFTW_K748510748) * (tre0_6_0 + tre0_7_0)) - (((FFTW_REAL) FFTW_K354604887) * (tre0_3_0 + tre0_10_0)) - (((FFTW_REAL) FFTW_K970941817) * (tre0_2_0 + tre0_11_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K992708874) * (tim0_12_0 - tim0_1_0)) + (((FFTW_REAL) FFTW_K239315664) * (tim0_11_0 - tim0_2_0)) + (((FFTW_REAL) FFTW_K935016242) * (tim0_3_0 - tim0_10_0)) + (((FFTW_REAL) FFTW_K464723172) * (tim0_4_0 - tim0_9_0)) + (((FFTW_REAL) FFTW_K822983865) * (tim0_8_0 - tim0_5_0)) + (((FFTW_REAL) FFTW_K663122658) * (tim0_7_0 - tim0_6_0));
	  c_re(out[3 * ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[10 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K120536680) * (tim0_1_0 + tim0_12_0)) + (((FFTW_REAL) FFTW_K885456025) * (tim0_4_0 + tim0_9_0)) + (((FFTW_REAL) FFTW_K568064746) * (tim0_5_0 + tim0_8_0)) - (((FFTW_REAL) FFTW_K748510748) * (tim0_6_0 + tim0_7_0)) - (((FFTW_REAL) FFTW_K354604887) * (tim0_3_0 + tim0_10_0)) - (((FFTW_REAL) FFTW_K970941817) * (tim0_2_0 + tim0_11_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K992708874) * (tre0_1_0 - tre0_12_0)) + (((FFTW_REAL) FFTW_K239315664) * (tre0_2_0 - tre0_11_0)) + (((FFTW_REAL) FFTW_K935016242) * (tre0_10_0 - tre0_3_0)) + (((FFTW_REAL) FFTW_K464723172) * (tre0_9_0 - tre0_4_0)) + (((FFTW_REAL) FFTW_K822983865) * (tre0_5_0 - tre0_8_0)) + (((FFTW_REAL) FFTW_K663122658) * (tre0_6_0 - tre0_7_0));
	  c_im(out[3 * ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[10 * ostride]) = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K885456025) * (tre0_3_0 + tre0_10_0)) + (((FFTW_REAL) FFTW_K120536680) * (tre0_4_0 + tre0_9_0)) + (((FFTW_REAL) FFTW_K568064746) * (tre0_6_0 + tre0_7_0)) - (((FFTW_REAL) FFTW_K970941817) * (tre0_5_0 + tre0_8_0)) - (((FFTW_REAL) FFTW_K748510748) * (tre0_2_0 + tre0_11_0)) - (((FFTW_REAL) FFTW_K354604887) * (tre0_1_0 + tre0_12_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K935016242) * (tim0_12_0 - tim0_1_0)) + (((FFTW_REAL) FFTW_K663122658) * (tim0_2_0 - tim0_11_0)) + (((FFTW_REAL) FFTW_K464723172) * (tim0_3_0 - tim0_10_0)) + (((FFTW_REAL) FFTW_K992708874) * (tim0_9_0 - tim0_4_0)) + (((FFTW_REAL) FFTW_K239315664) * (tim0_5_0 - tim0_8_0)) + (((FFTW_REAL) FFTW_K822983865) * (tim0_6_0 - tim0_7_0));
	  c_re(out[4 * ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[9 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K885456025) * (tim0_3_0 + tim0_10_0)) + (((FFTW_REAL) FFTW_K120536680) * (tim0_4_0 + tim0_9_0)) + (((FFTW_REAL) FFTW_K568064746) * (tim0_6_0 + tim0_7_0)) - (((FFTW_REAL) FFTW_K970941817) * (tim0_5_0 + tim0_8_0)) - (((FFTW_REAL) FFTW_K748510748) * (tim0_2_0 + tim0_11_0)) - (((FFTW_REAL) FFTW_K354604887) * (tim0_1_0 + tim0_12_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K935016242) * (tre0_1_0 - tre0_12_0)) + (((FFTW_REAL) FFTW_K663122658) * (tre0_11_0 - tre0_2_0)) + (((FFTW_REAL) FFTW_K464723172) * (tre0_10_0 - tre0_3_0)) + (((FFTW_REAL) FFTW_K992708874) * (tre0_4_0 - tre0_9_0)) + (((FFTW_REAL) FFTW_K239315664) * (tre0_8_0 - tre0_5_0)) + (((FFTW_REAL) FFTW_K822983865) * (tre0_7_0 - tre0_6_0));
	  c_im(out[4 * ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[9 * ostride]) = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K120536680) * (tre0_2_0 + tre0_11_0)) + (((FFTW_REAL) FFTW_K568064746) * (tre0_3_0 + tre0_10_0)) + (((FFTW_REAL) FFTW_K885456025) * (tre0_5_0 + tre0_8_0)) - (((FFTW_REAL) FFTW_K354604887) * (tre0_6_0 + tre0_7_0)) - (((FFTW_REAL) FFTW_K970941817) * (tre0_4_0 + tre0_9_0)) - (((FFTW_REAL) FFTW_K748510748) * (tre0_1_0 + tre0_12_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K663122658) * (tim0_12_0 - tim0_1_0)) + (((FFTW_REAL) FFTW_K992708874) * (tim0_2_0 - tim0_11_0)) + (((FFTW_REAL) FFTW_K822983865) * (tim0_10_0 - tim0_3_0)) + (((FFTW_REAL) FFTW_K239315664) * (tim0_4_0 - tim0_9_0)) + (((FFTW_REAL) FFTW_K464723172) * (tim0_5_0 - tim0_8_0)) + (((FFTW_REAL) FFTW_K935016242) * (tim0_7_0 - tim0_6_0));
	  c_re(out[5 * ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[8 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K120536680) * (tim0_2_0 + tim0_11_0)) + (((FFTW_REAL) FFTW_K568064746) * (tim0_3_0 + tim0_10_0)) + (((FFTW_REAL) FFTW_K885456025) * (tim0_5_0 + tim0_8_0)) - (((FFTW_REAL) FFTW_K354604887) * (tim0_6_0 + tim0_7_0)) - (((FFTW_REAL) FFTW_K970941817) * (tim0_4_0 + tim0_9_0)) - (((FFTW_REAL) FFTW_K748510748) * (tim0_1_0 + tim0_12_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K663122658) * (tre0_1_0 - tre0_12_0)) + (((FFTW_REAL) FFTW_K992708874) * (tre0_11_0 - tre0_2_0)) + (((FFTW_REAL) FFTW_K822983865) * (tre0_3_0 - tre0_10_0)) + (((FFTW_REAL) FFTW_K239315664) * (tre0_9_0 - tre0_4_0)) + (((FFTW_REAL) FFTW_K464723172) * (tre0_8_0 - tre0_5_0)) + (((FFTW_REAL) FFTW_K935016242) * (tre0_6_0 - tre0_7_0));
	  c_im(out[5 * ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[8 * ostride]) = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K885456025) * (tre0_2_0 + tre0_11_0)) + (((FFTW_REAL) FFTW_K568064746) * (tre0_4_0 + tre0_9_0)) + (((FFTW_REAL) FFTW_K120536680) * (tre0_6_0 + tre0_7_0)) - (((FFTW_REAL) FFTW_K354604887) * (tre0_5_0 + tre0_8_0)) - (((FFTW_REAL) FFTW_K748510748) * (tre0_3_0 + tre0_10_0)) - (((FFTW_REAL) FFTW_K970941817) * (tre0_1_0 + tre0_12_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K239315664) * (tim0_12_0 - tim0_1_0)) + (((FFTW_REAL) FFTW_K464723172) * (tim0_2_0 - tim0_11_0)) + (((FFTW_REAL) FFTW_K663122658) * (tim0_10_0 - tim0_3_0)) + (((FFTW_REAL) FFTW_K822983865) * (tim0_4_0 - tim0_9_0)) + (((FFTW_REAL) FFTW_K935016242) * (tim0_8_0 - tim0_5_0)) + (((FFTW_REAL) FFTW_K992708874) * (tim0_6_0 - tim0_7_0));
	  c_re(out[6 * ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[7 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K885456025) * (tim0_2_0 + tim0_11_0)) + (((FFTW_REAL) FFTW_K568064746) * (tim0_4_0 + tim0_9_0)) + (((FFTW_REAL) FFTW_K120536680) * (tim0_6_0 + tim0_7_0)) - (((FFTW_REAL) FFTW_K354604887) * (tim0_5_0 + tim0_8_0)) - (((FFTW_REAL) FFTW_K748510748) * (tim0_3_0 + tim0_10_0)) - (((FFTW_REAL) FFTW_K970941817) * (tim0_1_0 + tim0_12_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K239315664) * (tre0_1_0 - tre0_12_0)) + (((FFTW_REAL) FFTW_K464723172) * (tre0_11_0 - tre0_2_0)) + (((FFTW_REAL) FFTW_K663122658) * (tre0_3_0 - tre0_10_0)) + (((FFTW_REAL) FFTW_K822983865) * (tre0_9_0 - tre0_4_0)) + (((FFTW_REAL) FFTW_K935016242) * (tre0_5_0 - tre0_8_0)) + (((FFTW_REAL) FFTW_K992708874) * (tre0_7_0 - tre0_6_0));
	  c_im(out[6 * ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[7 * ostride]) = tim1_0_0 - tim1_1_0;
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 208 FP additions and 72 FP multiplications */

void fftwi_no_twiddle_14(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_0_1;
     FFTW_REAL tim0_0_1;
     FFTW_REAL tre0_0_2;
     FFTW_REAL tim0_0_2;
     FFTW_REAL tre0_0_3;
     FFTW_REAL tim0_0_3;
     FFTW_REAL tre0_0_4;
     FFTW_REAL tim0_0_4;
     FFTW_REAL tre0_0_5;
     FFTW_REAL tim0_0_5;
     FFTW_REAL tre0_0_6;
     FFTW_REAL tim0_0_6;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_1_1;
     FFTW_REAL tim0_1_1;
     FFTW_REAL tre0_1_2;
     FFTW_REAL tim0_1_2;
     FFTW_REAL tre0_1_3;
     FFTW_REAL tim0_1_3;
     FFTW_REAL tre0_1_4;
     FFTW_REAL tim0_1_4;
     FFTW_REAL tre0_1_5;
     FFTW_REAL tim0_1_5;
     FFTW_REAL tre0_1_6;
     FFTW_REAL tim0_1_6;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[0]);
	  tim1_0_0 = c_im(in[0]);
	  tre1_1_0 = c_re(in[7 * istride]);
	  tim1_1_0 = c_im(in[7 * istride]);
	  tre0_0_0 = tre1_0_0 + tre1_1_0;
	  tim0_0_0 = tim1_0_0 + tim1_1_0;
	  tre0_1_0 = tre1_0_0 - tre1_1_0;
	  tim0_1_0 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[2 * istride]);
	  tim1_0_0 = c_im(in[2 * istride]);
	  tre1_1_0 = c_re(in[9 * istride]);
	  tim1_1_0 = c_im(in[9 * istride]);
	  tre0_0_1 = tre1_0_0 + tre1_1_0;
	  tim0_0_1 = tim1_0_0 + tim1_1_0;
	  tre0_1_1 = tre1_0_0 - tre1_1_0;
	  tim0_1_1 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[4 * istride]);
	  tim1_0_0 = c_im(in[4 * istride]);
	  tre1_1_0 = c_re(in[11 * istride]);
	  tim1_1_0 = c_im(in[11 * istride]);
	  tre0_0_2 = tre1_0_0 + tre1_1_0;
	  tim0_0_2 = tim1_0_0 + tim1_1_0;
	  tre0_1_2 = tre1_0_0 - tre1_1_0;
	  tim0_1_2 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[6 * istride]);
	  tim1_0_0 = c_im(in[6 * istride]);
	  tre1_1_0 = c_re(in[13 * istride]);
	  tim1_1_0 = c_im(in[13 * istride]);
	  tre0_0_3 = tre1_0_0 + tre1_1_0;
	  tim0_0_3 = tim1_0_0 + tim1_1_0;
	  tre0_1_3 = tre1_0_0 - tre1_1_0;
	  tim0_1_3 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[8 * istride]);
	  tim1_0_0 = c_im(in[8 * istride]);
	  tre1_1_0 = c_re(in[istride]);
	  tim1_1_0 = c_im(in[istride]);
	  tre0_0_4 = tre1_0_0 + tre1_1_0;
	  tim0_0_4 = tim1_0_0 + tim1_1_0;
	  tre0_1_4 = tre1_0_0 - tre1_1_0;
	  tim0_1_4 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[10 * istride]);
	  tim1_0_0 = c_im(in[10 * istride]);
	  tre1_1_0 = c_re(in[3 * istride]);
	  tim1_1_0 = c_im(in[3 * istride]);
	  tre0_0_5 = tre1_0_0 + tre1_1_0;
	  tim0_0_5 = tim1_0_0 + tim1_1_0;
	  tre0_1_5 = tre1_0_0 - tre1_1_0;
	  tim0_1_5 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[12 * istride]);
	  tim1_0_0 = c_im(in[12 * istride]);
	  tre1_1_0 = c_re(in[5 * istride]);
	  tim1_1_0 = c_im(in[5 * istride]);
	  tre0_0_6 = tre1_0_0 + tre1_1_0;
	  tim0_0_6 = tim1_0_0 + tim1_1_0;
	  tre0_1_6 = tre1_0_0 - tre1_1_0;
	  tim0_1_6 = tim1_0_0 - tim1_1_0;
     }
     c_re(out[0]) = tre0_0_0 + tre0_0_1 + tre0_0_2 + tre0_0_3 + tre0_0_4 + tre0_0_5 + tre0_0_6;
     c_im(out[0]) = tim0_0_0 + tim0_0_1 + tim0_0_2 + tim0_0_3 + tim0_0_4 + tim0_0_5 + tim0_0_6;
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tre0_0_1 + tre0_0_6)) - (((FFTW_REAL) FFTW_K900968867) * (tre0_0_3 + tre0_0_4)) - (((FFTW_REAL) FFTW_K222520933) * (tre0_0_2 + tre0_0_5));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K781831482) * (tim0_0_6 - tim0_0_1)) + (((FFTW_REAL) FFTW_K974927912) * (tim0_0_5 - tim0_0_2)) + (((FFTW_REAL) FFTW_K433883739) * (tim0_0_4 - tim0_0_3));
	  c_re(out[8 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[6 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tim0_0_1 + tim0_0_6)) - (((FFTW_REAL) FFTW_K900968867) * (tim0_0_3 + tim0_0_4)) - (((FFTW_REAL) FFTW_K222520933) * (tim0_0_2 + tim0_0_5));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K781831482) * (tre0_0_1 - tre0_0_6)) + (((FFTW_REAL) FFTW_K974927912) * (tre0_0_2 - tre0_0_5)) + (((FFTW_REAL) FFTW_K433883739) * (tre0_0_3 - tre0_0_4));
	  c_im(out[8 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[6 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tre0_0_3 + tre0_0_4)) - (((FFTW_REAL) FFTW_K900968867) * (tre0_0_2 + tre0_0_5)) - (((FFTW_REAL) FFTW_K222520933) * (tre0_0_1 + tre0_0_6));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K974927912) * (tim0_0_6 - tim0_0_1)) + (((FFTW_REAL) FFTW_K433883739) * (tim0_0_2 - tim0_0_5)) + (((FFTW_REAL) FFTW_K781831482) * (tim0_0_3 - tim0_0_4));
	  c_re(out[2 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[12 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tim0_0_3 + tim0_0_4)) - (((FFTW_REAL) FFTW_K900968867) * (tim0_0_2 + tim0_0_5)) - (((FFTW_REAL) FFTW_K222520933) * (tim0_0_1 + tim0_0_6));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K974927912) * (tre0_0_1 - tre0_0_6)) + (((FFTW_REAL) FFTW_K433883739) * (tre0_0_5 - tre0_0_2)) + (((FFTW_REAL) FFTW_K781831482) * (tre0_0_4 - tre0_0_3));
	  c_im(out[2 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[12 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tre0_0_2 + tre0_0_5)) - (((FFTW_REAL) FFTW_K222520933) * (tre0_0_3 + tre0_0_4)) - (((FFTW_REAL) FFTW_K900968867) * (tre0_0_1 + tre0_0_6));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K433883739) * (tim0_0_6 - tim0_0_1)) + (((FFTW_REAL) FFTW_K781831482) * (tim0_0_2 - tim0_0_5)) + (((FFTW_REAL) FFTW_K974927912) * (tim0_0_4 - tim0_0_3));
	  c_re(out[10 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[4 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tim0_0_2 + tim0_0_5)) - (((FFTW_REAL) FFTW_K222520933) * (tim0_0_3 + tim0_0_4)) - (((FFTW_REAL) FFTW_K900968867) * (tim0_0_1 + tim0_0_6));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K433883739) * (tre0_0_1 - tre0_0_6)) + (((FFTW_REAL) FFTW_K781831482) * (tre0_0_5 - tre0_0_2)) + (((FFTW_REAL) FFTW_K974927912) * (tre0_0_3 - tre0_0_4));
	  c_im(out[10 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[4 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     c_re(out[7 * ostride]) = tre0_1_0 + tre0_1_1 + tre0_1_2 + tre0_1_3 + tre0_1_4 + tre0_1_5 + tre0_1_6;
     c_im(out[7 * ostride]) = tim0_1_0 + tim0_1_1 + tim0_1_2 + tim0_1_3 + tim0_1_4 + tim0_1_5 + tim0_1_6;
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_1_0 + (((FFTW_REAL) FFTW_K623489801) * (tre0_1_1 + tre0_1_6)) - (((FFTW_REAL) FFTW_K900968867) * (tre0_1_3 + tre0_1_4)) - (((FFTW_REAL) FFTW_K222520933) * (tre0_1_2 + tre0_1_5));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K781831482) * (tim0_1_6 - tim0_1_1)) + (((FFTW_REAL) FFTW_K974927912) * (tim0_1_5 - tim0_1_2)) + (((FFTW_REAL) FFTW_K433883739) * (tim0_1_4 - tim0_1_3));
	  c_re(out[ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[13 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_1_0 + (((FFTW_REAL) FFTW_K623489801) * (tim0_1_1 + tim0_1_6)) - (((FFTW_REAL) FFTW_K900968867) * (tim0_1_3 + tim0_1_4)) - (((FFTW_REAL) FFTW_K222520933) * (tim0_1_2 + tim0_1_5));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K781831482) * (tre0_1_1 - tre0_1_6)) + (((FFTW_REAL) FFTW_K974927912) * (tre0_1_2 - tre0_1_5)) + (((FFTW_REAL) FFTW_K433883739) * (tre0_1_3 - tre0_1_4));
	  c_im(out[ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[13 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_1_0 + (((FFTW_REAL) FFTW_K623489801) * (tre0_1_3 + tre0_1_4)) - (((FFTW_REAL) FFTW_K900968867) * (tre0_1_2 + tre0_1_5)) - (((FFTW_REAL) FFTW_K222520933) * (tre0_1_1 + tre0_1_6));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K974927912) * (tim0_1_6 - tim0_1_1)) + (((FFTW_REAL) FFTW_K433883739) * (tim0_1_2 - tim0_1_5)) + (((FFTW_REAL) FFTW_K781831482) * (tim0_1_3 - tim0_1_4));
	  c_re(out[9 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[5 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_1_0 + (((FFTW_REAL) FFTW_K623489801) * (tim0_1_3 + tim0_1_4)) - (((FFTW_REAL) FFTW_K900968867) * (tim0_1_2 + tim0_1_5)) - (((FFTW_REAL) FFTW_K222520933) * (tim0_1_1 + tim0_1_6));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K974927912) * (tre0_1_1 - tre0_1_6)) + (((FFTW_REAL) FFTW_K433883739) * (tre0_1_5 - tre0_1_2)) + (((FFTW_REAL) FFTW_K781831482) * (tre0_1_4 - tre0_1_3));
	  c_im(out[9 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[5 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_1_0 + (((FFTW_REAL) FFTW_K623489801) * (tre0_1_2 + tre0_1_5)) - (((FFTW_REAL) FFTW_K222520933) * (tre0_1_3 + tre0_1_4)) - (((FFTW_REAL) FFTW_K900968867) * (tre0_1_1 + tre0_1_6));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K433883739) * (tim0_1_6 - tim0_1_1)) + (((FFTW_REAL) FFTW_K781831482) * (tim0_1_2 - tim0_1_5)) + (((FFTW_REAL) FFTW_K974927912) * (tim0_1_4 - tim0_1_3));
	  c_re(out[3 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[11 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_1_0 + (((FFTW_REAL) FFTW_K623489801) * (tim0_1_2 + tim0_1_5)) - (((FFTW_REAL) FFTW_K222520933) * (tim0_1_3 + tim0_1_4)) - (((FFTW_REAL) FFTW_K900968867) * (tim0_1_1 + tim0_1_6));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K433883739) * (tre0_1_1 - tre0_1_6)) + (((FFTW_REAL) FFTW_K781831482) * (tre0_1_5 - tre0_1_2)) + (((FFTW_REAL) FFTW_K974927912) * (tre0_1_3 - tre0_1_4));
	  c_im(out[3 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[11 * ostride]) = tim2_0_0 - tim2_1_0;
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 202 FP additions and 68 FP multiplications */

void fftwi_no_twiddle_15(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_0_1;
     FFTW_REAL tim0_0_1;
     FFTW_REAL tre0_0_2;
     FFTW_REAL tim0_0_2;
     FFTW_REAL tre0_0_3;
     FFTW_REAL tim0_0_3;
     FFTW_REAL tre0_0_4;
     FFTW_REAL tim0_0_4;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_1_1;
     FFTW_REAL tim0_1_1;
     FFTW_REAL tre0_1_2;
     FFTW_REAL tim0_1_2;
     FFTW_REAL tre0_1_3;
     FFTW_REAL tim0_1_3;
     FFTW_REAL tre0_1_4;
     FFTW_REAL tim0_1_4;
     FFTW_REAL tre0_2_0;
     FFTW_REAL tim0_2_0;
     FFTW_REAL tre0_2_1;
     FFTW_REAL tim0_2_1;
     FFTW_REAL tre0_2_2;
     FFTW_REAL tim0_2_2;
     FFTW_REAL tre0_2_3;
     FFTW_REAL tim0_2_3;
     FFTW_REAL tre0_2_4;
     FFTW_REAL tim0_2_4;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_2_0;
	  FFTW_REAL tim1_2_0;
	  tre1_0_0 = c_re(in[0]);
	  tim1_0_0 = c_im(in[0]);
	  tre1_1_0 = c_re(in[5 * istride]);
	  tim1_1_0 = c_im(in[5 * istride]);
	  tre1_2_0 = c_re(in[10 * istride]);
	  tim1_2_0 = c_im(in[10 * istride]);
	  tre0_0_0 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	  tim0_0_0 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_2_0 - tim1_1_0);
	       tre0_1_0 = tre2_0_0 + tre2_1_0;
	       tre0_2_0 = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_1_0 - tre1_2_0);
	       tim0_1_0 = tim2_0_0 + tim2_1_0;
	       tim0_2_0 = tim2_0_0 - tim2_1_0;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_2_0;
	  FFTW_REAL tim1_2_0;
	  tre1_0_0 = c_re(in[3 * istride]);
	  tim1_0_0 = c_im(in[3 * istride]);
	  tre1_1_0 = c_re(in[8 * istride]);
	  tim1_1_0 = c_im(in[8 * istride]);
	  tre1_2_0 = c_re(in[13 * istride]);
	  tim1_2_0 = c_im(in[13 * istride]);
	  tre0_0_1 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	  tim0_0_1 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_2_0 - tim1_1_0);
	       tre0_1_1 = tre2_0_0 + tre2_1_0;
	       tre0_2_1 = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_1_0 - tre1_2_0);
	       tim0_1_1 = tim2_0_0 + tim2_1_0;
	       tim0_2_1 = tim2_0_0 - tim2_1_0;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_2_0;
	  FFTW_REAL tim1_2_0;
	  tre1_0_0 = c_re(in[6 * istride]);
	  tim1_0_0 = c_im(in[6 * istride]);
	  tre1_1_0 = c_re(in[11 * istride]);
	  tim1_1_0 = c_im(in[11 * istride]);
	  tre1_2_0 = c_re(in[istride]);
	  tim1_2_0 = c_im(in[istride]);
	  tre0_0_2 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	  tim0_0_2 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_2_0 - tim1_1_0);
	       tre0_1_2 = tre2_0_0 + tre2_1_0;
	       tre0_2_2 = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_1_0 - tre1_2_0);
	       tim0_1_2 = tim2_0_0 + tim2_1_0;
	       tim0_2_2 = tim2_0_0 - tim2_1_0;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_2_0;
	  FFTW_REAL tim1_2_0;
	  tre1_0_0 = c_re(in[9 * istride]);
	  tim1_0_0 = c_im(in[9 * istride]);
	  tre1_1_0 = c_re(in[14 * istride]);
	  tim1_1_0 = c_im(in[14 * istride]);
	  tre1_2_0 = c_re(in[4 * istride]);
	  tim1_2_0 = c_im(in[4 * istride]);
	  tre0_0_3 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	  tim0_0_3 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_2_0 - tim1_1_0);
	       tre0_1_3 = tre2_0_0 + tre2_1_0;
	       tre0_2_3 = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_1_0 - tre1_2_0);
	       tim0_1_3 = tim2_0_0 + tim2_1_0;
	       tim0_2_3 = tim2_0_0 - tim2_1_0;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_2_0;
	  FFTW_REAL tim1_2_0;
	  tre1_0_0 = c_re(in[12 * istride]);
	  tim1_0_0 = c_im(in[12 * istride]);
	  tre1_1_0 = c_re(in[2 * istride]);
	  tim1_1_0 = c_im(in[2 * istride]);
	  tre1_2_0 = c_re(in[7 * istride]);
	  tim1_2_0 = c_im(in[7 * istride]);
	  tre0_0_4 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	  tim0_0_4 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_2_0 - tim1_1_0);
	       tre0_1_4 = tre2_0_0 + tre2_1_0;
	       tre0_2_4 = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_1_0 - tre1_2_0);
	       tim0_1_4 = tim2_0_0 + tim2_1_0;
	       tim0_2_4 = tim2_0_0 - tim2_1_0;
	  }
     }
     c_re(out[0]) = tre0_0_0 + tre0_0_1 + tre0_0_2 + tre0_0_3 + tre0_0_4;
     c_im(out[0]) = tim0_0_0 + tim0_0_1 + tim0_0_2 + tim0_0_3 + tim0_0_4;
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_0_1 + tre0_0_4)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_0_2 + tre0_0_3));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tim0_0_4 - tim0_0_1)) + (((FFTW_REAL) FFTW_K587785252) * (tim0_0_3 - tim0_0_2));
	  c_re(out[6 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[9 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_0_1 + tim0_0_4)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_0_2 + tim0_0_3));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tre0_0_1 - tre0_0_4)) + (((FFTW_REAL) FFTW_K587785252) * (tre0_0_2 - tre0_0_3));
	  c_im(out[6 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[9 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_0_2 + tre0_0_3)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_0_1 + tre0_0_4));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tim0_0_4 - tim0_0_1)) + (((FFTW_REAL) FFTW_K951056516) * (tim0_0_2 - tim0_0_3));
	  c_re(out[12 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[3 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_0_2 + tim0_0_3)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_0_1 + tim0_0_4));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tre0_0_1 - tre0_0_4)) + (((FFTW_REAL) FFTW_K951056516) * (tre0_0_3 - tre0_0_2));
	  c_im(out[12 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[3 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     c_re(out[10 * ostride]) = tre0_1_0 + tre0_1_1 + tre0_1_2 + tre0_1_3 + tre0_1_4;
     c_im(out[10 * ostride]) = tim0_1_0 + tim0_1_1 + tim0_1_2 + tim0_1_3 + tim0_1_4;
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_1_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_1_1 + tre0_1_4)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_1_2 + tre0_1_3));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tim0_1_4 - tim0_1_1)) + (((FFTW_REAL) FFTW_K587785252) * (tim0_1_3 - tim0_1_2));
	  c_re(out[ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[4 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_1_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_1_1 + tim0_1_4)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_1_2 + tim0_1_3));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tre0_1_1 - tre0_1_4)) + (((FFTW_REAL) FFTW_K587785252) * (tre0_1_2 - tre0_1_3));
	  c_im(out[ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[4 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_1_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_1_2 + tre0_1_3)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_1_1 + tre0_1_4));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tim0_1_4 - tim0_1_1)) + (((FFTW_REAL) FFTW_K951056516) * (tim0_1_2 - tim0_1_3));
	  c_re(out[7 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[13 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_1_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_1_2 + tim0_1_3)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_1_1 + tim0_1_4));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tre0_1_1 - tre0_1_4)) + (((FFTW_REAL) FFTW_K951056516) * (tre0_1_3 - tre0_1_2));
	  c_im(out[7 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[13 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     c_re(out[5 * ostride]) = tre0_2_0 + tre0_2_1 + tre0_2_2 + tre0_2_3 + tre0_2_4;
     c_im(out[5 * ostride]) = tim0_2_0 + tim0_2_1 + tim0_2_2 + tim0_2_3 + tim0_2_4;
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_2_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_2_1 + tre0_2_4)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_2_2 + tre0_2_3));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tim0_2_4 - tim0_2_1)) + (((FFTW_REAL) FFTW_K587785252) * (tim0_2_3 - tim0_2_2));
	  c_re(out[11 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[14 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_2_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_2_1 + tim0_2_4)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_2_2 + tim0_2_3));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tre0_2_1 - tre0_2_4)) + (((FFTW_REAL) FFTW_K587785252) * (tre0_2_2 - tre0_2_3));
	  c_im(out[11 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[14 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_2_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_2_2 + tre0_2_3)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_2_1 + tre0_2_4));
	  tre2_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tim0_2_4 - tim0_2_1)) + (((FFTW_REAL) FFTW_K951056516) * (tim0_2_2 - tim0_2_3));
	  c_re(out[2 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[8 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_2_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_2_2 + tim0_2_3)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_2_1 + tim0_2_4));
	  tim2_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tre0_2_1 - tre0_2_4)) + (((FFTW_REAL) FFTW_K951056516) * (tre0_2_3 - tre0_2_2));
	  c_im(out[2 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[8 * ostride]) = tim2_0_0 - tim2_1_0;
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 144 FP additions and 24 FP multiplications */

void fftwi_no_twiddle_16(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_0_1;
     FFTW_REAL tim0_0_1;
     FFTW_REAL tre0_0_2;
     FFTW_REAL tim0_0_2;
     FFTW_REAL tre0_0_3;
     FFTW_REAL tim0_0_3;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_1_1;
     FFTW_REAL tim0_1_1;
     FFTW_REAL tre0_1_2;
     FFTW_REAL tim0_1_2;
     FFTW_REAL tre0_1_3;
     FFTW_REAL tim0_1_3;
     FFTW_REAL tre0_2_0;
     FFTW_REAL tim0_2_0;
     FFTW_REAL tre0_2_1;
     FFTW_REAL tim0_2_1;
     FFTW_REAL tre0_2_2;
     FFTW_REAL tim0_2_2;
     FFTW_REAL tre0_2_3;
     FFTW_REAL tim0_2_3;
     FFTW_REAL tre0_3_0;
     FFTW_REAL tim0_3_0;
     FFTW_REAL tre0_3_1;
     FFTW_REAL tim0_3_1;
     FFTW_REAL tre0_3_2;
     FFTW_REAL tim0_3_2;
     FFTW_REAL tre0_3_3;
     FFTW_REAL tim0_3_3;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[0]);
	       tim2_0_0 = c_im(in[0]);
	       tre2_1_0 = c_re(in[8 * istride]);
	       tim2_1_0 = c_im(in[8 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[4 * istride]);
	       tim2_0_0 = c_im(in[4 * istride]);
	       tre2_1_0 = c_re(in[12 * istride]);
	       tim2_1_0 = c_im(in[12 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  tre0_0_0 = tre1_0_0 + tre1_0_1;
	  tim0_0_0 = tim1_0_0 + tim1_0_1;
	  tre0_2_0 = tre1_0_0 - tre1_0_1;
	  tim0_2_0 = tim1_0_0 - tim1_0_1;
	  tre0_1_0 = tre1_1_0 - tim1_1_1;
	  tim0_1_0 = tim1_1_0 + tre1_1_1;
	  tre0_3_0 = tre1_1_0 + tim1_1_1;
	  tim0_3_0 = tim1_1_0 - tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[istride]);
	       tim2_0_0 = c_im(in[istride]);
	       tre2_1_0 = c_re(in[9 * istride]);
	       tim2_1_0 = c_im(in[9 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[5 * istride]);
	       tim2_0_0 = c_im(in[5 * istride]);
	       tre2_1_0 = c_re(in[13 * istride]);
	       tim2_1_0 = c_im(in[13 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  tre0_0_1 = tre1_0_0 + tre1_0_1;
	  tim0_0_1 = tim1_0_0 + tim1_0_1;
	  tre0_2_1 = tre1_0_0 - tre1_0_1;
	  tim0_2_1 = tim1_0_0 - tim1_0_1;
	  tre0_1_1 = tre1_1_0 - tim1_1_1;
	  tim0_1_1 = tim1_1_0 + tre1_1_1;
	  tre0_3_1 = tre1_1_0 + tim1_1_1;
	  tim0_3_1 = tim1_1_0 - tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[2 * istride]);
	       tim2_0_0 = c_im(in[2 * istride]);
	       tre2_1_0 = c_re(in[10 * istride]);
	       tim2_1_0 = c_im(in[10 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[6 * istride]);
	       tim2_0_0 = c_im(in[6 * istride]);
	       tre2_1_0 = c_re(in[14 * istride]);
	       tim2_1_0 = c_im(in[14 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  tre0_0_2 = tre1_0_0 + tre1_0_1;
	  tim0_0_2 = tim1_0_0 + tim1_0_1;
	  tre0_2_2 = tre1_0_0 - tre1_0_1;
	  tim0_2_2 = tim1_0_0 - tim1_0_1;
	  tre0_1_2 = tre1_1_0 - tim1_1_1;
	  tim0_1_2 = tim1_1_0 + tre1_1_1;
	  tre0_3_2 = tre1_1_0 + tim1_1_1;
	  tim0_3_2 = tim1_1_0 - tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[3 * istride]);
	       tim2_0_0 = c_im(in[3 * istride]);
	       tre2_1_0 = c_re(in[11 * istride]);
	       tim2_1_0 = c_im(in[11 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[7 * istride]);
	       tim2_0_0 = c_im(in[7 * istride]);
	       tre2_1_0 = c_re(in[15 * istride]);
	       tim2_1_0 = c_im(in[15 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  tre0_0_3 = tre1_0_0 + tre1_0_1;
	  tim0_0_3 = tim1_0_0 + tim1_0_1;
	  tre0_2_3 = tre1_0_0 - tre1_0_1;
	  tim0_2_3 = tim1_0_0 - tim1_0_1;
	  tre0_1_3 = tre1_1_0 - tim1_1_1;
	  tim0_1_3 = tim1_1_0 + tre1_1_1;
	  tre0_3_3 = tre1_1_0 + tim1_1_1;
	  tim0_3_3 = tim1_1_0 - tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  tre1_0_0 = tre0_0_0 + tre0_0_2;
	  tim1_0_0 = tim0_0_0 + tim0_0_2;
	  tre1_1_0 = tre0_0_0 - tre0_0_2;
	  tim1_1_0 = tim0_0_0 - tim0_0_2;
	  tre1_0_1 = tre0_0_1 + tre0_0_3;
	  tim1_0_1 = tim0_0_1 + tim0_0_3;
	  tre1_1_1 = tre0_0_1 - tre0_0_3;
	  tim1_1_1 = tim0_0_1 - tim0_0_3;
	  c_re(out[0]) = tre1_0_0 + tre1_0_1;
	  c_im(out[0]) = tim1_0_0 + tim1_0_1;
	  c_re(out[8 * ostride]) = tre1_0_0 - tre1_0_1;
	  c_im(out[8 * ostride]) = tim1_0_0 - tim1_0_1;
	  c_re(out[4 * ostride]) = tre1_1_0 - tim1_1_1;
	  c_im(out[4 * ostride]) = tim1_1_0 + tre1_1_1;
	  c_re(out[12 * ostride]) = tre1_1_0 + tim1_1_1;
	  c_im(out[12 * ostride]) = tim1_1_0 - tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  {
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_1_2 - tim0_1_2);
	       tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_1_2 + tre0_1_2);
	       tre1_0_0 = tre0_1_0 + tre2_1_0;
	       tim1_0_0 = tim0_1_0 + tim2_1_0;
	       tre1_1_0 = tre0_1_0 - tre2_1_0;
	       tim1_1_0 = tim0_1_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_1_1) - (((FFTW_REAL) FFTW_K382683432) * tim0_1_1);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_1_1) + (((FFTW_REAL) FFTW_K382683432) * tre0_1_1);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_1_3) - (((FFTW_REAL) FFTW_K923879532) * tim0_1_3);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_1_3) + (((FFTW_REAL) FFTW_K923879532) * tre0_1_3);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  c_re(out[ostride]) = tre1_0_0 + tre1_0_1;
	  c_im(out[ostride]) = tim1_0_0 + tim1_0_1;
	  c_re(out[9 * ostride]) = tre1_0_0 - tre1_0_1;
	  c_im(out[9 * ostride]) = tim1_0_0 - tim1_0_1;
	  c_re(out[5 * ostride]) = tre1_1_0 - tim1_1_1;
	  c_im(out[5 * ostride]) = tim1_1_0 + tre1_1_1;
	  c_re(out[13 * ostride]) = tre1_1_0 + tim1_1_1;
	  c_im(out[13 * ostride]) = tim1_1_0 - tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  tre1_0_0 = tre0_2_0 - tim0_2_2;
	  tim1_0_0 = tim0_2_0 + tre0_2_2;
	  tre1_1_0 = tre0_2_0 + tim0_2_2;
	  tim1_1_0 = tim0_2_0 - tre0_2_2;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_2_1 - tim0_2_1);
	       tim2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_2_1 + tre0_2_1);
	       tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_2_3 + tim0_2_3);
	       tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_2_3 - tim0_2_3);
	       tre1_0_1 = tre2_0_0 - tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 + tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  c_re(out[2 * ostride]) = tre1_0_0 + tre1_0_1;
	  c_im(out[2 * ostride]) = tim1_0_0 + tim1_0_1;
	  c_re(out[10 * ostride]) = tre1_0_0 - tre1_0_1;
	  c_im(out[10 * ostride]) = tim1_0_0 - tim1_0_1;
	  c_re(out[6 * ostride]) = tre1_1_0 - tim1_1_1;
	  c_im(out[6 * ostride]) = tim1_1_0 + tre1_1_1;
	  c_re(out[14 * ostride]) = tre1_1_0 + tim1_1_1;
	  c_im(out[14 * ostride]) = tim1_1_0 - tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  {
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_3_2 + tim0_3_2);
	       tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_3_2 - tim0_3_2);
	       tre1_0_0 = tre0_3_0 - tre2_1_0;
	       tim1_0_0 = tim0_3_0 + tim2_1_0;
	       tre1_1_0 = tre0_3_0 + tre2_1_0;
	       tim1_1_0 = tim0_3_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_3_1) - (((FFTW_REAL) FFTW_K923879532) * tim0_3_1);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_3_1) + (((FFTW_REAL) FFTW_K923879532) * tre0_3_1);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_3_3) - (((FFTW_REAL) FFTW_K923879532) * tre0_3_3);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_3_3) + (((FFTW_REAL) FFTW_K382683432) * tre0_3_3);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 - tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 + tim2_1_0;
	  }
	  c_re(out[3 * ostride]) = tre1_0_0 + tre1_0_1;
	  c_im(out[3 * ostride]) = tim1_0_0 + tim1_0_1;
	  c_re(out[11 * ostride]) = tre1_0_0 - tre1_0_1;
	  c_im(out[11 * ostride]) = tim1_0_0 - tim1_0_1;
	  c_re(out[7 * ostride]) = tre1_1_0 - tim1_1_1;
	  c_im(out[7 * ostride]) = tim1_1_0 + tre1_1_1;
	  c_re(out[15 * ostride]) = tre1_1_0 + tim1_1_1;
	  c_im(out[15 * ostride]) = tim1_1_0 - tre1_1_1;
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 4 FP additions and 0 FP multiplications */

void fftwi_no_twiddle_2(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     tre0_0_0 = c_re(in[0]);
     tim0_0_0 = c_im(in[0]);
     tre0_1_0 = c_re(in[istride]);
     tim0_1_0 = c_im(in[istride]);
     c_re(out[0]) = tre0_0_0 + tre0_1_0;
     c_im(out[0]) = tim0_0_0 + tim0_1_0;
     c_re(out[ostride]) = tre0_0_0 - tre0_1_0;
     c_im(out[ostride]) = tim0_0_0 - tim0_1_0;
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 14 FP additions and 4 FP multiplications */

void fftwi_no_twiddle_3(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_2_0;
     FFTW_REAL tim0_2_0;
     tre0_0_0 = c_re(in[0]);
     tim0_0_0 = c_im(in[0]);
     tre0_1_0 = c_re(in[istride]);
     tim0_1_0 = c_im(in[istride]);
     tre0_2_0 = c_re(in[2 * istride]);
     tim0_2_0 = c_im(in[2 * istride]);
     c_re(out[0]) = tre0_0_0 + tre0_1_0 + tre0_2_0;
     c_im(out[0]) = tim0_0_0 + tim0_1_0 + tim0_2_0;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre0_1_0 + tre0_2_0));
	  tre1_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim0_2_0 - tim0_1_0);
	  c_re(out[ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[2 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim0_1_0 + tim0_2_0));
	  tim1_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre0_1_0 - tre0_2_0);
	  c_im(out[ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[2 * ostride]) = tim1_0_0 - tim1_1_0;
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 376 FP additions and 88 FP multiplications */

void fftwi_no_twiddle_32(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_0_1;
     FFTW_REAL tim0_0_1;
     FFTW_REAL tre0_0_2;
     FFTW_REAL tim0_0_2;
     FFTW_REAL tre0_0_3;
     FFTW_REAL tim0_0_3;
     FFTW_REAL tre0_0_4;
     FFTW_REAL tim0_0_4;
     FFTW_REAL tre0_0_5;
     FFTW_REAL tim0_0_5;
     FFTW_REAL tre0_0_6;
     FFTW_REAL tim0_0_6;
     FFTW_REAL tre0_0_7;
     FFTW_REAL tim0_0_7;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_1_1;
     FFTW_REAL tim0_1_1;
     FFTW_REAL tre0_1_2;
     FFTW_REAL tim0_1_2;
     FFTW_REAL tre0_1_3;
     FFTW_REAL tim0_1_3;
     FFTW_REAL tre0_1_4;
     FFTW_REAL tim0_1_4;
     FFTW_REAL tre0_1_5;
     FFTW_REAL tim0_1_5;
     FFTW_REAL tre0_1_6;
     FFTW_REAL tim0_1_6;
     FFTW_REAL tre0_1_7;
     FFTW_REAL tim0_1_7;
     FFTW_REAL tre0_2_0;
     FFTW_REAL tim0_2_0;
     FFTW_REAL tre0_2_1;
     FFTW_REAL tim0_2_1;
     FFTW_REAL tre0_2_2;
     FFTW_REAL tim0_2_2;
     FFTW_REAL tre0_2_3;
     FFTW_REAL tim0_2_3;
     FFTW_REAL tre0_2_4;
     FFTW_REAL tim0_2_4;
     FFTW_REAL tre0_2_5;
     FFTW_REAL tim0_2_5;
     FFTW_REAL tre0_2_6;
     FFTW_REAL tim0_2_6;
     FFTW_REAL tre0_2_7;
     FFTW_REAL tim0_2_7;
     FFTW_REAL tre0_3_0;
     FFTW_REAL tim0_3_0;
     FFTW_REAL tre0_3_1;
     FFTW_REAL tim0_3_1;
     FFTW_REAL tre0_3_2;
     FFTW_REAL tim0_3_2;
     FFTW_REAL tre0_3_3;
     FFTW_REAL tim0_3_3;
     FFTW_REAL tre0_3_4;
     FFTW_REAL tim0_3_4;
     FFTW_REAL tre0_3_5;
     FFTW_REAL tim0_3_5;
     FFTW_REAL tre0_3_6;
     FFTW_REAL tim0_3_6;
     FFTW_REAL tre0_3_7;
     FFTW_REAL tim0_3_7;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[0]);
	       tim2_0_0 = c_im(in[0]);
	       tre2_1_0 = c_re(in[16 * istride]);
	       tim2_1_0 = c_im(in[16 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[8 * istride]);
	       tim2_0_0 = c_im(in[8 * istride]);
	       tre2_1_0 = c_re(in[24 * istride]);
	       tim2_1_0 = c_im(in[24 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  tre0_0_0 = tre1_0_0 + tre1_0_1;
	  tim0_0_0 = tim1_0_0 + tim1_0_1;
	  tre0_2_0 = tre1_0_0 - tre1_0_1;
	  tim0_2_0 = tim1_0_0 - tim1_0_1;
	  tre0_1_0 = tre1_1_0 - tim1_1_1;
	  tim0_1_0 = tim1_1_0 + tre1_1_1;
	  tre0_3_0 = tre1_1_0 + tim1_1_1;
	  tim0_3_0 = tim1_1_0 - tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[istride]);
	       tim2_0_0 = c_im(in[istride]);
	       tre2_1_0 = c_re(in[17 * istride]);
	       tim2_1_0 = c_im(in[17 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[9 * istride]);
	       tim2_0_0 = c_im(in[9 * istride]);
	       tre2_1_0 = c_re(in[25 * istride]);
	       tim2_1_0 = c_im(in[25 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  tre0_0_1 = tre1_0_0 + tre1_0_1;
	  tim0_0_1 = tim1_0_0 + tim1_0_1;
	  tre0_2_1 = tre1_0_0 - tre1_0_1;
	  tim0_2_1 = tim1_0_0 - tim1_0_1;
	  tre0_1_1 = tre1_1_0 - tim1_1_1;
	  tim0_1_1 = tim1_1_0 + tre1_1_1;
	  tre0_3_1 = tre1_1_0 + tim1_1_1;
	  tim0_3_1 = tim1_1_0 - tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[2 * istride]);
	       tim2_0_0 = c_im(in[2 * istride]);
	       tre2_1_0 = c_re(in[18 * istride]);
	       tim2_1_0 = c_im(in[18 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[10 * istride]);
	       tim2_0_0 = c_im(in[10 * istride]);
	       tre2_1_0 = c_re(in[26 * istride]);
	       tim2_1_0 = c_im(in[26 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  tre0_0_2 = tre1_0_0 + tre1_0_1;
	  tim0_0_2 = tim1_0_0 + tim1_0_1;
	  tre0_2_2 = tre1_0_0 - tre1_0_1;
	  tim0_2_2 = tim1_0_0 - tim1_0_1;
	  tre0_1_2 = tre1_1_0 - tim1_1_1;
	  tim0_1_2 = tim1_1_0 + tre1_1_1;
	  tre0_3_2 = tre1_1_0 + tim1_1_1;
	  tim0_3_2 = tim1_1_0 - tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[3 * istride]);
	       tim2_0_0 = c_im(in[3 * istride]);
	       tre2_1_0 = c_re(in[19 * istride]);
	       tim2_1_0 = c_im(in[19 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[11 * istride]);
	       tim2_0_0 = c_im(in[11 * istride]);
	       tre2_1_0 = c_re(in[27 * istride]);
	       tim2_1_0 = c_im(in[27 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  tre0_0_3 = tre1_0_0 + tre1_0_1;
	  tim0_0_3 = tim1_0_0 + tim1_0_1;
	  tre0_2_3 = tre1_0_0 - tre1_0_1;
	  tim0_2_3 = tim1_0_0 - tim1_0_1;
	  tre0_1_3 = tre1_1_0 - tim1_1_1;
	  tim0_1_3 = tim1_1_0 + tre1_1_1;
	  tre0_3_3 = tre1_1_0 + tim1_1_1;
	  tim0_3_3 = tim1_1_0 - tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[4 * istride]);
	       tim2_0_0 = c_im(in[4 * istride]);
	       tre2_1_0 = c_re(in[20 * istride]);
	       tim2_1_0 = c_im(in[20 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[12 * istride]);
	       tim2_0_0 = c_im(in[12 * istride]);
	       tre2_1_0 = c_re(in[28 * istride]);
	       tim2_1_0 = c_im(in[28 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  tre0_0_4 = tre1_0_0 + tre1_0_1;
	  tim0_0_4 = tim1_0_0 + tim1_0_1;
	  tre0_2_4 = tre1_0_0 - tre1_0_1;
	  tim0_2_4 = tim1_0_0 - tim1_0_1;
	  tre0_1_4 = tre1_1_0 - tim1_1_1;
	  tim0_1_4 = tim1_1_0 + tre1_1_1;
	  tre0_3_4 = tre1_1_0 + tim1_1_1;
	  tim0_3_4 = tim1_1_0 - tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[5 * istride]);
	       tim2_0_0 = c_im(in[5 * istride]);
	       tre2_1_0 = c_re(in[21 * istride]);
	       tim2_1_0 = c_im(in[21 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[13 * istride]);
	       tim2_0_0 = c_im(in[13 * istride]);
	       tre2_1_0 = c_re(in[29 * istride]);
	       tim2_1_0 = c_im(in[29 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  tre0_0_5 = tre1_0_0 + tre1_0_1;
	  tim0_0_5 = tim1_0_0 + tim1_0_1;
	  tre0_2_5 = tre1_0_0 - tre1_0_1;
	  tim0_2_5 = tim1_0_0 - tim1_0_1;
	  tre0_1_5 = tre1_1_0 - tim1_1_1;
	  tim0_1_5 = tim1_1_0 + tre1_1_1;
	  tre0_3_5 = tre1_1_0 + tim1_1_1;
	  tim0_3_5 = tim1_1_0 - tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[6 * istride]);
	       tim2_0_0 = c_im(in[6 * istride]);
	       tre2_1_0 = c_re(in[22 * istride]);
	       tim2_1_0 = c_im(in[22 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[14 * istride]);
	       tim2_0_0 = c_im(in[14 * istride]);
	       tre2_1_0 = c_re(in[30 * istride]);
	       tim2_1_0 = c_im(in[30 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  tre0_0_6 = tre1_0_0 + tre1_0_1;
	  tim0_0_6 = tim1_0_0 + tim1_0_1;
	  tre0_2_6 = tre1_0_0 - tre1_0_1;
	  tim0_2_6 = tim1_0_0 - tim1_0_1;
	  tre0_1_6 = tre1_1_0 - tim1_1_1;
	  tim0_1_6 = tim1_1_0 + tre1_1_1;
	  tre0_3_6 = tre1_1_0 + tim1_1_1;
	  tim0_3_6 = tim1_1_0 - tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[7 * istride]);
	       tim2_0_0 = c_im(in[7 * istride]);
	       tre2_1_0 = c_re(in[23 * istride]);
	       tim2_1_0 = c_im(in[23 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[15 * istride]);
	       tim2_0_0 = c_im(in[15 * istride]);
	       tre2_1_0 = c_re(in[31 * istride]);
	       tim2_1_0 = c_im(in[31 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  tre0_0_7 = tre1_0_0 + tre1_0_1;
	  tim0_0_7 = tim1_0_0 + tim1_0_1;
	  tre0_2_7 = tre1_0_0 - tre1_0_1;
	  tim0_2_7 = tim1_0_0 - tim1_0_1;
	  tre0_1_7 = tre1_1_0 - tim1_1_1;
	  tim0_1_7 = tim1_1_0 + tre1_1_1;
	  tre0_3_7 = tre1_1_0 + tim1_1_1;
	  tim0_3_7 = tim1_1_0 - tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  tre1_0_0 = tre0_0_0 + tre0_0_4;
	  tim1_0_0 = tim0_0_0 + tim0_0_4;
	  tre1_1_0 = tre0_0_0 - tre0_0_4;
	  tim1_1_0 = tim0_0_0 - tim0_0_4;
	  tre1_0_1 = tre0_0_1 + tre0_0_5;
	  tim1_0_1 = tim0_0_1 + tim0_0_5;
	  tre1_1_1 = tre0_0_1 - tre0_0_5;
	  tim1_1_1 = tim0_0_1 - tim0_0_5;
	  tre1_0_2 = tre0_0_2 + tre0_0_6;
	  tim1_0_2 = tim0_0_2 + tim0_0_6;
	  tre1_1_2 = tre0_0_2 - tre0_0_6;
	  tim1_1_2 = tim0_0_2 - tim0_0_6;
	  tre1_0_3 = tre0_0_3 + tre0_0_7;
	  tim1_0_3 = tim0_0_3 + tim0_0_7;
	  tre1_1_3 = tre0_0_3 - tre0_0_7;
	  tim1_1_3 = tim0_0_3 - tim0_0_7;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       c_re(out[0]) = tre2_0_0 + tre2_0_1;
	       c_im(out[0]) = tim2_0_0 + tim2_0_1;
	       c_re(out[16 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[16 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[8 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[8 * ostride]) = tim2_1_0 + tre2_1_1;
	       c_re(out[24 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[24 * ostride]) = tim2_1_0 - tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 - tim1_1_2;
	       tim2_0_0 = tim1_1_0 + tre1_1_2;
	       tre2_1_0 = tre1_1_0 + tim1_1_2;
	       tim2_1_0 = tim1_1_0 - tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
		    tre2_0_1 = tre3_0_0 - tre3_1_0;
		    tim2_0_1 = tim3_0_0 + tim3_1_0;
		    tre2_1_1 = tre3_0_0 + tre3_1_0;
		    tim2_1_1 = tim3_0_0 - tim3_1_0;
	       }
	       c_re(out[4 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[4 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[20 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[20 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[12 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[12 * ostride]) = tim2_1_0 + tre2_1_1;
	       c_re(out[28 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[28 * ostride]) = tim2_1_0 - tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_1_4 - tim0_1_4);
	       tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_1_4 + tre0_1_4);
	       tre1_0_0 = tre0_1_0 + tre2_1_0;
	       tim1_0_0 = tim0_1_0 + tim2_1_0;
	       tre1_1_0 = tre0_1_0 - tre2_1_0;
	       tim1_1_0 = tim0_1_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tre0_1_1) - (((FFTW_REAL) FFTW_K195090322) * tim0_1_1);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tim0_1_1) + (((FFTW_REAL) FFTW_K195090322) * tre0_1_1);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K555570233) * tre0_1_5) - (((FFTW_REAL) FFTW_K831469612) * tim0_1_5);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K555570233) * tim0_1_5) + (((FFTW_REAL) FFTW_K831469612) * tre0_1_5);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_1_2) - (((FFTW_REAL) FFTW_K382683432) * tim0_1_2);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_1_2) + (((FFTW_REAL) FFTW_K382683432) * tre0_1_2);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_1_6) - (((FFTW_REAL) FFTW_K923879532) * tim0_1_6);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_1_6) + (((FFTW_REAL) FFTW_K923879532) * tre0_1_6);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_1_3) - (((FFTW_REAL) FFTW_K555570233) * tim0_1_3);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_1_3) + (((FFTW_REAL) FFTW_K555570233) * tre0_1_3);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tre0_1_7) - (((FFTW_REAL) FFTW_K980785280) * tim0_1_7);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tim0_1_7) + (((FFTW_REAL) FFTW_K980785280) * tre0_1_7);
	       tre1_0_3 = tre2_0_0 + tre2_1_0;
	       tim1_0_3 = tim2_0_0 + tim2_1_0;
	       tre1_1_3 = tre2_0_0 - tre2_1_0;
	       tim1_1_3 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       c_re(out[ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[17 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[17 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[9 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[9 * ostride]) = tim2_1_0 + tre2_1_1;
	       c_re(out[25 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[25 * ostride]) = tim2_1_0 - tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 - tim1_1_2;
	       tim2_0_0 = tim1_1_0 + tre1_1_2;
	       tre2_1_0 = tre1_1_0 + tim1_1_2;
	       tim2_1_0 = tim1_1_0 - tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
		    tre2_0_1 = tre3_0_0 - tre3_1_0;
		    tim2_0_1 = tim3_0_0 + tim3_1_0;
		    tre2_1_1 = tre3_0_0 + tre3_1_0;
		    tim2_1_1 = tim3_0_0 - tim3_1_0;
	       }
	       c_re(out[5 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[5 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[21 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[21 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[13 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[13 * ostride]) = tim2_1_0 + tre2_1_1;
	       c_re(out[29 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[29 * ostride]) = tim2_1_0 - tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  tre1_0_0 = tre0_2_0 - tim0_2_4;
	  tim1_0_0 = tim0_2_0 + tre0_2_4;
	  tre1_1_0 = tre0_2_0 + tim0_2_4;
	  tim1_1_0 = tim0_2_0 - tre0_2_4;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_2_1) - (((FFTW_REAL) FFTW_K382683432) * tim0_2_1);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_2_1) + (((FFTW_REAL) FFTW_K382683432) * tre0_2_1);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_2_5) + (((FFTW_REAL) FFTW_K923879532) * tim0_2_5);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_2_5) - (((FFTW_REAL) FFTW_K382683432) * tim0_2_5);
	       tre1_0_1 = tre2_0_0 - tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 + tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_2_2 - tim0_2_2);
	       tim2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_2_2 + tre0_2_2);
	       tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_2_6 + tim0_2_6);
	       tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_2_6 - tim0_2_6);
	       tre1_0_2 = tre2_0_0 - tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 + tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_2_3) - (((FFTW_REAL) FFTW_K923879532) * tim0_2_3);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_2_3) + (((FFTW_REAL) FFTW_K923879532) * tre0_2_3);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_2_7) + (((FFTW_REAL) FFTW_K382683432) * tim0_2_7);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_2_7) - (((FFTW_REAL) FFTW_K923879532) * tim0_2_7);
	       tre1_0_3 = tre2_0_0 - tre2_1_0;
	       tim1_0_3 = tim2_0_0 + tim2_1_0;
	       tre1_1_3 = tre2_0_0 + tre2_1_0;
	       tim1_1_3 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       c_re(out[2 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[2 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[18 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[18 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[10 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[10 * ostride]) = tim2_1_0 + tre2_1_1;
	       c_re(out[26 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[26 * ostride]) = tim2_1_0 - tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 - tim1_1_2;
	       tim2_0_0 = tim1_1_0 + tre1_1_2;
	       tre2_1_0 = tre1_1_0 + tim1_1_2;
	       tim2_1_0 = tim1_1_0 - tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
		    tre2_0_1 = tre3_0_0 - tre3_1_0;
		    tim2_0_1 = tim3_0_0 + tim3_1_0;
		    tre2_1_1 = tre3_0_0 + tre3_1_0;
		    tim2_1_1 = tim3_0_0 - tim3_1_0;
	       }
	       c_re(out[6 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[6 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[22 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[22 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[14 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[14 * ostride]) = tim2_1_0 + tre2_1_1;
	       c_re(out[30 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[30 * ostride]) = tim2_1_0 - tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_3_4 + tim0_3_4);
	       tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_3_4 - tim0_3_4);
	       tre1_0_0 = tre0_3_0 - tre2_1_0;
	       tim1_0_0 = tim0_3_0 + tim2_1_0;
	       tre1_1_0 = tre0_3_0 + tre2_1_0;
	       tim1_1_0 = tim0_3_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_3_1) - (((FFTW_REAL) FFTW_K555570233) * tim0_3_1);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_3_1) + (((FFTW_REAL) FFTW_K555570233) * tre0_3_1);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K980785280) * tre0_3_5) + (((FFTW_REAL) FFTW_K195090322) * tim0_3_5);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tre0_3_5) - (((FFTW_REAL) FFTW_K980785280) * tim0_3_5);
	       tre1_0_1 = tre2_0_0 - tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 + tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_3_2) - (((FFTW_REAL) FFTW_K923879532) * tim0_3_2);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_3_2) + (((FFTW_REAL) FFTW_K923879532) * tre0_3_2);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_3_6) - (((FFTW_REAL) FFTW_K923879532) * tre0_3_6);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_3_6) + (((FFTW_REAL) FFTW_K382683432) * tre0_3_6);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 - tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 + tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K195090322) * tre0_3_3) + (((FFTW_REAL) FFTW_K980785280) * tim0_3_3);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tre0_3_3) - (((FFTW_REAL) FFTW_K195090322) * tim0_3_3);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_3_7) - (((FFTW_REAL) FFTW_K555570233) * tre0_3_7);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K555570233) * tim0_3_7) + (((FFTW_REAL) FFTW_K831469612) * tre0_3_7);
	       tre1_0_3 = tre2_1_0 - tre2_0_0;
	       tim1_0_3 = tim2_0_0 - tim2_1_0;
	       tre1_1_3 = (-(tre2_0_0 + tre2_1_0));
	       tim1_1_3 = tim2_0_0 + tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       c_re(out[3 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[3 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[19 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[19 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[11 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[11 * ostride]) = tim2_1_0 + tre2_1_1;
	       c_re(out[27 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[27 * ostride]) = tim2_1_0 - tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 - tim1_1_2;
	       tim2_0_0 = tim1_1_0 + tre1_1_2;
	       tre2_1_0 = tre1_1_0 + tim1_1_2;
	       tim2_1_0 = tim1_1_0 - tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
		    tre2_0_1 = tre3_0_0 - tre3_1_0;
		    tim2_0_1 = tim3_0_0 + tim3_1_0;
		    tre2_1_1 = tre3_0_0 + tre3_1_0;
		    tim2_1_1 = tim3_0_0 - tim3_1_0;
	       }
	       c_re(out[7 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[7 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[23 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[23 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[15 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[15 * ostride]) = tim2_1_0 + tre2_1_1;
	       c_re(out[31 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[31 * ostride]) = tim2_1_0 - tre2_1_1;
	  }
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 16 FP additions and 0 FP multiplications */

void fftwi_no_twiddle_4(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_0_1;
     FFTW_REAL tim0_0_1;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_1_1;
     FFTW_REAL tim0_1_1;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[0]);
	  tim1_0_0 = c_im(in[0]);
	  tre1_1_0 = c_re(in[2 * istride]);
	  tim1_1_0 = c_im(in[2 * istride]);
	  tre0_0_0 = tre1_0_0 + tre1_1_0;
	  tim0_0_0 = tim1_0_0 + tim1_1_0;
	  tre0_1_0 = tre1_0_0 - tre1_1_0;
	  tim0_1_0 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[istride]);
	  tim1_0_0 = c_im(in[istride]);
	  tre1_1_0 = c_re(in[3 * istride]);
	  tim1_1_0 = c_im(in[3 * istride]);
	  tre0_0_1 = tre1_0_0 + tre1_1_0;
	  tim0_0_1 = tim1_0_0 + tim1_1_0;
	  tre0_1_1 = tre1_0_0 - tre1_1_0;
	  tim0_1_1 = tim1_0_0 - tim1_1_0;
     }
     c_re(out[0]) = tre0_0_0 + tre0_0_1;
     c_im(out[0]) = tim0_0_0 + tim0_0_1;
     c_re(out[2 * ostride]) = tre0_0_0 - tre0_0_1;
     c_im(out[2 * ostride]) = tim0_0_0 - tim0_0_1;
     c_re(out[ostride]) = tre0_1_0 - tim0_1_1;
     c_im(out[ostride]) = tim0_1_0 + tre0_1_1;
     c_re(out[3 * ostride]) = tre0_1_0 + tim0_1_1;
     c_im(out[3 * ostride]) = tim0_1_0 - tre0_1_1;
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 44 FP additions and 16 FP multiplications */

void fftwi_no_twiddle_5(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_2_0;
     FFTW_REAL tim0_2_0;
     FFTW_REAL tre0_3_0;
     FFTW_REAL tim0_3_0;
     FFTW_REAL tre0_4_0;
     FFTW_REAL tim0_4_0;
     tre0_0_0 = c_re(in[0]);
     tim0_0_0 = c_im(in[0]);
     tre0_1_0 = c_re(in[istride]);
     tim0_1_0 = c_im(in[istride]);
     tre0_2_0 = c_re(in[2 * istride]);
     tim0_2_0 = c_im(in[2 * istride]);
     tre0_3_0 = c_re(in[3 * istride]);
     tim0_3_0 = c_im(in[3 * istride]);
     tre0_4_0 = c_re(in[4 * istride]);
     tim0_4_0 = c_im(in[4 * istride]);
     c_re(out[0]) = tre0_0_0 + tre0_1_0 + tre0_2_0 + tre0_3_0 + tre0_4_0;
     c_im(out[0]) = tim0_0_0 + tim0_1_0 + tim0_2_0 + tim0_3_0 + tim0_4_0;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_1_0 + tre0_4_0)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_2_0 + tre0_3_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tim0_4_0 - tim0_1_0)) + (((FFTW_REAL) FFTW_K587785252) * (tim0_3_0 - tim0_2_0));
	  c_re(out[ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[4 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_1_0 + tim0_4_0)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_2_0 + tim0_3_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tre0_1_0 - tre0_4_0)) + (((FFTW_REAL) FFTW_K587785252) * (tre0_2_0 - tre0_3_0));
	  c_im(out[ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[4 * ostride]) = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_2_0 + tre0_3_0)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_1_0 + tre0_4_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tim0_4_0 - tim0_1_0)) + (((FFTW_REAL) FFTW_K951056516) * (tim0_2_0 - tim0_3_0));
	  c_re(out[2 * ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[3 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_2_0 + tim0_3_0)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_1_0 + tim0_4_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tre0_1_0 - tre0_4_0)) + (((FFTW_REAL) FFTW_K951056516) * (tre0_3_0 - tre0_2_0));
	  c_im(out[2 * ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[3 * ostride]) = tim1_0_0 - tim1_1_0;
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 40 FP additions and 8 FP multiplications */

void fftwi_no_twiddle_6(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_0_1;
     FFTW_REAL tim0_0_1;
     FFTW_REAL tre0_0_2;
     FFTW_REAL tim0_0_2;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_1_1;
     FFTW_REAL tim0_1_1;
     FFTW_REAL tre0_1_2;
     FFTW_REAL tim0_1_2;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[0]);
	  tim1_0_0 = c_im(in[0]);
	  tre1_1_0 = c_re(in[3 * istride]);
	  tim1_1_0 = c_im(in[3 * istride]);
	  tre0_0_0 = tre1_0_0 + tre1_1_0;
	  tim0_0_0 = tim1_0_0 + tim1_1_0;
	  tre0_1_0 = tre1_0_0 - tre1_1_0;
	  tim0_1_0 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[2 * istride]);
	  tim1_0_0 = c_im(in[2 * istride]);
	  tre1_1_0 = c_re(in[5 * istride]);
	  tim1_1_0 = c_im(in[5 * istride]);
	  tre0_0_1 = tre1_0_0 + tre1_1_0;
	  tim0_0_1 = tim1_0_0 + tim1_1_0;
	  tre0_1_1 = tre1_0_0 - tre1_1_0;
	  tim0_1_1 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[4 * istride]);
	  tim1_0_0 = c_im(in[4 * istride]);
	  tre1_1_0 = c_re(in[istride]);
	  tim1_1_0 = c_im(in[istride]);
	  tre0_0_2 = tre1_0_0 + tre1_1_0;
	  tim0_0_2 = tim1_0_0 + tim1_1_0;
	  tre0_1_2 = tre1_0_0 - tre1_1_0;
	  tim0_1_2 = tim1_0_0 - tim1_1_0;
     }
     c_re(out[0]) = tre0_0_0 + tre0_0_1 + tre0_0_2;
     c_im(out[0]) = tim0_0_0 + tim0_0_1 + tim0_0_2;
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre0_0_1 + tre0_0_2));
	  tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim0_0_2 - tim0_0_1);
	  c_re(out[4 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[2 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim0_0_1 + tim0_0_2));
	  tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre0_0_1 - tre0_0_2);
	  c_im(out[4 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[2 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     c_re(out[3 * ostride]) = tre0_1_0 + tre0_1_1 + tre0_1_2;
     c_im(out[3 * ostride]) = tim0_1_0 + tim0_1_1 + tim0_1_2;
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_1_0 - (((FFTW_REAL) FFTW_K499999999) * (tre0_1_1 + tre0_1_2));
	  tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim0_1_2 - tim0_1_1);
	  c_re(out[ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[5 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_1_0 - (((FFTW_REAL) FFTW_K499999999) * (tim0_1_1 + tim0_1_2));
	  tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre0_1_1 - tre0_1_2);
	  c_im(out[ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[5 * ostride]) = tim2_0_0 - tim2_1_0;
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 928 FP additions and 248 FP multiplications */

void fftwi_no_twiddle_64(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_0_1;
     FFTW_REAL tim0_0_1;
     FFTW_REAL tre0_0_2;
     FFTW_REAL tim0_0_2;
     FFTW_REAL tre0_0_3;
     FFTW_REAL tim0_0_3;
     FFTW_REAL tre0_0_4;
     FFTW_REAL tim0_0_4;
     FFTW_REAL tre0_0_5;
     FFTW_REAL tim0_0_5;
     FFTW_REAL tre0_0_6;
     FFTW_REAL tim0_0_6;
     FFTW_REAL tre0_0_7;
     FFTW_REAL tim0_0_7;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_1_1;
     FFTW_REAL tim0_1_1;
     FFTW_REAL tre0_1_2;
     FFTW_REAL tim0_1_2;
     FFTW_REAL tre0_1_3;
     FFTW_REAL tim0_1_3;
     FFTW_REAL tre0_1_4;
     FFTW_REAL tim0_1_4;
     FFTW_REAL tre0_1_5;
     FFTW_REAL tim0_1_5;
     FFTW_REAL tre0_1_6;
     FFTW_REAL tim0_1_6;
     FFTW_REAL tre0_1_7;
     FFTW_REAL tim0_1_7;
     FFTW_REAL tre0_2_0;
     FFTW_REAL tim0_2_0;
     FFTW_REAL tre0_2_1;
     FFTW_REAL tim0_2_1;
     FFTW_REAL tre0_2_2;
     FFTW_REAL tim0_2_2;
     FFTW_REAL tre0_2_3;
     FFTW_REAL tim0_2_3;
     FFTW_REAL tre0_2_4;
     FFTW_REAL tim0_2_4;
     FFTW_REAL tre0_2_5;
     FFTW_REAL tim0_2_5;
     FFTW_REAL tre0_2_6;
     FFTW_REAL tim0_2_6;
     FFTW_REAL tre0_2_7;
     FFTW_REAL tim0_2_7;
     FFTW_REAL tre0_3_0;
     FFTW_REAL tim0_3_0;
     FFTW_REAL tre0_3_1;
     FFTW_REAL tim0_3_1;
     FFTW_REAL tre0_3_2;
     FFTW_REAL tim0_3_2;
     FFTW_REAL tre0_3_3;
     FFTW_REAL tim0_3_3;
     FFTW_REAL tre0_3_4;
     FFTW_REAL tim0_3_4;
     FFTW_REAL tre0_3_5;
     FFTW_REAL tim0_3_5;
     FFTW_REAL tre0_3_6;
     FFTW_REAL tim0_3_6;
     FFTW_REAL tre0_3_7;
     FFTW_REAL tim0_3_7;
     FFTW_REAL tre0_4_0;
     FFTW_REAL tim0_4_0;
     FFTW_REAL tre0_4_1;
     FFTW_REAL tim0_4_1;
     FFTW_REAL tre0_4_2;
     FFTW_REAL tim0_4_2;
     FFTW_REAL tre0_4_3;
     FFTW_REAL tim0_4_3;
     FFTW_REAL tre0_4_4;
     FFTW_REAL tim0_4_4;
     FFTW_REAL tre0_4_5;
     FFTW_REAL tim0_4_5;
     FFTW_REAL tre0_4_6;
     FFTW_REAL tim0_4_6;
     FFTW_REAL tre0_4_7;
     FFTW_REAL tim0_4_7;
     FFTW_REAL tre0_5_0;
     FFTW_REAL tim0_5_0;
     FFTW_REAL tre0_5_1;
     FFTW_REAL tim0_5_1;
     FFTW_REAL tre0_5_2;
     FFTW_REAL tim0_5_2;
     FFTW_REAL tre0_5_3;
     FFTW_REAL tim0_5_3;
     FFTW_REAL tre0_5_4;
     FFTW_REAL tim0_5_4;
     FFTW_REAL tre0_5_5;
     FFTW_REAL tim0_5_5;
     FFTW_REAL tre0_5_6;
     FFTW_REAL tim0_5_6;
     FFTW_REAL tre0_5_7;
     FFTW_REAL tim0_5_7;
     FFTW_REAL tre0_6_0;
     FFTW_REAL tim0_6_0;
     FFTW_REAL tre0_6_1;
     FFTW_REAL tim0_6_1;
     FFTW_REAL tre0_6_2;
     FFTW_REAL tim0_6_2;
     FFTW_REAL tre0_6_3;
     FFTW_REAL tim0_6_3;
     FFTW_REAL tre0_6_4;
     FFTW_REAL tim0_6_4;
     FFTW_REAL tre0_6_5;
     FFTW_REAL tim0_6_5;
     FFTW_REAL tre0_6_6;
     FFTW_REAL tim0_6_6;
     FFTW_REAL tre0_6_7;
     FFTW_REAL tim0_6_7;
     FFTW_REAL tre0_7_0;
     FFTW_REAL tim0_7_0;
     FFTW_REAL tre0_7_1;
     FFTW_REAL tim0_7_1;
     FFTW_REAL tre0_7_2;
     FFTW_REAL tim0_7_2;
     FFTW_REAL tre0_7_3;
     FFTW_REAL tim0_7_3;
     FFTW_REAL tre0_7_4;
     FFTW_REAL tim0_7_4;
     FFTW_REAL tre0_7_5;
     FFTW_REAL tim0_7_5;
     FFTW_REAL tre0_7_6;
     FFTW_REAL tim0_7_6;
     FFTW_REAL tre0_7_7;
     FFTW_REAL tim0_7_7;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[0]);
	       tim2_0_0 = c_im(in[0]);
	       tre2_1_0 = c_re(in[32 * istride]);
	       tim2_1_0 = c_im(in[32 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[8 * istride]);
	       tim2_0_0 = c_im(in[8 * istride]);
	       tre2_1_0 = c_re(in[40 * istride]);
	       tim2_1_0 = c_im(in[40 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[16 * istride]);
	       tim2_0_0 = c_im(in[16 * istride]);
	       tre2_1_0 = c_re(in[48 * istride]);
	       tim2_1_0 = c_im(in[48 * istride]);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[24 * istride]);
	       tim2_0_0 = c_im(in[24 * istride]);
	       tre2_1_0 = c_re(in[56 * istride]);
	       tim2_1_0 = c_im(in[56 * istride]);
	       tre1_0_3 = tre2_0_0 + tre2_1_0;
	       tim1_0_3 = tim2_0_0 + tim2_1_0;
	       tre1_1_3 = tre2_0_0 - tre2_1_0;
	       tim1_1_3 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       tre0_0_0 = tre2_0_0 + tre2_0_1;
	       tim0_0_0 = tim2_0_0 + tim2_0_1;
	       tre0_4_0 = tre2_0_0 - tre2_0_1;
	       tim0_4_0 = tim2_0_0 - tim2_0_1;
	       tre0_2_0 = tre2_1_0 - tim2_1_1;
	       tim0_2_0 = tim2_1_0 + tre2_1_1;
	       tre0_6_0 = tre2_1_0 + tim2_1_1;
	       tim0_6_0 = tim2_1_0 - tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 - tim1_1_2;
	       tim2_0_0 = tim1_1_0 + tre1_1_2;
	       tre2_1_0 = tre1_1_0 + tim1_1_2;
	       tim2_1_0 = tim1_1_0 - tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
		    tre2_0_1 = tre3_0_0 - tre3_1_0;
		    tim2_0_1 = tim3_0_0 + tim3_1_0;
		    tre2_1_1 = tre3_0_0 + tre3_1_0;
		    tim2_1_1 = tim3_0_0 - tim3_1_0;
	       }
	       tre0_1_0 = tre2_0_0 + tre2_0_1;
	       tim0_1_0 = tim2_0_0 + tim2_0_1;
	       tre0_5_0 = tre2_0_0 - tre2_0_1;
	       tim0_5_0 = tim2_0_0 - tim2_0_1;
	       tre0_3_0 = tre2_1_0 - tim2_1_1;
	       tim0_3_0 = tim2_1_0 + tre2_1_1;
	       tre0_7_0 = tre2_1_0 + tim2_1_1;
	       tim0_7_0 = tim2_1_0 - tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[istride]);
	       tim2_0_0 = c_im(in[istride]);
	       tre2_1_0 = c_re(in[33 * istride]);
	       tim2_1_0 = c_im(in[33 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[9 * istride]);
	       tim2_0_0 = c_im(in[9 * istride]);
	       tre2_1_0 = c_re(in[41 * istride]);
	       tim2_1_0 = c_im(in[41 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[17 * istride]);
	       tim2_0_0 = c_im(in[17 * istride]);
	       tre2_1_0 = c_re(in[49 * istride]);
	       tim2_1_0 = c_im(in[49 * istride]);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[25 * istride]);
	       tim2_0_0 = c_im(in[25 * istride]);
	       tre2_1_0 = c_re(in[57 * istride]);
	       tim2_1_0 = c_im(in[57 * istride]);
	       tre1_0_3 = tre2_0_0 + tre2_1_0;
	       tim1_0_3 = tim2_0_0 + tim2_1_0;
	       tre1_1_3 = tre2_0_0 - tre2_1_0;
	       tim1_1_3 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       tre0_0_1 = tre2_0_0 + tre2_0_1;
	       tim0_0_1 = tim2_0_0 + tim2_0_1;
	       tre0_4_1 = tre2_0_0 - tre2_0_1;
	       tim0_4_1 = tim2_0_0 - tim2_0_1;
	       tre0_2_1 = tre2_1_0 - tim2_1_1;
	       tim0_2_1 = tim2_1_0 + tre2_1_1;
	       tre0_6_1 = tre2_1_0 + tim2_1_1;
	       tim0_6_1 = tim2_1_0 - tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 - tim1_1_2;
	       tim2_0_0 = tim1_1_0 + tre1_1_2;
	       tre2_1_0 = tre1_1_0 + tim1_1_2;
	       tim2_1_0 = tim1_1_0 - tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
		    tre2_0_1 = tre3_0_0 - tre3_1_0;
		    tim2_0_1 = tim3_0_0 + tim3_1_0;
		    tre2_1_1 = tre3_0_0 + tre3_1_0;
		    tim2_1_1 = tim3_0_0 - tim3_1_0;
	       }
	       tre0_1_1 = tre2_0_0 + tre2_0_1;
	       tim0_1_1 = tim2_0_0 + tim2_0_1;
	       tre0_5_1 = tre2_0_0 - tre2_0_1;
	       tim0_5_1 = tim2_0_0 - tim2_0_1;
	       tre0_3_1 = tre2_1_0 - tim2_1_1;
	       tim0_3_1 = tim2_1_0 + tre2_1_1;
	       tre0_7_1 = tre2_1_0 + tim2_1_1;
	       tim0_7_1 = tim2_1_0 - tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[2 * istride]);
	       tim2_0_0 = c_im(in[2 * istride]);
	       tre2_1_0 = c_re(in[34 * istride]);
	       tim2_1_0 = c_im(in[34 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[10 * istride]);
	       tim2_0_0 = c_im(in[10 * istride]);
	       tre2_1_0 = c_re(in[42 * istride]);
	       tim2_1_0 = c_im(in[42 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[18 * istride]);
	       tim2_0_0 = c_im(in[18 * istride]);
	       tre2_1_0 = c_re(in[50 * istride]);
	       tim2_1_0 = c_im(in[50 * istride]);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[26 * istride]);
	       tim2_0_0 = c_im(in[26 * istride]);
	       tre2_1_0 = c_re(in[58 * istride]);
	       tim2_1_0 = c_im(in[58 * istride]);
	       tre1_0_3 = tre2_0_0 + tre2_1_0;
	       tim1_0_3 = tim2_0_0 + tim2_1_0;
	       tre1_1_3 = tre2_0_0 - tre2_1_0;
	       tim1_1_3 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       tre0_0_2 = tre2_0_0 + tre2_0_1;
	       tim0_0_2 = tim2_0_0 + tim2_0_1;
	       tre0_4_2 = tre2_0_0 - tre2_0_1;
	       tim0_4_2 = tim2_0_0 - tim2_0_1;
	       tre0_2_2 = tre2_1_0 - tim2_1_1;
	       tim0_2_2 = tim2_1_0 + tre2_1_1;
	       tre0_6_2 = tre2_1_0 + tim2_1_1;
	       tim0_6_2 = tim2_1_0 - tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 - tim1_1_2;
	       tim2_0_0 = tim1_1_0 + tre1_1_2;
	       tre2_1_0 = tre1_1_0 + tim1_1_2;
	       tim2_1_0 = tim1_1_0 - tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
		    tre2_0_1 = tre3_0_0 - tre3_1_0;
		    tim2_0_1 = tim3_0_0 + tim3_1_0;
		    tre2_1_1 = tre3_0_0 + tre3_1_0;
		    tim2_1_1 = tim3_0_0 - tim3_1_0;
	       }
	       tre0_1_2 = tre2_0_0 + tre2_0_1;
	       tim0_1_2 = tim2_0_0 + tim2_0_1;
	       tre0_5_2 = tre2_0_0 - tre2_0_1;
	       tim0_5_2 = tim2_0_0 - tim2_0_1;
	       tre0_3_2 = tre2_1_0 - tim2_1_1;
	       tim0_3_2 = tim2_1_0 + tre2_1_1;
	       tre0_7_2 = tre2_1_0 + tim2_1_1;
	       tim0_7_2 = tim2_1_0 - tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[3 * istride]);
	       tim2_0_0 = c_im(in[3 * istride]);
	       tre2_1_0 = c_re(in[35 * istride]);
	       tim2_1_0 = c_im(in[35 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[11 * istride]);
	       tim2_0_0 = c_im(in[11 * istride]);
	       tre2_1_0 = c_re(in[43 * istride]);
	       tim2_1_0 = c_im(in[43 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[19 * istride]);
	       tim2_0_0 = c_im(in[19 * istride]);
	       tre2_1_0 = c_re(in[51 * istride]);
	       tim2_1_0 = c_im(in[51 * istride]);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[27 * istride]);
	       tim2_0_0 = c_im(in[27 * istride]);
	       tre2_1_0 = c_re(in[59 * istride]);
	       tim2_1_0 = c_im(in[59 * istride]);
	       tre1_0_3 = tre2_0_0 + tre2_1_0;
	       tim1_0_3 = tim2_0_0 + tim2_1_0;
	       tre1_1_3 = tre2_0_0 - tre2_1_0;
	       tim1_1_3 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       tre0_0_3 = tre2_0_0 + tre2_0_1;
	       tim0_0_3 = tim2_0_0 + tim2_0_1;
	       tre0_4_3 = tre2_0_0 - tre2_0_1;
	       tim0_4_3 = tim2_0_0 - tim2_0_1;
	       tre0_2_3 = tre2_1_0 - tim2_1_1;
	       tim0_2_3 = tim2_1_0 + tre2_1_1;
	       tre0_6_3 = tre2_1_0 + tim2_1_1;
	       tim0_6_3 = tim2_1_0 - tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 - tim1_1_2;
	       tim2_0_0 = tim1_1_0 + tre1_1_2;
	       tre2_1_0 = tre1_1_0 + tim1_1_2;
	       tim2_1_0 = tim1_1_0 - tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
		    tre2_0_1 = tre3_0_0 - tre3_1_0;
		    tim2_0_1 = tim3_0_0 + tim3_1_0;
		    tre2_1_1 = tre3_0_0 + tre3_1_0;
		    tim2_1_1 = tim3_0_0 - tim3_1_0;
	       }
	       tre0_1_3 = tre2_0_0 + tre2_0_1;
	       tim0_1_3 = tim2_0_0 + tim2_0_1;
	       tre0_5_3 = tre2_0_0 - tre2_0_1;
	       tim0_5_3 = tim2_0_0 - tim2_0_1;
	       tre0_3_3 = tre2_1_0 - tim2_1_1;
	       tim0_3_3 = tim2_1_0 + tre2_1_1;
	       tre0_7_3 = tre2_1_0 + tim2_1_1;
	       tim0_7_3 = tim2_1_0 - tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[4 * istride]);
	       tim2_0_0 = c_im(in[4 * istride]);
	       tre2_1_0 = c_re(in[36 * istride]);
	       tim2_1_0 = c_im(in[36 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[12 * istride]);
	       tim2_0_0 = c_im(in[12 * istride]);
	       tre2_1_0 = c_re(in[44 * istride]);
	       tim2_1_0 = c_im(in[44 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[20 * istride]);
	       tim2_0_0 = c_im(in[20 * istride]);
	       tre2_1_0 = c_re(in[52 * istride]);
	       tim2_1_0 = c_im(in[52 * istride]);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[28 * istride]);
	       tim2_0_0 = c_im(in[28 * istride]);
	       tre2_1_0 = c_re(in[60 * istride]);
	       tim2_1_0 = c_im(in[60 * istride]);
	       tre1_0_3 = tre2_0_0 + tre2_1_0;
	       tim1_0_3 = tim2_0_0 + tim2_1_0;
	       tre1_1_3 = tre2_0_0 - tre2_1_0;
	       tim1_1_3 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       tre0_0_4 = tre2_0_0 + tre2_0_1;
	       tim0_0_4 = tim2_0_0 + tim2_0_1;
	       tre0_4_4 = tre2_0_0 - tre2_0_1;
	       tim0_4_4 = tim2_0_0 - tim2_0_1;
	       tre0_2_4 = tre2_1_0 - tim2_1_1;
	       tim0_2_4 = tim2_1_0 + tre2_1_1;
	       tre0_6_4 = tre2_1_0 + tim2_1_1;
	       tim0_6_4 = tim2_1_0 - tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 - tim1_1_2;
	       tim2_0_0 = tim1_1_0 + tre1_1_2;
	       tre2_1_0 = tre1_1_0 + tim1_1_2;
	       tim2_1_0 = tim1_1_0 - tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
		    tre2_0_1 = tre3_0_0 - tre3_1_0;
		    tim2_0_1 = tim3_0_0 + tim3_1_0;
		    tre2_1_1 = tre3_0_0 + tre3_1_0;
		    tim2_1_1 = tim3_0_0 - tim3_1_0;
	       }
	       tre0_1_4 = tre2_0_0 + tre2_0_1;
	       tim0_1_4 = tim2_0_0 + tim2_0_1;
	       tre0_5_4 = tre2_0_0 - tre2_0_1;
	       tim0_5_4 = tim2_0_0 - tim2_0_1;
	       tre0_3_4 = tre2_1_0 - tim2_1_1;
	       tim0_3_4 = tim2_1_0 + tre2_1_1;
	       tre0_7_4 = tre2_1_0 + tim2_1_1;
	       tim0_7_4 = tim2_1_0 - tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[5 * istride]);
	       tim2_0_0 = c_im(in[5 * istride]);
	       tre2_1_0 = c_re(in[37 * istride]);
	       tim2_1_0 = c_im(in[37 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[13 * istride]);
	       tim2_0_0 = c_im(in[13 * istride]);
	       tre2_1_0 = c_re(in[45 * istride]);
	       tim2_1_0 = c_im(in[45 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[21 * istride]);
	       tim2_0_0 = c_im(in[21 * istride]);
	       tre2_1_0 = c_re(in[53 * istride]);
	       tim2_1_0 = c_im(in[53 * istride]);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[29 * istride]);
	       tim2_0_0 = c_im(in[29 * istride]);
	       tre2_1_0 = c_re(in[61 * istride]);
	       tim2_1_0 = c_im(in[61 * istride]);
	       tre1_0_3 = tre2_0_0 + tre2_1_0;
	       tim1_0_3 = tim2_0_0 + tim2_1_0;
	       tre1_1_3 = tre2_0_0 - tre2_1_0;
	       tim1_1_3 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       tre0_0_5 = tre2_0_0 + tre2_0_1;
	       tim0_0_5 = tim2_0_0 + tim2_0_1;
	       tre0_4_5 = tre2_0_0 - tre2_0_1;
	       tim0_4_5 = tim2_0_0 - tim2_0_1;
	       tre0_2_5 = tre2_1_0 - tim2_1_1;
	       tim0_2_5 = tim2_1_0 + tre2_1_1;
	       tre0_6_5 = tre2_1_0 + tim2_1_1;
	       tim0_6_5 = tim2_1_0 - tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 - tim1_1_2;
	       tim2_0_0 = tim1_1_0 + tre1_1_2;
	       tre2_1_0 = tre1_1_0 + tim1_1_2;
	       tim2_1_0 = tim1_1_0 - tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
		    tre2_0_1 = tre3_0_0 - tre3_1_0;
		    tim2_0_1 = tim3_0_0 + tim3_1_0;
		    tre2_1_1 = tre3_0_0 + tre3_1_0;
		    tim2_1_1 = tim3_0_0 - tim3_1_0;
	       }
	       tre0_1_5 = tre2_0_0 + tre2_0_1;
	       tim0_1_5 = tim2_0_0 + tim2_0_1;
	       tre0_5_5 = tre2_0_0 - tre2_0_1;
	       tim0_5_5 = tim2_0_0 - tim2_0_1;
	       tre0_3_5 = tre2_1_0 - tim2_1_1;
	       tim0_3_5 = tim2_1_0 + tre2_1_1;
	       tre0_7_5 = tre2_1_0 + tim2_1_1;
	       tim0_7_5 = tim2_1_0 - tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[6 * istride]);
	       tim2_0_0 = c_im(in[6 * istride]);
	       tre2_1_0 = c_re(in[38 * istride]);
	       tim2_1_0 = c_im(in[38 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[14 * istride]);
	       tim2_0_0 = c_im(in[14 * istride]);
	       tre2_1_0 = c_re(in[46 * istride]);
	       tim2_1_0 = c_im(in[46 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[22 * istride]);
	       tim2_0_0 = c_im(in[22 * istride]);
	       tre2_1_0 = c_re(in[54 * istride]);
	       tim2_1_0 = c_im(in[54 * istride]);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[30 * istride]);
	       tim2_0_0 = c_im(in[30 * istride]);
	       tre2_1_0 = c_re(in[62 * istride]);
	       tim2_1_0 = c_im(in[62 * istride]);
	       tre1_0_3 = tre2_0_0 + tre2_1_0;
	       tim1_0_3 = tim2_0_0 + tim2_1_0;
	       tre1_1_3 = tre2_0_0 - tre2_1_0;
	       tim1_1_3 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       tre0_0_6 = tre2_0_0 + tre2_0_1;
	       tim0_0_6 = tim2_0_0 + tim2_0_1;
	       tre0_4_6 = tre2_0_0 - tre2_0_1;
	       tim0_4_6 = tim2_0_0 - tim2_0_1;
	       tre0_2_6 = tre2_1_0 - tim2_1_1;
	       tim0_2_6 = tim2_1_0 + tre2_1_1;
	       tre0_6_6 = tre2_1_0 + tim2_1_1;
	       tim0_6_6 = tim2_1_0 - tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 - tim1_1_2;
	       tim2_0_0 = tim1_1_0 + tre1_1_2;
	       tre2_1_0 = tre1_1_0 + tim1_1_2;
	       tim2_1_0 = tim1_1_0 - tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
		    tre2_0_1 = tre3_0_0 - tre3_1_0;
		    tim2_0_1 = tim3_0_0 + tim3_1_0;
		    tre2_1_1 = tre3_0_0 + tre3_1_0;
		    tim2_1_1 = tim3_0_0 - tim3_1_0;
	       }
	       tre0_1_6 = tre2_0_0 + tre2_0_1;
	       tim0_1_6 = tim2_0_0 + tim2_0_1;
	       tre0_5_6 = tre2_0_0 - tre2_0_1;
	       tim0_5_6 = tim2_0_0 - tim2_0_1;
	       tre0_3_6 = tre2_1_0 - tim2_1_1;
	       tim0_3_6 = tim2_1_0 + tre2_1_1;
	       tre0_7_6 = tre2_1_0 + tim2_1_1;
	       tim0_7_6 = tim2_1_0 - tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[7 * istride]);
	       tim2_0_0 = c_im(in[7 * istride]);
	       tre2_1_0 = c_re(in[39 * istride]);
	       tim2_1_0 = c_im(in[39 * istride]);
	       tre1_0_0 = tre2_0_0 + tre2_1_0;
	       tim1_0_0 = tim2_0_0 + tim2_1_0;
	       tre1_1_0 = tre2_0_0 - tre2_1_0;
	       tim1_1_0 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[15 * istride]);
	       tim2_0_0 = c_im(in[15 * istride]);
	       tre2_1_0 = c_re(in[47 * istride]);
	       tim2_1_0 = c_im(in[47 * istride]);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[23 * istride]);
	       tim2_0_0 = c_im(in[23 * istride]);
	       tre2_1_0 = c_re(in[55 * istride]);
	       tim2_1_0 = c_im(in[55 * istride]);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = c_re(in[31 * istride]);
	       tim2_0_0 = c_im(in[31 * istride]);
	       tre2_1_0 = c_re(in[63 * istride]);
	       tim2_1_0 = c_im(in[63 * istride]);
	       tre1_0_3 = tre2_0_0 + tre2_1_0;
	       tim1_0_3 = tim2_0_0 + tim2_1_0;
	       tre1_1_3 = tre2_0_0 - tre2_1_0;
	       tim1_1_3 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       tre0_0_7 = tre2_0_0 + tre2_0_1;
	       tim0_0_7 = tim2_0_0 + tim2_0_1;
	       tre0_4_7 = tre2_0_0 - tre2_0_1;
	       tim0_4_7 = tim2_0_0 - tim2_0_1;
	       tre0_2_7 = tre2_1_0 - tim2_1_1;
	       tim0_2_7 = tim2_1_0 + tre2_1_1;
	       tre0_6_7 = tre2_1_0 + tim2_1_1;
	       tim0_6_7 = tim2_1_0 - tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 - tim1_1_2;
	       tim2_0_0 = tim1_1_0 + tre1_1_2;
	       tre2_1_0 = tre1_1_0 + tim1_1_2;
	       tim2_1_0 = tim1_1_0 - tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
		    tre2_0_1 = tre3_0_0 - tre3_1_0;
		    tim2_0_1 = tim3_0_0 + tim3_1_0;
		    tre2_1_1 = tre3_0_0 + tre3_1_0;
		    tim2_1_1 = tim3_0_0 - tim3_1_0;
	       }
	       tre0_1_7 = tre2_0_0 + tre2_0_1;
	       tim0_1_7 = tim2_0_0 + tim2_0_1;
	       tre0_5_7 = tre2_0_0 - tre2_0_1;
	       tim0_5_7 = tim2_0_0 - tim2_0_1;
	       tre0_3_7 = tre2_1_0 - tim2_1_1;
	       tim0_3_7 = tim2_1_0 + tre2_1_1;
	       tre0_7_7 = tre2_1_0 + tim2_1_1;
	       tim0_7_7 = tim2_1_0 - tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  tre1_0_0 = tre0_0_0 + tre0_0_4;
	  tim1_0_0 = tim0_0_0 + tim0_0_4;
	  tre1_1_0 = tre0_0_0 - tre0_0_4;
	  tim1_1_0 = tim0_0_0 - tim0_0_4;
	  tre1_0_1 = tre0_0_1 + tre0_0_5;
	  tim1_0_1 = tim0_0_1 + tim0_0_5;
	  tre1_1_1 = tre0_0_1 - tre0_0_5;
	  tim1_1_1 = tim0_0_1 - tim0_0_5;
	  tre1_0_2 = tre0_0_2 + tre0_0_6;
	  tim1_0_2 = tim0_0_2 + tim0_0_6;
	  tre1_1_2 = tre0_0_2 - tre0_0_6;
	  tim1_1_2 = tim0_0_2 - tim0_0_6;
	  tre1_0_3 = tre0_0_3 + tre0_0_7;
	  tim1_0_3 = tim0_0_3 + tim0_0_7;
	  tre1_1_3 = tre0_0_3 - tre0_0_7;
	  tim1_1_3 = tim0_0_3 - tim0_0_7;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       c_re(out[0]) = tre2_0_0 + tre2_0_1;
	       c_im(out[0]) = tim2_0_0 + tim2_0_1;
	       c_re(out[32 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[32 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[16 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[16 * ostride]) = tim2_1_0 + tre2_1_1;
	       c_re(out[48 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[48 * ostride]) = tim2_1_0 - tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 - tim1_1_2;
	       tim2_0_0 = tim1_1_0 + tre1_1_2;
	       tre2_1_0 = tre1_1_0 + tim1_1_2;
	       tim2_1_0 = tim1_1_0 - tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
		    tre2_0_1 = tre3_0_0 - tre3_1_0;
		    tim2_0_1 = tim3_0_0 + tim3_1_0;
		    tre2_1_1 = tre3_0_0 + tre3_1_0;
		    tim2_1_1 = tim3_0_0 - tim3_1_0;
	       }
	       c_re(out[8 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[8 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[40 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[40 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[24 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[24 * ostride]) = tim2_1_0 + tre2_1_1;
	       c_re(out[56 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[56 * ostride]) = tim2_1_0 - tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_1_4) - (((FFTW_REAL) FFTW_K382683432) * tim0_1_4);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_1_4) + (((FFTW_REAL) FFTW_K382683432) * tre0_1_4);
	       tre1_0_0 = tre0_1_0 + tre2_1_0;
	       tim1_0_0 = tim0_1_0 + tim2_1_0;
	       tre1_1_0 = tre0_1_0 - tre2_1_0;
	       tim1_1_0 = tim0_1_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K995184726) * tre0_1_1) - (((FFTW_REAL) FFTW_K098017140) * tim0_1_1);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K995184726) * tim0_1_1) + (((FFTW_REAL) FFTW_K098017140) * tre0_1_1);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K881921264) * tre0_1_5) - (((FFTW_REAL) FFTW_K471396736) * tim0_1_5);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K881921264) * tim0_1_5) + (((FFTW_REAL) FFTW_K471396736) * tre0_1_5);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tre0_1_2) - (((FFTW_REAL) FFTW_K195090322) * tim0_1_2);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tim0_1_2) + (((FFTW_REAL) FFTW_K195090322) * tre0_1_2);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_1_6) - (((FFTW_REAL) FFTW_K555570233) * tim0_1_6);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_1_6) + (((FFTW_REAL) FFTW_K555570233) * tre0_1_6);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K956940335) * tre0_1_3) - (((FFTW_REAL) FFTW_K290284677) * tim0_1_3);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K956940335) * tim0_1_3) + (((FFTW_REAL) FFTW_K290284677) * tre0_1_3);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K773010453) * tre0_1_7) - (((FFTW_REAL) FFTW_K634393284) * tim0_1_7);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K773010453) * tim0_1_7) + (((FFTW_REAL) FFTW_K634393284) * tre0_1_7);
	       tre1_0_3 = tre2_0_0 + tre2_1_0;
	       tim1_0_3 = tim2_0_0 + tim2_1_0;
	       tre1_1_3 = tre2_0_0 - tre2_1_0;
	       tim1_1_3 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       c_re(out[ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[33 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[33 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[17 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[17 * ostride]) = tim2_1_0 + tre2_1_1;
	       c_re(out[49 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[49 * ostride]) = tim2_1_0 - tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 - tim1_1_2;
	       tim2_0_0 = tim1_1_0 + tre1_1_2;
	       tre2_1_0 = tre1_1_0 + tim1_1_2;
	       tim2_1_0 = tim1_1_0 - tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
		    tre2_0_1 = tre3_0_0 - tre3_1_0;
		    tim2_0_1 = tim3_0_0 + tim3_1_0;
		    tre2_1_1 = tre3_0_0 + tre3_1_0;
		    tim2_1_1 = tim3_0_0 - tim3_1_0;
	       }
	       c_re(out[9 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[9 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[41 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[41 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[25 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[25 * ostride]) = tim2_1_0 + tre2_1_1;
	       c_re(out[57 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[57 * ostride]) = tim2_1_0 - tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_2_4 - tim0_2_4);
	       tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_2_4 + tre0_2_4);
	       tre1_0_0 = tre0_2_0 + tre2_1_0;
	       tim1_0_0 = tim0_2_0 + tim2_1_0;
	       tre1_1_0 = tre0_2_0 - tre2_1_0;
	       tim1_1_0 = tim0_2_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tre0_2_1) - (((FFTW_REAL) FFTW_K195090322) * tim0_2_1);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tim0_2_1) + (((FFTW_REAL) FFTW_K195090322) * tre0_2_1);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K555570233) * tre0_2_5) - (((FFTW_REAL) FFTW_K831469612) * tim0_2_5);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K555570233) * tim0_2_5) + (((FFTW_REAL) FFTW_K831469612) * tre0_2_5);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_2_2) - (((FFTW_REAL) FFTW_K382683432) * tim0_2_2);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_2_2) + (((FFTW_REAL) FFTW_K382683432) * tre0_2_2);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_2_6) - (((FFTW_REAL) FFTW_K923879532) * tim0_2_6);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_2_6) + (((FFTW_REAL) FFTW_K923879532) * tre0_2_6);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_2_3) - (((FFTW_REAL) FFTW_K555570233) * tim0_2_3);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_2_3) + (((FFTW_REAL) FFTW_K555570233) * tre0_2_3);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tre0_2_7) - (((FFTW_REAL) FFTW_K980785280) * tim0_2_7);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tim0_2_7) + (((FFTW_REAL) FFTW_K980785280) * tre0_2_7);
	       tre1_0_3 = tre2_0_0 + tre2_1_0;
	       tim1_0_3 = tim2_0_0 + tim2_1_0;
	       tre1_1_3 = tre2_0_0 - tre2_1_0;
	       tim1_1_3 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       c_re(out[2 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[2 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[34 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[34 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[18 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[18 * ostride]) = tim2_1_0 + tre2_1_1;
	       c_re(out[50 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[50 * ostride]) = tim2_1_0 - tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 - tim1_1_2;
	       tim2_0_0 = tim1_1_0 + tre1_1_2;
	       tre2_1_0 = tre1_1_0 + tim1_1_2;
	       tim2_1_0 = tim1_1_0 - tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
		    tre2_0_1 = tre3_0_0 - tre3_1_0;
		    tim2_0_1 = tim3_0_0 + tim3_1_0;
		    tre2_1_1 = tre3_0_0 + tre3_1_0;
		    tim2_1_1 = tim3_0_0 - tim3_1_0;
	       }
	       c_re(out[10 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[10 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[42 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[42 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[26 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[26 * ostride]) = tim2_1_0 + tre2_1_1;
	       c_re(out[58 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[58 * ostride]) = tim2_1_0 - tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_3_4) - (((FFTW_REAL) FFTW_K923879532) * tim0_3_4);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_3_4) + (((FFTW_REAL) FFTW_K923879532) * tre0_3_4);
	       tre1_0_0 = tre0_3_0 + tre2_1_0;
	       tim1_0_0 = tim0_3_0 + tim2_1_0;
	       tre1_1_0 = tre0_3_0 - tre2_1_0;
	       tim1_1_0 = tim0_3_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K956940335) * tre0_3_1) - (((FFTW_REAL) FFTW_K290284677) * tim0_3_1);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K956940335) * tim0_3_1) + (((FFTW_REAL) FFTW_K290284677) * tre0_3_1);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K098017140) * tre0_3_5) - (((FFTW_REAL) FFTW_K995184726) * tim0_3_5);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K098017140) * tim0_3_5) + (((FFTW_REAL) FFTW_K995184726) * tre0_3_5);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_3_2) - (((FFTW_REAL) FFTW_K555570233) * tim0_3_2);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_3_2) + (((FFTW_REAL) FFTW_K555570233) * tre0_3_2);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tre0_3_6) + (((FFTW_REAL) FFTW_K980785280) * tim0_3_6);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K980785280) * tre0_3_6) - (((FFTW_REAL) FFTW_K195090322) * tim0_3_6);
	       tre1_0_2 = tre2_0_0 - tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 + tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K634393284) * tre0_3_3) - (((FFTW_REAL) FFTW_K773010453) * tim0_3_3);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K634393284) * tim0_3_3) + (((FFTW_REAL) FFTW_K773010453) * tre0_3_3);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K471396736) * tre0_3_7) + (((FFTW_REAL) FFTW_K881921264) * tim0_3_7);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K881921264) * tre0_3_7) - (((FFTW_REAL) FFTW_K471396736) * tim0_3_7);
	       tre1_0_3 = tre2_0_0 - tre2_1_0;
	       tim1_0_3 = tim2_0_0 + tim2_1_0;
	       tre1_1_3 = tre2_0_0 + tre2_1_0;
	       tim1_1_3 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       c_re(out[3 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[3 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[35 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[35 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[19 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[19 * ostride]) = tim2_1_0 + tre2_1_1;
	       c_re(out[51 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[51 * ostride]) = tim2_1_0 - tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 - tim1_1_2;
	       tim2_0_0 = tim1_1_0 + tre1_1_2;
	       tre2_1_0 = tre1_1_0 + tim1_1_2;
	       tim2_1_0 = tim1_1_0 - tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
		    tre2_0_1 = tre3_0_0 - tre3_1_0;
		    tim2_0_1 = tim3_0_0 + tim3_1_0;
		    tre2_1_1 = tre3_0_0 + tre3_1_0;
		    tim2_1_1 = tim3_0_0 - tim3_1_0;
	       }
	       c_re(out[11 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[11 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[43 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[43 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[27 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[27 * ostride]) = tim2_1_0 + tre2_1_1;
	       c_re(out[59 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[59 * ostride]) = tim2_1_0 - tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  tre1_0_0 = tre0_4_0 - tim0_4_4;
	  tim1_0_0 = tim0_4_0 + tre0_4_4;
	  tre1_1_0 = tre0_4_0 + tim0_4_4;
	  tim1_1_0 = tim0_4_0 - tre0_4_4;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_4_1) - (((FFTW_REAL) FFTW_K382683432) * tim0_4_1);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_4_1) + (((FFTW_REAL) FFTW_K382683432) * tre0_4_1);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_4_5) + (((FFTW_REAL) FFTW_K923879532) * tim0_4_5);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_4_5) - (((FFTW_REAL) FFTW_K382683432) * tim0_4_5);
	       tre1_0_1 = tre2_0_0 - tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 + tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_4_2 - tim0_4_2);
	       tim2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_4_2 + tre0_4_2);
	       tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_4_6 + tim0_4_6);
	       tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_4_6 - tim0_4_6);
	       tre1_0_2 = tre2_0_0 - tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 + tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_4_3) - (((FFTW_REAL) FFTW_K923879532) * tim0_4_3);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_4_3) + (((FFTW_REAL) FFTW_K923879532) * tre0_4_3);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_4_7) + (((FFTW_REAL) FFTW_K382683432) * tim0_4_7);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_4_7) - (((FFTW_REAL) FFTW_K923879532) * tim0_4_7);
	       tre1_0_3 = tre2_0_0 - tre2_1_0;
	       tim1_0_3 = tim2_0_0 + tim2_1_0;
	       tre1_1_3 = tre2_0_0 + tre2_1_0;
	       tim1_1_3 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       c_re(out[4 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[4 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[36 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[36 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[20 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[20 * ostride]) = tim2_1_0 + tre2_1_1;
	       c_re(out[52 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[52 * ostride]) = tim2_1_0 - tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 - tim1_1_2;
	       tim2_0_0 = tim1_1_0 + tre1_1_2;
	       tre2_1_0 = tre1_1_0 + tim1_1_2;
	       tim2_1_0 = tim1_1_0 - tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
		    tre2_0_1 = tre3_0_0 - tre3_1_0;
		    tim2_0_1 = tim3_0_0 + tim3_1_0;
		    tre2_1_1 = tre3_0_0 + tre3_1_0;
		    tim2_1_1 = tim3_0_0 - tim3_1_0;
	       }
	       c_re(out[12 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[12 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[44 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[44 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[28 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[28 * ostride]) = tim2_1_0 + tre2_1_1;
	       c_re(out[60 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[60 * ostride]) = tim2_1_0 - tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_5_4) + (((FFTW_REAL) FFTW_K923879532) * tim0_5_4);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_5_4) - (((FFTW_REAL) FFTW_K382683432) * tim0_5_4);
	       tre1_0_0 = tre0_5_0 - tre2_1_0;
	       tim1_0_0 = tim0_5_0 + tim2_1_0;
	       tre1_1_0 = tre0_5_0 + tre2_1_0;
	       tim1_1_0 = tim0_5_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K881921264) * tre0_5_1) - (((FFTW_REAL) FFTW_K471396736) * tim0_5_1);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K881921264) * tim0_5_1) + (((FFTW_REAL) FFTW_K471396736) * tre0_5_1);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K773010453) * tre0_5_5) + (((FFTW_REAL) FFTW_K634393284) * tim0_5_5);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K634393284) * tre0_5_5) - (((FFTW_REAL) FFTW_K773010453) * tim0_5_5);
	       tre1_0_1 = tre2_0_0 - tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 + tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K555570233) * tre0_5_2) - (((FFTW_REAL) FFTW_K831469612) * tim0_5_2);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K555570233) * tim0_5_2) + (((FFTW_REAL) FFTW_K831469612) * tre0_5_2);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K980785280) * tre0_5_6) + (((FFTW_REAL) FFTW_K195090322) * tim0_5_6);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tre0_5_6) - (((FFTW_REAL) FFTW_K980785280) * tim0_5_6);
	       tre1_0_2 = tre2_0_0 - tre2_1_0;
	       tim1_0_2 = tim2_0_0 + tim2_1_0;
	       tre1_1_2 = tre2_0_0 + tre2_1_0;
	       tim1_1_2 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K098017140) * tre0_5_3) - (((FFTW_REAL) FFTW_K995184726) * tim0_5_3);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K098017140) * tim0_5_3) + (((FFTW_REAL) FFTW_K995184726) * tre0_5_3);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K290284677) * tim0_5_7) - (((FFTW_REAL) FFTW_K956940335) * tre0_5_7);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K956940335) * tim0_5_7) + (((FFTW_REAL) FFTW_K290284677) * tre0_5_7);
	       tre1_0_3 = tre2_0_0 + tre2_1_0;
	       tim1_0_3 = tim2_0_0 - tim2_1_0;
	       tre1_1_3 = tre2_0_0 - tre2_1_0;
	       tim1_1_3 = tim2_0_0 + tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       c_re(out[5 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[5 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[37 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[37 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[21 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[21 * ostride]) = tim2_1_0 + tre2_1_1;
	       c_re(out[53 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[53 * ostride]) = tim2_1_0 - tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 - tim1_1_2;
	       tim2_0_0 = tim1_1_0 + tre1_1_2;
	       tre2_1_0 = tre1_1_0 + tim1_1_2;
	       tim2_1_0 = tim1_1_0 - tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
		    tre2_0_1 = tre3_0_0 - tre3_1_0;
		    tim2_0_1 = tim3_0_0 + tim3_1_0;
		    tre2_1_1 = tre3_0_0 + tre3_1_0;
		    tim2_1_1 = tim3_0_0 - tim3_1_0;
	       }
	       c_re(out[13 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[13 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[45 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[45 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[29 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[29 * ostride]) = tim2_1_0 + tre2_1_1;
	       c_re(out[61 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[61 * ostride]) = tim2_1_0 - tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_6_4 + tim0_6_4);
	       tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_6_4 - tim0_6_4);
	       tre1_0_0 = tre0_6_0 - tre2_1_0;
	       tim1_0_0 = tim0_6_0 + tim2_1_0;
	       tre1_1_0 = tre0_6_0 + tre2_1_0;
	       tim1_1_0 = tim0_6_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_6_1) - (((FFTW_REAL) FFTW_K555570233) * tim0_6_1);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_6_1) + (((FFTW_REAL) FFTW_K555570233) * tre0_6_1);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K980785280) * tre0_6_5) + (((FFTW_REAL) FFTW_K195090322) * tim0_6_5);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tre0_6_5) - (((FFTW_REAL) FFTW_K980785280) * tim0_6_5);
	       tre1_0_1 = tre2_0_0 - tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 + tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_6_2) - (((FFTW_REAL) FFTW_K923879532) * tim0_6_2);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_6_2) + (((FFTW_REAL) FFTW_K923879532) * tre0_6_2);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_6_6) - (((FFTW_REAL) FFTW_K923879532) * tre0_6_6);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_6_6) + (((FFTW_REAL) FFTW_K382683432) * tre0_6_6);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 - tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 + tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K195090322) * tre0_6_3) + (((FFTW_REAL) FFTW_K980785280) * tim0_6_3);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tre0_6_3) - (((FFTW_REAL) FFTW_K195090322) * tim0_6_3);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_6_7) - (((FFTW_REAL) FFTW_K555570233) * tre0_6_7);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K555570233) * tim0_6_7) + (((FFTW_REAL) FFTW_K831469612) * tre0_6_7);
	       tre1_0_3 = tre2_1_0 - tre2_0_0;
	       tim1_0_3 = tim2_0_0 - tim2_1_0;
	       tre1_1_3 = (-(tre2_0_0 + tre2_1_0));
	       tim1_1_3 = tim2_0_0 + tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       c_re(out[6 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[6 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[38 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[38 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[22 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[22 * ostride]) = tim2_1_0 + tre2_1_1;
	       c_re(out[54 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[54 * ostride]) = tim2_1_0 - tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 - tim1_1_2;
	       tim2_0_0 = tim1_1_0 + tre1_1_2;
	       tre2_1_0 = tre1_1_0 + tim1_1_2;
	       tim2_1_0 = tim1_1_0 - tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
		    tre2_0_1 = tre3_0_0 - tre3_1_0;
		    tim2_0_1 = tim3_0_0 + tim3_1_0;
		    tre2_1_1 = tre3_0_0 + tre3_1_0;
		    tim2_1_1 = tim3_0_0 - tim3_1_0;
	       }
	       c_re(out[14 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[14 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[46 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[46 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[30 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[30 * ostride]) = tim2_1_0 + tre2_1_1;
	       c_re(out[62 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[62 * ostride]) = tim2_1_0 - tre2_1_1;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_0_2;
	  FFTW_REAL tim1_0_2;
	  FFTW_REAL tre1_0_3;
	  FFTW_REAL tim1_0_3;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  FFTW_REAL tre1_1_2;
	  FFTW_REAL tim1_1_2;
	  FFTW_REAL tre1_1_3;
	  FFTW_REAL tim1_1_3;
	  {
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_7_4) + (((FFTW_REAL) FFTW_K382683432) * tim0_7_4);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_7_4) - (((FFTW_REAL) FFTW_K923879532) * tim0_7_4);
	       tre1_0_0 = tre0_7_0 - tre2_1_0;
	       tim1_0_0 = tim0_7_0 + tim2_1_0;
	       tre1_1_0 = tre0_7_0 + tre2_1_0;
	       tim1_1_0 = tim0_7_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K773010453) * tre0_7_1) - (((FFTW_REAL) FFTW_K634393284) * tim0_7_1);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K773010453) * tim0_7_1) + (((FFTW_REAL) FFTW_K634393284) * tre0_7_1);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K290284677) * tim0_7_5) - (((FFTW_REAL) FFTW_K956940335) * tre0_7_5);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K956940335) * tim0_7_5) + (((FFTW_REAL) FFTW_K290284677) * tre0_7_5);
	       tre1_0_1 = tre2_0_0 + tre2_1_0;
	       tim1_0_1 = tim2_0_0 - tim2_1_0;
	       tre1_1_1 = tre2_0_0 - tre2_1_0;
	       tim1_1_1 = tim2_0_0 + tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K195090322) * tre0_7_2) - (((FFTW_REAL) FFTW_K980785280) * tim0_7_2);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K195090322) * tim0_7_2) + (((FFTW_REAL) FFTW_K980785280) * tre0_7_2);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_7_6) - (((FFTW_REAL) FFTW_K555570233) * tre0_7_6);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K555570233) * tim0_7_6) + (((FFTW_REAL) FFTW_K831469612) * tre0_7_6);
	       tre1_0_2 = tre2_0_0 + tre2_1_0;
	       tim1_0_2 = tim2_0_0 - tim2_1_0;
	       tre1_1_2 = tre2_0_0 - tre2_1_0;
	       tim1_1_2 = tim2_0_0 + tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = (((FFTW_REAL) FFTW_K471396736) * tre0_7_3) + (((FFTW_REAL) FFTW_K881921264) * tim0_7_3);
	       tim2_0_0 = (((FFTW_REAL) FFTW_K881921264) * tre0_7_3) - (((FFTW_REAL) FFTW_K471396736) * tim0_7_3);
	       tre2_1_0 = (((FFTW_REAL) FFTW_K098017140) * tre0_7_7) + (((FFTW_REAL) FFTW_K995184726) * tim0_7_7);
	       tim2_1_0 = (((FFTW_REAL) FFTW_K098017140) * tim0_7_7) - (((FFTW_REAL) FFTW_K995184726) * tre0_7_7);
	       tre1_0_3 = tre2_1_0 - tre2_0_0;
	       tim1_0_3 = tim2_0_0 + tim2_1_0;
	       tre1_1_3 = (-(tre2_0_0 + tre2_1_0));
	       tim1_1_3 = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_0_0 + tre1_0_2;
	       tim2_0_0 = tim1_0_0 + tim1_0_2;
	       tre2_1_0 = tre1_0_0 - tre1_0_2;
	       tim2_1_0 = tim1_0_0 - tim1_0_2;
	       tre2_0_1 = tre1_0_1 + tre1_0_3;
	       tim2_0_1 = tim1_0_1 + tim1_0_3;
	       tre2_1_1 = tre1_0_1 - tre1_0_3;
	       tim2_1_1 = tim1_0_1 - tim1_0_3;
	       c_re(out[7 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[7 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[39 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[39 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[23 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[23 * ostride]) = tim2_1_0 + tre2_1_1;
	       c_re(out[55 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[55 * ostride]) = tim2_1_0 - tre2_1_1;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_0_1;
	       FFTW_REAL tim2_0_1;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       FFTW_REAL tre2_1_1;
	       FFTW_REAL tim2_1_1;
	       tre2_0_0 = tre1_1_0 - tim1_1_2;
	       tim2_0_0 = tim1_1_0 + tre1_1_2;
	       tre2_1_0 = tre1_1_0 + tim1_1_2;
	       tim2_1_0 = tim1_1_0 - tre1_1_2;
	       {
		    FFTW_REAL tre3_0_0;
		    FFTW_REAL tim3_0_0;
		    FFTW_REAL tre3_1_0;
		    FFTW_REAL tim3_1_0;
		    tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
		    tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
		    tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
		    tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
		    tre2_0_1 = tre3_0_0 - tre3_1_0;
		    tim2_0_1 = tim3_0_0 + tim3_1_0;
		    tre2_1_1 = tre3_0_0 + tre3_1_0;
		    tim2_1_1 = tim3_0_0 - tim3_1_0;
	       }
	       c_re(out[15 * ostride]) = tre2_0_0 + tre2_0_1;
	       c_im(out[15 * ostride]) = tim2_0_0 + tim2_0_1;
	       c_re(out[47 * ostride]) = tre2_0_0 - tre2_0_1;
	       c_im(out[47 * ostride]) = tim2_0_0 - tim2_0_1;
	       c_re(out[31 * ostride]) = tre2_1_0 - tim2_1_1;
	       c_im(out[31 * ostride]) = tim2_1_0 + tre2_1_1;
	       c_re(out[63 * ostride]) = tre2_1_0 + tim2_1_1;
	       c_im(out[63 * ostride]) = tim2_1_0 - tre2_1_1;
	  }
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 90 FP additions and 36 FP multiplications */

void fftwi_no_twiddle_7(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_2_0;
     FFTW_REAL tim0_2_0;
     FFTW_REAL tre0_3_0;
     FFTW_REAL tim0_3_0;
     FFTW_REAL tre0_4_0;
     FFTW_REAL tim0_4_0;
     FFTW_REAL tre0_5_0;
     FFTW_REAL tim0_5_0;
     FFTW_REAL tre0_6_0;
     FFTW_REAL tim0_6_0;
     tre0_0_0 = c_re(in[0]);
     tim0_0_0 = c_im(in[0]);
     tre0_1_0 = c_re(in[istride]);
     tim0_1_0 = c_im(in[istride]);
     tre0_2_0 = c_re(in[2 * istride]);
     tim0_2_0 = c_im(in[2 * istride]);
     tre0_3_0 = c_re(in[3 * istride]);
     tim0_3_0 = c_im(in[3 * istride]);
     tre0_4_0 = c_re(in[4 * istride]);
     tim0_4_0 = c_im(in[4 * istride]);
     tre0_5_0 = c_re(in[5 * istride]);
     tim0_5_0 = c_im(in[5 * istride]);
     tre0_6_0 = c_re(in[6 * istride]);
     tim0_6_0 = c_im(in[6 * istride]);
     c_re(out[0]) = tre0_0_0 + tre0_1_0 + tre0_2_0 + tre0_3_0 + tre0_4_0 + tre0_5_0 + tre0_6_0;
     c_im(out[0]) = tim0_0_0 + tim0_1_0 + tim0_2_0 + tim0_3_0 + tim0_4_0 + tim0_5_0 + tim0_6_0;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tre0_1_0 + tre0_6_0)) - (((FFTW_REAL) FFTW_K900968867) * (tre0_3_0 + tre0_4_0)) - (((FFTW_REAL) FFTW_K222520933) * (tre0_2_0 + tre0_5_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K781831482) * (tim0_6_0 - tim0_1_0)) + (((FFTW_REAL) FFTW_K974927912) * (tim0_5_0 - tim0_2_0)) + (((FFTW_REAL) FFTW_K433883739) * (tim0_4_0 - tim0_3_0));
	  c_re(out[ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[6 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tim0_1_0 + tim0_6_0)) - (((FFTW_REAL) FFTW_K900968867) * (tim0_3_0 + tim0_4_0)) - (((FFTW_REAL) FFTW_K222520933) * (tim0_2_0 + tim0_5_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K781831482) * (tre0_1_0 - tre0_6_0)) + (((FFTW_REAL) FFTW_K974927912) * (tre0_2_0 - tre0_5_0)) + (((FFTW_REAL) FFTW_K433883739) * (tre0_3_0 - tre0_4_0));
	  c_im(out[ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[6 * ostride]) = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tre0_3_0 + tre0_4_0)) - (((FFTW_REAL) FFTW_K900968867) * (tre0_2_0 + tre0_5_0)) - (((FFTW_REAL) FFTW_K222520933) * (tre0_1_0 + tre0_6_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K974927912) * (tim0_6_0 - tim0_1_0)) + (((FFTW_REAL) FFTW_K433883739) * (tim0_2_0 - tim0_5_0)) + (((FFTW_REAL) FFTW_K781831482) * (tim0_3_0 - tim0_4_0));
	  c_re(out[2 * ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[5 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tim0_3_0 + tim0_4_0)) - (((FFTW_REAL) FFTW_K900968867) * (tim0_2_0 + tim0_5_0)) - (((FFTW_REAL) FFTW_K222520933) * (tim0_1_0 + tim0_6_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K974927912) * (tre0_1_0 - tre0_6_0)) + (((FFTW_REAL) FFTW_K433883739) * (tre0_5_0 - tre0_2_0)) + (((FFTW_REAL) FFTW_K781831482) * (tre0_4_0 - tre0_3_0));
	  c_im(out[2 * ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[5 * ostride]) = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tre1_1_0;
	  tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tre0_2_0 + tre0_5_0)) - (((FFTW_REAL) FFTW_K222520933) * (tre0_3_0 + tre0_4_0)) - (((FFTW_REAL) FFTW_K900968867) * (tre0_1_0 + tre0_6_0));
	  tre1_1_0 = (((FFTW_REAL) FFTW_K433883739) * (tim0_6_0 - tim0_1_0)) + (((FFTW_REAL) FFTW_K781831482) * (tim0_2_0 - tim0_5_0)) + (((FFTW_REAL) FFTW_K974927912) * (tim0_4_0 - tim0_3_0));
	  c_re(out[3 * ostride]) = tre1_0_0 + tre1_1_0;
	  c_re(out[4 * ostride]) = tre1_0_0 - tre1_1_0;
     }
     {
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tim1_1_0;
	  tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tim0_2_0 + tim0_5_0)) - (((FFTW_REAL) FFTW_K222520933) * (tim0_3_0 + tim0_4_0)) - (((FFTW_REAL) FFTW_K900968867) * (tim0_1_0 + tim0_6_0));
	  tim1_1_0 = (((FFTW_REAL) FFTW_K433883739) * (tre0_1_0 - tre0_6_0)) + (((FFTW_REAL) FFTW_K781831482) * (tre0_5_0 - tre0_2_0)) + (((FFTW_REAL) FFTW_K974927912) * (tre0_3_0 - tre0_4_0));
	  c_im(out[3 * ostride]) = tim1_0_0 + tim1_1_0;
	  c_im(out[4 * ostride]) = tim1_0_0 - tim1_1_0;
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 52 FP additions and 4 FP multiplications */

void fftwi_no_twiddle_8(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_0_1;
     FFTW_REAL tim0_0_1;
     FFTW_REAL tre0_0_2;
     FFTW_REAL tim0_0_2;
     FFTW_REAL tre0_0_3;
     FFTW_REAL tim0_0_3;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_1_1;
     FFTW_REAL tim0_1_1;
     FFTW_REAL tre0_1_2;
     FFTW_REAL tim0_1_2;
     FFTW_REAL tre0_1_3;
     FFTW_REAL tim0_1_3;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[0]);
	  tim1_0_0 = c_im(in[0]);
	  tre1_1_0 = c_re(in[4 * istride]);
	  tim1_1_0 = c_im(in[4 * istride]);
	  tre0_0_0 = tre1_0_0 + tre1_1_0;
	  tim0_0_0 = tim1_0_0 + tim1_1_0;
	  tre0_1_0 = tre1_0_0 - tre1_1_0;
	  tim0_1_0 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[istride]);
	  tim1_0_0 = c_im(in[istride]);
	  tre1_1_0 = c_re(in[5 * istride]);
	  tim1_1_0 = c_im(in[5 * istride]);
	  tre0_0_1 = tre1_0_0 + tre1_1_0;
	  tim0_0_1 = tim1_0_0 + tim1_1_0;
	  tre0_1_1 = tre1_0_0 - tre1_1_0;
	  tim0_1_1 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[2 * istride]);
	  tim1_0_0 = c_im(in[2 * istride]);
	  tre1_1_0 = c_re(in[6 * istride]);
	  tim1_1_0 = c_im(in[6 * istride]);
	  tre0_0_2 = tre1_0_0 + tre1_1_0;
	  tim0_0_2 = tim1_0_0 + tim1_1_0;
	  tre0_1_2 = tre1_0_0 - tre1_1_0;
	  tim0_1_2 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  tre1_0_0 = c_re(in[3 * istride]);
	  tim1_0_0 = c_im(in[3 * istride]);
	  tre1_1_0 = c_re(in[7 * istride]);
	  tim1_1_0 = c_im(in[7 * istride]);
	  tre0_0_3 = tre1_0_0 + tre1_1_0;
	  tim0_0_3 = tim1_0_0 + tim1_1_0;
	  tre0_1_3 = tre1_0_0 - tre1_1_0;
	  tim0_1_3 = tim1_0_0 - tim1_1_0;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  tre1_0_0 = tre0_0_0 + tre0_0_2;
	  tim1_0_0 = tim0_0_0 + tim0_0_2;
	  tre1_1_0 = tre0_0_0 - tre0_0_2;
	  tim1_1_0 = tim0_0_0 - tim0_0_2;
	  tre1_0_1 = tre0_0_1 + tre0_0_3;
	  tim1_0_1 = tim0_0_1 + tim0_0_3;
	  tre1_1_1 = tre0_0_1 - tre0_0_3;
	  tim1_1_1 = tim0_0_1 - tim0_0_3;
	  c_re(out[0]) = tre1_0_0 + tre1_0_1;
	  c_im(out[0]) = tim1_0_0 + tim1_0_1;
	  c_re(out[4 * ostride]) = tre1_0_0 - tre1_0_1;
	  c_im(out[4 * ostride]) = tim1_0_0 - tim1_0_1;
	  c_re(out[2 * ostride]) = tre1_1_0 - tim1_1_1;
	  c_im(out[2 * ostride]) = tim1_1_0 + tre1_1_1;
	  c_re(out[6 * ostride]) = tre1_1_0 + tim1_1_1;
	  c_im(out[6 * ostride]) = tim1_1_0 - tre1_1_1;
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_0_1;
	  FFTW_REAL tim1_0_1;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_1_1;
	  FFTW_REAL tim1_1_1;
	  tre1_0_0 = tre0_1_0 - tim0_1_2;
	  tim1_0_0 = tim0_1_0 + tre0_1_2;
	  tre1_1_0 = tre0_1_0 + tim0_1_2;
	  tim1_1_0 = tim0_1_0 - tre0_1_2;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tre2_1_0;
	       FFTW_REAL tim2_1_0;
	       tre2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_1_1 - tim0_1_1);
	       tim2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_1_1 + tre0_1_1);
	       tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_1_3 + tim0_1_3);
	       tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_1_3 - tim0_1_3);
	       tre1_0_1 = tre2_0_0 - tre2_1_0;
	       tim1_0_1 = tim2_0_0 + tim2_1_0;
	       tre1_1_1 = tre2_0_0 + tre2_1_0;
	       tim1_1_1 = tim2_0_0 - tim2_1_0;
	  }
	  c_re(out[ostride]) = tre1_0_0 + tre1_0_1;
	  c_im(out[ostride]) = tim1_0_0 + tim1_0_1;
	  c_re(out[5 * ostride]) = tre1_0_0 - tre1_0_1;
	  c_im(out[5 * ostride]) = tim1_0_0 - tim1_0_1;
	  c_re(out[3 * ostride]) = tre1_1_0 - tim1_1_1;
	  c_im(out[3 * ostride]) = tim1_1_0 + tre1_1_1;
	  c_re(out[7 * ostride]) = tre1_1_0 + tim1_1_1;
	  c_im(out[7 * ostride]) = tim1_1_0 - tre1_1_1;
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 92 FP additions and 40 FP multiplications */

void fftwi_no_twiddle_9(const FFTW_COMPLEX *in, FFTW_COMPLEX *out, int istride, int ostride)
{
     FFTW_REAL tre0_0_0;
     FFTW_REAL tim0_0_0;
     FFTW_REAL tre0_0_1;
     FFTW_REAL tim0_0_1;
     FFTW_REAL tre0_0_2;
     FFTW_REAL tim0_0_2;
     FFTW_REAL tre0_1_0;
     FFTW_REAL tim0_1_0;
     FFTW_REAL tre0_1_1;
     FFTW_REAL tim0_1_1;
     FFTW_REAL tre0_1_2;
     FFTW_REAL tim0_1_2;
     FFTW_REAL tre0_2_0;
     FFTW_REAL tim0_2_0;
     FFTW_REAL tre0_2_1;
     FFTW_REAL tim0_2_1;
     FFTW_REAL tre0_2_2;
     FFTW_REAL tim0_2_2;
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_2_0;
	  FFTW_REAL tim1_2_0;
	  tre1_0_0 = c_re(in[0]);
	  tim1_0_0 = c_im(in[0]);
	  tre1_1_0 = c_re(in[3 * istride]);
	  tim1_1_0 = c_im(in[3 * istride]);
	  tre1_2_0 = c_re(in[6 * istride]);
	  tim1_2_0 = c_im(in[6 * istride]);
	  tre0_0_0 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	  tim0_0_0 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_2_0 - tim1_1_0);
	       tre0_1_0 = tre2_0_0 + tre2_1_0;
	       tre0_2_0 = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_1_0 - tre1_2_0);
	       tim0_1_0 = tim2_0_0 + tim2_1_0;
	       tim0_2_0 = tim2_0_0 - tim2_1_0;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_2_0;
	  FFTW_REAL tim1_2_0;
	  tre1_0_0 = c_re(in[istride]);
	  tim1_0_0 = c_im(in[istride]);
	  tre1_1_0 = c_re(in[4 * istride]);
	  tim1_1_0 = c_im(in[4 * istride]);
	  tre1_2_0 = c_re(in[7 * istride]);
	  tim1_2_0 = c_im(in[7 * istride]);
	  tre0_0_1 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	  tim0_0_1 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_2_0 - tim1_1_0);
	       tre0_1_1 = tre2_0_0 + tre2_1_0;
	       tre0_2_1 = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_1_0 - tre1_2_0);
	       tim0_1_1 = tim2_0_0 + tim2_1_0;
	       tim0_2_1 = tim2_0_0 - tim2_1_0;
	  }
     }
     {
	  FFTW_REAL tre1_0_0;
	  FFTW_REAL tim1_0_0;
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_2_0;
	  FFTW_REAL tim1_2_0;
	  tre1_0_0 = c_re(in[2 * istride]);
	  tim1_0_0 = c_im(in[2 * istride]);
	  tre1_1_0 = c_re(in[5 * istride]);
	  tim1_1_0 = c_im(in[5 * istride]);
	  tre1_2_0 = c_re(in[8 * istride]);
	  tim1_2_0 = c_im(in[8 * istride]);
	  tre0_0_2 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	  tim0_0_2 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_2_0 - tim1_1_0);
	       tre0_1_2 = tre2_0_0 + tre2_1_0;
	       tre0_2_2 = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_1_0 - tre1_2_0);
	       tim0_1_2 = tim2_0_0 + tim2_1_0;
	       tim0_2_2 = tim2_0_0 - tim2_1_0;
	  }
     }
     c_re(out[0]) = tre0_0_0 + tre0_0_1 + tre0_0_2;
     c_im(out[0]) = tim0_0_0 + tim0_0_1 + tim0_0_2;
     {
	  FFTW_REAL tre2_0_0;
	  FFTW_REAL tre2_1_0;
	  tre2_0_0 = tre0_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre0_0_1 + tre0_0_2));
	  tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim0_0_2 - tim0_0_1);
	  c_re(out[3 * ostride]) = tre2_0_0 + tre2_1_0;
	  c_re(out[6 * ostride]) = tre2_0_0 - tre2_1_0;
     }
     {
	  FFTW_REAL tim2_0_0;
	  FFTW_REAL tim2_1_0;
	  tim2_0_0 = tim0_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim0_0_1 + tim0_0_2));
	  tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre0_0_1 - tre0_0_2);
	  c_im(out[3 * ostride]) = tim2_0_0 + tim2_1_0;
	  c_im(out[6 * ostride]) = tim2_0_0 - tim2_1_0;
     }
     {
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_2_0;
	  FFTW_REAL tim1_2_0;
	  tre1_1_0 = (((FFTW_REAL) FFTW_K766044443) * tre0_1_1) - (((FFTW_REAL) FFTW_K642787609) * tim0_1_1);
	  tim1_1_0 = (((FFTW_REAL) FFTW_K766044443) * tim0_1_1) + (((FFTW_REAL) FFTW_K642787609) * tre0_1_1);
	  tre1_2_0 = (((FFTW_REAL) FFTW_K173648177) * tre0_1_2) - (((FFTW_REAL) FFTW_K984807753) * tim0_1_2);
	  tim1_2_0 = (((FFTW_REAL) FFTW_K173648177) * tim0_1_2) + (((FFTW_REAL) FFTW_K984807753) * tre0_1_2);
	  c_re(out[ostride]) = tre0_1_0 + tre1_1_0 + tre1_2_0;
	  c_im(out[ostride]) = tim0_1_0 + tim1_1_0 + tim1_2_0;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre0_1_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_2_0 - tim1_1_0);
	       c_re(out[4 * ostride]) = tre2_0_0 + tre2_1_0;
	       c_re(out[7 * ostride]) = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim0_1_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_1_0 - tre1_2_0);
	       c_im(out[4 * ostride]) = tim2_0_0 + tim2_1_0;
	       c_im(out[7 * ostride]) = tim2_0_0 - tim2_1_0;
	  }
     }
     {
	  FFTW_REAL tre1_1_0;
	  FFTW_REAL tim1_1_0;
	  FFTW_REAL tre1_2_0;
	  FFTW_REAL tim1_2_0;
	  tre1_1_0 = (((FFTW_REAL) FFTW_K173648177) * tre0_2_1) - (((FFTW_REAL) FFTW_K984807753) * tim0_2_1);
	  tim1_1_0 = (((FFTW_REAL) FFTW_K173648177) * tim0_2_1) + (((FFTW_REAL) FFTW_K984807753) * tre0_2_1);
	  tre1_2_0 = (((FFTW_REAL) FFTW_K939692620) * tre0_2_2) + (((FFTW_REAL) FFTW_K342020143) * tim0_2_2);
	  tim1_2_0 = (((FFTW_REAL) FFTW_K342020143) * tre0_2_2) - (((FFTW_REAL) FFTW_K939692620) * tim0_2_2);
	  c_re(out[2 * ostride]) = tre0_2_0 + tre1_1_0 - tre1_2_0;
	  c_im(out[2 * ostride]) = tim0_2_0 + tim1_1_0 + tim1_2_0;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre0_2_0 + (((FFTW_REAL) FFTW_K499999999) * (tre1_2_0 - tre1_1_0));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_2_0 - tim1_1_0);
	       c_re(out[5 * ostride]) = tre2_0_0 + tre2_1_0;
	       c_re(out[8 * ostride]) = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim0_2_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_1_0 + tre1_2_0);
	       c_im(out[5 * ostride]) = tim2_0_0 + tim2_1_0;
	       c_im(out[8 * ostride]) = tim2_0_0 - tim2_1_0;
	  }
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 126 FP additions and 68 FP multiplications */

void fftw_twiddle_10(FFTW_COMPLEX *A, const FFTW_COMPLEX *W, int stride, int m, int dist)
{
     int i;
     COMPLEX *inout;
     inout = A;
     for (i = 0; i < m; i = i + 1, inout = inout + dist, W = W + 9) {
	  FFTW_REAL tre0_0_0;
	  FFTW_REAL tim0_0_0;
	  FFTW_REAL tre0_0_1;
	  FFTW_REAL tim0_0_1;
	  FFTW_REAL tre0_0_2;
	  FFTW_REAL tim0_0_2;
	  FFTW_REAL tre0_0_3;
	  FFTW_REAL tim0_0_3;
	  FFTW_REAL tre0_0_4;
	  FFTW_REAL tim0_0_4;
	  FFTW_REAL tre0_1_0;
	  FFTW_REAL tim0_1_0;
	  FFTW_REAL tre0_1_1;
	  FFTW_REAL tim0_1_1;
	  FFTW_REAL tre0_1_2;
	  FFTW_REAL tim0_1_2;
	  FFTW_REAL tre0_1_3;
	  FFTW_REAL tim0_1_3;
	  FFTW_REAL tre0_1_4;
	  FFTW_REAL tim0_1_4;
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       tre1_0_0 = c_re(inout[0]);
	       tim1_0_0 = c_im(inout[0]);
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[5 * stride]);
		    ti = c_im(inout[5 * stride]);
		    twr = c_re(W[4]);
		    twi = c_im(W[4]);
		    tre1_1_0 = (tr * twr) - (ti * twi);
		    tim1_1_0 = (tr * twi) + (ti * twr);
	       }
	       tre0_0_0 = tre1_0_0 + tre1_1_0;
	       tim0_0_0 = tim1_0_0 + tim1_1_0;
	       tre0_1_0 = tre1_0_0 - tre1_1_0;
	       tim0_1_0 = tim1_0_0 - tim1_1_0;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[2 * stride]);
		    ti = c_im(inout[2 * stride]);
		    twr = c_re(W[1]);
		    twi = c_im(W[1]);
		    tre1_0_0 = (tr * twr) - (ti * twi);
		    tim1_0_0 = (tr * twi) + (ti * twr);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[7 * stride]);
		    ti = c_im(inout[7 * stride]);
		    twr = c_re(W[6]);
		    twi = c_im(W[6]);
		    tre1_1_0 = (tr * twr) - (ti * twi);
		    tim1_1_0 = (tr * twi) + (ti * twr);
	       }
	       tre0_0_1 = tre1_0_0 + tre1_1_0;
	       tim0_0_1 = tim1_0_0 + tim1_1_0;
	       tre0_1_1 = tre1_0_0 - tre1_1_0;
	       tim0_1_1 = tim1_0_0 - tim1_1_0;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[4 * stride]);
		    ti = c_im(inout[4 * stride]);
		    twr = c_re(W[3]);
		    twi = c_im(W[3]);
		    tre1_0_0 = (tr * twr) - (ti * twi);
		    tim1_0_0 = (tr * twi) + (ti * twr);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[9 * stride]);
		    ti = c_im(inout[9 * stride]);
		    twr = c_re(W[8]);
		    twi = c_im(W[8]);
		    tre1_1_0 = (tr * twr) - (ti * twi);
		    tim1_1_0 = (tr * twi) + (ti * twr);
	       }
	       tre0_0_2 = tre1_0_0 + tre1_1_0;
	       tim0_0_2 = tim1_0_0 + tim1_1_0;
	       tre0_1_2 = tre1_0_0 - tre1_1_0;
	       tim0_1_2 = tim1_0_0 - tim1_1_0;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[6 * stride]);
		    ti = c_im(inout[6 * stride]);
		    twr = c_re(W[5]);
		    twi = c_im(W[5]);
		    tre1_0_0 = (tr * twr) - (ti * twi);
		    tim1_0_0 = (tr * twi) + (ti * twr);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[stride]);
		    ti = c_im(inout[stride]);
		    twr = c_re(W[0]);
		    twi = c_im(W[0]);
		    tre1_1_0 = (tr * twr) - (ti * twi);
		    tim1_1_0 = (tr * twi) + (ti * twr);
	       }
	       tre0_0_3 = tre1_0_0 + tre1_1_0;
	       tim0_0_3 = tim1_0_0 + tim1_1_0;
	       tre0_1_3 = tre1_0_0 - tre1_1_0;
	       tim0_1_3 = tim1_0_0 - tim1_1_0;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[8 * stride]);
		    ti = c_im(inout[8 * stride]);
		    twr = c_re(W[7]);
		    twi = c_im(W[7]);
		    tre1_0_0 = (tr * twr) - (ti * twi);
		    tim1_0_0 = (tr * twi) + (ti * twr);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[3 * stride]);
		    ti = c_im(inout[3 * stride]);
		    twr = c_re(W[2]);
		    twi = c_im(W[2]);
		    tre1_1_0 = (tr * twr) - (ti * twi);
		    tim1_1_0 = (tr * twi) + (ti * twr);
	       }
	       tre0_0_4 = tre1_0_0 + tre1_1_0;
	       tim0_0_4 = tim1_0_0 + tim1_1_0;
	       tre0_1_4 = tre1_0_0 - tre1_1_0;
	       tim0_1_4 = tim1_0_0 - tim1_1_0;
	  }
	  c_re(inout[0]) = tre0_0_0 + tre0_0_1 + tre0_0_2 + tre0_0_3 + tre0_0_4;
	  c_im(inout[0]) = tim0_0_0 + tim0_0_1 + tim0_0_2 + tim0_0_3 + tim0_0_4;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_0_1 + tre0_0_4)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_0_2 + tre0_0_3));
	       tre2_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tim0_0_1 - tim0_0_4)) + (((FFTW_REAL) FFTW_K587785252) * (tim0_0_2 - tim0_0_3));
	       c_re(inout[6 * stride]) = tre2_0_0 + tre2_1_0;
	       c_re(inout[4 * stride]) = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_0_1 + tim0_0_4)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_0_2 + tim0_0_3));
	       tim2_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tre0_0_4 - tre0_0_1)) + (((FFTW_REAL) FFTW_K587785252) * (tre0_0_3 - tre0_0_2));
	       c_im(inout[6 * stride]) = tim2_0_0 + tim2_1_0;
	       c_im(inout[4 * stride]) = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_0_2 + tre0_0_3)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_0_1 + tre0_0_4));
	       tre2_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tim0_0_1 - tim0_0_4)) + (((FFTW_REAL) FFTW_K951056516) * (tim0_0_3 - tim0_0_2));
	       c_re(inout[2 * stride]) = tre2_0_0 + tre2_1_0;
	       c_re(inout[8 * stride]) = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_0_2 + tim0_0_3)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_0_1 + tim0_0_4));
	       tim2_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tre0_0_4 - tre0_0_1)) + (((FFTW_REAL) FFTW_K951056516) * (tre0_0_2 - tre0_0_3));
	       c_im(inout[2 * stride]) = tim2_0_0 + tim2_1_0;
	       c_im(inout[8 * stride]) = tim2_0_0 - tim2_1_0;
	  }
	  c_re(inout[5 * stride]) = tre0_1_0 + tre0_1_1 + tre0_1_2 + tre0_1_3 + tre0_1_4;
	  c_im(inout[5 * stride]) = tim0_1_0 + tim0_1_1 + tim0_1_2 + tim0_1_3 + tim0_1_4;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre0_1_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_1_1 + tre0_1_4)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_1_2 + tre0_1_3));
	       tre2_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tim0_1_1 - tim0_1_4)) + (((FFTW_REAL) FFTW_K587785252) * (tim0_1_2 - tim0_1_3));
	       c_re(inout[stride]) = tre2_0_0 + tre2_1_0;
	       c_re(inout[9 * stride]) = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim0_1_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_1_1 + tim0_1_4)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_1_2 + tim0_1_3));
	       tim2_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tre0_1_4 - tre0_1_1)) + (((FFTW_REAL) FFTW_K587785252) * (tre0_1_3 - tre0_1_2));
	       c_im(inout[stride]) = tim2_0_0 + tim2_1_0;
	       c_im(inout[9 * stride]) = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre0_1_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_1_2 + tre0_1_3)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_1_1 + tre0_1_4));
	       tre2_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tim0_1_1 - tim0_1_4)) + (((FFTW_REAL) FFTW_K951056516) * (tim0_1_3 - tim0_1_2));
	       c_re(inout[7 * stride]) = tre2_0_0 + tre2_1_0;
	       c_re(inout[3 * stride]) = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim0_1_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_1_2 + tim0_1_3)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_1_1 + tim0_1_4));
	       tim2_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tre0_1_4 - tre0_1_1)) + (((FFTW_REAL) FFTW_K951056516) * (tre0_1_2 - tre0_1_3));
	       c_im(inout[7 * stride]) = tim2_0_0 + tim2_1_0;
	       c_im(inout[3 * stride]) = tim2_0_0 - tim2_1_0;
	  }
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 174 FP additions and 84 FP multiplications */

void fftw_twiddle_16(FFTW_COMPLEX *A, const FFTW_COMPLEX *W, int stride, int m, int dist)
{
     int i;
     COMPLEX *inout;
     inout = A;
     for (i = 0; i < m; i = i + 1, inout = inout + dist, W = W + 15) {
	  FFTW_REAL tre0_0_0;
	  FFTW_REAL tim0_0_0;
	  FFTW_REAL tre0_0_1;
	  FFTW_REAL tim0_0_1;
	  FFTW_REAL tre0_0_2;
	  FFTW_REAL tim0_0_2;
	  FFTW_REAL tre0_0_3;
	  FFTW_REAL tim0_0_3;
	  FFTW_REAL tre0_1_0;
	  FFTW_REAL tim0_1_0;
	  FFTW_REAL tre0_1_1;
	  FFTW_REAL tim0_1_1;
	  FFTW_REAL tre0_1_2;
	  FFTW_REAL tim0_1_2;
	  FFTW_REAL tre0_1_3;
	  FFTW_REAL tim0_1_3;
	  FFTW_REAL tre0_2_0;
	  FFTW_REAL tim0_2_0;
	  FFTW_REAL tre0_2_1;
	  FFTW_REAL tim0_2_1;
	  FFTW_REAL tre0_2_2;
	  FFTW_REAL tim0_2_2;
	  FFTW_REAL tre0_2_3;
	  FFTW_REAL tim0_2_3;
	  FFTW_REAL tre0_3_0;
	  FFTW_REAL tim0_3_0;
	  FFTW_REAL tre0_3_1;
	  FFTW_REAL tim0_3_1;
	  FFTW_REAL tre0_3_2;
	  FFTW_REAL tim0_3_2;
	  FFTW_REAL tre0_3_3;
	  FFTW_REAL tim0_3_3;
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = c_re(inout[0]);
		    tim2_0_0 = c_im(inout[0]);
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[8 * stride]);
			 ti = c_im(inout[8 * stride]);
			 twr = c_re(W[7]);
			 twi = c_im(W[7]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[4 * stride]);
			 ti = c_im(inout[4 * stride]);
			 twr = c_re(W[3]);
			 twi = c_im(W[3]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[12 * stride]);
			 ti = c_im(inout[12 * stride]);
			 twr = c_re(W[11]);
			 twi = c_im(W[11]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       tre0_0_0 = tre1_0_0 + tre1_0_1;
	       tim0_0_0 = tim1_0_0 + tim1_0_1;
	       tre0_2_0 = tre1_0_0 - tre1_0_1;
	       tim0_2_0 = tim1_0_0 - tim1_0_1;
	       tre0_1_0 = tre1_1_0 + tim1_1_1;
	       tim0_1_0 = tim1_1_0 - tre1_1_1;
	       tre0_3_0 = tre1_1_0 - tim1_1_1;
	       tim0_3_0 = tim1_1_0 + tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[stride]);
			 ti = c_im(inout[stride]);
			 twr = c_re(W[0]);
			 twi = c_im(W[0]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[9 * stride]);
			 ti = c_im(inout[9 * stride]);
			 twr = c_re(W[8]);
			 twi = c_im(W[8]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[5 * stride]);
			 ti = c_im(inout[5 * stride]);
			 twr = c_re(W[4]);
			 twi = c_im(W[4]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[13 * stride]);
			 ti = c_im(inout[13 * stride]);
			 twr = c_re(W[12]);
			 twi = c_im(W[12]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       tre0_0_1 = tre1_0_0 + tre1_0_1;
	       tim0_0_1 = tim1_0_0 + tim1_0_1;
	       tre0_2_1 = tre1_0_0 - tre1_0_1;
	       tim0_2_1 = tim1_0_0 - tim1_0_1;
	       tre0_1_1 = tre1_1_0 + tim1_1_1;
	       tim0_1_1 = tim1_1_0 - tre1_1_1;
	       tre0_3_1 = tre1_1_0 - tim1_1_1;
	       tim0_3_1 = tim1_1_0 + tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[2 * stride]);
			 ti = c_im(inout[2 * stride]);
			 twr = c_re(W[1]);
			 twi = c_im(W[1]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[10 * stride]);
			 ti = c_im(inout[10 * stride]);
			 twr = c_re(W[9]);
			 twi = c_im(W[9]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[6 * stride]);
			 ti = c_im(inout[6 * stride]);
			 twr = c_re(W[5]);
			 twi = c_im(W[5]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[14 * stride]);
			 ti = c_im(inout[14 * stride]);
			 twr = c_re(W[13]);
			 twi = c_im(W[13]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       tre0_0_2 = tre1_0_0 + tre1_0_1;
	       tim0_0_2 = tim1_0_0 + tim1_0_1;
	       tre0_2_2 = tre1_0_0 - tre1_0_1;
	       tim0_2_2 = tim1_0_0 - tim1_0_1;
	       tre0_1_2 = tre1_1_0 + tim1_1_1;
	       tim0_1_2 = tim1_1_0 - tre1_1_1;
	       tre0_3_2 = tre1_1_0 - tim1_1_1;
	       tim0_3_2 = tim1_1_0 + tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[3 * stride]);
			 ti = c_im(inout[3 * stride]);
			 twr = c_re(W[2]);
			 twi = c_im(W[2]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[11 * stride]);
			 ti = c_im(inout[11 * stride]);
			 twr = c_re(W[10]);
			 twi = c_im(W[10]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[7 * stride]);
			 ti = c_im(inout[7 * stride]);
			 twr = c_re(W[6]);
			 twi = c_im(W[6]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[15 * stride]);
			 ti = c_im(inout[15 * stride]);
			 twr = c_re(W[14]);
			 twi = c_im(W[14]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       tre0_0_3 = tre1_0_0 + tre1_0_1;
	       tim0_0_3 = tim1_0_0 + tim1_0_1;
	       tre0_2_3 = tre1_0_0 - tre1_0_1;
	       tim0_2_3 = tim1_0_0 - tim1_0_1;
	       tre0_1_3 = tre1_1_0 + tim1_1_1;
	       tim0_1_3 = tim1_1_0 - tre1_1_1;
	       tre0_3_3 = tre1_1_0 - tim1_1_1;
	       tim0_3_3 = tim1_1_0 + tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       tre1_0_0 = tre0_0_0 + tre0_0_2;
	       tim1_0_0 = tim0_0_0 + tim0_0_2;
	       tre1_1_0 = tre0_0_0 - tre0_0_2;
	       tim1_1_0 = tim0_0_0 - tim0_0_2;
	       tre1_0_1 = tre0_0_1 + tre0_0_3;
	       tim1_0_1 = tim0_0_1 + tim0_0_3;
	       tre1_1_1 = tre0_0_1 - tre0_0_3;
	       tim1_1_1 = tim0_0_1 - tim0_0_3;
	       c_re(inout[0]) = tre1_0_0 + tre1_0_1;
	       c_im(inout[0]) = tim1_0_0 + tim1_0_1;
	       c_re(inout[8 * stride]) = tre1_0_0 - tre1_0_1;
	       c_im(inout[8 * stride]) = tim1_0_0 - tim1_0_1;
	       c_re(inout[4 * stride]) = tre1_1_0 + tim1_1_1;
	       c_im(inout[4 * stride]) = tim1_1_0 - tre1_1_1;
	       c_re(inout[12 * stride]) = tre1_1_0 - tim1_1_1;
	       c_im(inout[12 * stride]) = tim1_1_0 + tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       {
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_1_2 + tim0_1_2);
		    tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_1_2 - tre0_1_2);
		    tre1_0_0 = tre0_1_0 + tre2_1_0;
		    tim1_0_0 = tim0_1_0 + tim2_1_0;
		    tre1_1_0 = tre0_1_0 - tre2_1_0;
		    tim1_1_0 = tim0_1_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_1_1) + (((FFTW_REAL) FFTW_K382683432) * tim0_1_1);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_1_1) - (((FFTW_REAL) FFTW_K382683432) * tre0_1_1);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_1_3) + (((FFTW_REAL) FFTW_K923879532) * tim0_1_3);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_1_3) - (((FFTW_REAL) FFTW_K923879532) * tre0_1_3);
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       c_re(inout[stride]) = tre1_0_0 + tre1_0_1;
	       c_im(inout[stride]) = tim1_0_0 + tim1_0_1;
	       c_re(inout[9 * stride]) = tre1_0_0 - tre1_0_1;
	       c_im(inout[9 * stride]) = tim1_0_0 - tim1_0_1;
	       c_re(inout[5 * stride]) = tre1_1_0 + tim1_1_1;
	       c_im(inout[5 * stride]) = tim1_1_0 - tre1_1_1;
	       c_re(inout[13 * stride]) = tre1_1_0 - tim1_1_1;
	       c_im(inout[13 * stride]) = tim1_1_0 + tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       tre1_0_0 = tre0_2_0 + tim0_2_2;
	       tim1_0_0 = tim0_2_0 - tre0_2_2;
	       tre1_1_0 = tre0_2_0 - tim0_2_2;
	       tim1_1_0 = tim0_2_0 + tre0_2_2;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_2_1 + tim0_2_1);
		    tim2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_2_1 - tre0_2_1);
		    tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_2_3 - tre0_2_3);
		    tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_2_3 + tre0_2_3);
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 - tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 + tim2_1_0;
	       }
	       c_re(inout[2 * stride]) = tre1_0_0 + tre1_0_1;
	       c_im(inout[2 * stride]) = tim1_0_0 + tim1_0_1;
	       c_re(inout[10 * stride]) = tre1_0_0 - tre1_0_1;
	       c_im(inout[10 * stride]) = tim1_0_0 - tim1_0_1;
	       c_re(inout[6 * stride]) = tre1_1_0 + tim1_1_1;
	       c_im(inout[6 * stride]) = tim1_1_0 - tre1_1_1;
	       c_re(inout[14 * stride]) = tre1_1_0 - tim1_1_1;
	       c_im(inout[14 * stride]) = tim1_1_0 + tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       {
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_3_2 - tre0_3_2);
		    tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_3_2 + tre0_3_2);
		    tre1_0_0 = tre0_3_0 + tre2_1_0;
		    tim1_0_0 = tim0_3_0 - tim2_1_0;
		    tre1_1_0 = tre0_3_0 - tre2_1_0;
		    tim1_1_0 = tim0_3_0 + tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_3_1) + (((FFTW_REAL) FFTW_K923879532) * tim0_3_1);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_3_1) - (((FFTW_REAL) FFTW_K923879532) * tre0_3_1);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_3_3) + (((FFTW_REAL) FFTW_K382683432) * tim0_3_3);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_3_3) - (((FFTW_REAL) FFTW_K923879532) * tim0_3_3);
		    tre1_0_1 = tre2_0_0 - tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 + tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       c_re(inout[3 * stride]) = tre1_0_0 + tre1_0_1;
	       c_im(inout[3 * stride]) = tim1_0_0 + tim1_0_1;
	       c_re(inout[11 * stride]) = tre1_0_0 - tre1_0_1;
	       c_im(inout[11 * stride]) = tim1_0_0 - tim1_0_1;
	       c_re(inout[7 * stride]) = tre1_1_0 + tim1_1_1;
	       c_im(inout[7 * stride]) = tim1_1_0 - tre1_1_1;
	       c_re(inout[15 * stride]) = tre1_1_0 - tim1_1_1;
	       c_im(inout[15 * stride]) = tim1_1_0 + tre1_1_1;
	  }
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 6 FP additions and 4 FP multiplications */

void fftw_twiddle_2(FFTW_COMPLEX *A, const FFTW_COMPLEX *W, int stride, int m, int dist)
{
     int i;
     COMPLEX *inout;
     inout = A;
     for (i = 0; i < m; i = i + 1, inout = inout + dist, W = W + 1) {
	  FFTW_REAL tre0_0_0;
	  FFTW_REAL tim0_0_0;
	  FFTW_REAL tre0_1_0;
	  FFTW_REAL tim0_1_0;
	  tre0_0_0 = c_re(inout[0]);
	  tim0_0_0 = c_im(inout[0]);
	  {
	       FFTW_REAL tr;
	       FFTW_REAL ti;
	       FFTW_REAL twr;
	       FFTW_REAL twi;
	       tr = c_re(inout[stride]);
	       ti = c_im(inout[stride]);
	       twr = c_re(W[0]);
	       twi = c_im(W[0]);
	       tre0_1_0 = (tr * twr) - (ti * twi);
	       tim0_1_0 = (tr * twi) + (ti * twr);
	  }
	  c_re(inout[0]) = tre0_0_0 + tre0_1_0;
	  c_im(inout[0]) = tim0_0_0 + tim0_1_0;
	  c_re(inout[stride]) = tre0_0_0 - tre0_1_0;
	  c_im(inout[stride]) = tim0_0_0 - tim0_1_0;
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 18 FP additions and 12 FP multiplications */

void fftw_twiddle_3(FFTW_COMPLEX *A, const FFTW_COMPLEX *W, int stride, int m, int dist)
{
     int i;
     COMPLEX *inout;
     inout = A;
     for (i = 0; i < m; i = i + 1, inout = inout + dist, W = W + 2) {
	  FFTW_REAL tre0_0_0;
	  FFTW_REAL tim0_0_0;
	  FFTW_REAL tre0_1_0;
	  FFTW_REAL tim0_1_0;
	  FFTW_REAL tre0_2_0;
	  FFTW_REAL tim0_2_0;
	  tre0_0_0 = c_re(inout[0]);
	  tim0_0_0 = c_im(inout[0]);
	  {
	       FFTW_REAL tr;
	       FFTW_REAL ti;
	       FFTW_REAL twr;
	       FFTW_REAL twi;
	       tr = c_re(inout[stride]);
	       ti = c_im(inout[stride]);
	       twr = c_re(W[0]);
	       twi = c_im(W[0]);
	       tre0_1_0 = (tr * twr) - (ti * twi);
	       tim0_1_0 = (tr * twi) + (ti * twr);
	  }
	  {
	       FFTW_REAL tr;
	       FFTW_REAL ti;
	       FFTW_REAL twr;
	       FFTW_REAL twi;
	       tr = c_re(inout[2 * stride]);
	       ti = c_im(inout[2 * stride]);
	       twr = c_re(W[1]);
	       twi = c_im(W[1]);
	       tre0_2_0 = (tr * twr) - (ti * twi);
	       tim0_2_0 = (tr * twi) + (ti * twr);
	  }
	  c_re(inout[0]) = tre0_0_0 + tre0_1_0 + tre0_2_0;
	  c_im(inout[0]) = tim0_0_0 + tim0_1_0 + tim0_2_0;
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tre1_1_0;
	       tre1_0_0 = tre0_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre0_1_0 + tre0_2_0));
	       tre1_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim0_1_0 - tim0_2_0);
	       c_re(inout[stride]) = tre1_0_0 + tre1_1_0;
	       c_re(inout[2 * stride]) = tre1_0_0 - tre1_1_0;
	  }
	  {
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tim1_1_0;
	       tim1_0_0 = tim0_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim0_1_0 + tim0_2_0));
	       tim1_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre0_2_0 - tre0_1_0);
	       c_im(inout[stride]) = tim1_0_0 + tim1_1_0;
	       c_im(inout[2 * stride]) = tim1_0_0 - tim1_1_0;
	  }
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 438 FP additions and 212 FP multiplications */

void fftw_twiddle_32(FFTW_COMPLEX *A, const FFTW_COMPLEX *W, int stride, int m, int dist)
{
     int i;
     COMPLEX *inout;
     inout = A;
     for (i = 0; i < m; i = i + 1, inout = inout + dist, W = W + 31) {
	  FFTW_REAL tre0_0_0;
	  FFTW_REAL tim0_0_0;
	  FFTW_REAL tre0_0_1;
	  FFTW_REAL tim0_0_1;
	  FFTW_REAL tre0_0_2;
	  FFTW_REAL tim0_0_2;
	  FFTW_REAL tre0_0_3;
	  FFTW_REAL tim0_0_3;
	  FFTW_REAL tre0_0_4;
	  FFTW_REAL tim0_0_4;
	  FFTW_REAL tre0_0_5;
	  FFTW_REAL tim0_0_5;
	  FFTW_REAL tre0_0_6;
	  FFTW_REAL tim0_0_6;
	  FFTW_REAL tre0_0_7;
	  FFTW_REAL tim0_0_7;
	  FFTW_REAL tre0_1_0;
	  FFTW_REAL tim0_1_0;
	  FFTW_REAL tre0_1_1;
	  FFTW_REAL tim0_1_1;
	  FFTW_REAL tre0_1_2;
	  FFTW_REAL tim0_1_2;
	  FFTW_REAL tre0_1_3;
	  FFTW_REAL tim0_1_3;
	  FFTW_REAL tre0_1_4;
	  FFTW_REAL tim0_1_4;
	  FFTW_REAL tre0_1_5;
	  FFTW_REAL tim0_1_5;
	  FFTW_REAL tre0_1_6;
	  FFTW_REAL tim0_1_6;
	  FFTW_REAL tre0_1_7;
	  FFTW_REAL tim0_1_7;
	  FFTW_REAL tre0_2_0;
	  FFTW_REAL tim0_2_0;
	  FFTW_REAL tre0_2_1;
	  FFTW_REAL tim0_2_1;
	  FFTW_REAL tre0_2_2;
	  FFTW_REAL tim0_2_2;
	  FFTW_REAL tre0_2_3;
	  FFTW_REAL tim0_2_3;
	  FFTW_REAL tre0_2_4;
	  FFTW_REAL tim0_2_4;
	  FFTW_REAL tre0_2_5;
	  FFTW_REAL tim0_2_5;
	  FFTW_REAL tre0_2_6;
	  FFTW_REAL tim0_2_6;
	  FFTW_REAL tre0_2_7;
	  FFTW_REAL tim0_2_7;
	  FFTW_REAL tre0_3_0;
	  FFTW_REAL tim0_3_0;
	  FFTW_REAL tre0_3_1;
	  FFTW_REAL tim0_3_1;
	  FFTW_REAL tre0_3_2;
	  FFTW_REAL tim0_3_2;
	  FFTW_REAL tre0_3_3;
	  FFTW_REAL tim0_3_3;
	  FFTW_REAL tre0_3_4;
	  FFTW_REAL tim0_3_4;
	  FFTW_REAL tre0_3_5;
	  FFTW_REAL tim0_3_5;
	  FFTW_REAL tre0_3_6;
	  FFTW_REAL tim0_3_6;
	  FFTW_REAL tre0_3_7;
	  FFTW_REAL tim0_3_7;
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = c_re(inout[0]);
		    tim2_0_0 = c_im(inout[0]);
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[16 * stride]);
			 ti = c_im(inout[16 * stride]);
			 twr = c_re(W[15]);
			 twi = c_im(W[15]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[8 * stride]);
			 ti = c_im(inout[8 * stride]);
			 twr = c_re(W[7]);
			 twi = c_im(W[7]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[24 * stride]);
			 ti = c_im(inout[24 * stride]);
			 twr = c_re(W[23]);
			 twi = c_im(W[23]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       tre0_0_0 = tre1_0_0 + tre1_0_1;
	       tim0_0_0 = tim1_0_0 + tim1_0_1;
	       tre0_2_0 = tre1_0_0 - tre1_0_1;
	       tim0_2_0 = tim1_0_0 - tim1_0_1;
	       tre0_1_0 = tre1_1_0 + tim1_1_1;
	       tim0_1_0 = tim1_1_0 - tre1_1_1;
	       tre0_3_0 = tre1_1_0 - tim1_1_1;
	       tim0_3_0 = tim1_1_0 + tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[stride]);
			 ti = c_im(inout[stride]);
			 twr = c_re(W[0]);
			 twi = c_im(W[0]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[17 * stride]);
			 ti = c_im(inout[17 * stride]);
			 twr = c_re(W[16]);
			 twi = c_im(W[16]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[9 * stride]);
			 ti = c_im(inout[9 * stride]);
			 twr = c_re(W[8]);
			 twi = c_im(W[8]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[25 * stride]);
			 ti = c_im(inout[25 * stride]);
			 twr = c_re(W[24]);
			 twi = c_im(W[24]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       tre0_0_1 = tre1_0_0 + tre1_0_1;
	       tim0_0_1 = tim1_0_0 + tim1_0_1;
	       tre0_2_1 = tre1_0_0 - tre1_0_1;
	       tim0_2_1 = tim1_0_0 - tim1_0_1;
	       tre0_1_1 = tre1_1_0 + tim1_1_1;
	       tim0_1_1 = tim1_1_0 - tre1_1_1;
	       tre0_3_1 = tre1_1_0 - tim1_1_1;
	       tim0_3_1 = tim1_1_0 + tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[2 * stride]);
			 ti = c_im(inout[2 * stride]);
			 twr = c_re(W[1]);
			 twi = c_im(W[1]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[18 * stride]);
			 ti = c_im(inout[18 * stride]);
			 twr = c_re(W[17]);
			 twi = c_im(W[17]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[10 * stride]);
			 ti = c_im(inout[10 * stride]);
			 twr = c_re(W[9]);
			 twi = c_im(W[9]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[26 * stride]);
			 ti = c_im(inout[26 * stride]);
			 twr = c_re(W[25]);
			 twi = c_im(W[25]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       tre0_0_2 = tre1_0_0 + tre1_0_1;
	       tim0_0_2 = tim1_0_0 + tim1_0_1;
	       tre0_2_2 = tre1_0_0 - tre1_0_1;
	       tim0_2_2 = tim1_0_0 - tim1_0_1;
	       tre0_1_2 = tre1_1_0 + tim1_1_1;
	       tim0_1_2 = tim1_1_0 - tre1_1_1;
	       tre0_3_2 = tre1_1_0 - tim1_1_1;
	       tim0_3_2 = tim1_1_0 + tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[3 * stride]);
			 ti = c_im(inout[3 * stride]);
			 twr = c_re(W[2]);
			 twi = c_im(W[2]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[19 * stride]);
			 ti = c_im(inout[19 * stride]);
			 twr = c_re(W[18]);
			 twi = c_im(W[18]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[11 * stride]);
			 ti = c_im(inout[11 * stride]);
			 twr = c_re(W[10]);
			 twi = c_im(W[10]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[27 * stride]);
			 ti = c_im(inout[27 * stride]);
			 twr = c_re(W[26]);
			 twi = c_im(W[26]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       tre0_0_3 = tre1_0_0 + tre1_0_1;
	       tim0_0_3 = tim1_0_0 + tim1_0_1;
	       tre0_2_3 = tre1_0_0 - tre1_0_1;
	       tim0_2_3 = tim1_0_0 - tim1_0_1;
	       tre0_1_3 = tre1_1_0 + tim1_1_1;
	       tim0_1_3 = tim1_1_0 - tre1_1_1;
	       tre0_3_3 = tre1_1_0 - tim1_1_1;
	       tim0_3_3 = tim1_1_0 + tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[4 * stride]);
			 ti = c_im(inout[4 * stride]);
			 twr = c_re(W[3]);
			 twi = c_im(W[3]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[20 * stride]);
			 ti = c_im(inout[20 * stride]);
			 twr = c_re(W[19]);
			 twi = c_im(W[19]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[12 * stride]);
			 ti = c_im(inout[12 * stride]);
			 twr = c_re(W[11]);
			 twi = c_im(W[11]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[28 * stride]);
			 ti = c_im(inout[28 * stride]);
			 twr = c_re(W[27]);
			 twi = c_im(W[27]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       tre0_0_4 = tre1_0_0 + tre1_0_1;
	       tim0_0_4 = tim1_0_0 + tim1_0_1;
	       tre0_2_4 = tre1_0_0 - tre1_0_1;
	       tim0_2_4 = tim1_0_0 - tim1_0_1;
	       tre0_1_4 = tre1_1_0 + tim1_1_1;
	       tim0_1_4 = tim1_1_0 - tre1_1_1;
	       tre0_3_4 = tre1_1_0 - tim1_1_1;
	       tim0_3_4 = tim1_1_0 + tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[5 * stride]);
			 ti = c_im(inout[5 * stride]);
			 twr = c_re(W[4]);
			 twi = c_im(W[4]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[21 * stride]);
			 ti = c_im(inout[21 * stride]);
			 twr = c_re(W[20]);
			 twi = c_im(W[20]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[13 * stride]);
			 ti = c_im(inout[13 * stride]);
			 twr = c_re(W[12]);
			 twi = c_im(W[12]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[29 * stride]);
			 ti = c_im(inout[29 * stride]);
			 twr = c_re(W[28]);
			 twi = c_im(W[28]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       tre0_0_5 = tre1_0_0 + tre1_0_1;
	       tim0_0_5 = tim1_0_0 + tim1_0_1;
	       tre0_2_5 = tre1_0_0 - tre1_0_1;
	       tim0_2_5 = tim1_0_0 - tim1_0_1;
	       tre0_1_5 = tre1_1_0 + tim1_1_1;
	       tim0_1_5 = tim1_1_0 - tre1_1_1;
	       tre0_3_5 = tre1_1_0 - tim1_1_1;
	       tim0_3_5 = tim1_1_0 + tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[6 * stride]);
			 ti = c_im(inout[6 * stride]);
			 twr = c_re(W[5]);
			 twi = c_im(W[5]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[22 * stride]);
			 ti = c_im(inout[22 * stride]);
			 twr = c_re(W[21]);
			 twi = c_im(W[21]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[14 * stride]);
			 ti = c_im(inout[14 * stride]);
			 twr = c_re(W[13]);
			 twi = c_im(W[13]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[30 * stride]);
			 ti = c_im(inout[30 * stride]);
			 twr = c_re(W[29]);
			 twi = c_im(W[29]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       tre0_0_6 = tre1_0_0 + tre1_0_1;
	       tim0_0_6 = tim1_0_0 + tim1_0_1;
	       tre0_2_6 = tre1_0_0 - tre1_0_1;
	       tim0_2_6 = tim1_0_0 - tim1_0_1;
	       tre0_1_6 = tre1_1_0 + tim1_1_1;
	       tim0_1_6 = tim1_1_0 - tre1_1_1;
	       tre0_3_6 = tre1_1_0 - tim1_1_1;
	       tim0_3_6 = tim1_1_0 + tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[7 * stride]);
			 ti = c_im(inout[7 * stride]);
			 twr = c_re(W[6]);
			 twi = c_im(W[6]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[23 * stride]);
			 ti = c_im(inout[23 * stride]);
			 twr = c_re(W[22]);
			 twi = c_im(W[22]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[15 * stride]);
			 ti = c_im(inout[15 * stride]);
			 twr = c_re(W[14]);
			 twi = c_im(W[14]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[31 * stride]);
			 ti = c_im(inout[31 * stride]);
			 twr = c_re(W[30]);
			 twi = c_im(W[30]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       tre0_0_7 = tre1_0_0 + tre1_0_1;
	       tim0_0_7 = tim1_0_0 + tim1_0_1;
	       tre0_2_7 = tre1_0_0 - tre1_0_1;
	       tim0_2_7 = tim1_0_0 - tim1_0_1;
	       tre0_1_7 = tre1_1_0 + tim1_1_1;
	       tim0_1_7 = tim1_1_0 - tre1_1_1;
	       tre0_3_7 = tre1_1_0 - tim1_1_1;
	       tim0_3_7 = tim1_1_0 + tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       tre1_0_0 = tre0_0_0 + tre0_0_4;
	       tim1_0_0 = tim0_0_0 + tim0_0_4;
	       tre1_1_0 = tre0_0_0 - tre0_0_4;
	       tim1_1_0 = tim0_0_0 - tim0_0_4;
	       tre1_0_1 = tre0_0_1 + tre0_0_5;
	       tim1_0_1 = tim0_0_1 + tim0_0_5;
	       tre1_1_1 = tre0_0_1 - tre0_0_5;
	       tim1_1_1 = tim0_0_1 - tim0_0_5;
	       tre1_0_2 = tre0_0_2 + tre0_0_6;
	       tim1_0_2 = tim0_0_2 + tim0_0_6;
	       tre1_1_2 = tre0_0_2 - tre0_0_6;
	       tim1_1_2 = tim0_0_2 - tim0_0_6;
	       tre1_0_3 = tre0_0_3 + tre0_0_7;
	       tim1_0_3 = tim0_0_3 + tim0_0_7;
	       tre1_1_3 = tre0_0_3 - tre0_0_7;
	       tim1_1_3 = tim0_0_3 - tim0_0_7;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    c_re(inout[0]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[0]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[16 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[16 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[8 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[8 * stride]) = tim2_1_0 - tre2_1_1;
		    c_re(inout[24 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[24 * stride]) = tim2_1_0 + tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 + tim1_1_2;
		    tim2_0_0 = tim1_1_0 - tre1_1_2;
		    tre2_1_0 = tre1_1_0 - tim1_1_2;
		    tim2_1_0 = tim1_1_0 + tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
			 tre2_0_1 = tre3_0_0 + tre3_1_0;
			 tim2_0_1 = tim3_0_0 - tim3_1_0;
			 tre2_1_1 = tre3_0_0 - tre3_1_0;
			 tim2_1_1 = tim3_0_0 + tim3_1_0;
		    }
		    c_re(inout[4 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[4 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[20 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[20 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[12 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[12 * stride]) = tim2_1_0 - tre2_1_1;
		    c_re(inout[28 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[28 * stride]) = tim2_1_0 + tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_1_4 + tim0_1_4);
		    tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_1_4 - tre0_1_4);
		    tre1_0_0 = tre0_1_0 + tre2_1_0;
		    tim1_0_0 = tim0_1_0 + tim2_1_0;
		    tre1_1_0 = tre0_1_0 - tre2_1_0;
		    tim1_1_0 = tim0_1_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tre0_1_1) + (((FFTW_REAL) FFTW_K195090322) * tim0_1_1);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tim0_1_1) - (((FFTW_REAL) FFTW_K195090322) * tre0_1_1);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K555570233) * tre0_1_5) + (((FFTW_REAL) FFTW_K831469612) * tim0_1_5);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K555570233) * tim0_1_5) - (((FFTW_REAL) FFTW_K831469612) * tre0_1_5);
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_1_2) + (((FFTW_REAL) FFTW_K382683432) * tim0_1_2);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_1_2) - (((FFTW_REAL) FFTW_K382683432) * tre0_1_2);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_1_6) + (((FFTW_REAL) FFTW_K923879532) * tim0_1_6);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_1_6) - (((FFTW_REAL) FFTW_K923879532) * tre0_1_6);
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_1_3) + (((FFTW_REAL) FFTW_K555570233) * tim0_1_3);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_1_3) - (((FFTW_REAL) FFTW_K555570233) * tre0_1_3);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tre0_1_7) + (((FFTW_REAL) FFTW_K980785280) * tim0_1_7);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tim0_1_7) - (((FFTW_REAL) FFTW_K980785280) * tre0_1_7);
		    tre1_0_3 = tre2_0_0 + tre2_1_0;
		    tim1_0_3 = tim2_0_0 + tim2_1_0;
		    tre1_1_3 = tre2_0_0 - tre2_1_0;
		    tim1_1_3 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    c_re(inout[stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[17 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[17 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[9 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[9 * stride]) = tim2_1_0 - tre2_1_1;
		    c_re(inout[25 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[25 * stride]) = tim2_1_0 + tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 + tim1_1_2;
		    tim2_0_0 = tim1_1_0 - tre1_1_2;
		    tre2_1_0 = tre1_1_0 - tim1_1_2;
		    tim2_1_0 = tim1_1_0 + tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
			 tre2_0_1 = tre3_0_0 + tre3_1_0;
			 tim2_0_1 = tim3_0_0 - tim3_1_0;
			 tre2_1_1 = tre3_0_0 - tre3_1_0;
			 tim2_1_1 = tim3_0_0 + tim3_1_0;
		    }
		    c_re(inout[5 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[5 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[21 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[21 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[13 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[13 * stride]) = tim2_1_0 - tre2_1_1;
		    c_re(inout[29 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[29 * stride]) = tim2_1_0 + tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       tre1_0_0 = tre0_2_0 + tim0_2_4;
	       tim1_0_0 = tim0_2_0 - tre0_2_4;
	       tre1_1_0 = tre0_2_0 - tim0_2_4;
	       tim1_1_0 = tim0_2_0 + tre0_2_4;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_2_1) + (((FFTW_REAL) FFTW_K382683432) * tim0_2_1);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_2_1) - (((FFTW_REAL) FFTW_K382683432) * tre0_2_1);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_2_5) - (((FFTW_REAL) FFTW_K382683432) * tre0_2_5);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_2_5) + (((FFTW_REAL) FFTW_K923879532) * tre0_2_5);
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 - tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 + tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_2_2 + tim0_2_2);
		    tim2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_2_2 - tre0_2_2);
		    tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_2_6 - tre0_2_6);
		    tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_2_6 + tre0_2_6);
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 - tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 + tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_2_3) + (((FFTW_REAL) FFTW_K923879532) * tim0_2_3);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_2_3) - (((FFTW_REAL) FFTW_K923879532) * tre0_2_3);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_2_7) - (((FFTW_REAL) FFTW_K923879532) * tre0_2_7);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_2_7) + (((FFTW_REAL) FFTW_K382683432) * tre0_2_7);
		    tre1_0_3 = tre2_0_0 + tre2_1_0;
		    tim1_0_3 = tim2_0_0 - tim2_1_0;
		    tre1_1_3 = tre2_0_0 - tre2_1_0;
		    tim1_1_3 = tim2_0_0 + tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    c_re(inout[2 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[2 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[18 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[18 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[10 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[10 * stride]) = tim2_1_0 - tre2_1_1;
		    c_re(inout[26 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[26 * stride]) = tim2_1_0 + tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 + tim1_1_2;
		    tim2_0_0 = tim1_1_0 - tre1_1_2;
		    tre2_1_0 = tre1_1_0 - tim1_1_2;
		    tim2_1_0 = tim1_1_0 + tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
			 tre2_0_1 = tre3_0_0 + tre3_1_0;
			 tim2_0_1 = tim3_0_0 - tim3_1_0;
			 tre2_1_1 = tre3_0_0 - tre3_1_0;
			 tim2_1_1 = tim3_0_0 + tim3_1_0;
		    }
		    c_re(inout[6 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[6 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[22 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[22 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[14 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[14 * stride]) = tim2_1_0 - tre2_1_1;
		    c_re(inout[30 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[30 * stride]) = tim2_1_0 + tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_3_4 - tre0_3_4);
		    tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_3_4 + tre0_3_4);
		    tre1_0_0 = tre0_3_0 + tre2_1_0;
		    tim1_0_0 = tim0_3_0 - tim2_1_0;
		    tre1_1_0 = tre0_3_0 - tre2_1_0;
		    tim1_1_0 = tim0_3_0 + tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_3_1) + (((FFTW_REAL) FFTW_K555570233) * tim0_3_1);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_3_1) - (((FFTW_REAL) FFTW_K555570233) * tre0_3_1);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tim0_3_5) - (((FFTW_REAL) FFTW_K980785280) * tre0_3_5);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K980785280) * tim0_3_5) + (((FFTW_REAL) FFTW_K195090322) * tre0_3_5);
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 - tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 + tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_3_2) + (((FFTW_REAL) FFTW_K923879532) * tim0_3_2);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_3_2) - (((FFTW_REAL) FFTW_K923879532) * tre0_3_2);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_3_6) + (((FFTW_REAL) FFTW_K382683432) * tim0_3_6);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_3_6) - (((FFTW_REAL) FFTW_K923879532) * tim0_3_6);
		    tre1_0_2 = tre2_0_0 - tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 + tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tim0_3_3) - (((FFTW_REAL) FFTW_K195090322) * tre0_3_3);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K195090322) * tim0_3_3) + (((FFTW_REAL) FFTW_K980785280) * tre0_3_3);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K555570233) * tre0_3_7) + (((FFTW_REAL) FFTW_K831469612) * tim0_3_7);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_3_7) - (((FFTW_REAL) FFTW_K555570233) * tim0_3_7);
		    tre1_0_3 = tre2_0_0 - tre2_1_0;
		    tim1_0_3 = tim2_1_0 - tim2_0_0;
		    tre1_1_3 = tre2_0_0 + tre2_1_0;
		    tim1_1_3 = (-(tim2_0_0 + tim2_1_0));
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    c_re(inout[3 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[3 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[19 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[19 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[11 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[11 * stride]) = tim2_1_0 - tre2_1_1;
		    c_re(inout[27 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[27 * stride]) = tim2_1_0 + tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 + tim1_1_2;
		    tim2_0_0 = tim1_1_0 - tre1_1_2;
		    tre2_1_0 = tre1_1_0 - tim1_1_2;
		    tim2_1_0 = tim1_1_0 + tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
			 tre2_0_1 = tre3_0_0 + tre3_1_0;
			 tim2_0_1 = tim3_0_0 - tim3_1_0;
			 tre2_1_1 = tre3_0_0 - tre3_1_0;
			 tim2_1_1 = tim3_0_0 + tim3_1_0;
		    }
		    c_re(inout[7 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[7 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[23 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[23 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[15 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[15 * stride]) = tim2_1_0 - tre2_1_1;
		    c_re(inout[31 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[31 * stride]) = tim2_1_0 + tre2_1_1;
	       }
	  }
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 22 FP additions and 12 FP multiplications */

void fftw_twiddle_4(FFTW_COMPLEX *A, const FFTW_COMPLEX *W, int stride, int m, int dist)
{
     int i;
     COMPLEX *inout;
     inout = A;
     for (i = 0; i < m; i = i + 1, inout = inout + dist, W = W + 3) {
	  FFTW_REAL tre0_0_0;
	  FFTW_REAL tim0_0_0;
	  FFTW_REAL tre0_0_1;
	  FFTW_REAL tim0_0_1;
	  FFTW_REAL tre0_1_0;
	  FFTW_REAL tim0_1_0;
	  FFTW_REAL tre0_1_1;
	  FFTW_REAL tim0_1_1;
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       tre1_0_0 = c_re(inout[0]);
	       tim1_0_0 = c_im(inout[0]);
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[2 * stride]);
		    ti = c_im(inout[2 * stride]);
		    twr = c_re(W[1]);
		    twi = c_im(W[1]);
		    tre1_1_0 = (tr * twr) - (ti * twi);
		    tim1_1_0 = (tr * twi) + (ti * twr);
	       }
	       tre0_0_0 = tre1_0_0 + tre1_1_0;
	       tim0_0_0 = tim1_0_0 + tim1_1_0;
	       tre0_1_0 = tre1_0_0 - tre1_1_0;
	       tim0_1_0 = tim1_0_0 - tim1_1_0;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[stride]);
		    ti = c_im(inout[stride]);
		    twr = c_re(W[0]);
		    twi = c_im(W[0]);
		    tre1_0_0 = (tr * twr) - (ti * twi);
		    tim1_0_0 = (tr * twi) + (ti * twr);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[3 * stride]);
		    ti = c_im(inout[3 * stride]);
		    twr = c_re(W[2]);
		    twi = c_im(W[2]);
		    tre1_1_0 = (tr * twr) - (ti * twi);
		    tim1_1_0 = (tr * twi) + (ti * twr);
	       }
	       tre0_0_1 = tre1_0_0 + tre1_1_0;
	       tim0_0_1 = tim1_0_0 + tim1_1_0;
	       tre0_1_1 = tre1_0_0 - tre1_1_0;
	       tim0_1_1 = tim1_0_0 - tim1_1_0;
	  }
	  c_re(inout[0]) = tre0_0_0 + tre0_0_1;
	  c_im(inout[0]) = tim0_0_0 + tim0_0_1;
	  c_re(inout[2 * stride]) = tre0_0_0 - tre0_0_1;
	  c_im(inout[2 * stride]) = tim0_0_0 - tim0_0_1;
	  c_re(inout[stride]) = tre0_1_0 + tim0_1_1;
	  c_im(inout[stride]) = tim0_1_0 - tre0_1_1;
	  c_re(inout[3 * stride]) = tre0_1_0 - tim0_1_1;
	  c_im(inout[3 * stride]) = tim0_1_0 + tre0_1_1;
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 52 FP additions and 32 FP multiplications */

void fftw_twiddle_5(FFTW_COMPLEX *A, const FFTW_COMPLEX *W, int stride, int m, int dist)
{
     int i;
     COMPLEX *inout;
     inout = A;
     for (i = 0; i < m; i = i + 1, inout = inout + dist, W = W + 4) {
	  FFTW_REAL tre0_0_0;
	  FFTW_REAL tim0_0_0;
	  FFTW_REAL tre0_1_0;
	  FFTW_REAL tim0_1_0;
	  FFTW_REAL tre0_2_0;
	  FFTW_REAL tim0_2_0;
	  FFTW_REAL tre0_3_0;
	  FFTW_REAL tim0_3_0;
	  FFTW_REAL tre0_4_0;
	  FFTW_REAL tim0_4_0;
	  tre0_0_0 = c_re(inout[0]);
	  tim0_0_0 = c_im(inout[0]);
	  {
	       FFTW_REAL tr;
	       FFTW_REAL ti;
	       FFTW_REAL twr;
	       FFTW_REAL twi;
	       tr = c_re(inout[stride]);
	       ti = c_im(inout[stride]);
	       twr = c_re(W[0]);
	       twi = c_im(W[0]);
	       tre0_1_0 = (tr * twr) - (ti * twi);
	       tim0_1_0 = (tr * twi) + (ti * twr);
	  }
	  {
	       FFTW_REAL tr;
	       FFTW_REAL ti;
	       FFTW_REAL twr;
	       FFTW_REAL twi;
	       tr = c_re(inout[2 * stride]);
	       ti = c_im(inout[2 * stride]);
	       twr = c_re(W[1]);
	       twi = c_im(W[1]);
	       tre0_2_0 = (tr * twr) - (ti * twi);
	       tim0_2_0 = (tr * twi) + (ti * twr);
	  }
	  {
	       FFTW_REAL tr;
	       FFTW_REAL ti;
	       FFTW_REAL twr;
	       FFTW_REAL twi;
	       tr = c_re(inout[3 * stride]);
	       ti = c_im(inout[3 * stride]);
	       twr = c_re(W[2]);
	       twi = c_im(W[2]);
	       tre0_3_0 = (tr * twr) - (ti * twi);
	       tim0_3_0 = (tr * twi) + (ti * twr);
	  }
	  {
	       FFTW_REAL tr;
	       FFTW_REAL ti;
	       FFTW_REAL twr;
	       FFTW_REAL twi;
	       tr = c_re(inout[4 * stride]);
	       ti = c_im(inout[4 * stride]);
	       twr = c_re(W[3]);
	       twi = c_im(W[3]);
	       tre0_4_0 = (tr * twr) - (ti * twi);
	       tim0_4_0 = (tr * twi) + (ti * twr);
	  }
	  c_re(inout[0]) = tre0_0_0 + tre0_1_0 + tre0_2_0 + tre0_3_0 + tre0_4_0;
	  c_im(inout[0]) = tim0_0_0 + tim0_1_0 + tim0_2_0 + tim0_3_0 + tim0_4_0;
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tre1_1_0;
	       tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_1_0 + tre0_4_0)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_2_0 + tre0_3_0));
	       tre1_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tim0_1_0 - tim0_4_0)) + (((FFTW_REAL) FFTW_K587785252) * (tim0_2_0 - tim0_3_0));
	       c_re(inout[stride]) = tre1_0_0 + tre1_1_0;
	       c_re(inout[4 * stride]) = tre1_0_0 - tre1_1_0;
	  }
	  {
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tim1_1_0;
	       tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_1_0 + tim0_4_0)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_2_0 + tim0_3_0));
	       tim1_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tre0_4_0 - tre0_1_0)) + (((FFTW_REAL) FFTW_K587785252) * (tre0_3_0 - tre0_2_0));
	       c_im(inout[stride]) = tim1_0_0 + tim1_1_0;
	       c_im(inout[4 * stride]) = tim1_0_0 - tim1_1_0;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tre1_1_0;
	       tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_2_0 + tre0_3_0)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_1_0 + tre0_4_0));
	       tre1_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tim0_1_0 - tim0_4_0)) + (((FFTW_REAL) FFTW_K951056516) * (tim0_3_0 - tim0_2_0));
	       c_re(inout[2 * stride]) = tre1_0_0 + tre1_1_0;
	       c_re(inout[3 * stride]) = tre1_0_0 - tre1_1_0;
	  }
	  {
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tim1_1_0;
	       tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_2_0 + tim0_3_0)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_1_0 + tim0_4_0));
	       tim1_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tre0_4_0 - tre0_1_0)) + (((FFTW_REAL) FFTW_K951056516) * (tre0_2_0 - tre0_3_0));
	       c_im(inout[2 * stride]) = tim1_0_0 + tim1_1_0;
	       c_im(inout[3 * stride]) = tim1_0_0 - tim1_1_0;
	  }
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 50 FP additions and 28 FP multiplications */

void fftw_twiddle_6(FFTW_COMPLEX *A, const FFTW_COMPLEX *W, int stride, int m, int dist)
{
     int i;
     COMPLEX *inout;
     inout = A;
     for (i = 0; i < m; i = i + 1, inout = inout + dist, W = W + 5) {
	  FFTW_REAL tre0_0_0;
	  FFTW_REAL tim0_0_0;
	  FFTW_REAL tre0_0_1;
	  FFTW_REAL tim0_0_1;
	  FFTW_REAL tre0_0_2;
	  FFTW_REAL tim0_0_2;
	  FFTW_REAL tre0_1_0;
	  FFTW_REAL tim0_1_0;
	  FFTW_REAL tre0_1_1;
	  FFTW_REAL tim0_1_1;
	  FFTW_REAL tre0_1_2;
	  FFTW_REAL tim0_1_2;
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       tre1_0_0 = c_re(inout[0]);
	       tim1_0_0 = c_im(inout[0]);
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[3 * stride]);
		    ti = c_im(inout[3 * stride]);
		    twr = c_re(W[2]);
		    twi = c_im(W[2]);
		    tre1_1_0 = (tr * twr) - (ti * twi);
		    tim1_1_0 = (tr * twi) + (ti * twr);
	       }
	       tre0_0_0 = tre1_0_0 + tre1_1_0;
	       tim0_0_0 = tim1_0_0 + tim1_1_0;
	       tre0_1_0 = tre1_0_0 - tre1_1_0;
	       tim0_1_0 = tim1_0_0 - tim1_1_0;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[2 * stride]);
		    ti = c_im(inout[2 * stride]);
		    twr = c_re(W[1]);
		    twi = c_im(W[1]);
		    tre1_0_0 = (tr * twr) - (ti * twi);
		    tim1_0_0 = (tr * twi) + (ti * twr);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[5 * stride]);
		    ti = c_im(inout[5 * stride]);
		    twr = c_re(W[4]);
		    twi = c_im(W[4]);
		    tre1_1_0 = (tr * twr) - (ti * twi);
		    tim1_1_0 = (tr * twi) + (ti * twr);
	       }
	       tre0_0_1 = tre1_0_0 + tre1_1_0;
	       tim0_0_1 = tim1_0_0 + tim1_1_0;
	       tre0_1_1 = tre1_0_0 - tre1_1_0;
	       tim0_1_1 = tim1_0_0 - tim1_1_0;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[4 * stride]);
		    ti = c_im(inout[4 * stride]);
		    twr = c_re(W[3]);
		    twi = c_im(W[3]);
		    tre1_0_0 = (tr * twr) - (ti * twi);
		    tim1_0_0 = (tr * twi) + (ti * twr);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[stride]);
		    ti = c_im(inout[stride]);
		    twr = c_re(W[0]);
		    twi = c_im(W[0]);
		    tre1_1_0 = (tr * twr) - (ti * twi);
		    tim1_1_0 = (tr * twi) + (ti * twr);
	       }
	       tre0_0_2 = tre1_0_0 + tre1_1_0;
	       tim0_0_2 = tim1_0_0 + tim1_1_0;
	       tre0_1_2 = tre1_0_0 - tre1_1_0;
	       tim0_1_2 = tim1_0_0 - tim1_1_0;
	  }
	  c_re(inout[0]) = tre0_0_0 + tre0_0_1 + tre0_0_2;
	  c_im(inout[0]) = tim0_0_0 + tim0_0_1 + tim0_0_2;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre0_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre0_0_1 + tre0_0_2));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim0_0_1 - tim0_0_2);
	       c_re(inout[4 * stride]) = tre2_0_0 + tre2_1_0;
	       c_re(inout[2 * stride]) = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim0_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim0_0_1 + tim0_0_2));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre0_0_2 - tre0_0_1);
	       c_im(inout[4 * stride]) = tim2_0_0 + tim2_1_0;
	       c_im(inout[2 * stride]) = tim2_0_0 - tim2_1_0;
	  }
	  c_re(inout[3 * stride]) = tre0_1_0 + tre0_1_1 + tre0_1_2;
	  c_im(inout[3 * stride]) = tim0_1_0 + tim0_1_1 + tim0_1_2;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre0_1_0 - (((FFTW_REAL) FFTW_K499999999) * (tre0_1_1 + tre0_1_2));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim0_1_1 - tim0_1_2);
	       c_re(inout[stride]) = tre2_0_0 + tre2_1_0;
	       c_re(inout[5 * stride]) = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim0_1_0 - (((FFTW_REAL) FFTW_K499999999) * (tim0_1_1 + tim0_1_2));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre0_1_2 - tre0_1_1);
	       c_im(inout[stride]) = tim2_0_0 + tim2_1_0;
	       c_im(inout[5 * stride]) = tim2_0_0 - tim2_1_0;
	  }
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 1054 FP additions and 500 FP multiplications */

void fftw_twiddle_64(FFTW_COMPLEX *A, const FFTW_COMPLEX *W, int stride, int m, int dist)
{
     int i;
     COMPLEX *inout;
     inout = A;
     for (i = 0; i < m; i = i + 1, inout = inout + dist, W = W + 63) {
	  FFTW_REAL tre0_0_0;
	  FFTW_REAL tim0_0_0;
	  FFTW_REAL tre0_0_1;
	  FFTW_REAL tim0_0_1;
	  FFTW_REAL tre0_0_2;
	  FFTW_REAL tim0_0_2;
	  FFTW_REAL tre0_0_3;
	  FFTW_REAL tim0_0_3;
	  FFTW_REAL tre0_0_4;
	  FFTW_REAL tim0_0_4;
	  FFTW_REAL tre0_0_5;
	  FFTW_REAL tim0_0_5;
	  FFTW_REAL tre0_0_6;
	  FFTW_REAL tim0_0_6;
	  FFTW_REAL tre0_0_7;
	  FFTW_REAL tim0_0_7;
	  FFTW_REAL tre0_1_0;
	  FFTW_REAL tim0_1_0;
	  FFTW_REAL tre0_1_1;
	  FFTW_REAL tim0_1_1;
	  FFTW_REAL tre0_1_2;
	  FFTW_REAL tim0_1_2;
	  FFTW_REAL tre0_1_3;
	  FFTW_REAL tim0_1_3;
	  FFTW_REAL tre0_1_4;
	  FFTW_REAL tim0_1_4;
	  FFTW_REAL tre0_1_5;
	  FFTW_REAL tim0_1_5;
	  FFTW_REAL tre0_1_6;
	  FFTW_REAL tim0_1_6;
	  FFTW_REAL tre0_1_7;
	  FFTW_REAL tim0_1_7;
	  FFTW_REAL tre0_2_0;
	  FFTW_REAL tim0_2_0;
	  FFTW_REAL tre0_2_1;
	  FFTW_REAL tim0_2_1;
	  FFTW_REAL tre0_2_2;
	  FFTW_REAL tim0_2_2;
	  FFTW_REAL tre0_2_3;
	  FFTW_REAL tim0_2_3;
	  FFTW_REAL tre0_2_4;
	  FFTW_REAL tim0_2_4;
	  FFTW_REAL tre0_2_5;
	  FFTW_REAL tim0_2_5;
	  FFTW_REAL tre0_2_6;
	  FFTW_REAL tim0_2_6;
	  FFTW_REAL tre0_2_7;
	  FFTW_REAL tim0_2_7;
	  FFTW_REAL tre0_3_0;
	  FFTW_REAL tim0_3_0;
	  FFTW_REAL tre0_3_1;
	  FFTW_REAL tim0_3_1;
	  FFTW_REAL tre0_3_2;
	  FFTW_REAL tim0_3_2;
	  FFTW_REAL tre0_3_3;
	  FFTW_REAL tim0_3_3;
	  FFTW_REAL tre0_3_4;
	  FFTW_REAL tim0_3_4;
	  FFTW_REAL tre0_3_5;
	  FFTW_REAL tim0_3_5;
	  FFTW_REAL tre0_3_6;
	  FFTW_REAL tim0_3_6;
	  FFTW_REAL tre0_3_7;
	  FFTW_REAL tim0_3_7;
	  FFTW_REAL tre0_4_0;
	  FFTW_REAL tim0_4_0;
	  FFTW_REAL tre0_4_1;
	  FFTW_REAL tim0_4_1;
	  FFTW_REAL tre0_4_2;
	  FFTW_REAL tim0_4_2;
	  FFTW_REAL tre0_4_3;
	  FFTW_REAL tim0_4_3;
	  FFTW_REAL tre0_4_4;
	  FFTW_REAL tim0_4_4;
	  FFTW_REAL tre0_4_5;
	  FFTW_REAL tim0_4_5;
	  FFTW_REAL tre0_4_6;
	  FFTW_REAL tim0_4_6;
	  FFTW_REAL tre0_4_7;
	  FFTW_REAL tim0_4_7;
	  FFTW_REAL tre0_5_0;
	  FFTW_REAL tim0_5_0;
	  FFTW_REAL tre0_5_1;
	  FFTW_REAL tim0_5_1;
	  FFTW_REAL tre0_5_2;
	  FFTW_REAL tim0_5_2;
	  FFTW_REAL tre0_5_3;
	  FFTW_REAL tim0_5_3;
	  FFTW_REAL tre0_5_4;
	  FFTW_REAL tim0_5_4;
	  FFTW_REAL tre0_5_5;
	  FFTW_REAL tim0_5_5;
	  FFTW_REAL tre0_5_6;
	  FFTW_REAL tim0_5_6;
	  FFTW_REAL tre0_5_7;
	  FFTW_REAL tim0_5_7;
	  FFTW_REAL tre0_6_0;
	  FFTW_REAL tim0_6_0;
	  FFTW_REAL tre0_6_1;
	  FFTW_REAL tim0_6_1;
	  FFTW_REAL tre0_6_2;
	  FFTW_REAL tim0_6_2;
	  FFTW_REAL tre0_6_3;
	  FFTW_REAL tim0_6_3;
	  FFTW_REAL tre0_6_4;
	  FFTW_REAL tim0_6_4;
	  FFTW_REAL tre0_6_5;
	  FFTW_REAL tim0_6_5;
	  FFTW_REAL tre0_6_6;
	  FFTW_REAL tim0_6_6;
	  FFTW_REAL tre0_6_7;
	  FFTW_REAL tim0_6_7;
	  FFTW_REAL tre0_7_0;
	  FFTW_REAL tim0_7_0;
	  FFTW_REAL tre0_7_1;
	  FFTW_REAL tim0_7_1;
	  FFTW_REAL tre0_7_2;
	  FFTW_REAL tim0_7_2;
	  FFTW_REAL tre0_7_3;
	  FFTW_REAL tim0_7_3;
	  FFTW_REAL tre0_7_4;
	  FFTW_REAL tim0_7_4;
	  FFTW_REAL tre0_7_5;
	  FFTW_REAL tim0_7_5;
	  FFTW_REAL tre0_7_6;
	  FFTW_REAL tim0_7_6;
	  FFTW_REAL tre0_7_7;
	  FFTW_REAL tim0_7_7;
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = c_re(inout[0]);
		    tim2_0_0 = c_im(inout[0]);
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[32 * stride]);
			 ti = c_im(inout[32 * stride]);
			 twr = c_re(W[31]);
			 twi = c_im(W[31]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[8 * stride]);
			 ti = c_im(inout[8 * stride]);
			 twr = c_re(W[7]);
			 twi = c_im(W[7]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[40 * stride]);
			 ti = c_im(inout[40 * stride]);
			 twr = c_re(W[39]);
			 twi = c_im(W[39]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[16 * stride]);
			 ti = c_im(inout[16 * stride]);
			 twr = c_re(W[15]);
			 twi = c_im(W[15]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[48 * stride]);
			 ti = c_im(inout[48 * stride]);
			 twr = c_re(W[47]);
			 twi = c_im(W[47]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[24 * stride]);
			 ti = c_im(inout[24 * stride]);
			 twr = c_re(W[23]);
			 twi = c_im(W[23]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[56 * stride]);
			 ti = c_im(inout[56 * stride]);
			 twr = c_re(W[55]);
			 twi = c_im(W[55]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_3 = tre2_0_0 + tre2_1_0;
		    tim1_0_3 = tim2_0_0 + tim2_1_0;
		    tre1_1_3 = tre2_0_0 - tre2_1_0;
		    tim1_1_3 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    tre0_0_0 = tre2_0_0 + tre2_0_1;
		    tim0_0_0 = tim2_0_0 + tim2_0_1;
		    tre0_4_0 = tre2_0_0 - tre2_0_1;
		    tim0_4_0 = tim2_0_0 - tim2_0_1;
		    tre0_2_0 = tre2_1_0 + tim2_1_1;
		    tim0_2_0 = tim2_1_0 - tre2_1_1;
		    tre0_6_0 = tre2_1_0 - tim2_1_1;
		    tim0_6_0 = tim2_1_0 + tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 + tim1_1_2;
		    tim2_0_0 = tim1_1_0 - tre1_1_2;
		    tre2_1_0 = tre1_1_0 - tim1_1_2;
		    tim2_1_0 = tim1_1_0 + tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
			 tre2_0_1 = tre3_0_0 + tre3_1_0;
			 tim2_0_1 = tim3_0_0 - tim3_1_0;
			 tre2_1_1 = tre3_0_0 - tre3_1_0;
			 tim2_1_1 = tim3_0_0 + tim3_1_0;
		    }
		    tre0_1_0 = tre2_0_0 + tre2_0_1;
		    tim0_1_0 = tim2_0_0 + tim2_0_1;
		    tre0_5_0 = tre2_0_0 - tre2_0_1;
		    tim0_5_0 = tim2_0_0 - tim2_0_1;
		    tre0_3_0 = tre2_1_0 + tim2_1_1;
		    tim0_3_0 = tim2_1_0 - tre2_1_1;
		    tre0_7_0 = tre2_1_0 - tim2_1_1;
		    tim0_7_0 = tim2_1_0 + tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[stride]);
			 ti = c_im(inout[stride]);
			 twr = c_re(W[0]);
			 twi = c_im(W[0]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[33 * stride]);
			 ti = c_im(inout[33 * stride]);
			 twr = c_re(W[32]);
			 twi = c_im(W[32]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[9 * stride]);
			 ti = c_im(inout[9 * stride]);
			 twr = c_re(W[8]);
			 twi = c_im(W[8]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[41 * stride]);
			 ti = c_im(inout[41 * stride]);
			 twr = c_re(W[40]);
			 twi = c_im(W[40]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[17 * stride]);
			 ti = c_im(inout[17 * stride]);
			 twr = c_re(W[16]);
			 twi = c_im(W[16]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[49 * stride]);
			 ti = c_im(inout[49 * stride]);
			 twr = c_re(W[48]);
			 twi = c_im(W[48]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[25 * stride]);
			 ti = c_im(inout[25 * stride]);
			 twr = c_re(W[24]);
			 twi = c_im(W[24]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[57 * stride]);
			 ti = c_im(inout[57 * stride]);
			 twr = c_re(W[56]);
			 twi = c_im(W[56]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_3 = tre2_0_0 + tre2_1_0;
		    tim1_0_3 = tim2_0_0 + tim2_1_0;
		    tre1_1_3 = tre2_0_0 - tre2_1_0;
		    tim1_1_3 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    tre0_0_1 = tre2_0_0 + tre2_0_1;
		    tim0_0_1 = tim2_0_0 + tim2_0_1;
		    tre0_4_1 = tre2_0_0 - tre2_0_1;
		    tim0_4_1 = tim2_0_0 - tim2_0_1;
		    tre0_2_1 = tre2_1_0 + tim2_1_1;
		    tim0_2_1 = tim2_1_0 - tre2_1_1;
		    tre0_6_1 = tre2_1_0 - tim2_1_1;
		    tim0_6_1 = tim2_1_0 + tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 + tim1_1_2;
		    tim2_0_0 = tim1_1_0 - tre1_1_2;
		    tre2_1_0 = tre1_1_0 - tim1_1_2;
		    tim2_1_0 = tim1_1_0 + tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
			 tre2_0_1 = tre3_0_0 + tre3_1_0;
			 tim2_0_1 = tim3_0_0 - tim3_1_0;
			 tre2_1_1 = tre3_0_0 - tre3_1_0;
			 tim2_1_1 = tim3_0_0 + tim3_1_0;
		    }
		    tre0_1_1 = tre2_0_0 + tre2_0_1;
		    tim0_1_1 = tim2_0_0 + tim2_0_1;
		    tre0_5_1 = tre2_0_0 - tre2_0_1;
		    tim0_5_1 = tim2_0_0 - tim2_0_1;
		    tre0_3_1 = tre2_1_0 + tim2_1_1;
		    tim0_3_1 = tim2_1_0 - tre2_1_1;
		    tre0_7_1 = tre2_1_0 - tim2_1_1;
		    tim0_7_1 = tim2_1_0 + tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[2 * stride]);
			 ti = c_im(inout[2 * stride]);
			 twr = c_re(W[1]);
			 twi = c_im(W[1]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[34 * stride]);
			 ti = c_im(inout[34 * stride]);
			 twr = c_re(W[33]);
			 twi = c_im(W[33]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[10 * stride]);
			 ti = c_im(inout[10 * stride]);
			 twr = c_re(W[9]);
			 twi = c_im(W[9]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[42 * stride]);
			 ti = c_im(inout[42 * stride]);
			 twr = c_re(W[41]);
			 twi = c_im(W[41]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[18 * stride]);
			 ti = c_im(inout[18 * stride]);
			 twr = c_re(W[17]);
			 twi = c_im(W[17]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[50 * stride]);
			 ti = c_im(inout[50 * stride]);
			 twr = c_re(W[49]);
			 twi = c_im(W[49]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[26 * stride]);
			 ti = c_im(inout[26 * stride]);
			 twr = c_re(W[25]);
			 twi = c_im(W[25]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[58 * stride]);
			 ti = c_im(inout[58 * stride]);
			 twr = c_re(W[57]);
			 twi = c_im(W[57]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_3 = tre2_0_0 + tre2_1_0;
		    tim1_0_3 = tim2_0_0 + tim2_1_0;
		    tre1_1_3 = tre2_0_0 - tre2_1_0;
		    tim1_1_3 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    tre0_0_2 = tre2_0_0 + tre2_0_1;
		    tim0_0_2 = tim2_0_0 + tim2_0_1;
		    tre0_4_2 = tre2_0_0 - tre2_0_1;
		    tim0_4_2 = tim2_0_0 - tim2_0_1;
		    tre0_2_2 = tre2_1_0 + tim2_1_1;
		    tim0_2_2 = tim2_1_0 - tre2_1_1;
		    tre0_6_2 = tre2_1_0 - tim2_1_1;
		    tim0_6_2 = tim2_1_0 + tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 + tim1_1_2;
		    tim2_0_0 = tim1_1_0 - tre1_1_2;
		    tre2_1_0 = tre1_1_0 - tim1_1_2;
		    tim2_1_0 = tim1_1_0 + tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
			 tre2_0_1 = tre3_0_0 + tre3_1_0;
			 tim2_0_1 = tim3_0_0 - tim3_1_0;
			 tre2_1_1 = tre3_0_0 - tre3_1_0;
			 tim2_1_1 = tim3_0_0 + tim3_1_0;
		    }
		    tre0_1_2 = tre2_0_0 + tre2_0_1;
		    tim0_1_2 = tim2_0_0 + tim2_0_1;
		    tre0_5_2 = tre2_0_0 - tre2_0_1;
		    tim0_5_2 = tim2_0_0 - tim2_0_1;
		    tre0_3_2 = tre2_1_0 + tim2_1_1;
		    tim0_3_2 = tim2_1_0 - tre2_1_1;
		    tre0_7_2 = tre2_1_0 - tim2_1_1;
		    tim0_7_2 = tim2_1_0 + tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[3 * stride]);
			 ti = c_im(inout[3 * stride]);
			 twr = c_re(W[2]);
			 twi = c_im(W[2]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[35 * stride]);
			 ti = c_im(inout[35 * stride]);
			 twr = c_re(W[34]);
			 twi = c_im(W[34]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[11 * stride]);
			 ti = c_im(inout[11 * stride]);
			 twr = c_re(W[10]);
			 twi = c_im(W[10]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[43 * stride]);
			 ti = c_im(inout[43 * stride]);
			 twr = c_re(W[42]);
			 twi = c_im(W[42]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[19 * stride]);
			 ti = c_im(inout[19 * stride]);
			 twr = c_re(W[18]);
			 twi = c_im(W[18]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[51 * stride]);
			 ti = c_im(inout[51 * stride]);
			 twr = c_re(W[50]);
			 twi = c_im(W[50]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[27 * stride]);
			 ti = c_im(inout[27 * stride]);
			 twr = c_re(W[26]);
			 twi = c_im(W[26]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[59 * stride]);
			 ti = c_im(inout[59 * stride]);
			 twr = c_re(W[58]);
			 twi = c_im(W[58]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_3 = tre2_0_0 + tre2_1_0;
		    tim1_0_3 = tim2_0_0 + tim2_1_0;
		    tre1_1_3 = tre2_0_0 - tre2_1_0;
		    tim1_1_3 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    tre0_0_3 = tre2_0_0 + tre2_0_1;
		    tim0_0_3 = tim2_0_0 + tim2_0_1;
		    tre0_4_3 = tre2_0_0 - tre2_0_1;
		    tim0_4_3 = tim2_0_0 - tim2_0_1;
		    tre0_2_3 = tre2_1_0 + tim2_1_1;
		    tim0_2_3 = tim2_1_0 - tre2_1_1;
		    tre0_6_3 = tre2_1_0 - tim2_1_1;
		    tim0_6_3 = tim2_1_0 + tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 + tim1_1_2;
		    tim2_0_0 = tim1_1_0 - tre1_1_2;
		    tre2_1_0 = tre1_1_0 - tim1_1_2;
		    tim2_1_0 = tim1_1_0 + tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
			 tre2_0_1 = tre3_0_0 + tre3_1_0;
			 tim2_0_1 = tim3_0_0 - tim3_1_0;
			 tre2_1_1 = tre3_0_0 - tre3_1_0;
			 tim2_1_1 = tim3_0_0 + tim3_1_0;
		    }
		    tre0_1_3 = tre2_0_0 + tre2_0_1;
		    tim0_1_3 = tim2_0_0 + tim2_0_1;
		    tre0_5_3 = tre2_0_0 - tre2_0_1;
		    tim0_5_3 = tim2_0_0 - tim2_0_1;
		    tre0_3_3 = tre2_1_0 + tim2_1_1;
		    tim0_3_3 = tim2_1_0 - tre2_1_1;
		    tre0_7_3 = tre2_1_0 - tim2_1_1;
		    tim0_7_3 = tim2_1_0 + tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[4 * stride]);
			 ti = c_im(inout[4 * stride]);
			 twr = c_re(W[3]);
			 twi = c_im(W[3]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[36 * stride]);
			 ti = c_im(inout[36 * stride]);
			 twr = c_re(W[35]);
			 twi = c_im(W[35]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[12 * stride]);
			 ti = c_im(inout[12 * stride]);
			 twr = c_re(W[11]);
			 twi = c_im(W[11]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[44 * stride]);
			 ti = c_im(inout[44 * stride]);
			 twr = c_re(W[43]);
			 twi = c_im(W[43]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[20 * stride]);
			 ti = c_im(inout[20 * stride]);
			 twr = c_re(W[19]);
			 twi = c_im(W[19]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[52 * stride]);
			 ti = c_im(inout[52 * stride]);
			 twr = c_re(W[51]);
			 twi = c_im(W[51]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[28 * stride]);
			 ti = c_im(inout[28 * stride]);
			 twr = c_re(W[27]);
			 twi = c_im(W[27]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[60 * stride]);
			 ti = c_im(inout[60 * stride]);
			 twr = c_re(W[59]);
			 twi = c_im(W[59]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_3 = tre2_0_0 + tre2_1_0;
		    tim1_0_3 = tim2_0_0 + tim2_1_0;
		    tre1_1_3 = tre2_0_0 - tre2_1_0;
		    tim1_1_3 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    tre0_0_4 = tre2_0_0 + tre2_0_1;
		    tim0_0_4 = tim2_0_0 + tim2_0_1;
		    tre0_4_4 = tre2_0_0 - tre2_0_1;
		    tim0_4_4 = tim2_0_0 - tim2_0_1;
		    tre0_2_4 = tre2_1_0 + tim2_1_1;
		    tim0_2_4 = tim2_1_0 - tre2_1_1;
		    tre0_6_4 = tre2_1_0 - tim2_1_1;
		    tim0_6_4 = tim2_1_0 + tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 + tim1_1_2;
		    tim2_0_0 = tim1_1_0 - tre1_1_2;
		    tre2_1_0 = tre1_1_0 - tim1_1_2;
		    tim2_1_0 = tim1_1_0 + tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
			 tre2_0_1 = tre3_0_0 + tre3_1_0;
			 tim2_0_1 = tim3_0_0 - tim3_1_0;
			 tre2_1_1 = tre3_0_0 - tre3_1_0;
			 tim2_1_1 = tim3_0_0 + tim3_1_0;
		    }
		    tre0_1_4 = tre2_0_0 + tre2_0_1;
		    tim0_1_4 = tim2_0_0 + tim2_0_1;
		    tre0_5_4 = tre2_0_0 - tre2_0_1;
		    tim0_5_4 = tim2_0_0 - tim2_0_1;
		    tre0_3_4 = tre2_1_0 + tim2_1_1;
		    tim0_3_4 = tim2_1_0 - tre2_1_1;
		    tre0_7_4 = tre2_1_0 - tim2_1_1;
		    tim0_7_4 = tim2_1_0 + tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[5 * stride]);
			 ti = c_im(inout[5 * stride]);
			 twr = c_re(W[4]);
			 twi = c_im(W[4]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[37 * stride]);
			 ti = c_im(inout[37 * stride]);
			 twr = c_re(W[36]);
			 twi = c_im(W[36]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[13 * stride]);
			 ti = c_im(inout[13 * stride]);
			 twr = c_re(W[12]);
			 twi = c_im(W[12]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[45 * stride]);
			 ti = c_im(inout[45 * stride]);
			 twr = c_re(W[44]);
			 twi = c_im(W[44]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[21 * stride]);
			 ti = c_im(inout[21 * stride]);
			 twr = c_re(W[20]);
			 twi = c_im(W[20]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[53 * stride]);
			 ti = c_im(inout[53 * stride]);
			 twr = c_re(W[52]);
			 twi = c_im(W[52]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[29 * stride]);
			 ti = c_im(inout[29 * stride]);
			 twr = c_re(W[28]);
			 twi = c_im(W[28]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[61 * stride]);
			 ti = c_im(inout[61 * stride]);
			 twr = c_re(W[60]);
			 twi = c_im(W[60]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_3 = tre2_0_0 + tre2_1_0;
		    tim1_0_3 = tim2_0_0 + tim2_1_0;
		    tre1_1_3 = tre2_0_0 - tre2_1_0;
		    tim1_1_3 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    tre0_0_5 = tre2_0_0 + tre2_0_1;
		    tim0_0_5 = tim2_0_0 + tim2_0_1;
		    tre0_4_5 = tre2_0_0 - tre2_0_1;
		    tim0_4_5 = tim2_0_0 - tim2_0_1;
		    tre0_2_5 = tre2_1_0 + tim2_1_1;
		    tim0_2_5 = tim2_1_0 - tre2_1_1;
		    tre0_6_5 = tre2_1_0 - tim2_1_1;
		    tim0_6_5 = tim2_1_0 + tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 + tim1_1_2;
		    tim2_0_0 = tim1_1_0 - tre1_1_2;
		    tre2_1_0 = tre1_1_0 - tim1_1_2;
		    tim2_1_0 = tim1_1_0 + tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
			 tre2_0_1 = tre3_0_0 + tre3_1_0;
			 tim2_0_1 = tim3_0_0 - tim3_1_0;
			 tre2_1_1 = tre3_0_0 - tre3_1_0;
			 tim2_1_1 = tim3_0_0 + tim3_1_0;
		    }
		    tre0_1_5 = tre2_0_0 + tre2_0_1;
		    tim0_1_5 = tim2_0_0 + tim2_0_1;
		    tre0_5_5 = tre2_0_0 - tre2_0_1;
		    tim0_5_5 = tim2_0_0 - tim2_0_1;
		    tre0_3_5 = tre2_1_0 + tim2_1_1;
		    tim0_3_5 = tim2_1_0 - tre2_1_1;
		    tre0_7_5 = tre2_1_0 - tim2_1_1;
		    tim0_7_5 = tim2_1_0 + tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[6 * stride]);
			 ti = c_im(inout[6 * stride]);
			 twr = c_re(W[5]);
			 twi = c_im(W[5]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[38 * stride]);
			 ti = c_im(inout[38 * stride]);
			 twr = c_re(W[37]);
			 twi = c_im(W[37]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[14 * stride]);
			 ti = c_im(inout[14 * stride]);
			 twr = c_re(W[13]);
			 twi = c_im(W[13]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[46 * stride]);
			 ti = c_im(inout[46 * stride]);
			 twr = c_re(W[45]);
			 twi = c_im(W[45]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[22 * stride]);
			 ti = c_im(inout[22 * stride]);
			 twr = c_re(W[21]);
			 twi = c_im(W[21]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[54 * stride]);
			 ti = c_im(inout[54 * stride]);
			 twr = c_re(W[53]);
			 twi = c_im(W[53]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[30 * stride]);
			 ti = c_im(inout[30 * stride]);
			 twr = c_re(W[29]);
			 twi = c_im(W[29]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[62 * stride]);
			 ti = c_im(inout[62 * stride]);
			 twr = c_re(W[61]);
			 twi = c_im(W[61]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_3 = tre2_0_0 + tre2_1_0;
		    tim1_0_3 = tim2_0_0 + tim2_1_0;
		    tre1_1_3 = tre2_0_0 - tre2_1_0;
		    tim1_1_3 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    tre0_0_6 = tre2_0_0 + tre2_0_1;
		    tim0_0_6 = tim2_0_0 + tim2_0_1;
		    tre0_4_6 = tre2_0_0 - tre2_0_1;
		    tim0_4_6 = tim2_0_0 - tim2_0_1;
		    tre0_2_6 = tre2_1_0 + tim2_1_1;
		    tim0_2_6 = tim2_1_0 - tre2_1_1;
		    tre0_6_6 = tre2_1_0 - tim2_1_1;
		    tim0_6_6 = tim2_1_0 + tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 + tim1_1_2;
		    tim2_0_0 = tim1_1_0 - tre1_1_2;
		    tre2_1_0 = tre1_1_0 - tim1_1_2;
		    tim2_1_0 = tim1_1_0 + tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
			 tre2_0_1 = tre3_0_0 + tre3_1_0;
			 tim2_0_1 = tim3_0_0 - tim3_1_0;
			 tre2_1_1 = tre3_0_0 - tre3_1_0;
			 tim2_1_1 = tim3_0_0 + tim3_1_0;
		    }
		    tre0_1_6 = tre2_0_0 + tre2_0_1;
		    tim0_1_6 = tim2_0_0 + tim2_0_1;
		    tre0_5_6 = tre2_0_0 - tre2_0_1;
		    tim0_5_6 = tim2_0_0 - tim2_0_1;
		    tre0_3_6 = tre2_1_0 + tim2_1_1;
		    tim0_3_6 = tim2_1_0 - tre2_1_1;
		    tre0_7_6 = tre2_1_0 - tim2_1_1;
		    tim0_7_6 = tim2_1_0 + tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[7 * stride]);
			 ti = c_im(inout[7 * stride]);
			 twr = c_re(W[6]);
			 twi = c_im(W[6]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[39 * stride]);
			 ti = c_im(inout[39 * stride]);
			 twr = c_re(W[38]);
			 twi = c_im(W[38]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[15 * stride]);
			 ti = c_im(inout[15 * stride]);
			 twr = c_re(W[14]);
			 twi = c_im(W[14]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[47 * stride]);
			 ti = c_im(inout[47 * stride]);
			 twr = c_re(W[46]);
			 twi = c_im(W[46]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[23 * stride]);
			 ti = c_im(inout[23 * stride]);
			 twr = c_re(W[22]);
			 twi = c_im(W[22]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[55 * stride]);
			 ti = c_im(inout[55 * stride]);
			 twr = c_re(W[54]);
			 twi = c_im(W[54]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[31 * stride]);
			 ti = c_im(inout[31 * stride]);
			 twr = c_re(W[30]);
			 twi = c_im(W[30]);
			 tre2_0_0 = (tr * twr) - (ti * twi);
			 tim2_0_0 = (tr * twi) + (ti * twr);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[63 * stride]);
			 ti = c_im(inout[63 * stride]);
			 twr = c_re(W[62]);
			 twi = c_im(W[62]);
			 tre2_1_0 = (tr * twr) - (ti * twi);
			 tim2_1_0 = (tr * twi) + (ti * twr);
		    }
		    tre1_0_3 = tre2_0_0 + tre2_1_0;
		    tim1_0_3 = tim2_0_0 + tim2_1_0;
		    tre1_1_3 = tre2_0_0 - tre2_1_0;
		    tim1_1_3 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    tre0_0_7 = tre2_0_0 + tre2_0_1;
		    tim0_0_7 = tim2_0_0 + tim2_0_1;
		    tre0_4_7 = tre2_0_0 - tre2_0_1;
		    tim0_4_7 = tim2_0_0 - tim2_0_1;
		    tre0_2_7 = tre2_1_0 + tim2_1_1;
		    tim0_2_7 = tim2_1_0 - tre2_1_1;
		    tre0_6_7 = tre2_1_0 - tim2_1_1;
		    tim0_6_7 = tim2_1_0 + tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 + tim1_1_2;
		    tim2_0_0 = tim1_1_0 - tre1_1_2;
		    tre2_1_0 = tre1_1_0 - tim1_1_2;
		    tim2_1_0 = tim1_1_0 + tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
			 tre2_0_1 = tre3_0_0 + tre3_1_0;
			 tim2_0_1 = tim3_0_0 - tim3_1_0;
			 tre2_1_1 = tre3_0_0 - tre3_1_0;
			 tim2_1_1 = tim3_0_0 + tim3_1_0;
		    }
		    tre0_1_7 = tre2_0_0 + tre2_0_1;
		    tim0_1_7 = tim2_0_0 + tim2_0_1;
		    tre0_5_7 = tre2_0_0 - tre2_0_1;
		    tim0_5_7 = tim2_0_0 - tim2_0_1;
		    tre0_3_7 = tre2_1_0 + tim2_1_1;
		    tim0_3_7 = tim2_1_0 - tre2_1_1;
		    tre0_7_7 = tre2_1_0 - tim2_1_1;
		    tim0_7_7 = tim2_1_0 + tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       tre1_0_0 = tre0_0_0 + tre0_0_4;
	       tim1_0_0 = tim0_0_0 + tim0_0_4;
	       tre1_1_0 = tre0_0_0 - tre0_0_4;
	       tim1_1_0 = tim0_0_0 - tim0_0_4;
	       tre1_0_1 = tre0_0_1 + tre0_0_5;
	       tim1_0_1 = tim0_0_1 + tim0_0_5;
	       tre1_1_1 = tre0_0_1 - tre0_0_5;
	       tim1_1_1 = tim0_0_1 - tim0_0_5;
	       tre1_0_2 = tre0_0_2 + tre0_0_6;
	       tim1_0_2 = tim0_0_2 + tim0_0_6;
	       tre1_1_2 = tre0_0_2 - tre0_0_6;
	       tim1_1_2 = tim0_0_2 - tim0_0_6;
	       tre1_0_3 = tre0_0_3 + tre0_0_7;
	       tim1_0_3 = tim0_0_3 + tim0_0_7;
	       tre1_1_3 = tre0_0_3 - tre0_0_7;
	       tim1_1_3 = tim0_0_3 - tim0_0_7;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    c_re(inout[0]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[0]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[32 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[32 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[16 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[16 * stride]) = tim2_1_0 - tre2_1_1;
		    c_re(inout[48 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[48 * stride]) = tim2_1_0 + tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 + tim1_1_2;
		    tim2_0_0 = tim1_1_0 - tre1_1_2;
		    tre2_1_0 = tre1_1_0 - tim1_1_2;
		    tim2_1_0 = tim1_1_0 + tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
			 tre2_0_1 = tre3_0_0 + tre3_1_0;
			 tim2_0_1 = tim3_0_0 - tim3_1_0;
			 tre2_1_1 = tre3_0_0 - tre3_1_0;
			 tim2_1_1 = tim3_0_0 + tim3_1_0;
		    }
		    c_re(inout[8 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[8 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[40 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[40 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[24 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[24 * stride]) = tim2_1_0 - tre2_1_1;
		    c_re(inout[56 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[56 * stride]) = tim2_1_0 + tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_1_4) + (((FFTW_REAL) FFTW_K382683432) * tim0_1_4);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_1_4) - (((FFTW_REAL) FFTW_K382683432) * tre0_1_4);
		    tre1_0_0 = tre0_1_0 + tre2_1_0;
		    tim1_0_0 = tim0_1_0 + tim2_1_0;
		    tre1_1_0 = tre0_1_0 - tre2_1_0;
		    tim1_1_0 = tim0_1_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K995184726) * tre0_1_1) + (((FFTW_REAL) FFTW_K098017140) * tim0_1_1);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K995184726) * tim0_1_1) - (((FFTW_REAL) FFTW_K098017140) * tre0_1_1);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K881921264) * tre0_1_5) + (((FFTW_REAL) FFTW_K471396736) * tim0_1_5);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K881921264) * tim0_1_5) - (((FFTW_REAL) FFTW_K471396736) * tre0_1_5);
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tre0_1_2) + (((FFTW_REAL) FFTW_K195090322) * tim0_1_2);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tim0_1_2) - (((FFTW_REAL) FFTW_K195090322) * tre0_1_2);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_1_6) + (((FFTW_REAL) FFTW_K555570233) * tim0_1_6);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_1_6) - (((FFTW_REAL) FFTW_K555570233) * tre0_1_6);
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K956940335) * tre0_1_3) + (((FFTW_REAL) FFTW_K290284677) * tim0_1_3);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K956940335) * tim0_1_3) - (((FFTW_REAL) FFTW_K290284677) * tre0_1_3);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K773010453) * tre0_1_7) + (((FFTW_REAL) FFTW_K634393284) * tim0_1_7);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K773010453) * tim0_1_7) - (((FFTW_REAL) FFTW_K634393284) * tre0_1_7);
		    tre1_0_3 = tre2_0_0 + tre2_1_0;
		    tim1_0_3 = tim2_0_0 + tim2_1_0;
		    tre1_1_3 = tre2_0_0 - tre2_1_0;
		    tim1_1_3 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    c_re(inout[stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[33 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[33 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[17 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[17 * stride]) = tim2_1_0 - tre2_1_1;
		    c_re(inout[49 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[49 * stride]) = tim2_1_0 + tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 + tim1_1_2;
		    tim2_0_0 = tim1_1_0 - tre1_1_2;
		    tre2_1_0 = tre1_1_0 - tim1_1_2;
		    tim2_1_0 = tim1_1_0 + tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
			 tre2_0_1 = tre3_0_0 + tre3_1_0;
			 tim2_0_1 = tim3_0_0 - tim3_1_0;
			 tre2_1_1 = tre3_0_0 - tre3_1_0;
			 tim2_1_1 = tim3_0_0 + tim3_1_0;
		    }
		    c_re(inout[9 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[9 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[41 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[41 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[25 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[25 * stride]) = tim2_1_0 - tre2_1_1;
		    c_re(inout[57 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[57 * stride]) = tim2_1_0 + tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_2_4 + tim0_2_4);
		    tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_2_4 - tre0_2_4);
		    tre1_0_0 = tre0_2_0 + tre2_1_0;
		    tim1_0_0 = tim0_2_0 + tim2_1_0;
		    tre1_1_0 = tre0_2_0 - tre2_1_0;
		    tim1_1_0 = tim0_2_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tre0_2_1) + (((FFTW_REAL) FFTW_K195090322) * tim0_2_1);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tim0_2_1) - (((FFTW_REAL) FFTW_K195090322) * tre0_2_1);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K555570233) * tre0_2_5) + (((FFTW_REAL) FFTW_K831469612) * tim0_2_5);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K555570233) * tim0_2_5) - (((FFTW_REAL) FFTW_K831469612) * tre0_2_5);
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_2_2) + (((FFTW_REAL) FFTW_K382683432) * tim0_2_2);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_2_2) - (((FFTW_REAL) FFTW_K382683432) * tre0_2_2);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_2_6) + (((FFTW_REAL) FFTW_K923879532) * tim0_2_6);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_2_6) - (((FFTW_REAL) FFTW_K923879532) * tre0_2_6);
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_2_3) + (((FFTW_REAL) FFTW_K555570233) * tim0_2_3);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_2_3) - (((FFTW_REAL) FFTW_K555570233) * tre0_2_3);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tre0_2_7) + (((FFTW_REAL) FFTW_K980785280) * tim0_2_7);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tim0_2_7) - (((FFTW_REAL) FFTW_K980785280) * tre0_2_7);
		    tre1_0_3 = tre2_0_0 + tre2_1_0;
		    tim1_0_3 = tim2_0_0 + tim2_1_0;
		    tre1_1_3 = tre2_0_0 - tre2_1_0;
		    tim1_1_3 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    c_re(inout[2 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[2 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[34 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[34 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[18 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[18 * stride]) = tim2_1_0 - tre2_1_1;
		    c_re(inout[50 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[50 * stride]) = tim2_1_0 + tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 + tim1_1_2;
		    tim2_0_0 = tim1_1_0 - tre1_1_2;
		    tre2_1_0 = tre1_1_0 - tim1_1_2;
		    tim2_1_0 = tim1_1_0 + tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
			 tre2_0_1 = tre3_0_0 + tre3_1_0;
			 tim2_0_1 = tim3_0_0 - tim3_1_0;
			 tre2_1_1 = tre3_0_0 - tre3_1_0;
			 tim2_1_1 = tim3_0_0 + tim3_1_0;
		    }
		    c_re(inout[10 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[10 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[42 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[42 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[26 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[26 * stride]) = tim2_1_0 - tre2_1_1;
		    c_re(inout[58 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[58 * stride]) = tim2_1_0 + tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_3_4) + (((FFTW_REAL) FFTW_K923879532) * tim0_3_4);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_3_4) - (((FFTW_REAL) FFTW_K923879532) * tre0_3_4);
		    tre1_0_0 = tre0_3_0 + tre2_1_0;
		    tim1_0_0 = tim0_3_0 + tim2_1_0;
		    tre1_1_0 = tre0_3_0 - tre2_1_0;
		    tim1_1_0 = tim0_3_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K956940335) * tre0_3_1) + (((FFTW_REAL) FFTW_K290284677) * tim0_3_1);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K956940335) * tim0_3_1) - (((FFTW_REAL) FFTW_K290284677) * tre0_3_1);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K098017140) * tre0_3_5) + (((FFTW_REAL) FFTW_K995184726) * tim0_3_5);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K098017140) * tim0_3_5) - (((FFTW_REAL) FFTW_K995184726) * tre0_3_5);
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_3_2) + (((FFTW_REAL) FFTW_K555570233) * tim0_3_2);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_3_2) - (((FFTW_REAL) FFTW_K555570233) * tre0_3_2);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K980785280) * tim0_3_6) - (((FFTW_REAL) FFTW_K195090322) * tre0_3_6);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tim0_3_6) + (((FFTW_REAL) FFTW_K980785280) * tre0_3_6);
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 - tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 + tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K634393284) * tre0_3_3) + (((FFTW_REAL) FFTW_K773010453) * tim0_3_3);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K634393284) * tim0_3_3) - (((FFTW_REAL) FFTW_K773010453) * tre0_3_3);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K881921264) * tim0_3_7) - (((FFTW_REAL) FFTW_K471396736) * tre0_3_7);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K471396736) * tim0_3_7) + (((FFTW_REAL) FFTW_K881921264) * tre0_3_7);
		    tre1_0_3 = tre2_0_0 + tre2_1_0;
		    tim1_0_3 = tim2_0_0 - tim2_1_0;
		    tre1_1_3 = tre2_0_0 - tre2_1_0;
		    tim1_1_3 = tim2_0_0 + tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    c_re(inout[3 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[3 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[35 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[35 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[19 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[19 * stride]) = tim2_1_0 - tre2_1_1;
		    c_re(inout[51 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[51 * stride]) = tim2_1_0 + tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 + tim1_1_2;
		    tim2_0_0 = tim1_1_0 - tre1_1_2;
		    tre2_1_0 = tre1_1_0 - tim1_1_2;
		    tim2_1_0 = tim1_1_0 + tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
			 tre2_0_1 = tre3_0_0 + tre3_1_0;
			 tim2_0_1 = tim3_0_0 - tim3_1_0;
			 tre2_1_1 = tre3_0_0 - tre3_1_0;
			 tim2_1_1 = tim3_0_0 + tim3_1_0;
		    }
		    c_re(inout[11 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[11 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[43 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[43 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[27 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[27 * stride]) = tim2_1_0 - tre2_1_1;
		    c_re(inout[59 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[59 * stride]) = tim2_1_0 + tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       tre1_0_0 = tre0_4_0 + tim0_4_4;
	       tim1_0_0 = tim0_4_0 - tre0_4_4;
	       tre1_1_0 = tre0_4_0 - tim0_4_4;
	       tim1_1_0 = tim0_4_0 + tre0_4_4;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_4_1) + (((FFTW_REAL) FFTW_K382683432) * tim0_4_1);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_4_1) - (((FFTW_REAL) FFTW_K382683432) * tre0_4_1);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_4_5) - (((FFTW_REAL) FFTW_K382683432) * tre0_4_5);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_4_5) + (((FFTW_REAL) FFTW_K923879532) * tre0_4_5);
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 - tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 + tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_4_2 + tim0_4_2);
		    tim2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_4_2 - tre0_4_2);
		    tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_4_6 - tre0_4_6);
		    tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_4_6 + tre0_4_6);
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 - tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 + tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_4_3) + (((FFTW_REAL) FFTW_K923879532) * tim0_4_3);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_4_3) - (((FFTW_REAL) FFTW_K923879532) * tre0_4_3);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_4_7) - (((FFTW_REAL) FFTW_K923879532) * tre0_4_7);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_4_7) + (((FFTW_REAL) FFTW_K382683432) * tre0_4_7);
		    tre1_0_3 = tre2_0_0 + tre2_1_0;
		    tim1_0_3 = tim2_0_0 - tim2_1_0;
		    tre1_1_3 = tre2_0_0 - tre2_1_0;
		    tim1_1_3 = tim2_0_0 + tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    c_re(inout[4 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[4 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[36 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[36 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[20 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[20 * stride]) = tim2_1_0 - tre2_1_1;
		    c_re(inout[52 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[52 * stride]) = tim2_1_0 + tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 + tim1_1_2;
		    tim2_0_0 = tim1_1_0 - tre1_1_2;
		    tre2_1_0 = tre1_1_0 - tim1_1_2;
		    tim2_1_0 = tim1_1_0 + tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
			 tre2_0_1 = tre3_0_0 + tre3_1_0;
			 tim2_0_1 = tim3_0_0 - tim3_1_0;
			 tre2_1_1 = tre3_0_0 - tre3_1_0;
			 tim2_1_1 = tim3_0_0 + tim3_1_0;
		    }
		    c_re(inout[12 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[12 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[44 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[44 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[28 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[28 * stride]) = tim2_1_0 - tre2_1_1;
		    c_re(inout[60 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[60 * stride]) = tim2_1_0 + tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_5_4) - (((FFTW_REAL) FFTW_K382683432) * tre0_5_4);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_5_4) + (((FFTW_REAL) FFTW_K923879532) * tre0_5_4);
		    tre1_0_0 = tre0_5_0 + tre2_1_0;
		    tim1_0_0 = tim0_5_0 - tim2_1_0;
		    tre1_1_0 = tre0_5_0 - tre2_1_0;
		    tim1_1_0 = tim0_5_0 + tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K881921264) * tre0_5_1) + (((FFTW_REAL) FFTW_K471396736) * tim0_5_1);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K881921264) * tim0_5_1) - (((FFTW_REAL) FFTW_K471396736) * tre0_5_1);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K634393284) * tim0_5_5) - (((FFTW_REAL) FFTW_K773010453) * tre0_5_5);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K773010453) * tim0_5_5) + (((FFTW_REAL) FFTW_K634393284) * tre0_5_5);
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 - tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 + tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K555570233) * tre0_5_2) + (((FFTW_REAL) FFTW_K831469612) * tim0_5_2);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K555570233) * tim0_5_2) - (((FFTW_REAL) FFTW_K831469612) * tre0_5_2);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tim0_5_6) - (((FFTW_REAL) FFTW_K980785280) * tre0_5_6);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K980785280) * tim0_5_6) + (((FFTW_REAL) FFTW_K195090322) * tre0_5_6);
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 - tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 + tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K098017140) * tre0_5_3) + (((FFTW_REAL) FFTW_K995184726) * tim0_5_3);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K098017140) * tim0_5_3) - (((FFTW_REAL) FFTW_K995184726) * tre0_5_3);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K956940335) * tre0_5_7) + (((FFTW_REAL) FFTW_K290284677) * tim0_5_7);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K290284677) * tre0_5_7) - (((FFTW_REAL) FFTW_K956940335) * tim0_5_7);
		    tre1_0_3 = tre2_0_0 - tre2_1_0;
		    tim1_0_3 = tim2_0_0 + tim2_1_0;
		    tre1_1_3 = tre2_0_0 + tre2_1_0;
		    tim1_1_3 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    c_re(inout[5 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[5 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[37 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[37 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[21 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[21 * stride]) = tim2_1_0 - tre2_1_1;
		    c_re(inout[53 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[53 * stride]) = tim2_1_0 + tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 + tim1_1_2;
		    tim2_0_0 = tim1_1_0 - tre1_1_2;
		    tre2_1_0 = tre1_1_0 - tim1_1_2;
		    tim2_1_0 = tim1_1_0 + tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
			 tre2_0_1 = tre3_0_0 + tre3_1_0;
			 tim2_0_1 = tim3_0_0 - tim3_1_0;
			 tre2_1_1 = tre3_0_0 - tre3_1_0;
			 tim2_1_1 = tim3_0_0 + tim3_1_0;
		    }
		    c_re(inout[13 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[13 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[45 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[45 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[29 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[29 * stride]) = tim2_1_0 - tre2_1_1;
		    c_re(inout[61 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[61 * stride]) = tim2_1_0 + tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_6_4 - tre0_6_4);
		    tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_6_4 + tre0_6_4);
		    tre1_0_0 = tre0_6_0 + tre2_1_0;
		    tim1_0_0 = tim0_6_0 - tim2_1_0;
		    tre1_1_0 = tre0_6_0 - tre2_1_0;
		    tim1_1_0 = tim0_6_0 + tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_6_1) + (((FFTW_REAL) FFTW_K555570233) * tim0_6_1);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_6_1) - (((FFTW_REAL) FFTW_K555570233) * tre0_6_1);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tim0_6_5) - (((FFTW_REAL) FFTW_K980785280) * tre0_6_5);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K980785280) * tim0_6_5) + (((FFTW_REAL) FFTW_K195090322) * tre0_6_5);
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 - tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 + tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_6_2) + (((FFTW_REAL) FFTW_K923879532) * tim0_6_2);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_6_2) - (((FFTW_REAL) FFTW_K923879532) * tre0_6_2);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_6_6) + (((FFTW_REAL) FFTW_K382683432) * tim0_6_6);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_6_6) - (((FFTW_REAL) FFTW_K923879532) * tim0_6_6);
		    tre1_0_2 = tre2_0_0 - tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 + tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tim0_6_3) - (((FFTW_REAL) FFTW_K195090322) * tre0_6_3);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K195090322) * tim0_6_3) + (((FFTW_REAL) FFTW_K980785280) * tre0_6_3);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K555570233) * tre0_6_7) + (((FFTW_REAL) FFTW_K831469612) * tim0_6_7);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_6_7) - (((FFTW_REAL) FFTW_K555570233) * tim0_6_7);
		    tre1_0_3 = tre2_0_0 - tre2_1_0;
		    tim1_0_3 = tim2_1_0 - tim2_0_0;
		    tre1_1_3 = tre2_0_0 + tre2_1_0;
		    tim1_1_3 = (-(tim2_0_0 + tim2_1_0));
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    c_re(inout[6 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[6 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[38 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[38 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[22 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[22 * stride]) = tim2_1_0 - tre2_1_1;
		    c_re(inout[54 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[54 * stride]) = tim2_1_0 + tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 + tim1_1_2;
		    tim2_0_0 = tim1_1_0 - tre1_1_2;
		    tre2_1_0 = tre1_1_0 - tim1_1_2;
		    tim2_1_0 = tim1_1_0 + tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
			 tre2_0_1 = tre3_0_0 + tre3_1_0;
			 tim2_0_1 = tim3_0_0 - tim3_1_0;
			 tre2_1_1 = tre3_0_0 - tre3_1_0;
			 tim2_1_1 = tim3_0_0 + tim3_1_0;
		    }
		    c_re(inout[14 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[14 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[46 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[46 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[30 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[30 * stride]) = tim2_1_0 - tre2_1_1;
		    c_re(inout[62 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[62 * stride]) = tim2_1_0 + tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_7_4) - (((FFTW_REAL) FFTW_K923879532) * tre0_7_4);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_7_4) + (((FFTW_REAL) FFTW_K382683432) * tre0_7_4);
		    tre1_0_0 = tre0_7_0 + tre2_1_0;
		    tim1_0_0 = tim0_7_0 - tim2_1_0;
		    tre1_1_0 = tre0_7_0 - tre2_1_0;
		    tim1_1_0 = tim0_7_0 + tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K773010453) * tre0_7_1) + (((FFTW_REAL) FFTW_K634393284) * tim0_7_1);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K773010453) * tim0_7_1) - (((FFTW_REAL) FFTW_K634393284) * tre0_7_1);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K956940335) * tre0_7_5) + (((FFTW_REAL) FFTW_K290284677) * tim0_7_5);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K290284677) * tre0_7_5) - (((FFTW_REAL) FFTW_K956940335) * tim0_7_5);
		    tre1_0_1 = tre2_0_0 - tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 + tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K195090322) * tre0_7_2) + (((FFTW_REAL) FFTW_K980785280) * tim0_7_2);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K195090322) * tim0_7_2) - (((FFTW_REAL) FFTW_K980785280) * tre0_7_2);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K555570233) * tre0_7_6) + (((FFTW_REAL) FFTW_K831469612) * tim0_7_6);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_7_6) - (((FFTW_REAL) FFTW_K555570233) * tim0_7_6);
		    tre1_0_2 = tre2_0_0 - tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 + tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K881921264) * tim0_7_3) - (((FFTW_REAL) FFTW_K471396736) * tre0_7_3);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K471396736) * tim0_7_3) + (((FFTW_REAL) FFTW_K881921264) * tre0_7_3);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K098017140) * tre0_7_7) - (((FFTW_REAL) FFTW_K995184726) * tim0_7_7);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K098017140) * tim0_7_7) + (((FFTW_REAL) FFTW_K995184726) * tre0_7_7);
		    tre1_0_3 = tre2_0_0 + tre2_1_0;
		    tim1_0_3 = tim2_1_0 - tim2_0_0;
		    tre1_1_3 = tre2_0_0 - tre2_1_0;
		    tim1_1_3 = (-(tim2_0_0 + tim2_1_0));
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    c_re(inout[7 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[7 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[39 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[39 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[23 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[23 * stride]) = tim2_1_0 - tre2_1_1;
		    c_re(inout[55 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[55 * stride]) = tim2_1_0 + tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 + tim1_1_2;
		    tim2_0_0 = tim1_1_0 - tre1_1_2;
		    tre2_1_0 = tre1_1_0 - tim1_1_2;
		    tim2_1_0 = tim1_1_0 + tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 + tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 - tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 - tre1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_3 + tre1_1_3);
			 tre2_0_1 = tre3_0_0 + tre3_1_0;
			 tim2_0_1 = tim3_0_0 - tim3_1_0;
			 tre2_1_1 = tre3_0_0 - tre3_1_0;
			 tim2_1_1 = tim3_0_0 + tim3_1_0;
		    }
		    c_re(inout[15 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[15 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[47 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[47 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[31 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[31 * stride]) = tim2_1_0 - tre2_1_1;
		    c_re(inout[63 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[63 * stride]) = tim2_1_0 + tre2_1_1;
	       }
	  }
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 102 FP additions and 60 FP multiplications */

void fftw_twiddle_7(FFTW_COMPLEX *A, const FFTW_COMPLEX *W, int stride, int m, int dist)
{
     int i;
     COMPLEX *inout;
     inout = A;
     for (i = 0; i < m; i = i + 1, inout = inout + dist, W = W + 6) {
	  FFTW_REAL tre0_0_0;
	  FFTW_REAL tim0_0_0;
	  FFTW_REAL tre0_1_0;
	  FFTW_REAL tim0_1_0;
	  FFTW_REAL tre0_2_0;
	  FFTW_REAL tim0_2_0;
	  FFTW_REAL tre0_3_0;
	  FFTW_REAL tim0_3_0;
	  FFTW_REAL tre0_4_0;
	  FFTW_REAL tim0_4_0;
	  FFTW_REAL tre0_5_0;
	  FFTW_REAL tim0_5_0;
	  FFTW_REAL tre0_6_0;
	  FFTW_REAL tim0_6_0;
	  tre0_0_0 = c_re(inout[0]);
	  tim0_0_0 = c_im(inout[0]);
	  {
	       FFTW_REAL tr;
	       FFTW_REAL ti;
	       FFTW_REAL twr;
	       FFTW_REAL twi;
	       tr = c_re(inout[stride]);
	       ti = c_im(inout[stride]);
	       twr = c_re(W[0]);
	       twi = c_im(W[0]);
	       tre0_1_0 = (tr * twr) - (ti * twi);
	       tim0_1_0 = (tr * twi) + (ti * twr);
	  }
	  {
	       FFTW_REAL tr;
	       FFTW_REAL ti;
	       FFTW_REAL twr;
	       FFTW_REAL twi;
	       tr = c_re(inout[2 * stride]);
	       ti = c_im(inout[2 * stride]);
	       twr = c_re(W[1]);
	       twi = c_im(W[1]);
	       tre0_2_0 = (tr * twr) - (ti * twi);
	       tim0_2_0 = (tr * twi) + (ti * twr);
	  }
	  {
	       FFTW_REAL tr;
	       FFTW_REAL ti;
	       FFTW_REAL twr;
	       FFTW_REAL twi;
	       tr = c_re(inout[3 * stride]);
	       ti = c_im(inout[3 * stride]);
	       twr = c_re(W[2]);
	       twi = c_im(W[2]);
	       tre0_3_0 = (tr * twr) - (ti * twi);
	       tim0_3_0 = (tr * twi) + (ti * twr);
	  }
	  {
	       FFTW_REAL tr;
	       FFTW_REAL ti;
	       FFTW_REAL twr;
	       FFTW_REAL twi;
	       tr = c_re(inout[4 * stride]);
	       ti = c_im(inout[4 * stride]);
	       twr = c_re(W[3]);
	       twi = c_im(W[3]);
	       tre0_4_0 = (tr * twr) - (ti * twi);
	       tim0_4_0 = (tr * twi) + (ti * twr);
	  }
	  {
	       FFTW_REAL tr;
	       FFTW_REAL ti;
	       FFTW_REAL twr;
	       FFTW_REAL twi;
	       tr = c_re(inout[5 * stride]);
	       ti = c_im(inout[5 * stride]);
	       twr = c_re(W[4]);
	       twi = c_im(W[4]);
	       tre0_5_0 = (tr * twr) - (ti * twi);
	       tim0_5_0 = (tr * twi) + (ti * twr);
	  }
	  {
	       FFTW_REAL tr;
	       FFTW_REAL ti;
	       FFTW_REAL twr;
	       FFTW_REAL twi;
	       tr = c_re(inout[6 * stride]);
	       ti = c_im(inout[6 * stride]);
	       twr = c_re(W[5]);
	       twi = c_im(W[5]);
	       tre0_6_0 = (tr * twr) - (ti * twi);
	       tim0_6_0 = (tr * twi) + (ti * twr);
	  }
	  c_re(inout[0]) = tre0_0_0 + tre0_1_0 + tre0_2_0 + tre0_3_0 + tre0_4_0 + tre0_5_0 + tre0_6_0;
	  c_im(inout[0]) = tim0_0_0 + tim0_1_0 + tim0_2_0 + tim0_3_0 + tim0_4_0 + tim0_5_0 + tim0_6_0;
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tre1_1_0;
	       tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tre0_1_0 + tre0_6_0)) - (((FFTW_REAL) FFTW_K900968867) * (tre0_3_0 + tre0_4_0)) - (((FFTW_REAL) FFTW_K222520933) * (tre0_2_0 + tre0_5_0));
	       tre1_1_0 = (((FFTW_REAL) FFTW_K781831482) * (tim0_1_0 - tim0_6_0)) + (((FFTW_REAL) FFTW_K974927912) * (tim0_2_0 - tim0_5_0)) + (((FFTW_REAL) FFTW_K433883739) * (tim0_3_0 - tim0_4_0));
	       c_re(inout[stride]) = tre1_0_0 + tre1_1_0;
	       c_re(inout[6 * stride]) = tre1_0_0 - tre1_1_0;
	  }
	  {
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tim1_1_0;
	       tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tim0_1_0 + tim0_6_0)) - (((FFTW_REAL) FFTW_K900968867) * (tim0_3_0 + tim0_4_0)) - (((FFTW_REAL) FFTW_K222520933) * (tim0_2_0 + tim0_5_0));
	       tim1_1_0 = (((FFTW_REAL) FFTW_K781831482) * (tre0_6_0 - tre0_1_0)) + (((FFTW_REAL) FFTW_K974927912) * (tre0_5_0 - tre0_2_0)) + (((FFTW_REAL) FFTW_K433883739) * (tre0_4_0 - tre0_3_0));
	       c_im(inout[stride]) = tim1_0_0 + tim1_1_0;
	       c_im(inout[6 * stride]) = tim1_0_0 - tim1_1_0;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tre1_1_0;
	       tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tre0_3_0 + tre0_4_0)) - (((FFTW_REAL) FFTW_K900968867) * (tre0_2_0 + tre0_5_0)) - (((FFTW_REAL) FFTW_K222520933) * (tre0_1_0 + tre0_6_0));
	       tre1_1_0 = (((FFTW_REAL) FFTW_K974927912) * (tim0_1_0 - tim0_6_0)) + (((FFTW_REAL) FFTW_K433883739) * (tim0_5_0 - tim0_2_0)) + (((FFTW_REAL) FFTW_K781831482) * (tim0_4_0 - tim0_3_0));
	       c_re(inout[2 * stride]) = tre1_0_0 + tre1_1_0;
	       c_re(inout[5 * stride]) = tre1_0_0 - tre1_1_0;
	  }
	  {
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tim1_1_0;
	       tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tim0_3_0 + tim0_4_0)) - (((FFTW_REAL) FFTW_K900968867) * (tim0_2_0 + tim0_5_0)) - (((FFTW_REAL) FFTW_K222520933) * (tim0_1_0 + tim0_6_0));
	       tim1_1_0 = (((FFTW_REAL) FFTW_K974927912) * (tre0_6_0 - tre0_1_0)) + (((FFTW_REAL) FFTW_K433883739) * (tre0_2_0 - tre0_5_0)) + (((FFTW_REAL) FFTW_K781831482) * (tre0_3_0 - tre0_4_0));
	       c_im(inout[2 * stride]) = tim1_0_0 + tim1_1_0;
	       c_im(inout[5 * stride]) = tim1_0_0 - tim1_1_0;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tre1_1_0;
	       tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tre0_2_0 + tre0_5_0)) - (((FFTW_REAL) FFTW_K222520933) * (tre0_3_0 + tre0_4_0)) - (((FFTW_REAL) FFTW_K900968867) * (tre0_1_0 + tre0_6_0));
	       tre1_1_0 = (((FFTW_REAL) FFTW_K433883739) * (tim0_1_0 - tim0_6_0)) + (((FFTW_REAL) FFTW_K781831482) * (tim0_5_0 - tim0_2_0)) + (((FFTW_REAL) FFTW_K974927912) * (tim0_3_0 - tim0_4_0));
	       c_re(inout[3 * stride]) = tre1_0_0 + tre1_1_0;
	       c_re(inout[4 * stride]) = tre1_0_0 - tre1_1_0;
	  }
	  {
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tim1_1_0;
	       tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tim0_2_0 + tim0_5_0)) - (((FFTW_REAL) FFTW_K222520933) * (tim0_3_0 + tim0_4_0)) - (((FFTW_REAL) FFTW_K900968867) * (tim0_1_0 + tim0_6_0));
	       tim1_1_0 = (((FFTW_REAL) FFTW_K433883739) * (tre0_6_0 - tre0_1_0)) + (((FFTW_REAL) FFTW_K781831482) * (tre0_2_0 - tre0_5_0)) + (((FFTW_REAL) FFTW_K974927912) * (tre0_4_0 - tre0_3_0));
	       c_im(inout[3 * stride]) = tim1_0_0 + tim1_1_0;
	       c_im(inout[4 * stride]) = tim1_0_0 - tim1_1_0;
	  }
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 66 FP additions and 32 FP multiplications */

void fftw_twiddle_8(FFTW_COMPLEX *A, const FFTW_COMPLEX *W, int stride, int m, int dist)
{
     int i;
     COMPLEX *inout;
     inout = A;
     for (i = 0; i < m; i = i + 1, inout = inout + dist, W = W + 7) {
	  FFTW_REAL tre0_0_0;
	  FFTW_REAL tim0_0_0;
	  FFTW_REAL tre0_0_1;
	  FFTW_REAL tim0_0_1;
	  FFTW_REAL tre0_0_2;
	  FFTW_REAL tim0_0_2;
	  FFTW_REAL tre0_0_3;
	  FFTW_REAL tim0_0_3;
	  FFTW_REAL tre0_1_0;
	  FFTW_REAL tim0_1_0;
	  FFTW_REAL tre0_1_1;
	  FFTW_REAL tim0_1_1;
	  FFTW_REAL tre0_1_2;
	  FFTW_REAL tim0_1_2;
	  FFTW_REAL tre0_1_3;
	  FFTW_REAL tim0_1_3;
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       tre1_0_0 = c_re(inout[0]);
	       tim1_0_0 = c_im(inout[0]);
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[4 * stride]);
		    ti = c_im(inout[4 * stride]);
		    twr = c_re(W[3]);
		    twi = c_im(W[3]);
		    tre1_1_0 = (tr * twr) - (ti * twi);
		    tim1_1_0 = (tr * twi) + (ti * twr);
	       }
	       tre0_0_0 = tre1_0_0 + tre1_1_0;
	       tim0_0_0 = tim1_0_0 + tim1_1_0;
	       tre0_1_0 = tre1_0_0 - tre1_1_0;
	       tim0_1_0 = tim1_0_0 - tim1_1_0;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[stride]);
		    ti = c_im(inout[stride]);
		    twr = c_re(W[0]);
		    twi = c_im(W[0]);
		    tre1_0_0 = (tr * twr) - (ti * twi);
		    tim1_0_0 = (tr * twi) + (ti * twr);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[5 * stride]);
		    ti = c_im(inout[5 * stride]);
		    twr = c_re(W[4]);
		    twi = c_im(W[4]);
		    tre1_1_0 = (tr * twr) - (ti * twi);
		    tim1_1_0 = (tr * twi) + (ti * twr);
	       }
	       tre0_0_1 = tre1_0_0 + tre1_1_0;
	       tim0_0_1 = tim1_0_0 + tim1_1_0;
	       tre0_1_1 = tre1_0_0 - tre1_1_0;
	       tim0_1_1 = tim1_0_0 - tim1_1_0;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[2 * stride]);
		    ti = c_im(inout[2 * stride]);
		    twr = c_re(W[1]);
		    twi = c_im(W[1]);
		    tre1_0_0 = (tr * twr) - (ti * twi);
		    tim1_0_0 = (tr * twi) + (ti * twr);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[6 * stride]);
		    ti = c_im(inout[6 * stride]);
		    twr = c_re(W[5]);
		    twi = c_im(W[5]);
		    tre1_1_0 = (tr * twr) - (ti * twi);
		    tim1_1_0 = (tr * twi) + (ti * twr);
	       }
	       tre0_0_2 = tre1_0_0 + tre1_1_0;
	       tim0_0_2 = tim1_0_0 + tim1_1_0;
	       tre0_1_2 = tre1_0_0 - tre1_1_0;
	       tim0_1_2 = tim1_0_0 - tim1_1_0;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[3 * stride]);
		    ti = c_im(inout[3 * stride]);
		    twr = c_re(W[2]);
		    twi = c_im(W[2]);
		    tre1_0_0 = (tr * twr) - (ti * twi);
		    tim1_0_0 = (tr * twi) + (ti * twr);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[7 * stride]);
		    ti = c_im(inout[7 * stride]);
		    twr = c_re(W[6]);
		    twi = c_im(W[6]);
		    tre1_1_0 = (tr * twr) - (ti * twi);
		    tim1_1_0 = (tr * twi) + (ti * twr);
	       }
	       tre0_0_3 = tre1_0_0 + tre1_1_0;
	       tim0_0_3 = tim1_0_0 + tim1_1_0;
	       tre0_1_3 = tre1_0_0 - tre1_1_0;
	       tim0_1_3 = tim1_0_0 - tim1_1_0;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       tre1_0_0 = tre0_0_0 + tre0_0_2;
	       tim1_0_0 = tim0_0_0 + tim0_0_2;
	       tre1_1_0 = tre0_0_0 - tre0_0_2;
	       tim1_1_0 = tim0_0_0 - tim0_0_2;
	       tre1_0_1 = tre0_0_1 + tre0_0_3;
	       tim1_0_1 = tim0_0_1 + tim0_0_3;
	       tre1_1_1 = tre0_0_1 - tre0_0_3;
	       tim1_1_1 = tim0_0_1 - tim0_0_3;
	       c_re(inout[0]) = tre1_0_0 + tre1_0_1;
	       c_im(inout[0]) = tim1_0_0 + tim1_0_1;
	       c_re(inout[4 * stride]) = tre1_0_0 - tre1_0_1;
	       c_im(inout[4 * stride]) = tim1_0_0 - tim1_0_1;
	       c_re(inout[2 * stride]) = tre1_1_0 + tim1_1_1;
	       c_im(inout[2 * stride]) = tim1_1_0 - tre1_1_1;
	       c_re(inout[6 * stride]) = tre1_1_0 - tim1_1_1;
	       c_im(inout[6 * stride]) = tim1_1_0 + tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       tre1_0_0 = tre0_1_0 + tim0_1_2;
	       tim1_0_0 = tim0_1_0 - tre0_1_2;
	       tre1_1_0 = tre0_1_0 - tim0_1_2;
	       tim1_1_0 = tim0_1_0 + tre0_1_2;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_1_1 + tim0_1_1);
		    tim2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_1_1 - tre0_1_1);
		    tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_1_3 - tre0_1_3);
		    tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_1_3 + tre0_1_3);
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 - tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 + tim2_1_0;
	       }
	       c_re(inout[stride]) = tre1_0_0 + tre1_0_1;
	       c_im(inout[stride]) = tim1_0_0 + tim1_0_1;
	       c_re(inout[5 * stride]) = tre1_0_0 - tre1_0_1;
	       c_im(inout[5 * stride]) = tim1_0_0 - tim1_0_1;
	       c_re(inout[3 * stride]) = tre1_1_0 + tim1_1_1;
	       c_im(inout[3 * stride]) = tim1_1_0 - tre1_1_1;
	       c_re(inout[7 * stride]) = tre1_1_0 - tim1_1_1;
	       c_im(inout[7 * stride]) = tim1_1_0 + tre1_1_1;
	  }
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 108 FP additions and 72 FP multiplications */

void fftw_twiddle_9(FFTW_COMPLEX *A, const FFTW_COMPLEX *W, int stride, int m, int dist)
{
     int i;
     COMPLEX *inout;
     inout = A;
     for (i = 0; i < m; i = i + 1, inout = inout + dist, W = W + 8) {
	  FFTW_REAL tre0_0_0;
	  FFTW_REAL tim0_0_0;
	  FFTW_REAL tre0_0_1;
	  FFTW_REAL tim0_0_1;
	  FFTW_REAL tre0_0_2;
	  FFTW_REAL tim0_0_2;
	  FFTW_REAL tre0_1_0;
	  FFTW_REAL tim0_1_0;
	  FFTW_REAL tre0_1_1;
	  FFTW_REAL tim0_1_1;
	  FFTW_REAL tre0_1_2;
	  FFTW_REAL tim0_1_2;
	  FFTW_REAL tre0_2_0;
	  FFTW_REAL tim0_2_0;
	  FFTW_REAL tre0_2_1;
	  FFTW_REAL tim0_2_1;
	  FFTW_REAL tre0_2_2;
	  FFTW_REAL tim0_2_2;
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_2_0;
	       FFTW_REAL tim1_2_0;
	       tre1_0_0 = c_re(inout[0]);
	       tim1_0_0 = c_im(inout[0]);
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[3 * stride]);
		    ti = c_im(inout[3 * stride]);
		    twr = c_re(W[2]);
		    twi = c_im(W[2]);
		    tre1_1_0 = (tr * twr) - (ti * twi);
		    tim1_1_0 = (tr * twi) + (ti * twr);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[6 * stride]);
		    ti = c_im(inout[6 * stride]);
		    twr = c_re(W[5]);
		    twi = c_im(W[5]);
		    tre1_2_0 = (tr * twr) - (ti * twi);
		    tim1_2_0 = (tr * twi) + (ti * twr);
	       }
	       tre0_0_0 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	       tim0_0_0 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tre2_1_0;
		    tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
		    tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_1_0 - tim1_2_0);
		    tre0_1_0 = tre2_0_0 + tre2_1_0;
		    tre0_2_0 = tre2_0_0 - tre2_1_0;
	       }
	       {
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tim2_1_0;
		    tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
		    tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_2_0 - tre1_1_0);
		    tim0_1_0 = tim2_0_0 + tim2_1_0;
		    tim0_2_0 = tim2_0_0 - tim2_1_0;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_2_0;
	       FFTW_REAL tim1_2_0;
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[stride]);
		    ti = c_im(inout[stride]);
		    twr = c_re(W[0]);
		    twi = c_im(W[0]);
		    tre1_0_0 = (tr * twr) - (ti * twi);
		    tim1_0_0 = (tr * twi) + (ti * twr);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[4 * stride]);
		    ti = c_im(inout[4 * stride]);
		    twr = c_re(W[3]);
		    twi = c_im(W[3]);
		    tre1_1_0 = (tr * twr) - (ti * twi);
		    tim1_1_0 = (tr * twi) + (ti * twr);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[7 * stride]);
		    ti = c_im(inout[7 * stride]);
		    twr = c_re(W[6]);
		    twi = c_im(W[6]);
		    tre1_2_0 = (tr * twr) - (ti * twi);
		    tim1_2_0 = (tr * twi) + (ti * twr);
	       }
	       tre0_0_1 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	       tim0_0_1 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tre2_1_0;
		    tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
		    tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_1_0 - tim1_2_0);
		    tre0_1_1 = tre2_0_0 + tre2_1_0;
		    tre0_2_1 = tre2_0_0 - tre2_1_0;
	       }
	       {
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tim2_1_0;
		    tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
		    tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_2_0 - tre1_1_0);
		    tim0_1_1 = tim2_0_0 + tim2_1_0;
		    tim0_2_1 = tim2_0_0 - tim2_1_0;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_2_0;
	       FFTW_REAL tim1_2_0;
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[2 * stride]);
		    ti = c_im(inout[2 * stride]);
		    twr = c_re(W[1]);
		    twi = c_im(W[1]);
		    tre1_0_0 = (tr * twr) - (ti * twi);
		    tim1_0_0 = (tr * twi) + (ti * twr);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[5 * stride]);
		    ti = c_im(inout[5 * stride]);
		    twr = c_re(W[4]);
		    twi = c_im(W[4]);
		    tre1_1_0 = (tr * twr) - (ti * twi);
		    tim1_1_0 = (tr * twi) + (ti * twr);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[8 * stride]);
		    ti = c_im(inout[8 * stride]);
		    twr = c_re(W[7]);
		    twi = c_im(W[7]);
		    tre1_2_0 = (tr * twr) - (ti * twi);
		    tim1_2_0 = (tr * twi) + (ti * twr);
	       }
	       tre0_0_2 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	       tim0_0_2 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tre2_1_0;
		    tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
		    tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_1_0 - tim1_2_0);
		    tre0_1_2 = tre2_0_0 + tre2_1_0;
		    tre0_2_2 = tre2_0_0 - tre2_1_0;
	       }
	       {
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tim2_1_0;
		    tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
		    tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_2_0 - tre1_1_0);
		    tim0_1_2 = tim2_0_0 + tim2_1_0;
		    tim0_2_2 = tim2_0_0 - tim2_1_0;
	       }
	  }
	  c_re(inout[0]) = tre0_0_0 + tre0_0_1 + tre0_0_2;
	  c_im(inout[0]) = tim0_0_0 + tim0_0_1 + tim0_0_2;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre0_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre0_0_1 + tre0_0_2));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim0_0_1 - tim0_0_2);
	       c_re(inout[3 * stride]) = tre2_0_0 + tre2_1_0;
	       c_re(inout[6 * stride]) = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim0_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim0_0_1 + tim0_0_2));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre0_0_2 - tre0_0_1);
	       c_im(inout[3 * stride]) = tim2_0_0 + tim2_1_0;
	       c_im(inout[6 * stride]) = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_2_0;
	       FFTW_REAL tim1_2_0;
	       tre1_1_0 = (((FFTW_REAL) FFTW_K766044443) * tre0_1_1) + (((FFTW_REAL) FFTW_K642787609) * tim0_1_1);
	       tim1_1_0 = (((FFTW_REAL) FFTW_K766044443) * tim0_1_1) - (((FFTW_REAL) FFTW_K642787609) * tre0_1_1);
	       tre1_2_0 = (((FFTW_REAL) FFTW_K173648177) * tre0_1_2) + (((FFTW_REAL) FFTW_K984807753) * tim0_1_2);
	       tim1_2_0 = (((FFTW_REAL) FFTW_K173648177) * tim0_1_2) - (((FFTW_REAL) FFTW_K984807753) * tre0_1_2);
	       c_re(inout[stride]) = tre0_1_0 + tre1_1_0 + tre1_2_0;
	       c_im(inout[stride]) = tim0_1_0 + tim1_1_0 + tim1_2_0;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tre2_1_0;
		    tre2_0_0 = tre0_1_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
		    tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_1_0 - tim1_2_0);
		    c_re(inout[4 * stride]) = tre2_0_0 + tre2_1_0;
		    c_re(inout[7 * stride]) = tre2_0_0 - tre2_1_0;
	       }
	       {
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tim2_1_0;
		    tim2_0_0 = tim0_1_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
		    tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_2_0 - tre1_1_0);
		    c_im(inout[4 * stride]) = tim2_0_0 + tim2_1_0;
		    c_im(inout[7 * stride]) = tim2_0_0 - tim2_1_0;
	       }
	  }
	  {
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_2_0;
	       FFTW_REAL tim1_2_0;
	       tre1_1_0 = (((FFTW_REAL) FFTW_K173648177) * tre0_2_1) + (((FFTW_REAL) FFTW_K984807753) * tim0_2_1);
	       tim1_1_0 = (((FFTW_REAL) FFTW_K173648177) * tim0_2_1) - (((FFTW_REAL) FFTW_K984807753) * tre0_2_1);
	       tre1_2_0 = (((FFTW_REAL) FFTW_K342020143) * tim0_2_2) - (((FFTW_REAL) FFTW_K939692620) * tre0_2_2);
	       tim1_2_0 = (((FFTW_REAL) FFTW_K939692620) * tim0_2_2) + (((FFTW_REAL) FFTW_K342020143) * tre0_2_2);
	       c_re(inout[2 * stride]) = tre0_2_0 + tre1_1_0 + tre1_2_0;
	       c_im(inout[2 * stride]) = tim0_2_0 + tim1_1_0 - tim1_2_0;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tre2_1_0;
		    tre2_0_0 = tre0_2_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
		    tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_1_0 + tim1_2_0);
		    c_re(inout[5 * stride]) = tre2_0_0 + tre2_1_0;
		    c_re(inout[8 * stride]) = tre2_0_0 - tre2_1_0;
	       }
	       {
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tim2_1_0;
		    tim2_0_0 = tim0_2_0 + (((FFTW_REAL) FFTW_K499999999) * (tim1_2_0 - tim1_1_0));
		    tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_2_0 - tre1_1_0);
		    c_im(inout[5 * stride]) = tim2_0_0 + tim2_1_0;
		    c_im(inout[8 * stride]) = tim2_0_0 - tim2_1_0;
	       }
	  }
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 126 FP additions and 68 FP multiplications */

void fftwi_twiddle_10(FFTW_COMPLEX *A, const FFTW_COMPLEX *W, int stride, int m, int dist)
{
     int i;
     COMPLEX *inout;
     inout = A;
     for (i = 0; i < m; i = i + 1, inout = inout + dist, W = W + 9) {
	  FFTW_REAL tre0_0_0;
	  FFTW_REAL tim0_0_0;
	  FFTW_REAL tre0_0_1;
	  FFTW_REAL tim0_0_1;
	  FFTW_REAL tre0_0_2;
	  FFTW_REAL tim0_0_2;
	  FFTW_REAL tre0_0_3;
	  FFTW_REAL tim0_0_3;
	  FFTW_REAL tre0_0_4;
	  FFTW_REAL tim0_0_4;
	  FFTW_REAL tre0_1_0;
	  FFTW_REAL tim0_1_0;
	  FFTW_REAL tre0_1_1;
	  FFTW_REAL tim0_1_1;
	  FFTW_REAL tre0_1_2;
	  FFTW_REAL tim0_1_2;
	  FFTW_REAL tre0_1_3;
	  FFTW_REAL tim0_1_3;
	  FFTW_REAL tre0_1_4;
	  FFTW_REAL tim0_1_4;
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       tre1_0_0 = c_re(inout[0]);
	       tim1_0_0 = c_im(inout[0]);
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[5 * stride]);
		    ti = c_im(inout[5 * stride]);
		    twr = c_re(W[4]);
		    twi = c_im(W[4]);
		    tre1_1_0 = (tr * twr) + (ti * twi);
		    tim1_1_0 = (ti * twr) - (tr * twi);
	       }
	       tre0_0_0 = tre1_0_0 + tre1_1_0;
	       tim0_0_0 = tim1_0_0 + tim1_1_0;
	       tre0_1_0 = tre1_0_0 - tre1_1_0;
	       tim0_1_0 = tim1_0_0 - tim1_1_0;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[2 * stride]);
		    ti = c_im(inout[2 * stride]);
		    twr = c_re(W[1]);
		    twi = c_im(W[1]);
		    tre1_0_0 = (tr * twr) + (ti * twi);
		    tim1_0_0 = (ti * twr) - (tr * twi);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[7 * stride]);
		    ti = c_im(inout[7 * stride]);
		    twr = c_re(W[6]);
		    twi = c_im(W[6]);
		    tre1_1_0 = (tr * twr) + (ti * twi);
		    tim1_1_0 = (ti * twr) - (tr * twi);
	       }
	       tre0_0_1 = tre1_0_0 + tre1_1_0;
	       tim0_0_1 = tim1_0_0 + tim1_1_0;
	       tre0_1_1 = tre1_0_0 - tre1_1_0;
	       tim0_1_1 = tim1_0_0 - tim1_1_0;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[4 * stride]);
		    ti = c_im(inout[4 * stride]);
		    twr = c_re(W[3]);
		    twi = c_im(W[3]);
		    tre1_0_0 = (tr * twr) + (ti * twi);
		    tim1_0_0 = (ti * twr) - (tr * twi);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[9 * stride]);
		    ti = c_im(inout[9 * stride]);
		    twr = c_re(W[8]);
		    twi = c_im(W[8]);
		    tre1_1_0 = (tr * twr) + (ti * twi);
		    tim1_1_0 = (ti * twr) - (tr * twi);
	       }
	       tre0_0_2 = tre1_0_0 + tre1_1_0;
	       tim0_0_2 = tim1_0_0 + tim1_1_0;
	       tre0_1_2 = tre1_0_0 - tre1_1_0;
	       tim0_1_2 = tim1_0_0 - tim1_1_0;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[6 * stride]);
		    ti = c_im(inout[6 * stride]);
		    twr = c_re(W[5]);
		    twi = c_im(W[5]);
		    tre1_0_0 = (tr * twr) + (ti * twi);
		    tim1_0_0 = (ti * twr) - (tr * twi);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[stride]);
		    ti = c_im(inout[stride]);
		    twr = c_re(W[0]);
		    twi = c_im(W[0]);
		    tre1_1_0 = (tr * twr) + (ti * twi);
		    tim1_1_0 = (ti * twr) - (tr * twi);
	       }
	       tre0_0_3 = tre1_0_0 + tre1_1_0;
	       tim0_0_3 = tim1_0_0 + tim1_1_0;
	       tre0_1_3 = tre1_0_0 - tre1_1_0;
	       tim0_1_3 = tim1_0_0 - tim1_1_0;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[8 * stride]);
		    ti = c_im(inout[8 * stride]);
		    twr = c_re(W[7]);
		    twi = c_im(W[7]);
		    tre1_0_0 = (tr * twr) + (ti * twi);
		    tim1_0_0 = (ti * twr) - (tr * twi);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[3 * stride]);
		    ti = c_im(inout[3 * stride]);
		    twr = c_re(W[2]);
		    twi = c_im(W[2]);
		    tre1_1_0 = (tr * twr) + (ti * twi);
		    tim1_1_0 = (ti * twr) - (tr * twi);
	       }
	       tre0_0_4 = tre1_0_0 + tre1_1_0;
	       tim0_0_4 = tim1_0_0 + tim1_1_0;
	       tre0_1_4 = tre1_0_0 - tre1_1_0;
	       tim0_1_4 = tim1_0_0 - tim1_1_0;
	  }
	  c_re(inout[0]) = tre0_0_0 + tre0_0_1 + tre0_0_2 + tre0_0_3 + tre0_0_4;
	  c_im(inout[0]) = tim0_0_0 + tim0_0_1 + tim0_0_2 + tim0_0_3 + tim0_0_4;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_0_1 + tre0_0_4)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_0_2 + tre0_0_3));
	       tre2_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tim0_0_4 - tim0_0_1)) + (((FFTW_REAL) FFTW_K587785252) * (tim0_0_3 - tim0_0_2));
	       c_re(inout[6 * stride]) = tre2_0_0 + tre2_1_0;
	       c_re(inout[4 * stride]) = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_0_1 + tim0_0_4)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_0_2 + tim0_0_3));
	       tim2_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tre0_0_1 - tre0_0_4)) + (((FFTW_REAL) FFTW_K587785252) * (tre0_0_2 - tre0_0_3));
	       c_im(inout[6 * stride]) = tim2_0_0 + tim2_1_0;
	       c_im(inout[4 * stride]) = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_0_2 + tre0_0_3)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_0_1 + tre0_0_4));
	       tre2_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tim0_0_4 - tim0_0_1)) + (((FFTW_REAL) FFTW_K951056516) * (tim0_0_2 - tim0_0_3));
	       c_re(inout[2 * stride]) = tre2_0_0 + tre2_1_0;
	       c_re(inout[8 * stride]) = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_0_2 + tim0_0_3)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_0_1 + tim0_0_4));
	       tim2_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tre0_0_1 - tre0_0_4)) + (((FFTW_REAL) FFTW_K951056516) * (tre0_0_3 - tre0_0_2));
	       c_im(inout[2 * stride]) = tim2_0_0 + tim2_1_0;
	       c_im(inout[8 * stride]) = tim2_0_0 - tim2_1_0;
	  }
	  c_re(inout[5 * stride]) = tre0_1_0 + tre0_1_1 + tre0_1_2 + tre0_1_3 + tre0_1_4;
	  c_im(inout[5 * stride]) = tim0_1_0 + tim0_1_1 + tim0_1_2 + tim0_1_3 + tim0_1_4;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre0_1_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_1_1 + tre0_1_4)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_1_2 + tre0_1_3));
	       tre2_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tim0_1_4 - tim0_1_1)) + (((FFTW_REAL) FFTW_K587785252) * (tim0_1_3 - tim0_1_2));
	       c_re(inout[stride]) = tre2_0_0 + tre2_1_0;
	       c_re(inout[9 * stride]) = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim0_1_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_1_1 + tim0_1_4)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_1_2 + tim0_1_3));
	       tim2_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tre0_1_1 - tre0_1_4)) + (((FFTW_REAL) FFTW_K587785252) * (tre0_1_2 - tre0_1_3));
	       c_im(inout[stride]) = tim2_0_0 + tim2_1_0;
	       c_im(inout[9 * stride]) = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre0_1_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_1_2 + tre0_1_3)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_1_1 + tre0_1_4));
	       tre2_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tim0_1_4 - tim0_1_1)) + (((FFTW_REAL) FFTW_K951056516) * (tim0_1_2 - tim0_1_3));
	       c_re(inout[7 * stride]) = tre2_0_0 + tre2_1_0;
	       c_re(inout[3 * stride]) = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim0_1_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_1_2 + tim0_1_3)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_1_1 + tim0_1_4));
	       tim2_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tre0_1_1 - tre0_1_4)) + (((FFTW_REAL) FFTW_K951056516) * (tre0_1_3 - tre0_1_2));
	       c_im(inout[7 * stride]) = tim2_0_0 + tim2_1_0;
	       c_im(inout[3 * stride]) = tim2_0_0 - tim2_1_0;
	  }
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 174 FP additions and 84 FP multiplications */

void fftwi_twiddle_16(FFTW_COMPLEX *A, const FFTW_COMPLEX *W, int stride, int m, int dist)
{
     int i;
     COMPLEX *inout;
     inout = A;
     for (i = 0; i < m; i = i + 1, inout = inout + dist, W = W + 15) {
	  FFTW_REAL tre0_0_0;
	  FFTW_REAL tim0_0_0;
	  FFTW_REAL tre0_0_1;
	  FFTW_REAL tim0_0_1;
	  FFTW_REAL tre0_0_2;
	  FFTW_REAL tim0_0_2;
	  FFTW_REAL tre0_0_3;
	  FFTW_REAL tim0_0_3;
	  FFTW_REAL tre0_1_0;
	  FFTW_REAL tim0_1_0;
	  FFTW_REAL tre0_1_1;
	  FFTW_REAL tim0_1_1;
	  FFTW_REAL tre0_1_2;
	  FFTW_REAL tim0_1_2;
	  FFTW_REAL tre0_1_3;
	  FFTW_REAL tim0_1_3;
	  FFTW_REAL tre0_2_0;
	  FFTW_REAL tim0_2_0;
	  FFTW_REAL tre0_2_1;
	  FFTW_REAL tim0_2_1;
	  FFTW_REAL tre0_2_2;
	  FFTW_REAL tim0_2_2;
	  FFTW_REAL tre0_2_3;
	  FFTW_REAL tim0_2_3;
	  FFTW_REAL tre0_3_0;
	  FFTW_REAL tim0_3_0;
	  FFTW_REAL tre0_3_1;
	  FFTW_REAL tim0_3_1;
	  FFTW_REAL tre0_3_2;
	  FFTW_REAL tim0_3_2;
	  FFTW_REAL tre0_3_3;
	  FFTW_REAL tim0_3_3;
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = c_re(inout[0]);
		    tim2_0_0 = c_im(inout[0]);
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[8 * stride]);
			 ti = c_im(inout[8 * stride]);
			 twr = c_re(W[7]);
			 twi = c_im(W[7]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[4 * stride]);
			 ti = c_im(inout[4 * stride]);
			 twr = c_re(W[3]);
			 twi = c_im(W[3]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[12 * stride]);
			 ti = c_im(inout[12 * stride]);
			 twr = c_re(W[11]);
			 twi = c_im(W[11]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       tre0_0_0 = tre1_0_0 + tre1_0_1;
	       tim0_0_0 = tim1_0_0 + tim1_0_1;
	       tre0_2_0 = tre1_0_0 - tre1_0_1;
	       tim0_2_0 = tim1_0_0 - tim1_0_1;
	       tre0_1_0 = tre1_1_0 - tim1_1_1;
	       tim0_1_0 = tim1_1_0 + tre1_1_1;
	       tre0_3_0 = tre1_1_0 + tim1_1_1;
	       tim0_3_0 = tim1_1_0 - tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[stride]);
			 ti = c_im(inout[stride]);
			 twr = c_re(W[0]);
			 twi = c_im(W[0]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[9 * stride]);
			 ti = c_im(inout[9 * stride]);
			 twr = c_re(W[8]);
			 twi = c_im(W[8]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[5 * stride]);
			 ti = c_im(inout[5 * stride]);
			 twr = c_re(W[4]);
			 twi = c_im(W[4]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[13 * stride]);
			 ti = c_im(inout[13 * stride]);
			 twr = c_re(W[12]);
			 twi = c_im(W[12]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       tre0_0_1 = tre1_0_0 + tre1_0_1;
	       tim0_0_1 = tim1_0_0 + tim1_0_1;
	       tre0_2_1 = tre1_0_0 - tre1_0_1;
	       tim0_2_1 = tim1_0_0 - tim1_0_1;
	       tre0_1_1 = tre1_1_0 - tim1_1_1;
	       tim0_1_1 = tim1_1_0 + tre1_1_1;
	       tre0_3_1 = tre1_1_0 + tim1_1_1;
	       tim0_3_1 = tim1_1_0 - tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[2 * stride]);
			 ti = c_im(inout[2 * stride]);
			 twr = c_re(W[1]);
			 twi = c_im(W[1]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[10 * stride]);
			 ti = c_im(inout[10 * stride]);
			 twr = c_re(W[9]);
			 twi = c_im(W[9]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[6 * stride]);
			 ti = c_im(inout[6 * stride]);
			 twr = c_re(W[5]);
			 twi = c_im(W[5]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[14 * stride]);
			 ti = c_im(inout[14 * stride]);
			 twr = c_re(W[13]);
			 twi = c_im(W[13]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       tre0_0_2 = tre1_0_0 + tre1_0_1;
	       tim0_0_2 = tim1_0_0 + tim1_0_1;
	       tre0_2_2 = tre1_0_0 - tre1_0_1;
	       tim0_2_2 = tim1_0_0 - tim1_0_1;
	       tre0_1_2 = tre1_1_0 - tim1_1_1;
	       tim0_1_2 = tim1_1_0 + tre1_1_1;
	       tre0_3_2 = tre1_1_0 + tim1_1_1;
	       tim0_3_2 = tim1_1_0 - tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[3 * stride]);
			 ti = c_im(inout[3 * stride]);
			 twr = c_re(W[2]);
			 twi = c_im(W[2]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[11 * stride]);
			 ti = c_im(inout[11 * stride]);
			 twr = c_re(W[10]);
			 twi = c_im(W[10]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[7 * stride]);
			 ti = c_im(inout[7 * stride]);
			 twr = c_re(W[6]);
			 twi = c_im(W[6]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[15 * stride]);
			 ti = c_im(inout[15 * stride]);
			 twr = c_re(W[14]);
			 twi = c_im(W[14]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       tre0_0_3 = tre1_0_0 + tre1_0_1;
	       tim0_0_3 = tim1_0_0 + tim1_0_1;
	       tre0_2_3 = tre1_0_0 - tre1_0_1;
	       tim0_2_3 = tim1_0_0 - tim1_0_1;
	       tre0_1_3 = tre1_1_0 - tim1_1_1;
	       tim0_1_3 = tim1_1_0 + tre1_1_1;
	       tre0_3_3 = tre1_1_0 + tim1_1_1;
	       tim0_3_3 = tim1_1_0 - tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       tre1_0_0 = tre0_0_0 + tre0_0_2;
	       tim1_0_0 = tim0_0_0 + tim0_0_2;
	       tre1_1_0 = tre0_0_0 - tre0_0_2;
	       tim1_1_0 = tim0_0_0 - tim0_0_2;
	       tre1_0_1 = tre0_0_1 + tre0_0_3;
	       tim1_0_1 = tim0_0_1 + tim0_0_3;
	       tre1_1_1 = tre0_0_1 - tre0_0_3;
	       tim1_1_1 = tim0_0_1 - tim0_0_3;
	       c_re(inout[0]) = tre1_0_0 + tre1_0_1;
	       c_im(inout[0]) = tim1_0_0 + tim1_0_1;
	       c_re(inout[8 * stride]) = tre1_0_0 - tre1_0_1;
	       c_im(inout[8 * stride]) = tim1_0_0 - tim1_0_1;
	       c_re(inout[4 * stride]) = tre1_1_0 - tim1_1_1;
	       c_im(inout[4 * stride]) = tim1_1_0 + tre1_1_1;
	       c_re(inout[12 * stride]) = tre1_1_0 + tim1_1_1;
	       c_im(inout[12 * stride]) = tim1_1_0 - tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       {
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_1_2 - tim0_1_2);
		    tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_1_2 + tre0_1_2);
		    tre1_0_0 = tre0_1_0 + tre2_1_0;
		    tim1_0_0 = tim0_1_0 + tim2_1_0;
		    tre1_1_0 = tre0_1_0 - tre2_1_0;
		    tim1_1_0 = tim0_1_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_1_1) - (((FFTW_REAL) FFTW_K382683432) * tim0_1_1);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_1_1) + (((FFTW_REAL) FFTW_K382683432) * tre0_1_1);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_1_3) - (((FFTW_REAL) FFTW_K923879532) * tim0_1_3);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_1_3) + (((FFTW_REAL) FFTW_K923879532) * tre0_1_3);
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       c_re(inout[stride]) = tre1_0_0 + tre1_0_1;
	       c_im(inout[stride]) = tim1_0_0 + tim1_0_1;
	       c_re(inout[9 * stride]) = tre1_0_0 - tre1_0_1;
	       c_im(inout[9 * stride]) = tim1_0_0 - tim1_0_1;
	       c_re(inout[5 * stride]) = tre1_1_0 - tim1_1_1;
	       c_im(inout[5 * stride]) = tim1_1_0 + tre1_1_1;
	       c_re(inout[13 * stride]) = tre1_1_0 + tim1_1_1;
	       c_im(inout[13 * stride]) = tim1_1_0 - tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       tre1_0_0 = tre0_2_0 - tim0_2_2;
	       tim1_0_0 = tim0_2_0 + tre0_2_2;
	       tre1_1_0 = tre0_2_0 + tim0_2_2;
	       tim1_1_0 = tim0_2_0 - tre0_2_2;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_2_1 - tim0_2_1);
		    tim2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_2_1 + tre0_2_1);
		    tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_2_3 + tim0_2_3);
		    tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_2_3 - tim0_2_3);
		    tre1_0_1 = tre2_0_0 - tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 + tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       c_re(inout[2 * stride]) = tre1_0_0 + tre1_0_1;
	       c_im(inout[2 * stride]) = tim1_0_0 + tim1_0_1;
	       c_re(inout[10 * stride]) = tre1_0_0 - tre1_0_1;
	       c_im(inout[10 * stride]) = tim1_0_0 - tim1_0_1;
	       c_re(inout[6 * stride]) = tre1_1_0 - tim1_1_1;
	       c_im(inout[6 * stride]) = tim1_1_0 + tre1_1_1;
	       c_re(inout[14 * stride]) = tre1_1_0 + tim1_1_1;
	       c_im(inout[14 * stride]) = tim1_1_0 - tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       {
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_3_2 + tim0_3_2);
		    tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_3_2 - tim0_3_2);
		    tre1_0_0 = tre0_3_0 - tre2_1_0;
		    tim1_0_0 = tim0_3_0 + tim2_1_0;
		    tre1_1_0 = tre0_3_0 + tre2_1_0;
		    tim1_1_0 = tim0_3_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_3_1) - (((FFTW_REAL) FFTW_K923879532) * tim0_3_1);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_3_1) + (((FFTW_REAL) FFTW_K923879532) * tre0_3_1);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_3_3) - (((FFTW_REAL) FFTW_K923879532) * tre0_3_3);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_3_3) + (((FFTW_REAL) FFTW_K382683432) * tre0_3_3);
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 - tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 + tim2_1_0;
	       }
	       c_re(inout[3 * stride]) = tre1_0_0 + tre1_0_1;
	       c_im(inout[3 * stride]) = tim1_0_0 + tim1_0_1;
	       c_re(inout[11 * stride]) = tre1_0_0 - tre1_0_1;
	       c_im(inout[11 * stride]) = tim1_0_0 - tim1_0_1;
	       c_re(inout[7 * stride]) = tre1_1_0 - tim1_1_1;
	       c_im(inout[7 * stride]) = tim1_1_0 + tre1_1_1;
	       c_re(inout[15 * stride]) = tre1_1_0 + tim1_1_1;
	       c_im(inout[15 * stride]) = tim1_1_0 - tre1_1_1;
	  }
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 6 FP additions and 4 FP multiplications */

void fftwi_twiddle_2(FFTW_COMPLEX *A, const FFTW_COMPLEX *W, int stride, int m, int dist)
{
     int i;
     COMPLEX *inout;
     inout = A;
     for (i = 0; i < m; i = i + 1, inout = inout + dist, W = W + 1) {
	  FFTW_REAL tre0_0_0;
	  FFTW_REAL tim0_0_0;
	  FFTW_REAL tre0_1_0;
	  FFTW_REAL tim0_1_0;
	  tre0_0_0 = c_re(inout[0]);
	  tim0_0_0 = c_im(inout[0]);
	  {
	       FFTW_REAL tr;
	       FFTW_REAL ti;
	       FFTW_REAL twr;
	       FFTW_REAL twi;
	       tr = c_re(inout[stride]);
	       ti = c_im(inout[stride]);
	       twr = c_re(W[0]);
	       twi = c_im(W[0]);
	       tre0_1_0 = (tr * twr) + (ti * twi);
	       tim0_1_0 = (ti * twr) - (tr * twi);
	  }
	  c_re(inout[0]) = tre0_0_0 + tre0_1_0;
	  c_im(inout[0]) = tim0_0_0 + tim0_1_0;
	  c_re(inout[stride]) = tre0_0_0 - tre0_1_0;
	  c_im(inout[stride]) = tim0_0_0 - tim0_1_0;
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 18 FP additions and 12 FP multiplications */

void fftwi_twiddle_3(FFTW_COMPLEX *A, const FFTW_COMPLEX *W, int stride, int m, int dist)
{
     int i;
     COMPLEX *inout;
     inout = A;
     for (i = 0; i < m; i = i + 1, inout = inout + dist, W = W + 2) {
	  FFTW_REAL tre0_0_0;
	  FFTW_REAL tim0_0_0;
	  FFTW_REAL tre0_1_0;
	  FFTW_REAL tim0_1_0;
	  FFTW_REAL tre0_2_0;
	  FFTW_REAL tim0_2_0;
	  tre0_0_0 = c_re(inout[0]);
	  tim0_0_0 = c_im(inout[0]);
	  {
	       FFTW_REAL tr;
	       FFTW_REAL ti;
	       FFTW_REAL twr;
	       FFTW_REAL twi;
	       tr = c_re(inout[stride]);
	       ti = c_im(inout[stride]);
	       twr = c_re(W[0]);
	       twi = c_im(W[0]);
	       tre0_1_0 = (tr * twr) + (ti * twi);
	       tim0_1_0 = (ti * twr) - (tr * twi);
	  }
	  {
	       FFTW_REAL tr;
	       FFTW_REAL ti;
	       FFTW_REAL twr;
	       FFTW_REAL twi;
	       tr = c_re(inout[2 * stride]);
	       ti = c_im(inout[2 * stride]);
	       twr = c_re(W[1]);
	       twi = c_im(W[1]);
	       tre0_2_0 = (tr * twr) + (ti * twi);
	       tim0_2_0 = (ti * twr) - (tr * twi);
	  }
	  c_re(inout[0]) = tre0_0_0 + tre0_1_0 + tre0_2_0;
	  c_im(inout[0]) = tim0_0_0 + tim0_1_0 + tim0_2_0;
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tre1_1_0;
	       tre1_0_0 = tre0_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre0_1_0 + tre0_2_0));
	       tre1_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim0_2_0 - tim0_1_0);
	       c_re(inout[stride]) = tre1_0_0 + tre1_1_0;
	       c_re(inout[2 * stride]) = tre1_0_0 - tre1_1_0;
	  }
	  {
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tim1_1_0;
	       tim1_0_0 = tim0_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim0_1_0 + tim0_2_0));
	       tim1_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre0_1_0 - tre0_2_0);
	       c_im(inout[stride]) = tim1_0_0 + tim1_1_0;
	       c_im(inout[2 * stride]) = tim1_0_0 - tim1_1_0;
	  }
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 438 FP additions and 212 FP multiplications */

void fftwi_twiddle_32(FFTW_COMPLEX *A, const FFTW_COMPLEX *W, int stride, int m, int dist)
{
     int i;
     COMPLEX *inout;
     inout = A;
     for (i = 0; i < m; i = i + 1, inout = inout + dist, W = W + 31) {
	  FFTW_REAL tre0_0_0;
	  FFTW_REAL tim0_0_0;
	  FFTW_REAL tre0_0_1;
	  FFTW_REAL tim0_0_1;
	  FFTW_REAL tre0_0_2;
	  FFTW_REAL tim0_0_2;
	  FFTW_REAL tre0_0_3;
	  FFTW_REAL tim0_0_3;
	  FFTW_REAL tre0_0_4;
	  FFTW_REAL tim0_0_4;
	  FFTW_REAL tre0_0_5;
	  FFTW_REAL tim0_0_5;
	  FFTW_REAL tre0_0_6;
	  FFTW_REAL tim0_0_6;
	  FFTW_REAL tre0_0_7;
	  FFTW_REAL tim0_0_7;
	  FFTW_REAL tre0_1_0;
	  FFTW_REAL tim0_1_0;
	  FFTW_REAL tre0_1_1;
	  FFTW_REAL tim0_1_1;
	  FFTW_REAL tre0_1_2;
	  FFTW_REAL tim0_1_2;
	  FFTW_REAL tre0_1_3;
	  FFTW_REAL tim0_1_3;
	  FFTW_REAL tre0_1_4;
	  FFTW_REAL tim0_1_4;
	  FFTW_REAL tre0_1_5;
	  FFTW_REAL tim0_1_5;
	  FFTW_REAL tre0_1_6;
	  FFTW_REAL tim0_1_6;
	  FFTW_REAL tre0_1_7;
	  FFTW_REAL tim0_1_7;
	  FFTW_REAL tre0_2_0;
	  FFTW_REAL tim0_2_0;
	  FFTW_REAL tre0_2_1;
	  FFTW_REAL tim0_2_1;
	  FFTW_REAL tre0_2_2;
	  FFTW_REAL tim0_2_2;
	  FFTW_REAL tre0_2_3;
	  FFTW_REAL tim0_2_3;
	  FFTW_REAL tre0_2_4;
	  FFTW_REAL tim0_2_4;
	  FFTW_REAL tre0_2_5;
	  FFTW_REAL tim0_2_5;
	  FFTW_REAL tre0_2_6;
	  FFTW_REAL tim0_2_6;
	  FFTW_REAL tre0_2_7;
	  FFTW_REAL tim0_2_7;
	  FFTW_REAL tre0_3_0;
	  FFTW_REAL tim0_3_0;
	  FFTW_REAL tre0_3_1;
	  FFTW_REAL tim0_3_1;
	  FFTW_REAL tre0_3_2;
	  FFTW_REAL tim0_3_2;
	  FFTW_REAL tre0_3_3;
	  FFTW_REAL tim0_3_3;
	  FFTW_REAL tre0_3_4;
	  FFTW_REAL tim0_3_4;
	  FFTW_REAL tre0_3_5;
	  FFTW_REAL tim0_3_5;
	  FFTW_REAL tre0_3_6;
	  FFTW_REAL tim0_3_6;
	  FFTW_REAL tre0_3_7;
	  FFTW_REAL tim0_3_7;
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = c_re(inout[0]);
		    tim2_0_0 = c_im(inout[0]);
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[16 * stride]);
			 ti = c_im(inout[16 * stride]);
			 twr = c_re(W[15]);
			 twi = c_im(W[15]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[8 * stride]);
			 ti = c_im(inout[8 * stride]);
			 twr = c_re(W[7]);
			 twi = c_im(W[7]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[24 * stride]);
			 ti = c_im(inout[24 * stride]);
			 twr = c_re(W[23]);
			 twi = c_im(W[23]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       tre0_0_0 = tre1_0_0 + tre1_0_1;
	       tim0_0_0 = tim1_0_0 + tim1_0_1;
	       tre0_2_0 = tre1_0_0 - tre1_0_1;
	       tim0_2_0 = tim1_0_0 - tim1_0_1;
	       tre0_1_0 = tre1_1_0 - tim1_1_1;
	       tim0_1_0 = tim1_1_0 + tre1_1_1;
	       tre0_3_0 = tre1_1_0 + tim1_1_1;
	       tim0_3_0 = tim1_1_0 - tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[stride]);
			 ti = c_im(inout[stride]);
			 twr = c_re(W[0]);
			 twi = c_im(W[0]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[17 * stride]);
			 ti = c_im(inout[17 * stride]);
			 twr = c_re(W[16]);
			 twi = c_im(W[16]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[9 * stride]);
			 ti = c_im(inout[9 * stride]);
			 twr = c_re(W[8]);
			 twi = c_im(W[8]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[25 * stride]);
			 ti = c_im(inout[25 * stride]);
			 twr = c_re(W[24]);
			 twi = c_im(W[24]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       tre0_0_1 = tre1_0_0 + tre1_0_1;
	       tim0_0_1 = tim1_0_0 + tim1_0_1;
	       tre0_2_1 = tre1_0_0 - tre1_0_1;
	       tim0_2_1 = tim1_0_0 - tim1_0_1;
	       tre0_1_1 = tre1_1_0 - tim1_1_1;
	       tim0_1_1 = tim1_1_0 + tre1_1_1;
	       tre0_3_1 = tre1_1_0 + tim1_1_1;
	       tim0_3_1 = tim1_1_0 - tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[2 * stride]);
			 ti = c_im(inout[2 * stride]);
			 twr = c_re(W[1]);
			 twi = c_im(W[1]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[18 * stride]);
			 ti = c_im(inout[18 * stride]);
			 twr = c_re(W[17]);
			 twi = c_im(W[17]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[10 * stride]);
			 ti = c_im(inout[10 * stride]);
			 twr = c_re(W[9]);
			 twi = c_im(W[9]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[26 * stride]);
			 ti = c_im(inout[26 * stride]);
			 twr = c_re(W[25]);
			 twi = c_im(W[25]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       tre0_0_2 = tre1_0_0 + tre1_0_1;
	       tim0_0_2 = tim1_0_0 + tim1_0_1;
	       tre0_2_2 = tre1_0_0 - tre1_0_1;
	       tim0_2_2 = tim1_0_0 - tim1_0_1;
	       tre0_1_2 = tre1_1_0 - tim1_1_1;
	       tim0_1_2 = tim1_1_0 + tre1_1_1;
	       tre0_3_2 = tre1_1_0 + tim1_1_1;
	       tim0_3_2 = tim1_1_0 - tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[3 * stride]);
			 ti = c_im(inout[3 * stride]);
			 twr = c_re(W[2]);
			 twi = c_im(W[2]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[19 * stride]);
			 ti = c_im(inout[19 * stride]);
			 twr = c_re(W[18]);
			 twi = c_im(W[18]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[11 * stride]);
			 ti = c_im(inout[11 * stride]);
			 twr = c_re(W[10]);
			 twi = c_im(W[10]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[27 * stride]);
			 ti = c_im(inout[27 * stride]);
			 twr = c_re(W[26]);
			 twi = c_im(W[26]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       tre0_0_3 = tre1_0_0 + tre1_0_1;
	       tim0_0_3 = tim1_0_0 + tim1_0_1;
	       tre0_2_3 = tre1_0_0 - tre1_0_1;
	       tim0_2_3 = tim1_0_0 - tim1_0_1;
	       tre0_1_3 = tre1_1_0 - tim1_1_1;
	       tim0_1_3 = tim1_1_0 + tre1_1_1;
	       tre0_3_3 = tre1_1_0 + tim1_1_1;
	       tim0_3_3 = tim1_1_0 - tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[4 * stride]);
			 ti = c_im(inout[4 * stride]);
			 twr = c_re(W[3]);
			 twi = c_im(W[3]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[20 * stride]);
			 ti = c_im(inout[20 * stride]);
			 twr = c_re(W[19]);
			 twi = c_im(W[19]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[12 * stride]);
			 ti = c_im(inout[12 * stride]);
			 twr = c_re(W[11]);
			 twi = c_im(W[11]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[28 * stride]);
			 ti = c_im(inout[28 * stride]);
			 twr = c_re(W[27]);
			 twi = c_im(W[27]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       tre0_0_4 = tre1_0_0 + tre1_0_1;
	       tim0_0_4 = tim1_0_0 + tim1_0_1;
	       tre0_2_4 = tre1_0_0 - tre1_0_1;
	       tim0_2_4 = tim1_0_0 - tim1_0_1;
	       tre0_1_4 = tre1_1_0 - tim1_1_1;
	       tim0_1_4 = tim1_1_0 + tre1_1_1;
	       tre0_3_4 = tre1_1_0 + tim1_1_1;
	       tim0_3_4 = tim1_1_0 - tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[5 * stride]);
			 ti = c_im(inout[5 * stride]);
			 twr = c_re(W[4]);
			 twi = c_im(W[4]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[21 * stride]);
			 ti = c_im(inout[21 * stride]);
			 twr = c_re(W[20]);
			 twi = c_im(W[20]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[13 * stride]);
			 ti = c_im(inout[13 * stride]);
			 twr = c_re(W[12]);
			 twi = c_im(W[12]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[29 * stride]);
			 ti = c_im(inout[29 * stride]);
			 twr = c_re(W[28]);
			 twi = c_im(W[28]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       tre0_0_5 = tre1_0_0 + tre1_0_1;
	       tim0_0_5 = tim1_0_0 + tim1_0_1;
	       tre0_2_5 = tre1_0_0 - tre1_0_1;
	       tim0_2_5 = tim1_0_0 - tim1_0_1;
	       tre0_1_5 = tre1_1_0 - tim1_1_1;
	       tim0_1_5 = tim1_1_0 + tre1_1_1;
	       tre0_3_5 = tre1_1_0 + tim1_1_1;
	       tim0_3_5 = tim1_1_0 - tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[6 * stride]);
			 ti = c_im(inout[6 * stride]);
			 twr = c_re(W[5]);
			 twi = c_im(W[5]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[22 * stride]);
			 ti = c_im(inout[22 * stride]);
			 twr = c_re(W[21]);
			 twi = c_im(W[21]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[14 * stride]);
			 ti = c_im(inout[14 * stride]);
			 twr = c_re(W[13]);
			 twi = c_im(W[13]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[30 * stride]);
			 ti = c_im(inout[30 * stride]);
			 twr = c_re(W[29]);
			 twi = c_im(W[29]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       tre0_0_6 = tre1_0_0 + tre1_0_1;
	       tim0_0_6 = tim1_0_0 + tim1_0_1;
	       tre0_2_6 = tre1_0_0 - tre1_0_1;
	       tim0_2_6 = tim1_0_0 - tim1_0_1;
	       tre0_1_6 = tre1_1_0 - tim1_1_1;
	       tim0_1_6 = tim1_1_0 + tre1_1_1;
	       tre0_3_6 = tre1_1_0 + tim1_1_1;
	       tim0_3_6 = tim1_1_0 - tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[7 * stride]);
			 ti = c_im(inout[7 * stride]);
			 twr = c_re(W[6]);
			 twi = c_im(W[6]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[23 * stride]);
			 ti = c_im(inout[23 * stride]);
			 twr = c_re(W[22]);
			 twi = c_im(W[22]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[15 * stride]);
			 ti = c_im(inout[15 * stride]);
			 twr = c_re(W[14]);
			 twi = c_im(W[14]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[31 * stride]);
			 ti = c_im(inout[31 * stride]);
			 twr = c_re(W[30]);
			 twi = c_im(W[30]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       tre0_0_7 = tre1_0_0 + tre1_0_1;
	       tim0_0_7 = tim1_0_0 + tim1_0_1;
	       tre0_2_7 = tre1_0_0 - tre1_0_1;
	       tim0_2_7 = tim1_0_0 - tim1_0_1;
	       tre0_1_7 = tre1_1_0 - tim1_1_1;
	       tim0_1_7 = tim1_1_0 + tre1_1_1;
	       tre0_3_7 = tre1_1_0 + tim1_1_1;
	       tim0_3_7 = tim1_1_0 - tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       tre1_0_0 = tre0_0_0 + tre0_0_4;
	       tim1_0_0 = tim0_0_0 + tim0_0_4;
	       tre1_1_0 = tre0_0_0 - tre0_0_4;
	       tim1_1_0 = tim0_0_0 - tim0_0_4;
	       tre1_0_1 = tre0_0_1 + tre0_0_5;
	       tim1_0_1 = tim0_0_1 + tim0_0_5;
	       tre1_1_1 = tre0_0_1 - tre0_0_5;
	       tim1_1_1 = tim0_0_1 - tim0_0_5;
	       tre1_0_2 = tre0_0_2 + tre0_0_6;
	       tim1_0_2 = tim0_0_2 + tim0_0_6;
	       tre1_1_2 = tre0_0_2 - tre0_0_6;
	       tim1_1_2 = tim0_0_2 - tim0_0_6;
	       tre1_0_3 = tre0_0_3 + tre0_0_7;
	       tim1_0_3 = tim0_0_3 + tim0_0_7;
	       tre1_1_3 = tre0_0_3 - tre0_0_7;
	       tim1_1_3 = tim0_0_3 - tim0_0_7;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    c_re(inout[0]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[0]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[16 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[16 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[8 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[8 * stride]) = tim2_1_0 + tre2_1_1;
		    c_re(inout[24 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[24 * stride]) = tim2_1_0 - tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 - tim1_1_2;
		    tim2_0_0 = tim1_1_0 + tre1_1_2;
		    tre2_1_0 = tre1_1_0 + tim1_1_2;
		    tim2_1_0 = tim1_1_0 - tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
			 tre2_0_1 = tre3_0_0 - tre3_1_0;
			 tim2_0_1 = tim3_0_0 + tim3_1_0;
			 tre2_1_1 = tre3_0_0 + tre3_1_0;
			 tim2_1_1 = tim3_0_0 - tim3_1_0;
		    }
		    c_re(inout[4 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[4 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[20 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[20 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[12 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[12 * stride]) = tim2_1_0 + tre2_1_1;
		    c_re(inout[28 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[28 * stride]) = tim2_1_0 - tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_1_4 - tim0_1_4);
		    tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_1_4 + tre0_1_4);
		    tre1_0_0 = tre0_1_0 + tre2_1_0;
		    tim1_0_0 = tim0_1_0 + tim2_1_0;
		    tre1_1_0 = tre0_1_0 - tre2_1_0;
		    tim1_1_0 = tim0_1_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tre0_1_1) - (((FFTW_REAL) FFTW_K195090322) * tim0_1_1);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tim0_1_1) + (((FFTW_REAL) FFTW_K195090322) * tre0_1_1);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K555570233) * tre0_1_5) - (((FFTW_REAL) FFTW_K831469612) * tim0_1_5);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K555570233) * tim0_1_5) + (((FFTW_REAL) FFTW_K831469612) * tre0_1_5);
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_1_2) - (((FFTW_REAL) FFTW_K382683432) * tim0_1_2);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_1_2) + (((FFTW_REAL) FFTW_K382683432) * tre0_1_2);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_1_6) - (((FFTW_REAL) FFTW_K923879532) * tim0_1_6);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_1_6) + (((FFTW_REAL) FFTW_K923879532) * tre0_1_6);
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_1_3) - (((FFTW_REAL) FFTW_K555570233) * tim0_1_3);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_1_3) + (((FFTW_REAL) FFTW_K555570233) * tre0_1_3);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tre0_1_7) - (((FFTW_REAL) FFTW_K980785280) * tim0_1_7);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tim0_1_7) + (((FFTW_REAL) FFTW_K980785280) * tre0_1_7);
		    tre1_0_3 = tre2_0_0 + tre2_1_0;
		    tim1_0_3 = tim2_0_0 + tim2_1_0;
		    tre1_1_3 = tre2_0_0 - tre2_1_0;
		    tim1_1_3 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    c_re(inout[stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[17 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[17 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[9 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[9 * stride]) = tim2_1_0 + tre2_1_1;
		    c_re(inout[25 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[25 * stride]) = tim2_1_0 - tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 - tim1_1_2;
		    tim2_0_0 = tim1_1_0 + tre1_1_2;
		    tre2_1_0 = tre1_1_0 + tim1_1_2;
		    tim2_1_0 = tim1_1_0 - tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
			 tre2_0_1 = tre3_0_0 - tre3_1_0;
			 tim2_0_1 = tim3_0_0 + tim3_1_0;
			 tre2_1_1 = tre3_0_0 + tre3_1_0;
			 tim2_1_1 = tim3_0_0 - tim3_1_0;
		    }
		    c_re(inout[5 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[5 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[21 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[21 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[13 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[13 * stride]) = tim2_1_0 + tre2_1_1;
		    c_re(inout[29 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[29 * stride]) = tim2_1_0 - tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       tre1_0_0 = tre0_2_0 - tim0_2_4;
	       tim1_0_0 = tim0_2_0 + tre0_2_4;
	       tre1_1_0 = tre0_2_0 + tim0_2_4;
	       tim1_1_0 = tim0_2_0 - tre0_2_4;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_2_1) - (((FFTW_REAL) FFTW_K382683432) * tim0_2_1);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_2_1) + (((FFTW_REAL) FFTW_K382683432) * tre0_2_1);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_2_5) + (((FFTW_REAL) FFTW_K923879532) * tim0_2_5);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_2_5) - (((FFTW_REAL) FFTW_K382683432) * tim0_2_5);
		    tre1_0_1 = tre2_0_0 - tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 + tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_2_2 - tim0_2_2);
		    tim2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_2_2 + tre0_2_2);
		    tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_2_6 + tim0_2_6);
		    tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_2_6 - tim0_2_6);
		    tre1_0_2 = tre2_0_0 - tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 + tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_2_3) - (((FFTW_REAL) FFTW_K923879532) * tim0_2_3);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_2_3) + (((FFTW_REAL) FFTW_K923879532) * tre0_2_3);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_2_7) + (((FFTW_REAL) FFTW_K382683432) * tim0_2_7);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_2_7) - (((FFTW_REAL) FFTW_K923879532) * tim0_2_7);
		    tre1_0_3 = tre2_0_0 - tre2_1_0;
		    tim1_0_3 = tim2_0_0 + tim2_1_0;
		    tre1_1_3 = tre2_0_0 + tre2_1_0;
		    tim1_1_3 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    c_re(inout[2 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[2 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[18 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[18 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[10 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[10 * stride]) = tim2_1_0 + tre2_1_1;
		    c_re(inout[26 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[26 * stride]) = tim2_1_0 - tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 - tim1_1_2;
		    tim2_0_0 = tim1_1_0 + tre1_1_2;
		    tre2_1_0 = tre1_1_0 + tim1_1_2;
		    tim2_1_0 = tim1_1_0 - tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
			 tre2_0_1 = tre3_0_0 - tre3_1_0;
			 tim2_0_1 = tim3_0_0 + tim3_1_0;
			 tre2_1_1 = tre3_0_0 + tre3_1_0;
			 tim2_1_1 = tim3_0_0 - tim3_1_0;
		    }
		    c_re(inout[6 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[6 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[22 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[22 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[14 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[14 * stride]) = tim2_1_0 + tre2_1_1;
		    c_re(inout[30 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[30 * stride]) = tim2_1_0 - tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_3_4 + tim0_3_4);
		    tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_3_4 - tim0_3_4);
		    tre1_0_0 = tre0_3_0 - tre2_1_0;
		    tim1_0_0 = tim0_3_0 + tim2_1_0;
		    tre1_1_0 = tre0_3_0 + tre2_1_0;
		    tim1_1_0 = tim0_3_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_3_1) - (((FFTW_REAL) FFTW_K555570233) * tim0_3_1);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_3_1) + (((FFTW_REAL) FFTW_K555570233) * tre0_3_1);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K980785280) * tre0_3_5) + (((FFTW_REAL) FFTW_K195090322) * tim0_3_5);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tre0_3_5) - (((FFTW_REAL) FFTW_K980785280) * tim0_3_5);
		    tre1_0_1 = tre2_0_0 - tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 + tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_3_2) - (((FFTW_REAL) FFTW_K923879532) * tim0_3_2);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_3_2) + (((FFTW_REAL) FFTW_K923879532) * tre0_3_2);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_3_6) - (((FFTW_REAL) FFTW_K923879532) * tre0_3_6);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_3_6) + (((FFTW_REAL) FFTW_K382683432) * tre0_3_6);
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 - tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 + tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K195090322) * tre0_3_3) + (((FFTW_REAL) FFTW_K980785280) * tim0_3_3);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tre0_3_3) - (((FFTW_REAL) FFTW_K195090322) * tim0_3_3);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_3_7) - (((FFTW_REAL) FFTW_K555570233) * tre0_3_7);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K555570233) * tim0_3_7) + (((FFTW_REAL) FFTW_K831469612) * tre0_3_7);
		    tre1_0_3 = tre2_1_0 - tre2_0_0;
		    tim1_0_3 = tim2_0_0 - tim2_1_0;
		    tre1_1_3 = (-(tre2_0_0 + tre2_1_0));
		    tim1_1_3 = tim2_0_0 + tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    c_re(inout[3 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[3 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[19 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[19 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[11 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[11 * stride]) = tim2_1_0 + tre2_1_1;
		    c_re(inout[27 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[27 * stride]) = tim2_1_0 - tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 - tim1_1_2;
		    tim2_0_0 = tim1_1_0 + tre1_1_2;
		    tre2_1_0 = tre1_1_0 + tim1_1_2;
		    tim2_1_0 = tim1_1_0 - tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
			 tre2_0_1 = tre3_0_0 - tre3_1_0;
			 tim2_0_1 = tim3_0_0 + tim3_1_0;
			 tre2_1_1 = tre3_0_0 + tre3_1_0;
			 tim2_1_1 = tim3_0_0 - tim3_1_0;
		    }
		    c_re(inout[7 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[7 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[23 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[23 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[15 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[15 * stride]) = tim2_1_0 + tre2_1_1;
		    c_re(inout[31 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[31 * stride]) = tim2_1_0 - tre2_1_1;
	       }
	  }
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 22 FP additions and 12 FP multiplications */

void fftwi_twiddle_4(FFTW_COMPLEX *A, const FFTW_COMPLEX *W, int stride, int m, int dist)
{
     int i;
     COMPLEX *inout;
     inout = A;
     for (i = 0; i < m; i = i + 1, inout = inout + dist, W = W + 3) {
	  FFTW_REAL tre0_0_0;
	  FFTW_REAL tim0_0_0;
	  FFTW_REAL tre0_0_1;
	  FFTW_REAL tim0_0_1;
	  FFTW_REAL tre0_1_0;
	  FFTW_REAL tim0_1_0;
	  FFTW_REAL tre0_1_1;
	  FFTW_REAL tim0_1_1;
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       tre1_0_0 = c_re(inout[0]);
	       tim1_0_0 = c_im(inout[0]);
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[2 * stride]);
		    ti = c_im(inout[2 * stride]);
		    twr = c_re(W[1]);
		    twi = c_im(W[1]);
		    tre1_1_0 = (tr * twr) + (ti * twi);
		    tim1_1_0 = (ti * twr) - (tr * twi);
	       }
	       tre0_0_0 = tre1_0_0 + tre1_1_0;
	       tim0_0_0 = tim1_0_0 + tim1_1_0;
	       tre0_1_0 = tre1_0_0 - tre1_1_0;
	       tim0_1_0 = tim1_0_0 - tim1_1_0;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[stride]);
		    ti = c_im(inout[stride]);
		    twr = c_re(W[0]);
		    twi = c_im(W[0]);
		    tre1_0_0 = (tr * twr) + (ti * twi);
		    tim1_0_0 = (ti * twr) - (tr * twi);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[3 * stride]);
		    ti = c_im(inout[3 * stride]);
		    twr = c_re(W[2]);
		    twi = c_im(W[2]);
		    tre1_1_0 = (tr * twr) + (ti * twi);
		    tim1_1_0 = (ti * twr) - (tr * twi);
	       }
	       tre0_0_1 = tre1_0_0 + tre1_1_0;
	       tim0_0_1 = tim1_0_0 + tim1_1_0;
	       tre0_1_1 = tre1_0_0 - tre1_1_0;
	       tim0_1_1 = tim1_0_0 - tim1_1_0;
	  }
	  c_re(inout[0]) = tre0_0_0 + tre0_0_1;
	  c_im(inout[0]) = tim0_0_0 + tim0_0_1;
	  c_re(inout[2 * stride]) = tre0_0_0 - tre0_0_1;
	  c_im(inout[2 * stride]) = tim0_0_0 - tim0_0_1;
	  c_re(inout[stride]) = tre0_1_0 - tim0_1_1;
	  c_im(inout[stride]) = tim0_1_0 + tre0_1_1;
	  c_re(inout[3 * stride]) = tre0_1_0 + tim0_1_1;
	  c_im(inout[3 * stride]) = tim0_1_0 - tre0_1_1;
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 52 FP additions and 32 FP multiplications */

void fftwi_twiddle_5(FFTW_COMPLEX *A, const FFTW_COMPLEX *W, int stride, int m, int dist)
{
     int i;
     COMPLEX *inout;
     inout = A;
     for (i = 0; i < m; i = i + 1, inout = inout + dist, W = W + 4) {
	  FFTW_REAL tre0_0_0;
	  FFTW_REAL tim0_0_0;
	  FFTW_REAL tre0_1_0;
	  FFTW_REAL tim0_1_0;
	  FFTW_REAL tre0_2_0;
	  FFTW_REAL tim0_2_0;
	  FFTW_REAL tre0_3_0;
	  FFTW_REAL tim0_3_0;
	  FFTW_REAL tre0_4_0;
	  FFTW_REAL tim0_4_0;
	  tre0_0_0 = c_re(inout[0]);
	  tim0_0_0 = c_im(inout[0]);
	  {
	       FFTW_REAL tr;
	       FFTW_REAL ti;
	       FFTW_REAL twr;
	       FFTW_REAL twi;
	       tr = c_re(inout[stride]);
	       ti = c_im(inout[stride]);
	       twr = c_re(W[0]);
	       twi = c_im(W[0]);
	       tre0_1_0 = (tr * twr) + (ti * twi);
	       tim0_1_0 = (ti * twr) - (tr * twi);
	  }
	  {
	       FFTW_REAL tr;
	       FFTW_REAL ti;
	       FFTW_REAL twr;
	       FFTW_REAL twi;
	       tr = c_re(inout[2 * stride]);
	       ti = c_im(inout[2 * stride]);
	       twr = c_re(W[1]);
	       twi = c_im(W[1]);
	       tre0_2_0 = (tr * twr) + (ti * twi);
	       tim0_2_0 = (ti * twr) - (tr * twi);
	  }
	  {
	       FFTW_REAL tr;
	       FFTW_REAL ti;
	       FFTW_REAL twr;
	       FFTW_REAL twi;
	       tr = c_re(inout[3 * stride]);
	       ti = c_im(inout[3 * stride]);
	       twr = c_re(W[2]);
	       twi = c_im(W[2]);
	       tre0_3_0 = (tr * twr) + (ti * twi);
	       tim0_3_0 = (ti * twr) - (tr * twi);
	  }
	  {
	       FFTW_REAL tr;
	       FFTW_REAL ti;
	       FFTW_REAL twr;
	       FFTW_REAL twi;
	       tr = c_re(inout[4 * stride]);
	       ti = c_im(inout[4 * stride]);
	       twr = c_re(W[3]);
	       twi = c_im(W[3]);
	       tre0_4_0 = (tr * twr) + (ti * twi);
	       tim0_4_0 = (ti * twr) - (tr * twi);
	  }
	  c_re(inout[0]) = tre0_0_0 + tre0_1_0 + tre0_2_0 + tre0_3_0 + tre0_4_0;
	  c_im(inout[0]) = tim0_0_0 + tim0_1_0 + tim0_2_0 + tim0_3_0 + tim0_4_0;
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tre1_1_0;
	       tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_1_0 + tre0_4_0)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_2_0 + tre0_3_0));
	       tre1_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tim0_4_0 - tim0_1_0)) + (((FFTW_REAL) FFTW_K587785252) * (tim0_3_0 - tim0_2_0));
	       c_re(inout[stride]) = tre1_0_0 + tre1_1_0;
	       c_re(inout[4 * stride]) = tre1_0_0 - tre1_1_0;
	  }
	  {
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tim1_1_0;
	       tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_1_0 + tim0_4_0)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_2_0 + tim0_3_0));
	       tim1_1_0 = (((FFTW_REAL) FFTW_K951056516) * (tre0_1_0 - tre0_4_0)) + (((FFTW_REAL) FFTW_K587785252) * (tre0_2_0 - tre0_3_0));
	       c_im(inout[stride]) = tim1_0_0 + tim1_1_0;
	       c_im(inout[4 * stride]) = tim1_0_0 - tim1_1_0;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tre1_1_0;
	       tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tre0_2_0 + tre0_3_0)) - (((FFTW_REAL) FFTW_K809016994) * (tre0_1_0 + tre0_4_0));
	       tre1_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tim0_4_0 - tim0_1_0)) + (((FFTW_REAL) FFTW_K951056516) * (tim0_2_0 - tim0_3_0));
	       c_re(inout[2 * stride]) = tre1_0_0 + tre1_1_0;
	       c_re(inout[3 * stride]) = tre1_0_0 - tre1_1_0;
	  }
	  {
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tim1_1_0;
	       tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K309016994) * (tim0_2_0 + tim0_3_0)) - (((FFTW_REAL) FFTW_K809016994) * (tim0_1_0 + tim0_4_0));
	       tim1_1_0 = (((FFTW_REAL) FFTW_K587785252) * (tre0_1_0 - tre0_4_0)) + (((FFTW_REAL) FFTW_K951056516) * (tre0_3_0 - tre0_2_0));
	       c_im(inout[2 * stride]) = tim1_0_0 + tim1_1_0;
	       c_im(inout[3 * stride]) = tim1_0_0 - tim1_1_0;
	  }
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 50 FP additions and 28 FP multiplications */

void fftwi_twiddle_6(FFTW_COMPLEX *A, const FFTW_COMPLEX *W, int stride, int m, int dist)
{
     int i;
     COMPLEX *inout;
     inout = A;
     for (i = 0; i < m; i = i + 1, inout = inout + dist, W = W + 5) {
	  FFTW_REAL tre0_0_0;
	  FFTW_REAL tim0_0_0;
	  FFTW_REAL tre0_0_1;
	  FFTW_REAL tim0_0_1;
	  FFTW_REAL tre0_0_2;
	  FFTW_REAL tim0_0_2;
	  FFTW_REAL tre0_1_0;
	  FFTW_REAL tim0_1_0;
	  FFTW_REAL tre0_1_1;
	  FFTW_REAL tim0_1_1;
	  FFTW_REAL tre0_1_2;
	  FFTW_REAL tim0_1_2;
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       tre1_0_0 = c_re(inout[0]);
	       tim1_0_0 = c_im(inout[0]);
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[3 * stride]);
		    ti = c_im(inout[3 * stride]);
		    twr = c_re(W[2]);
		    twi = c_im(W[2]);
		    tre1_1_0 = (tr * twr) + (ti * twi);
		    tim1_1_0 = (ti * twr) - (tr * twi);
	       }
	       tre0_0_0 = tre1_0_0 + tre1_1_0;
	       tim0_0_0 = tim1_0_0 + tim1_1_0;
	       tre0_1_0 = tre1_0_0 - tre1_1_0;
	       tim0_1_0 = tim1_0_0 - tim1_1_0;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[2 * stride]);
		    ti = c_im(inout[2 * stride]);
		    twr = c_re(W[1]);
		    twi = c_im(W[1]);
		    tre1_0_0 = (tr * twr) + (ti * twi);
		    tim1_0_0 = (ti * twr) - (tr * twi);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[5 * stride]);
		    ti = c_im(inout[5 * stride]);
		    twr = c_re(W[4]);
		    twi = c_im(W[4]);
		    tre1_1_0 = (tr * twr) + (ti * twi);
		    tim1_1_0 = (ti * twr) - (tr * twi);
	       }
	       tre0_0_1 = tre1_0_0 + tre1_1_0;
	       tim0_0_1 = tim1_0_0 + tim1_1_0;
	       tre0_1_1 = tre1_0_0 - tre1_1_0;
	       tim0_1_1 = tim1_0_0 - tim1_1_0;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[4 * stride]);
		    ti = c_im(inout[4 * stride]);
		    twr = c_re(W[3]);
		    twi = c_im(W[3]);
		    tre1_0_0 = (tr * twr) + (ti * twi);
		    tim1_0_0 = (ti * twr) - (tr * twi);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[stride]);
		    ti = c_im(inout[stride]);
		    twr = c_re(W[0]);
		    twi = c_im(W[0]);
		    tre1_1_0 = (tr * twr) + (ti * twi);
		    tim1_1_0 = (ti * twr) - (tr * twi);
	       }
	       tre0_0_2 = tre1_0_0 + tre1_1_0;
	       tim0_0_2 = tim1_0_0 + tim1_1_0;
	       tre0_1_2 = tre1_0_0 - tre1_1_0;
	       tim0_1_2 = tim1_0_0 - tim1_1_0;
	  }
	  c_re(inout[0]) = tre0_0_0 + tre0_0_1 + tre0_0_2;
	  c_im(inout[0]) = tim0_0_0 + tim0_0_1 + tim0_0_2;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre0_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre0_0_1 + tre0_0_2));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim0_0_2 - tim0_0_1);
	       c_re(inout[4 * stride]) = tre2_0_0 + tre2_1_0;
	       c_re(inout[2 * stride]) = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim0_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim0_0_1 + tim0_0_2));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre0_0_1 - tre0_0_2);
	       c_im(inout[4 * stride]) = tim2_0_0 + tim2_1_0;
	       c_im(inout[2 * stride]) = tim2_0_0 - tim2_1_0;
	  }
	  c_re(inout[3 * stride]) = tre0_1_0 + tre0_1_1 + tre0_1_2;
	  c_im(inout[3 * stride]) = tim0_1_0 + tim0_1_1 + tim0_1_2;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre0_1_0 - (((FFTW_REAL) FFTW_K499999999) * (tre0_1_1 + tre0_1_2));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim0_1_2 - tim0_1_1);
	       c_re(inout[stride]) = tre2_0_0 + tre2_1_0;
	       c_re(inout[5 * stride]) = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim0_1_0 - (((FFTW_REAL) FFTW_K499999999) * (tim0_1_1 + tim0_1_2));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre0_1_1 - tre0_1_2);
	       c_im(inout[stride]) = tim2_0_0 + tim2_1_0;
	       c_im(inout[5 * stride]) = tim2_0_0 - tim2_1_0;
	  }
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 1054 FP additions and 500 FP multiplications */

void fftwi_twiddle_64(FFTW_COMPLEX *A, const FFTW_COMPLEX *W, int stride, int m, int dist)
{
     int i;
     COMPLEX *inout;
     inout = A;
     for (i = 0; i < m; i = i + 1, inout = inout + dist, W = W + 63) {
	  FFTW_REAL tre0_0_0;
	  FFTW_REAL tim0_0_0;
	  FFTW_REAL tre0_0_1;
	  FFTW_REAL tim0_0_1;
	  FFTW_REAL tre0_0_2;
	  FFTW_REAL tim0_0_2;
	  FFTW_REAL tre0_0_3;
	  FFTW_REAL tim0_0_3;
	  FFTW_REAL tre0_0_4;
	  FFTW_REAL tim0_0_4;
	  FFTW_REAL tre0_0_5;
	  FFTW_REAL tim0_0_5;
	  FFTW_REAL tre0_0_6;
	  FFTW_REAL tim0_0_6;
	  FFTW_REAL tre0_0_7;
	  FFTW_REAL tim0_0_7;
	  FFTW_REAL tre0_1_0;
	  FFTW_REAL tim0_1_0;
	  FFTW_REAL tre0_1_1;
	  FFTW_REAL tim0_1_1;
	  FFTW_REAL tre0_1_2;
	  FFTW_REAL tim0_1_2;
	  FFTW_REAL tre0_1_3;
	  FFTW_REAL tim0_1_3;
	  FFTW_REAL tre0_1_4;
	  FFTW_REAL tim0_1_4;
	  FFTW_REAL tre0_1_5;
	  FFTW_REAL tim0_1_5;
	  FFTW_REAL tre0_1_6;
	  FFTW_REAL tim0_1_6;
	  FFTW_REAL tre0_1_7;
	  FFTW_REAL tim0_1_7;
	  FFTW_REAL tre0_2_0;
	  FFTW_REAL tim0_2_0;
	  FFTW_REAL tre0_2_1;
	  FFTW_REAL tim0_2_1;
	  FFTW_REAL tre0_2_2;
	  FFTW_REAL tim0_2_2;
	  FFTW_REAL tre0_2_3;
	  FFTW_REAL tim0_2_3;
	  FFTW_REAL tre0_2_4;
	  FFTW_REAL tim0_2_4;
	  FFTW_REAL tre0_2_5;
	  FFTW_REAL tim0_2_5;
	  FFTW_REAL tre0_2_6;
	  FFTW_REAL tim0_2_6;
	  FFTW_REAL tre0_2_7;
	  FFTW_REAL tim0_2_7;
	  FFTW_REAL tre0_3_0;
	  FFTW_REAL tim0_3_0;
	  FFTW_REAL tre0_3_1;
	  FFTW_REAL tim0_3_1;
	  FFTW_REAL tre0_3_2;
	  FFTW_REAL tim0_3_2;
	  FFTW_REAL tre0_3_3;
	  FFTW_REAL tim0_3_3;
	  FFTW_REAL tre0_3_4;
	  FFTW_REAL tim0_3_4;
	  FFTW_REAL tre0_3_5;
	  FFTW_REAL tim0_3_5;
	  FFTW_REAL tre0_3_6;
	  FFTW_REAL tim0_3_6;
	  FFTW_REAL tre0_3_7;
	  FFTW_REAL tim0_3_7;
	  FFTW_REAL tre0_4_0;
	  FFTW_REAL tim0_4_0;
	  FFTW_REAL tre0_4_1;
	  FFTW_REAL tim0_4_1;
	  FFTW_REAL tre0_4_2;
	  FFTW_REAL tim0_4_2;
	  FFTW_REAL tre0_4_3;
	  FFTW_REAL tim0_4_3;
	  FFTW_REAL tre0_4_4;
	  FFTW_REAL tim0_4_4;
	  FFTW_REAL tre0_4_5;
	  FFTW_REAL tim0_4_5;
	  FFTW_REAL tre0_4_6;
	  FFTW_REAL tim0_4_6;
	  FFTW_REAL tre0_4_7;
	  FFTW_REAL tim0_4_7;
	  FFTW_REAL tre0_5_0;
	  FFTW_REAL tim0_5_0;
	  FFTW_REAL tre0_5_1;
	  FFTW_REAL tim0_5_1;
	  FFTW_REAL tre0_5_2;
	  FFTW_REAL tim0_5_2;
	  FFTW_REAL tre0_5_3;
	  FFTW_REAL tim0_5_3;
	  FFTW_REAL tre0_5_4;
	  FFTW_REAL tim0_5_4;
	  FFTW_REAL tre0_5_5;
	  FFTW_REAL tim0_5_5;
	  FFTW_REAL tre0_5_6;
	  FFTW_REAL tim0_5_6;
	  FFTW_REAL tre0_5_7;
	  FFTW_REAL tim0_5_7;
	  FFTW_REAL tre0_6_0;
	  FFTW_REAL tim0_6_0;
	  FFTW_REAL tre0_6_1;
	  FFTW_REAL tim0_6_1;
	  FFTW_REAL tre0_6_2;
	  FFTW_REAL tim0_6_2;
	  FFTW_REAL tre0_6_3;
	  FFTW_REAL tim0_6_3;
	  FFTW_REAL tre0_6_4;
	  FFTW_REAL tim0_6_4;
	  FFTW_REAL tre0_6_5;
	  FFTW_REAL tim0_6_5;
	  FFTW_REAL tre0_6_6;
	  FFTW_REAL tim0_6_6;
	  FFTW_REAL tre0_6_7;
	  FFTW_REAL tim0_6_7;
	  FFTW_REAL tre0_7_0;
	  FFTW_REAL tim0_7_0;
	  FFTW_REAL tre0_7_1;
	  FFTW_REAL tim0_7_1;
	  FFTW_REAL tre0_7_2;
	  FFTW_REAL tim0_7_2;
	  FFTW_REAL tre0_7_3;
	  FFTW_REAL tim0_7_3;
	  FFTW_REAL tre0_7_4;
	  FFTW_REAL tim0_7_4;
	  FFTW_REAL tre0_7_5;
	  FFTW_REAL tim0_7_5;
	  FFTW_REAL tre0_7_6;
	  FFTW_REAL tim0_7_6;
	  FFTW_REAL tre0_7_7;
	  FFTW_REAL tim0_7_7;
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = c_re(inout[0]);
		    tim2_0_0 = c_im(inout[0]);
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[32 * stride]);
			 ti = c_im(inout[32 * stride]);
			 twr = c_re(W[31]);
			 twi = c_im(W[31]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[8 * stride]);
			 ti = c_im(inout[8 * stride]);
			 twr = c_re(W[7]);
			 twi = c_im(W[7]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[40 * stride]);
			 ti = c_im(inout[40 * stride]);
			 twr = c_re(W[39]);
			 twi = c_im(W[39]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[16 * stride]);
			 ti = c_im(inout[16 * stride]);
			 twr = c_re(W[15]);
			 twi = c_im(W[15]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[48 * stride]);
			 ti = c_im(inout[48 * stride]);
			 twr = c_re(W[47]);
			 twi = c_im(W[47]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[24 * stride]);
			 ti = c_im(inout[24 * stride]);
			 twr = c_re(W[23]);
			 twi = c_im(W[23]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[56 * stride]);
			 ti = c_im(inout[56 * stride]);
			 twr = c_re(W[55]);
			 twi = c_im(W[55]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_3 = tre2_0_0 + tre2_1_0;
		    tim1_0_3 = tim2_0_0 + tim2_1_0;
		    tre1_1_3 = tre2_0_0 - tre2_1_0;
		    tim1_1_3 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    tre0_0_0 = tre2_0_0 + tre2_0_1;
		    tim0_0_0 = tim2_0_0 + tim2_0_1;
		    tre0_4_0 = tre2_0_0 - tre2_0_1;
		    tim0_4_0 = tim2_0_0 - tim2_0_1;
		    tre0_2_0 = tre2_1_0 - tim2_1_1;
		    tim0_2_0 = tim2_1_0 + tre2_1_1;
		    tre0_6_0 = tre2_1_0 + tim2_1_1;
		    tim0_6_0 = tim2_1_0 - tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 - tim1_1_2;
		    tim2_0_0 = tim1_1_0 + tre1_1_2;
		    tre2_1_0 = tre1_1_0 + tim1_1_2;
		    tim2_1_0 = tim1_1_0 - tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
			 tre2_0_1 = tre3_0_0 - tre3_1_0;
			 tim2_0_1 = tim3_0_0 + tim3_1_0;
			 tre2_1_1 = tre3_0_0 + tre3_1_0;
			 tim2_1_1 = tim3_0_0 - tim3_1_0;
		    }
		    tre0_1_0 = tre2_0_0 + tre2_0_1;
		    tim0_1_0 = tim2_0_0 + tim2_0_1;
		    tre0_5_0 = tre2_0_0 - tre2_0_1;
		    tim0_5_0 = tim2_0_0 - tim2_0_1;
		    tre0_3_0 = tre2_1_0 - tim2_1_1;
		    tim0_3_0 = tim2_1_0 + tre2_1_1;
		    tre0_7_0 = tre2_1_0 + tim2_1_1;
		    tim0_7_0 = tim2_1_0 - tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[stride]);
			 ti = c_im(inout[stride]);
			 twr = c_re(W[0]);
			 twi = c_im(W[0]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[33 * stride]);
			 ti = c_im(inout[33 * stride]);
			 twr = c_re(W[32]);
			 twi = c_im(W[32]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[9 * stride]);
			 ti = c_im(inout[9 * stride]);
			 twr = c_re(W[8]);
			 twi = c_im(W[8]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[41 * stride]);
			 ti = c_im(inout[41 * stride]);
			 twr = c_re(W[40]);
			 twi = c_im(W[40]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[17 * stride]);
			 ti = c_im(inout[17 * stride]);
			 twr = c_re(W[16]);
			 twi = c_im(W[16]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[49 * stride]);
			 ti = c_im(inout[49 * stride]);
			 twr = c_re(W[48]);
			 twi = c_im(W[48]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[25 * stride]);
			 ti = c_im(inout[25 * stride]);
			 twr = c_re(W[24]);
			 twi = c_im(W[24]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[57 * stride]);
			 ti = c_im(inout[57 * stride]);
			 twr = c_re(W[56]);
			 twi = c_im(W[56]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_3 = tre2_0_0 + tre2_1_0;
		    tim1_0_3 = tim2_0_0 + tim2_1_0;
		    tre1_1_3 = tre2_0_0 - tre2_1_0;
		    tim1_1_3 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    tre0_0_1 = tre2_0_0 + tre2_0_1;
		    tim0_0_1 = tim2_0_0 + tim2_0_1;
		    tre0_4_1 = tre2_0_0 - tre2_0_1;
		    tim0_4_1 = tim2_0_0 - tim2_0_1;
		    tre0_2_1 = tre2_1_0 - tim2_1_1;
		    tim0_2_1 = tim2_1_0 + tre2_1_1;
		    tre0_6_1 = tre2_1_0 + tim2_1_1;
		    tim0_6_1 = tim2_1_0 - tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 - tim1_1_2;
		    tim2_0_0 = tim1_1_0 + tre1_1_2;
		    tre2_1_0 = tre1_1_0 + tim1_1_2;
		    tim2_1_0 = tim1_1_0 - tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
			 tre2_0_1 = tre3_0_0 - tre3_1_0;
			 tim2_0_1 = tim3_0_0 + tim3_1_0;
			 tre2_1_1 = tre3_0_0 + tre3_1_0;
			 tim2_1_1 = tim3_0_0 - tim3_1_0;
		    }
		    tre0_1_1 = tre2_0_0 + tre2_0_1;
		    tim0_1_1 = tim2_0_0 + tim2_0_1;
		    tre0_5_1 = tre2_0_0 - tre2_0_1;
		    tim0_5_1 = tim2_0_0 - tim2_0_1;
		    tre0_3_1 = tre2_1_0 - tim2_1_1;
		    tim0_3_1 = tim2_1_0 + tre2_1_1;
		    tre0_7_1 = tre2_1_0 + tim2_1_1;
		    tim0_7_1 = tim2_1_0 - tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[2 * stride]);
			 ti = c_im(inout[2 * stride]);
			 twr = c_re(W[1]);
			 twi = c_im(W[1]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[34 * stride]);
			 ti = c_im(inout[34 * stride]);
			 twr = c_re(W[33]);
			 twi = c_im(W[33]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[10 * stride]);
			 ti = c_im(inout[10 * stride]);
			 twr = c_re(W[9]);
			 twi = c_im(W[9]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[42 * stride]);
			 ti = c_im(inout[42 * stride]);
			 twr = c_re(W[41]);
			 twi = c_im(W[41]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[18 * stride]);
			 ti = c_im(inout[18 * stride]);
			 twr = c_re(W[17]);
			 twi = c_im(W[17]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[50 * stride]);
			 ti = c_im(inout[50 * stride]);
			 twr = c_re(W[49]);
			 twi = c_im(W[49]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[26 * stride]);
			 ti = c_im(inout[26 * stride]);
			 twr = c_re(W[25]);
			 twi = c_im(W[25]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[58 * stride]);
			 ti = c_im(inout[58 * stride]);
			 twr = c_re(W[57]);
			 twi = c_im(W[57]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_3 = tre2_0_0 + tre2_1_0;
		    tim1_0_3 = tim2_0_0 + tim2_1_0;
		    tre1_1_3 = tre2_0_0 - tre2_1_0;
		    tim1_1_3 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    tre0_0_2 = tre2_0_0 + tre2_0_1;
		    tim0_0_2 = tim2_0_0 + tim2_0_1;
		    tre0_4_2 = tre2_0_0 - tre2_0_1;
		    tim0_4_2 = tim2_0_0 - tim2_0_1;
		    tre0_2_2 = tre2_1_0 - tim2_1_1;
		    tim0_2_2 = tim2_1_0 + tre2_1_1;
		    tre0_6_2 = tre2_1_0 + tim2_1_1;
		    tim0_6_2 = tim2_1_0 - tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 - tim1_1_2;
		    tim2_0_0 = tim1_1_0 + tre1_1_2;
		    tre2_1_0 = tre1_1_0 + tim1_1_2;
		    tim2_1_0 = tim1_1_0 - tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
			 tre2_0_1 = tre3_0_0 - tre3_1_0;
			 tim2_0_1 = tim3_0_0 + tim3_1_0;
			 tre2_1_1 = tre3_0_0 + tre3_1_0;
			 tim2_1_1 = tim3_0_0 - tim3_1_0;
		    }
		    tre0_1_2 = tre2_0_0 + tre2_0_1;
		    tim0_1_2 = tim2_0_0 + tim2_0_1;
		    tre0_5_2 = tre2_0_0 - tre2_0_1;
		    tim0_5_2 = tim2_0_0 - tim2_0_1;
		    tre0_3_2 = tre2_1_0 - tim2_1_1;
		    tim0_3_2 = tim2_1_0 + tre2_1_1;
		    tre0_7_2 = tre2_1_0 + tim2_1_1;
		    tim0_7_2 = tim2_1_0 - tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[3 * stride]);
			 ti = c_im(inout[3 * stride]);
			 twr = c_re(W[2]);
			 twi = c_im(W[2]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[35 * stride]);
			 ti = c_im(inout[35 * stride]);
			 twr = c_re(W[34]);
			 twi = c_im(W[34]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[11 * stride]);
			 ti = c_im(inout[11 * stride]);
			 twr = c_re(W[10]);
			 twi = c_im(W[10]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[43 * stride]);
			 ti = c_im(inout[43 * stride]);
			 twr = c_re(W[42]);
			 twi = c_im(W[42]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[19 * stride]);
			 ti = c_im(inout[19 * stride]);
			 twr = c_re(W[18]);
			 twi = c_im(W[18]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[51 * stride]);
			 ti = c_im(inout[51 * stride]);
			 twr = c_re(W[50]);
			 twi = c_im(W[50]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[27 * stride]);
			 ti = c_im(inout[27 * stride]);
			 twr = c_re(W[26]);
			 twi = c_im(W[26]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[59 * stride]);
			 ti = c_im(inout[59 * stride]);
			 twr = c_re(W[58]);
			 twi = c_im(W[58]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_3 = tre2_0_0 + tre2_1_0;
		    tim1_0_3 = tim2_0_0 + tim2_1_0;
		    tre1_1_3 = tre2_0_0 - tre2_1_0;
		    tim1_1_3 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    tre0_0_3 = tre2_0_0 + tre2_0_1;
		    tim0_0_3 = tim2_0_0 + tim2_0_1;
		    tre0_4_3 = tre2_0_0 - tre2_0_1;
		    tim0_4_3 = tim2_0_0 - tim2_0_1;
		    tre0_2_3 = tre2_1_0 - tim2_1_1;
		    tim0_2_3 = tim2_1_0 + tre2_1_1;
		    tre0_6_3 = tre2_1_0 + tim2_1_1;
		    tim0_6_3 = tim2_1_0 - tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 - tim1_1_2;
		    tim2_0_0 = tim1_1_0 + tre1_1_2;
		    tre2_1_0 = tre1_1_0 + tim1_1_2;
		    tim2_1_0 = tim1_1_0 - tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
			 tre2_0_1 = tre3_0_0 - tre3_1_0;
			 tim2_0_1 = tim3_0_0 + tim3_1_0;
			 tre2_1_1 = tre3_0_0 + tre3_1_0;
			 tim2_1_1 = tim3_0_0 - tim3_1_0;
		    }
		    tre0_1_3 = tre2_0_0 + tre2_0_1;
		    tim0_1_3 = tim2_0_0 + tim2_0_1;
		    tre0_5_3 = tre2_0_0 - tre2_0_1;
		    tim0_5_3 = tim2_0_0 - tim2_0_1;
		    tre0_3_3 = tre2_1_0 - tim2_1_1;
		    tim0_3_3 = tim2_1_0 + tre2_1_1;
		    tre0_7_3 = tre2_1_0 + tim2_1_1;
		    tim0_7_3 = tim2_1_0 - tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[4 * stride]);
			 ti = c_im(inout[4 * stride]);
			 twr = c_re(W[3]);
			 twi = c_im(W[3]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[36 * stride]);
			 ti = c_im(inout[36 * stride]);
			 twr = c_re(W[35]);
			 twi = c_im(W[35]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[12 * stride]);
			 ti = c_im(inout[12 * stride]);
			 twr = c_re(W[11]);
			 twi = c_im(W[11]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[44 * stride]);
			 ti = c_im(inout[44 * stride]);
			 twr = c_re(W[43]);
			 twi = c_im(W[43]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[20 * stride]);
			 ti = c_im(inout[20 * stride]);
			 twr = c_re(W[19]);
			 twi = c_im(W[19]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[52 * stride]);
			 ti = c_im(inout[52 * stride]);
			 twr = c_re(W[51]);
			 twi = c_im(W[51]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[28 * stride]);
			 ti = c_im(inout[28 * stride]);
			 twr = c_re(W[27]);
			 twi = c_im(W[27]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[60 * stride]);
			 ti = c_im(inout[60 * stride]);
			 twr = c_re(W[59]);
			 twi = c_im(W[59]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_3 = tre2_0_0 + tre2_1_0;
		    tim1_0_3 = tim2_0_0 + tim2_1_0;
		    tre1_1_3 = tre2_0_0 - tre2_1_0;
		    tim1_1_3 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    tre0_0_4 = tre2_0_0 + tre2_0_1;
		    tim0_0_4 = tim2_0_0 + tim2_0_1;
		    tre0_4_4 = tre2_0_0 - tre2_0_1;
		    tim0_4_4 = tim2_0_0 - tim2_0_1;
		    tre0_2_4 = tre2_1_0 - tim2_1_1;
		    tim0_2_4 = tim2_1_0 + tre2_1_1;
		    tre0_6_4 = tre2_1_0 + tim2_1_1;
		    tim0_6_4 = tim2_1_0 - tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 - tim1_1_2;
		    tim2_0_0 = tim1_1_0 + tre1_1_2;
		    tre2_1_0 = tre1_1_0 + tim1_1_2;
		    tim2_1_0 = tim1_1_0 - tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
			 tre2_0_1 = tre3_0_0 - tre3_1_0;
			 tim2_0_1 = tim3_0_0 + tim3_1_0;
			 tre2_1_1 = tre3_0_0 + tre3_1_0;
			 tim2_1_1 = tim3_0_0 - tim3_1_0;
		    }
		    tre0_1_4 = tre2_0_0 + tre2_0_1;
		    tim0_1_4 = tim2_0_0 + tim2_0_1;
		    tre0_5_4 = tre2_0_0 - tre2_0_1;
		    tim0_5_4 = tim2_0_0 - tim2_0_1;
		    tre0_3_4 = tre2_1_0 - tim2_1_1;
		    tim0_3_4 = tim2_1_0 + tre2_1_1;
		    tre0_7_4 = tre2_1_0 + tim2_1_1;
		    tim0_7_4 = tim2_1_0 - tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[5 * stride]);
			 ti = c_im(inout[5 * stride]);
			 twr = c_re(W[4]);
			 twi = c_im(W[4]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[37 * stride]);
			 ti = c_im(inout[37 * stride]);
			 twr = c_re(W[36]);
			 twi = c_im(W[36]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[13 * stride]);
			 ti = c_im(inout[13 * stride]);
			 twr = c_re(W[12]);
			 twi = c_im(W[12]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[45 * stride]);
			 ti = c_im(inout[45 * stride]);
			 twr = c_re(W[44]);
			 twi = c_im(W[44]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[21 * stride]);
			 ti = c_im(inout[21 * stride]);
			 twr = c_re(W[20]);
			 twi = c_im(W[20]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[53 * stride]);
			 ti = c_im(inout[53 * stride]);
			 twr = c_re(W[52]);
			 twi = c_im(W[52]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[29 * stride]);
			 ti = c_im(inout[29 * stride]);
			 twr = c_re(W[28]);
			 twi = c_im(W[28]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[61 * stride]);
			 ti = c_im(inout[61 * stride]);
			 twr = c_re(W[60]);
			 twi = c_im(W[60]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_3 = tre2_0_0 + tre2_1_0;
		    tim1_0_3 = tim2_0_0 + tim2_1_0;
		    tre1_1_3 = tre2_0_0 - tre2_1_0;
		    tim1_1_3 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    tre0_0_5 = tre2_0_0 + tre2_0_1;
		    tim0_0_5 = tim2_0_0 + tim2_0_1;
		    tre0_4_5 = tre2_0_0 - tre2_0_1;
		    tim0_4_5 = tim2_0_0 - tim2_0_1;
		    tre0_2_5 = tre2_1_0 - tim2_1_1;
		    tim0_2_5 = tim2_1_0 + tre2_1_1;
		    tre0_6_5 = tre2_1_0 + tim2_1_1;
		    tim0_6_5 = tim2_1_0 - tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 - tim1_1_2;
		    tim2_0_0 = tim1_1_0 + tre1_1_2;
		    tre2_1_0 = tre1_1_0 + tim1_1_2;
		    tim2_1_0 = tim1_1_0 - tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
			 tre2_0_1 = tre3_0_0 - tre3_1_0;
			 tim2_0_1 = tim3_0_0 + tim3_1_0;
			 tre2_1_1 = tre3_0_0 + tre3_1_0;
			 tim2_1_1 = tim3_0_0 - tim3_1_0;
		    }
		    tre0_1_5 = tre2_0_0 + tre2_0_1;
		    tim0_1_5 = tim2_0_0 + tim2_0_1;
		    tre0_5_5 = tre2_0_0 - tre2_0_1;
		    tim0_5_5 = tim2_0_0 - tim2_0_1;
		    tre0_3_5 = tre2_1_0 - tim2_1_1;
		    tim0_3_5 = tim2_1_0 + tre2_1_1;
		    tre0_7_5 = tre2_1_0 + tim2_1_1;
		    tim0_7_5 = tim2_1_0 - tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[6 * stride]);
			 ti = c_im(inout[6 * stride]);
			 twr = c_re(W[5]);
			 twi = c_im(W[5]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[38 * stride]);
			 ti = c_im(inout[38 * stride]);
			 twr = c_re(W[37]);
			 twi = c_im(W[37]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[14 * stride]);
			 ti = c_im(inout[14 * stride]);
			 twr = c_re(W[13]);
			 twi = c_im(W[13]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[46 * stride]);
			 ti = c_im(inout[46 * stride]);
			 twr = c_re(W[45]);
			 twi = c_im(W[45]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[22 * stride]);
			 ti = c_im(inout[22 * stride]);
			 twr = c_re(W[21]);
			 twi = c_im(W[21]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[54 * stride]);
			 ti = c_im(inout[54 * stride]);
			 twr = c_re(W[53]);
			 twi = c_im(W[53]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[30 * stride]);
			 ti = c_im(inout[30 * stride]);
			 twr = c_re(W[29]);
			 twi = c_im(W[29]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[62 * stride]);
			 ti = c_im(inout[62 * stride]);
			 twr = c_re(W[61]);
			 twi = c_im(W[61]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_3 = tre2_0_0 + tre2_1_0;
		    tim1_0_3 = tim2_0_0 + tim2_1_0;
		    tre1_1_3 = tre2_0_0 - tre2_1_0;
		    tim1_1_3 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    tre0_0_6 = tre2_0_0 + tre2_0_1;
		    tim0_0_6 = tim2_0_0 + tim2_0_1;
		    tre0_4_6 = tre2_0_0 - tre2_0_1;
		    tim0_4_6 = tim2_0_0 - tim2_0_1;
		    tre0_2_6 = tre2_1_0 - tim2_1_1;
		    tim0_2_6 = tim2_1_0 + tre2_1_1;
		    tre0_6_6 = tre2_1_0 + tim2_1_1;
		    tim0_6_6 = tim2_1_0 - tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 - tim1_1_2;
		    tim2_0_0 = tim1_1_0 + tre1_1_2;
		    tre2_1_0 = tre1_1_0 + tim1_1_2;
		    tim2_1_0 = tim1_1_0 - tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
			 tre2_0_1 = tre3_0_0 - tre3_1_0;
			 tim2_0_1 = tim3_0_0 + tim3_1_0;
			 tre2_1_1 = tre3_0_0 + tre3_1_0;
			 tim2_1_1 = tim3_0_0 - tim3_1_0;
		    }
		    tre0_1_6 = tre2_0_0 + tre2_0_1;
		    tim0_1_6 = tim2_0_0 + tim2_0_1;
		    tre0_5_6 = tre2_0_0 - tre2_0_1;
		    tim0_5_6 = tim2_0_0 - tim2_0_1;
		    tre0_3_6 = tre2_1_0 - tim2_1_1;
		    tim0_3_6 = tim2_1_0 + tre2_1_1;
		    tre0_7_6 = tre2_1_0 + tim2_1_1;
		    tim0_7_6 = tim2_1_0 - tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[7 * stride]);
			 ti = c_im(inout[7 * stride]);
			 twr = c_re(W[6]);
			 twi = c_im(W[6]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[39 * stride]);
			 ti = c_im(inout[39 * stride]);
			 twr = c_re(W[38]);
			 twi = c_im(W[38]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_0 = tre2_0_0 + tre2_1_0;
		    tim1_0_0 = tim2_0_0 + tim2_1_0;
		    tre1_1_0 = tre2_0_0 - tre2_1_0;
		    tim1_1_0 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[15 * stride]);
			 ti = c_im(inout[15 * stride]);
			 twr = c_re(W[14]);
			 twi = c_im(W[14]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[47 * stride]);
			 ti = c_im(inout[47 * stride]);
			 twr = c_re(W[46]);
			 twi = c_im(W[46]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[23 * stride]);
			 ti = c_im(inout[23 * stride]);
			 twr = c_re(W[22]);
			 twi = c_im(W[22]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[55 * stride]);
			 ti = c_im(inout[55 * stride]);
			 twr = c_re(W[54]);
			 twi = c_im(W[54]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[31 * stride]);
			 ti = c_im(inout[31 * stride]);
			 twr = c_re(W[30]);
			 twi = c_im(W[30]);
			 tre2_0_0 = (tr * twr) + (ti * twi);
			 tim2_0_0 = (ti * twr) - (tr * twi);
		    }
		    {
			 FFTW_REAL tr;
			 FFTW_REAL ti;
			 FFTW_REAL twr;
			 FFTW_REAL twi;
			 tr = c_re(inout[63 * stride]);
			 ti = c_im(inout[63 * stride]);
			 twr = c_re(W[62]);
			 twi = c_im(W[62]);
			 tre2_1_0 = (tr * twr) + (ti * twi);
			 tim2_1_0 = (ti * twr) - (tr * twi);
		    }
		    tre1_0_3 = tre2_0_0 + tre2_1_0;
		    tim1_0_3 = tim2_0_0 + tim2_1_0;
		    tre1_1_3 = tre2_0_0 - tre2_1_0;
		    tim1_1_3 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    tre0_0_7 = tre2_0_0 + tre2_0_1;
		    tim0_0_7 = tim2_0_0 + tim2_0_1;
		    tre0_4_7 = tre2_0_0 - tre2_0_1;
		    tim0_4_7 = tim2_0_0 - tim2_0_1;
		    tre0_2_7 = tre2_1_0 - tim2_1_1;
		    tim0_2_7 = tim2_1_0 + tre2_1_1;
		    tre0_6_7 = tre2_1_0 + tim2_1_1;
		    tim0_6_7 = tim2_1_0 - tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 - tim1_1_2;
		    tim2_0_0 = tim1_1_0 + tre1_1_2;
		    tre2_1_0 = tre1_1_0 + tim1_1_2;
		    tim2_1_0 = tim1_1_0 - tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
			 tre2_0_1 = tre3_0_0 - tre3_1_0;
			 tim2_0_1 = tim3_0_0 + tim3_1_0;
			 tre2_1_1 = tre3_0_0 + tre3_1_0;
			 tim2_1_1 = tim3_0_0 - tim3_1_0;
		    }
		    tre0_1_7 = tre2_0_0 + tre2_0_1;
		    tim0_1_7 = tim2_0_0 + tim2_0_1;
		    tre0_5_7 = tre2_0_0 - tre2_0_1;
		    tim0_5_7 = tim2_0_0 - tim2_0_1;
		    tre0_3_7 = tre2_1_0 - tim2_1_1;
		    tim0_3_7 = tim2_1_0 + tre2_1_1;
		    tre0_7_7 = tre2_1_0 + tim2_1_1;
		    tim0_7_7 = tim2_1_0 - tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       tre1_0_0 = tre0_0_0 + tre0_0_4;
	       tim1_0_0 = tim0_0_0 + tim0_0_4;
	       tre1_1_0 = tre0_0_0 - tre0_0_4;
	       tim1_1_0 = tim0_0_0 - tim0_0_4;
	       tre1_0_1 = tre0_0_1 + tre0_0_5;
	       tim1_0_1 = tim0_0_1 + tim0_0_5;
	       tre1_1_1 = tre0_0_1 - tre0_0_5;
	       tim1_1_1 = tim0_0_1 - tim0_0_5;
	       tre1_0_2 = tre0_0_2 + tre0_0_6;
	       tim1_0_2 = tim0_0_2 + tim0_0_6;
	       tre1_1_2 = tre0_0_2 - tre0_0_6;
	       tim1_1_2 = tim0_0_2 - tim0_0_6;
	       tre1_0_3 = tre0_0_3 + tre0_0_7;
	       tim1_0_3 = tim0_0_3 + tim0_0_7;
	       tre1_1_3 = tre0_0_3 - tre0_0_7;
	       tim1_1_3 = tim0_0_3 - tim0_0_7;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    c_re(inout[0]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[0]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[32 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[32 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[16 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[16 * stride]) = tim2_1_0 + tre2_1_1;
		    c_re(inout[48 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[48 * stride]) = tim2_1_0 - tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 - tim1_1_2;
		    tim2_0_0 = tim1_1_0 + tre1_1_2;
		    tre2_1_0 = tre1_1_0 + tim1_1_2;
		    tim2_1_0 = tim1_1_0 - tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
			 tre2_0_1 = tre3_0_0 - tre3_1_0;
			 tim2_0_1 = tim3_0_0 + tim3_1_0;
			 tre2_1_1 = tre3_0_0 + tre3_1_0;
			 tim2_1_1 = tim3_0_0 - tim3_1_0;
		    }
		    c_re(inout[8 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[8 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[40 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[40 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[24 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[24 * stride]) = tim2_1_0 + tre2_1_1;
		    c_re(inout[56 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[56 * stride]) = tim2_1_0 - tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_1_4) - (((FFTW_REAL) FFTW_K382683432) * tim0_1_4);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_1_4) + (((FFTW_REAL) FFTW_K382683432) * tre0_1_4);
		    tre1_0_0 = tre0_1_0 + tre2_1_0;
		    tim1_0_0 = tim0_1_0 + tim2_1_0;
		    tre1_1_0 = tre0_1_0 - tre2_1_0;
		    tim1_1_0 = tim0_1_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K995184726) * tre0_1_1) - (((FFTW_REAL) FFTW_K098017140) * tim0_1_1);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K995184726) * tim0_1_1) + (((FFTW_REAL) FFTW_K098017140) * tre0_1_1);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K881921264) * tre0_1_5) - (((FFTW_REAL) FFTW_K471396736) * tim0_1_5);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K881921264) * tim0_1_5) + (((FFTW_REAL) FFTW_K471396736) * tre0_1_5);
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tre0_1_2) - (((FFTW_REAL) FFTW_K195090322) * tim0_1_2);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tim0_1_2) + (((FFTW_REAL) FFTW_K195090322) * tre0_1_2);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_1_6) - (((FFTW_REAL) FFTW_K555570233) * tim0_1_6);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_1_6) + (((FFTW_REAL) FFTW_K555570233) * tre0_1_6);
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K956940335) * tre0_1_3) - (((FFTW_REAL) FFTW_K290284677) * tim0_1_3);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K956940335) * tim0_1_3) + (((FFTW_REAL) FFTW_K290284677) * tre0_1_3);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K773010453) * tre0_1_7) - (((FFTW_REAL) FFTW_K634393284) * tim0_1_7);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K773010453) * tim0_1_7) + (((FFTW_REAL) FFTW_K634393284) * tre0_1_7);
		    tre1_0_3 = tre2_0_0 + tre2_1_0;
		    tim1_0_3 = tim2_0_0 + tim2_1_0;
		    tre1_1_3 = tre2_0_0 - tre2_1_0;
		    tim1_1_3 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    c_re(inout[stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[33 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[33 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[17 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[17 * stride]) = tim2_1_0 + tre2_1_1;
		    c_re(inout[49 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[49 * stride]) = tim2_1_0 - tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 - tim1_1_2;
		    tim2_0_0 = tim1_1_0 + tre1_1_2;
		    tre2_1_0 = tre1_1_0 + tim1_1_2;
		    tim2_1_0 = tim1_1_0 - tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
			 tre2_0_1 = tre3_0_0 - tre3_1_0;
			 tim2_0_1 = tim3_0_0 + tim3_1_0;
			 tre2_1_1 = tre3_0_0 + tre3_1_0;
			 tim2_1_1 = tim3_0_0 - tim3_1_0;
		    }
		    c_re(inout[9 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[9 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[41 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[41 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[25 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[25 * stride]) = tim2_1_0 + tre2_1_1;
		    c_re(inout[57 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[57 * stride]) = tim2_1_0 - tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_2_4 - tim0_2_4);
		    tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_2_4 + tre0_2_4);
		    tre1_0_0 = tre0_2_0 + tre2_1_0;
		    tim1_0_0 = tim0_2_0 + tim2_1_0;
		    tre1_1_0 = tre0_2_0 - tre2_1_0;
		    tim1_1_0 = tim0_2_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tre0_2_1) - (((FFTW_REAL) FFTW_K195090322) * tim0_2_1);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tim0_2_1) + (((FFTW_REAL) FFTW_K195090322) * tre0_2_1);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K555570233) * tre0_2_5) - (((FFTW_REAL) FFTW_K831469612) * tim0_2_5);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K555570233) * tim0_2_5) + (((FFTW_REAL) FFTW_K831469612) * tre0_2_5);
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_2_2) - (((FFTW_REAL) FFTW_K382683432) * tim0_2_2);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_2_2) + (((FFTW_REAL) FFTW_K382683432) * tre0_2_2);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_2_6) - (((FFTW_REAL) FFTW_K923879532) * tim0_2_6);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_2_6) + (((FFTW_REAL) FFTW_K923879532) * tre0_2_6);
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_2_3) - (((FFTW_REAL) FFTW_K555570233) * tim0_2_3);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_2_3) + (((FFTW_REAL) FFTW_K555570233) * tre0_2_3);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tre0_2_7) - (((FFTW_REAL) FFTW_K980785280) * tim0_2_7);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tim0_2_7) + (((FFTW_REAL) FFTW_K980785280) * tre0_2_7);
		    tre1_0_3 = tre2_0_0 + tre2_1_0;
		    tim1_0_3 = tim2_0_0 + tim2_1_0;
		    tre1_1_3 = tre2_0_0 - tre2_1_0;
		    tim1_1_3 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    c_re(inout[2 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[2 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[34 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[34 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[18 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[18 * stride]) = tim2_1_0 + tre2_1_1;
		    c_re(inout[50 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[50 * stride]) = tim2_1_0 - tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 - tim1_1_2;
		    tim2_0_0 = tim1_1_0 + tre1_1_2;
		    tre2_1_0 = tre1_1_0 + tim1_1_2;
		    tim2_1_0 = tim1_1_0 - tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
			 tre2_0_1 = tre3_0_0 - tre3_1_0;
			 tim2_0_1 = tim3_0_0 + tim3_1_0;
			 tre2_1_1 = tre3_0_0 + tre3_1_0;
			 tim2_1_1 = tim3_0_0 - tim3_1_0;
		    }
		    c_re(inout[10 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[10 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[42 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[42 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[26 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[26 * stride]) = tim2_1_0 + tre2_1_1;
		    c_re(inout[58 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[58 * stride]) = tim2_1_0 - tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_3_4) - (((FFTW_REAL) FFTW_K923879532) * tim0_3_4);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_3_4) + (((FFTW_REAL) FFTW_K923879532) * tre0_3_4);
		    tre1_0_0 = tre0_3_0 + tre2_1_0;
		    tim1_0_0 = tim0_3_0 + tim2_1_0;
		    tre1_1_0 = tre0_3_0 - tre2_1_0;
		    tim1_1_0 = tim0_3_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K956940335) * tre0_3_1) - (((FFTW_REAL) FFTW_K290284677) * tim0_3_1);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K956940335) * tim0_3_1) + (((FFTW_REAL) FFTW_K290284677) * tre0_3_1);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K098017140) * tre0_3_5) - (((FFTW_REAL) FFTW_K995184726) * tim0_3_5);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K098017140) * tim0_3_5) + (((FFTW_REAL) FFTW_K995184726) * tre0_3_5);
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_3_2) - (((FFTW_REAL) FFTW_K555570233) * tim0_3_2);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_3_2) + (((FFTW_REAL) FFTW_K555570233) * tre0_3_2);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tre0_3_6) + (((FFTW_REAL) FFTW_K980785280) * tim0_3_6);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K980785280) * tre0_3_6) - (((FFTW_REAL) FFTW_K195090322) * tim0_3_6);
		    tre1_0_2 = tre2_0_0 - tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 + tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K634393284) * tre0_3_3) - (((FFTW_REAL) FFTW_K773010453) * tim0_3_3);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K634393284) * tim0_3_3) + (((FFTW_REAL) FFTW_K773010453) * tre0_3_3);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K471396736) * tre0_3_7) + (((FFTW_REAL) FFTW_K881921264) * tim0_3_7);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K881921264) * tre0_3_7) - (((FFTW_REAL) FFTW_K471396736) * tim0_3_7);
		    tre1_0_3 = tre2_0_0 - tre2_1_0;
		    tim1_0_3 = tim2_0_0 + tim2_1_0;
		    tre1_1_3 = tre2_0_0 + tre2_1_0;
		    tim1_1_3 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    c_re(inout[3 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[3 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[35 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[35 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[19 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[19 * stride]) = tim2_1_0 + tre2_1_1;
		    c_re(inout[51 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[51 * stride]) = tim2_1_0 - tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 - tim1_1_2;
		    tim2_0_0 = tim1_1_0 + tre1_1_2;
		    tre2_1_0 = tre1_1_0 + tim1_1_2;
		    tim2_1_0 = tim1_1_0 - tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
			 tre2_0_1 = tre3_0_0 - tre3_1_0;
			 tim2_0_1 = tim3_0_0 + tim3_1_0;
			 tre2_1_1 = tre3_0_0 + tre3_1_0;
			 tim2_1_1 = tim3_0_0 - tim3_1_0;
		    }
		    c_re(inout[11 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[11 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[43 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[43 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[27 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[27 * stride]) = tim2_1_0 + tre2_1_1;
		    c_re(inout[59 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[59 * stride]) = tim2_1_0 - tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       tre1_0_0 = tre0_4_0 - tim0_4_4;
	       tim1_0_0 = tim0_4_0 + tre0_4_4;
	       tre1_1_0 = tre0_4_0 + tim0_4_4;
	       tim1_1_0 = tim0_4_0 - tre0_4_4;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_4_1) - (((FFTW_REAL) FFTW_K382683432) * tim0_4_1);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_4_1) + (((FFTW_REAL) FFTW_K382683432) * tre0_4_1);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_4_5) + (((FFTW_REAL) FFTW_K923879532) * tim0_4_5);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_4_5) - (((FFTW_REAL) FFTW_K382683432) * tim0_4_5);
		    tre1_0_1 = tre2_0_0 - tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 + tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_4_2 - tim0_4_2);
		    tim2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_4_2 + tre0_4_2);
		    tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_4_6 + tim0_4_6);
		    tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_4_6 - tim0_4_6);
		    tre1_0_2 = tre2_0_0 - tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 + tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_4_3) - (((FFTW_REAL) FFTW_K923879532) * tim0_4_3);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_4_3) + (((FFTW_REAL) FFTW_K923879532) * tre0_4_3);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_4_7) + (((FFTW_REAL) FFTW_K382683432) * tim0_4_7);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_4_7) - (((FFTW_REAL) FFTW_K923879532) * tim0_4_7);
		    tre1_0_3 = tre2_0_0 - tre2_1_0;
		    tim1_0_3 = tim2_0_0 + tim2_1_0;
		    tre1_1_3 = tre2_0_0 + tre2_1_0;
		    tim1_1_3 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    c_re(inout[4 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[4 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[36 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[36 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[20 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[20 * stride]) = tim2_1_0 + tre2_1_1;
		    c_re(inout[52 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[52 * stride]) = tim2_1_0 - tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 - tim1_1_2;
		    tim2_0_0 = tim1_1_0 + tre1_1_2;
		    tre2_1_0 = tre1_1_0 + tim1_1_2;
		    tim2_1_0 = tim1_1_0 - tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
			 tre2_0_1 = tre3_0_0 - tre3_1_0;
			 tim2_0_1 = tim3_0_0 + tim3_1_0;
			 tre2_1_1 = tre3_0_0 + tre3_1_0;
			 tim2_1_1 = tim3_0_0 - tim3_1_0;
		    }
		    c_re(inout[12 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[12 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[44 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[44 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[28 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[28 * stride]) = tim2_1_0 + tre2_1_1;
		    c_re(inout[60 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[60 * stride]) = tim2_1_0 - tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_5_4) + (((FFTW_REAL) FFTW_K923879532) * tim0_5_4);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_5_4) - (((FFTW_REAL) FFTW_K382683432) * tim0_5_4);
		    tre1_0_0 = tre0_5_0 - tre2_1_0;
		    tim1_0_0 = tim0_5_0 + tim2_1_0;
		    tre1_1_0 = tre0_5_0 + tre2_1_0;
		    tim1_1_0 = tim0_5_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K881921264) * tre0_5_1) - (((FFTW_REAL) FFTW_K471396736) * tim0_5_1);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K881921264) * tim0_5_1) + (((FFTW_REAL) FFTW_K471396736) * tre0_5_1);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K773010453) * tre0_5_5) + (((FFTW_REAL) FFTW_K634393284) * tim0_5_5);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K634393284) * tre0_5_5) - (((FFTW_REAL) FFTW_K773010453) * tim0_5_5);
		    tre1_0_1 = tre2_0_0 - tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 + tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K555570233) * tre0_5_2) - (((FFTW_REAL) FFTW_K831469612) * tim0_5_2);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K555570233) * tim0_5_2) + (((FFTW_REAL) FFTW_K831469612) * tre0_5_2);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K980785280) * tre0_5_6) + (((FFTW_REAL) FFTW_K195090322) * tim0_5_6);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tre0_5_6) - (((FFTW_REAL) FFTW_K980785280) * tim0_5_6);
		    tre1_0_2 = tre2_0_0 - tre2_1_0;
		    tim1_0_2 = tim2_0_0 + tim2_1_0;
		    tre1_1_2 = tre2_0_0 + tre2_1_0;
		    tim1_1_2 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K098017140) * tre0_5_3) - (((FFTW_REAL) FFTW_K995184726) * tim0_5_3);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K098017140) * tim0_5_3) + (((FFTW_REAL) FFTW_K995184726) * tre0_5_3);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K290284677) * tim0_5_7) - (((FFTW_REAL) FFTW_K956940335) * tre0_5_7);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K956940335) * tim0_5_7) + (((FFTW_REAL) FFTW_K290284677) * tre0_5_7);
		    tre1_0_3 = tre2_0_0 + tre2_1_0;
		    tim1_0_3 = tim2_0_0 - tim2_1_0;
		    tre1_1_3 = tre2_0_0 - tre2_1_0;
		    tim1_1_3 = tim2_0_0 + tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    c_re(inout[5 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[5 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[37 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[37 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[21 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[21 * stride]) = tim2_1_0 + tre2_1_1;
		    c_re(inout[53 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[53 * stride]) = tim2_1_0 - tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 - tim1_1_2;
		    tim2_0_0 = tim1_1_0 + tre1_1_2;
		    tre2_1_0 = tre1_1_0 + tim1_1_2;
		    tim2_1_0 = tim1_1_0 - tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
			 tre2_0_1 = tre3_0_0 - tre3_1_0;
			 tim2_0_1 = tim3_0_0 + tim3_1_0;
			 tre2_1_1 = tre3_0_0 + tre3_1_0;
			 tim2_1_1 = tim3_0_0 - tim3_1_0;
		    }
		    c_re(inout[13 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[13 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[45 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[45 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[29 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[29 * stride]) = tim2_1_0 + tre2_1_1;
		    c_re(inout[61 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[61 * stride]) = tim2_1_0 - tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_6_4 + tim0_6_4);
		    tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_6_4 - tim0_6_4);
		    tre1_0_0 = tre0_6_0 - tre2_1_0;
		    tim1_0_0 = tim0_6_0 + tim2_1_0;
		    tre1_1_0 = tre0_6_0 + tre2_1_0;
		    tim1_1_0 = tim0_6_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tre0_6_1) - (((FFTW_REAL) FFTW_K555570233) * tim0_6_1);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_6_1) + (((FFTW_REAL) FFTW_K555570233) * tre0_6_1);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K980785280) * tre0_6_5) + (((FFTW_REAL) FFTW_K195090322) * tim0_6_5);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K195090322) * tre0_6_5) - (((FFTW_REAL) FFTW_K980785280) * tim0_6_5);
		    tre1_0_1 = tre2_0_0 - tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 + tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_6_2) - (((FFTW_REAL) FFTW_K923879532) * tim0_6_2);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_6_2) + (((FFTW_REAL) FFTW_K923879532) * tre0_6_2);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tim0_6_6) - (((FFTW_REAL) FFTW_K923879532) * tre0_6_6);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tim0_6_6) + (((FFTW_REAL) FFTW_K382683432) * tre0_6_6);
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 - tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 + tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K195090322) * tre0_6_3) + (((FFTW_REAL) FFTW_K980785280) * tim0_6_3);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K980785280) * tre0_6_3) - (((FFTW_REAL) FFTW_K195090322) * tim0_6_3);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_6_7) - (((FFTW_REAL) FFTW_K555570233) * tre0_6_7);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K555570233) * tim0_6_7) + (((FFTW_REAL) FFTW_K831469612) * tre0_6_7);
		    tre1_0_3 = tre2_1_0 - tre2_0_0;
		    tim1_0_3 = tim2_0_0 - tim2_1_0;
		    tre1_1_3 = (-(tre2_0_0 + tre2_1_0));
		    tim1_1_3 = tim2_0_0 + tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    c_re(inout[6 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[6 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[38 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[38 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[22 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[22 * stride]) = tim2_1_0 + tre2_1_1;
		    c_re(inout[54 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[54 * stride]) = tim2_1_0 - tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 - tim1_1_2;
		    tim2_0_0 = tim1_1_0 + tre1_1_2;
		    tre2_1_0 = tre1_1_0 + tim1_1_2;
		    tim2_1_0 = tim1_1_0 - tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
			 tre2_0_1 = tre3_0_0 - tre3_1_0;
			 tim2_0_1 = tim3_0_0 + tim3_1_0;
			 tre2_1_1 = tre3_0_0 + tre3_1_0;
			 tim2_1_1 = tim3_0_0 - tim3_1_0;
		    }
		    c_re(inout[14 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[14 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[46 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[46 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[30 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[30 * stride]) = tim2_1_0 + tre2_1_1;
		    c_re(inout[62 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[62 * stride]) = tim2_1_0 - tre2_1_1;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_0_2;
	       FFTW_REAL tim1_0_2;
	       FFTW_REAL tre1_0_3;
	       FFTW_REAL tim1_0_3;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       FFTW_REAL tre1_1_2;
	       FFTW_REAL tim1_1_2;
	       FFTW_REAL tre1_1_3;
	       FFTW_REAL tim1_1_3;
	       {
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_1_0 = (((FFTW_REAL) FFTW_K923879532) * tre0_7_4) + (((FFTW_REAL) FFTW_K382683432) * tim0_7_4);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K382683432) * tre0_7_4) - (((FFTW_REAL) FFTW_K923879532) * tim0_7_4);
		    tre1_0_0 = tre0_7_0 - tre2_1_0;
		    tim1_0_0 = tim0_7_0 + tim2_1_0;
		    tre1_1_0 = tre0_7_0 + tre2_1_0;
		    tim1_1_0 = tim0_7_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K773010453) * tre0_7_1) - (((FFTW_REAL) FFTW_K634393284) * tim0_7_1);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K773010453) * tim0_7_1) + (((FFTW_REAL) FFTW_K634393284) * tre0_7_1);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K290284677) * tim0_7_5) - (((FFTW_REAL) FFTW_K956940335) * tre0_7_5);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K956940335) * tim0_7_5) + (((FFTW_REAL) FFTW_K290284677) * tre0_7_5);
		    tre1_0_1 = tre2_0_0 + tre2_1_0;
		    tim1_0_1 = tim2_0_0 - tim2_1_0;
		    tre1_1_1 = tre2_0_0 - tre2_1_0;
		    tim1_1_1 = tim2_0_0 + tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K195090322) * tre0_7_2) - (((FFTW_REAL) FFTW_K980785280) * tim0_7_2);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K195090322) * tim0_7_2) + (((FFTW_REAL) FFTW_K980785280) * tre0_7_2);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K831469612) * tim0_7_6) - (((FFTW_REAL) FFTW_K555570233) * tre0_7_6);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K555570233) * tim0_7_6) + (((FFTW_REAL) FFTW_K831469612) * tre0_7_6);
		    tre1_0_2 = tre2_0_0 + tre2_1_0;
		    tim1_0_2 = tim2_0_0 - tim2_1_0;
		    tre1_1_2 = tre2_0_0 - tre2_1_0;
		    tim1_1_2 = tim2_0_0 + tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = (((FFTW_REAL) FFTW_K471396736) * tre0_7_3) + (((FFTW_REAL) FFTW_K881921264) * tim0_7_3);
		    tim2_0_0 = (((FFTW_REAL) FFTW_K881921264) * tre0_7_3) - (((FFTW_REAL) FFTW_K471396736) * tim0_7_3);
		    tre2_1_0 = (((FFTW_REAL) FFTW_K098017140) * tre0_7_7) + (((FFTW_REAL) FFTW_K995184726) * tim0_7_7);
		    tim2_1_0 = (((FFTW_REAL) FFTW_K098017140) * tim0_7_7) - (((FFTW_REAL) FFTW_K995184726) * tre0_7_7);
		    tre1_0_3 = tre2_1_0 - tre2_0_0;
		    tim1_0_3 = tim2_0_0 + tim2_1_0;
		    tre1_1_3 = (-(tre2_0_0 + tre2_1_0));
		    tim1_1_3 = tim2_0_0 - tim2_1_0;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_0_0 + tre1_0_2;
		    tim2_0_0 = tim1_0_0 + tim1_0_2;
		    tre2_1_0 = tre1_0_0 - tre1_0_2;
		    tim2_1_0 = tim1_0_0 - tim1_0_2;
		    tre2_0_1 = tre1_0_1 + tre1_0_3;
		    tim2_0_1 = tim1_0_1 + tim1_0_3;
		    tre2_1_1 = tre1_0_1 - tre1_0_3;
		    tim2_1_1 = tim1_0_1 - tim1_0_3;
		    c_re(inout[7 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[7 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[39 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[39 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[23 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[23 * stride]) = tim2_1_0 + tre2_1_1;
		    c_re(inout[55 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[55 * stride]) = tim2_1_0 - tre2_1_1;
	       }
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_0_1;
		    FFTW_REAL tim2_0_1;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    FFTW_REAL tre2_1_1;
		    FFTW_REAL tim2_1_1;
		    tre2_0_0 = tre1_1_0 - tim1_1_2;
		    tim2_0_0 = tim1_1_0 + tre1_1_2;
		    tre2_1_0 = tre1_1_0 + tim1_1_2;
		    tim2_1_0 = tim1_1_0 - tre1_1_2;
		    {
			 FFTW_REAL tre3_0_0;
			 FFTW_REAL tim3_0_0;
			 FFTW_REAL tre3_1_0;
			 FFTW_REAL tim3_1_0;
			 tre3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_1 - tim1_1_1);
			 tim3_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim1_1_1 + tre1_1_1);
			 tre3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 + tim1_1_3);
			 tim3_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre1_1_3 - tim1_1_3);
			 tre2_0_1 = tre3_0_0 - tre3_1_0;
			 tim2_0_1 = tim3_0_0 + tim3_1_0;
			 tre2_1_1 = tre3_0_0 + tre3_1_0;
			 tim2_1_1 = tim3_0_0 - tim3_1_0;
		    }
		    c_re(inout[15 * stride]) = tre2_0_0 + tre2_0_1;
		    c_im(inout[15 * stride]) = tim2_0_0 + tim2_0_1;
		    c_re(inout[47 * stride]) = tre2_0_0 - tre2_0_1;
		    c_im(inout[47 * stride]) = tim2_0_0 - tim2_0_1;
		    c_re(inout[31 * stride]) = tre2_1_0 - tim2_1_1;
		    c_im(inout[31 * stride]) = tim2_1_0 + tre2_1_1;
		    c_re(inout[63 * stride]) = tre2_1_0 + tim2_1_1;
		    c_im(inout[63 * stride]) = tim2_1_0 - tre2_1_1;
	       }
	  }
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 102 FP additions and 60 FP multiplications */

void fftwi_twiddle_7(FFTW_COMPLEX *A, const FFTW_COMPLEX *W, int stride, int m, int dist)
{
     int i;
     COMPLEX *inout;
     inout = A;
     for (i = 0; i < m; i = i + 1, inout = inout + dist, W = W + 6) {
	  FFTW_REAL tre0_0_0;
	  FFTW_REAL tim0_0_0;
	  FFTW_REAL tre0_1_0;
	  FFTW_REAL tim0_1_0;
	  FFTW_REAL tre0_2_0;
	  FFTW_REAL tim0_2_0;
	  FFTW_REAL tre0_3_0;
	  FFTW_REAL tim0_3_0;
	  FFTW_REAL tre0_4_0;
	  FFTW_REAL tim0_4_0;
	  FFTW_REAL tre0_5_0;
	  FFTW_REAL tim0_5_0;
	  FFTW_REAL tre0_6_0;
	  FFTW_REAL tim0_6_0;
	  tre0_0_0 = c_re(inout[0]);
	  tim0_0_0 = c_im(inout[0]);
	  {
	       FFTW_REAL tr;
	       FFTW_REAL ti;
	       FFTW_REAL twr;
	       FFTW_REAL twi;
	       tr = c_re(inout[stride]);
	       ti = c_im(inout[stride]);
	       twr = c_re(W[0]);
	       twi = c_im(W[0]);
	       tre0_1_0 = (tr * twr) + (ti * twi);
	       tim0_1_0 = (ti * twr) - (tr * twi);
	  }
	  {
	       FFTW_REAL tr;
	       FFTW_REAL ti;
	       FFTW_REAL twr;
	       FFTW_REAL twi;
	       tr = c_re(inout[2 * stride]);
	       ti = c_im(inout[2 * stride]);
	       twr = c_re(W[1]);
	       twi = c_im(W[1]);
	       tre0_2_0 = (tr * twr) + (ti * twi);
	       tim0_2_0 = (ti * twr) - (tr * twi);
	  }
	  {
	       FFTW_REAL tr;
	       FFTW_REAL ti;
	       FFTW_REAL twr;
	       FFTW_REAL twi;
	       tr = c_re(inout[3 * stride]);
	       ti = c_im(inout[3 * stride]);
	       twr = c_re(W[2]);
	       twi = c_im(W[2]);
	       tre0_3_0 = (tr * twr) + (ti * twi);
	       tim0_3_0 = (ti * twr) - (tr * twi);
	  }
	  {
	       FFTW_REAL tr;
	       FFTW_REAL ti;
	       FFTW_REAL twr;
	       FFTW_REAL twi;
	       tr = c_re(inout[4 * stride]);
	       ti = c_im(inout[4 * stride]);
	       twr = c_re(W[3]);
	       twi = c_im(W[3]);
	       tre0_4_0 = (tr * twr) + (ti * twi);
	       tim0_4_0 = (ti * twr) - (tr * twi);
	  }
	  {
	       FFTW_REAL tr;
	       FFTW_REAL ti;
	       FFTW_REAL twr;
	       FFTW_REAL twi;
	       tr = c_re(inout[5 * stride]);
	       ti = c_im(inout[5 * stride]);
	       twr = c_re(W[4]);
	       twi = c_im(W[4]);
	       tre0_5_0 = (tr * twr) + (ti * twi);
	       tim0_5_0 = (ti * twr) - (tr * twi);
	  }
	  {
	       FFTW_REAL tr;
	       FFTW_REAL ti;
	       FFTW_REAL twr;
	       FFTW_REAL twi;
	       tr = c_re(inout[6 * stride]);
	       ti = c_im(inout[6 * stride]);
	       twr = c_re(W[5]);
	       twi = c_im(W[5]);
	       tre0_6_0 = (tr * twr) + (ti * twi);
	       tim0_6_0 = (ti * twr) - (tr * twi);
	  }
	  c_re(inout[0]) = tre0_0_0 + tre0_1_0 + tre0_2_0 + tre0_3_0 + tre0_4_0 + tre0_5_0 + tre0_6_0;
	  c_im(inout[0]) = tim0_0_0 + tim0_1_0 + tim0_2_0 + tim0_3_0 + tim0_4_0 + tim0_5_0 + tim0_6_0;
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tre1_1_0;
	       tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tre0_1_0 + tre0_6_0)) - (((FFTW_REAL) FFTW_K900968867) * (tre0_3_0 + tre0_4_0)) - (((FFTW_REAL) FFTW_K222520933) * (tre0_2_0 + tre0_5_0));
	       tre1_1_0 = (((FFTW_REAL) FFTW_K781831482) * (tim0_6_0 - tim0_1_0)) + (((FFTW_REAL) FFTW_K974927912) * (tim0_5_0 - tim0_2_0)) + (((FFTW_REAL) FFTW_K433883739) * (tim0_4_0 - tim0_3_0));
	       c_re(inout[stride]) = tre1_0_0 + tre1_1_0;
	       c_re(inout[6 * stride]) = tre1_0_0 - tre1_1_0;
	  }
	  {
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tim1_1_0;
	       tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tim0_1_0 + tim0_6_0)) - (((FFTW_REAL) FFTW_K900968867) * (tim0_3_0 + tim0_4_0)) - (((FFTW_REAL) FFTW_K222520933) * (tim0_2_0 + tim0_5_0));
	       tim1_1_0 = (((FFTW_REAL) FFTW_K781831482) * (tre0_1_0 - tre0_6_0)) + (((FFTW_REAL) FFTW_K974927912) * (tre0_2_0 - tre0_5_0)) + (((FFTW_REAL) FFTW_K433883739) * (tre0_3_0 - tre0_4_0));
	       c_im(inout[stride]) = tim1_0_0 + tim1_1_0;
	       c_im(inout[6 * stride]) = tim1_0_0 - tim1_1_0;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tre1_1_0;
	       tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tre0_3_0 + tre0_4_0)) - (((FFTW_REAL) FFTW_K900968867) * (tre0_2_0 + tre0_5_0)) - (((FFTW_REAL) FFTW_K222520933) * (tre0_1_0 + tre0_6_0));
	       tre1_1_0 = (((FFTW_REAL) FFTW_K974927912) * (tim0_6_0 - tim0_1_0)) + (((FFTW_REAL) FFTW_K433883739) * (tim0_2_0 - tim0_5_0)) + (((FFTW_REAL) FFTW_K781831482) * (tim0_3_0 - tim0_4_0));
	       c_re(inout[2 * stride]) = tre1_0_0 + tre1_1_0;
	       c_re(inout[5 * stride]) = tre1_0_0 - tre1_1_0;
	  }
	  {
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tim1_1_0;
	       tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tim0_3_0 + tim0_4_0)) - (((FFTW_REAL) FFTW_K900968867) * (tim0_2_0 + tim0_5_0)) - (((FFTW_REAL) FFTW_K222520933) * (tim0_1_0 + tim0_6_0));
	       tim1_1_0 = (((FFTW_REAL) FFTW_K974927912) * (tre0_1_0 - tre0_6_0)) + (((FFTW_REAL) FFTW_K433883739) * (tre0_5_0 - tre0_2_0)) + (((FFTW_REAL) FFTW_K781831482) * (tre0_4_0 - tre0_3_0));
	       c_im(inout[2 * stride]) = tim1_0_0 + tim1_1_0;
	       c_im(inout[5 * stride]) = tim1_0_0 - tim1_1_0;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tre1_1_0;
	       tre1_0_0 = tre0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tre0_2_0 + tre0_5_0)) - (((FFTW_REAL) FFTW_K222520933) * (tre0_3_0 + tre0_4_0)) - (((FFTW_REAL) FFTW_K900968867) * (tre0_1_0 + tre0_6_0));
	       tre1_1_0 = (((FFTW_REAL) FFTW_K433883739) * (tim0_6_0 - tim0_1_0)) + (((FFTW_REAL) FFTW_K781831482) * (tim0_2_0 - tim0_5_0)) + (((FFTW_REAL) FFTW_K974927912) * (tim0_4_0 - tim0_3_0));
	       c_re(inout[3 * stride]) = tre1_0_0 + tre1_1_0;
	       c_re(inout[4 * stride]) = tre1_0_0 - tre1_1_0;
	  }
	  {
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tim1_1_0;
	       tim1_0_0 = tim0_0_0 + (((FFTW_REAL) FFTW_K623489801) * (tim0_2_0 + tim0_5_0)) - (((FFTW_REAL) FFTW_K222520933) * (tim0_3_0 + tim0_4_0)) - (((FFTW_REAL) FFTW_K900968867) * (tim0_1_0 + tim0_6_0));
	       tim1_1_0 = (((FFTW_REAL) FFTW_K433883739) * (tre0_1_0 - tre0_6_0)) + (((FFTW_REAL) FFTW_K781831482) * (tre0_5_0 - tre0_2_0)) + (((FFTW_REAL) FFTW_K974927912) * (tre0_3_0 - tre0_4_0));
	       c_im(inout[3 * stride]) = tim1_0_0 + tim1_1_0;
	       c_im(inout[4 * stride]) = tim1_0_0 - tim1_1_0;
	  }
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 66 FP additions and 32 FP multiplications */

void fftwi_twiddle_8(FFTW_COMPLEX *A, const FFTW_COMPLEX *W, int stride, int m, int dist)
{
     int i;
     COMPLEX *inout;
     inout = A;
     for (i = 0; i < m; i = i + 1, inout = inout + dist, W = W + 7) {
	  FFTW_REAL tre0_0_0;
	  FFTW_REAL tim0_0_0;
	  FFTW_REAL tre0_0_1;
	  FFTW_REAL tim0_0_1;
	  FFTW_REAL tre0_0_2;
	  FFTW_REAL tim0_0_2;
	  FFTW_REAL tre0_0_3;
	  FFTW_REAL tim0_0_3;
	  FFTW_REAL tre0_1_0;
	  FFTW_REAL tim0_1_0;
	  FFTW_REAL tre0_1_1;
	  FFTW_REAL tim0_1_1;
	  FFTW_REAL tre0_1_2;
	  FFTW_REAL tim0_1_2;
	  FFTW_REAL tre0_1_3;
	  FFTW_REAL tim0_1_3;
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       tre1_0_0 = c_re(inout[0]);
	       tim1_0_0 = c_im(inout[0]);
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[4 * stride]);
		    ti = c_im(inout[4 * stride]);
		    twr = c_re(W[3]);
		    twi = c_im(W[3]);
		    tre1_1_0 = (tr * twr) + (ti * twi);
		    tim1_1_0 = (ti * twr) - (tr * twi);
	       }
	       tre0_0_0 = tre1_0_0 + tre1_1_0;
	       tim0_0_0 = tim1_0_0 + tim1_1_0;
	       tre0_1_0 = tre1_0_0 - tre1_1_0;
	       tim0_1_0 = tim1_0_0 - tim1_1_0;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[stride]);
		    ti = c_im(inout[stride]);
		    twr = c_re(W[0]);
		    twi = c_im(W[0]);
		    tre1_0_0 = (tr * twr) + (ti * twi);
		    tim1_0_0 = (ti * twr) - (tr * twi);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[5 * stride]);
		    ti = c_im(inout[5 * stride]);
		    twr = c_re(W[4]);
		    twi = c_im(W[4]);
		    tre1_1_0 = (tr * twr) + (ti * twi);
		    tim1_1_0 = (ti * twr) - (tr * twi);
	       }
	       tre0_0_1 = tre1_0_0 + tre1_1_0;
	       tim0_0_1 = tim1_0_0 + tim1_1_0;
	       tre0_1_1 = tre1_0_0 - tre1_1_0;
	       tim0_1_1 = tim1_0_0 - tim1_1_0;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[2 * stride]);
		    ti = c_im(inout[2 * stride]);
		    twr = c_re(W[1]);
		    twi = c_im(W[1]);
		    tre1_0_0 = (tr * twr) + (ti * twi);
		    tim1_0_0 = (ti * twr) - (tr * twi);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[6 * stride]);
		    ti = c_im(inout[6 * stride]);
		    twr = c_re(W[5]);
		    twi = c_im(W[5]);
		    tre1_1_0 = (tr * twr) + (ti * twi);
		    tim1_1_0 = (ti * twr) - (tr * twi);
	       }
	       tre0_0_2 = tre1_0_0 + tre1_1_0;
	       tim0_0_2 = tim1_0_0 + tim1_1_0;
	       tre0_1_2 = tre1_0_0 - tre1_1_0;
	       tim0_1_2 = tim1_0_0 - tim1_1_0;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[3 * stride]);
		    ti = c_im(inout[3 * stride]);
		    twr = c_re(W[2]);
		    twi = c_im(W[2]);
		    tre1_0_0 = (tr * twr) + (ti * twi);
		    tim1_0_0 = (ti * twr) - (tr * twi);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[7 * stride]);
		    ti = c_im(inout[7 * stride]);
		    twr = c_re(W[6]);
		    twi = c_im(W[6]);
		    tre1_1_0 = (tr * twr) + (ti * twi);
		    tim1_1_0 = (ti * twr) - (tr * twi);
	       }
	       tre0_0_3 = tre1_0_0 + tre1_1_0;
	       tim0_0_3 = tim1_0_0 + tim1_1_0;
	       tre0_1_3 = tre1_0_0 - tre1_1_0;
	       tim0_1_3 = tim1_0_0 - tim1_1_0;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       tre1_0_0 = tre0_0_0 + tre0_0_2;
	       tim1_0_0 = tim0_0_0 + tim0_0_2;
	       tre1_1_0 = tre0_0_0 - tre0_0_2;
	       tim1_1_0 = tim0_0_0 - tim0_0_2;
	       tre1_0_1 = tre0_0_1 + tre0_0_3;
	       tim1_0_1 = tim0_0_1 + tim0_0_3;
	       tre1_1_1 = tre0_0_1 - tre0_0_3;
	       tim1_1_1 = tim0_0_1 - tim0_0_3;
	       c_re(inout[0]) = tre1_0_0 + tre1_0_1;
	       c_im(inout[0]) = tim1_0_0 + tim1_0_1;
	       c_re(inout[4 * stride]) = tre1_0_0 - tre1_0_1;
	       c_im(inout[4 * stride]) = tim1_0_0 - tim1_0_1;
	       c_re(inout[2 * stride]) = tre1_1_0 - tim1_1_1;
	       c_im(inout[2 * stride]) = tim1_1_0 + tre1_1_1;
	       c_re(inout[6 * stride]) = tre1_1_0 + tim1_1_1;
	       c_im(inout[6 * stride]) = tim1_1_0 - tre1_1_1;
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_0_1;
	       FFTW_REAL tim1_0_1;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_1_1;
	       FFTW_REAL tim1_1_1;
	       tre1_0_0 = tre0_1_0 - tim0_1_2;
	       tim1_0_0 = tim0_1_0 + tre0_1_2;
	       tre1_1_0 = tre0_1_0 + tim0_1_2;
	       tim1_1_0 = tim0_1_0 - tre0_1_2;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tre2_1_0;
		    FFTW_REAL tim2_1_0;
		    tre2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_1_1 - tim0_1_1);
		    tim2_0_0 = ((FFTW_REAL) FFTW_K707106781) * (tim0_1_1 + tre0_1_1);
		    tre2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_1_3 + tim0_1_3);
		    tim2_1_0 = ((FFTW_REAL) FFTW_K707106781) * (tre0_1_3 - tim0_1_3);
		    tre1_0_1 = tre2_0_0 - tre2_1_0;
		    tim1_0_1 = tim2_0_0 + tim2_1_0;
		    tre1_1_1 = tre2_0_0 + tre2_1_0;
		    tim1_1_1 = tim2_0_0 - tim2_1_0;
	       }
	       c_re(inout[stride]) = tre1_0_0 + tre1_0_1;
	       c_im(inout[stride]) = tim1_0_0 + tim1_0_1;
	       c_re(inout[5 * stride]) = tre1_0_0 - tre1_0_1;
	       c_im(inout[5 * stride]) = tim1_0_0 - tim1_0_1;
	       c_re(inout[3 * stride]) = tre1_1_0 - tim1_1_1;
	       c_im(inout[3 * stride]) = tim1_1_0 + tre1_1_1;
	       c_re(inout[7 * stride]) = tre1_1_0 + tim1_1_1;
	       c_im(inout[7 * stride]) = tim1_1_0 - tre1_1_1;
	  }
     }
}
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

/* This file has been automatically generated --- DO NOT EDIT */

#include "fftw.h"
#include "konst.h"

/* Generated by $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

/* This function contains 108 FP additions and 72 FP multiplications */

void fftwi_twiddle_9(FFTW_COMPLEX *A, const FFTW_COMPLEX *W, int stride, int m, int dist)
{
     int i;
     COMPLEX *inout;
     inout = A;
     for (i = 0; i < m; i = i + 1, inout = inout + dist, W = W + 8) {
	  FFTW_REAL tre0_0_0;
	  FFTW_REAL tim0_0_0;
	  FFTW_REAL tre0_0_1;
	  FFTW_REAL tim0_0_1;
	  FFTW_REAL tre0_0_2;
	  FFTW_REAL tim0_0_2;
	  FFTW_REAL tre0_1_0;
	  FFTW_REAL tim0_1_0;
	  FFTW_REAL tre0_1_1;
	  FFTW_REAL tim0_1_1;
	  FFTW_REAL tre0_1_2;
	  FFTW_REAL tim0_1_2;
	  FFTW_REAL tre0_2_0;
	  FFTW_REAL tim0_2_0;
	  FFTW_REAL tre0_2_1;
	  FFTW_REAL tim0_2_1;
	  FFTW_REAL tre0_2_2;
	  FFTW_REAL tim0_2_2;
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_2_0;
	       FFTW_REAL tim1_2_0;
	       tre1_0_0 = c_re(inout[0]);
	       tim1_0_0 = c_im(inout[0]);
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[3 * stride]);
		    ti = c_im(inout[3 * stride]);
		    twr = c_re(W[2]);
		    twi = c_im(W[2]);
		    tre1_1_0 = (tr * twr) + (ti * twi);
		    tim1_1_0 = (ti * twr) - (tr * twi);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[6 * stride]);
		    ti = c_im(inout[6 * stride]);
		    twr = c_re(W[5]);
		    twi = c_im(W[5]);
		    tre1_2_0 = (tr * twr) + (ti * twi);
		    tim1_2_0 = (ti * twr) - (tr * twi);
	       }
	       tre0_0_0 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	       tim0_0_0 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tre2_1_0;
		    tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
		    tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_2_0 - tim1_1_0);
		    tre0_1_0 = tre2_0_0 + tre2_1_0;
		    tre0_2_0 = tre2_0_0 - tre2_1_0;
	       }
	       {
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tim2_1_0;
		    tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
		    tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_1_0 - tre1_2_0);
		    tim0_1_0 = tim2_0_0 + tim2_1_0;
		    tim0_2_0 = tim2_0_0 - tim2_1_0;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_2_0;
	       FFTW_REAL tim1_2_0;
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[stride]);
		    ti = c_im(inout[stride]);
		    twr = c_re(W[0]);
		    twi = c_im(W[0]);
		    tre1_0_0 = (tr * twr) + (ti * twi);
		    tim1_0_0 = (ti * twr) - (tr * twi);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[4 * stride]);
		    ti = c_im(inout[4 * stride]);
		    twr = c_re(W[3]);
		    twi = c_im(W[3]);
		    tre1_1_0 = (tr * twr) + (ti * twi);
		    tim1_1_0 = (ti * twr) - (tr * twi);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[7 * stride]);
		    ti = c_im(inout[7 * stride]);
		    twr = c_re(W[6]);
		    twi = c_im(W[6]);
		    tre1_2_0 = (tr * twr) + (ti * twi);
		    tim1_2_0 = (ti * twr) - (tr * twi);
	       }
	       tre0_0_1 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	       tim0_0_1 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tre2_1_0;
		    tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
		    tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_2_0 - tim1_1_0);
		    tre0_1_1 = tre2_0_0 + tre2_1_0;
		    tre0_2_1 = tre2_0_0 - tre2_1_0;
	       }
	       {
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tim2_1_0;
		    tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
		    tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_1_0 - tre1_2_0);
		    tim0_1_1 = tim2_0_0 + tim2_1_0;
		    tim0_2_1 = tim2_0_0 - tim2_1_0;
	       }
	  }
	  {
	       FFTW_REAL tre1_0_0;
	       FFTW_REAL tim1_0_0;
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_2_0;
	       FFTW_REAL tim1_2_0;
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[2 * stride]);
		    ti = c_im(inout[2 * stride]);
		    twr = c_re(W[1]);
		    twi = c_im(W[1]);
		    tre1_0_0 = (tr * twr) + (ti * twi);
		    tim1_0_0 = (ti * twr) - (tr * twi);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[5 * stride]);
		    ti = c_im(inout[5 * stride]);
		    twr = c_re(W[4]);
		    twi = c_im(W[4]);
		    tre1_1_0 = (tr * twr) + (ti * twi);
		    tim1_1_0 = (ti * twr) - (tr * twi);
	       }
	       {
		    FFTW_REAL tr;
		    FFTW_REAL ti;
		    FFTW_REAL twr;
		    FFTW_REAL twi;
		    tr = c_re(inout[8 * stride]);
		    ti = c_im(inout[8 * stride]);
		    twr = c_re(W[7]);
		    twi = c_im(W[7]);
		    tre1_2_0 = (tr * twr) + (ti * twi);
		    tim1_2_0 = (ti * twr) - (tr * twi);
	       }
	       tre0_0_2 = tre1_0_0 + tre1_1_0 + tre1_2_0;
	       tim0_0_2 = tim1_0_0 + tim1_1_0 + tim1_2_0;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tre2_1_0;
		    tre2_0_0 = tre1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
		    tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_2_0 - tim1_1_0);
		    tre0_1_2 = tre2_0_0 + tre2_1_0;
		    tre0_2_2 = tre2_0_0 - tre2_1_0;
	       }
	       {
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tim2_1_0;
		    tim2_0_0 = tim1_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
		    tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_1_0 - tre1_2_0);
		    tim0_1_2 = tim2_0_0 + tim2_1_0;
		    tim0_2_2 = tim2_0_0 - tim2_1_0;
	       }
	  }
	  c_re(inout[0]) = tre0_0_0 + tre0_0_1 + tre0_0_2;
	  c_im(inout[0]) = tim0_0_0 + tim0_0_1 + tim0_0_2;
	  {
	       FFTW_REAL tre2_0_0;
	       FFTW_REAL tre2_1_0;
	       tre2_0_0 = tre0_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tre0_0_1 + tre0_0_2));
	       tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim0_0_2 - tim0_0_1);
	       c_re(inout[3 * stride]) = tre2_0_0 + tre2_1_0;
	       c_re(inout[6 * stride]) = tre2_0_0 - tre2_1_0;
	  }
	  {
	       FFTW_REAL tim2_0_0;
	       FFTW_REAL tim2_1_0;
	       tim2_0_0 = tim0_0_0 - (((FFTW_REAL) FFTW_K499999999) * (tim0_0_1 + tim0_0_2));
	       tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre0_0_1 - tre0_0_2);
	       c_im(inout[3 * stride]) = tim2_0_0 + tim2_1_0;
	       c_im(inout[6 * stride]) = tim2_0_0 - tim2_1_0;
	  }
	  {
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_2_0;
	       FFTW_REAL tim1_2_0;
	       tre1_1_0 = (((FFTW_REAL) FFTW_K766044443) * tre0_1_1) - (((FFTW_REAL) FFTW_K642787609) * tim0_1_1);
	       tim1_1_0 = (((FFTW_REAL) FFTW_K766044443) * tim0_1_1) + (((FFTW_REAL) FFTW_K642787609) * tre0_1_1);
	       tre1_2_0 = (((FFTW_REAL) FFTW_K173648177) * tre0_1_2) - (((FFTW_REAL) FFTW_K984807753) * tim0_1_2);
	       tim1_2_0 = (((FFTW_REAL) FFTW_K173648177) * tim0_1_2) + (((FFTW_REAL) FFTW_K984807753) * tre0_1_2);
	       c_re(inout[stride]) = tre0_1_0 + tre1_1_0 + tre1_2_0;
	       c_im(inout[stride]) = tim0_1_0 + tim1_1_0 + tim1_2_0;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tre2_1_0;
		    tre2_0_0 = tre0_1_0 - (((FFTW_REAL) FFTW_K499999999) * (tre1_1_0 + tre1_2_0));
		    tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_2_0 - tim1_1_0);
		    c_re(inout[4 * stride]) = tre2_0_0 + tre2_1_0;
		    c_re(inout[7 * stride]) = tre2_0_0 - tre2_1_0;
	       }
	       {
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tim2_1_0;
		    tim2_0_0 = tim0_1_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
		    tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_1_0 - tre1_2_0);
		    c_im(inout[4 * stride]) = tim2_0_0 + tim2_1_0;
		    c_im(inout[7 * stride]) = tim2_0_0 - tim2_1_0;
	       }
	  }
	  {
	       FFTW_REAL tre1_1_0;
	       FFTW_REAL tim1_1_0;
	       FFTW_REAL tre1_2_0;
	       FFTW_REAL tim1_2_0;
	       tre1_1_0 = (((FFTW_REAL) FFTW_K173648177) * tre0_2_1) - (((FFTW_REAL) FFTW_K984807753) * tim0_2_1);
	       tim1_1_0 = (((FFTW_REAL) FFTW_K173648177) * tim0_2_1) + (((FFTW_REAL) FFTW_K984807753) * tre0_2_1);
	       tre1_2_0 = (((FFTW_REAL) FFTW_K939692620) * tre0_2_2) + (((FFTW_REAL) FFTW_K342020143) * tim0_2_2);
	       tim1_2_0 = (((FFTW_REAL) FFTW_K342020143) * tre0_2_2) - (((FFTW_REAL) FFTW_K939692620) * tim0_2_2);
	       c_re(inout[2 * stride]) = tre0_2_0 + tre1_1_0 - tre1_2_0;
	       c_im(inout[2 * stride]) = tim0_2_0 + tim1_1_0 + tim1_2_0;
	       {
		    FFTW_REAL tre2_0_0;
		    FFTW_REAL tre2_1_0;
		    tre2_0_0 = tre0_2_0 + (((FFTW_REAL) FFTW_K499999999) * (tre1_2_0 - tre1_1_0));
		    tre2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tim1_2_0 - tim1_1_0);
		    c_re(inout[5 * stride]) = tre2_0_0 + tre2_1_0;
		    c_re(inout[8 * stride]) = tre2_0_0 - tre2_1_0;
	       }
	       {
		    FFTW_REAL tim2_0_0;
		    FFTW_REAL tim2_1_0;
		    tim2_0_0 = tim0_2_0 - (((FFTW_REAL) FFTW_K499999999) * (tim1_1_0 + tim1_2_0));
		    tim2_1_0 = ((FFTW_REAL) FFTW_K866025403) * (tre1_1_0 + tre1_2_0);
		    c_im(inout[5 * stride]) = tim2_0_0 + tim2_1_0;
		    c_im(inout[8 * stride]) = tim2_0_0 - tim2_1_0;
	       }
	  }
     }
}
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

/*
 * generic.c -- "generic" solvers.  They work for all
 * n (and are slow)
 */
#include "fftw.h"
#include <math.h>
#include <stdlib.h>

void fftw_twiddle_generic(FFTW_COMPLEX *A, const FFTW_COMPLEX *W,
			  int m, int r, int n, int stride)
{
     int i, j, k;
     const FFTW_COMPLEX *jp;
     FFTW_COMPLEX *kp;
     FFTW_COMPLEX *tmp = (FFTW_COMPLEX *)
     fftw_malloc(r * sizeof(FFTW_COMPLEX));

     for (i = 0; i < m; ++i) {
	  for (k = 0, kp = tmp; k < r; ++k, kp++) {
	       FFTW_REAL r0, i0, rt, it, rw, iw;
	       int l1 = i + m * k;
	       int l0;

	       r0 = i0 = 0.0;
	       for (j = 0, jp = A + i * stride, l0 = 0; j < r; ++j,
		    jp += m * stride) {
		    rw = c_re(W[l0]);
		    iw = c_im(W[l0]);
		    rt = c_re(*jp);
		    it = c_im(*jp);
		    r0 += rt * rw - it * iw;
		    i0 += rt * iw + it * rw;
		    l0 += l1;
		    if (l0 > n)
			 l0 -= n;
	       }
	       c_re(*kp) = r0;
	       c_im(*kp) = i0;
	  }
	  for (k = 0, kp = A + i * stride; k < r; ++k, kp += m * stride)
	       *kp = tmp[k];
     }

     fftw_free(tmp);
}

void fftwi_twiddle_generic(FFTW_COMPLEX *A, const FFTW_COMPLEX *W,
			   int m, int r, int n, int stride)
{
     int i, j, k;
     const FFTW_COMPLEX *jp;
     FFTW_COMPLEX *kp;
     FFTW_COMPLEX *tmp = (FFTW_COMPLEX *)
     fftw_malloc(r * sizeof(FFTW_COMPLEX));

     for (i = 0; i < m; ++i) {
	  for (k = 0, kp = tmp; k < r; ++k, kp++) {
	       FFTW_REAL r0, i0, rt, it, rw, iw;
	       int l1 = i + m * k;
	       int l0;

	       r0 = i0 = 0.0;
	       for (j = 0, jp = A + i * stride, l0 = 0; j < r; ++j,
		    jp += m * stride) {
		    rw = c_re(W[l0]);
		    iw = c_im(W[l0]);
		    rt = c_re(*jp);
		    it = c_im(*jp);
		    r0 += rt * rw + it * iw;
		    i0 += it * rw - rt * iw;
		    l0 += l1;
		    if (l0 > n)
			 l0 -= n;
	       }
	       c_re(*kp) = r0;
	       c_im(*kp) = i0;
	  }
	  for (k = 0, kp = A + i * stride; k < r; ++k, kp += m * stride)
	       *kp = tmp[k];
     }

     fftw_free(tmp);
}
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

/*
 * malloc.c -- memory allocation related functions
 */

/* $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */
#if defined FFTW_USING_CILK
#include <cilk.h>
#include <cilk-compat.h>
#endif

#include "fftw.h"
#include <stdio.h>
#include <stdlib.h>

int fftw_malloc_cnt = 0;
void *(*fftw_malloc_hook) (size_t n) = (void *(*)(size_t n)) 0;
void (*fftw_free_hook) (void *p) = (void (*)(void *p)) 0;

#define FFTW_MALLOC_DEBUG 0
/* sorry for this debugging hack ... */
#define COMMA ,

#if FFTW_MALLOC_DEBUG
#define WHEN_DEBUG(a) a

/*
 * debugging malloc/free.  Initialize every malloced and freed area to
 * random values, just to make sure we are not using uninitialized
 * pointers.  Also check for writes past the ends of allocated blocks,
 * and a couple of other things.
 *
 * This code is a quick and dirty hack -- use at your own risk.
 */

int fftw_malloc_total = 0;

#define MAGIC 0xABadCafe
#define PAD_FACTOR 2
#define TWOINTS (2 * sizeof(int))

#define VERBOSE_ALLOCATION 0

#if VERBOSE_ALLOCATION
#define WHEN_VERBOSE(a) a
#else
#define WHEN_VERBOSE(a) 
#endif

void *fftw_malloc(size_t n)
{
     char *p;
     int i;

     WHEN_VERBOSE({
	  printf("FFTW_MALLOC %d\n",n);
	  fflush(stdout);
     })

     if (n == 0)
	  fftw_die("Tried to allocate a block of zero size!\n");

     fftw_malloc_total += n;

     p = (char *) malloc(PAD_FACTOR*n + TWOINTS);
     if (!p)
	  fftw_die("fftw_malloc: out of memory\n");

     /* store the size in a known position */
     ((int *) p)[0] = n;
     ((int *) p)[1] = MAGIC;
     for (i = 0; i < PAD_FACTOR*n; ++i)
	  p[i + TWOINTS] = (char) (i ^ 0xDEADBEEF);

     ++fftw_malloc_cnt;

     /* skip the size we stored previously */
     return (void *) (p + TWOINTS);
}

void fftw_free(void *p)
{
     char *q = ((char *) p) - TWOINTS;

     if (!p)
	  fftw_die("fftw_free: tried to free NULL pointer!\n");

     if (!q)
	  fftw_die("fftw_free: tried to free NULL+TWOINTS pointer!\n");

     {
	  int n = ((int *) q)[0];
	  int magic = ((int *) q)[1];
	  int i;
	  
	  WHEN_VERBOSE({
	       printf("FFTW_FREE %d\n",n);
	       fflush(stdout);
	  })
	  
	  if (n == 0)
	       fftw_die("Tried to free a freed pointer!\n");
	  *((int *) q) = 0; /* set to zero to detect duplicate free's */
	  
	  if (magic != MAGIC)
	       fftw_die("Wrong magic in fftw_free()!\n");	       
	  ((int *) q)[1] = ~MAGIC;

	  if (n < 0)
	       fftw_die("Tried to free block with corrupt size descriptor!\n");
	  
	  fftw_malloc_total -= n;
	  
	  if (fftw_malloc_total < 0)
	       fftw_die("fftw_malloc_total went negative!\n");
	  
	  /* check for writing past end of array: */
	  for (i = n; i < PAD_FACTOR*n; ++i)
	       if (q[i+TWOINTS] != (char) (i ^ 0xDEADBEEF)) {
		    fprintf(stderr, "Byte %d past end of array has changed!\n",
			    i - n + 1);
		    fftw_die("Array bounds overwritten!\n");
	       }
	  
	  for (i = 0; i < PAD_FACTOR*n; ++i)
	       q[i + TWOINTS] = (char) (i ^ 0xBEEFDEAD);

	  --fftw_malloc_cnt;
	  free(q);
     }
}

#else				/* production version, no hacks */
#define WHEN_DEBUG(a) 

void *fftw_malloc(size_t n)
{
     void *p;

     if (fftw_malloc_hook)
	  return fftw_malloc_hook(n);

     if (n == 0)
	  n = 1;

     p = malloc(n);

     if (!p)
	  fftw_die("fftw_malloc: out of memory\n");

     return p;
}

void fftw_free(void *p)
{
     if (p) {
	  if (fftw_free_hook) {
	       fftw_free_hook(p);
	       return;
	  }
	  free(p);
     }
}

#endif

/* die when fatal errors occur */
void fftw_die(char *s)
{
     fprintf(stderr, "%s", s);
     exit(1);
}

/* check for memory leaks when debugging */
void fftw_check_memory_leaks(void)
{
     extern int fftw_node_cnt, fftw_plan_cnt, fftw_twiddle_size;

     if (WHEN_DEBUG(fftw_malloc_cnt ||)
	 WHEN_DEBUG(fftw_malloc_total ||)
	 fftw_node_cnt || fftw_plan_cnt || fftw_twiddle_size) {
	  fprintf(stderr,
		  "MEMORY LEAK!!!\n"
		  WHEN_DEBUG("fftw_malloc = %d")
		  " node=%d plan=%d twiddle=%d\n"
		  WHEN_DEBUG("fftw_malloc_total = %d\n"), 
		  WHEN_DEBUG(fftw_malloc_cnt COMMA)
		  fftw_node_cnt, fftw_plan_cnt, fftw_twiddle_size
		  WHEN_DEBUG(COMMA fftw_malloc_total));
	  exit(1);
     }
}

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

/* $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */
#include "fftw.h"
#include <math.h>

/*
 * Naive O(n^2) algorithm, used for testing purposes
 */
void fftw_naive(int n, FFTW_COMPLEX *in, FFTW_COMPLEX *out)
{
     int i, j;
     FFTW_COMPLEX sum;
     FFTW_COMPLEX w;
     FFTW_REAL pi = 3.1415926535897932384626434;

     for (j = 0; j < n; ++j) {
	  c_re(sum) = c_im(sum) = 0.0;
	  for (i = 0; i < n; ++i) {
	       c_re(w) = cos((2.0 * pi * (i * j % n)) / n);
	       c_im(w) = -sin((2.0 * pi * (i * j % n)) / n);
	       c_re(sum) += c_re(in[i]) * c_re(w) - c_im(in[i]) * c_im(w);
	       c_im(sum) += c_im(in[i]) * c_re(w) + c_re(in[i]) * c_im(w);
	  }
	  out[j] = sum;
     }
     return;
}

/*
 * Naive O(n^2) algorithm, for the inverse.
 */
void fftwi_naive(int n, FFTW_COMPLEX *in, FFTW_COMPLEX *out)
{
     int i, j;
     FFTW_COMPLEX sum;
     FFTW_COMPLEX w;
     FFTW_REAL pi = 3.1415926535897932384626434;

     for (j = 0; j < n; ++j) {
	  c_re(sum) = c_im(sum) = 0.0;
	  for (i = 0; i < n; ++i) {
	       c_re(w) = cos((2.0 * pi * (i * j % n)) / n);
	       c_im(w) = sin((2.0 * pi * (i * j % n)) / n);
	       c_re(sum) += c_re(in[i]) * c_re(w) - c_im(in[i]) * c_im(w);
	       c_im(sum) += c_im(in[i]) * c_re(w) + c_re(in[i]) * c_im(w);
	  }
	  out[j] = sum;
     }
     return;
}
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

/*
 * planner.c -- find the optimal plan
 */

/* $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */
#if defined FFTW_USING_CILK
#include <cilk.h>
#include <cilk-compat.h>
#endif

#include "fftw.h"
#include <stdlib.h>
#include <stdio.h>

int fftw_node_cnt = 0;
int fftw_plan_cnt = 0;

#define NOTW_OPTIMAL_SIZE 32
#define TWIDDLE_OPTIMAL_SIZE 12

/* wisdom prototypes */
extern int fftw_wisdom_lookup(int n, int flags, fftw_direction dir,
			    enum fftw_node_type *type,
			    int *signature, int replace_p);
extern void fftw_wisdom_add(int n, int flags, fftw_direction dir,
			  enum fftw_node_type type,
			  int signature);

/* constructors --- I wish I had ML */
static fftw_plan_node *make_node(void)
{
     fftw_plan_node *p = (fftw_plan_node *)
     fftw_malloc(sizeof(fftw_plan_node));
     p->refcnt = 0;
     fftw_node_cnt++;
     return p;
}

static void use_node(fftw_plan_node *p)
{
     ++p->refcnt;
}

static fftw_plan_node *make_node_notw(int size, notw_codelet *codelet)
{
     fftw_plan_node *p = make_node();

     p->type = FFTW_NOTW;
     p->nodeu.notw.size = size;
     p->nodeu.notw.codelet = codelet;
     return p;
}

static fftw_plan_node *make_node_twiddle(int n, int size, twiddle_codelet *codelet,
					 fftw_plan_node *recurse,
					 int flags)
{
     fftw_plan_node *p = make_node();

     p->type = FFTW_TWIDDLE;
     p->nodeu.twiddle.size = size;
     p->nodeu.twiddle.codelet = codelet;
     p->nodeu.twiddle.recurse = recurse;
     use_node(recurse);
     if (flags & FFTW_MEASURE)
	  p->nodeu.twiddle.tw = fftw_create_twiddle(n, size, n / size);
     else
	  p->nodeu.twiddle.tw = 0;
     return p;
}

static fftw_plan_node *make_node_generic(int n, int size,
					 generic_codelet *codelet,
					 fftw_plan_node *recurse,
					 int flags)
{
     fftw_plan_node *p = make_node();

     p->type = FFTW_GENERIC;
     p->nodeu.generic.size = size;
     p->nodeu.generic.codelet = codelet;
     p->nodeu.generic.recurse = recurse;
     use_node(recurse);

     if (flags & FFTW_MEASURE)
	  p->nodeu.generic.tw = fftw_create_twiddle(n, 2, n);
     else
	  p->nodeu.generic.tw = 0;
     return p;
}

static void destroy_tree(fftw_plan_node *p)
{
     if (p) {
	  --p->refcnt;
	  if (p->refcnt == 0) {
	       switch (p->type) {
		   case FFTW_NOTW:
			break;

		   case FFTW_TWIDDLE:
			if (p->nodeu.twiddle.tw)
			     fftw_destroy_twiddle(p->nodeu.twiddle.tw);
			destroy_tree(p->nodeu.twiddle.recurse);
			break;

		   case FFTW_GENERIC:
			if (p->nodeu.generic.tw)
			     fftw_destroy_twiddle(p->nodeu.generic.tw);
			destroy_tree(p->nodeu.generic.recurse);
			break;
	       }

	       fftw_free(p);
	       fftw_node_cnt--;
	  }
     }
}

/* create a plan with twiddle factors, and other bells and whistles */
static fftw_plan make_plan(int n, fftw_direction dir,
			   fftw_plan_node *root, int flags,
			   enum fftw_node_type wisdom_type,
			   int wisdom_signature)
{
     fftw_plan p = (fftw_plan) fftw_malloc(sizeof(struct fftw_plan_struct));

     p->n = n;
     p->dir = dir;
     p->flags = flags;
     use_node(root);
     p->root = root;
     p->cost = 0.0;
     p->wisdom_type = wisdom_type;
     p->wisdom_signature = wisdom_signature;
     p->next = (fftw_plan) 0;
     p->refcnt = 0;
     fftw_plan_cnt++;
     return p;
}

/*
 * complete with twiddle factors (because nodes don't have
 * them when FFTW_ESTIMATE is set)
 */
static void complete_twiddle(fftw_plan_node *p, int n)
{
     int r;
     switch (p->type) {
	 case FFTW_NOTW:
	      break;

	 case FFTW_TWIDDLE:
	      r = p->nodeu.twiddle.size;
	      if (!p->nodeu.twiddle.tw)
		   p->nodeu.twiddle.tw = fftw_create_twiddle(n, r, n / r);
	      complete_twiddle(p->nodeu.twiddle.recurse, n / r);
	      break;

	 case FFTW_GENERIC:
	      r = p->nodeu.generic.size;
	      if (!p->nodeu.generic.tw)
		   p->nodeu.generic.tw = fftw_create_twiddle(n, 2, n);
	      complete_twiddle(p->nodeu.generic.recurse, n / r);
	      break;
     }
}

static void use_plan(fftw_plan p)
{
     ++p->refcnt;
}

static void destroy_plan(fftw_plan p)
{
     --p->refcnt;

     if (p->refcnt == 0) {
	  destroy_tree(p->root);
	  fftw_plan_cnt--;
	  fftw_free(p);
     }
}

/* end of constructors */

/* management of plan tables */
static void make_empty_table(fftw_plan *table)
{
     *table = (fftw_plan) 0;
}

static void insert(fftw_plan *table, fftw_plan this_plan, int n)
{
     use_plan(this_plan);
     this_plan->n = n;
     this_plan->next = *table;
     *table = this_plan;
}

static fftw_plan lookup(fftw_plan *table, int n, int flags)
{
     fftw_plan p;

     for (p = *table; p &&
	  ((p->n != n) || (p->flags != flags)); p = p->next);

     return p;
}

static void destroy_table(fftw_plan *table)
{
     fftw_plan p, q;

     for (p = *table; p; p = q) {
	  q = p->next;
	  destroy_plan(p);
     }
}

static double estimate_node(fftw_plan_node *p)
{
     int k;

     switch (p->type) {
	 case FFTW_NOTW:
	      k = p->nodeu.notw.size;
	      return 1.0 + 0.1 * (k - NOTW_OPTIMAL_SIZE) *
		  (k - NOTW_OPTIMAL_SIZE);

	 case FFTW_TWIDDLE:
	      k = p->nodeu.twiddle.size;
	      return 1.0 + 0.1 * (k - TWIDDLE_OPTIMAL_SIZE) *
		  (k - TWIDDLE_OPTIMAL_SIZE)
		  + estimate_node(p->nodeu.twiddle.recurse);

	 case FFTW_GENERIC:
	      k = p->nodeu.generic.size;
	      return 10.0 + k * k
		  + estimate_node(p->nodeu.generic.recurse);
     }
     return 1.0E20;
}

/* auxiliary functions */
static void compute_cost(fftw_plan plan)
{
     if (plan->flags & FFTW_MEASURE)
	  plan->cost = fftw_measure_runtime(plan);
     else {
	  double c;
	  c = plan->n * estimate_node(plan->root);
	  plan->cost = c;
     }
}

/* pick the better of two plans and destroy the other one. */
static fftw_plan pick_better(fftw_plan p1, fftw_plan p2)
{
     if (!p1)
	  return p2;

     if (!p2)
	  return p1;

     if (p1->cost > p2->cost) {
	  destroy_plan(p1);
	  return p2;
     } else {
	  destroy_plan(p2);
	  return p1;
     }
}

/* find the smallest prime factor of n */
static int factor(int n)
{
     int r;

     /* try 2 */
     if ((n & 1) == 0)
	  return 2;

     /* try odd numbers up to sqrt(n) */
     for (r = 3; r * r <= n; r += 2)
	  if (n % r == 0)
	       return r;

     /* n is prime */
     return n;
}

/* 
 * Some macrology for the planner.  If you have to write
 * the same line of code twice, there must be some bug.
 */
#define NOTW_ITERATOR(p, dir)                                \
      config_notw *p =                                       \
	  p = (dir == FFTW_FORWARD ?                         \
	       fftw_config_notw : fftwi_config_notw)

#define TWIDDLE_ITERATOR(p, dir)                             \
      config_twiddle *p =                                    \
	  p = (dir == FFTW_FORWARD ?                         \
	       fftw_config_twiddle : fftwi_config_twiddle);

#define FORALL_NOTW(p)             \
	 for (; p->size; ++p) 

#define FORALL_TWIDDLE(p)          \
	 for (; p->size; ++p) 

/******************************************
 *      Recursive planner                 *
 ******************************************/
fftw_plan planner(fftw_plan *table, int n, fftw_direction dir, int flags);

/*
 * the planner consists of two parts: one that tries to
 * use accumulated wisdom, and one that does not.
 * A small driver invokes both parts in sequence
 */

/* planner with wisdom: look up the codelet suggested by the wisdom */
fftw_plan planner_wisdom(fftw_plan *table, int n,
			 fftw_direction dir, int flags)
{
     fftw_plan best = (fftw_plan) 0;
     fftw_plan_node *node;
     int have_wisdom;
     enum fftw_node_type wisdom_type;
     int wisdom_signature;

     /* see if we remember any wisdom for this case */
     have_wisdom = fftw_wisdom_lookup(n, flags, dir, 
				      &wisdom_type, &wisdom_signature, 0);

     if (!have_wisdom)
	  return best;

     if (wisdom_type == FFTW_NOTW) {
	  NOTW_ITERATOR(p, dir);
	       
	  FORALL_NOTW(p) {
	       /* see if wisdom applies */
	       if (wisdom_signature == p->signature &&
		   p->size == n) {
		    node = make_node_notw(n, p->codelet);
		    best = make_plan(n, dir, node, flags,
				     FFTW_NOTW, p->signature);
		    use_plan(best);
		    return best;
	       }
	  }
     }
	  
     if (wisdom_type == FFTW_TWIDDLE) {
	  TWIDDLE_ITERATOR(p, dir);
	       
	  FORALL_TWIDDLE(p) {
	       /* see if wisdom applies */
	       if (wisdom_signature == p->signature &&
		   (n % p->size) == 0) {
		    fftw_plan r = planner(table, n / p->size, dir, flags);
		    node = make_node_twiddle(n, p->size, p->codelet,
					     r->root, flags);
		    best = make_plan(n, dir, node, flags,
				     FFTW_TWIDDLE, p->signature);
		    use_plan(best);
		    destroy_plan(r);
		    return best;
	       }
	  }
     }

     /* 
      * BUG (or: TODO)  Can we have generic wisdom? This is probably
      * an academic question
      */

     return best;
}

/*
 * planner with no wisdom: try all combinations and pick
 * the best
 */
fftw_plan planner_normal(fftw_plan *table, int n, fftw_direction dir,
			   int flags)
{
     fftw_plan best = (fftw_plan) 0;
     fftw_plan newplan;
     fftw_plan_node *node;

     /* see if we have any codelet that solves the problem */
     {
	  NOTW_ITERATOR(p, dir);
	       
	  FORALL_NOTW(p) {
	       if (p->size == n) {
		    node = make_node_notw(n, p->codelet);
		    newplan = make_plan(n, dir, node, flags,
					FFTW_NOTW, p->signature);
		    use_plan(newplan);
		    compute_cost(newplan);
		    best = pick_better(newplan, best);
	       }
	  }
     }

     /* Then, try all available twiddle codelets */
     {
	  TWIDDLE_ITERATOR(p, dir);
	       
	  FORALL_TWIDDLE(p) {
	       if ((n % p->size) == 0 &&
		   (!best || n != p->size)) {
		    fftw_plan r = planner(table, n / p->size, dir, flags);
		    node = make_node_twiddle(n, p->size, p->codelet,
					     r->root, flags);
		    newplan = make_plan(n, dir, node, flags,
					FFTW_TWIDDLE, p->signature);
		    use_plan(newplan);
		    destroy_plan(r);
		    compute_cost(newplan);
		    best = pick_better(newplan, best);
	       }
	  }
     }

     /* 
      * if no plan has been found so far, resort to generic codelets 
      */
     if (!best) {
	  generic_codelet *codelet = (dir == FFTW_FORWARD ?
			   fftw_twiddle_generic : fftwi_twiddle_generic);
	  int size = factor(n);
	  fftw_plan r = planner(table, n / size, dir, flags);

	  node = make_node_generic(n, size, codelet, r->root, flags);
	  newplan = make_plan(n, dir, node, flags, FFTW_GENERIC, 0);
	  use_plan(newplan);
	  destroy_plan(r);
	  compute_cost(newplan);
	  best = pick_better(newplan, best);
     }

     return best;
}

fftw_plan planner(fftw_plan *table, int n, fftw_direction dir,
			   int flags)
{
     fftw_plan best = (fftw_plan) 0;

     /* see if plan has already been computed */
     best = lookup(table, n, flags);
     if (best) {
	  use_plan(best);
	  return best;
     }

     /* try a wise plan */
     best = planner_wisdom(table, n, dir, flags);

     if (!best) {
	  /* No wisdom.  Plan normally. */
	  best = planner_normal(table, n, dir, flags);
     }

     if (best) {
	  insert(table, best, n);

	  /* remember the wisdom */
	  fftw_wisdom_add(n, flags, dir, best->wisdom_type,
			  best->wisdom_signature);
     }

     return best;
}

fftw_plan fftw_create_plan(int n, fftw_direction dir, int flags)
{
     fftw_plan table;
     fftw_plan p1;

     /* validate parameters */
     if (n <= 0)
	  return (fftw_plan) 0;

     if ((dir != FFTW_FORWARD) && (dir != FFTW_BACKWARD))
	  return (fftw_plan) 0;

     make_empty_table(&table);
     p1 = planner(&table, n, dir, flags);
     destroy_table(&table);

     complete_twiddle(p1->root, n);
     return p1;
}

void fftw_destroy_plan(fftw_plan plan)
{
     destroy_plan(plan);
}

static void print_node(FILE * f, fftw_plan_node *p, int indent)
{
     if (p) {
	  switch (p->type) {
	      case FFTW_NOTW:
		   fprintf(f, "%*sFFTW_NOTW %d\n", indent, "",
			   p->nodeu.notw.size);
		   break;
	      case FFTW_TWIDDLE:
		   fprintf(f, "%*sFFTW_TWIDDLE %d\n", indent, "",
			   p->nodeu.twiddle.size);
		   print_node(f, p->nodeu.twiddle.recurse, indent);
		   break;
	      case FFTW_GENERIC:
		   fprintf(f, "%*sFFTW_GENERIC %d\n", indent, "",
			   p->nodeu.generic.size);
		   print_node(f, p->nodeu.generic.recurse, indent);
		   break;
	  }
     }
}

void fftw_fprint_plan(FILE * f, fftw_plan p)
{
     fprintf(f, "plan: (cost = %e)\n", p->cost);
     print_node(f, p->root, 0);
}

void fftw_print_plan(fftw_plan p)
{
     fftw_fprint_plan(stdout, p);
}
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

/*
 * timer.c -- this file measures the execution time of 
 *            ffts.  This information is used by the planner.
 */

/* $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */

#include <time.h>
#include "fftw.h"
#include <math.h>
#include <stdlib.h>

/*
 * The timer keeps doubling the number of iterations
 * until the program runs for more than FFTW_TIME_MIN
 */
double fftw_measure_runtime(fftw_plan plan)
{
     FFTW_COMPLEX *in, *out;
     fftw_time begin, end;
     double t;
     int i, iter;
     int n;

     n = plan->n;

     iter = 1;

retry:
     in = (FFTW_COMPLEX *) fftw_malloc(n * sizeof(FFTW_COMPLEX));
     out = (FFTW_COMPLEX *) fftw_malloc(n * sizeof(FFTW_COMPLEX));

     begin = fftw_get_time();
     for (i = 0; i < iter; ++i) {
	  int j;

	  /* generate random inputs */
	  for (j = 0; j < n; ++j) {
	       c_re(in[j]) = 1.0;
	       c_im(in[j]) = 32.432;
	  }

	  fftw(plan, 1, in, 1, 0, out, 1, 0);
     }
     end = fftw_get_time();

     t = fftw_time_to_sec(fftw_time_diff(end,begin));

     fftw_free(in);
     fftw_free(out);

     if (t < FFTW_TIME_MIN) {
	  iter *= 2;
	  /* 
	   * See D. E. Knuth, Structured Programming with GOTO Statements,
	   * Computing Surveys (6), December 1974, for a justification
	   * of this `goto' in the `n + 1/2' loop.
	   */
	  goto retry;
     }

     return t / (double)iter;
}

#if defined(MAC) || defined(macintosh)

/* Use Macintosh Time Manager to get the time: */

#pragma only_std_keywords off  /* make sure compiler (CW) recognizes the pascal
				  keywords that are in Timer.h */

#include <Timer.h>

#pragma only_std_keywords reset

fftw_time get_Mac_microseconds(void)
{
     fftw_time t;
     UnsignedWide microsec;	/* 
				 * microsec.lo and microsec.hi are
				 * unsigned long's, and are the two parts
				 * of a 64 bit unsigned integer 
				 */

     Microseconds(&microsec);	/* get time in microseconds */

     /* store lo and hi words into our structure: */
     t.lo = microsec.lo; t.hi = microsec.hi;

     return t;
}

fftw_time fftw_time_diff(fftw_time t1, fftw_time t2)
/* This function takes the difference t1 - t2 of two 64 bit
   integers, represented by the 32 bit lo and hi words.
   if t1 < t2, returns 0. */
{
     fftw_time diff;

     if (t1.hi < t2.hi) { /* something is wrong...t1 < t2! */
	  diff.hi = diff.lo = 0;
	  return diff;
     }
     else
	  diff.hi = t1.hi - t2.hi;

     if (t1.lo < t2.lo) {
	  if (diff.hi > 0)
	       diff.hi -= 1; /* carry */
	  else { /* something is wrong...t1 < t2! */
	       diff.hi = diff.lo = 0;
	       return diff;
	  }
     }
     
     diff.lo = t1.lo - t2.lo;

     return diff;
}

#endif

#if defined __WIN32__
#include <windows.h>

static LARGE_INTEGER gFreq;
static int gHaveHiResTimer = 0;
static int gFirstTime = 1;

unsigned long GetPerfTime(void)
{
     LARGE_INTEGER lCounter;

     if (gFirstTime) {
	  gFirstTime = 0;

	  if (QueryPerformanceFrequency(&gFreq)) {
	       gHaveHiResTimer = 1;
	  }
     }
     if (gHaveHiResTimer) {
	  QueryPerformanceCounter(&lCounter);
	  return lCounter.u.LowPart;
     } else {
#if defined(__QK_USER__)
          return (unsigned long) (dclock() * 1000000.0L)
#else
	  return (unsigned long) clock();
#endif
     }
}

double GetPerfSec(double pTime)
{
     if (gHaveHiResTimer) {
	  return pTime / gFreq.u.LowPart;	/* assumes HighPart==0 */

     } else {
	  return pTime / CLOCKS_PER_SEC;
     }
}

#endif				/* __WIN32__ */

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

/*
 * twiddle.c -- compute twiddle factors
 * These are the twiddle factors for *direct* fft.  Flip sign to get
 * the inverse
 */

/* $Id: fftw.c,v 1.3 2010-01-26 14:06:59 giannozz Exp $ */
#if defined FFTW_USING_CILK
#include <cilk.h>
#include <cilk-compat.h>
#endif

#include "fftw.h"
#include <math.h>
#include <stdlib.h>

#define FFTW_K2PI 6.2831853071795864769252867665590057683943387987502

/*
 * compute the W coefficients (that is, powers of the root of 1)
 * and store them into an array.
 */
static void fftw_compute_twiddle(int n, int r, int m, FFTW_COMPLEX *W)
{
     double twoPiOverN;
     int i, j;

     twoPiOverN = FFTW_K2PI / (double) n;
     for (i = 0; i < m; ++i)
	  for (j = 1; j < r; ++j) {
	       int k = i * (r - 1) + (j - 1);
	       c_re(W[k]) = cos(twoPiOverN * (double) i * (double) j);
	       c_im(W[k]) = -sin(twoPiOverN * (double) i * (double) j);
	  }
}

/*
 * these routines implement a simple reference-count-based 
 * management of twiddle structures
 */
static fftw_twiddle *twlist = (fftw_twiddle *) 0;
int fftw_twiddle_size = 0;	/* total allocated size, for debugging */

fftw_twiddle *fftw_create_twiddle(int n, int r, int m)
{
     fftw_twiddle *tw;
     FFTW_COMPLEX *W;

     /* lookup for this n in the twiddle list */
     for (tw = twlist; tw; tw = tw->next)
	  if (tw->n == n && tw->r == r && tw->m == m) {
	       ++tw->refcnt;
	       return tw;
	  }
     /* not found --- allocate a new struct twiddle */
     tw = (fftw_twiddle *) fftw_malloc(sizeof(fftw_twiddle));
     W = (FFTW_COMPLEX *) fftw_malloc(m * (r - 1) * sizeof(FFTW_COMPLEX));
     fftw_twiddle_size += n;

     tw->n = n;
     tw->r = r;
     tw->m = m;
     tw->twarray = W;
     tw->refcnt = 1;
     fftw_compute_twiddle(n, r, m, W);

     /* enqueue the new struct */
     tw->next = twlist;
     twlist = tw;

     return tw;
}

void fftw_destroy_twiddle(fftw_twiddle * tw)
{
     fftw_twiddle **p;
     --tw->refcnt;

     if (tw->refcnt == 0) {
	  /* remove from the list of known twiddle factors */
	  for (p = &twlist; p; p = &((*p)->next))
	       if (*p == tw) {
		    *p = tw->next;
		    fftw_twiddle_size -= tw->n;
		    fftw_free(tw->twarray);
		    fftw_free(tw);
		    return;
	       }
	  fftw_die("BUG in fftw_destroy_twiddle\n");
     }
}
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

/*
 * wisdom.c -- manage the wisdom
 */

#include "fftw.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

struct wisdom {
     int n;
     int flags;
     fftw_direction dir;
     enum fftw_node_type type;	/* this is the wisdom */
     int signature;		/* this is the wisdom */
     struct wisdom *next;
};

/* list of wisdom */
static struct wisdom *wisdom_list = (struct wisdom *) 0;

int fftw_wisdom_lookup(int n, int flags, fftw_direction dir,
		     enum fftw_node_type *type,
		     int *signature, int replacep)
{
     struct wisdom *p;

     if (!(flags & FFTW_USE_WISDOM))
	  return 0;		/* simply ignore if wisdom is disabled */

     flags |= FFTW_MEASURE; /* always use (only) wisdom from measurements */

     for (p = wisdom_list; p; p = p->next) {
	  if (p->n == n && p->flags == flags && p->dir == dir) {
	       /* found wisdom */
	       if (replacep) {
		    /* replace old wisdom with new */
		    p->type = *type;
		    p->signature = *signature;
	       } else {
		    *type = p->type;
		    *signature = p->signature;
	       }
	       return 1;
	  }
     }

     return 0;
}

void fftw_wisdom_add(int n, int flags, fftw_direction dir,
		   enum fftw_node_type type,
		   int signature)
{
     struct wisdom *p;

     if (!(flags & FFTW_USE_WISDOM))
	  return;		/* simply ignore if wisdom is disabled */

     if (!(flags & FFTW_MEASURE))
	  return;  /* only measurements produce wisdom */

     if (fftw_wisdom_lookup(n, flags, dir, &type, &signature, 1))
	  return;		/* wisdom overwrote old wisdom */

     p = (struct wisdom *) fftw_malloc(sizeof(struct wisdom));

     p->n = n;
     p->flags = flags;
     p->dir = dir;
     p->type = type;
     p->signature = signature;

     /* remember this wisdom */
     p->next = wisdom_list;
     wisdom_list = p;
}

void fftw_forget_wisdom(void)
{
     while (wisdom_list) {
	  struct wisdom *p;

	  p = wisdom_list;
	  wisdom_list = wisdom_list->next;
	  fftw_free(p);
     }
}

/*
 * user-visible routines, to convert wisdom into strings etc.
 */
#define WISDOM_FORMAT_VERSION "FFTW-1.2"

static void (*emit)(char c, void *data);

static void emit_string(char *s, void *data)
{
     while (*s) 
	  emit(*s++, data);
}

static void emit_int(int n, void *data)
{
     char buf[128];

     sprintf(buf, "%d", n);
     emit_string(buf, data);
}

/* dump wisdom in lisp-like format */
void fftw_export_wisdom(void (*emitter)(char c, void *), void *data)
{
     struct wisdom *p;

     /* install the output handler */
     emit = emitter;

     emit('(',data);
     emit_string(WISDOM_FORMAT_VERSION,data);

     for (p = wisdom_list; p; p = p->next) {
	  emit(' ',data);	/* separator to make the output nicer */
	  emit('(',data);
	  emit_int((int) p->n, data);
	  emit(' ',data);
	  emit_int((int) p->flags, data);
	  emit(' ',data);
	  emit_int((int) p->dir, data);
	  emit(' ',data);
	  emit_int((int) p->type, data);
	  emit(' ',data);
	  emit_int((int) p->signature, data);
	  emit(')',data);
     }
     emit(')',data);
}

/* input part */
static int next_char;
static int (*get_input)(void *data);
static fftw_status input_error;

static void read_char(void *data)
{
     next_char = get_input(data);
     if (next_char == 0 ||
	 next_char == EOF)
	  input_error = FFTW_FAILURE;
}

/* skip blanks, newlines, tabs, etc */
static void eat_blanks(void *data)
{
     while (isspace(next_char))
	  read_char(data);
}

static int read_int(void *data)
{
     int sign = 1;
     int n = 0;

     eat_blanks(data);
     if (next_char == '-') {
	  sign = -1;
	  read_char(data);
	  eat_blanks(data);
     }

     if (!isdigit(next_char)) {
	  /* error, no digit */
	  input_error = FFTW_FAILURE;
	  return 0;
     }

     while (isdigit(next_char)) {
	  n = n * 10 + (next_char - '0');
	  read_char(data);
     }

     return sign * n;
}

#define EXPECT(c)                     \
{				      \
     eat_blanks(data);		      \
     if (input_error == FFTW_FAILURE || \
         next_char != c)	      \
	  return FFTW_FAILURE;	      \
     read_char(data);		      \
}
				      
#define EXPECT_INT(n)                                 \
{				                      \
     n = read_int(data);	                      \
     if (input_error == FFTW_FAILURE)                 \
	  return FFTW_FAILURE;		              \
}				      
				      
#define EXPECT_STRING(s)             \
{                                    \
     char *s1 = s;		     \
     while (*s1) {		     \
	  EXPECT(*s1);		     \
	  ++s1;			     \
     }				     \
}       			      
                                      
fftw_status fftw_import_wisdom(int (*g)(void *), void *data)
{
     int n;
     int flags;
     fftw_direction dir;
     enum fftw_node_type type;
     int signature;

     get_input = g;
     input_error = FFTW_SUCCESS;

     read_char(data);

     eat_blanks(data);
     EXPECT('(');
     eat_blanks(data);
     EXPECT_STRING(WISDOM_FORMAT_VERSION);
     eat_blanks(data);

     while (next_char != ')') {
	  EXPECT('(');
	  EXPECT_INT(n);
	  EXPECT_INT(flags);
	  EXPECT_INT(dir);
	  EXPECT_INT(type);
	  EXPECT_INT(signature);
	  eat_blanks(data);
	  EXPECT(')');

	  /* the wisdom has been read properly. Add it */
	  fftw_wisdom_add(n, flags, dir, type, signature);

	  /* prepare for next morsel of wisdom */
	  eat_blanks(data);
     }

     return FFTW_SUCCESS;
}
