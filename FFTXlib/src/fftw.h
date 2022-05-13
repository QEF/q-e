#ifndef FFTW_H
#define FFTW_H

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

/* ----------- */

#ifdef __cplusplus
extern "C" {
#endif				/* __cplusplus */

#define c_re(c)  ((c).re)
#define c_im(c)  ((c).im)

typedef enum {
     FFTW_FORWARD = -1, FFTW_BACKWARD = 1
} fftw_direction;

/*********************************************
 * Success or failure status
 *********************************************/

typedef enum {
     FFTW_SUCCESS = 0, FFTW_FAILURE = -1
} fftw_status;


enum fftw_node_type {
	FFTW_NOTW, FFTW_TWIDDLE, FFTW_GENERIC
};

/* flags for the planner */
#define  FFTW_ESTIMATE (0)
#define  FFTW_MEASURE  (1)

#define FFTW_IN_PLACE (8)
#define FFTW_USE_WISDOM (16)

void fftw_die(char* s);

#ifdef __cplusplus
}                               /* extern "C" */
#endif				/* __cplusplus */

#endif				/* FFTW_H */

