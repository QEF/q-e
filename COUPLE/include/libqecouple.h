/*
 * Copyright (C) 2013 Quantum ESPRESSO group
 * This file is distributed under the terms of the
 * GNU General Public License. See the file `License'
 * in the root directory of the present distribution,
 * or http://www.gnu.org/copyleft/gpl.txt .
 */

/* C/C++ interface to the codes of the Quantum ESPRESSO package */

#ifndef QE_LIBCOUPLE_H
#define QE_LIBCOUPLE_H

/* API version of the COUPLE library C interface.
 * Increment, if incompatible changes are made to the API. */

#define QE_LIBCOUPLE_API_VERSION 1

#ifdef __cplusplus
extern "C" {
#endif

/* interface to pw.x */
/* launch a pw.x-like calculation */
void c2libpwscf(int lib_comm, int nimage, int npot, int npool, int ntaskgroup,
                int nband, int ndiag, int *exit_status, char *input_file);

/* interface to cp.x */
/* launch a cp.x-like calculation */
void c2libcpv(int lib_comm, int nimage, int npot, int npool, int ntaskgroup,
              int nband, int ndiag, int *exit_status, char *input_file);

/* accessing the qmmm.f90 module */
/* pass in the inter program communicator */
void c2qmmm_mpi_config(int qmmm_mode, int inter_comm, int verb, int inter_rank);
    
#ifdef __cplusplus
}
#endif

#endif /* QE_LIBCOUPLE_H */
