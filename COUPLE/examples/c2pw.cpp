//
// Copyright (C) 2013 Quantum ESPRESSO group
// This file is distributed under the terms of the
// GNU General Public License. See the file `License'
// in the root directory of the present distribution,
// or http://www.gnu.org/copyleft/gpl.txt .
//

#include <mpi.h>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include "libqecouple.h"

// ... Test program for Q-E library interface
int main(int argc, char **argv)
{
    int retval, pw_comm, ncpu, key, me;
    char input[81] = { ' ', '\0' };
    MPI_Comm new_comm;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&ncpu);
    MPI_Comm_rank(MPI_COMM_WORLD,&me);

    // parse command line flags.
    int i=1;
    int nimage=1, npots=1, npools=1, ntg=1, nband=1, ndiag=1, nres=0;
    while (i < argc-1) {
        if (strncmp("-i",argv[i],2) == 0) {
            ++i;
            strncpy(input, argv[i], 80);
            input[80] = '\0';
            ++i;
            continue;
        }
        if (strncmp("-ni",argv[i],3) == 0) {
            ++i;
            nimage=std::atoi(argv[i]);
            ++i;
            continue;
        }
        if (strncmp("-npot",argv[i],5) == 0) {
            ++i;
            npots=std::atoi(argv[i]);
            ++i;
            continue;
        }
        if ((strncmp("-nk",argv[i],3) == 0)
            || (strncmp("-npoo",argv[i],5) == 0)) {
            ++i;
            ndiag=std::atoi(argv[i]);
            ++i;
            continue;
        }
        if (strncmp("-nt",argv[i],3) == 0) {
            ++i;
            ntg=std::atoi(argv[i]);
            ++i;
            continue;
        }
        if (strncmp("-nb",argv[i],3) == 0) {
            ++i;
            nband=std::atoi(argv[i]);
            ++i;
            continue;
        }
        if ((strncmp("-nd",argv[i],3) == 0)
            || (strncmp("-no",argv[i],3) == 0)
            || (strcmp("-nproc_diag",argv[i]) == 0)
            || (strcmp("-nproc_ortho",argv[i]) == 0)) {
            ++i;
            ndiag=std::atoi(argv[i]);
            ++i;
            continue;
        }
        if (strncmp("-nr",argv[i],3) == 0) {
            ++i;
            nres=std::atoi(argv[i]);
            ++i;
            continue;
        }
        std::cerr << "usage: " << argv[0] << " -flag1 <value1> -flag2 <value2>\n"
                  << std::endl;
        return -1;
    }

    if (i != argc) {
        std::cerr << "usage: " << argv[0] << " -flag1 <value1> -flag2 <value2>\n"
                  << std::endl;
        return -1;
    }

    // Create new C-style communicator and convert to Fortran
    key = MPI_UNDEFINED;
    if (me < (ncpu - nres))
        key = 1;
    
    MPI_Comm_split(MPI_COMM_WORLD, key, me, &new_comm);

    if (new_comm != MPI_COMM_NULL) {
        pw_comm = MPI_Comm_c2f(new_comm);
        // call Q-E
        c2libpwscf(pw_comm,nimage,npots,npools,ntg,nband,ndiag,&retval,input);
  
        std::cout << " rank " << me << " return value is: " << retval << std::endl;
    } else {
        std::cout << " rank " << me << " of " << ncpu -1 << " is reserved" << std::endl;
        retval = 0;
    }
    

    MPI_Finalize();
    return retval;
}
