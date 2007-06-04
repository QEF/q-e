/*
  Copyright (C) 2007 PWSCF group
  This file is distributed under the terms of the
  GNU General Public License. See the file `License'
  in the root directory of the present distribution,
  or http://www.gnu.org/copyleft/gpl.txt .
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include "c_defs.h"

void F77_FUNC(remove_stack_limit,REMOVE_STACK_LIMIT) (void) {
 
  struct rlimit rlim = { RLIM_INFINITY, RLIM_INFINITY };

  if ( setrlimit(RLIMIT_STACK, &rlim) == -1 ) {
    perror("Cannot set stack size to infinity");
    /* exit(1); */
  }
}
