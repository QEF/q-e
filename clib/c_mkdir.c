/*
  Copyright (C) 2003 PWSCF group
  This file is distributed under the terms of the
  GNU General Public License. See the file `License'
  in the root directory of the present distribution,
  or http://www.gnu.org/copyleft/gpl.txt .
*/

#include <stdio.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "cp.h"

int C_MKDIR( char * dirname )
{
   mode_t mode = 0777 ;
   return mkdir( dirname, mode );
}
