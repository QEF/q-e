/* ///////////////////////////////////////////////////////////////////////////
// @file    secd.c
// @author  Michael Holst
// @brief   Timing routines.
// @version $Id: secd.c,v 1.1 2008-06-11 14:24:03 degironc Exp $
// @attention
// @verbatim
//
// PMG -- Parallel algebraic MultiGrid
// Copyright (c) 1994-2006.  Michael Holst.
//
// Michael Holst <mholst@math.ucsd.edu>
// University of California, San Diego
// Department of Mathematics, 5739 AP&M
// 9500 Gilman Drive, Dept. 0112
// La Jolla, CA 92093-0112 USA                                                  
// http://math.ucsd.edu/~mholst
//
// This file is part of PMG.
//
// PMG is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// PMG is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with PMG; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
//
// Linking PMG statically or dynamically with other modules is making a
// combined work based on PMG. Thus, the terms and conditions of the GNU
// General Public License cover the whole combination.
// 
// SPECIAL GPL EXCEPTION
// In addition, as a special exception, the copyright holders of PMG
// give you permission to combine the PMG program with free software
// programs and libraries that are released under the GNU LGPL or with
// code included in releases of ISIM, PMV, PyMOL, SMOL, VMD, and Vision.
// Such combined software may be linked with PMG and redistributed together 
// in original or modified form as mere aggregation without requirement that 
// the entire work be under the scope of the GNU General Public License.
// This special exception permission is also extended to any software listed
// in the SPECIAL GPL EXCEPTION clauses by the FEtk and APBS libraries.
// 
// Note that people who make modified versions of PMG are not obligated
// to grant this special exception for their modified versions; it is
// their choice whether to do so. The GNU General Public License gives
// permission to release a modified version without this exception; this
// exception also makes it possible to release a modified version which
// carries forward this exception.
//
// @endverbatim
// //////////////////////////////////////////////////////////////////////// */

/* ******************************************************************* */
/* purpose:                                                            */
/*                                                                     */
/*    Returns the time in seconds used by the process.                 */
/*    The system call is "getrusage", which returns a                  */
/*    struct containing various pieces of resource usage information.  */
/*                                                                     */
/* Machine:                                                            */
/*                                                                     */
/*    most UNIX machines                                               */
/*                                                                     */
/* author:  michael holst                                              */
/* ******************************************************************* */
#include <stdio.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/resource.h>
double tsecnd_()
{
        static struct rusage temp;
        long foo, foo1;

        getrusage(RUSAGE_SELF,&temp)    ;                /* get clock */
        foo     = temp.ru_utime.tv_sec  ;                /* seconds */
        foo1    = temp.ru_utime.tv_usec ;                /* uSecs */
        return ((double)foo + (double)foo1*0.000001) ;   /* milliseconds */
}
