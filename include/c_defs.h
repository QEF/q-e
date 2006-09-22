/*
Copyright (C) 2002-2006 Quantum-ESPRESSO group
This file is distributed under the terms of the
GNU General Public License. See the file `License'
in the root directory of the present distribution,
or http://www.gnu.org/copyleft/gpl.txt .
*/

/* Machine-dependent Fortran to C calling convention
   redefine C symbols so that Fortran finds them

   C routine 'name' is called as 
      F77_FUNC('name','NAME')
   if 'name' does not conatin an underscore; by
      F77_FUNC_('name','NAME')
   if it does . The following defines F77_FUNC and F77_FUNC_

   XLF (Aix, Mac OS-X), HP-UX:
      lowercase, with no added underscores 
   G95, EKOPath, Alpha Linux:
      lowercase, with one added underscore if the name does
      not contain underscores, with two if it does
   Most other cases: 
      lowercase, with one added underscore 

*/

#if defined (__ABSOFT)

/* convert to capital letters, no added underscores */

#define F77_FUNC(name,NAME) NAME
#define F77_FUNC_(name,NAME) NAME

#elif ( defined (__XLF) || defined (__HP) )

/* convert to lowercase letters, no added underscores */

#define F77_FUNC(name,NAME) name
#define F77_FUNC_(name,NAME) name

#elif ( defined (__G95) || defined (__EKO) || (defined __ALPHA && defined __LINUX64) )

/* convert to lowercase letters, one added underscore if none in the name,
   two added underscores if one is already present in the name */

#define F77_FUNC(name,NAME) name ## _
#define F77_FUNC_(name,NAME) name ## __

#else

/* convert to lowercase letters, one added underscore */

#define F77_FUNC(name,NAME) name ## _
#define F77_FUNC_(name,NAME) name ## _

#endif
