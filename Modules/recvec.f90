!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!=----------------------------------------------------------------------------=!
   MODULE reciprocal_vectors
!=----------------------------------------------------------------------------=!

     USE kinds

     IMPLICIT NONE
     SAVE

     ! ...   declare module-scope variables

     ! ...   G vectors less than the wave function cut-off ( ecutwfc )
     INTEGER ngwt   ! global number of G vectors
     INTEGER ngwl   ! local number of G vectors
     INTEGER ngwlx  ! maximum local number of G vectors

     ! ...   G vectors less than the potential cut-off ( ecutrho )
     INTEGER ngmt   ! same as above
     INTEGER ngml
     INTEGER ngmlx

     ! ...   G vectors less than the smooth grid cut-off ( ? )
     INTEGER ngst   ! same as above
     INTEGER ngsl
     INTEGER ngslx

     LOGICAL gzero   ! .TRUE. if the first G vectors on this processor is
                     ! the null G vector ( i.e. |G| == 0 )
     INTEGER gstart  ! index of the first G vectors whose module is greather
                     ! than 0 . 
                     ! gstart = 2 when gzero == .TRUE., gstart = 1 otherwise 

!=----------------------------------------------------------------------------=!
   END MODULE reciprocal_vectors
!=----------------------------------------------------------------------------=!
