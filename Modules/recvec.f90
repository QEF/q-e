!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------------=!
   MODULE gvecw
!=----------------------------------------------------------------------------=!
     USE kinds

     IMPLICIT NONE
     SAVE

     ! ...   G vectors less than the wave function cut-off ( ecutwfc )
     INTEGER ngw    ! local number of G vectors
     INTEGER ngwt   ! in parallel execution global number of G vectors,
                    ! in serial execution this is equal to ngw
     INTEGER ngwl   ! number of G-vector shells up to ngw
     INTEGER ngwx   ! maximum local number of G vectors
     INTEGER ng0    ! first G-vector with nonzero modulus
                    ! needed in the parallel case (G=0 is on one node only!)

     REAL(dbl) :: ecutw

!=----------------------------------------------------------------------------=!
   END MODULE gvecw
!=----------------------------------------------------------------------------=!

!=----------------------------------------------------------------------------=!
   MODULE gvecp
!=----------------------------------------------------------------------------=!
     USE kinds

     IMPLICIT NONE
     SAVE

     ! ...   G vectors less than the potential cut-off ( ecutrho )
     INTEGER ngm    ! local number of G vectors
     INTEGER ngmt   ! in parallel execution global number of G vectors,
                    ! in serial execution this is equal to ngm
     INTEGER ngml   ! number of G-vector shells up to ngw
     INTEGER ngmx   ! maximum local number of G vectors

     REAL(dbl) :: ecutp

!=----------------------------------------------------------------------------=!
   END MODULE gvecp
!=----------------------------------------------------------------------------=!

!=----------------------------------------------------------------------------=!
   MODULE gvecs
!=----------------------------------------------------------------------------=!
     USE kinds

     IMPLICIT NONE
     SAVE

     ! ...   G vectors less than the smooth grid cut-off ( ? )
     INTEGER ngs    ! local number of G vectors
     INTEGER ngst   ! in parallel execution global number of G vectors,
                    ! in serial execution this is equal to ngw
     INTEGER ngsl   ! number of G-vector shells up to ngw
     INTEGER ngsx   ! maximum local number of G vectors

     INTEGER, ALLOCATABLE :: nps(:), nms(:)

     REAL(dbl) :: ecuts

!=----------------------------------------------------------------------------=!
   END MODULE gvecs
!=----------------------------------------------------------------------------=!

!=----------------------------------------------------------------------------=!
   MODULE gvecb
!=----------------------------------------------------------------------------=!
     USE kinds

     IMPLICIT NONE
     SAVE

     ! ...   G vectors less than the box grid cut-off ( ? )
     INTEGER ngb    ! local number of G vectors
     INTEGER ngbt   ! in parallel execution global number of G vectors,
                    ! in serial execution this is equal to ngw
     INTEGER ngbl   ! number of G-vector shells up to ngw
     INTEGER ngbx   ! maximum local number of G vectors

     REAL(dbl), ALLOCATABLE :: gb(:), gxb(:,:), gxnb(:,:), glb(:)
     INTEGER, ALLOCATABLE :: npb(:), nmb(:), iglb(:), in1pb(:), in2pb(:), in3pb(:)

     REAL(dbl) :: ecutb

!=----------------------------------------------------------------------------=!
   END MODULE gvecb
!=----------------------------------------------------------------------------=!


!=----------------------------------------------------------------------------=!
   MODULE reciprocal_vectors
!=----------------------------------------------------------------------------=!

     USE kinds
     USE gvecp
     USE gvecb
     USE gvecs
     USE gvecw

     IMPLICIT NONE
     SAVE

     ! ...   declare module-scope variables

     LOGICAL gzero   ! .TRUE. if the first G vectors on this processor is
                     ! the null G vector ( i.e. |G| == 0 )
     INTEGER gstart  ! index of the first G vectors whose module is greather
                     ! than 0 . 
                     ! gstart = 2 when gzero == .TRUE., gstart = 1 otherwise 

     REAL(dbl) :: tpiba, tpiba2

!=----------------------------------------------------------------------------=!
   END MODULE reciprocal_vectors
!=----------------------------------------------------------------------------=!
