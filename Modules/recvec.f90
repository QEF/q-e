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
     INTEGER :: ngw  = 0  ! local number of G vectors
     INTEGER :: ngwt = 0  ! in parallel execution global number of G vectors,
                       ! in serial execution this is equal to ngw
     INTEGER :: ngwl = 0  ! number of G-vector shells up to ngw
     INTEGER :: ngwx = 0  ! maximum local number of G vectors
     INTEGER :: ng0  = 0  ! first G-vector with nonzero modulus
                       ! needed in the parallel case (G=0 is on one node only!)

     REAL(dbl) :: ecutw = 0.0d0
     REAL(dbl) :: gcutw = 0.0d0

     !   values for costant cut-off computations

     REAL(dbl) :: ecfix = 0.0d0     ! value of the constant cut-off
     REAL(dbl) :: gcfix = 0.0d0     ! ecfix / tpiba**2
     REAL(dbl) :: ecutz = 0.0d0     ! height of the penalty function (above ecfix)
     REAL(dbl) :: gcutz = 0.0d0
     REAL(dbl) :: ecsig = 0.0d0     ! spread of the penalty function around ecfix
     REAL(dbl) :: gcsig = 0.0d0
     LOGICAL   :: tecfix = .FALSE.  ! .TRUE. if constant cut-off is in use

     ! augmented cut-off for k-point calculation

     REAL(dbl) :: ekcut = 0.0d0  
     REAL(dbl) :: gkcut = 0.0d0
    

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
     INTEGER :: ngm  = 0  ! local number of G vectors
     INTEGER :: ngmt = 0  ! in parallel execution global number of G vectors,
                       ! in serial execution this is equal to ngm
     INTEGER :: ngml = 0  ! number of G-vector shells up to ngw
     INTEGER :: ngmx = 0  ! maximum local number of G vectors

     REAL(dbl) :: ecutp = 0.0d0
     REAL(dbl) :: gcutp = 0.0d0

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
     INTEGER :: ngs  = 0  ! local number of G vectors
     INTEGER :: ngst = 0  ! in parallel execution global number of G vectors,
                       ! in serial execution this is equal to ngw
     INTEGER :: ngsl = 0  ! number of G-vector shells up to ngw
     INTEGER :: ngsx = 0  ! maximum local number of G vectors

     INTEGER, ALLOCATABLE :: nps(:), nms(:)

     REAL(dbl) :: ecuts = 0.0d0
     REAL(dbl) :: gcuts = 0.0d0

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
     INTEGER :: ngb  = 0  ! local number of G vectors
     INTEGER :: ngbt = 0  ! in parallel execution global number of G vectors,
                       ! in serial execution this is equal to ngw
     INTEGER :: ngbl = 0  ! number of G-vector shells up to ngw
     INTEGER :: ngbx = 0  ! maximum local number of G vectors

     REAL(dbl), ALLOCATABLE :: gb(:), gxb(:,:), gxnb(:,:), glb(:)
     INTEGER, ALLOCATABLE :: npb(:), nmb(:), iglb(:), in1pb(:), in2pb(:), in3pb(:)

     REAL(dbl) :: ecutb = 0.0d0
     REAL(dbl) :: gcutb = 0.0d0

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

     LOGICAL :: gzero  = .TRUE.   ! .TRUE. if the first G vectors on this processor is
                                  ! the null G vector ( i.e. |G| == 0 )
     INTEGER :: gstart = 2        ! index of the first G vectors whose module is greather
                                  ! than 0 . 
                                  ! gstart = 2 when gzero == .TRUE., gstart = 1 otherwise 

     REAL(dbl) :: tpiba  = 0.0d0
     REAL(dbl) :: tpiba2 = 0.0d0


     !     G^2 in increasing order (in units of tpiba2=(2pi/a)^2)
     !
     REAL(dbl), ALLOCATABLE, TARGET :: g(:) 

     !     shells of G^2
     !
     REAL(dbl), ALLOCATABLE, TARGET :: gl(:) 

     !     G-vectors cartesian components ( units tpiba =(2pi/a)  )
     !
     REAL(dbl), ALLOCATABLE, TARGET :: gx(:,:) 

     !     g2_g    = all G^2 in increasing order, replicated on all procs
     !
     REAL(dbl), ALLOCATABLE, TARGET :: g2_g(:)

     !     mill_g  = miller index of G vecs (increasing order), replicated on all procs
     !
     INTEGER, ALLOCATABLE, TARGET :: mill_g(:,:)

     !     mill_l  = miller index of G vecs local to the processors
     !
     INTEGER, ALLOCATABLE, TARGET :: mill_l(:,:)

     !     ig_l2g  = "l2g" means local to global, this array convert a local
     !               G-vector index into the global index, in other words
     !               the index of the G-v. in the overall array of G-vectors
     !
     INTEGER, ALLOCATABLE, TARGET :: ig_l2g(:)

     !     igl = index of the g-vector shells
     !
     INTEGER, ALLOCATABLE, TARGET :: igl(:)

     !     bi  = base vector used to generate the reciprocal space
     !
     REAL(dbl) :: bi1(3) = (/ 0.0d0, 0.0d0, 0.0d0 /)
     REAL(dbl) :: bi2(3) = (/ 0.0d0, 0.0d0, 0.0d0 /)
     REAL(dbl) :: bi3(3) = (/ 0.0d0, 0.0d0, 0.0d0 /)

!=----------------------------------------------------------------------------=!
   END MODULE reciprocal_vectors
!=----------------------------------------------------------------------------=!
