!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"

!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------
!  Car-Parrinello Parallel Program
!  Carlo Cavazzoni - Gerardo Ballabio
!  SISSA, Trieste, Italy - 1997-99
!  Last modified: Fri Oct  8 15:08:56 MDT; 1999
!  ----------------------------------------------
!  BEGIN manual

!=----------------------------------------------------------------------------=!
      MODULE reciprocal_space_mesh
!=----------------------------------------------------------------------------=!

!  (describe briefly what this module does...)
!  ----------------------------------------------
!  routines in this module:
!  SUBROUTINE gmeshset(ngw,ng)
!  INTEGER FUNCTION owner_of_gvec(ig)
!  SUBROUTINE newg(gv,kp,htm1)
!  ----------------------------------------------
!  END manual

! ...   declare included modules

        USE kinds

        USE cell_base, ONLY: tpiba, tpiba2

        IMPLICIT NONE
        SAVE

        !! ...    quantities related to k+G vectors

        REAL(dbl), ALLOCATABLE :: gkcutz_l(:,:) ! smooth cutoff factor index: G vector
                                                ! first index: G vector
                                                ! second index: k point

        REAL(dbl), ALLOCATABLE :: gkmask_l(:,:) ! cutoff mask
        REAL(dbl), ALLOCATABLE :: gk_l(:,:)     ! length squared of k+G
        REAL(dbl), ALLOCATABLE :: gkx_l(:,:,:)  ! components of k+G
                                                ! first index: G vector
                                                ! second index: x,y,z
                                                ! third index: k point


        PRIVATE

        PUBLIC :: recvecs_units, newgk, gmeshinfo, gindex_closeup
        PUBLIC :: gkmask_l, gkcutz_l, gkx_l, gk_l

! ...   end of module-scope declarations
!  ----------------------------------------------

!=----------------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------------=!


   SUBROUTINE recvecs_units( alat )
     USE constants, ONLY: tpi
     REAL(dbl) :: alat
     IF( alat == 0.0d0 ) &
       CALL errore(' recvecs_units ', ' alat is zero ! ', 1 )
     tpiba = tpi / alat
     tpiba2 = tpiba ** 2
     RETURN
   END SUBROUTINE


!  ----------------------------------------------



   SUBROUTINE gmeshinfo( )
          
!  (describe briefly what this routine does...)
!  ----------------------------------------------

! ...     declare modules
          USE mp_global, ONLY: nproc, mpime, group
          USE io_global, ONLY: ionode, ionode_id, stdout
          USE mp, ONLY: mp_max, mp_gather
          USE brillouin, ONLY: kp
          USE reciprocal_vectors, only: &
              ngw_g  => ngwt,   &
              ngw_l  => ngw ,   &
              ngw_lx => ngwx,   &
              ng_g   => ngmt,   &
              ng_l   => ngm ,   &
              ng_lx  => ngmx,   &
              ig_l2g, &
              mill_l

          IMPLICIT NONE

          INTEGER :: ip, ng_snd(6), ng_rcv(6,nproc)

! ... end of declarations
!  ----------------------------------------------

! ...     diagnostics
          ng_snd(1) = ng_g
          ng_snd(2) = ng_l
          ng_snd(3) = ng_lx
          ng_snd(4) = ngw_g
          ng_snd(5) = ngw_l
          ng_snd(6) = ngw_lx
          CALL mp_gather(ng_snd, ng_rcv, ionode_id, group)
          IF(ionode) THEN
            WRITE( stdout,*)
            WRITE( stdout,*) '  Reciprocal Space Mesh'
            WRITE( stdout,*) '  ---------------------'
            WRITE( stdout,1000)
            DO ip = 1, nproc
              WRITE( stdout,1010) ip, ng_rcv(1,ip), ng_rcv(2,ip), ng_rcv(3,ip), &
                ng_rcv(4,ip), ng_rcv(5,ip), ng_rcv(6,ip)
            END DO
          END IF

1000      FORMAT(16X,'Large Mesh Number of G',15X,'Small Mesh Number of G' &
          ,/,'   PE       Global       Local   Max Local' &
                 ,'       Global       Local   Max Local' )
 
1010      FORMAT( I5,1X,3I12,1X,3I12 )

          RETURN 

   END SUBROUTINE  

        
!  ----------------------------------------------


   INTEGER FUNCTION owner_of_gvec(mill)

     USE stick_base, ONLY: stown => sticks_owner

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: mill(:)

     owner_of_gvec = stown( mill(1), mill(2) )

     RETURN
   END FUNCTION owner_of_gvec


!  ----------------------------------------------


      SUBROUTINE newgk( kp, htm1 )

!  this routine computes the squared modulus, and Ciartesian components 
!  of G vectors, from their Miller indices and current cell shape. 
!  G vectors are expressed in units of 2*pi/alat.
!  In generic-k-points calculations, squared modulus and components of
!  k+G vectors are computed too, together with cutoff masks
!  ----------------------------------------------

! ... declare modules
      USE gvecw,              ONLY: ngw, ecfix, ecsig, ecutz, gcutw
      USE gvecp,              ONLY: ngm
      USE reciprocal_vectors, only: g, gx
      USE cell_base,          ONLY: tpiba2
      USE brillouin,          ONLY: kpoints

      IMPLICIT NONE

! ... declare subroutine argumentsz
      TYPE (kpoints), INTENT(IN)  :: kp
      REAL(dbl), INTENT(IN) :: htm1(3,3)
      REAL(dbl) :: xk(3, SIZE(kp%xk, 2) )

! ... declare external function
      REAL(dbl) :: erf, erfc
      EXTERNAL erf, erfc

! ... declare other variables
      INTEGER   :: ig, i, j, k, ik, nkp
      INTEGER   :: isign

! ... end of declarations
!  ----------------------------------------------

      nkp = kp%nkpt
      IF( .NOT. ALLOCATED( gk_l ) )       ALLOCATE( gk_l( ngw, nkp ) )
      IF( .NOT. ALLOCATED( gkx_l ) )      ALLOCATE( gkx_l( 3, ngw, nkp ) )
      IF( .NOT. ALLOCATED( gkcutz_l ) )   ALLOCATE( gkcutz_l( ngw, nkp ) )
      IF( .NOT. ALLOCATED( gkmask_l ) )   ALLOCATE( gkmask_l  ( ngw, nkp ) )
      IF( SIZE( gk_l ) /= ( ngw * nkp ) ) THEN
        DEALLOCATE( gk_l )
        DEALLOCATE( gkx_l )
        DEALLOCATE( gkcutz_l )
        DEALLOCATE( gkmask_l )
        ALLOCATE( gk_l( ngw, nkp ) )
        ALLOCATE( gkx_l( 3, ngw, nkp ) )
        ALLOCATE( gkcutz_l( ngw, nkp ) )
        ALLOCATE( gkmask_l( ngw, nkp ) )
      END IF

      IF( kp%scheme == 'gamma' ) THEN

        DO ik = 1, kp%nkpt
          gk_l( 1:ngw ,ik)       = g( 1:ngw )
          gkx_l( 1:3, 1:ngw, ik) = gx( 1:3, 1:ngw )
        END DO

      ELSE

        DO ik = 1, kp%nkpt

          xk(:,ik) = kp%xk(:,ik)

! ...     compute components of G
          DO ig = 1, ngw

! ...       compute components of G+k

            gkx_l(1,ig,ik) = gx(1,ig) + xk(1,ik)
            gkx_l(2,ig,ik) = gx(2,ig) + xk(2,ik)
            gkx_l(3,ig,ik) = gx(3,ig) + xk(3,ik)

! ...       compute squared length

            gk_l(ig,ik) = gkx_l(1,ig,ik)**2 + gkx_l(2,ig,ik)**2 + gkx_l(3,ig,ik)**2

            !  compute cutoff mask for k+G

            IF( gk_l(ig, ik) <= gcutw ) THEN
              gkmask_l(ig, ik) = 1.0d0
            ELSE 
              gkmask_l(ig, ik) = 0.0d0
            END IF

          END DO
        END DO
      END IF

      DO ik = 1, kp%nkpt
        IF( ecutz > 0.0d0 ) THEN
          DO ig = 1, ngw
            ! ... compute smooth cutoff G+k vectors
            gkcutz_l(ig,ik) = erf( ( tpiba2 * gk_l(ig,ik) - ecfix ) / ecsig )
            gkcutz_l(ig,ik) = gk_l(ig,ik) + ecutz / tpiba2 * ( 1.0d0 + gkcutz_l(ig,ik) )
          END DO
        ELSE
          gkcutz_l(1:ngw,ik) = gk_l(1:ngw,ik)
        END IF
      END DO

      RETURN
      END SUBROUTINE newgk

!=----------------------------------------------------------------------------=!

      SUBROUTINE gindex_closeup
        IMPLICIT NONE
      RETURN
      END SUBROUTINE
     

!=----------------------------------------------------------------------------=!
      END MODULE reciprocal_space_mesh
!=----------------------------------------------------------------------------=!

