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
!  SUBROUTINE gindexset(gv,nr1,nr2,nr3)
!  INTEGER FUNCTION owner_of_gvec(ig)
!  SUBROUTINE newg(gv,kp,htm1)
!  ----------------------------------------------
!  END manual

! ...   declare included modules

        USE kinds
        USE cp_types

        USE cell_base, ONLY: tpiba, tpiba2

        IMPLICIT NONE
        SAVE

        PRIVATE

        PUBLIC :: recvecs_units, newg, gindexset, gmeshinfo, gindex_closeup

! ...   end of module-scope declarations
!  ----------------------------------------------

!=----------------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------------=!

!  subroutines


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
!  ----------------------------------------------
        SUBROUTINE gindexset(gv, mill, b1, b2, b3)

!  (describe briefly what this routine does...)
!  ----------------------------------------------

! ...     declare modules
          USE mp_global, ONLY: mpime
          USE io_global, ONLY: stdout
          USE stick, ONLY: dfftp, fft_dlay_descriptor
          USE stick_base, ONLY: stown => sticks_owner
          USE reciprocal_vectors, only: &
              ngw_g  => ngwt,   &
              ngw_l  => ngw ,   &
              ngw_lx => ngwx,   &
              ng_g   => ngmt,   &
              ng_l   => ngm ,   &
              ng_lx  => ngmx,   &
              gstart, &
              gzero,  &
              ig_l2g, &
              mill_l, g, gx
          

          IMPLICIT NONE
       
! ...     declare subroutine arguments
          TYPE (recvecs) gv
          INTEGER, INTENT(IN) :: mill(:,:)
          REAL(dbl), INTENT(IN) :: b1(3), b2(3), b3(3) 

! ...     declare other variables
          INTEGER ig, ng, i, j, k, is

! ... end of declarations
!  ----------------------------------------------

          ng    = 0

          gv%gzero  = gzero
          gv%gstart = gstart

          gv%bi1 = b1
          gv%bi2 = b2
          gv%bi3 = b3
          gv%b1 = b1
          gv%b2 = b2
          gv%b3 = b3

          gv%hg_l => g
          gv%ig   => ig_l2g( 1:ng_l )
          gv%mill => mill_l
          gv%gx_l => gx

          IF( ng_l /= gv%ng_l ) THEN
            WRITE( stdout,*) ' MSG: gv%ng_l = ',gv%ng_l,' ng = ', ng_l 
            CALL errore(' gindexset ',' inconsistent ng ', 2)
          END IF

          ng = ng_l

          DO ig = 1, ng

            CALL miller2nxh( gv%mill(1,ig), gv%mill(2,ig), gv%mill(3,ig), &
              i, j, k, dfftp%nr1, dfftp%nr2, dfftp%nr3 )
            is = dfftp%isind( i + ( j - 1 ) * dfftp%nr1x )
            gv%ig2st( 1, ig ) = k
            gv%ig2st( 2, ig ) = is

            IF( is == 0 ) &
              CALL errore( ' gindexset ', ' wrong stick index ', 1 )

            IF( gv%hspace ) THEN
              CALL miller2indx( gv%mill(1,ig), gv%mill(2,ig), gv%mill(3,ig), &
                i, j, k, dfftp%nr1, dfftp%nr2, dfftp%nr3 )
              is = dfftp%isind( i + ( j - 1 ) * dfftp%nr1x )
              gv%ig2st( 3, ig ) = k
              gv%ig2st( 4, ig ) = is
            END IF

            IF( is == 0 ) &
              CALL errore( ' gindexset ', ' wrong stick index ', 2 )
            
          END DO

          IF( gv%gzero .and. ( .not. ( gv%gstart == 2 ) ) ) THEN
            CALL errore(' gindexset ',' gzero and gstart are inconsistent ', 3)
          END IF
          IF( gv%gstart < 1 .or. gv%gstart > 2 ) THEN
            CALL errore(' gindexset ',' gstart out of range ', 4)
          END IF

 200      FORMAT(3I5)

          RETURN
        END SUBROUTINE gindexset
        
!  ----------------------------------------------
!  ----------------------------------------------
        INTEGER FUNCTION owner_of_gvec(mill)

!  (describe briefly what this routine does...)
!  ----------------------------------------------

! ...     declare modules
          USE stick_base, ONLY: stown => sticks_owner

          IMPLICIT NONE

! ...     declare function arguments
          INTEGER, INTENT(IN) :: mill(:)

! ... end of declarations
!  ----------------------------------------------

          owner_of_gvec = stown( mill(1), mill(2) )

          RETURN
        END FUNCTION owner_of_gvec

!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE newg(gv,kp,htm1)

!  this routine computes the squared modulus, Cartesian components and
!  smooth cutoff masks of G vectors, from their Miller indices and the
!  current cell shape. G vectors are expressed in units of 2*pi/alat.
!  In generic-k-points calculations, squared modulus and components of
!  k+G vectors are computed too, together with cutoff masks
!  ----------------------------------------------

! ... declare modules
      USE gvecw, ONLY: tecfix, gcfix, gcsig, gcutz, gcutw
      USE cell_module, ONLY: alat
      USE brillouin, ONLY: kpoints
        USE reciprocal_vectors, only: &
              ngw_g  => ngwt,   &
              ngw_l  => ngw ,   &
              ngw_lx => ngwx,   &
              ng_g   => ngmt,   &
              ng_l   => ngm ,   &
              ng_lx  => ngmx,   &
              gstart, &
              gzero,  &
              ig_l2g, &
              mill_l

      IMPLICIT NONE

! ... declare subroutine argumentsz
      TYPE (recvecs), INTENT(INOUT) :: gv
      TYPE (kpoints), INTENT(IN)  :: kp
      REAL(dbl), INTENT(IN) :: htm1(3,3)
      REAL(dbl) :: xk(3, SIZE(kp%xk, 2) )

! ... declare external function
      REAL(dbl) :: erf, erfc
      EXTERNAL erf, erfc

! ... declare other variables
      INTEGER ig,i,j,k,ik
      INTEGER isign
      REAL(dbl) :: b1(3), b2(3), b3(3)

! ... end of declarations
!  ----------------------------------------------

      b1 = htm1(:,1)
      b2 = htm1(:,2)
      b3 = htm1(:,3)
      gv%b1 = htm1(:,1)
      gv%b2 = htm1(:,2)
      gv%b3 = htm1(:,3)

      DO ig = 1, gv%ng_l

        i = gv%mill(1,ig)
        j = gv%mill(2,ig)
        k = gv%mill(3,ig)

! ...   compute components of G
        gv%gx_l(1,ig) = b1(1)*i + b2(1)*j + b3(1)*k
        gv%gx_l(2,ig) = b1(2)*i + b2(2)*j + b3(2)*k
        gv%gx_l(3,ig) = b1(3)*i + b2(3)*j + b3(3)*k

        gv%gx_l(1,ig) = gv%gx_l(1,ig) * alat
        gv%gx_l(2,ig) = gv%gx_l(2,ig) * alat
        gv%gx_l(3,ig) = gv%gx_l(3,ig) * alat

! ...   compute squared length
        gv%hg_l(ig) = gv%gx_l(1,ig)**2 + gv%gx_l(2,ig)**2 + gv%gx_l(3,ig)**2

      END DO

      IF( kp%scheme == 'gamma' ) THEN

        DO ik = 1, kp%nkpt
          gv%khg_l( 1:SIZE(gv%khg_l,1) ,ik)      = gv%hg_l( 1:SIZE(gv%khg_l,1) )
          gv%kgx_l( 1:3, 1:SIZE(gv%khg_l,1), ik) = gv%gx_l( 1:3, 1:SIZE(gv%khg_l,1) )
        END DO

      ELSE

        DO ik = 1, kp%nkpt

! ...     Bring k-points in the unscaled reciprocal space
          xk(:,ik) = kp%xk(:,ik)
          !xk(1,ik) = htm1(1,1)*kp%xk(1,ik) + htm1(1,2)*kp%xk(2,ik) + htm1(1,3)*kp%xk(3,ik)
          !xk(1,ik) = xk(1,ik) * alat
          !xk(2,ik) = htm1(2,1)*kp%xk(1,ik) + htm1(2,2)*kp%xk(2,ik) + htm1(2,3)*kp%xk(3,ik)
          !xk(2,ik) = xk(2,ik) * alat
          !xk(3,ik) = htm1(3,1)*kp%xk(1,ik) + htm1(3,2)*kp%xk(2,ik) + htm1(3,3)*kp%xk(3,ik)
          !xk(3,ik) = xk(3,ik) * alat
          !WRITE( stdout,*) alat, xk(1,ik), xk(2,ik), xk(3,ik)

! ...     compute components of G
          DO ig = 1, SIZE( gv%khg_l, 1 )

! ...       compute components of G+k
            gv%kgx_l(1,ig,ik) = gv%gx_l(1,ig) + xk(1,ik)
            gv%kgx_l(2,ig,ik) = gv%gx_l(2,ig) + xk(2,ik)
            gv%kgx_l(3,ig,ik) = gv%gx_l(3,ig) + xk(3,ik)

! ...       compute squared length
            gv%khg_l(ig,ik) = gv%kgx_l(1,ig,ik)**2 + gv%kgx_l(2,ig,ik)**2 + gv%kgx_l(3,ig,ik)**2

! ...       compute cutoff mask for k+G
            IF( gv%khg_l(ig, ik) .LE. gcutw) THEN
              gv%kg_mask_l(ig, ik) = 1.0d0
            ELSE 
              gv%kg_mask_l(ig, ik) = 0.0d0
            END IF

          END DO
        END DO
      END IF

      IF(tecfix) THEN
        DO ik = 1, kp%nkpt
          DO ig = 1, SIZE( gv%khgcutz_l, 1 ) 
! ...       compute smooth cutoff G+k vectors
            gv%khgcutz_l(ig,ik) = erf((gv%khg_l(ig,ik) - gcfix)/gcsig)
            gv%khgcutz_l(ig,ik) = gv%khg_l(ig,ik) + gcutz * ( 1.0d0 + gv%khgcutz_l(ig,ik))
          END DO
        END DO
      END IF


      RETURN
      END SUBROUTINE newg


      SUBROUTINE gindex_closeup
        IMPLICIT NONE
      RETURN
      END SUBROUTINE
     

!=----------------------------------------------------------------------------=!
      END MODULE reciprocal_space_mesh
!=----------------------------------------------------------------------------=!

