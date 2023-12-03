!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE upf_spinorb
  !
  !! Variables needed for calculations with spin-orbit
  !
  USE upf_kinds,   ONLY : DP
  USE upf_params,  ONLY : lmaxx, lqmax 
  !
  !! FIXME: rot_ylm could be dynamically allocated
  !
  IMPLICIT NONE
  SAVE

  LOGICAL :: is_spinorbit
  !! if .TRUE. this is a spin-orbit calculation
  !! internal flag, set by uspp_allocate, to be used only in upflib
  COMPLEX (DP) :: rot_ylm(lqmax,lqmax)
  !! transform real spherical harmonics into complex ones
  COMPLEX (DP), ALLOCATABLE :: fcoef(:,:,:,:,:)
  !! function needed to account for spinors.
  !
  INTERFACE transform_qq_so
     MODULE PROCEDURE transform_qqr_so, transform_qqc_so
  END INTERFACE transform_qq_so
CONTAINS
  !
  !----------------------------------------------------------------------
  SUBROUTINE transform_qqc_so( qq, qq_so )
  !----------------------------------------------------------------------
  ! For Berry-phase / macroscopic electric-field calculations (qq complex)
  !  
  USE uspp_param,   ONLY : nsp, upf, nhm, nh
  !
  COMPLEX(DP), INTENT(in)  :: qq(nhm, nhm, nsp)
  COMPLEX(DP), INTENT(out) :: qq_so(nhm, nhm, 4, nsp)
  !
  INTEGER :: nt, ih, jh, kh, lh, ijs, is1, is2, is
  !
  qq_so=(0.0_DP, 0.0_DP)
  DO nt = 1, nsp
    IF ( upf(nt)%tvanp ) THEN
      IF (upf(nt)%has_so) THEN
        DO ih=1,nh(nt)
          DO jh=1,nh(nt)
            DO kh=1,nh(nt)
              DO lh=1,nh(nt)
                ijs=0
                DO is1=1,2
                  DO is2=1,2
                    ijs=ijs+1
                    DO is=1,2
                      qq_so(kh,lh,ijs,nt) = qq_so(kh,lh,ijs,nt)    &
                           + qq(ih,jh,nt) * fcoef(kh,ih,is1,is,nt) &
                                          * fcoef(jh,lh,is,is2,nt)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO ih = 1, nh (nt)
          DO jh = ih, nh (nt)
             IF (is_spinorbit) THEN
                 qq_so (ih, jh, 1, nt) = qq (ih, jh, nt) 
                 qq_so (jh, ih, 1, nt) = qq_so (ih, jh, 1, nt)
                 qq_so (ih, jh, 4, nt) = qq_so (ih, jh, 1, nt)
                 qq_so (jh, ih, 4, nt) = qq_so (ih, jh, 4, nt)
             ENDIF
          ENDDO
        ENDDO
      ENDIF
    ENDIF
  ENDDO
  !
END SUBROUTINE transform_qqc_so
  !
  !----------------------------------------------------------------------
  SUBROUTINE transform_qqr_so( qq, qq_so )
  !----------------------------------------------------------------------
  !
  ! For non-Berry-phase/macroscopic electric-field calculations (qq real)
  !
  USE uspp_param,   ONLY : nsp, upf, nhm, nh
  !
  REAL(DP), INTENT(in)  :: qq(nhm, nhm, nsp)
  COMPLEX(DP), INTENT(out) :: qq_so(nhm, nhm, 4, nsp)
  !
  INTEGER :: nt, ih, jh, kh, lh, ijs, is1, is2, is
  !
  qq_so=(0.0_DP, 0.0_DP)
  DO nt = 1, nsp
    IF ( upf(nt)%tvanp ) THEN
      IF (upf(nt)%has_so) THEN
        DO ih=1,nh(nt)
          DO jh=1,nh(nt)
            DO kh=1,nh(nt)
              DO lh=1,nh(nt)
                ijs=0
                DO is1=1,2
                  DO is2=1,2
                    ijs=ijs+1
                    DO is=1,2
                      qq_so(kh,lh,ijs,nt) = qq_so(kh,lh,ijs,nt)    &
                           + qq(ih,jh,nt) * fcoef(kh,ih,is1,is,nt) &
                                          * fcoef(jh,lh,is,is2,nt)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO ih = 1, nh (nt)
          DO jh = ih, nh (nt)
             IF (is_spinorbit) THEN
                 qq_so (ih, jh, 1, nt) = qq (ih, jh, nt) 
                 qq_so (jh, ih, 1, nt) = qq_so (ih, jh, 1, nt)
                 qq_so (ih, jh, 4, nt) = qq_so (ih, jh, 1, nt)
                 qq_so (jh, ih, 4, nt) = qq_so (ih, jh, 4, nt)
             ENDIF
          ENDDO
        ENDDO
      ENDIF
    ENDIF
  ENDDO

END SUBROUTINE transform_qqr_so

END MODULE upf_spinorb

