!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------
!  Car-Parrinello Parallel Program
!  Carlo Cavazzoni - Gerardo Ballabio
!  SISSA, Trieste, Italy - 1997-99
!  Last modified: Sun Nov 14 07:54:47 MET 1999
!  ----------------------------------------------
!  BEGIN manual

      MODULE cp_types

!  this module contains the definitions of several TYPE structures,
!  together with their allocation/deallocation routines
!  ----------------------------------------------
!  routines in this module:
!  SUBROUTINE allocate_recvecs(igs,ngl,ngg,ngwl,ngwg,gv, nk, tk)
!  SUBROUTINE deallocate_recvecs(gv)
!  SUBROUTINE allocate_pseudo(ps,nsp,ng,ngw,lnlx,ngh,tcc)
!  SUBROUTINE deallocate_pseudo(ps)
!  SUBROUTINE allocate_phfac(eigr,nr1,nr2,nr3,nsp,nat,ngw)
!  SUBROUTINE deallocate_phfac(eigr)
!  ----------------------------------------------
!  END manual

        USE kinds
        USE pseudo_types, ONLY: pseudo_ncpp

        IMPLICIT NONE
        SAVE

        PRIVATE

!  BEGIN manual
!  TYPE DEFINITIONS

!! ...  Structure Factor
        TYPE structure_factor
          COMPLEX(dbl), POINTER :: s(:)  ! S(G) = sum_I exp(i G dot R_I)
        END TYPE

!! ...  phase factors exp(i G dot R_I)
!! ...  G = reciprocal lattice vectors
!! ...  R_I = ionic positions
        TYPE phase_factors
          COMPLEX(dbl), POINTER :: xyz(:,:)  ! exp(i G dot R_I)
          COMPLEX(dbl), POINTER :: x(:,:)    ! exp(i G_x x_I)
          COMPLEX(dbl), POINTER :: y(:,:)    ! exp(i G_y y_I)
          COMPLEX(dbl), POINTER :: z(:,:)    ! exp(i G_z z_I)
                                           !   first index: G vector
                                           !   second index: ion
        END TYPE phase_factors


!! ...  pseudopotential
        TYPE pseudo
!! ...    local part
          REAL(dbl), POINTER :: vps(:,:)     ! form factors
          REAL(dbl), POINTER :: dvps(:,:)    ! cell derivatives of form factors
          REAL(dbl), POINTER :: rhops(:,:)   ! Ewald pseudocharges form factors
!! ...    nonlocal part
          REAL(dbl), POINTER :: wsg(:,:)     ! inverse of Kleinman-Bylander
                                          !   denominators
                                          ! <Y phi | V | phi Y>**(-1)
                                          !   first index: orbital
                                          !   second index: atomic species
          REAL(dbl), POINTER :: wnl(:,:,:,:) ! Kleinman-Bylander products
                                          ! <Y phi | V | exp(i(k+G) dot r)>
                                          !   first index: G vector
                                          !   second index: orbital
                                          !   third index: atomic species
                                          !   fourth index: k point
!! ...    core corrections
          LOGICAL, POINTER :: tnlcc (:)
          REAL(dbl), POINTER :: rhoc1(:,:)  ! correction to pseudopotential
          REAL(dbl), POINTER :: rhocp(:,:)  ! cell derivative

          TYPE (pseudo_ncpp), POINTER :: ap(:)

        END TYPE pseudo


!! ...  arrays of reciprocal space vectors (G and k+G)
        TYPE recvecs
          LOGICAL gzero
          LOGICAL hspace
          INTEGER gstart ! 2 if hg_l(1) = 0, 1 otherwise 
          INTEGER ng_l   ! local number of G vectors
          INTEGER ng_g   ! global number of G vectors

          INTEGER ngw_l
          INTEGER ngw_g

!! ...    quantities related to G vectors
          REAL(dbl), POINTER :: hg_l(:)        ! length squared of G vectors

          REAL(dbl), POINTER :: gx_l(:,:)   ! components of G vectors 
                                            !   first index: x,y,z
                                            !   second index: G vector

          INTEGER, POINTER :: mill(:,:)     ! indices for FFT
                                            !   first index: x,y,z
                                            !   second index: G vector

          INTEGER, POINTER :: ig(:)     ! Global indices of the G vectors
                                        !   index: G vector

          REAL(dbl), POINTER :: khgcutz_l(:,:) ! smooth cutoff factor index: G vector
                                               !   first index: G vector
                                               !   second index: k point

!! ...    quantities related to k+G vectors

          REAL(dbl), POINTER :: kg_mask_l(:,:) ! cutoff mask
          REAL(dbl), POINTER :: khg_l(:,:)     ! length squared of k+G
          REAL(dbl), POINTER :: kgx_l(:,:,:)   ! components of k+G
                                            !   first index: G vector
                                            !   second index: x,y,z
                                            !   third index: k point
          REAL(dbl) :: bi1(3), bi2(3), bi3(3)  ! Initial reciprocal lattice base (used to count the g's)
          REAL(dbl) :: b1(3), b2(3), b3(3)  ! Actual reciprocal lattice base (used to determine g2 )
        END TYPE recvecs

!  ----------------------------------------------
!  END manual

        PUBLIC :: pseudo, recvecs, phase_factors, pseudo_ncpp
        PUBLIC :: allocate_pseudo, deallocate_pseudo
        PUBLIC :: allocate_recvecs, deallocate_rvecs
        PUBLIC :: allocate_phfac, deallocate_phfac

!  end of module-scope declarations
!  ----------------------------------------------

      CONTAINS

!  subroutines


        SUBROUTINE allocate_recvecs(gv,ngl,ngg,ngwl,ngwg, tk,nk)

          INTEGER, INTENT(IN) :: ngl,ngg,ngwl,ngwg,nk
          LOGICAL, INTENT(IN) :: tk
          TYPE (recvecs), INTENT(OUT) :: gv
          INTEGER :: ierr

            gv%gstart = -100000
            gv%gzero  = .TRUE.

            gv%ng_l = ngl
            gv%ng_g = ngg
            gv%ngw_l = ngwl
            gv%ngw_g = ngwg

            !ALLOCATE(gv%hg_l(ngl), STAT=ierr)
            !IF( ierr /= 0 ) CALL errore(' allocate_recvecs ', ' allocating %hg_l ', ierr )
            NULLIFY( gv%hg_l )

            !ALLOCATE(gv%gx_l(ngl,3), STAT=ierr)
            !IF( ierr /= 0 ) CALL errore(' allocate_recvecs ', ' allocating %gx_l ', ierr )
            NULLIFY( gv%gx_l )

            IF (tk) THEN
              gv%hspace = .FALSE.
            ELSE 
              gv%hspace = .TRUE.
            END IF

            IF (tk) THEN
              ALLOCATE(gv%kg_mask_l(ngwl,nk), STAT=ierr)
            ELSE 
              ALLOCATE(gv%kg_mask_l(1,1), STAT=ierr)
            END IF
            IF( ierr /= 0 ) CALL errore(' allocate_recvecs ', ' allocating %kg_mask_l ', ierr )

            ALLOCATE(gv%kgx_l(3,ngwl,nk), STAT=ierr)
            IF( ierr /= 0 ) CALL errore(' allocate_recvecs ', ' allocating %kgx_l ', ierr )
            ALLOCATE(gv%khg_l(ngwl,nk), STAT=ierr)
            IF( ierr /= 0 ) CALL errore(' allocate_recvecs ', ' allocating %khg_l ', ierr )
            ALLOCATE(gv%khgcutz_l(ngwl,nk), STAT=ierr)
            IF( ierr /= 0 ) CALL errore(' allocate_recvecs ', ' allocating %khgcutz_l ', ierr )

            NULLIFY( gv%mill )
            NULLIFY( gv%ig )

            !ALLOCATE(gv%mill(3,ngl), STAT=ierr)
            !IF( ierr /= 0 ) CALL errore(' allocate_recvecs ', ' allocating %mill ', ierr )
            !ALLOCATE(gv%ig(ngl), STAT=ierr)
            !IF( ierr /= 0 ) CALL errore(' allocate_recvecs ', ' allocating %ig ', ierr )

          RETURN
        END SUBROUTINE allocate_recvecs

!  ----------------------------------------------
!  ----------------------------------------------

        SUBROUTINE deallocate_rvecs(gv)
          TYPE (recvecs), INTENT(INOUT) :: gv
          INTEGER :: ierr
          gv%gstart = 0
          gv%ng_l = 0
          gv%ng_g = 0
          gv%ngw_l = 0
          gv%ngw_g = 0

          IF(ASSOCIATED(gv%hg_l)) THEN
            NULLIFY( gv%hg_l )
          END IF
          IF(ASSOCIATED(gv%gx_l)) THEN
            NULLIFY( gv%gx_l )
          END IF
          IF(ASSOCIATED(gv%kg_mask_l)) THEN
            DEALLOCATE(gv%kg_mask_l, STAT=ierr)
            IF( ierr /= 0 ) CALL errore(' deallocate_rvecs ', ' deallocating %kg_mask_l ', ierr )
          END IF
          IF(ASSOCIATED(gv%kgx_l)) THEN
            DEALLOCATE(gv%kgx_l, STAT=ierr)
            IF( ierr /= 0 ) CALL errore(' deallocate_rvecs ', ' deallocating %kgx_l ', ierr )
          END IF
          IF(ASSOCIATED(gv%khg_l)) THEN
            DEALLOCATE(gv%khg_l, STAT=ierr)
            IF( ierr /= 0 ) CALL errore(' deallocate_rvecs ', ' deallocating %khg_l ', ierr )
          END IF
          IF(ASSOCIATED(gv%khgcutz_l)) THEN
            DEALLOCATE(gv%khgcutz_l, STAT=ierr)
            IF( ierr /= 0 ) CALL errore(' deallocate_rvecs ', ' deallocating %khgcutz_l ', ierr )
          END IF
          IF(ASSOCIATED(gv%mill)) THEN
            NULLIFY( gv%mill )
          END IF
          IF(ASSOCIATED(gv%ig)) THEN
            NULLIFY( gv%ig )
          END IF
          RETURN
        END SUBROUTINE deallocate_rvecs

!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE allocate_pseudo(ps,nsp,ng,ngw,nk,lnlx,ngh,tcc)

        TYPE (pseudo) ps
        INTEGER, INTENT(IN) :: ng,nsp,lnlx,ngh,ngw,nk
        LOGICAL, INTENT(IN) :: tcc(:)
        INTEGER :: ierr

        ALLOCATE(ps%vps(ng,nsp), STAT=ierr)
        IF( ierr /= 0 ) CALL errore(' allocate_pseudo ', ' allocating %vps ', ierr )
        ALLOCATE(ps%dvps(ng,nsp), STAT=ierr)
        IF( ierr /= 0 ) CALL errore(' allocate_pseudo ', ' allocating %dvps ', ierr )
        ALLOCATE(ps%rhops(ng,nsp), STAT=ierr)
        IF( ierr /= 0 ) CALL errore(' allocate_pseudo ', ' allocating %rhops ', ierr )
        ALLOCATE(ps%wnl(ngw,lnlx,nsp,nk), STAT=ierr)
        IF( ierr /= 0 ) CALL errore(' allocate_pseudo ', ' allocating %wnl ', ierr )
        ALLOCATE(ps%wsg(ngh,nsp), STAT=ierr)
        IF( ierr /= 0 ) CALL errore(' allocate_pseudo ', ' allocating %wsg ', ierr )

        NULLIFY(ps%ap)

        IF(ANY(tcc)) THEN 
          ALLOCATE(ps%tnlcc(nsp), STAT=ierr)
          IF( ierr /= 0 ) CALL errore(' allocate_pseudo ', ' allocating %tnlcc ', ierr )
          ALLOCATE(ps%rhoc1(ng,nsp), STAT=ierr)
          IF( ierr /= 0 ) CALL errore(' allocate_pseudo ', ' allocating %rhoc1 ', ierr )
          ALLOCATE(ps%rhocp(ng,nsp), STAT=ierr)
          IF( ierr /= 0 ) CALL errore(' allocate_pseudo ', ' allocating %rhocp ', ierr )
          ps%tnlcc = tcc
        ELSE
          ALLOCATE(ps%tnlcc(1), STAT=ierr)
          IF( ierr /= 0 ) CALL errore(' allocate_pseudo ', ' allocating %tnlcc ', ierr )
          ALLOCATE(ps%rhoc1(1,1), STAT=ierr)
          IF( ierr /= 0 ) CALL errore(' allocate_pseudo ', ' allocating %rhoc1 ', ierr )
          ALLOCATE(ps%rhocp(1,1), STAT=ierr)
          IF( ierr /= 0 ) CALL errore(' allocate_pseudo ', ' allocating %rhocp ', ierr )
          ps%tnlcc = .FALSE.
        END IF

        RETURN
      END SUBROUTINE allocate_pseudo

!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE deallocate_pseudo(ps)
        TYPE (pseudo) ps
        INTEGER :: ierr

        IF(ASSOCIATED(ps%vps)) THEN
          DEALLOCATE(ps%vps, STAT=ierr)
          IF( ierr /= 0 ) CALL errore(' deallocate_pseudo ', ' deallocating %vps ', ierr )
        END IF
        IF(ASSOCIATED(ps%dvps)) THEN
          DEALLOCATE(ps%dvps, STAT=ierr)
          IF( ierr /= 0 ) CALL errore(' deallocate_pseudo ', ' deallocating %dvps ', ierr )
        END IF
        IF(ASSOCIATED(ps%rhops)) THEN
          DEALLOCATE(ps%rhops, STAT=ierr)
          IF( ierr /= 0 ) CALL errore(' deallocate_pseudo ', ' deallocating %rhops ', ierr )
        END IF
        IF(ASSOCIATED(ps%wnl)) THEN
          DEALLOCATE(ps%wnl, STAT=ierr)
          IF( ierr /= 0 ) CALL errore(' deallocate_pseudo ', ' deallocating %wnl ', ierr )
        END IF
        IF(ASSOCIATED(ps%wsg)) THEN
          DEALLOCATE(ps%wsg, STAT=ierr)
          IF( ierr /= 0 ) CALL errore(' deallocate_pseudo ', ' deallocating %wsg ', ierr )
        END IF
        IF(ASSOCIATED(ps%tnlcc)) THEN
          DEALLOCATE(ps%tnlcc, STAT=ierr)
          IF( ierr /= 0 ) CALL errore(' deallocate_pseudo ', ' deallocating %tnlcc ', ierr )
        END IF
        IF(ASSOCIATED(ps%rhoc1)) THEN
          DEALLOCATE(ps%rhoc1, STAT=ierr)
          IF( ierr /= 0 ) CALL errore(' deallocate_pseudo ', ' deallocating %rhoc1 ', ierr )
        END IF
        IF(ASSOCIATED(ps%rhocp)) THEN
          DEALLOCATE(ps%rhocp, STAT=ierr)
          IF( ierr /= 0 ) CALL errore(' deallocate_pseudo ', ' deallocating %rhocp ', ierr )
        END IF

        RETURN
      END SUBROUTINE deallocate_pseudo

!  ----------------------------------------------
!  ----------------------------------------------

      SUBROUTINE allocate_phfac(eigr,nr1,nr2,nr3,nsp,nat,ngw,ng)

        INTEGER, INTENT(IN) :: nr1,nr2,nr3,nsp,nat,ngw,ng
        TYPE (phase_factors) :: eigr

        INTEGER :: is, ierr

        ALLOCATE(eigr%xyz(ngw,nat), STAT=ierr)
        IF( ierr /= 0 ) CALL errore(' allocate_phfac ', ' allocating %xyz ', ierr )
        ALLOCATE(eigr%x(-nr1:nr1,nat), STAT=ierr)
        IF( ierr /= 0 ) CALL errore(' allocate_phfac ', ' allocating %x ', ierr )
        ALLOCATE(eigr%y(-nr2:nr2,nat), STAT=ierr)
        IF( ierr /= 0 ) CALL errore(' allocate_phfac ', ' allocating %y ', ierr )
        ALLOCATE(eigr%z(-nr3:nr3,nat), STAT=ierr)
        IF( ierr /= 0 ) CALL errore(' allocate_phfac ', ' allocating %z ', ierr )
        RETURN
      END SUBROUTINE allocate_phfac

!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE deallocate_phfac(eigr)

        TYPE (phase_factors) eigr

        INTEGER :: is, ierr

        IF(ASSOCIATED(eigr%xyz)) THEN
          DEALLOCATE(eigr%xyz, STAT=ierr)
          IF( ierr /= 0 ) CALL errore(' deallocate_phfac ', ' deallocating %xyz ', ierr )
        END IF
        IF(ASSOCIATED(eigr%x)) THEN
          DEALLOCATE(eigr%x, STAT=ierr)
          IF( ierr /= 0 ) CALL errore(' deallocate_phfac ', ' deallocating %x ', ierr )
        END IF
        IF(ASSOCIATED(eigr%y)) THEN
          DEALLOCATE(eigr%y, STAT=ierr)
          IF( ierr /= 0 ) CALL errore(' deallocate_phfac ', ' deallocating %y ', ierr )
        END IF
        IF(ASSOCIATED(eigr%z)) THEN
          DEALLOCATE(eigr%z, STAT=ierr)
          IF( ierr /= 0 ) CALL errore(' deallocate_phfac ', ' deallocating %z ', ierr )
        END IF
        RETURN
      END SUBROUTINE deallocate_phfac

!  ----------------------------------------------
!  ----------------------------------------------

      END MODULE cp_types

