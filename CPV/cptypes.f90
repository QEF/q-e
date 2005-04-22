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
!  ----------------------------------------------
!  END manual

        USE kinds
        USE pseudo_types, ONLY: pseudo_ncpp

        IMPLICIT NONE
        SAVE

        PRIVATE

!  BEGIN manual
!  TYPE DEFINITIONS

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


!  ----------------------------------------------
!  END manual

        PUBLIC :: pseudo, pseudo_ncpp
        PUBLIC :: allocate_pseudo, deallocate_pseudo

!  end of module-scope declarations
!  ----------------------------------------------

      CONTAINS

!  subroutines


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

      END MODULE cp_types

