!
! Copyright (C) 2002-2005 FPMD-CPV groups
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
      SUBROUTINE allocate_pseudo(ps,nsp,ngw,lnlx,ngh,nk)

        TYPE (pseudo) ps
        INTEGER, INTENT(IN) :: nsp,lnlx,ngh,ngw,nk
        INTEGER :: ierr

        ALLOCATE(ps%wnl(ngw,lnlx,nsp,nk), STAT=ierr)
        IF( ierr /= 0 ) CALL errore(' allocate_pseudo ', ' allocating %wnl ', ierr )
        ALLOCATE(ps%wsg(ngh,nsp), STAT=ierr)
        IF( ierr /= 0 ) CALL errore(' allocate_pseudo ', ' allocating %wsg ', ierr )

        NULLIFY(ps%ap)

        RETURN
      END SUBROUTINE allocate_pseudo

!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE deallocate_pseudo(ps)
        TYPE (pseudo) ps
        INTEGER :: ierr

        IF(ASSOCIATED(ps%wnl)) THEN
          DEALLOCATE(ps%wnl, STAT=ierr)
          IF( ierr /= 0 ) CALL errore(' deallocate_pseudo ', ' deallocating %wnl ', ierr )
        END IF
        IF(ASSOCIATED(ps%wsg)) THEN
          DEALLOCATE(ps%wsg, STAT=ierr)
          IF( ierr /= 0 ) CALL errore(' deallocate_pseudo ', ' deallocating %wsg ', ierr )
        END IF

        RETURN
      END SUBROUTINE deallocate_pseudo

!  ----------------------------------------------
!  ----------------------------------------------

      END MODULE cp_types

