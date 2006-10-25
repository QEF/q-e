!
! Copyright (C) 2006 QUANTUM-ESPRESSO
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE write_states(nrzp,n2d,norb,norbf,nchan,nrx,nry,ounit, &
                                                  funz0,vec,norm_flag)
!
! Using the basis functions obtained by scatter_forw (funz0), and
! the coefficients obtained in compbs or in transmit (vec) it computes
! the states on a three dimensional grid and writes them in files
! readable by the post-processing programs. 
! If norm_flag is .true. the states are normalized in such a way
! that their square modulus integrated over all volume is one.
!
#include "f_defs.h"
 USE kinds, ONLY : DP
 USE noncollin_module, ONLY : npol
 USE cond,      ONLY : ngper, newbg, nl_2d
 USE pwcom,     ONLY : omega, ibrav, celldm, gcutm, dual, ecutwfc, at
 USE ions_base, ONLY : ityp, zv, tau, atm

 IMPLICIT NONE
 INTEGER :: nrzp,n2d,norb,norbf,nchan,nrx,nry,ounit 
 COMPLEX(DP) :: funz0(n2d, 2*n2d+norbf*npol, nrzp), vec(2*n2d+npol*norb,nchan)
 LOGICAL :: norm_flag
 INTEGER :: ik, ig, ir, mu, ichan, ios
 REAL(DP) :: DDOT
 COMPLEX(DP), PARAMETER :: one=(1.d0,0.d0), zero=(0.d0,0.d0)
 CHARACTER(LEN=50) :: filename
 CHARACTER(LEN=75) :: title_here

 INTEGER :: ix, jx
 COMPLEX(DP), ALLOCATABLE :: funr(:), fung(:), kfunz(:,:,:)

 REAL(DP), ALLOCATABLE :: aux_plot(:), norm(:)
!
!    Calculate the functions in all z points
!
 ALLOCATE(funr(nrx*nry))
 ALLOCATE(fung(npol*ngper))
 ALLOCATE(kfunz(n2d,nchan,nrzp))
 DO ik=1,nrzp
    CALL ZGEMM('n', 'n', n2d, nchan, 2*n2d+npol*norb, one, funz0(1,1,ik), &
                n2d, vec, 2*n2d+npol*norb, zero, kfunz(1,1,ik), n2d)
 END DO
 IF (norm_flag) THEN
    ALLOCATE (norm(nchan))
    norm=0.d0
    DO ichan=1,nchan
       DO ik=1,nrzp
          norm(ichan) = norm(ichan) + &
                     DDOT(2*n2d,kfunz(1,ichan,ik),1,kfunz(1,ichan,ik),1)
       END DO
    END DO
    norm=norm*omega/nrzp
    norm=1.d0/SQRT(norm)
    DO ik=1,nrzp
       DO ichan=1,nchan
          CALL DSCAL(2*n2d,norm(ichan),kfunz(1,ichan,ik),1)
       END DO 
    END DO
    DEALLOCATE (norm)
 END IF
!
!    Transform and save all propagating functions. One function after the other.
!
 ALLOCATE(aux_plot(nrx*nry*nrzp))
 DO ichan=1,nchan
    ounit = 34
    IF (ichan>9) THEN
      write(filename,'("wfc_n.",i2)') ichan
    ELSE
      write(filename,'("wfc_n.",i1)') ichan
    ENDIF
    filename=TRIM(filename)
    OPEN (UNIT=ounit, FILE=filename, FORM='formatted', &
                           STATUS='unknown', ERR=100, IOSTAT=ios)
100 CALL errore('write_states','opening file'//filename,ABS(ios))

    aux_plot=(0.d0,0.d0)
    DO ik=1,nrzp
       fung=(0.d0,0.d0)
       DO ig=1,npol*ngper
          DO mu=1,n2d
             fung(ig)=fung(ig) + kfunz(mu,ichan,ik) * newbg(ig,mu)
          END DO
       END DO

       funr=(0.d0,0.d0)
       DO ig=1,ngper*npol
          funr(nl_2d(ig))=fung(ig)
       END DO

       CALL cft3(funr,nrx,nry,1,nrx,nry,1,1)

       DO ix=1,nrx
          DO jx=1,nry
             ir=ix + (jx - 1) * nrx + (ik - 1) * nrx * nry
             ig=ix+(jx-1)*nrx
             aux_plot(ir)=REAL(funr(ig))**2+AIMAG(funr(ig))**2
          END DO
       END DO
    END DO
    CLOSE(ounit)

    title_here='written by pwcond'

    CALL plot_io (TRIM(filename), title_here, nrx, nry, nrzp,       &
                  nrx, nry, nrzp, 0, 0, ibrav, celldm, at, gcutm,   &
                  dual, ecutwfc, 7, atm, ityp, zv, tau, aux_plot, + 1)
 END DO
 DEALLOCATE(aux_plot)
 DEALLOCATE(kfunz)
 DEALLOCATE(funr)
 DEALLOCATE(fung)

 RETURN
END SUBROUTINE write_states
