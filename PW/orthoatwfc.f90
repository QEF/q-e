!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
subroutine orthoatwfc
  !-----------------------------------------------------------------------
  !
  ! This routine is meant to orthogonalize all the atomic wfcs. This is
  ! useful when we want to compute the occupation of the atomic orbitals
  ! in order to make lda+U calculations
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : nchix
  USE io_global,  ONLY : stdout
  USE io_files,   ONLY : iunat, nwordatwfc, iunigk
  USE ions_base,  ONLY : nat
  USE basis,      ONLY : natomwfc
  USE klist,      ONLY : nks, xk
  USE ldaU,       ONLY : swfcatom, U_projection
  USE wvfct,      ONLY : npwx, npw, igk, gamma_only
  USE uspp,       ONLY : nkb, vkb
  USE becmod,     ONLY : becp, rbecp
#ifdef __PARA
  USE para
#endif
  ! 
  IMPLICIT NONE
  !
  !
  integer :: ik, ibnd, info, i, j, k, na, nb, nt, isym, n, ntemp, m, &
       l, lm, ltot, ntot
  integer, allocatable ::  nwfc (:), ml (:)
  ! the k point under consideration
  ! counter on bands
  ! number of wfc of atom na
  ! corrispondence m <-> natomwfc
  integer :: nchi_, lchi_ (nchix)
  logical :: exst
  real(kind=DP) :: t0, scnds
  ! cpu time spent
  logical :: orthogonalize_wfc
     
  complex(kind=DP) :: temp, t (5)
  complex(kind=DP) , allocatable :: wfcatom (:,:), work (:,:), overlap (:,:)
  real(kind=DP) , allocatable :: e (:)

  t0 = scnds ()
  
  allocate (wfcatom( npwx, natomwfc))    
  allocate (overlap( natomwfc , natomwfc))    
  allocate (work   ( natomwfc , natomwfc))    
  allocate (e      ( natomwfc))    
  allocate (ml     ( natomwfc))    
  allocate (nwfc   ( nat))    

  if (U_projection=="file") then
     WRITE( stdout,*) 'LDA+U Projector read from file '
     return
  end if

  if (U_projection=="atomic") then
     orthogonalize_wfc = .false.
     WRITE( stdout,*) 'Atomic wfc used for LDA+U Projector are NOT orthogonalized'
  else if (U_projection=="ortho-atomic") then
     orthogonalize_wfc = .true.
     WRITE( stdout,*) 'Atomic wfc used for LDA+U Projector are orthogonalized'
     if (gamma_only) then
        WRITE( stdout,*) 'Gamma-only calculation for this case not implemented'
        stop
     end if
  else
     write( stdout,*) "U_projection_type =", U_projection
     call errore ("orthoatwfc"," this U_projection_type is not valid",1)
  end if

  ! Allocate the array becp = <beta|wfcatom>
  if ( gamma_only ) then 
     allocate (rbecp (nkb,natomwfc)) 
  else
     allocate ( becp (nkb,natomwfc)) 
  end if
  
  if (nks > 1) rewind (iunigk)
  
  do ik = 1, nks
     
     if (nks > 1) read (iunigk) npw, igk
     
     overlap(:,:) = (0.d0,0.d0)
     work(:,:) = (0.d0,0.d0)
     
     call atomic_wfc (ik, wfcatom)
     call init_us_2 (npw, igk, xk (1, ik), vkb)
     
     if ( gamma_only ) then 
        call pw_gemm ('Y', nkb, natomwfc, npw, vkb, npwx, &
             wfcatom, npwx, rbecp, nkb) 
     else
        call ccalbec (nkb, npwx, npw, natomwfc, becp, vkb, wfcatom)
     endif

     call s_psi (npwx, npw, natomwfc, wfcatom, swfcatom)

     if (orthogonalize_wfc) then
        !
        ! calculate overlap matrix
        !
        call ZGEMM ('c', 'n', natomwfc, natomwfc, npw, (1.d0, 0.d0) , &
             wfcatom, npwx, swfcatom, npwx, (0.d0, 0.d0) , overlap, natomwfc)
#ifdef __PARA
        call reduce (2 * natomwfc * natomwfc, overlap)
#endif
        !
        ! find O^-.5
        !
        call cdiagh (natomwfc, overlap, natomwfc, e, work)
        do i = 1, natomwfc
           e (i) = 1.d0 / dsqrt (e (i) )
        enddo
        do i = 1, natomwfc
           do j = i, natomwfc
              temp = (0.d0, 0.d0)
              do k = 1, natomwfc
                 temp = temp + e (k) * work (j, k) * conjg (work (i, k) )
              enddo
              overlap (i, j) = temp
              if (j.ne.i) overlap (j, i) = conjg (temp)
           enddo
        enddo
        !
        ! trasform atomic orbitals O^-.5 psi
        !
        do i = 1, npw
           work(:,1) = (0.d0,0.d0)
           call ZGEMV ('n', natomwfc, natomwfc, (1.d0, 0.d0) , overlap, &
                natomwfc, swfcatom (i, 1) , npwx, (0.d0, 0.d0) , work, 1)
           call ZCOPY (natomwfc, work, 1, swfcatom (i, 1), npwx)
        enddo
        
     end if
     
     call davcio (swfcatom, nwordatwfc, iunat, ik, 1)
     
  enddo
  deallocate (nwfc)
  deallocate (ml)
  deallocate (overlap)
  deallocate (work)
  deallocate (e)
  deallocate (wfcatom)
  if ( gamma_only ) then 
     deallocate (rbecp) 
  else
     deallocate ( becp) 
  end if
  !
  return
     
END SUBROUTINE orthoatwfc


