!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------

subroutine orthoatwfc  
  !-----------------------------------------------------------------------
  !
  ! This routine is meant to orthogonalize all the atomic wfcs. This is
  ! useful when we want to compute the occupation of the atomic orbitals
  ! in order to make lda+U calculations
  !
#include "machine.h"
  use pwcom  
  use becmod  
  use allocate
#ifdef PARA
  use para
#endif
  implicit none
  integer :: ik, ibnd, info, i, j, k, na, nb, nt, isym, n, ntemp, m, &
       l, lm, ltot, ntot
  integer, pointer ::  nwfc (:), ml (:)
  ! the k point under consideration
  ! counter on bands
  ! number of wfc of atom na
  ! corrispondence m <-> natomwfc
  integer :: nchi_, lchi_ (nchix)  
  logical :: exst  
  real(kind=DP) :: t0, scnds  
  ! cpu time spent
  logical :: orthogonalize_wfc

  complex(kind=DP) :: ZDOTC, temp, t (5)
  complex(kind=DP) , pointer :: wfcatom (:,:), work (:,:), overlap (:,:)
  real(kind=DP) , pointer :: e (:)  


  t0 = scnds ()  

  call mallocate(wfcatom, npwx, natomwfc)  
  call mallocate(overlap, natomwfc , natomwfc)  
  call mallocate(work   , natomwfc , natomwfc)  
  call mallocate(e      , natomwfc)  
  call mallocate(ml     , natomwfc)  
  call mallocate(nwfc   , nat)  

  orthogonalize_wfc = .false.
  if (orthogonalize_wfc) then
     write (6,*) 'Atomic wfc used in LDA+U are orthogonalized'
  else
     write (6,*) 'Atomic wfc used in LDA+U are NOT orthogonalized'
  end if

  if (nks.gt.1) rewind (iunigk)  

  do ik = 1, nks  

     if (nks.gt.1) read (iunigk) npw, igk  

     overlap(:,:) = (0.d0,0.d0)
     work(:,:) = (0.d0,0.d0)

     call atomic_wfc (ik, wfcatom)  
     call init_us_2 (npw, igk, xk (1, ik), vkb)  

     call ccalbec (nkb, npwx, npw, natomwfc, becp, vkb, wfcatom)  

     call s_psi (npwx, npw, natomwfc, wfcatom, swfcatom)  

     if (orthogonalize_wfc) then
        !
        ! calculate overlap matrix
        !
        call ZGEMM ('c', 'n', natomwfc, natomwfc, npw, (1.d0, 0.d0) , &
             wfcatom, npwx, swfcatom, npwx, (0.d0, 0.d0) , overlap, natomwfc)
#ifdef PARA
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
           call setv (2 * natomwfc, 0.d0, work, 1)  
           call ZGEMV ('n', natomwfc, natomwfc, (1.d0, 0.d0) , overlap, &
                natomwfc, swfcatom (i, 1) , npwx, (0.d0, 0.d0) , work, 1)
           call ZCOPY (natomwfc, work, 1, swfcatom (i, 1), npwx)  
        enddo

     end if

     call davcio (swfcatom, nwordatwfc, iunat, ik, 1)  

  enddo
  call mfree (nwfc)  
  call mfree (ml)  
  call mfree (overlap)  
  call mfree (work)  
  call mfree (e)  
  call mfree (wfcatom)  
  return  

end subroutine orthoatwfc

