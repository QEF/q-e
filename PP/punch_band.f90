!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine do_bands (nodenumber)
  !-----------------------------------------------------------------------
  use pwcom
  use becmod
  use io

  implicit none
  character (len=3)  :: nodenumber
  character (len=14) :: filband
  integer :: ios
  namelist / inputpp / tmp_dir, prefix, filband
  !
  nd_nmbr = nodenumber
  !
  !   set default values for variables in namelist
  !
  prefix = ' '
  tmp_dir = './'
  filband = ' '
  !
  read (5, inputpp, err = 200, iostat = ios)
200 call errore ('projwave', 'reading inputpp namelist', abs (ios) )
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  call read_file
  call openfil
  call init_us_1
  !
  call punch_band (filband)
  !
  return
end subroutine do_bands
!
!-----------------------------------------------------------------------
subroutine punch_band (filband)
  !-----------------------------------------------------------------------
  !
  !    This routine writes the band energies on a file. The routine orders
  !    the eigenvalues using the overlap of the eigenvectors to give
  !    an estimate crossing and anticrossing of the bands. This simplified
  !    method works in many, but not in all the cases.
  !
  !
#include "machine.h"

  use pwcom
  use becmod

  implicit none
  character (len=14) :: filband
  real(kind=DP) :: proold, modulo
  ! the best overlap product
  ! the x coordinate in k space
  complex(kind=DP) :: ZDOTC,  pro, cgracsc
  ! scalar product function
  ! the product of wavefunctions
  ! scalar product with the S matrix

  complex(kind=DP), allocatable :: psiold (:,:), old (:), actual (:), &
       becpold (:,:)
  ! the space used to save the eigenf
  ! the old testing wavefunction
  ! the testing wavefunction
  ! products of wavefunctions and bet

  integer :: ibnd, jbnd, ik, ikb, ig, npwold, ios
  ! if 1 the band has been already us
  ! counter on bands
  ! counter on bands
  ! counter on k points
  ! counter of beta functions
  ! counter on g vectors
  ! save the igk
  ! save the number of plane waves
  ! index of changes
  ! index of changes
  ! used to control I/O status
  integer, allocatable :: ok (:), igkold (:), il (:), ilm (:)


  if (filband.eq.' ') return
  iunpun = 18
  open (unit = iunpun, file = filband, status = 'unknown', form = &
       'formatted', err = 100, iostat = ios)
100 call errore ('punch_band', 'Opening filband file', abs (ios) )
  rewind (iunpun)
  allocate (psiold( npwx, nbnd))    
  allocate (old(   ngm))    
  allocate (actual(ngm))    
  allocate (becpold (nkb , nbnd))    
  allocate (igkold( npwx))    
  allocate (ok ( nbnd))    
  allocate (il ( nbnd))    
  allocate (ilm( nbnd))    

  do ik = 1, nks
     !
     !    prepare the indices of this k point
     !
     call gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, &
          igk, g2kin)
     !
     !   reads the eigenfunctions
     !
     call davcio (evc, nwordwfc, iunwfc, ik, - 1)
     if (ik.eq.1) then
        !
        !     The first eigenfunctions are saved as they are
        !
        do ibnd = 1, nbnd
           !
           !     the order is the increasing energy in this case
           !
           il (ibnd) = ibnd
           do ig = 1, npw
              psiold (ig, ibnd) = evc (ig, ibnd)
           enddo
        enddo
        do ig = 1, npw
           igkold (ig) = igk (ig)
        enddo
        npwold = npw
        modulo = wk (ik) * nks / 2.d0
        write(iunpun, '(14x,3f7.4)') xk(1,ik),xk(2,ik),xk(3,ik)
        write (iunpun, '(10f8.3)')  (et (il (ibnd) , ik) &
             * rytoev, ibnd = 1, nbnd)
        !
        !   The bec function for this k point are computed
        !
        call init_us_2 (npw, igk, xk (1, ik), vkb)

        call ccalbec (nkb, npwx, npw, nbnd, becpold, vkb, evc)
     else
        !
        !    here we are at a generic step, not the first, no eigenfunction
        !    has been already chosen
        !
        do ibnd = 1, nbnd
           ok (ibnd) = 0
        enddo
        call init_us_2 (npw, igk, xk (1, ik), vkb)

        call ccalbec (nkb, npwx, npw, nbnd, becp, vkb, evc)
        do ibnd = 1, nbnd
           call setv (2 * ngm, 0.d0, old, 1)
           do ig = 1, npwold
              old (igkold (ig) ) = psiold (ig, ibnd)
           enddo
           proold = 0.d0
           do jbnd = 1, nbnd
              if (ok (jbnd) .eq.0) then
                 call setv (2 * ngm, 0.d0, actual, 1)
                 do ig = 1, npw
                    actual (igk (ig) ) = evc (ig, jbnd)
                 enddo
                 pro = cgracsc (nkb, becp (1, jbnd), becpold (1, ibnd), &
                      nhm, ntyp, nh, qq, nat, ityp, ngm, actual, old, tvanp)
                 if (abs (pro) .gt.proold) then
                    ilm (ibnd) = jbnd
                    proold = abs (pro)
                 endif
              endif
           enddo
           ok (ilm (ibnd) ) = 1
        enddo
        !
        !   Now the order of the new eigenfunctions has been established,
        !   prepare the next k point
        !
        do ibnd = 1, nbnd
           il (ibnd) = ilm (ibnd)
           do ig = 1, npw
              psiold (ig, ibnd) = evc (ig, il (ibnd) )
           enddo
           !
           !   copy the becp in the becpold
           !
           do ikb = 1, nkb
              becpold (ikb, ibnd) = becp (ikb, il (ibnd) )
           enddo
        enddo
        do ig = 1, npw
           igkold (ig) = igk (ig)
        enddo
        npwold = npw
        !
        !     When a band calculation is performed the weight of the k point is
        !     used as the coordinate in k space
        !
        modulo = wk (ik) * nks / 2.d0

        write(iunpun, '(14x,3f7.4)') xk(1,ik),xk(2,ik),xk(3,ik)
        write (iunpun, '(10f8.3)') (et (il (ibnd) , ik) &
             * rytoev, ibnd = 1, nbnd)
     endif

  enddo
  deallocate(ilm)
  deallocate(il)
  deallocate(ok)
  deallocate(igkold)
  deallocate(becpold)
  deallocate(actual)
  deallocate(old)
  deallocate(psiold)

  close (iunpun)
  return
end subroutine punch_band
