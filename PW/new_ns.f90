!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------

subroutine new_ns
  !-----------------------------------------------------------------------
  !
  ! This routine computes the new value for ns (the occupation numbers of
  ! ortogonalized atomic wfcs).
  ! These quantities are defined as follows: ns_{I,s,m1,m2} = \sum_{k,v}
  ! f_{kv} <\fi^{at}_{I,m1}|\psi_{k,v,s}><\psi_{k,v,s}|\fi^{at}_{I,m2}>
  !
#include "machine.h"
  use pwcom
  use io
#ifdef __PARA
  use para
#endif
  implicit none
  integer :: ik, ibnd, is, i, na, nb, nt, isym, n, counter, m1, m2, &
       m0, m00, l, ldim
  integer, allocatable ::  offset (:)
  ! counter on k points
  !    "    "  bands
  !    "    "  spins
  ! offset of d electrons of atom d
  ! in the natomwfc ordering
  real(kind=DP) , allocatable :: nr (:,:,:,:)
  real(kind=DP) ::  t0, scnds
  ! cpu time spent

  complex(kind=DP) :: ZDOTC
  complex(kind=DP) , allocatable :: proj(:,:)

  real(kind=DP) :: psum

  t0 = scnds ()  
  ldim = 2 * Hubbard_lmax + 1
  allocate( offset(nat), proj(natomwfc,nbnd), nr(nat,nspin,ldim,ldim) )  
  !
  ! D_Sl for l=1, l=2 and l=3 are already initialized, for l=0 D_S0 is 1
  !
  counter = 0  
  do na = 1, nat  
     nt = ityp (na)  
     do n = 1, nchi (nt)  
        if (oc (n, nt) .gt.0.d0.or..not.newpseudo (nt) ) then  
           l = lchi (n, nt)  
           if (l.eq.Hubbard_l(nt)) offset (na) = counter  
           counter = counter + 2 * l + 1  
        endif
     enddo

  enddo

  if (counter.ne.natomwfc) call errore ('new_ns', 'nstart<>counter', 1)
  nr    (:,:,:,:) = 0.d0
  nsnew (:,:,:,:) = 0.d0
  !
  !    we start a loop on k points
  !

  if (nks.gt.1) rewind (iunigk)

  do ik = 1, nks
     if (lsda) current_spin = isk(ik)
     if (nks.gt.1) read (iunigk) npw, igk
     call davcio (evc, nwordwfc, iunwfc, ik, - 1)

     call davcio (swfcatom, nwordatwfc, iunat, ik, - 1)
     !
     ! make the projection
     !
     do ibnd = 1, nbnd
        do i = 1, natomwfc
           proj (i, ibnd) = ZDOTC (npw, swfcatom (1, i), 1, evc (1, ibnd), 1)
        enddo
     enddo
#ifdef __PARA
     call reduce (2 * natomwfc * nbnd, proj)
#endif
     !
     ! compute the occupation numbers (the quantities n(m1,m2)) of the
     ! atomic orbitals
     !
     do na = 1, nat  
        nt = ityp (na)  
        if (Hubbard_U(nt).ne.0.d0 .or. Hubbard_alpha(nt).ne.0.d0) then  
           do m1 = 1, 2 * Hubbard_l(nt) + 1  
              do m2 = m1, 2 * Hubbard_l(nt) + 1
                 do ibnd = 1, nbnd  
                    nr(na,current_spin,m1,m2) = nr(na,current_spin,m1,m2) + &
                            wg(ibnd,ik) * DREAL( proj(offset(na)+m2,ibnd) * &
                                           conjg(proj(offset(na)+m1,ibnd)) )
                 enddo
              enddo
           enddo
        endif

     enddo
     ! on k-points

  enddo
#ifdef __PARA
  call poolreduce (nat * nspin * ldim * ldim, nr)  
#endif
  !
  ! impose hermiticity of n_{m1,m2}
  !
  do na = 1, nat  
     nt = ityp(na)
     do is = 1, nspin  
        do m1 = 1, 2 * Hubbard_l(nt) + 1
           do m2 = m1 + 1, 2 * Hubbard_l(nt) + 1  
              nr (na, is, m2, m1) = nr (na, is, m1, m2)  
           enddo
        enddo
     enddo
  enddo

  ! symmetryze the quantities nr -> nsnew
  do na = 1, nat  
     nt = ityp (na)  
     if (Hubbard_U(nt).ne.0.d0 .or. Hubbard_alpha(nt).ne.0.d0) then  
        do is = 1, nspin  
           do m1 = 1, 2 * Hubbard_l(nt) + 1  
              do m2 = 1, 2 * Hubbard_l(nt) + 1  
                 do isym = 1, nsym  
                    nb = irt (isym, na)  
                    do m0 = 1, 2 * Hubbard_l(nt) + 1  
                       do m00 = 1, 2 * Hubbard_l(nt) + 1  
                          if (Hubbard_l(nt).eq.0) then
                             nsnew(na,is,m1,m2) = nsnew(na,is,m1,m2) +  &
                                   nr(nb,is,m0,m00) / nsym
                          else if (Hubbard_l(nt).eq.1) then
                             nsnew(na,is,m1,m2) = nsnew(na,is,m1,m2) +  &
                                   d1(m0 ,m1,isym) * nr(nb,is,m0,m00) * &
                                   d1(m00,m2,isym) / nsym
                          else if (Hubbard_l(nt).eq.2) then
                             nsnew(na,is,m1,m2) = nsnew(na,is,m1,m2) +  &
                                   d2(m0 ,m1,isym) * nr(nb,is,m0,m00) * &
                                   d2(m00,m2,isym) / nsym
                          else if (Hubbard_l(nt).eq.3) then
                             nsnew(na,is,m1,m2) = nsnew(na,is,m1,m2) +  &
                                   d3(m0 ,m1,isym) * nr(nb,is,m0,m00) * &
                                   d3(m00,m2,isym) / nsym
                          else
                             call errore ('new_ns', &
                                         'angular momentum not implemented', &
                                          abs(Hubbard_l(nt)) )
                          end if
                       enddo
                    enddo
                 enddo
              enddo
           enddo
        enddo
     endif
  enddo

  ! Now we make the matrix ns(m1,m2) strictly hermitean
  do na = 1, nat  
     nt = ityp (na)  
     if (Hubbard_U(nt).ne.0.d0 .or. Hubbard_alpha(nt).ne.0.d0) then  
        do is = 1, nspin  
           do m1 = 1, 2 * Hubbard_l(nt) + 1  
              do m2 = m1, 2 * Hubbard_l(nt) + 1  
                 psum = abs ( nsnew(na,is,m1,m2) - nsnew(na,is,m2,m1) )  
                 if (psum.gt.1.d-10) then  
                    write (6, * ) na, is, m1, m2  
                    write (6, * ) nsnew (na, is, m1, m2)  
                    write (6, * ) nsnew (na, is, m2, m1)  
                    call errore ('new_ns', 'non hermitean matrix', 1)  
                 else  
                    nsnew(na,is,m1,m2) = 0.5d0 * (nsnew(na,is,m1,m2) + &
                                                  nsnew(na,is,m2,m1) )
                    nsnew(na,is,m2,m1) = nsnew(na,is,m1,m2)
                 endif
              enddo
           enddo
        enddo
     endif
  enddo
  !
  ! Now the contribution to the total energy is computed. The corrections
  ! needed to obtain a variational expression are already included
  !
  eth = 0.d0  
  do na = 1, nat  
     nt = ityp (na)  
     if (Hubbard_U(nt).ne.0.d0 .or. Hubbard_alpha(nt).ne.0.d0) then  
        do is = 1, nspin  
           do m1 = 1, 2 * Hubbard_l(nt) + 1  
              do m2 = 1, 2 * Hubbard_l(nt) + 1  
                 eth = eth + Hubbard_U(nt) * nsnew(na,is,m1,m2) * &
                          (ns(na,is,m2,m1) - nsnew(na,is,m2,m1) * 0.5d0)
              enddo
           enddo
        enddo
     endif
  enddo
  deallocate ( offset, proj, nr )
  if (nspin.eq.1) eth = 2.d0 * eth

  return

end subroutine new_ns
