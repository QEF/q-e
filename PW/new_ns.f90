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
  use allocate
#ifdef PARA
  use para
#endif
  implicit none
  integer :: ik, ibnd, is, i, na, nb, nt, isym, n, counter, m1, m2, &
       m0, m00, l
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
  allocate( offset(nat), proj(natomwfc,nbnd), nr(nat,nspin,5,5) )  
  !
  ! D_Sl for l=1 and l=2 are already initialized, for l=0 D_S0 is 1
  !
  counter = 0  
  do na = 1, nat  
     nt = ityp (na)  
     do n = 1, nchi (nt)  
        if (oc (n, nt) .gt.0.d0.or..not.newpseudo (nt) ) then  
           l = lchi (n, nt)  
           if (l.eq.2) offset (na) = counter  
           counter = counter + 2 * l + 1  
        endif
     enddo

  enddo

  if (counter.ne.natomwfc) call error ('new_ns', 'nstart<>counter', 1)
  nr    (:,:,:,:) = 0.d0
  nsnew (:,:,:,:) = 0.d0
  !
  !    we start a loop on k points
  !

  if (nks.gt.1) rewind (iunigk)  

  do ik = 1, nks  
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
#ifdef PARA
     call reduce (2 * natomwfc * nbnd, proj)  
#endif
     !
     ! compute the occupation numbers (the quantities n(m1,m2)) of the
     ! atomic orbitals
     !
     do na = 1, nat  
        nt = ityp (na)  
        if (Hubbard_U(nt).ne.0.d0 .or. Hubbard_alpha(nt).ne.0.d0) then  
           do m1 = 1, 5  
              do m2 = m1, 5  
                 do ibnd = 1, nbnd  
                    nr(na,isk(ik),m1,m2) = nr(na,isk(ik),m1,m2) + &
                            wg(ibnd,ik) * DREAL( proj(offset(na)+m2,ibnd) * &
                                           conjg(proj(offset(na)+m1,ibnd)) )
                 enddo
              enddo
           enddo
        endif

     enddo
     ! on k-points

  enddo
#ifdef PARA
  call poolreduce (nat * nspin * 25, nr)  
#endif
  !
  ! impose hermiticity of n_{m1,m2}
  !
  do na = 1, nat  
     do is = 1, nspin  
        do m1 = 1, 5  
           do m2 = m1 + 1, 5  
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
           do m1 = 1, 5  
              do m2 = 1, 5  
                 do isym = 1, nsym  
                    nb = irt (isym, na)  
                    do m0 = 1, 5  
                       do m00 = 1, 5  
                          nsnew(na,is,m1,m2) = nsnew(na,is,m1,m2) +  &
                                   d2(m1,m0 ,isym) * nr(nb,is,m0,m00) * &
                                   d2(m2,m00,isym) / nsym
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
           do m1 = 1, 5  
              do m2 = m1, 5  
                 psum = abs ( nsnew(na,is,m1,m2) - nsnew(na,is,m2,m1) )  
                 if (psum.gt.1.d-10) then  
                    write (6, * ) na, is, m1, m2  
                    write (6, * ) nsnew (na, is, m1, m2)  
                    write (6, * ) nsnew (na, is, m2, m1)  
                    call error ('new_ns', 'non hermitean matrix', 1)  
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
           do m1 = 1, 5  
              do m2 = 1, 5  
                 eth = eth + Hubbard_U(nt) * nsnew(na,is,m1,m2) * &
                          (ns(na,is,m2,m1) - nsnew(na,is,m2,m1) * 0.5d0)
              enddo
           enddo
        enddo
     endif

  enddo
  deallocate ( offset, proj, nr )
  return  

end subroutine new_ns

