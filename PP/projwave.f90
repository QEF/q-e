!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine do_projwfc (nodenumber)
  !-----------------------------------------------------------------------
  !
  ! projects wavefunctions onto orthogonalized atomic wavefunctions
  ! calculates Lowdin charges, spilling parameter
  ! input: namelist "&inputpp", with variables
  ! prefix      prefix of input files saved by program pwscf
  ! tmp_dir     temporary directory where files resides
  ! filproj     output file containing the results
  !
  use pwcom
  use io

  implicit none
  character (len=3)  :: nodenumber
  character (len=14) :: filproj
  integer :: ios
  namelist / inputpp / tmp_dir, prefix, filproj
  !
  nd_nmbr = nodenumber
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  tmp_dir = './'
  filproj = ' '
  !
  read (5, inputpp, err = 200, iostat = ios)
200 call error ('projwave', 'reading inputpp namelist', abs (ios) )
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  call read_file
  call openfil
  !
  call projwave (filproj)
  !
  return
end subroutine do_projwfc

!-----------------------------------------------------------------------
subroutine projwave (filproj)
  !-----------------------------------------------------------------------
  !
#include "machine.h"
  use pwcom
  use becmod
  use io
#ifdef PARA
  use para
#endif
  implicit none
  character (len=14) :: filproj
  !
  type wfc_label
     integer na, n, l, m
  end type wfc_label
  type(wfc_label), allocatable :: nlmchi(:)
  !
  integer :: ik, ibnd, i, j, k, na, nb, nt, isym, n,  m, m1, l, lm, nwfc,&
       nwfc1, lmax_wfc
  logical :: exst
  real(kind=DP) :: psum, totcharge
  real(kind=DP), allocatable :: e (:), proj (:,:,:), charges(:,:)
  complex(kind=DP), allocatable :: wfcatom (:,:),  overlap (:,:), &
       work (:,:), work1(:), proj0(:,:)
  integer, allocatable :: index(:)
  !
  !
  if (filproj.eq.' ') return
  write (6, '(/5x,"Calling projwave .... ", &
       &            /5x,"Projections are written on file ",a)') filproj

  !
  allocate(swfcatom (npwx , natomwfc ) )
  allocate(wfcatom (npwx, natomwfc) )
  allocate(proj (natomwfc, nbnd, nkstot) )
  allocate(overlap (natomwfc, natomwfc) )
  allocate(e (natomwfc) )

  proj   = 0.d0
  overlap= (0.d0,0.d0)
  !
  ! initialize D_Sl for l=1 and l=2, for l=0 D_S0 is 1
  !
  call d_matrix (d1, d2)
  !
  ! fill structure nlmchi
  !
  allocate (nlmchi(natomwfc))
  nwfc=0
  lmax_wfc = 0
  do na = 1, nat
     nt = ityp (na)
     do n = 1, nchi (nt)
        if (oc (n, nt) .gt.0.d0.or..not.newpseudo (nt) ) then
           l = lchi (n, nt)
           lmax_wfc = max (lmax_wfc, l )
           do m = 1, 2 * l + 1
              nwfc=nwfc+1
              nlmchi(nwfc)%na = na
              nlmchi(nwfc)%n  =  n
              nlmchi(nwfc)%l  =  l
              nlmchi(nwfc)%m  =  m
           enddo
        endif
     enddo
  enddo

  if (lmax_wfc.gt.2) call error ('projwave', 'l > 2 not yet implemented', 1)
  if (nwfc.ne.natomwfc) call error ('projwave', 'wrong # of atomic wfcs?', 1)
  !
  !    loop on k points
  !
  call init_us_1
  !
  do ik = 1, nks
     npw = npwx
     call gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
     call davcio (evc, nwordwfc, iunwfc, ik, - 1)

     call atomic_wfc (ik, wfcatom)
     call init_us_2 (npw, igk, xk (1, ik), vkb)

     call ccalbec (nkb, npwx, npw, natomwfc, becp, vkb, wfcatom)

     call s_psi (npwx, npw, natomwfc, wfcatom, swfcatom)
     !
     ! wfcatom = |phi_i> , swfcatom = \hat S |phi_i>
     ! calculate overlap matrix O_ij = <phi_i|\hat S|\phi_j>
     !
     call ZGEMM ('c', 'n', natomwfc, natomwfc, npw, (1.d0, 0.d0) , &
          wfcatom, npwx, swfcatom, npwx, (0.d0, 0.d0) , overlap, natomwfc)
#ifdef PARA
     call reduce (2 * natomwfc * natomwfc, overlap)
#endif
     !
     ! calculate O^{-1/2}
     !
     allocate(work (natomwfc, natomwfc) )
     call cdiagh (natomwfc, overlap, natomwfc, e, work)
     do i = 1, natomwfc
        e (i) = 1.d0 / dsqrt (e (i) )
     enddo
     do i = 1, natomwfc
        do j = i, natomwfc
           overlap (i, j) = (0.d0, 0.d0)
           do k = 1, natomwfc
              overlap (i, j) = overlap (i, j) + e (k) * work (j, k) * conjg (work (i, k) )
           enddo
           if (j.ne.i) overlap (j, i) = conjg (overlap (i, j))
        enddo
     enddo
     deallocate (work)
     !
     ! calculate wfcatom = O^{-1/2} \hat S | phi>
     !
     call ZGEMM ('n', 't', npw, natomwfc, natomwfc, (1.d0, 0.d0) , &
          swfcatom, npwx,  overlap, natomwfc, (0.d0, 0.d0), wfcatom, npwx)
     !
     ! make the projection <psi_i| O^{-1/2} \hat S | phi_j>
     !
     allocate(proj0(natomwfc,nbnd) )
     call ZGEMM ('c', 'n', natomwfc, nbnd, npw, (1.d0, 0.d0) , &
          wfcatom, npwx, evc, npwx, (0.d0, 0.d0) , proj0, natomwfc)
#ifdef PARA
     call reduce (2 * natomwfc * nbnd, proj0)
#endif
     !
     ! symmetrize the projections
     !
     allocate(work1 (nbnd) )
     do nwfc = 1, natomwfc
        !
        !  atomic wavefunction nwfc is on atom na
        !
        na= nlmchi(nwfc)%na
        n = nlmchi(nwfc)%n
        l = nlmchi(nwfc)%l
        m = nlmchi(nwfc)%m
        !
        do isym = 1, nsym
           nb = irt (isym, na)
           do nwfc1 =1, natomwfc
              if (nlmchi(nwfc1)%na.eq. nb             .and. &
                   nlmchi(nwfc1)%n .eq. nlmchi(nwfc)%n .and. &
                   nlmchi(nwfc1)%l .eq. nlmchi(nwfc)%l .and. &
                   nlmchi(nwfc1)%m .eq. 1 ) go to 10
           end do
           call error('projwave','cannot symmetrize',1)
10         nwfc1=nwfc1-1
           !
           !  nwfc1 is the first rotated atomic wfc corresponding to nwfc
           !
           if (l.eq.0) then
              work1(:) = proj0 (nwfc1 + 1,:)
           else if (l.eq.1) then
              work1(:) = 0.d0
              do m1 = 1, 3
                 work1(:) = work1(:) + d1 (m, m1, isym) * &
                      proj0 (nwfc1 + m1,:)
              enddo
           else if (l.eq.2) then
              work1(:) = 0.d0
              do m1 = 1, 5
                 work1(:) = work1(:) + d2 (m, m1, isym) * &
                      proj0 (nwfc1 + m1,:)
              enddo
           endif
           do ibnd = 1, nbnd
              proj (nwfc, ibnd, ik) = proj (nwfc, ibnd, ik) + &
                   work1(ibnd) * conjg (work1(ibnd)) / nsym
           enddo
        enddo
     enddo
     deallocate (work1)
     deallocate (proj0 )
     ! on k-points
  enddo
#ifdef PARA
  !
  !   recover the vector proj over the pools
  !
  call poolrecover (et, nbndx, nkstot, nks)
  call poolrecover (proj, nbnd * natomwfc, nkstot, nks)
  !
  if (me.eq.1.and.mypool.eq.1) then
#endif
     !
     ! write on the output file
     !
     call seqopn (4, filproj, 'formatted', exst)
     write(4,'(/"Projection on atomic states:"/)')
     do nwfc = 1, natomwfc
        write(4,'(5x,"state #",i3,": atom ",i3," (",a3,"), wfc ",i2, &
             &       " (l=",i1," m=",i2,")")') &
             nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), &
             nlmchi(nwfc)%n, nlmchi(nwfc)%l, nlmchi(nwfc)%m
     end do
     !
     allocate(index (natomwfc) )
     do ik = 1, nkstot
        write (4, '(/" k = ",3f14.10)') (xk (i, ik) , i = 1, 3)
        do ibnd = 1, nbnd
           write (4, '(5x,"e = ",f14.10," eV")') et (ibnd, ik) * rytoev
           !
           ! sort projections by magnitude, in decreasing order
           !
           do nwfc = 1, natomwfc
              index (nwfc) = 0
              e (nwfc) = - proj (nwfc, ibnd, ik)
           end do
           call hpsort (natomwfc, e, index)
           !
           !  only projections that are larger than 0.001 are written
           !
           do nwfc = 1, natomwfc
              e (nwfc) = - e(nwfc)
              if ( abs (e(nwfc)).lt.0.001 ) go to 20
           end do
           nwfc = natomwfc + 1
20         nwfc = nwfc -1
           !
           ! fancy (?!?) formatting
           !
           write (4, '(5x,"psi = ",5(f5.3,"*[#",i3,"]+"))') &
               (e (i), index(i), i = 1, min(5,nwfc))
           do j = 1, (nwfc-1)/5
              write (4, '(10x,"+",5(f5.3,"*[#",i3,"]+"))') &
                (e (i), index(i), i = 5*j+1, min(5*(j+1),nwfc))
           end do
           psum = 0.d0
           do nwfc = 1, natomwfc
              psum = psum + proj (nwfc, ibnd, ik)
           end do
           write (4, '(4x,"|psi|^2 = ",f5.3)') psum
           !
        enddo
     enddo
     deallocate (index)
     !
     ! estimate partial charges (Loewdin) on each atom
     !
     allocate ( charges (nat, 0:lmax_wfc ) )
     charges = 0.0
     do ik = 1, nkstot
        do ibnd = 1, nbnd
           do nwfc = 1, natomwfc
              na= nlmchi(nwfc)%na
              l = nlmchi(nwfc)%l
              charges(na,l) = charges(na,l) + wg (ibnd,ik) * &
                   proj (nwfc, ibnd, ik)
           enddo
        end do
     end do
     !
     write (4, '(/"Lowdin Charges: "/)')
     !
     psum = 0.0
     do na = 1, nat
        totcharge = 0.d0
        do l = 0, lmax_wfc
           totcharge = totcharge + charges(na,l)
        end do
        psum = psum + totcharge
        write (4, '(5x,"Atom # ",i3,": total charge = ",f8.4, &
      &                ", s, p, d = ",4f8.4    )') &
             na, totcharge, ( charges(na,l), l= 0,lmax_wfc)
     end do
     psum = psum / nelec
     write (4, '(5x,"Spilling Parameter: ",f8.4)') 1.0 - psum
     !
     ! Sanchez-Portal et al., Sol. State Commun.  95, 685 (1995).
     ! The spilling parameter measures the ability of the basis provided by
     ! the pseudo-atomic wfc to represent the PW eigenstates,
     ! by measuring how much of the subspace of the Hamiltonian
     ! eigenstates falls outside the subspace spanned by the atomic basis
     !
     close (unit=4)
     deallocate (charges)
#ifdef PARA
  endif
#endif
  deallocate (nlmchi)
  deallocate (e)
  deallocate (overlap)
  deallocate (proj)
  deallocate (wfcatom)
  deallocate (swfcatom)

  return

end subroutine projwave
