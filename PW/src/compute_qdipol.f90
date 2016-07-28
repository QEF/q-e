!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE compute_qdipol(dpqq)
  !
  ! This routine computes the term dpqq, i.e. the dipole moment of the
  ! augmentation charge. The output is given on cartesian coordinates
  !
  USE kinds, only: DP
  USE constants, ONLY: fpi
  USE atom, ONLY: rgrid
  USE ions_base, ONLY: ntyp => nsp
  USE uspp, only: nhtol, nhtolm, indv, nlx, ap
  USE uspp_param, only: upf, nbetam, nh, nhm

  implicit none

  REAL(DP) :: dpqq( nhm, nhm, 3, ntyp)
  real(DP), allocatable :: qrad2(:,:,:), qtot(:,:,:), aux(:)
  real(DP) :: fact
  integer :: nt, l, ir, nb, mb, ijv, ilast, ipol, ih, ivl, jh, jvl, lp, ndm

  call start_clock('cmpt_qdipol')
  ndm = MAXVAL ( upf(1:ntyp)%kkbeta )
  allocate (qrad2( nbetam , nbetam, ntyp))    
  allocate (aux( ndm))    
  allocate (qtot( ndm, nbetam, nbetam))    

  qrad2(:,:,:)=0.d0
  dpqq=0.d0

  do nt = 1, ntyp
     if ( upf(nt)%tvanp ) then
        l=1
!
!   Only l=1 terms enter in the dipole of Q
!
        do nb = 1, upf(nt)%nbeta
           do mb = nb, upf(nt)%nbeta
              ijv = mb * (mb-1) /2 + nb
              if ( ( l >= abs(upf(nt)%lll(nb) - upf(nt)%lll(mb)) ) .and. &
                   ( l <=     upf(nt)%lll(nb) + upf(nt)%lll(mb)  ) .and. &
                   (mod (l+upf(nt)%lll(nb)+upf(nt)%lll(mb), 2) == 0) ) then
                 qtot(1:upf(nt)%kkbeta,nb,mb) = upf(nt)%qfuncl(1:upf(nt)%kkbeta,ijv,l)
              endif
           enddo
        enddo
        do nb=1, upf(nt)%nbeta
           !
           !    the Q are symmetric with respect to indices
           !
           do mb=nb, upf(nt)%nbeta
              if ( ( l >= abs(upf(nt)%lll(nb) - upf(nt)%lll(mb)) ) .and. &
                   ( l <=     upf(nt)%lll(nb) + upf(nt)%lll(mb)  ) .and. &
                   (mod (l+upf(nt)%lll(nb)+upf(nt)%lll(mb), 2) == 0) ) then
                 do ir = 1, upf(nt)%kkbeta
                    aux(ir)=rgrid(nt)%r(ir)*qtot(ir, nb, mb)
                 enddo
                 call simpson ( upf(nt)%kkbeta, aux, rgrid(nt)%rab, &
                                qrad2(nb,mb,nt) )
              endif
           enddo
        enddo
     endif
     ! ntyp
  enddo
  
  do ipol = 1,3
     fact=-sqrt(fpi/3.d0)
     if (ipol.eq.1) lp=3
     if (ipol.eq.2) lp=4
     if (ipol.eq.3) then
        lp=2
        fact=-fact
     endif
     do nt = 1,ntyp
        if ( upf(nt)%tvanp ) then
           do ih = 1, nh(nt)
              ivl = nhtolm(ih, nt)
              mb = indv(ih, nt)
              do jh = ih, nh (nt)
                 jvl = nhtolm(jh, nt)
                 nb=indv(jh,nt)
                 if (ivl > nlx) call errore('compute_qdipol',' ivl > nlx', ivl)
                 if (jvl > nlx) call errore('compute_qdipol',' jvl > nlx', jvl)
                 if (nb > nbetam) &
                    call errore('compute_qdipol',' nb out of bounds', nb)
                 if (mb > nbetam) &
                    call errore('compute_qdipol',' mb out of bounds', mb)
                 if (mb > nb) call errore('compute_qdipol',' mb > nb', 1)
                 dpqq(ih,jh,ipol,nt)=fact*ap(lp,ivl,jvl)*qrad2(mb,nb,nt)
                 dpqq(jh,ih,ipol,nt)=dpqq(ih,jh,ipol,nt)
                 ! WRITE( stdout,'(3i5,2f15.9)') ih,jh,ipol,dpqq(ih,jh,ipol,nt)
              enddo
           enddo
        endif
     enddo
  enddo
  deallocate(qtot)
  deallocate(aux)
  deallocate(qrad2)
  call stop_clock('cmpt_qdipol')

  return
end subroutine compute_qdipol
