!
! Copyright (C) 2005 Paolo Umari
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
subroutine h_epsi_her_apply(lda, n,nbande, psi, hpsi, pdir, e_field)
  !-----------------------------------------------------------------------
  !
  ! this subroutine applies w_k+w_k* on psi, 
  ! (as in Souza et al.  PRB B 69, 085106 (2004))
  ! the output is put into hpsi
  !
  ! evcel must contain the wavefunctions from previous iteration
  ! spin polarized systems supported only with fixed occupations

  USE noncollin_module,     ONLY : noncolin, npol
  USE kinds,    ONLY : DP
  USE spin_orb, ONLY: lspinorb
  USE us
  USE wvfct,    ONLY : npwx, nbnd, ik => current_k
  USE ldaU,     ONLY : lda_plus_u
  USE lsda_mod, ONLY : current_spin, nspin
  USE scf,      ONLY : vrs  
  USE gvect
  USE uspp
  USE uspp_param, ONLY: nh, nhm, nbetam
  USE bp
  USE klist
  USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
  USE ions_base, ONLY: ityp, tau, nat,ntyp => nsp
  USE constants, ONLY : e2, pi, tpi, fpi
  USE fixed_occ
  USE io_global, ONLY : stdout
  USE becmod,    ONLY : calbec,bec_type,allocate_bec_type,deallocate_bec_type
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,        ONLY : mp_sum
  !
  implicit none
  INTEGER, INTENT(in) :: pdir!direction on which the polarization is calculated
  REAL(DP) :: e_field!electric field along pdir    

  !
  INTEGER :: lda !leading dimension
  INTEGER ::  n! total number of wavefunctions 
  INTEGER :: nbande!number of wavefunctions to be calculated 
  
  COMPLEX(DP) :: psi (lda*npol, nbande ), hpsi (lda*npol,nbande)


  COMPLEX(DP), EXTERNAL :: zdotc
  
  COMPLEX(DP), ALLOCATABLE  :: evct(:,:)!temporary array
  COMPLEX(DP) :: ps(nkb,nbnd*npol)


  TYPE(bec_type) :: becp0
  INTEGER :: nkbtona(nkb)
   INTEGER :: nkbtonh(nkb)
  COMPLEX(DP) :: sca, sca1, pref
  INTEGER :: npw, nb,mb, jkb, nhjkb, na, np, nhjkbm,jkb1,i,j,iv
  INTEGER :: jkb_bp,nt,ig, ijkb0,ibnd,jh,ih,ikb
  REAL(dp) :: eps
  COMPLEX(kind=DP), ALLOCATABLE :: sca_mat(:,:),sca_mat1(:,:)
  COMPLEX(kind=DP) :: pref0(4)

  !  --- Define a small number ---
  eps=0.000001d0
  if(ABS(e_field)<eps) return
  call start_clock('h_epsi_apply')

  ALLOCATE( evct(npwx*npol,nbnd))
  call allocate_bec_type(nkb,nbnd,becp0)
  npw = ngk(ik) 
  if(okvan) then
!  --- Initialize arrays ---
     jkb_bp=0
      DO nt=1,ntyp
         DO na=1,nat
            IF (ityp(na)== nt) THEN
               DO i=1, nh(nt)
                  jkb_bp=jkb_bp+1
                  nkbtona(jkb_bp) = na
                  nkbtonh(jkb_bp) = i
               END DO
            END IF
         END DO
      END DO
      CALL calbec ( npw, vkb, psi, becp0, nbande )
  endif

  allocate(sca_mat(nbnd,nbande),sca_mat1(nbnd,nbande))
  call ZGEMM('C','N',nbnd,nbande,npw,(1.d0,0.d0),evcel,npwx*npol,psi,npwx*npol,(0.d0,0.d0),sca_mat,nbnd)
  if(noncolin) then
     call ZGEMM('C','N',nbnd,nbande,npw,(1.d0,0.d0),evcel(npwx+1,1),npwx*npol,psi(npwx+1,1),npwx*npol,(1.d0,0.d0),sca_mat,nbnd)
  endif
  call mp_sum( sca_mat, intra_bgrp_comm )
  if(okvan) then
     call start_clock('h_eps_van2')
      
!apply w_k 

      
        do nb=1,nbande
           DO jkb=1,nkb
              nhjkb = nkbtonh(jkb)
              na = nkbtona(jkb)
              np = ityp(na)
              nhjkbm = nh(np)
              jkb1 = jkb - nhjkb
              pref0=(0.d0,0.d0)
              DO j = 1,nhjkbm
                 ! bec_evcel is relative to ik
                 if(lspinorb) then
                    pref0(1) = pref0(1)+becp0%nc(jkb1+j,1,nb) &
                         *qq_so(nhjkb,j,1,np)
                    pref0(2) = pref0(2)+becp0%nc(jkb1+j,2,nb) &
                         *qq_so(nhjkb,j,2,np)
                    pref0(3) = pref0(3)+becp0%nc(jkb1+j,1,nb) &
                         *qq_so(nhjkb,j,3,np)
                    pref0(4) = pref0(4)+becp0%nc(jkb1+j,2,nb) &
                         *qq_so(nhjkb,j,4,np)

                 else
                    pref0(1) = pref0(1)+becp0%k(jkb1+j,nb) &
                         *qq(nhjkb,j,np)
                 endif
              END DO
              DO mb=1,nbnd
                 if(lspinorb) then
                    pref=(0.d0,0.d0)
                    pref = pref+CONJG(bec_evcel%nc(jkb,1,mb))*pref0(1)
                    pref = pref+CONJG(bec_evcel%nc(jkb,1,mb))*pref0(2)
                    pref = pref+CONJG(bec_evcel%nc(jkb,2,mb))*pref0(3)
                    pref = pref+CONJG(bec_evcel%nc(jkb,2,mb))*pref0(4)

                 else
                    pref = CONJG(bec_evcel%k(jkb,mb))*pref0(1)
                 endif
                 sca_mat(mb,nb)=sca_mat(mb,nb)+pref
              END DO


           ENDDO
        ENDDO

     call stop_clock('h_eps_van2')
  end if
  call ZGEMM('N','N',npw,nbande,nbnd,fact_hepsi(ik,pdir),evcelm(1,1,pdir),npwx*npol,&
       &sca_mat,nbnd,(1.d0,0.d0),hpsi,npwx*npol)
  call ZGEMM('N','N',npw,nbande,nbnd,-fact_hepsi(ik,pdir),evcelp(1,1,pdir),npwx*npol,&
       &sca_mat,nbnd,(1.d0,0.d0),hpsi,npwx*npol)
  if (noncolin) then
     call ZGEMM('N','N',npw,nbande,nbnd,fact_hepsi(ik,pdir),evcelm(1+npwx,1,pdir),npwx*npol,&
          &sca_mat,nbnd,(1.d0,0.d0),hpsi(1+npwx,1),npwx*npol)
     call ZGEMM('N','N',npw,nbande,nbnd,-fact_hepsi(ik,pdir),evcelp(1+npwx,1,pdir),npwx*npol,&
          &sca_mat,nbnd,(1.d0,0.d0),hpsi(1+npwx,1),npwx*npol)
  endif

!apply w_k*
  if(.not.okvan) then
     do nb=1,nbande
        do mb=1,nbnd!index on states of evcel        
           sca = zdotc(npw,evcelm(1,mb,pdir),1,psi(1,nb),1)
           IF (noncolin) sca = sca+zdotc(npw,evcelm(1+npwx,mb,pdir),1,psi(1+npwx,nb),1)
           sca1 = zdotc(npw,evcelp(1,mb,pdir),1,psi(1,nb),1)
           IF (noncolin) sca1 = sca1 + zdotc(npw,evcelp(1+npwx,mb,pdir),1,psi(1+npwx,nb),1)
           call mp_sum( sca, intra_bgrp_comm )
           call mp_sum(  sca1, intra_bgrp_comm )
           
           do ig=1,npw

              hpsi(ig,nb) = hpsi(ig,nb) + &
                   &     CONJG(fact_hepsi(ik,pdir))*evcel(ig,mb)*(sca-sca1)
              IF (noncolin) hpsi(ig+npwx,nb) = hpsi(ig+npwx,nb) + &
                   & CONJG(fact_hepsi(ik,pdir))*evcel(ig+npwx,mb)*(sca-sca1)
           enddo
        enddo
     end do
  else ! US case
    
          
     call start_clock('h_eps_ap_van')
! copy evcel into evct
     do iv=1,nbnd
        do ig=1,npwx*npol
           evct(ig,iv)=evcel(ig,iv)
        enddo
     enddo
!  calculate S|evct>
     call start_clock('h_eps_van2')
     ps (:,:) = (0.d0, 0.d0)
     ijkb0 = 0
     do nt = 1, ntyp
        do na = 1, nat
           if (ityp (na) == nt) then
              do ibnd = 1, nbnd
                 do jh = 1, nh (nt)
                    jkb = ijkb0 + jh
                    do ih = 1, nh (nt)
                       ikb = ijkb0 + ih
                       if(lspinorb) then
                          ps (ikb, (ibnd-1)*npol+1) = ps (ikb, (ibnd-1)*npol+1) + &
                               qq_so(ih,jh,1,nt)* bec_evcel%nc(jkb,1,ibnd)
                          ps (ikb, (ibnd-1)*npol+1) = ps (ikb, (ibnd-1)*npol+1) + &
                               qq_so(ih,jh,2,nt)* bec_evcel%nc(jkb,2,ibnd)
                          ps (ikb, (ibnd-1)*npol+2) = ps (ikb, (ibnd-1)*npol+2) + &
                               qq_so(ih,jh,3,nt)* bec_evcel%nc(jkb,1,ibnd)
                          ps (ikb, (ibnd-1)*npol+2) = ps (ikb, (ibnd-1)*npol+2) + &
                               qq_so(ih,jh,4,nt)* bec_evcel%nc(jkb,2,ibnd)
                          

                       else
                          ps (ikb, ibnd) = ps (ikb, ibnd) + &
                               qq(ih,jh,nt)* bec_evcel%k(jkb,ibnd)
                       endif
                    enddo
                 enddo
              enddo
              ijkb0 = ijkb0 + nh (nt)
           endif
        enddo
     enddo
     call stop_clock('h_eps_van2')
     call ZGEMM ('N', 'N', npw, nbnd*npol , nkb, (1.d0, 0.d0) , vkb, &!vkb is relative to the last ik read
          npwx, ps, nkb, (1.d0, 0.d0) , evct, npwx)
!!!
     call ZGEMM('C','N',nbnd,nbande,npw,(1.d0,0.d0),evcelm(1,1,pdir),npwx*npol,psi,npwx*npol,(0.d0,0.d0),sca_mat,nbnd)
     if(noncolin) then
        call ZGEMM('C','N',nbnd,nbande,npw,(1.d0,0.d0),evcelm(npwx+1,1,pdir),npwx*npol,&
             &psi(npwx+1,1),npwx*npol,(1.d0,0.d0),sca_mat,nbnd)
     endif
     call mp_sum( sca_mat, intra_bgrp_comm )
     call ZGEMM('C','N',nbnd,nbande,npw,(1.d0,0.d0),evcelp(1,1,pdir),npwx*npol,psi,npwx*npol,(0.d0,0.d0),sca_mat1,nbnd)
     if(noncolin) then
        call ZGEMM('C','N',nbnd,nbande,npw,(1.d0,0.d0),evcelp(npwx+1,1,pdir),npwx*npol,&
             &psi(npwx+1,1),npwx*npol,(1.d0,0.d0),sca_mat1,nbnd)
     endif
     call mp_sum( sca_mat1, intra_bgrp_comm )
!!!!!

     sca_mat(1:nbnd,1:nbande)=sca_mat(1:nbnd,1:nbande)-sca_mat1(1:nbnd,1:nbande)
     call ZGEMM('N','N',npw,nbande,nbnd,dconjg(fact_hepsi(ik,pdir)),evct(1,1),npwx*npol,&
          &sca_mat,nbnd,(1.d0,0.d0),hpsi,npwx*npol)
     if (noncolin) then
        call ZGEMM('N','N',npw,nbande,nbnd,dconjg(fact_hepsi(ik,pdir)),evct(1+npwx,1),npwx*npol,&
             &sca_mat,nbnd,(1.d0,0.d0),hpsi(1+npwx,1),npwx*npol)
     endif

     call stop_clock('h_eps_ap_van')
  END if

  DEALLOCATE( evct)

  call deallocate_bec_type(becp0)

  call stop_clock('h_epsi_apply')

  deallocate(sca_mat)
  deallocate(sca_mat1)
!  --
!------------------------------------------------------------------------------!
   return
 END SUBROUTINE h_epsi_her_apply
