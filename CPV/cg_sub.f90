!
! Copyright (C) 2002 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!=======================================================================
!
   subroutine runcg_uspp( nfi, tfirst, tlast, eigr, bec, irb, eigrb, &
      rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac, fion, ema0bg, becdr, &
      lambdap, lambda  )

      use kinds, only: dp
      use control_flags, only: iprint, thdyn, tpre, tbuff, iprsta, trhor, &
            tfor, tvlocw, trhow, taurdr, tprnfor
      use control_flags, only: ndr, ndw, nbeg, nomore, tsde, tortho, tnosee, &
            tnosep, trane, tranp, tsdp, tcp, tcap, ampre, amprp, tnoseh

      use atom, only: nlcc
      use core, only: nlcc_any
!---ensemble-DFT
      use energies, only: eht, epseu, exc, etot, eself, enl, ekin,          &
     &                    atot, entropy, egrand
      use electrons_base, only: f, nspin, nel, iupdwn, nupdwn, nudx, nelt, &
                                nx => nbspx, n => nbsp, ispin => fspin

      use ensemble_dft, only: tens, tgrand, ninner, ismear, etemp, ef,       &
     &                tdynz, tdynf, zmass, fmass, fricz, fricf, z0, c0diag,  &
                      becdiag, fmat0, becdrdiag, becm, bec0, fion2, atot0,   &
                      etot0, h0c0, c0hc0, epsi0, e0, dval, z1, f1, dfmat, fmat1, &
                      ef1, enocc, f0, fmatx, fx, zaux, zx, ex, zxt, atot1, etot1, &
                      dedx1, dentdx1, dx, dadx1, faux, eqc, eqa, atotmin, xmin, &
                      entropy2, f2, etot2, eqb, compute_entropy2, compute_entropy_der, &
                      compute_entropy
!---
      use gvecp, only: ngm
      use gvecs, only: ngs
      use gvecb, only: ngb
      use gvecw, only: ngw
      use reciprocal_vectors, only: ng0 => gstart
      use cvan, only: nvb, ish
      use ions_base, only: na, nat, pmass, nax, nsp, rcmax
      use grid_dimensions, only: nnr => nnrx, nr1, nr2, nr3
      use cell_base, only: ainv, a1, a2, a3
      use cell_base, only: omega, alat
      use cell_base, only: h, hold, deth, wmass, tpiba2
      use smooth_grid_dimensions, only: nnrsx, nr1s, nr2s, nr3s
      use smallbox_grid_dimensions, only: nnrb => nnrbx, nr1b, nr2b, nr3b
      use local_pseudo, only: vps, rhops
      use io_global, ONLY: io_global_start, stdout, ionode, ionode_id
      use mp_global, ONLY: mp_global_start
      use para_mod
      use dener
      use derho
      use cdvan
      use stre
      use parameters, only: nacx, natx, nsx, nbndxx
      use constants, only: pi, factem
      use io_files, only: psfile, pseudo_dir
      use io_files, only: outdir

      use uspp, only : nhsa=> nkb, betae => vkb, rhovan => becsum, deeq
      use uspp_param, only: nh
      use cg_module, only: ltresh, itercg, etotnew, etotold, tcutoff, &
          restartcg, passof, passov, passop, ene_ok, numok, maxiter, &
          enever, etresh, ene0, hpsi, gi, hi, esse, essenew, dene0, spasso, &
          ene1, passo, iter3, enesti, ninner_ef, emme
      use ions_positions, only: tau0
      use wavefunctions_module, only: c0, cm, phi => cp
      use efield_module, only: tefield, evalue, ctable, qmat, detq, ipolp, &
            berry_energy, ctabin, gqq, gqqm, df, pberryel
      use mp, only: mp_sum,mp_bcast
!
      implicit none
!
      integer :: nfi
      logical :: tfirst , tlast
      complex(8) :: eigr(ngw,nat)
      real(8) :: bec(nhsa,n)
      real(8) :: becdr(nhsa,n,3)
      integer irb(3,nat)
      complex(8) :: eigrb(ngb,nat)
      real(8) :: rhor(nnr,nspin)
      real(8) :: rhog(ngm,nspin)
      real(8) :: rhos(nnrsx,nspin)
      real(8) :: rhoc(nnr)
      complex(8) :: ei1(-nr1:nr1,nat)
      complex(8) :: ei2(-nr2:nr2,nat)
      complex(8) :: ei3(-nr3:nr3,nat)
      complex(8) :: sfac( ngs, nsp )
      real(8) :: fion(3,natx)
      real(8) :: ema0bg(ngw)
      real(8) :: lambdap(nx,nx)
      real(8) :: lambda(nx,nx)
!
!
      integer :: i, j, ig, k, is, ia, iv, jv, il
      integer :: inl, jnl, niter, istart, nss
      real(8) :: enb, enbi, x
      complex(8) :: c2( ngw )
      complex(8) :: c3( ngw )
      real(8) :: gamma, entmp, sta
!
!


     fion2=0.d0

      if(ionode) open(37,file='convergence.dat',status='unknown')!for debug and tuning purposes
      if(tfirst.and.ionode) write(stdout,*) 'PERFORMING CONJUGATE GRADIENT MINIMIZATION OF EL. STATES'

      call prefor(eigr,betae) 

      ltresh    = .false.
      itercg    = 1
      etotold   = 1.d8
      tcutoff   = .false.
      restartcg = .true.
      passof = passop
      ene_ok = .false.


      !orthonormalize c0

      call calbec(1,nsp,eigr,c0,bec)

      call gram(betae,bec,c0)

      call calbec(1,nsp,eigr,c0,bec)

      !calculates phi for pcdaga

      call calphiid(c0,bec,betae,phi)
         
      !set index on number of converged iterations

      numok = 0

      do while ( itercg .lt. maxiter .and. (.not.ltresh) )


        ENERGY_CHECK: if(.not. ene_ok ) then
          call calbec(1,nsp,eigr,c0,bec)
          if(.not.tens) then
            call rhoofr(nfi,c0,irb,eigrb,bec,rhovan,rhor,rhog,rhos,enl,ekin)
          else
            !     calculation of the rotated quantities
            call rotate(z0,c0(:,:,1,1),bec,c0diag,becdiag)
            !     calculation of rho corresponding to the rotated wavefunctions
            call rhoofr(nfi,c0diag,irb,eigrb,becdiag                         &
                     &                    ,rhovan,rhor,rhog,rhos,enl,ekin)
          endif

          !calculates the potential
          !
          !     put core charge (if present) in rhoc(r)
          !
          if (nlcc_any) call set_cc(irb,eigrb,rhoc)

          !
          !---ensemble-DFT

          if (.not.tens) then

            call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,             &
                  &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
            etotnew=etot

          else

            call compute_entropy2( entropy, f, n, nspin )
           
            call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,                &
                     &            ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
            etotnew=etot+entropy

          end if

          if(tefield  ) then!just in this case calculates elfield stuff at zeo field-->to be bettered
            
             call berry_energy( enb, enbi, bec, c0(:,:,1,1), fion )
             etot=etot+enb+enbi
          endif
        else

          etot=enever
          if(.not.tens) then 
             etotnew=etot
          else
             etotnew=etot+entropy
          endif
          ene_ok=.false.

        end if ENERGY_CHECK
        if(ionode) write(37,*)itercg, etotnew,pberryel!for debug and tuning purposes


        

        if(abs(etotnew-etotold).lt.etresh) then
           numok=numok+1
        else 
           numok=0
        endif

        if(numok.ge.4) then
           ltresh=.true.
        endif



        etotold=etotnew
        ene0=etot

        ! cal_emme style of orhogonalization              
        !               call cal_emme(c0,bec,emme, 1)

        !update d

        call newd(rhor,irb,eigrb,rhovan,fion)


        call prefor(eigr,betae)!ATTENZIONE

        do i=1,n,2
          call dforce(bec,betae,i,c0(1,i,1,1),c0(1,i+1,1,1),c2,c3,rhos)
          if(tefield .and. (evalue.ne.0.d0)) then
            call dforceb(c0, i, betae, ipolp, bec ,ctabin(1,1,ipolp), gqq, gqqm, qmat, deeq, df)
            c2(1:ngw)=c2(1:ngw)+evalue*df(1:ngw)
            call dforceb(c0, i+1, betae, ipolp, bec ,ctabin(1,1,ipolp), gqq, gqqm, qmat, deeq, df)
            c3(1:ngw)=c3(1:ngw)+evalue*df(1:ngw)
          endif
          hpsi(1:ngw,  i)=c2(1:ngw)
          hpsi(1:ngw,i+1)=c3(1:ngw)
          if (ng0.eq.2) then
            hpsi(1,  i)=CMPLX(DBLE(hpsi(1,  i)), 0.d0)
            hpsi(1,i+1)=CMPLX(DBLE(hpsi(1,i+1)), 0.d0)
          end if
        enddo

        gi(1:ngw,1:n) = hpsi(1:ngw,1:n)
               
        ! cal_emme strategy             
        !               call pc_emmedaga(c0,phi,gi,emme)

        call pcdaga2(c0,phi,gi)

        DO i = 1, n
          gi(1:ngw,i) = gi(1:ngw,i) * ema0bg(1:ngw)
        END DO
        
        call calcmt(f,z0,fmat0)

        call calbec(1,nsp,eigr,gi,becm)

        call pcdaga2(c0,phi,gi)
        call pcdaga2(c0,phi,hpsi) 

        !cal_emme strategy
        !              call pc_emmedaga(c0,phi,gi,emme)
        !              call pc_emmedaga(c0,phi,hpsi,emme)

        !case of first iteration

        if(itercg==1.or.(mod(itercg,20).eq.1).or.restartcg) then

          restartcg=.false.
          passof=passop
          hi(1:ngw,1:n)=gi(1:ngw,1:n)

          !calculates esse for the second iteration

          gamma=0.d0
          if(.not.tens) then
            call calbec(1,nsp,eigr,gi,becm)
            do i=1,n        
              do ig=1,ngw
                gamma=gamma+2*DBLE(CONJG(gi(ig,i))*gi(ig,i))
              enddo
              if (ng0.eq.2) then
                gamma=gamma-DBLE(CONJG(gi(1,i))*gi(1,i))
              endif
            enddo
            call mp_sum(gamma)
             
 !           if (nvb.gt.0) then
 !              do is=1,nvb
 !                 do iv=1,nh(is)
 !                    do jv=1,nh(is)
 !                       do ia=1,na(is)
 !                          inl=ish(is)+(iv-1)*na(is)+ia
 !                          jnl=ish(is)+(jv-1)*na(is)+ia
 !                          !      gamma=gamma+ qq(iv,jv,is)*becm(inl,i)*becm(jnl,i)  
 !                       end do
 !                    end do
 !                 end do
 !              end do
 !           endif

          else

            do is=1,nspin
               nss=nupdwn(is)
               istart=iupdwn(is)
               do i=1,nss
                  do j=1,nss
                     do ig=1,ngw
                       gamma=gamma+2*DBLE(CONJG(gi(ig,i+istart-1))*gi(ig,j+istart-1))*fmat0(j,i,is)
                     enddo
                     if (ng0.eq.2) then
                       gamma=gamma-DBLE(CONJG(gi(1,i+istart-1))*gi(1,j+istart-1))*fmat0(j,i,is)
                     endif
                  enddo
                enddo
            enddo
            call mp_sum(gamma)

          endif

          esse=gamma

        else

          !find direction hi for general case 
          !calculates gamma for general case, not using Polak Ribiere

          gamma=0.d0

          if(.not.tens) then

            call calbec(1,nsp,eigr,gi,becm)
            do i=1,n        
               do ig=1,ngw
                  gamma=gamma+2*DBLE(CONJG(gi(ig,i))*gi(ig,i))
               enddo
               if (ng0.eq.2) then
                  gamma=gamma-DBLE(CONJG(gi(1,i))*gi(1,i))
               endif
            enddo
                
            call mp_sum(gamma)
             
!            if (nvb.gt.0) then
!             do is=1,nvb
!                do iv=1,nh(is)
!                   do jv=1,nh(is)      
!                      do ia=1,na(is)
!                         inl=ish(is)+(iv-1)*na(is)+ia
!                         jnl=ish(is)+(jv-1)*na(is)+ia
!                         ! gamma=gamma+ qq(iv,jv,is)*becm(inl,i)*becm(jnl,i)  
!                      end do
!                   end do
!                end do
!             end do
!            endif

          else
           do is=1,nspin
               nss=nupdwn(is)
               istart=iupdwn(is)
               do i=1,nss
                  do j=1,nss
                     do ig=1,ngw
                       gamma=gamma+2*DBLE(CONJG(gi(ig,i+istart-1))*gi(ig,j+istart-1))*fmat0(j,i,is)
                     enddo
                     if (ng0.eq.2) then
                       gamma=gamma-DBLE(CONJG(gi(1,i+istart-1))*gi(1,j+istart-1))*fmat0(j,i,is)
                     endif
                  enddo
                enddo
            enddo

            call mp_sum(gamma)

          endif

          essenew=gamma
          gamma=gamma/esse
          esse=essenew

          hi(1:ngw,1:n)=gi(1:ngw,1:n)+gamma*hi(1:ngw,1:n)

        endif

        !find minimum along direction hi:

        !project hi on conduction sub-space

        call calbec(1,nsp,eigr,hi,bec0)
        call pc2(c0,bec,hi,bec0)
        ! cal_emme strategy
        !     call pc_emme(c0,bec,hi,bec0,emme)
        call calbec(1,nsp,eigr,hi,bec0)

        !do quadratic minimization
        !             
        !calculate derivative with respect to  lambda along direction hi

        dene0=0.
        if(.not.tens) then
          do i=1,n               
            do ig=1,ngw
              dene0=dene0-4.d0*DBLE(CONJG(hi(ig,i))*hpsi(ig,i))
            enddo
            if (ng0.eq.2) then
              dene0=dene0+2.d0*DBLE(CONJG(hi(1,i))*hpsi(1,i))
            endif
          end do
        else
          !in the ensamble case the derivative is Sum_ij (<hi|H|Psi_j>+ <Psi_i|H|hj>)*f_ji
          !     calculation of the kinetic energy x=xmin      
          call calcmt(f,z0,fmat0)
         do is=1,nspin
             nss=nupdwn(is)
             istart=iupdwn(is)
             do i=1,nss
                do j=1,nss
                   do ig=1,ngw
                    dene0=dene0-2.d0*DBLE(CONJG(hi(ig,i+istart-1))*hpsi(ig,j+istart-1))*fmat0(j,i,is)
                    dene0=dene0-2.d0*DBLE(CONJG(hpsi(ig,i+istart-1))*hi(ig,j+istart-1))*fmat0(j,i,is)
                   enddo
                   if (ng0.eq.2) then
                    dene0=dene0+DBLE(CONJG(hi(1,i+istart-1))*hpsi(1,j+istart-1))*fmat0(j,i,is)
                    dene0=dene0+DBLE(CONJG(hpsi(1,i+istart-1))*hi(1,j+istart-1))*fmat0(j,i,is)
                   end if
                  enddo
            enddo
          enddo
        endif

        call mp_sum(dene0)

        !if the derivative is positive, search along opposite direction

        if(dene0.gt.0.d0) then
          spasso=-1.D0
        else
          spasso=1.d0
        endif

        !calculates wave-functions on a point on direction hi

        cm(1:ngw,1:n,1,1)=c0(1:ngw,1:n,1,1)+spasso*passof*hi(1:ngw,1:n)


        !orthonormalize

        call calbec(1,nsp,eigr,cm,becm)
        call gram(betae,becm,cm)
        call calbec(1,nsp,eigr,cm,becm)
               
        !calculate energy
        if(.not.tens) then
          call rhoofr(nfi,cm,irb,eigrb,becm,rhovan,rhor,rhog,rhos,enl,ekin)
        else
          !     calculation of the rotated quantities
          call rotate(z0,cm(:,:,1,1),becm,c0diag,becdiag)
          !     calculation of rho corresponding to the rotated wavefunctions
          call rhoofr(nfi,c0diag,irb,eigrb,becdiag,rhovan,rhor,rhog,rhos,enl,ekin)
        endif

        !calculate potential
        !
        !     put core charge (if present) in rhoc(r)
        !
        if (nlcc_any) call set_cc(irb,eigrb,rhoc)
!
        if (.not.tens) then
          call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,             &
                      &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
        else
          call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,             &
                      &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
        end if

        if( tefield  ) then!to be bettered
          call berry_energy( enb, enbi, becm, cm(:,:,1,1), fion )
          etot=etot+enb+enbi
        endif
        ene1=etot
              
            
        !find the minimum

        call minparabola(ene0,spasso*dene0,ene1,passof,passo,enesti)

        if(iprsta.gt.1) write(6,*) ene0,dene0,ene1,passo, gamma, esse

        !set new step

        passov=passof
        passof=2.d0*passo
              
        !calculates wave-functions at minimum

        cm(1:ngw,1:n,1,1)=c0(1:ngw,1:n,1,1)+spasso*passo*hi(1:ngw,1:n)
        if(ng0.eq.2) then
          cm(1,:,1,1)=0.5d0*(cm(1,:,1,1)+CONJG(cm(1,:,1,1)))
        endif

        call calbec(1,nsp,eigr,cm,becm)
        call gram(betae,becm,cm)

        !test on energy: check the energy has really diminished

        call calbec(1,nsp,eigr,cm,becm)
        if(.not.tens) then
          call rhoofr(nfi,cm,irb,eigrb,becm,rhovan,rhor,rhog,rhos,enl,ekin)
        else
          !     calculation of the rotated quantities
          call rotate(z0,cm(:,:,1,1),becm,c0diag,becdiag)
          !     calculation of rho corresponding to the rotated wavefunctions
          call rhoofr(nfi,c0diag,irb,eigrb,becdiag,rhovan,rhor,rhog,rhos,enl,ekin)
        endif

        !calculates the potential
        !
        !     put core charge (if present) in rhoc(r)
        !
        if (nlcc_any) call set_cc(irb,eigrb,rhoc)
!
        if (.not.tens) then
          call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,             &
                       &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
        else
          call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,             &
                       &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
        end if
        if( tefield )  then!to be bettered
          call berry_energy( enb, enbi, becm, cm(:,:,1,1), fion )
          etot=etot+enb+enbi
        endif
        enever=etot
        if(ionode) write(37,'(a3,4f20.10)') 'CG1',ene0+entropy,ene1+entropy,enesti+entropy,enever+entropy
        if(ionode) write(37,'(a3,4f10.7)')  'CG2',spasso,passov,passo,(enever-ene0)/passo/dene0

        !check with  what supposed

        if(ionode) then
            if(iprsta.gt.1) then
                 write(stdout,*) 'cg_sub: estimate :'  , (enesti-enever)/(ene0-enever)
                 write(stdout,*) 'cg_sub: minmum   :'  , enever,passo,passov
             endif
        endif

        !if the energy has diminished with respect to  ene0 and ene1 , everything ok
        if( (enever.lt.ene0) .and. (enever.lt.ene1)) then
          c0(:,:,1,1)=cm(:,:,1,1)
          bec(:,:)=becm(:,:)
          ene_ok=.true.
          !if  ene1 << energy <  ene0; go to  ene1
        else if( (enever.ge.ene0).and.(ene0.gt.ene1)) then
          if(ionode) then
             write(stdout,*) 'cg_sub: missed minimum, case 2, iteration',itercg
          endif  
          c0(1:ngw,1:n,1,1)=c0(1:ngw,1:n,1,1)+spasso*passov*hi(1:ngw,1:n)
          restartcg=.true.!ATTENZIONE
          call calbec(1,nsp,eigr,c0,bec)
          call gram(betae,bec,c0)
          !if ene > ene0,en1 do a steepest descent step
          ene_ok=.false.
        else if((enever.ge.ene0).and.(ene0.le.ene1)) then
        if(ionode) then
             write(stdout,*) 'cg_sub: missed minimum, case 3, iteration',itercg
         endif

          iter3=0
          do while(enever.gt.ene0 .and. iter3.lt.4)
            iter3=iter3+1
            passov=passov*0.5d0
            cm(1:ngw,1:n,1,1)=c0(1:ngw,1:n,1,1)+spasso*passov*hi(1:ngw,1:n)
            ! chenge the searching direction
            spasso=spasso*(-1.d0)
            call calbec(1,nsp,eigr,cm,becm)
            call gram(betae,bec,cm)
            call calbec(1,nsp,eigr,cm,becm)
            if(.not.tens) then
              call rhoofr(nfi,cm,irb,eigrb,becm,rhovan,rhor,rhog,rhos,enl,ekin)
            else
              !     calculation of the rotated quantities
              call rotate(z0,cm(:,:,1,1),becm,c0diag,becdiag)
              !     calculation of rho corresponding to the rotated wavefunctions
              call rhoofr(nfi,c0diag,irb,eigrb,becdiag,rhovan,rhor,rhog,rhos,enl,ekin)
            endif
  
            !calculates the potential
            !
            !     put core charge (if present) in rhoc(r)
            !
            if (nlcc_any) call set_cc(irb,eigrb,rhoc)
  !
            if (.not.tens) then
              call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,             &
                        &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
                    
            else
              call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,             &
                        &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
            end if
            if( tefield)  then !to be bettered
              call berry_energy( enb, enbi, becm, cm(:,:,1,1), fion )
              etot=etot+enb+enbi
            endif
            enever=etot
           
          end do
  
          c0(:,:,1,1)=cm(:,:,1,1)
          restartcg=.true.
          ene_ok=.false.
        end if
  
        call calbec (1,nsp,eigr,c0,bec)
        !calculates phi for pc_daga
        call calphiid(c0,bec,betae,phi)
  
        !=======================================================================
        !
        !                 start of the inner loop
        !                 (Uij degrees of freedom)
        !
        !=======================================================================
  
            
        if(tens) then
          if(.not. ene_ok) then
            !     calculation of the array bec:
            call calbec (1,nsp,eigr,c0,bec)
  
                   
            call rotate(z0,c0(:,:,1,1),bec,c0diag,becdiag)
  
            call rhoofr (nfi,c0diag,irb,eigrb,becdiag,rhovan,rhor,rhog,rhos,enl,ekin)
  
            !     put core charge (if present) in rhoc(r)
            if (nlcc_any) call set_cc(irb,eigrb,rhoc)
  
            !     calculation of the potential 
  
            call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,                &
                 ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion2)
            !     calculation of the array deeq: 

            call compute_entropy2( entropy, f, n, nspin )

    
            !     deeq_i,lm = \int V_eff(r) q_i,lm(r) dr
          endif
  
          ene_ok=.false.
 !     calculation of ekinc
            call calcmt(f,z0,fmat0)
          call newd(rhor,irb,eigrb,rhovan,fion2)
    
          !     free energy at x=0
          atot0=etot+entropy    
          etot0=etot      
  
          !     start of the loop
          call prefor(eigr,betae)!ATTENZIONE
  
  
          do niter=1,ninner
            h0c0 = 0.0d0



            do i=1,n,2                      
              call dforce(bec,betae,i,c0(1,i,1,1), &
              &     c0(1,i+1,1,1),h0c0(1,i),h0c0(1,i+1),rhos)
            end do
                   
            do is=1,nspin
              nss=nupdwn(is)
              istart=iupdwn(is)
              do i=1,nss
                do k=1,nss
                  c0hc0(k,i,is)=0.d0
                  do ig=1,ngw
                    c0hc0(k,i,is)=c0hc0(k,i,is)- &
                 &   2.0*DBLE(CONJG(c0(ig,k+istart-1,1,1))*h0c0(ig,i+istart-1))
                  enddo
                  if (ng0.eq.2) then
                    c0hc0(k,i,is)=c0hc0(k,i,is)+&
                 &   DBLE(CONJG(c0(1,k+istart-1,1,1))*h0c0(1,i+istart-1))
                  endif
                end do
              end do
              call mp_sum(c0hc0(1:nss,1:nss,is))
            end do
 
!           do is=1,nspin
!              nss=nupdwn(is)
!              istart=iupdwn(is)
!              do i=1,nss
!                do k=i,nss
!                     c0hc0(k,i,is)=0.5d0*(c0hc0(k,i,is)+c0hc0(i,k,is))
!                     c0hc0(i,k,is)=c0hc0(k,i,is)
!                enddo
!             enddo
!           enddo 

 
  
            do is=1,nspin
              nss=nupdwn(is)
              epsi0(1:nss,1:nss,is)=c0hc0(1:nss,1:nss,is)!ATTENZIONE
            end do
   

          
            !    diagonalization of the matrix epsi0_ij
            e0 = 0.0d0
            do  is=1,nspin
              istart=iupdwn(is)
              nss=nupdwn(is) 
              if(ionode) then
                 call ddiag(nss,nss,epsi0(1,1,is),dval(1),z1(1,1,is),1)
              endif
              call mp_bcast(dval,ionode_id)
              call mp_bcast(z1(:,:,is),ionode_id)
             do i=1,nss
                e0(i+istart-1)=dval(i)
              enddo
            enddo
  

            !     calculation of the occupations and the fermi energy
            !     corresponding to the chosen ismear,etemp and nspin
  
           call efermi(nelt,n,etemp,1,f1,ef1,e0,enocc,ismear,nspin)
  
  
            !     calculation of the initial and final occupation matrices
            !     in the z0-rotated orbital basis
              
            call calcm(f1,z1,fmat1)
             
            !     calculation of dfmat
            do is=1,nspin
              nss=nupdwn(is)
              dfmat(1:nss,1:nss,is)=-fmat0(1:nss,1:nss,is)+fmat1(1:nss,1:nss,is)
            end do
                   
            f0(1:n)=f(1:n)
                   
            !     initialization when xmin is determined by sampling 
            do il=1,1
              !     this loop is useful to check that the sampling is correct
              x=1.*DBLE(il)
              do is=1,nspin
                nss=nupdwn(is)
                fmatx(1:nss,1:nss,is)=fmat0(1:nss,1:nss,is)+x*dfmat(1:nss,1:nss,is)
              end do
                      
              !     diagonalization of fmatx
              fx = 0.0d0
              do  is=1,nspin
                nss=nupdwn(is)
                istart=iupdwn(is)
                if(ionode) then
                   call ddiag(nss,nss,fmatx(1,1,is),dval(1),zaux(1,1,is),1)
                endif
                call mp_bcast(dval,ionode_id)
                call mp_bcast(zaux(:,:,is),ionode_id)
                do i=1,nss
                  faux(i+istart-1)=dval(i)
                enddo
              enddo
              !     reordering of the eigenvalues fx and eigenvectors zx
              do  is=1,nspin
                nss=nupdwn(is)
                istart=iupdwn(is)
                do i=1,nss
                  fx(i+istart-1)=faux(nss-i+istart)
                  do j=1,nss
                    zx(i,j,is)=zaux(i,nss-j+1,is)
                  end do
                enddo
              end do
              
              !     calculation of the entropy at x
              CALL compute_entropy2( entropy, fx, n, nspin )
  
              !     calculation of the entropy derivative at x
              CALL compute_entropy_der( ex, fx, n, nspin )
  
              !    update of f
              f(1:n)=fx(1:n)
  
              !    transposition of zx (recall zx^-1=zx^t)
              !    to obtain the rotation matrix at x
              do is=1,nspin
                nss=nupdwn(is)
                do i=1,nss
                  do k=1,nss
                    zxt(k,i,is)=zx(i,k,is)
                  end do
                end do
              end do
                         
              !     calculation of the quantities at x=1
              !     using the previously calculated rotation matrix 
              !     (similar to what has been done at x=0)
              call calcmt(f,zxt,fmatx)
             
              !     calculation of the rotated quantities for the calculation
              !     of the epsi0_ij matrix at x=1
              call rotate(zxt,c0(:,:,1,1),bec,c0diag,becdiag)   
              call rhoofr (nfi,c0diag,irb,eigrb,becdiag,rhovan,rhor,rhog,rhos,enl,ekin) 
              !     put core charge (if present) in rhoc(r)
              if (nlcc_any) call set_cc(irb,eigrb,rhoc)
  
              call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,                 &
                           &            ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion2)
              call newd(rhor,irb,eigrb,rhovan,fion2)
              call prefor(eigr,betae)
    
              !     free energy at x=1
              atot1=etot+entropy
              etot1=etot
              !if(ionode) write(37,'(a5,f4.2,3f15.10)') 'x :',x,etot,entropy,atot1
           end do
  
  
            !     calculation of c0hc0_ij at x=1
            call prefor(eigr,betae)!ATTENZIONE
            h0c0 = 0.0d0
            do i=1,n,2
              call dforce(bec,betae,i,c0(1,i,1,1),c0(1,i+1,1,1),&
                h0c0(1,i),h0c0(1,i+1),rhos)
            end do
  
            do is=1,nspin
              nss=nupdwn(is)
              istart=iupdwn(is)
              do i=1,nss
                do k=1,nss
                  c0hc0(k,i,is)=0.d0
                  do ig=1,ngw
                    c0hc0(k,i,is)=c0hc0(k,i,is)-&
                    2.0*DBLE(CONJG(c0(ig,k+istart-1,1,1))*h0c0(ig,i+istart-1))
                  enddo
                  if (ng0.eq.2) then
                    c0hc0(k,i,is)=c0hc0(k,i,is)+&
                    DBLE(CONJG(c0(1,k+istart-1,1,1))*h0c0(1,i+istart-1))
                  endif
                end do
              end do
              call mp_sum(c0hc0(1:nss,1:nss,is))
            end do
  
            do is=1,nspin
              nss=nupdwn(is)
              epsi0(1:nss,1:nss,is)=c0hc0(1:nss,1:nss,is)
            end do
                   
            !     calculation of da/dx(x=1)=da/df_ji*delta(f_ji) 
            !     recall: dtrS(f)/df_ij = S'(f)_ji
            !     The calculation of the d(-TS)/dx is done using 
            !       (ex)_j [(zt)_ji (dfmat)_ik (zt)_jk] 
            !     instead of [(zt)_jk (ex)_j (zt)_ji] (dfmat)_ik 
            !     because of the quantity (ex)_j is not well conditioned
                      
            dedx1=0.0
            dentdx1=0.0
            
            do is=1,nspin
              nss=nupdwn(is)
              istart=iupdwn(is)
              do i=1,nss
                 dx(i+istart-1)=0.
                 do k=1,nss
                    do j=1,nss
                       dx(i+istart-1)=dx(i+istart-1)-                              &
                            &         zxt(i,k,is)*fmat0(k,j,is)*zxt(i,j,is)
                    end do
                 end do
                 dx(i+istart-1)=dx(i+istart-1)+fx(i+istart-1)
              end do
            end do
  
            do is=1,nspin
              nss=nupdwn(is)
              istart=iupdwn(is)
              do i=1,nss
                 dentdx1=dentdx1-etemp*dx(i+istart-1)*ex(i+istart-1)
                 do k=1,nss
                    dedx1=dedx1+dfmat(i,k,is)*epsi0(k,i,is)
                 end do
              end do
            end do
            dadx1=dedx1+dentdx1
                   
            !     line minimization (using a second degree interpolation)
            !     (fermi-dirac distribution)
            !     The free energy curve is approximated as follows
            !     (a) E+Ek -> 2nd order polynomial P(x)=eqa*x**2+eqb*x+eqc
            !         such that P(1)=E+EK(1), P(0)=E+EK(0), P'(1)=(E+EK)'(1)
            !     (b) S -> S~(x)=\sum_i s(a_i*x**2+b_i*x+c_i)
            !         (where s(f)=-f*ln(f)-(1-f)*ln(1-f))
            !         such that S~(1)=S(1), S~(0)=S(0), S~'(1)=S'(1) 
            !         => S~=\sum_i s(f1_i+(x-1)*d_i+(-f1_i+d_i+f0_i)*(x-1)**2
            !            where d_i=zxt_i,k*dfmat_k,j*zxt_i,j
  
            eqc=etot0
            eqa=dedx1-etot1+etot0
            eqb=etot1-etot0-eqa
  
            !     sampling along the search direction to find xmin
            atotmin=atot0
            xmin=0.0
                   
            do il=0,2000
                      
              x=0.0005*DBLE(il)     
              entropy2=0.0
                      
              do is=1,nspin
                nss=nupdwn(is)
                istart=iupdwn(is)
                do i=1,nss
                  f2=fx(i+istart-1)+(x-1)*dx(i+istart-1) + &
                     (-fx(i+istart-1)+dx(i+istart-1)+f0(i+istart-1))*(x-1)**2
                  CALL compute_entropy( entmp, f2, nspin )
                  entropy2 = entropy2 + entmp
                end do
              end do
              etot2=eqa*x**2+eqb*x+eqc
!              write(*,*) 'Energie2: ',x,etot2+entropy2,entropy2
              if ((etot2+entropy2).lt.atotmin) then
                xmin=x
                atotmin=etot2+entropy2
              endif
  
            end do
 
  
            !                     eqc=atot0
            !                     eqa=dadx1-atot1+atot0
            !                     eqb=atot1-atot0-eqa
  
            !                     xmin=-eqb/(2.d0*eqa)
  
  
            !     if xmin=1, no need do recalculate the quantities
            if(ionode) write(37,'(a5,2f15.10)') 'XMIN', xmin,atotmin 
            if (xmin.eq.1.) goto 300
  
            !     calculation of the fmat at x=xmin
            !     this part can be optimized in the case where xmin=0
            do is=1,nspin
              nss=nupdwn(is)
              do i=1,nss
                do j=1,nss
                  fmatx(i,j,is)=fmat0(i,j,is)+xmin*dfmat(i,j,is)
                end do
              end do
            enddo
               
            !     diagonalization of fmat at x=xmin
            fx = 0.0d0
            do is=1,nspin
               nss=nupdwn(is)
               istart=iupdwn(is)
               if(ionode) call ddiag(nss,nss,fmatx(1,1,is),dval(1),zaux(1,1,is),1)  
               call mp_bcast(dval,ionode_id)
               call mp_bcast(zaux(:,:,is),ionode_id)
               do i=1,n
                  faux(i+istart-1)=dval(i)
               enddo
            enddo
            do  is=1,nspin
               nss=nupdwn(is)
               istart=iupdwn(is)
               do i=1,nss
                  fx(i+istart-1)=faux(nss-i+istart)
                  do j=1,nss
                     zx(i,j,is)=zaux(i,nss-j+1,is)
                  end do
               enddo
            end do
    
            !     calculation of the entropy at x=xmin
                CALL compute_entropy2( entropy, fx, n, nspin )
  
            !     update of f
            f(1:n)=fx(1:n)
  
            !     update of z0  
 300        continue
            do is=1,nspin
              nss=nupdwn(is)
              do i=1,nss
                do k=1,nss
                  z0(k,i,is)=zx(i,k,is)
                end do
              end do
            end do
  
            !     calculation of the rotated quantities





            call calcmt(f,z0,fmat0)
            call rotate(z0,c0(:,:,1,1),bec,c0diag,becdiag)

            call rhoofr (nfi,c0diag,irb,eigrb,becdiag,rhovan,rhor,rhog,rhos,enl,ekin) 
  
            !     put core charge (if present) in rhoc(r)
            if (nlcc_any) call set_cc(irb,eigrb,rhoc)
  
            call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,                &
                    &            ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion2)
            call newd(rhor,irb,eigrb,rhovan,fion2)
            CALL compute_entropy2( entropy, f, n, nspin )
            call prefor(eigr,betae)
            ene_ok=.true. !so does not calculate the energy again
    
            !     free energy at x=xmin
            atotmin=etot+entropy      
            if(ionode) write(37,'(a3,i2,2f15.10)') 'CI',niter,atot0,atotmin
  
  
            atot0=atotmin
            etot0=etot
            enever=etot
! set atot
            atot=atot0
            !     end of the loop
  
          end do!su ninnner
      
          !=======================================================================
          !                 end of the inner loop
          !=======================================================================
  
        endif !su tens
  
        itercg=itercg+1
  
      end do!on conjugate gradient iterations
      !calculates atomic forces and lambda
      call newd(rhor,irb,eigrb,rhovan,fion)
      if (.not.tens) then
        if (tfor .or. tprnfor) call nlfq(c0,eigr,bec,becdr,fion)
      else
        if (tfor .or. tprnfor) call nlfq(c0diag,eigr,becdiag,becdrdiag,fion)
      endif
  
        call prefor(eigr,betae)
        do i=1,n,2
          call dforce(bec,betae,i,c0(1,i,1,1),c0(1,i+1,1,1),c2,c3,rhos)
          if(tefield.and.(evalue .ne. 0.d0)) then
            call dforceb &
               (c0, i, betae, ipolp, bec ,ctabin(1,1,ipolp), gqq, gqqm, qmat, deeq, df)
            do ig=1,ngw
              c2(ig)=c2(ig)+evalue*df(ig)
            enddo
            call dforceb &
               (c0, i+1, betae, ipolp, bec ,ctabin(1,1,ipolp), gqq, gqqm, qmat, deeq, df)
            do ig=1,ngw
              c3(ig)=c3(ig)+evalue*df(ig)
            enddo
          endif
          do ig=1,ngw
            gi(ig,  i)=c2(ig)
            gi(ig,i+1)=c3(ig)
          end do
          if (ng0.eq.2) then
            gi(1,  i)=CMPLX(DBLE(gi(1,  i)),0.d0)
            gi(1,i+1)=CMPLX(DBLE(gi(1,i+1)),0.d0)
          end if
        enddo
        do i=1,n
          do j=i,n
            lambda(i,j)=0.d0
            do ig=1,ngw
              lambda(i,j)=lambda(i,j)-2.d0*DBLE(CONJG(c0(ig,i,1,1))*gi(ig,j))
            enddo
            if(ng0.eq.2) then
              lambda(i,j)=lambda(i,j)+DBLE(CONJG(c0(1,i,1,1))*gi(1,j))
            endif
            lambda(j,i)=lambda(i,j)
          enddo
        enddo
  
        call mp_sum(lambda)
  
        if(tens) then!in the ensemble case matrix labda must be multiplied with f
          do i=1,n
            do j=1,n
              lambdap(i,j)=0.d0
              do k=1,n
                lambdap(i,j)=lambdap(i,j)+lambda(i,k)*fmat0(k,j,1)
              end do
            end do
          enddo
          do i=1,n
            do j=1,n
              sta=lambda(i,j)
              lambda(i,j)=lambdap(i,j)
              lambdap(i,j)=sta
            enddo
          enddo
        call nlsm2(ngw,nhsa,n,eigr,c0(:,:,1,1),becdr,.true.)
        endif
        call nlfl(bec,becdr,lambda,fion)
          
        ! bforceion adds the force term due to electronic berry phase
        ! only in US-case
          
        if( tefield.and.(evalue .ne. 0.d0) ) then
           call bforceion(fion,tfor.or.tprnfor,ipolp, qmat,bec,becdr,gqq,evalue)

        endif

      if(ionode) close(37)!for debug and tuning purposes
END SUBROUTINE
