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
      use control_flags, only: iprint, thdyn, tpre, iprsta, &
            tfor, tvlocw, taurdr, tprnfor
      use control_flags, only: ndr, ndw, nbeg, nomore, tsde, tortho, tnosee, &
            tnosep, trane, tranp, tsdp, tcp, tcap, ampre, amprp, tnoseh

      use atom, only: nlcc
      use core, only: nlcc_any
!---ensemble-DFT
      use energies, only: eht, epseu, exc, etot, eself, enl, ekin,          &
     &                    atot, entropy, egrand
      use electrons_base, only: f, nspin, nel, iupdwn, nupdwn, nudx, nelt, &
                                nx => nbspx, n => nbsp, ispin

      use ensemble_dft, only: tens,   ef,  z0, c0diag,  &
                      becdiag, fmat0, e0, ismear
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
      use io_global,                ONLY : io_global_start, stdout, ionode, ionode_id
      use mp_global,                ONLY : intra_image_comm, np_ortho, me_ortho, ortho_comm
      use dener
      use derho
      use cdvan
      use stre
      use constants,                only : pi, au_gpa
      use io_files,                 only : psfile, pseudo_dir
      use io_files,                 only : outdir
      use uspp,                     only : nhsa=> nkb, nhsavb=> nkbus, betae => vkb, rhovan => becsum, deeq,qq
      use uspp_param,               only : nh
      use cg_module,                only : ltresh, itercg, etotnew, etotold, tcutoff, &
                                           restartcg, passof, passov, passop, ene_ok, numok, maxiter, &
                                           enever, conv_thr, ene0, esse, essenew, dene0, spasso, &
                                           ene1, passo, iter3, enesti, ninner_ef
      use ions_positions,           only : tau0
      use wavefunctions_module,     only : c0, cm, phi => cp
      use efield_module,            only : tefield, evalue, ctable, qmat, detq, ipolp, &
                                           berry_energy, ctabin, gqq, gqqm, df, pberryel, &
                                           tefield2, evalue2, ctable2, qmat2, detq2, ipolp2, &
                                           berry_energy2, ctabin2, gqq2, gqqm2, pberryel2
      use mp,                       only : mp_sum
      use cp_electronic_mass,       ONLY : emass_cutoff
      use orthogonalize_base,       ONLY : calphi
      use cp_interfaces,            ONLY : rhoofr, dforce
      USE cpr_subroutines,          ONLY : compute_stress
      USE printout_base,            ONLY : printout_stress
      USE cp_main_variables,        ONLY : nlax, collect_lambda, distribute_lambda, descla


!
      implicit none
!
      integer :: nfi
      logical :: tfirst , tlast
      complex(dp) :: eigr(ngw,nat)
      real(dp) :: bec(nhsa,n)
      real(dp) :: becdr(nhsa,n,3)
      integer irb(3,nat)
      complex(dp) :: eigrb(ngb,nat)
      real(dp) :: rhor(nnr,nspin)
      complex(dp) :: rhog(ngm,nspin)
      real(dp) :: rhos(nnrsx,nspin)
      real(dp) :: rhoc(nnr)
      complex(dp) :: ei1(-nr1:nr1,nat)
      complex(dp) :: ei2(-nr2:nr2,nat)
      complex(dp) :: ei3(-nr3:nr3,nat)
      complex(dp) :: sfac( ngs, nsp )
      real(dp) :: fion(3,nat)
      real(dp) :: ema0bg(ngw)
      real(dp) :: lambdap(nlax,nlax,nspin)
      real(dp) :: lambda(nlax,nlax,nspin)
!
!
      integer :: i, j, ig, k, is, iss,ia, iv, jv, il, ii, jj, kk
      integer :: inl, jnl, niter, istart, nss
      real(dp) :: enb, enbi, x
      complex(dp) :: c2( ngw )
      complex(dp) :: c3( ngw )
      real(dp) :: gamma, entmp, sta
      complex(dp),allocatable :: hpsi(:,:), hpsi0(:,:), gi(:,:), hi(:,:)
      real(DP), allocatable::               s_minus1(:,:)!factors for inverting US S matrix
      real(DP), allocatable::               k_minus1(:,:)!factors for inverting US preconditioning matrix 
      real(DP), allocatable :: lambda_repl(:,:) ! replicated copy of lambda
      real(DP), allocatable :: lambda_dist(:,:) ! replicated copy of lambda
      real(dp) :: sca, dumm(1)
      logical  :: newscheme, firstiter
      integer :: maxiter3
!
!
      real(kind=DP), allocatable :: bec0(:,:), becm(:,:), becdrdiag(:,:,:)



      allocate(bec0(nhsa,n),becm(nhsa,n), becdrdiag(nhsa,n,3))


      call start_clock('runcg_uspp')
      newscheme=.false.
      firstiter=.true.

      maxiter3=250


      if(ionode) open(37,file='convergence.dat',status='unknown')!for debug and tuning purposes
      if(tfirst.and.ionode) write(stdout,*) 'PERFORMING CONJUGATE GRADIENT MINIMIZATION OF EL. STATES'
      
!set tpa preconditioning

      call  emass_precond_tpa( ema0bg, tpiba2, emass_cutoff )
     
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

      call gram(betae,bec,nhsa,c0,ngw,n)

      !call calbec(1,nsp,eigr,c0,bec)

      !calculates phi for pcdaga

      ! call calphiid(c0,bec,betae,phi)
      CALL calphi( c0, SIZE(c0,1), bec, nhsa, betae, phi, n )

      !calculates the factors for S and K inversion in US case
      if(nvb.gt.0) then
         allocate( s_minus1(nhsavb,nhsavb))
         allocate( k_minus1(nhsavb,nhsavb))
        call  set_x_minus1(betae,s_minus1,dumm,.false.)
        call  set_x_minus1(betae,k_minus1,ema0bg,.true.)
      else
         allocate( s_minus1(1,1))
         allocate( k_minus1(1,1))
      endif  


!-------------verifica
!      do i=1,n
!      do ig=1,ngw
!         phi(ig,i,1,1)=phi(ig,i,1,1)+c0(ig,i,1,1)*(1.d0/ema0bg(ig)-1.d0)
!      enddo
!      enddo
!      !call calbec(1,nsp,eigr,phi,becm)
!      !call sminus1(phi,becm,betae)
!      call kminus1(phi,betae,ema0bg)
!      call calbec(1,nsp,eigr,phi,becm)
!      do i=1,n
!         do j=1,n
!             sca=0.d0
!                do ig=1,ngw
!                 sca=sca+2*DBLE(CONJG(phi(ig,i,1,1))*phi(ig,j,1,1))
!              enddo
!              if (ng0.eq.2) then
!                 sca=sca-DBLE(CONJG(phi(1,i,1,1))*phi(1,j,1,1))
!              endif


!           if (nvb.gt.0) then
!              do is=1,nvb
!                 do iv=1,nh(is)
!                    do jv=1,nh(is)
!                       do ia=1,na(is)
!                          inl=ish(is)+(iv-1)*na(is)+ia
!                          jnl=ish(is)+(jv-1)*na(is)+ia
!                          sca=sca+ qq(iv,jv,is)*becm(inl,i)*becm(jnl,j)
!                       end do
!                    end do
!                 end do
!              end do
!           endif
!            write(6,*) 'VERIFCA S :',i,j,sca
!         enddo
!       enddo






 
      !set index on number of converged iterations

      numok = 0

      allocate(hpsi(ngw,n),hpsi0(ngw,n),gi(ngw,n),hi(ngw,n))
      do while ( itercg .lt. maxiter .and. (.not.ltresh) )


        ENERGY_CHECK: if(.not. ene_ok ) then
          call calbec(1,nsp,eigr,c0,bec)
          if(.not.tens) then
            call rhoofr(nfi,c0(:,:),irb,eigrb,bec,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
          else

           if(newscheme.or.firstiter) then 
              if(ismear==2) then !fermi dirac smearing
                 call  inner_loop( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                      rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac,c0,bec,firstiter)
              else ! generalized smearings
                 call  inner_loop_cold( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                      rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac,c0,bec,firstiter)
              endif
               firstiter=.false.
           endif
            !     calculation of the rotated quantities
            call rotate(z0,c0(:,:),bec,c0diag,becdiag,.false.)
            !     calculation of rho corresponding to the rotated wavefunctions
            call rhoofr(nfi,c0diag,irb,eigrb,becdiag                         &
                     &                    ,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
          endif

          !calculates the potential
          !
          !     put core charge (if present) in rhoc(r)
          !
          if (nlcc_any) call set_cc(irb,eigrb,rhoc)

          !
          !---ensemble-DFT

          call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,             &
                 &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
          if (.not.tens) then
            etotnew=etot
          else
            etotnew=etot+entropy
          end if

          if(tefield  ) then!just in this case calculates elfield stuff at zeo field-->to be bettered
            
             call berry_energy( enb, enbi, bec, c0(:,:), fion )
             etot=etot+enb+enbi
          endif
          if(tefield2  ) then!just in this case calculates elfield stuff at zeo field-->to be bettered

             call berry_energy2( enb, enbi, bec, c0(:,:), fion )
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
        if(ionode) write(37,*)itercg, etotnew,pberryel,pberryel2!for debug and tuning purposes


        

        if(abs(etotnew-etotold).lt.conv_thr) then
           numok=numok+1
        else 
           numok=0
        endif

        if(numok.ge.4) then
           ltresh=.true.
        endif



        etotold=etotnew
        ene0=etot
        if(tens .and. newscheme) then
          ene0=ene0+entropy
        endif


        !update d

        call newd(rhor,irb,eigrb,rhovan,fion)


        call prefor(eigr,betae)!ATTENZIONE

        do i=1,n,2
          call dforce( i, bec, betae, c0,c2,c3,rhos, nnrsx, ispin,f,n,nspin)
          if(tefield .and. (evalue.ne.0.d0)) then
            call dforceb(c0, i, betae, ipolp, bec ,ctabin(1,1,ipolp), gqq, gqqm, qmat, deeq, df)
            c2(1:ngw)=c2(1:ngw)+evalue*df(1:ngw)
            call dforceb(c0, i+1, betae, ipolp, bec ,ctabin(1,1,ipolp), gqq, gqqm, qmat, deeq, df)
            c3(1:ngw)=c3(1:ngw)+evalue*df(1:ngw)
          endif
          if(tefield2 .and. (evalue2.ne.0.d0)) then
            call dforceb(c0, i, betae, ipolp2, bec ,ctabin2(1,1,ipolp2), gqq2, gqqm2, qmat2, deeq, df)
            c2(1:ngw)=c2(1:ngw)+evalue2*df(1:ngw)
            call dforceb(c0, i+1, betae, ipolp2, bec ,ctabin2(1,1,ipolp2), gqq2, gqqm2, qmat2, deeq, df)
            c3(1:ngw)=c3(1:ngw)+evalue2*df(1:ngw)
          endif

          hpsi(1:ngw,  i)=c2(1:ngw)
          if(i+1 <= n) then
            hpsi(1:ngw,i+1)=c3(1:ngw)
          endif
          if (ng0.eq.2) then
            hpsi(1,  i)=CMPLX(DBLE(hpsi(1,  i)), 0.d0)
            if(i+1 <= n) then
              hpsi(1,i+1)=CMPLX(DBLE(hpsi(1,i+1)), 0.d0)
            endif
          end if
        enddo

               
        call pcdaga2(c0,phi,hpsi)
               
        hpsi0(1:ngw,1:n)=hpsi(1:ngw,1:n)
        gi(1:ngw,1:n) = hpsi(1:ngw,1:n)
        
        call calbec(1,nsp,eigr,hpsi,becm)
        call xminus1(hpsi,betae,dumm,becm,s_minus1,.false.)
!        call sminus1(hpsi,becm,betae)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!look if the following two lines are really needed
        call calbec(1,nsp,eigr,hpsi,becm)
        call pc2(c0,bec,hpsi,becm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        call kminus1(gi,betae,ema0bg)
        call xminus1(gi,betae,ema0bg,becm,k_minus1,.true.)
        call calbec(1,nsp,eigr,gi,becm)
        call pc2(c0,bec,gi,becm)

        
        if(tens) call calcmt(f,z0,fmat0,.false.)

        call calbec(1,nsp,eigr,hpsi,bec0) 

!  calculates gamma
        gamma=0.d0
        
        if(.not.tens) then
           
           do i=1,n
              do ig=1,ngw
                 gamma=gamma+2*DBLE(CONJG(gi(ig,i))*hpsi(ig,i))
              enddo
              if (ng0.eq.2) then
                 gamma=gamma-DBLE(CONJG(gi(1,i))*hpsi(1,i))
              endif
           enddo
           
           call mp_sum( gamma, intra_image_comm )
           
           if (nvb.gt.0) then
              do is=1,nvb
                 do iv=1,nh(is)
                    do jv=1,nh(is)
                       do ia=1,na(is)
                          inl=ish(is)+(iv-1)*na(is)+ia
                          jnl=ish(is)+(jv-1)*na(is)+ia
                          gamma=gamma+ qq(iv,jv,is)*becm(inl,i)*bec0(jnl,i)
                       end do
                    end do
                 end do
              end do
           endif

        else

           do iss=1,nspin
              nss=nupdwn(iss)
              istart=iupdwn(iss)
              do i=1,nss
                 do j=1,nss
                    do ig=1,ngw
                       gamma=gamma+2*DBLE(CONJG(gi(ig,i+istart-1))*hpsi(ig,j+istart-1))*fmat0(j,i,iss)
                    enddo
                    if (ng0.eq.2) then
                       gamma=gamma-DBLE(CONJG(gi(1,i+istart-1))*hpsi(1,j+istart-1))*fmat0(j,i,iss)
                    endif
                 enddo
              enddo
           enddo
           call mp_sum( gamma, intra_image_comm )
           if(nvb.gt.0) then
              do iss=1,nspin
                 nss=nupdwn(iss)
                 istart=iupdwn(iss)
                 do i=1,nss
                    do j=1,nss
                       do is=1,nvb
                          do iv=1,nh(is)
                             do jv=1,nh(is)
                                do ia=1,na(is)
                                   inl=ish(is)+(iv-1)*na(is)+ia
                                   jnl=ish(is)+(jv-1)*na(is)+ia
                                   gamma=gamma+ qq(iv,jv,is)*becm(inl,i+istart-1)*bec0(jnl,j+istart-1)*fmat0(j,i,iss)
                                end do
                             end do
                          end do
                       enddo
                    enddo
                 enddo
              enddo
           endif
        endif


        !case of first iteration

        if(itercg==1.or.(mod(itercg,20).eq.1).or.restartcg) then

          restartcg=.false.
          passof=passop
          hi(1:ngw,1:n)=gi(1:ngw,1:n)!hi is the search direction
          esse=gamma


        else

          !find direction hi for general case 
          !calculates gamma for general case, not using Polak Ribiere
          
          essenew=gamma
          gamma=gamma/esse
          esse=essenew

          hi(1:ngw,1:n)=gi(1:ngw,1:n)+gamma*hi(1:ngw,1:n)

        endif
!note that hi, is saved  on gi, because we need it before projection on conduction states

        !find minimum along direction hi:

        !project hi on conduction sub-space

        call calbec(1,nsp,eigr,hi,bec0)
        call pc2(c0,bec,hi,bec0)
        

        !do quadratic minimization
        !             
        !calculate derivative with respect to  lambda along direction hi

        dene0=0.
        if(.not.tens) then
          do i=1,n               
            do ig=1,ngw
              dene0=dene0-4.d0*DBLE(CONJG(hi(ig,i))*hpsi0(ig,i))
            enddo
            if (ng0.eq.2) then
              dene0=dene0+2.d0*DBLE(CONJG(hi(1,i))*hpsi0(1,i))
            endif
          end do
        else
          !in the ensamble case the derivative is Sum_ij (<hi|H|Psi_j>+ <Psi_i|H|hj>)*f_ji
          !     calculation of the kinetic energy x=xmin      
         call calcmt(f,z0,fmat0,.false.)
         do is=1,nspin
            nss=nupdwn(is)
            istart=iupdwn(is)
            do i=1,nss
               do j=1,nss
                  do ig=1,ngw
                     dene0=dene0-2.d0*DBLE(CONJG(hi(ig,i+istart-1))*hpsi0(ig,j+istart-1))*fmat0(j,i,is)
                     dene0=dene0-2.d0*DBLE(CONJG(hpsi0(ig,i+istart-1))*hi(ig,j+istart-1))*fmat0(j,i,is)
                  enddo
                  if (ng0.eq.2) then
                     dene0=dene0+DBLE(CONJG(hi(1,i+istart-1))*hpsi0(1,j+istart-1))*fmat0(j,i,is)
                     dene0=dene0+DBLE(CONJG(hpsi0(1,i+istart-1))*hi(1,j+istart-1))*fmat0(j,i,is)
                  end if
               enddo
            enddo
         enddo
      endif

      call mp_sum( dene0, intra_image_comm )

        !if the derivative is positive, search along opposite direction
      if(dene0.gt.0.d0) then
         spasso=-1.D0
      else
         spasso=1.d0
      endif

        !calculates wave-functions on a point on direction hi

      cm(1:ngw,1:n)=c0(1:ngw,1:n)+spasso*passof*hi(1:ngw,1:n)


        !orthonormalize

      call calbec(1,nsp,eigr,cm,becm)
      call gram(betae,becm,nhsa,cm,ngw,n)
        !call calbec(1,nsp,eigr,cm,becm)
               
        !calculate energy
        if(.not.tens) then
          call rhoofr(nfi,cm(:,:),irb,eigrb,becm,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
        else
          if(newscheme) then 
             if(ismear==2) then
                call  inner_loop( nfi, tfirst, tlast, eigr,  irb, eigrb, &
           rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac,cm,becm,.false.  )
             else
                call  inner_loop_cold( nfi, tfirst, tlast, eigr,  irb, eigrb, &
           rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac,cm,becm,.false.  )  
             endif
          endif

          !     calculation of the rotated quantities
          call rotate(z0,cm(:,:),becm,c0diag,becdiag,.false.)
          !     calculation of rho corresponding to the rotated wavefunctions
          call rhoofr(nfi,c0diag,irb,eigrb,becdiag,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
        endif

        !calculate potential
        !
        !     put core charge (if present) in rhoc(r)
        !
        if (nlcc_any) call set_cc(irb,eigrb,rhoc)
!
        call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,             &
                      &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)

        if( tefield  ) then!to be bettered
          call berry_energy( enb, enbi, becm, cm(:,:), fion )
          etot=etot+enb+enbi
        endif
        if( tefield2  ) then!to be bettered
          call berry_energy2( enb, enbi, becm, cm(:,:), fion )
          etot=etot+enb+enbi
        endif

        ene1=etot
        if(tens.and.newscheme) then
          ene1=ene1+entropy
        endif
              
            
        !find the minimum

        call minparabola(ene0,spasso*dene0,ene1,passof,passo,enesti)

        if(iprsta.gt.1) write(6,*) ene0,dene0,ene1,passo, gamma, esse

        !set new step

        passov=passof
        passof=2.d0*passo
              
        !calculates wave-functions at minimum

        cm(1:ngw,1:n)=c0(1:ngw,1:n)+spasso*passo*hi(1:ngw,1:n)
        if(ng0.eq.2) then
          cm(1,:)=0.5d0*(cm(1,:)+CONJG(cm(1,:)))
        endif

        call calbec(1,nsp,eigr,cm,becm)
        call gram(betae,becm,nhsa,cm,ngw,n)

        !test on energy: check the energy has really diminished

        !call calbec(1,nsp,eigr,cm,becm)
        if(.not.tens) then
          call rhoofr(nfi,cm(:,:),irb,eigrb,becm,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
        else
          if(newscheme)  then
             if(ismear==2) then
                call  inner_loop( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                     rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac,cm,becm,.false.  )
             else
                call  inner_loop_cold( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                      rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac,cm,becm,.false.  )
             endif
          endif
          !     calculation of the rotated quantities
          call rotate(z0,cm(:,:),becm,c0diag,becdiag,.false.)
          !     calculation of rho corresponding to the rotated wavefunctions
          call rhoofr(nfi,c0diag,irb,eigrb,becdiag,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
        endif

        !calculates the potential
        !
        !     put core charge (if present) in rhoc(r)
        !
        if (nlcc_any) call set_cc(irb,eigrb,rhoc)
!
        call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,             &
                       &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)

        if( tefield )  then!to be bettered
          call berry_energy( enb, enbi, becm, cm(:,:), fion )
          etot=etot+enb+enbi
        endif
        if( tefield2 )  then!to be bettered
          call berry_energy2( enb, enbi, becm, cm(:,:), fion )
          etot=etot+enb+enbi
        endif

        enever=etot
        if(tens.and.newscheme) then
          enever=enever+entropy
        endif
        if(tens.and.newscheme) then
          if(ionode) write(37,'(a3,4f20.10)') 'CG1',ene0,ene1,enesti,enever
          if(ionode) write(37,'(a3,4f10.7)')  'CG2',spasso,passov,passo,(enever-ene0)/passo/dene0
        else
          if(ionode) write(37,'(a3,4f20.10)') 'CG1',ene0+entropy,ene1+entropy,enesti+entropy,enever+entropy
          if(ionode) write(37,'(a3,4f10.7)')  'CG2',spasso,passov,passo,(enever-ene0)/passo/dene0
        endif
        !check with  what supposed

        if(ionode) then
            if(iprsta.gt.1) then
                 write(stdout,*) 'cg_sub: estimate :'  , (enesti-enever)/(ene0-enever)
                 write(stdout,*) 'cg_sub: minmum   :'  , enever,passo,passov
             endif
        endif

        !if the energy has diminished with respect to  ene0 and ene1 , everything ok
        if( ((enever.lt.ene0) .and. (enever.lt.ene1)).or.(tefield.or.tefield2)) then
          c0(:,:)=cm(:,:)
          bec(:,:)=becm(:,:)
          ene_ok=.true.
        elseif( (enever.ge.ene1) .and. (enever.lt.ene0)) then
          if(ionode) then
             write(stdout,*) 'cg_sub: missed minimum, case 1, iteration',itercg
          endif
          c0(1:ngw,1:n)=c0(1:ngw,1:n)+spasso*passov*hi(1:ngw,1:n)
          restartcg=.true.
          call calbec(1,nsp,eigr,c0,bec)
          call gram(betae,bec,nhsa,c0,ngw,n)
          ene_ok=.false.
          !if  ene1 << energy <  ene0; go to  ene1
        else if( (enever.ge.ene0).and.(ene0.gt.ene1)) then
          if(ionode) then
             write(stdout,*) 'cg_sub: missed minimum, case 2, iteration',itercg
          endif  
          c0(1:ngw,1:n)=c0(1:ngw,1:n)+spasso*passov*hi(1:ngw,1:n)
          restartcg=.true.!ATTENZIONE
          call calbec(1,nsp,eigr,c0,bec)
          call gram(betae,bec,nhsa,c0,ngw,n)
          !if ene > ene0,en1 do a steepest descent step
          ene_ok=.false.
        else if((enever.ge.ene0).and.(ene0.le.ene1)) then
        if(ionode) then
             write(stdout,*) 'cg_sub: missed minimum, case 3, iteration',itercg
         endif

          iter3=0
          do while(enever.gt.ene0 .and. iter3.lt.maxiter3)
            iter3=iter3+1
            passov=passov*0.5d0
            cm(1:ngw,1:n)=c0(1:ngw,1:n)+spasso*passov*hi(1:ngw,1:n)
            ! chenge the searching direction
            spasso=spasso*(-1.d0)
            call calbec(1,nsp,eigr,cm,becm)
            call gram(betae,bec,nhsa,cm,ngw,n)
            call calbec(1,nsp,eigr,cm,becm)
            if(.not.tens) then
              call rhoofr(nfi,cm(:,:),irb,eigrb,becm,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
            else
              if(newscheme)  then
                 if(ismear==2) then
                    call  inner_loop( nfi, tfirst, tlast, eigr,  irb, eigrb, &
              rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac,cm,becm,.false.  )
                 else
                    call  inner_loop_cold( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                         rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac,cm,becm,.false.  )
                 endif
              endif
              !     calculation of the rotated quantities
              call rotate(z0,cm(:,:),becm,c0diag,becdiag,.false.)
              !     calculation of rho corresponding to the rotated wavefunctions
              call rhoofr(nfi,c0diag,irb,eigrb,becdiag,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
            endif
  
            !calculates the potential
            !
            !     put core charge (if present) in rhoc(r)
            !
            if (nlcc_any) call set_cc(irb,eigrb,rhoc)
  !
            call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,             &
                        &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)

            if( tefield)  then !to be bettered
              call berry_energy( enb, enbi, becm, cm(:,:), fion )
              etot=etot+enb+enbi
            endif
            if( tefield2)  then !to be bettered
              call berry_energy2( enb, enbi, becm, cm(:,:), fion )
              etot=etot+enb+enbi
            endif

            enever=etot
           if(tens.and.newscheme) then
             enever=enever+entropy
           endif

          end do
          if(iter3 == maxiter3) write(stdout,*) 'missed minimun: iter3 = maxiter3'
          c0(:,:)=cm(:,:)
          restartcg=.true.
          ene_ok=.false.
        end if
        
        if(tens.and.newscheme) enever=enever-entropy
 
        if(.not. ene_ok) call calbec (1,nsp,eigr,c0,bec)

        !calculates phi for pc_daga
        CALL calphi( c0, SIZE(c0,1), bec, nhsa, betae, phi, n )
  
        !=======================================================================
        !
        !                 start of the inner loop
        !                 (Uij degrees of freedom)
        !
        !=======================================================================
        if(tens.and. .not.newscheme) then
           if(ismear==2) then!fermi dirac
              call  inner_loop( nfi, tfirst, tlast, eigr,  irb, eigrb, &
           rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac,c0,bec,firstiter  )
           else!cold smearing and all others
               call  inner_loop_cold( nfi, tfirst, tlast, eigr,  irb, eigrb, &
           rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac,c0,bec,firstiter  )
            endif
!the following sets up the new energy
            enever=etot
         endif
      
          !=======================================================================
          !                 end of the inner loop
          !=======================================================================
  
  
        itercg=itercg+1

!   restore hi
!        hi(:,:)=gi(:,:) 


      end do!on conjugate gradient iterations
      !calculates atomic forces and lambda

       if(tpre) then!if pressure is need the following is written because of caldbec
          call  calbec(1,nsp,eigr,c0,bec)
          if(.not.tens) then
            call  caldbec( ngw, nhsa, n, 1, nsp, eigr, c0, dbec, .true. )
            call rhoofr(nfi,c0(:,:),irb,eigrb,bec,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
          else

            !     calculation of the rotated quantities
            call rotate(z0,c0(:,:),bec,c0diag,becdiag,.false.)
            !     calculation of rho corresponding to the rotated wavefunctions
            call  caldbec( ngw, nhsa, n, 1, nsp, eigr, c0diag, dbec, .true. )
            call rhoofr(nfi,c0diag,irb,eigrb,becdiag                         &
                     &                    ,rhovan,rhor,rhog,rhos,enl,denl,ekin,dekin6)
          endif

          !calculates the potential
          !
          !     put core charge (if present) in rhoc(r)
          !
          if (nlcc_any) call set_cc(irb,eigrb,rhoc)

          !
          !---ensemble-DFT

          call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,             &
                 &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)

   

     endif


     call calcmt(f,z0,fmat0,.false.)

      call newd(rhor,irb,eigrb,rhovan,fion)
      if (.not.tens) then
        if (tfor .or. tprnfor) call nlfq(c0,eigr,bec,becdr,fion)
      else
        if (tfor .or. tprnfor) call nlfq(c0diag,eigr,becdiag,becdrdiag,fion)
      endif
  
        call prefor(eigr,betae)
        do i=1,n,2
          call dforce(i,bec,betae,c0,c2,c3,rhos,nnrsx,ispin,f,n,nspin)
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
          if(tefield2.and.(evalue2 .ne. 0.d0)) then
            call dforceb &
               (c0, i, betae, ipolp2, bec ,ctabin2(1,1,ipolp2), gqq2, gqqm2, qmat2, deeq, df)
            do ig=1,ngw
              c2(ig)=c2(ig)+evalue2*df(ig)
            enddo
            call dforceb &
               (c0, i+1, betae, ipolp2, bec ,ctabin2(1,1,ipolp2), gqq2, gqqm2, qmat2, deeq, df)
            do ig=1,ngw
              c3(ig)=c3(ig)+evalue2*df(ig)
            enddo
          endif

          do ig=1,ngw
            gi(ig,  i)=c2(ig)
            if(i+1 <= n) then
              gi(ig,i+1)=c3(ig)
            endif
          end do
          if (ng0.eq.2) then
            gi(1,  i)=CMPLX(DBLE(gi(1,  i)),0.d0)
            if(i+1 <= n) then
              gi(1,i+1)=CMPLX(DBLE(gi(1,i+1)),0.d0)
            endif
          end if
        enddo

        ALLOCATE( lambda_repl( nudx, nudx ) )
        !
        do is = 1, nspin
           !
           nss = nupdwn(is)
           istart = iupdwn(is)
           lambda_repl = 0.d0
           !
           !
           do i = 1, nss
              do j = i, nss
                 ii = i + istart - 1
                 jj = j + istart - 1
                 do ig = 1, ngw
                    lambda_repl( i, j ) = lambda_repl( i, j ) - &
                       2.d0 * DBLE( CONJG( c0( ig, ii ) ) * gi( ig, jj) )
                 enddo
                 if( ng0 == 2 ) then
                    lambda_repl( i, j ) = lambda_repl( i, j ) + &
                       DBLE( CONJG( c0( 1, ii ) ) * gi( 1, jj ) )
                 endif
                 lambda_repl( j, i ) = lambda_repl( i, j )
              enddo
           enddo
           !
           CALL mp_sum( lambda_repl, intra_image_comm )
           !
           CALL distribute_lambda( lambda_repl, lambda( :, :, is ), descla( :, is ) )
           !
        end do

        DEALLOCATE( lambda_repl )
  
        if ( tens ) then
           !
           ! in the ensemble case matrix labda must be multiplied with f

           ALLOCATE( lambda_dist( nlax, nlax ) )
 
           do iss = 1, nspin
              !
              nss = nupdwn(iss)
              istart = iupdwn(iss)
              !
              lambdap(:,:,iss) = 0.0d0
              !
              CALL distribute_lambda( fmat0(:,:,iss), lambda_dist, descla( :, iss ) )
              !
              ! Perform lambdap = lambda * fmat0
              !
              CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, lambda(1,1,iss), nlax, lambda_dist, nlax, &
                                  0.0d0, lambdap(1,1,iss), nlax, descla(1,iss) )
              !
              lambda_dist      = lambda(:,:,iss)
              lambda(:,:,iss)  = lambdap(:,:,iss)
              lambdap(:,:,iss) = lambda_dist
              !
           end do
           !
           DEALLOCATE( lambda_dist )
           !
           call nlsm2(ngw,nhsa,n,eigr,c0(:,:),becdr,.true.)
           !
        endif
        !
  
        !
        call nlfl(bec,becdr,lambda,fion)
          
        ! bforceion adds the force term due to electronic berry phase
        ! only in US-case
          
        if( tefield.and.(evalue .ne. 0.d0) ) then
           call bforceion(fion,tfor.or.tprnfor,ipolp, qmat,bec,becdr,gqq,evalue)

        endif
        if( tefield2.and.(evalue2 .ne. 0.d0) ) then
           call bforceion(fion,tfor.or.tprnfor,ipolp2, qmat2,bec,becdr,gqq2,evalue2)
        endif

        deallocate(hpsi0,hpsi,gi,hi)
        deallocate( s_minus1,k_minus1)
       if(ionode) close(37)!for debug and tuning purposes
       call stop_clock('runcg_uspp')

       deallocate(bec0,becm,becdrdiag)

       return
     END SUBROUTINE runcg_uspp
