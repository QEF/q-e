!
! Copyright (C) 2002 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!
!=======================================================================
!
   subroutine runcg_uspp( nfi, tfirst, tlast, eigr, bec, irb, eigrb, &
      rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac, fion, ema0bg, becdr, &
      lambdap, lambda  )


      use control_flags, only: iprint, thdyn, tpre, tbuff, iprsta, trhor, &
            tfor, tvlocw, trhow, taurdr, tprnfor
      use control_flags, only: ndr, ndw, nbeg, nomore, tsde, tortho, tnosee, &
            tnosep, trane, tranp, tsdp, tcp, tcap, ampre, amprp, tnoseh

      use atom, only: nlcc
      use core, only: nlcc_any
      use core, only: deallocate_core
!---ensemble-DFT
      use energies, only: eht, epseu, exc, etot, eself, enl, ekin,          &
     &                    atot, entropy, egrand
      use electrons_base, only: f, nspin, nel, iupdwn, nupdwn, nudx, nelt, &
                                nx => nbspx, n => nbsp, ispin => fspin, &
                                deallocate_elct

      use ensemble_dft, only: tens, tgrand, ninner, ismear, etemp, ef,              &
     &                tdynz, tdynf, zmass, fmass, fricz, fricf, z0, c0diag, &
                      becdiag, fmat0, becdrdiag, becm, bec0, fion2, atot0, &
                      etot0, h0c0, c0hc0, epsi0, e0, dval, z1, f1, dfmat, fmat1, &
                      ef1, enocc, f0, fmatx, fx, zaux, zx, ex, zxt, atot1, etot1, &
                      dedx1, dentdx1, dx, dadx1, faux, eqc, eqa, atotmin, xmin, &
                      entropy2, f2, etot2, eqb, compute_entropy2, compute_entropy_der, &
                      compute_entropy
!---
      use gvec, only: tpiba2, ng
      use gvec, only: deallocate_gvec
      use gvecs, only: ngs
      use gvecb, only: ngb
      use gvecw, only: ngw
      use reciprocal_vectors, only: ng0 => gstart
      use cvan, only: ipp, nvb, ish
      use ions_base, only: na, nat, pmass, nax, nsp, rcmax
      use grid_dimensions, only: nnr => nnrx, nr1, nr2, nr3
      use cell_base, only: ainv, a1, a2, a3
      use cell_base, only: omega, alat
      use cell_base, only: h, hold, deth, wmass
      use smooth_grid_dimensions, only: nnrsx, nr1s, nr2s, nr3s
      use smallbox_grid_dimensions, only: nnrb => nnrbx, nr1b, nr2b, nr3b
      use pseu, only: vps, rhops
      use pseu, only: deallocate_pseu
      use io_global, ONLY: io_global_start, stdout, ionode
      use mp_global, ONLY: mp_global_start
      use mp, ONLY: mp_end
      use para_mod
      use dener
      use derho
      use dpseu
      use cdvan
      use stre
      use gvecw, only: ggp, agg => ecutz, sgg => ecsig, e0gg => ecfix
      use parameters, only: nacx, natx, nsx, nbndxx
      use constants, only: pi, factem
      use io_files, only: psfile, pseudo_dir
      use qgb_mod, only: deallocate_qgb_mod
      use dqgb_mod, only: deallocate_dqgb_mod
      use qradb_mod, only: deallocate_qradb_mod
      use dqrad_mod, only: deallocate_dqrad_mod
      use betax, only: deallocate_betax
      use input_parameters, only: outdir

      use uspp, only : nhsa=> nkb, betae => vkb, rhovan => becsum, deeq
      use uspp_param, only: nh
      use cg_module, only: ltresh, itercg, etotnew, etotold, tcutoff, &
          restartcg, passof, passov, passop, ene_ok, numok, maxiter, &
          enever, etresh, ene0, hpsi, gi, hi, esse, essenew, dene0, spasso, &
          ene1, passo, iter3, enesti, ninner_ef
      use ions_positions, only: tau0
      use wavefunctions_module, only: c0, cm, phi => cp
      use wavefunctions_module, only: deallocate_wavefunctions
      use efield_module, only: evalue, ctable, qmat, detq, ipolp, &
            berry_energy, ctabin, gqq, gqqm, df
      use mp, only: mp_sum
!
      implicit none
!
      integer :: nfi
      logical :: tfirst , tlast
      complex(kind=8) :: eigr(ngw,nax,nsp)
      real(kind=8) :: bec(nhsa,n)
      real(kind=8) :: becdr(nhsa,n,3)
      integer irb(3,natx,nsx)
      complex(kind=8) :: eigrb(ngb,nax,nsp)
      real(kind=8) :: rhor(nnr,nspin)
      real(kind=8) :: rhog(ng,nspin)
      real(kind=8) :: rhos(nnrsx,nspin)
      real(kind=8) :: rhoc(nnr)
      complex(kind=8) :: ei1(-nr1:nr1,nax,nsp)
      complex(kind=8) :: ei2(-nr2:nr2,nax,nsp)
      complex(kind=8) :: ei3(-nr3:nr3,nax,nsp)
      complex(kind=8) :: sfac( ngs, nsp )
      real(kind=8) :: fion(3,natx)
      real(kind=8) :: ema0bg(ngw)
      real(kind=8) :: lambdap(nx,nx)
      real(kind=8) :: lambda(nx,nx)
!
!
      integer :: i, j, ig, k, is, ia, iv, jv, il
      integer :: inl, jnl, niter, istart, nss
      real(kind=8) :: enb, enbi, x
      complex(kind=8) :: c2( ngw )
      complex(kind=8) :: c3( ngw )
      real(kind=8) :: gamma, entmp, sta
!
!
      if(tfirst) write(6,*) 'GRADIENTE CONIUGATO'

      ! in caso ismear=0 o =-1 mette tutte le f a posto:


      call prefor(eigr,betae) !basta na volta per tute in cg, verfica

      ltresh    = .false.
      itercg    = 1
      etotold   = 1.d8
      tcutoff   = .false.
      restartcg = .true.
      passof = passop
      ene_ok = .false.


      !ortonormalizza c0

      call calbec(1,nsp,eigr,c0,bec)

      call graham(betae,bec,c0)

      call calbec(1,nsp,eigr,c0,bec)

      !calcola phi per pcdaga

      call calphiid(c0,bec,betae,phi)
         
      !setta indice su numero passi convergiuti

      numok = 0

      do while ( itercg .lt. maxiter .and. (.not.ltresh) )

        if( tfirst ) write(6,*)  'Iterazione numero:', itercg

        !parte da c0, calcola bec, carica,potenziale
        !la struttura dipendente da posizioni atomiche e' gia' stata calcolata
        !ATTENZIONE estendere enever a caso metallico

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

          !calcula il potenzialo
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

          if(evalue .ne. 0.d0 ) then
            
             call berry_energy( enb, enbi, bec, c0(:,:,1,1), fion )
          
             enb=enb*evalue
             enbi=enbi*evalue
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

        !se prima iterazione diagonalizza stati

        !           if(itercg.eq.2) then
        !            do i=1,n
        !              do ig=1,ngw
        !                 c0(ig,i,1,1)=c0diag(ig,i)
        !              enddo
        !            enddo
        !            z0=id
        !            restartcg=.true.
        !           endif 


        !calcola el preconditioning dipendente da banda come in articolo

        if(abs(etotnew-etotold).lt.etresh) then
           numok=numok+1
        else 
           numok=0
        endif

        if(numok.ge.4) then
           ltresh=.true.
        endif

        !per calcolo stato eccitato rifa un giro
        !            else

        if(iprsta.gt.1) then
          if(etotnew.lt.etotold) then
            write(6,*) 'Energia     TOTALE  :',itercg, etotnew, etotnew-etotold
          else
            write(6,*) 'PORCO CAZZO TOTALE  :',itercg, etotnew
          endif
        endif

        !  if(abs(etotnew-etotold).lt.0.0001) tcutoff=.false.
        etotold=etotnew
        ene0=etot

        !  Non usa piu' emme, dovaria eser za bastanza preciso!!!             
        !               call cal_emme(c0,bec,emme, 1)
        !calcula nove d

        call newd(rhor,irb,eigrb,rhovan,deeq,fion)

        !calcula el gradiente al paso sucesivo, e calcula la soma energie de ks:

        call prefor(eigr,betae)!ATTENZIONE

        do i=1,n,2
          call dforce(bec,deeq,betae,i,c0(1,i,1,1),c0(1,i+1,1,1),c2,c3,rhos)
          if(evalue.ne.0.d0) then
            call dforceb(c0, i, betae, ipolp, bec ,ctabin(1,1,ipolp), gqq, gqqm, qmat, deeq, df)
            c2(1:ngw)=c2(1:ngw)+evalue*df(1:ngw)
            call dforceb(c0, i+1, betae, ipolp, bec ,ctabin(1,1,ipolp), gqq, gqqm, qmat, deeq, df)
            c3(1:ngw)=c3(1:ngw)+evalue*df(1:ngw)
          endif
          hpsi(1:ngw,  i)=c2(1:ngw)
          hpsi(1:ngw,i+1)=c3(1:ngw)
          if (ng0.eq.2) then
            hpsi(1,  i)=cmplx(real(hpsi(1,  i)))
            hpsi(1,i+1)=cmplx(real(hpsi(1,i+1)))
          end if
        enddo

        gi(1:ngw,1:n) = hpsi(1:ngw,1:n)
               
        !la riga sotto equivale a mettere il termine in lambda, vedi note              
        !               call pc_emmedaga(c0,phi,gi,emme)

        call pcdaga2(c0,phi,gi)

        DO i = 1, n
          gi(1:ngw,i) = gi(1:ngw,i) * ema0bg(1:ngw)
        END DO

        !se preconditioning su G calcola lambda per termine aggiunto in gradiente 
        !vedi Teter,et al. PRB 40,12255              
        !nel caso US ci sta S|psi_j>
        !               if(.not.tcutoff) then
        !                  do i=1,n
        !                     do j=i,n
        !                        lambda(i,j)=0.d0
        !                        do ig=1,ngw
        !                           lambda(i,j)=lambda(i,j)-2.d0*real(conjg(c0(ig,i,1,1))*hpsi(ig,j))
        !                        enddo
        !                        if(ng0.eq.2) then
        !                           lambda(i,j)=lambda(i,j)+real(conjg(c0(1,i,1,1))*hpsi(1,j))
        !                        endif
        !                        lambda(j,i)=lambda(i,j)
        !                     enddo
        !                  enddo
        !#ifdef __PARA
        !                  call reduce(nx*n,lambda)
        !#endif
        !               endif

        call calcmt(f,z0,fmat0)

        if(iprsta.gt.10) then!ATTENZIONE

          !stampa forza su ioni
          ! calculation of  contribution of the non-local part of the pseudopotential
          ! to the force on each ion
          if (.not.tens) then
            if (tfor .or. tprnfor) call nlfq(c0,deeq,eigr,bec,becdr,fion)
          else
            if (tfor .or. tprnfor) call nlfq(c0diag,deeq,eigr,becdiag,becdrdiag,fion)
            !EL PARAMETRO becdrdiag el xe in output, tuto bon...
          endif

          !aggiunge parte dipendente da lambda in caso US
          if(nvb.ge.1) then
            do i=1,n
              do j=1,n
                 lambda(i,j)=0.d0
                 do ig=1,ngw
                    lambda(i,j)=lambda(i,j)-2.d0*real(conjg(c0(ig,i,1,1))*gi(ig,j))
                 enddo
                 if(ng0.eq.2) then
                    lambda(i,j)=lambda(i,j)+real(conjg(c0(1,i,1,1))*gi(1,j))
                 endif
              enddo
            enddo

            call mp_sum(lambda)

            if(tens) then!caso us meti anca le f
               do i=1,n
                 do j=1,n
                   lambdap(i,j)=0.d0
                   do k=1,n
                     lambdap(i,j)=lambdap(i,j)+lambda(i,k)*fmat0(k,j,1)
                   end do
                 end do
               enddo
                   
               call nlsm2(eigr,c0(:,:,1,1),becdr)

            endif
            if(.not.tens) then
              call nlfl(bec,becdr,lambda,fion)
            else
              call nlfl(bec,becdr,lambdap,fion)
            endif

          end if


          do ia=1,nat
            write(6,*) 'F :', itercg,(fion(i,ia),i=1,3)
          enddo

        endif !su iprsta !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        call calbec(1,nsp,eigr,gi,becm)

        call pcdaga2(c0,phi,gi)
        call pcdaga2(c0,phi,hpsi) 

        !              call pc_emmedaga(c0,phi,gi,emme)!ATTENZIONE
        !              call pc_emmedaga(c0,phi,hpsi,emme)
        !caso prima iterazion

        if(itercg==1.or.(mod(itercg,20).eq.1).or.restartcg) then

          restartcg=.false.
          passof=passop
          hi(1:ngw,1:n)=gi(1:ngw,1:n)

          !calcula esse per la seconda interazion

          gamma=0.d0
          if(.not.tens) then
            call calbec(1,nsp,eigr,gi,becm)
            do i=1,n        
              do ig=1,ngw
                gamma=gamma+2*real(conjg(gi(ig,i))*gi(ig,i))
              enddo
              if (ng0.eq.2) then
                gamma=gamma-real(conjg(gi(1,i))*gi(1,i))
              endif
            enddo
            call mp_sum(gamma)
             
            if (nvb.gt.0) then
               do is=1,nvb
                  do iv=1,nh(is)
                     do jv=1,nh(is)
                        do ia=1,na(is)
                           inl=ish(is)+(iv-1)*na(is)+ia
                           jnl=ish(is)+(jv-1)*na(is)+ia
                           !      gamma=gamma+ qq(iv,jv,is)*becm(inl,i)*becm(jnl,i)  
                        end do
                     end do
                  end do
               end do
            endif

          else

            do i=1,n        
               do j=1,n
                  do ig=1,ngw
                     gamma=gamma+2*real(conjg(gi(ig,i))*gi(ig,j))*fmat0(j,i,1)
                  enddo
                  if (ng0.eq.2) then
                     gamma=gamma-real(conjg(gi(1,i))*gi(1,j))*fmat0(j,i,1)
                  endif
               enddo
            enddo

            call mp_sum(gamma)

          endif

          esse=gamma

        else

          !trova hi in caso generale
          !calcola gamma  caso generale no Polak Ribiere

          gamma=0.d0

          if(.not.tens) then

            call calbec(1,nsp,eigr,gi,becm)
            do i=1,n        
               do ig=1,ngw
                  gamma=gamma+2*real(conjg(gi(ig,i))*gi(ig,i))
               enddo
               if (ng0.eq.2) then
                  gamma=gamma-real(conjg(gi(1,i))*gi(1,i))
               endif
            enddo
                
            call mp_sum(gamma)
             
            if (nvb.gt.0) then
             do is=1,nvb
                do iv=1,nh(is)
                   do jv=1,nh(is)      
                      do ia=1,na(is)
                         inl=ish(is)+(iv-1)*na(is)+ia
                         jnl=ish(is)+(jv-1)*na(is)+ia
                         ! gamma=gamma+ qq(iv,jv,is)*becm(inl,i)*becm(jnl,i)  
                      end do
                   end do
                end do
             end do
            endif

          else

            do i=1,n        
              do j=1,n
                do ig=1,ngw
                  gamma=gamma+2*real(conjg(gi(ig,i))*gi(ig,j))*fmat0(j,i,1)
                enddo
                if (ng0.eq.2) then
                  gamma=gamma-real(conjg(gi(1,i))*gi(1,j))*fmat0(j,i,1)
                endif
              enddo
            enddo
            call mp_sum(gamma)

          endif

          essenew=gamma
          gamma=gamma/esse
          esse=essenew

          hi(1:ngw,1:n)=gi(1:ngw,1:n)+gamma*hi(1:ngw,1:n)

        endif

        !minimizza lungo hi:

        !proietta hi su spazio conduzione

        call calbec(1,nsp,eigr,hi,bec0)
        call pc2(c0,bec,hi,bec0)
        !     call pc_emme(c0,bec,hi,bec0,emme)
        call calbec(1,nsp,eigr,hi,bec0)

        !ora minimizzazione parabolica
        !             
        !calcola derivata rispetto a lambda su direzione hi

        dene0=0.
        if(.not.tens) then
          do i=1,n               
            do ig=1,ngw
              dene0=dene0-4.d0*real(conjg(hi(ig,i))*hpsi(ig,i))!ATTENZION iera gi
            enddo
            if (ng0.eq.2) then
              dene0=dene0+2.d0*real(conjg(hi(1,i))*hpsi(1,i))!ATTENZION iera gi
            endif
          end do
        else
          !nel caso metalico la derivata xe' Sum_ij (<hi|H|Psi_j>+ <Psi_i|H|hj>)*f_ji
          !     calculation of the kinetic energy x=xmin      
          call calcmt(f,z0,fmat0)
          do i=1,n
            do j=1,n
              do ig=1,ngw
                dene0=dene0-2.d0*real(conjg(hi(ig,i))*hpsi(ig,j))*fmat0(j,i,1)
                !ATTENZIONE solo caso nspin=1!!!!!
                dene0=dene0-2.d0*real(conjg(hpsi(ig,i))*hi(ig,j))*fmat0(j,i,1)
              enddo
              if (ng0.eq.2) then
                dene0=dene0+real(conjg(hi(1,i))*hpsi(1,j))*fmat0(j,i,1)
                dene0=dene0+real(conjg(hpsi(1,i))*hi(1,j))*fmat0(j,i,1)
              end if
            enddo
          enddo
        endif

        call mp_sum(dene0)
        !            if(tens.and.(nspin.eq.1)) dene0=dene0/2.d0 

        !se la darivata la xe positiva, zerca in direzion oposta:

        if(dene0.gt.0.d0) then
          spasso=-1.D0
        else
          spasso=1.d0
        endif

        !calcula fni de onda in punto un poco spostado

        cm(1:ngw,1:n,1,1)=c0(1:ngw,1:n,1,1)+spasso*passof*hi(1:ngw,1:n)

        ! ordina gli stati in base valore energia
        !              do i=1,n
        !                  e0(i)=lambda(i,i)
        !              enddo
        !              call  ordina(cm,e0)

        !le ortonormalizza

        call calbec(1,nsp,eigr,cm,becm)
        call graham(betae,becm,cm)
        !              call riordina(cm,e0) 
        call calbec(1,nsp,eigr,cm,becm)
               
        !calcula energia del passetto
        if(.not.tens) then
          call rhoofr(nfi,cm,irb,eigrb,becm,rhovan,rhor,rhog,rhos,enl,ekin)
        else
          !     calculation of the rotated quantities
          call rotate(z0,cm(:,:,1,1),becm,c0diag,becdiag)
          !     calculation of rho corresponding to the rotated wavefunctions
          call rhoofr(nfi,c0diag,irb,eigrb,becdiag,rhovan,rhor,rhog,rhos,enl,ekin)
        endif

        !calcula il potenzialo
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

        if( evalue .ne. 0.d0 ) then
          call berry_energy( enb, enbi, becm, cm(:,:,1,1), fion )
          etot=etot+enb+enbi
        endif
        ene1=etot
              
            
        !trova il minimo di parabola

        call minparabola(ene0,spasso*dene0,ene1,passof,passo,enesti)

        if(iprsta.gt.1) write(6,*) ene0,dene0,ene1,passo, gamma, esse

        !imposta nuovo passetto come doppiodistanza da minimo

        passov=passof
        passof=2.d0*passo
              
        !calcola f.ni donda al minimo e ortonormalizza
        !adesso xe c00

        cm(1:ngw,1:n,1,1)=c0(1:ngw,1:n,1,1)+spasso*passo*hi(1:ngw,1:n)
        if(ng0.eq.2) then
          cm(1,:,1,1)=0.5d0*(cm(1,:,1,1)+conjg(cm(1,:,1,1)))
        endif

        !               call ordina(cm,e0)               
        call calbec(1,nsp,eigr,cm,becm)
        call graham(betae,becm,cm)
        !              call riordina(cm,e0) 


        !calcola energia al minimo per vedere se e' veramente diminuita
        !siccome in caso metallico nella prossima iterazione zo, cambia
        !il test viene fatto adesso

        !parte da c0, calcola bec, carica,potenziale
        !la struttura dipendente da posizioni atomiche e' gia' stata calcolata

        call calbec(1,nsp,eigr,cm,becm)
        if(.not.tens) then
          call rhoofr(nfi,cm,irb,eigrb,becm,rhovan,rhor,rhog,rhos,enl,ekin)
        else
          !     calculation of the rotated quantities
          call rotate(z0,cm(:,:,1,1),becm,c0diag,becdiag)
          !     calculation of rho corresponding to the rotated wavefunctions
          call rhoofr(nfi,c0diag,irb,eigrb,becdiag,rhovan,rhor,rhog,rhos,enl,ekin)
        endif

        !calcula il potenzialo
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
        if( evalue .ne. 0.d0 ) then
          call berry_energy( enb, enbi, becm, cm(:,:,1,1), fion )
          etot=etot+enb+enbi
        endif
        enever=etot
        !confronto con quanto previsto
        if(iprsta.gt.1) then
          write(6,*) 'Confr :'  , (enesti-enever)/(ene0-enever)
          write(6,*) 'Enever.'  , enever,passo,passov           
        endif

        !se l'energia e' diminuita rispetto a ene0 ene1 , tutto ok
        if( (enever.lt.ene0) .and. (enever.lt.ene1)) then
          c0(:,:,1,1)=cm(:,:,1,1)
          bec(:,:)=becm(:,:)
          ene_ok=.true.
          !se l'energia si e' alzata ma ene1 << ene0 va in ene1
          if( tfirst) write(6,*)  'Tutto ok:', itercg
        else if( (enever.ge.ene0).and.(ene0.gt.ene1)) then
          if(iprsta.gt.1) write(6,*) 'CASO: 2'
          c0(1:ngw,1:n,1,1)=c0(1:ngw,1:n,1,1)+spasso*passov*hi(1:ngw,1:n)
          restartcg=.true.!ATTENZIONE
          !                  call ordina(c0,e0)
          call calbec(1,nsp,eigr,c0,bec)
          call graham(betae,bec,c0)
          !                  call riordina(c0,e0)
          !se anche ene1 e' piu grande di ene0 fa un passo di gradiente coniugato,
          !riducendo il passetto in scala 2
        else if((enever.ge.ene0).and.(ene0.le.ene1)) then
          if(iprsta.gt.1) write(6,*) 'CASO: 3'
          iter3=0
          do while(enever.gt.ene0 .and. iter3.lt.4)
            iter3=iter3+1
            passov=passov*0.5d0
            cm(1:ngw,1:n,1,1)=c0(1:ngw,1:n,1,1)+spasso*passov*hi(1:ngw,1:n)
            !cambia la direzione su cui cerca, se la derivata 
            !e' molto piccola la direzione e' indeterminata
            spasso=spasso*(-1.d0)
            !                   call ordina(cm,e0)
            call calbec(1,nsp,eigr,cm,becm)
            call graham(betae,bec,cm)
            call calbec(1,nsp,eigr,cm,becm)
            if(.not.tens) then
              call rhoofr(nfi,cm,irb,eigrb,becm,rhovan,rhor,rhog,rhos,enl,ekin)
            else
              !     calculation of the rotated quantities
              call rotate(z0,cm(:,:,1,1),becm,c0diag,becdiag)
              !     calculation of rho corresponding to the rotated wavefunctions
              call rhoofr(nfi,c0diag,irb,eigrb,becdiag,rhovan,rhor,rhog,rhos,enl,ekin)
            endif
  
            !calcula il potenzialo
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
            if( evalue .ne. 0.d0 ) then
              call berry_energy( enb, enbi, becm, cm(:,:,1,1), fion )
              etot=etot+enb+enbi
            endif
            enever=etot
  
          end do
  
          c0(:,:,1,1)=cm(:,:,1,1)
          restartcg=.true.
  
        end if
  
        call calbec (1,nsp,eigr,c0,bec)
        !calcula phi per l'operator Pc daga
        call calphiid(c0,bec,betae,phi)
  
        !calcolare energia , non calcolarla di nuovo in tens e all'inizio
        ! se e' minore di stato 0 di stato 1 ok, se e' minore
        ! di 0 maggiore di 1, a a 1 e ricomincia con sd,
        ! se e' maggiore di 1 e 1 e' maggiore di 0, prendi
        !passo piu' piccolo ancora stepest descend, fino a che trova minimo,
        !poi ricomincia con steepest descend                  
             
  
        !***ensemble-DFT
        !=======================================================================
        !
        !                 start of the inner loop
        !                 (Uij degrees of freedom)
        !
        !=======================================================================
  
        ! per lo stato eccitato ismear=-1, e poi similmente per tutti stati eccitati
        ! se ltresh e' = .true. fa un passo cambiando le f
            
        if(tens) then
  
          if(.not. ene_ok) then
  
            !     calculation of the array bec:
            call calbec (1,nsp,eigr,c0,bec)
  
            !     calculation of ekinc
            call calcmt(f,z0,fmat0)
                   
            call rotate(z0,c0(:,:,1,1),bec,c0diag,becdiag)
  
            call rhoofr (nfi,c0diag,irb,eigrb,becdiag,rhovan,rhor,rhog,rhos,enl,ekin)
  
            !     put core charge (if present) in rhoc(r)
            if (nlcc_any) call set_cc(irb,eigrb,rhoc)
  
            !     calculation of the potential 
  
            call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,                &
                 ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion2)
            !     calculation of the array deeq: 
            !     deeq_i,lm = \int V_eff(r) q_i,lm(r) dr
          endif
  
          ene_ok=.false.
          call newd(rhor,irb,eigrb,rhovan,deeq,fion2)
    
          !     free energy at x=0
          atot0=atot    
          etot0=etot      
  
          !     start of the loop
          call prefor(eigr,betae)!ATTENZIONE
  
          !a causa stato eccitato introduzione ninner_ef
  
          ninner_ef=ninner
          if(ismear.eq.-1) then
            if(ltresh) ninner_ef = 1
          endif
  
          do niter=1,ninner_ef
  
            h0c0 = 0.0d0
            do i=1,n,2                      
              call dforce(bec,deeq,betae,i,c0(1,i,1,1), &
                   c0(1,i+1,1,1),h0c0(1,i),h0c0(1,i+1),rhos)
            end do
                   
            do is=1,nspin
              nss=nupdwn(is)
              istart=iupdwn(is)
              do i=1,nss
                do k=1,nss
                  c0hc0(k,i,is)=0.d0
                  do ig=1,ngw
                    c0hc0(k,i,is)=c0hc0(k,i,is)- &
                    2.0*real(conjg(c0(ig,k+istart-1,1,1))*h0c0(ig,i+istart-1))
                  enddo
                  if (ng0.eq.2) then
                    c0hc0(k,i,is)=c0hc0(k,i,is)+&
                    real(conjg(c0(1,k+istart-1,1,1))*h0c0(1,i+istart-1))
                  endif
                end do
              end do
            end do
  
            call mp_sum(c0hc0)
  
            do is=1,nspin
              nss=nupdwn(is)
              epsi0(1:nss,1:nss,is)=c0hc0(1:nss,1:nss,is)!ATTENZIONE
            end do
                      
            !    diagonalization of the matrix epsi0_ij
            e0 = 0.0d0
            do  is=1,nspin
              istart=iupdwn(is)
              nss=nupdwn(is) 
              call ddiag(nss,nss,epsi0(1,1,is),dval(1),z1(1,1,is),1)
              do i=1,nss
                e0(i+istart-1)=dval(i)
              enddo
            enddo
                   
            !    printing of the eigenvalues in fort.101
            !              write(101,*)   '============='
            !              write(101,'(2i10)') nfi,niter
            !              write(101,*) 'Eigenvalues(x=0)'
            !              do is=1,nspin
            !                 nss=nupdwn(is)
            !                 write(101,*) 'spin=',is
            !                 write(101,*) (e0(j),j=1,nss)
            !              enddo
  
            !     calculation of the occupations and the fermi energy
            !     corresponding to the chosen ismear,etemp and nspin
  
            if(ismear.eq.-1) then
              if( ltresh) then 
                call efermi(nelt,n,etemp,1,f1,ef1,e0,enocc,-1,nspin)
              else
                call efermi(nelt,n,etemp,1,f1,ef1,e0,enocc,0,nspin)
              endif
            else
              call efermi(nelt,n,etemp,1,f1,ef1,e0,enocc,ismear,nspin)
            endif
  
            !    printing of the occupations in fort.101
            !              write(101,*) 'Fermi energy(x=0):',ef1 
            !              write(101,*) 'Occupations(x=0)'
            !              do is=1,nspin
            !                 nss=nupdwn(is)
            !                 write(101,*) 'spin=',is
            !                 write(101,*) (f1(j),j=1,nss)
            !              enddo
  
            !     calculation of the initial and final occupation matrices
            !     in the z0-rotated orbital basis
              
            call calcm(f1,z1,fmat1)
             
            !     calculation of dfmat
            do is=1,nspin
              nss=nupdwn(is)
              dfmat(1:nss,1:nss,is)=-fmat0(1:nss,1:nss,is)+fmat1(1:nss,1:nss,is)
            end do
                   
            !     printing of the loop index in fort.100
            ! write(100,'(2i10)') nfi,niter
            !     edition of f(i)
            f0(1:n)=f(1:n)
                   
            !     initialization when xmin is determined by sampling 
            do il=1,1
              !     this loop is useful to check that the sampling is correct
              x=1.0*il
              do is=1,nspin
                nss=nupdwn(is)
                fmatx(1:nss,1:nss,is)=fmat0(1:nss,1:nss,is)+x*dfmat(1:nss,1:nss,is)
              end do
                      
              !     diagonalization of fmatx
              fx = 0.0d0
              do  is=1,nspin
                nss=nupdwn(is)
                istart=iupdwn(is)
                call ddiag(nss,nss,fmatx(1,1,is),dval(1),zaux(1,1,is),1)
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
              call newd(rhor,irb,eigrb,rhovan,deeq,fion2)
              call prefor(eigr,betae)
    
              !     free energy at x=1
              atot1=atot
              etot1=etot
              ! write(100,'(3f12.6)') x,atot1,entropy   
  
            end do
  
            !write(100,*) "____________________"
  
            !     calculation of c0hc0_ij at x=1
            call prefor(eigr,betae)!ATTENZIONE
            h0c0 = 0.0d0
            do i=1,n,2
              call dforce(bec,deeq,betae,i,c0(1,i,1,1),c0(1,i+1,1,1),&
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
                    2.0*real(conjg(c0(ig,k+istart-1,1,1))*h0c0(ig,i+istart-1))
                  enddo
                  if (ng0.eq.2) then
                    c0hc0(k,i,is)=c0hc0(k,i,is)+&
                    real(conjg(c0(1,k+istart-1,1,1))*h0c0(1,i+istart-1))
                  endif
                end do
              end do
            end do
            call mp_sum(c0hc0)
  
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
                   
            do il=0,200
                      
              x=0.005*il     
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
              ! write(100,'(3f12.6)') x,etot2+entropy2,entropy2
              if ((etot2+entropy2).lt.atotmin) then
                xmin=x
                atotmin=etot2+entropy2
              endif
  
            end do
  
            !le righe qua soto le xe per far el minimo con fit parabolico
  
            !                     eqc=atot0
            !                     eqa=dadx1-atot1+atot0
            !                     eqb=atot1-atot0-eqa
  
            !                     xmin=-eqb/(2.d0*eqa)
  
            if(ismear.eq.-1) then
              if(ltresh) xmin=1
            endif
  
            !     if xmin=1, no need do recalculate the quantities
             
            if (xmin.eq.1) goto 300
  
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
               call ddiag(nss,nss,fmatx(1,1,is),dval(1),zaux(1,1,is),1)  
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
            call newd(rhor,irb,eigrb,rhovan,deeq,fion2)
            call prefor(eigr,betae)
            ene_ok=.true. !so does not calculate the energy again
    
            !     free energy at x=xmin
            atotmin=atot      
  
            !     output
            !    write(100,'(a35,f12.7)') 'Fermi energy =',ef1
            !    write(100,'(a35,6f12.7)') 'xmin,atot0,atotmin,atot1,datot=',      &
            !         &          xmin,atot0,atotmin,atot1,atotmin-atot0
            !    write(100,'(a35,3f12.7)') 'eqa,eqb,eqc=',eqa,eqb,eqc 
            !    write(100,'(a35,3f12.7)') 'dadx1,dedx1,dentdx1=',dadx1,dedx1,     &
            !         &                          dentdx1
  
            if(iprsta.gt.1) write(6,*) 'Ciclo :', itercg,atot0,atot1
            atot0=atotmin
            etot0=etot
            !     end of the loop
  
          end do!su ninnner
      
          !=======================================================================
          !                 end of the inner loop
          !=======================================================================
  
        endif !su tens
  
        !            endif!su raggiungimento treshold, eliminato per stato eccitato
        itercg=itercg+1
  
      end do!su iterazioni cg
  
      !calcola forze se richiesto
      !       if(tfor .or. tprnfor) then
      !calcola sempre le forze in modo da calcolare lambda che serve per gli autostati
      call newd(rhor,irb,eigrb,rhovan,deeq,fion)
      if (.not.tens) then
        if (tfor .or. tprnfor) call nlfq(c0,deeq,eigr,bec,becdr,fion)
      else
        if (tfor .or. tprnfor) call nlfq(c0diag,deeq,eigr,becdiag,becdrdiag,fion)
        !EL PARAMETRO becdrdiag el xe in output, tuto bon...
      endif
  
      !aggiunge parte dipendente da lambda in caso US
  
      if(nvb.ge.1) then
        call prefor(eigr,betae)!ATTENZIONE
        do i=1,n,2
          call dforce(bec,deeq,betae,i,c0(1,i,1,1),c0(1,i+1,1,1),c2,c3,rhos)
          if(evalue .ne. 0.d0) then
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
            gi(1,  i)=cmplx(real(gi(1,  i)))
            gi(1,i+1)=cmplx(real(gi(1,i+1)))
          end if
        enddo
        do i=1,n
          do j=i,n
            lambda(i,j)=0.d0
            do ig=1,ngw
              lambda(i,j)=lambda(i,j)-2.d0*real(conjg(c0(ig,i,1,1))*gi(ig,j))
            enddo
            if(ng0.eq.2) then
              lambda(i,j)=lambda(i,j)+real(conjg(c0(1,i,1,1))*gi(1,j))
            endif
            lambda(j,i)=lambda(i,j)
          enddo
        enddo
  
        call mp_sum(lambda)
  
        if(tens) then!caso us meti anca le f
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
          call nlsm2(eigr,c0(:,:,1,1),becdr)
        endif
        call nlfl(bec,becdr,lambda,fion)
          
        ! bforceion adds the force term due to electronic berry phase
        ! only in US-case
          
        call bforceion(fion,tfor,ipolp, qmat,bec,becdr,gqq,evalue)
  
      endif
      !        end if

END SUBROUTINE
