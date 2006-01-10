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
   subroutine inner_loop( nfi, tfirst, tlast, eigr,  irb, eigrb, &
      rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac,c0,bec  )

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
                                nx => nbspx, n => nbsp, ispin 

      use ensemble_dft, only: tens, tgrand, ninner, ismear, etemp, ef,       &
     &                tdynz, tdynf, zmass, fmass, fricz, fricf, z0, c0diag,  &
                      becdiag, fmat0,    fion2, atot0,   &
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
      use efield_module, only: tefield, evalue, ctable, qmat, detq, ipolp, &
            berry_energy, ctabin, gqq, gqqm, df, pberryel
      use mp, only: mp_sum,mp_bcast
!
      implicit none
!
      integer :: nfi
      logical :: tfirst , tlast
      complex(DP) :: eigr(ngw,nat)
      complex(DP) :: c0(ngw,n)
      real(DP) :: bec(nhsa,n)
      integer irb(3,nat)
      complex(DP) :: eigrb(ngb,nat)
      real(DP) :: rhor(nnr,nspin)
      real(DP) :: rhog(ngm,nspin)
      real(DP) :: rhos(nnrsx,nspin)
      real(DP) :: rhoc(nnr)
      complex(DP) :: ei1(-nr1:nr1,nat)
      complex(DP) :: ei2(-nr2:nr2,nat)
      complex(DP) :: ei3(-nr3:nr3,nat)
      complex(DP) :: sfac( ngs, nsp )
!
!
      integer :: i, j, ig, k, is, ia, iv, jv, il
      integer :: inl, jnl, niter, istart, nss
      real(DP) :: enb, enbi, x
      complex(DP) :: c2( ngw )
      complex(DP) :: c3( ngw )
      real(DP) :: gamma, entmp, sta
!
!


      fion2=0.d0



  
        !=======================================================================
        !
        !                 start of the inner loop
        !                 (Uij degrees of freedom)
        !
        !=======================================================================
  
            
          if(.not. ene_ok) then
            !     calculation of the array bec:
            call calbec (1,nsp,eigr,c0,bec)
  
                   
            call rotate(z0,c0(:,:),bec,c0diag,becdiag)
  
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
              call dforce(bec,betae,i,c0(1,i), &
              &     c0(1,i+1),h0c0(1,i),h0c0(1,i+1),rhos)
            end do
                   
            do is=1,nspin
              nss=nupdwn(is)
              istart=iupdwn(is)
              do i=1,nss
                do k=1,nss
                  c0hc0(k,i,is)=0.d0
                  do ig=1,ngw
                    c0hc0(k,i,is)=c0hc0(k,i,is)- &
                 &   2.0*DBLE(CONJG(c0(ig,k+istart-1))*h0c0(ig,i+istart-1))
                  enddo
                  if (ng0.eq.2) then
                    c0hc0(k,i,is)=c0hc0(k,i,is)+&
                 &   DBLE(CONJG(c0(1,k+istart-1))*h0c0(1,i+istart-1))
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
              call rotate(zxt,c0(:,:),bec,c0diag,becdiag)   
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
              call dforce(bec,betae,i,c0(1,i),c0(1,i+1),&
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
                    2.0*DBLE(CONJG(c0(ig,k+istart-1))*h0c0(ig,i+istart-1))
                  enddo
                  if (ng0.eq.2) then
                    c0hc0(k,i,is)=c0hc0(k,i,is)+&
                    DBLE(CONJG(c0(1,k+istart-1))*h0c0(1,i+istart-1))
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
            call rotate(z0,c0(:,:),bec,c0diag,becdiag)

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
  

END SUBROUTINE
