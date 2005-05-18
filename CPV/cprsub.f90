!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
     subroutine caldbec(nspmn,nspmx,eigr,c)
!-----------------------------------------------------------------------
!     this routine calculates array dbec, derivative of bec:
!
!        < psi_n | beta_i,i > = c_n(0) beta_i,i(0) +
!                 2 sum_g> re(c_n*(g) (-i)**l beta_i,i(g) e^-ig.r_i) 
!
!     with respect to cell parameters h
!
!     routine makes use of c(-g)=c*(g)  and  beta(-g)=beta*(g)
!
      use mp, only: mp_sum
      use mp_global, only: nproc
      use ions_base, only: na, nas => nax, nat
      use electrons_base, only: n => nbsp
      use gvecw, only: ngw
      use reciprocal_vectors, only: gstart
      use constants, only: pi, fpi
      use cvan, only: ish
      use uspp, only: nhtol, nhsa => nkb
      use uspp_param, only: nh, nhm
      use cdvan
!
      implicit none
      integer nspmn, nspmx
      complex(kind=8) c(ngw,n)
      real(kind=8) eigr(2,ngw,nat)
      complex(kind=8), allocatable :: wrk2(:,:)
!
      integer ig, is, iv, ia, l, ixr, ixi, inl, i, j, ii, isa
      real(kind=8) signre, signim, arg
!
      allocate( wrk2( ngw, nas ) )
!
      do j=1,3
         do i=1,3

            isa = 0
            do is = 1, nspmn - 1
              isa = isa + na(is)
            end do

            do is=nspmn,nspmx
               do iv=1,nh(is)
                  l=nhtol(iv,is)
                  if (l == 0) then
                     ixr = 1
                     ixi = 2
                     signre =  1.0
                     signim =  1.0
                  else if (l == 1) then
                     ixr = 2
                     ixi = 1
                     signre =  1.0
                     signim = -1.0
                  else if (l == 2) then
                     ixr = 1
                     ixi = 2
                     signre = -1.0
                     signim = -1.0
                  else if (l == 3) then
                     ixr = 2
                     ixi = 1
                     signre = -1.0
                     signim =  1.0
                  endif
                  !
                  do ia=1,na(is)
                     if (gstart == 2) then
!                   q = 0   component (with weight 1.0)
                        wrk2(1,ia)= cmplx(                             &
     &                  signre*dbeta(1,iv,is,i,j)*eigr(ixr,1,ia+isa),   &
     &                  signim*dbeta(1,iv,is,i,j)*eigr(ixi,1,ia+isa) )
!                   q > 0   components (with weight 2.0)
                     end if
                     do ig=gstart,ngw
                        arg = 2.0*dbeta(ig,iv,is,i,j)
                        wrk2(ig,ia) = cmplx(                           &
     &                         signre*arg*eigr(ixr,ig,ia+isa),          &
     &                         signim*arg*eigr(ixi,ig,ia+isa) )
                     end do
                  end do
                  inl=ish(is)+(iv-1)*na(is)+1
                  call MXMA(wrk2,2*ngw,1,c,1,2*ngw,dbec(inl,1,i,j),1,  &
     &                      nhsa,na(is),2*ngw,n)
               end do
               if( nproc > 1 ) then
                  inl=ish(is)+1
                  do ii=1,n
                     call mp_sum( dbec( inl : (inl + na(is)*nh(is) - 1), ii,i,j) )
                  end do
               end if
               isa = isa + na(is)
            end do
         end do
      end do

      deallocate( wrk2 )
!
      return
      end subroutine caldbec
!

!-----------------------------------------------------------------------
      subroutine formf(tfirst, eself)
!-----------------------------------------------------------------------
!computes (a) the self-energy eself of the ionic pseudocharges;
!         (b) the form factors of: (i) pseudopotential (vps),
!             (ii) ionic pseudocharge (rhops)
!         all quantities are returned in common /pseu/
!         also calculated the derivative of vps with respect to
!         g^2 (dvps)
! 
      use mp,            ONLY: mp_sum
      use control_flags, only: iprint, tpre, iprsta
      use io_global,     only: stdout
      use bhs,           only: rc1, rc2, wrc2, wrc1, rcl, al, bl, lloc
      use gvecs,         only: ngs
      use cell_base,     only: omega, tpiba2
      use ions_base,     only: rcmax, zv, nsp, na
      use cvan,          only: oldvan
      use local_pseudo,  only: vps, rhops, dvps, drhops
      use atom,          only: r, rab, mesh, numeric
      use uspp_param,    only: vloc_at
      use qrl_mod,       only: cmesh
      use pseudo_base,   only: compute_rhops, formfn, formfa, compute_eself
      use reciprocal_vectors, only: gstart, g
!
      implicit none
      logical      :: tfirst
      real(kind=8) :: eself
      !
      real(kind=8) :: vpsum, rhopsum
      integer      :: is

      call start_clock( 'formf' )
      !
      ! calculation of gaussian selfinteraction
      !
      eself = compute_eself( na, zv, rcmax, nsp )

      if( tfirst .or. iprsta .ge. 4 )then
         WRITE( stdout,1200) eself
      endif
 1200 format(2x,'formf: eself=',f10.5)
!
      do is=1,nsp

         if (numeric(is)) then

            call formfn( vps(:,is), dvps(:,is), r(:,is), rab(:,is), vloc_at(:,is), &
                         zv(is), rcmax(is), g, omega, tpiba2, cmesh(is), mesh(is), &
                         ngs, oldvan(is), tpre )

         else

            !     bhs pseudopotentials
            !
            call formfa( vps(:,is), dvps(:,is), rc1(is), rc2(is), wrc1(is), wrc2(is), &
                         rcl(:,is,lloc(is)), al(:,is,lloc(is)), bl(:,is,lloc(is)),    &
                         zv(is), rcmax(is), g, omega, tpiba2, ngs, gstart, tpre )

         end if


         !     fourier transform of local pp and gaussian nuclear charge
         !
         call compute_rhops( rhops(:,is), drhops(:,is), zv(is), rcmax(is), g,   &
                             omega, tpiba2, ngs, tpre )

         if(tfirst.or.(iprsta.ge.4))then
            vpsum = SUM( vps( 1:ngs, is ) )
            rhopsum = SUM( rhops( 1:ngs, is ) )
            call mp_sum( vpsum )
            call mp_sum( rhopsum )
            WRITE( stdout,1250) vps(1,is),rhops(1,is)
            WRITE( stdout,1300) vpsum,rhopsum
         endif
!
      end do
!
      call stop_clock( 'formf' )
!
 1250 format(2x,'formf:     vps(g=0)=',f12.7,'     rhops(g=0)=',f12.7)
 1300 format(2x,'formf: sum_g vps(g)=',f12.7,' sum_g rhops(g)=',f12.7)
!
      return
      end subroutine formf
!


!-----------------------------------------------------------------------
      subroutine newnlinit
!-----------------------------------------------------------------------
!     this routine calculates arrays beta, qradb, qq, qgb, rhocb
!     and derivatives w.r.t. cell parameters dbeta, dqrad 
!     See also comments in nlinit
!
      use control_flags, only: iprint, tpre, iprsta
      use io_global, only: stdout
      !use gvec
      use gvecw, only: ngw
      use reciprocal_vectors, only: g, gx
      use reciprocal_vectors, only: gstart
      use cell_base, only: omega, tpiba2, tpiba
      use cell_base, only: ainv
      use cvan, only: nvb
      use uspp, only: qq, nhtolm, beta
      use uspp_param, only: nh
      use core
      use constants, only: pi, fpi
      use ions_base, only: nsp
      use uspp_param, only: lmaxq, nqlc, lmaxkb, kkbeta, nbrx, nbeta
      use atom, only: nlcc, r, rab, mesh, rho_atc
      use qradb_mod
      use qgb_mod
      use gvecb
      use small_box,  only: omegab, tpibab
      use cdvan
      use dqrad_mod
      use dqgb_mod
      use betax
      use pseudo_base, only: compute_rhocg
!
      implicit none
      integer  is, l, lp, ig, ir, iv, jv, ijv, i,j, jj, ierr
      real(kind=8), allocatable:: fint(:), jl(:), dqradb(:,:,:,:,:)
      real(kind=8), allocatable:: ylmb(:,:), ylm(:,:), &
                                  dylmb(:,:,:,:), dylm(:,:,:,:)
      complex(kind=8), allocatable:: dqgbs(:,:,:)
      real(kind=8) xg, c, betagl, dbetagl, gg
!
! 
      allocate( ylmb( ngb, lmaxq*lmaxq ), STAT=ierr )
      IF( ierr  /= 0 ) &
        CALL errore(' newnlinit ', ' cannot allocate ylmb ', 1 )
!
      qradb(:,:,:,:,:) = 0.d0
      call ylmr2 (lmaxq*lmaxq, ngb, gxb, gb, ylmb)


!     ===============================================================
!     initialization for vanderbilt species
!     ===============================================================
      do is=1,nvb

!     ---------------------------------------------------------------
!     calculation of array qradb(igb,iv,jv,is)
!     ---------------------------------------------------------------
         if(iprsta.ge.4) WRITE( stdout,*)  '  qradb  '
         c=fpi/omegab
!
         do l=1,nqlc(is)
            do iv= 1,nbeta(is)
               do jv=iv,nbeta(is)
                  do ig=1,ngb
                     gg=gb(ig)*tpibab*tpibab/refg
                     jj=int(gg)+1
                     if(jj.ge.mmx) then
                        qradb(ig,iv,jv,l,is)=0.
                        qradb(ig,jv,iv,l,is)=qradb(ig,iv,jv,l,is)
                     else
                        qradb(ig,iv,jv,l,is)=                           &
     &                       c*qradx(jj+1,iv,jv,l,is)*(gg-real(jj-1))+  &
     &                       c*qradx(jj,iv,jv,l,is)*(real(jj)-gg)
                        qradb(ig,jv,iv,l,is)=qradb(ig,iv,jv,l,is)
                     endif
                  enddo
               enddo
            enddo
         enddo

!
!     ---------------------------------------------------------------
!     stocking of qgb(igb,ijv,is) and of qq(iv,jv,is)
!     ---------------------------------------------------------------
         ijv=0
         do iv= 1,nh(is)
            do jv=iv,nh(is)
!
!       compact indices because qgb is symmetric:
!       ivjv:  11 12 13 ... 22 23...
!       ijv :   1  2  3 ...  
!
               ijv=ijv+1
               call qvan2b(ngb,iv,jv,is,ylmb,qgb(1,ijv,is) )
!
               qq(iv,jv,is)=omegab*real(qgb(1,ijv,is))
               qq(jv,iv,is)=qq(iv,jv,is)
!
            end do
         end do

      end do
!
      if (tpre) then
!     ---------------------------------------------------------------
!     arrays required for stress calculation, variable-cell dynamics
!     ---------------------------------------------------------------
         allocate(dqradb(ngb,nbrx,nbrx,lmaxq,nsp))
         allocate(dylmb(ngb,lmaxq*lmaxq,3,3))
         allocate(dqgbs(ngb,3,3))
         dqrad(:,:,:,:,:,:,:) = 0.d0
         !
         call dylmr2_(lmaxq*lmaxq, ngb, gxb, gb, ainv, dylmb)
         !
         do is=1,nvb
            !
            do l=1,nqlc(is)
               do iv= 1,nbeta(is)
                  do jv=iv,nbeta(is) 
                    do ig=1,ngb
                        gg=gb(ig)*tpibab*tpibab/refg
                        jj=int(gg)+1
                        if(jj.ge.mmx) then
                           dqradb(ig,iv,jv,l,is) = 0.
                        else
                           dqradb(ig,iv,jv,l,is) =                      &
     &                       dqradx(jj+1,iv,jv,l,is)*(gg-real(jj-1))+   &
     &                       dqradx(jj,iv,jv,l,is)*(real(jj)-gg)
                        endif
                     enddo
                     do i=1,3
                        do j=1,3
                           dqrad(1,iv,jv,l,is,i,j) = &
                                -qradb(1,iv,jv,l,is) * ainv(j,i)
                           dqrad(1,jv,iv,l,is,i,j) = &
                                dqrad(1,iv,jv,l,is,i,j)
                           do ig=2,ngb
                              dqrad(ig,iv,jv,l,is,i,j) =                &
     &                          -qradb(ig,iv,jv,l,is)*ainv(j,i)         &
     &                          -c*dqradb(ig,iv,jv,l,is)*               &
     &                          gxb(i,ig)/gb(ig)*                       &
     &                          (gxb(1,ig)*ainv(j,1)+                   &
     &                           gxb(2,ig)*ainv(j,2)+                   &
     &                           gxb(3,ig)*ainv(j,3)) 
                              dqrad(ig,jv,iv,l,is,i,j) =                &
     &                          dqrad(ig,iv,jv,l,is,i,j)
                           enddo
                        enddo
                     enddo
                  end do
               enddo
            enddo
            !
            ijv=0
            !
            do iv= 1,nh(is)
               do jv=iv,nh(is)
                  !
                  !       compact indices because qgb is symmetric:
                  !       ivjv:  11 12 13 ... 22 23...
                  !       ijv :   1  2  3 ...  
                  !
                  ijv=ijv+1
                  call dqvan2b(ngb,iv,jv,is,ylmb,dylmb,dqgbs )
                  do i=1,3
                     do j=1,3
                        do ig=1,ngb
                           dqgb(ig,ijv,is,i,j)=dqgbs(ig,i,j)
                        enddo
                     enddo
                  enddo
               end do
            end do
         end do
         deallocate(dqgbs)
         deallocate(dylmb)
         deallocate(dqradb)
      end if
      deallocate(ylmb)
!
!     ===============================================================
!     initialization that is common to all species
!     ===============================================================
!
      allocate(ylm(ngw,(lmaxkb+1)**2))
      call ylmr2 ((lmaxkb+1)**2, ngw, gx, g, ylm)
!
      do is=1,nsp
!     ---------------------------------------------------------------
!     calculation of array  beta(ig,iv,is)
!     ---------------------------------------------------------------
         if(iprsta.ge.4) WRITE( stdout,*)  '  beta  '
         c=fpi/sqrt(omega)
         do iv=1,nh(is)
            lp=nhtolm(iv,is)
            do ig=1,ngw
               gg=g(ig)*tpiba*tpiba/refg
               jj=int(gg)+1
               betagl=betagx(jj+1,iv,is)*(gg-real(jj-1))+               &
     &              betagx(jj,iv,is)*(real(jj)-gg)
               beta(ig,iv,is)=c*ylm(ig,lp)*betagl
            end do
         end do
      end do
!
      if (tpre) then
!     ---------------------------------------------------------------
!     calculation of array dbeta required for stress, variable-cell
!     ---------------------------------------------------------------
         allocate(dylm(ngw,(lmaxkb+1)**2,3,3))
         !
         call dylmr2_((lmaxkb+1)**2, ngw, gx, g, ainv, dylm)
         !
         do is=1,nsp
            if(iprsta.ge.4) WRITE( stdout,*)  '  dbeta  '
            c=fpi/sqrt(omega)
            do iv=1,nh(is)
               lp=nhtolm(iv,is)
               betagl=betagx(1,iv,is)
               do i=1,3
                  do j=1,3
                     dbeta(1,iv,is,i,j)=-0.5*beta(1,iv,is)*ainv(j,i)    &
     &                                 +c*dylm(1,lp,i,j)*betagl
                  enddo
               enddo
               do ig=gstart,ngw
                  gg=g(ig)*tpiba*tpiba/refg
                  jj=int(gg)+1
                  betagl = betagx(jj+1,iv,is)*(gg-real(jj-1)) +         &
     &                     betagx(jj,iv,is)*(real(jj)-gg)
                  dbetagl= dbetagx(jj+1,iv,is)*(gg-real(jj-1)) +        &
     &                     dbetagx(jj,iv,is)*(real(jj)-gg)
                  do i=1,3
                     do j=1,3
                        dbeta(ig,iv,is,i,j)=                            &
     &                    -0.5*beta(ig,iv,is)*ainv(j,i)                 &
     &                    +c*dylm(ig,lp,i,j)*betagl                     &
     &                    -c*ylm (ig,lp)*dbetagl*gx(i,ig)/g(ig)         &
     &                    *(gx(1,ig)*ainv(j,1)+                         &
     &                      gx(2,ig)*ainv(j,2)+                         &
     &                      gx(3,ig)*ainv(j,3))
                     end do
                  end do
               end do
            end do
         end do
         !
         deallocate(dylm)
      end if
      !
      deallocate(ylm)
!     ---------------------------------------------------------------
!     non-linear core-correction   ( rhocb(ig,is) )
!     ---------------------------------------------------------------
      do is=1,nsp
         if(nlcc(is)) then
            CALL compute_rhocg( rhocb(:,is), rhocb(:,is), r(:,is), rab(:,is), &
                   rho_atc(:,is), gb, omegab, tpibab**2, kkbeta(is), ngb, 0 )
         endif
      end do
!
!
      return
      end subroutine newnlinit
!-----------------------------------------------------------------------
      subroutine nlfh(bec,dbec,lambda)
!-----------------------------------------------------------------------
!     contribution to the internal stress tensor due to the constraints
!
      use cvan, only: nvb, ish
      use uspp, only: nhsa => nkb, qq
      use uspp_param, only: nh, nhm
      use ions_base, only: na
      use electrons_base, only: nx => nbspx, n => nbsp
      use cell_base, only: omega, h
      use constants, only: pi, fpi
      use stre
!
      implicit none
      real(kind=8) bec(nhsa,n), dbec(nhsa,n,3,3), lambda(nx,nx)
!
      integer i, j, ii, jj, inl, iv, jv, ia, is
      real(kind=8) fpre(3,3), tmpbec(nhm,nx), tmpdh(nx,nhm), temp(nx,nx)
!
      fpre(:,:) = 0.d0
      do ii=1,3
         do jj=1,3
            do is=1,nvb
               do ia=1,na(is)
!
                  tmpbec(:, 1:n) = 0.d0
                  tmpdh (1:n, :) = 0.d0
!
                  do iv=1,nh(is)
                     do jv=1,nh(is)
                        inl=ish(is)+(jv-1)*na(is)+ia
                        if(abs(qq(iv,jv,is)).gt.1.e-5) then
                           do i=1,n
                              tmpbec(iv,i) = tmpbec(iv,i) +             &
     &                             qq(iv,jv,is)*bec(inl,i)
                           end do
                        endif
                     end do
                  end do
!
                  do iv=1,nh(is)
                     inl=ish(is)+(iv-1)*na(is)+ia
                     do i=1,n
                        tmpdh(i,iv)=dbec(inl,i,ii,jj)
                     end do
                  end do
!
                  if(nh(is).gt.0)then
                     temp(:, 1:n) = 0.d0
!
                     call MXMA                                          &
     &                    (tmpdh,1,nx,tmpbec,1,nhm,temp,1,nx,n,nh(is),n)
!
                     do j=1,n
                        do i=1,n
                           temp(i,j)=temp(i,j)*lambda(i,j)
                        end do
                     end do
!
                     fpre(ii,jj)=fpre(ii,jj)+2.*SUM(temp(1:n,1:n))
                  endif
!
               end do
            end do
         end do
      end do
      do i=1,3
         do j=1,3
            stress(i,j)=stress(i,j)+(fpre(i,1)*h(j,1)+                  &
     &           fpre(i,2)*h(j,2)+fpre(i,3)*h(j,3))/omega
         enddo
      enddo
!
      return
      end subroutine nlfh
!-----------------------------------------------------------------------
      subroutine nlinit
!-----------------------------------------------------------------------
!
!     this routine allocates and initalizes arrays beta, qradb, qq, qgb,
!     rhocb, and derivatives w.r.t. cell parameters dbeta, dqrad 
!
!       beta(ig,l,is) = 4pi/sqrt(omega) y^r(l,q^)
!                               int_0^inf dr r^2 j_l(qr) betar(l,is,r)
!
!       Note that beta(g)_lm,is = (-i)^l*beta(ig,l,is) (?)
!
!       qradb(ig,l,k,is) = 4pi/omega int_0^r dr r^2 j_l(qr) q(r,l,k,is)
!
!       qq_ij=int_0^r q_ij(r)=omega*qg(g=0)
!
!     beta and qradb are first calculated on a fixed linear grid in |G|
!     (betax, qradx) then calculated on the box grid by interpolation
!     (this is done in routine newnlinit)
!     
      use parameters, only: lmaxx
      use control_flags, only: iprint, tpre
      use io_global, only: stdout
      use gvecw, only: ngw
      use cvan, only: ish, nvb, oldvan
      use core, only: rhocb, nlcc_any
      use constants, only: pi, fpi
      use ions_base, only: na, nsp
      use uspp, only: aainit, beta, qq, dvan, nhtol, nhtolm, indv, &
           nhsa => nkb, nhsavb=>nkbus
      use uspp_param, only: kkbeta, qqq, nqlc, betar, nbrx, lmaxq, dion, &
           nbeta, lmaxkb, lll, nhm, nh, tvanp
      use qrl_mod, only: qrl, cmesh
      use atom, only: mesh, r, rab, nlcc, numeric
      use qradb_mod, only: qradb
      use qgb_mod, only: qgb
      use gvecb, only: ngb
      use cdvan, only: dbeta
      use dqrad_mod, only: dqrad
      use dqgb_mod, only: dqgb
      use betax, only: qradx, dqradx, refg, betagx, mmx, dbetagx
      use pseudopotential, only: pseudopotential_indexes, compute_dvan, &
            compute_betagx, compute_qradx

!
      implicit none
!
      integer  is, il, l, ir, iv, jv, lm, ind, ltmp, i0
      real(kind=8), allocatable:: fint(:), jl(:),  jltmp(:), djl(:),    &
     &              dfint(:)
      real(kind=8) xg, xrg, fac

      !
      !   initialize indexes
      !
      CALL pseudopotential_indexes()
      !
      !   initialize array ap
      !
      call aainit( lmaxkb + 1 )
!
      allocate(beta(ngw,nhm,nsp))
      allocate(qradb(ngb,nbrx,nbrx,lmaxq,nsp))
      allocate(qgb(ngb,nhm*(nhm+1)/2,nsp))
      allocate(qq(nhm,nhm,nsp))
      if (nlcc_any) allocate(rhocb(ngb,nsp))
!
      allocate(dqrad(ngb,nbrx,nbrx,lmaxq,nsp,3,3))
      allocate(dqgb(ngb,nhm*(nhm+1)/2,nsp,3,3))
      allocate(dbeta(ngw,nhm,nsp,3,3))
!
      qradb(:,:,:,:,:) = 0.d0
      qq  (:,:,:) =0.d0
      if(tpre) dqrad(:,:,:,:,:,:,:) = 0.d0
      !
      !     initialization for vanderbilt species
      !
      CALL compute_qradx( tpre )
      !    
      !     initialization that is common to all species
      !   
      WRITE( stdout, fmt="(//,3X,'Common initialization' )" )

      do is = 1, nsp
         WRITE( stdout, fmt="(/,3X,'Specie: ',I5)" ) is
         if ( .not. numeric(is) ) then
            fac=1.0
         else
            !     fac converts ry to hartree
            fac=0.5
         end if
         do iv = 1, nh(is)
            WRITE( stdout,901) iv, indv(iv,is), nhtol(iv,is)
         end do
 901     format(2x,i2,'  indv= ',i2,'   ang. mom= ',i2)
         !
         WRITE( stdout,*)
         WRITE( stdout,'(20x,a)') '    dion '
         do iv = 1, nbeta(is)
            WRITE( stdout,'(8f9.4)') ( fac * dion(iv,jv,is), jv = 1, nbeta(is) )
         end do
         !
      end do
      !
      !   calculation of array  betagx(ig,iv,is)
      !
      call compute_betagx( tpre )
      !
      !   calculate array  dvan(iv,jv,is)
      !
      call compute_dvan()
      !
      ! newnlinit stores qgb and qq, calculates arrays  beta  qradb  rhocb
      ! and derivatives wrt cell    dbeta dqrad
      !
      call newnlinit

      return
      end subroutine nlinit

!-------------------------------------------------------------------------
      subroutine qvan2b(ngy,iv,jv,is,ylm,qg)
!--------------------------------------------------------------------------
!     q(g,l,k) = sum_lm (-i)^l ap(lm,l,k) yr_lm(g^) qrad(g,l,l,k)
!
      use control_flags, only: iprint, tpre
      use qradb_mod
      use uspp, only: nlx, lpx, lpl, ap, indv, nhtolm
      use gvecb
      use cdvan
      use uspp_param, only: lmaxq
! 
      implicit none
      integer ngy, iv, jv, is
      complex(kind=8)   qg(ngb)
!
      integer ivs, jvs, ivl, jvl, i, ii, ij, l, lp, ig
      complex(kind=8) sig
      real(kind=8) :: ylm(ngb,lmaxq*lmaxq)
! 
!       iv  = 1..8     s_1 p_x1 p_z1 p_y1 s_2 p_x2 p_z2 p_y2
!       ivs = 1..4     s_1 s_2 p_1 p_2
!       ivl = 1..4     s p_x p_z p_y
! 
      ivs=indv(iv,is)
      jvs=indv(jv,is)
      ivl=nhtolm(iv,is)
      jvl=nhtolm(jv,is)
      if(ivl > nlx)  call errore(' qvan2b ',' ivl out of bounds ',ivl)
      if(jvl > nlx)  call errore(' qvan2b ',' jvl out of bounds ',jvl)
!
      qg(:) = (0.d0, 0.d0)
!
!     lpx = max number of allowed y_lm
!     lp  = composite lm to indentify them
!
      do i=1,lpx(ivl,jvl)
         lp=lpl(ivl,jvl,i)
         if (lp > lmaxq*lmaxq) call errore(' qvan2b ',' lp out of bounds ',lp)
!
!     extraction of angular momentum l from lp:  
!     l = int ( sqrt( float(l-1) + epsilon) ) + 1
!
         if (lp == 1) then
            l=1         
         else if ((lp >= 2) .and. (lp <= 4)) then
            l=2
         else if ((lp >= 5) .and. (lp <= 9)) then
            l=3
         else if ((lp >= 10).and.(lp <= 16)) then
            l=4
         else if ((lp >= 17).and.(lp <= 25)) then
            l=5
         else if ((lp >= 26).and.(lp <= 36)) then 
            l=6
         else if ((lp >= 37).and.(lp <= 49)) then 
            l=7
         else
            call errore(' qvan2b ',' not implemented ',lp)
         endif
!     
!       sig= (-i)^l
!
         sig=(0.,-1.)**(l-1)
         sig=sig*ap(lp,ivl,jvl)
         do ig=1,ngy
            qg(ig)=qg(ig)+sig*ylm(ig,lp)*qradb(ig,ivs,jvs,l,is)
         end do
      end do

      return
      end subroutine qvan2b

!-------------------------------------------------------------------------
      subroutine dqvan2b(ngy,iv,jv,is,ylm,dylm,dqg)
!--------------------------------------------------------------------------
!
!     dq(i,j) derivatives wrt to h(i,j) of q(g,l,k) calculated in qvan2b
!
      use control_flags, only: iprint, tpre
      use qradb_mod
      use uspp, only: nlx, lpx, lpl, ap, indv, nhtolm
      use gvecb
      use dqrad_mod
      use cdvan
      use uspp_param, only: lmaxq
! 
      implicit none
      integer ngy, iv, jv, is
      complex(kind=8) dqg(ngb,3,3)
!
      integer ivs, jvs, ivl, jvl, i, ii, ij, l, lp, ig
      complex(kind=8) sig
      real(kind=8) :: ylm(ngb,lmaxq*lmaxq), dylm(ngb,lmaxq*lmaxq,3,3)
! 
!       iv  = 1..8     s_1 p_x1 p_z1 p_y1 s_2 p_x2 p_z2 p_y2
!       ivs = 1..4     s_1 s_2 p_1 p_2
!       ivl = 1..4     s p_x p_z p_y
! 
      ivs=indv(iv,is)
      jvs=indv(jv,is)
      ivl=nhtolm(iv,is)
      jvl=nhtolm(jv,is)
      if(ivl > nlx)  call errore(' dqvan2b ',' ivl out of bounds ',ivl)
      if(jvl > nlx)  call errore(' dqvan2b ',' jvl out of bounds ',jvl)
!
      dqg(:,:,:) = (0.d0, 0.d0)
!
!     lpx = max number of allowed y_lm
!     lp  = composite lm to indentify them
!
      do i=1,lpx(ivl,jvl)
         lp=lpl(ivl,jvl,i)
         if (lp > lmaxq*lmaxq) call errore(' dqvan2b ',' lp out of bounds ',lp)
!
!     extraction of angular momentum l from lp:  
!     l = int ( sqrt( float(l-1) + epsilon) ) + 1
!
         if (lp == 1) then
            l=1         
         else if ((lp >= 2) .and. (lp <= 4)) then
            l=2
         else if ((lp >= 5) .and. (lp <= 9)) then
            l=3
         else if ((lp >= 10).and.(lp <= 16)) then
            l=4
         else if ((lp >= 17).and.(lp <= 25)) then
            l=5
         else if ((lp >= 26).and.(lp <= 36)) then 
            l=6
         else if ((lp >= 37).and.(lp <= 49)) then 
            l=7
         else
            call errore(' qvan2b ',' not implemented ',lp)
         endif
!     
!       sig= (-i)^l
!
         sig=(0.,-1.)**(l-1)
         sig=sig*ap(lp,ivl,jvl)
         do ij=1,3
            do ii=1,3
               do ig=1,ngy
                  dqg(ig,ii,ij) = dqg(ig,ii,ij) +  sig *                &
     &                    ( ylm(ig,lp) * dqrad(ig,ivs,jvs,l,is,ii,ij) + &
     &                     dylm(ig,lp,ii,ij)*qradb(ig,ivs,jvs,l,is)   )
               end do
            end do
         end do
      end do
!
      return
      end subroutine dqvan2b
!-----------------------------------------------------------------------
subroutine dylmr2_(nylm, ngy, g, gg, ainv, dylm)
  !-----------------------------------------------------------------------
  !
  ! temporary CP interface for PW routine dylmr2
  ! dylmr2  calculates d Y_{lm} /d G_ipol
  ! dylmr2_ calculates G_ipol \sum_k h^(-1)(jpol,k) (dY_{lm} /dG_k)
  !
  USE kinds
  implicit none
  !
  integer, intent(IN) :: nylm, ngy
  real(kind=DP), intent(IN) :: g (3, ngy), gg (ngy), ainv(3,3)
  real(kind=DP), intent(OUT) :: dylm (ngy, nylm, 3, 3)
  !
  integer :: ipol, jpol, lm
  real(kind=DP), allocatable :: dylmaux (:,:,:)
  !
  allocate ( dylmaux(ngy,nylm,3) )
  dylmaux(:,:,:) = 0.d0
  do ipol =1,3
     call dylmr2 (nylm, ngy, g, gg, dylmaux(1,1,ipol), ipol)
  enddo
  !
  do ipol =1,3
     do jpol =1,3
        do lm=1,nylm
           dylm (:,lm,ipol,jpol) = (dylmaux(:,lm,1) * ainv(jpol,1) + &
                                    dylmaux(:,lm,2) * ainv(jpol,2) + &
                                    dylmaux(:,lm,3) * ainv(jpol,3) ) &
                                    * g(ipol,:)
        end do
     end do
  end do
  deallocate ( dylmaux )
end subroutine dylmr2_
