program efg
  use kinds, only: DP
  use io_files, only: nd_nmbr,prefix, outdir, tmp_dir
  use parameters, only: ntypx
#ifdef __PARA 
  use para,       only : me 
  use mp, only: mp_bcast
#endif 

  implicit none
  character (len=80) :: filerec(ntypx)
  real(kind=DP) :: Q(ntypx)
  integer :: ios , ionode_id = 0

  namelist / inputpp / prefix, filerec, Q, outdir

  call start_postproc(nd_nmbr)
  !
  ! set default value
  !
  prefix = 'pwscf'
  outdir = './'
  Q=1.d0

#ifdef __PARA 
  if (me == 1)  then 
#endif 
  read (5, inputpp, err=200, iostat=ios)
200 call errore('efg.x', 'reading inputpp namelist', abs(ios))
  tmp_dir = trim(outdir)
#ifdef __PARA 
  end if
  ! 
  ! ... Broadcast variables 
  ! 
  CALL mp_bcast( prefix, ionode_id ) 
  CALL mp_bcast(tmp_dir, ionode_id ) 
  CALL mp_bcast(filerec, ionode_id )
  CALL mp_bcast(      Q, ionode_id )   
#endif 

  call read_file
  call openfil

#ifdef __PARA 
  if (me == 1)  then 
#endif 
  call read_recon(filerec)
#ifdef __PARA 
  end if
  ! 
  ! ... Broadcast variables 
  !
  CALL mp_bcast( paw_nbeta, ionode_id )
  CALL mp_bcast(     aephi, ionode_id )
  CALL mp_bcast(     psphi, ionode_id )
#endif 



  call do_efg(Q) 

  call stop_pp
  stop
end program efg

subroutine read_recon(filerec)

  use read_pseudo_module , only: scan_begin, scan_end
  use basis, only:ntyp
  use paw , only: aephi, psphi, paw_nbeta
  use atom, only: mesh 
  use kinds, only: DP
  use parameters, only : ntypx
  USE io_global,  ONLY : stdout
  implicit none

  character (len=80) :: filerec(ntypx)
  integer :: l,j,i,jtyp,kkphi,nbetam

  do jtyp=1,ntyp
     open(14,file=filerec(jtyp))
     call scan_begin(14,'PAW',.true.)
     read(14,*) paw_nbeta(jtyp)
     call scan_end(14,'PAW')
     close(14)
  enddo
     nbetam=maxval(paw_nbeta)
  allocate( psphi(ntyp,nbetam) )
  allocate( aephi(ntyp,nbetam) )


  recphi_read: do jtyp=1,ntyp
     open(14,file=filerec(jtyp))
     write (stdout,*) "N_AEwfc atom",jtyp,":",paw_nbeta(jtyp)
     recphi_loop: do i=1,paw_nbeta(jtyp)
        allocate(aephi(jtyp,i)%psi(maxval(mesh(1:ntyp))))
        aephi(jtyp,i)%label%nt=jtyp
        aephi(jtyp,i)%label%n=i
        call scan_begin(14,'REC',.false.)
        call scan_begin(14,'kkbeta',.false.)
        read(14,*)  kkphi
        call scan_end(14,'kkbeta')
        aephi(jtyp,i)%kkpsi=kkphi
        call scan_begin(14,'L',.false.)        
        read(14,*)  aephi(jtyp,i)%label%l
        call scan_end(14,'L')
        call scan_begin(14,'REC_AE',.false.)
        read(14,*) (aephi(jtyp,i)%psi(j),j=1,kkphi)
        call scan_end(14,'REC_AE')
        allocate (psphi(jtyp,i)%psi(maxval(mesh(1:ntyp))))
        psphi(jtyp,i)%label%nt=jtyp
        psphi(jtyp,i)%label%n=i
        psphi(jtyp,i)%label%l=aephi(jtyp,i)%label%l
        psphi(jtyp,i)%kkpsi=kkphi
        call scan_begin(14,'REC_PS',.false.)
        read(14,*) (psphi(jtyp,i)%psi(j),j=1,kkphi)
        call scan_end(14,'REC_PS')
        call scan_end(14,'REC')
     end do recphi_loop
     close(14)
  end do recphi_read
  
end subroutine read_recon

subroutine do_efg(Q)

  use io_files, only: nd_nmbr
  USE io_global,  ONLY : stdout
  use kinds ,only : DP 
  use parameters ,only: ntypx
  use constants, only: pi,tpi,fpi,ANGSTROM_AU,rytoev,ELECTRONVOLT_SI
  use scf, only: rho              !rho
  use gvect, only: nr1,nr2,nr3,nrx1,nrx1,nrx2,nrx3,nrxx,&
       g,gg,nl,gstart,ngm         !gvectors and parameters for the FFT
  use brilz , only: at,bg         !parameters of the cell
  use basis , only : nat, atm, tau, ityp         !coordinates of the atoms
  use symme , only: nsym, s, irt
  use pseud, only : zv !valence charge
  implicit none

  real(kind=DP) :: Q(ntypx), eta, Cq
  real(kind=DP) :: fac, trace, arg, e2
  integer :: alpha, beta, ig, na, i
  complex(kind=DP), allocatable:: aux(:)
  complex(kind=DP), allocatable:: efgg_el(:,:,:),efgr_el(:,:,:)
  complex(kind=DP), allocatable:: efg_io(:,:,:)
  real(kind=DP), allocatable:: zion(:), efg_corr_tens(:,:,:), efg(:,:,:)
  real(kind=DP):: efg_eig(3), v(3)
  complex(kind=DP) :: work(3,3), efg_vect(3,3)

  allocate(aux(nrxx))
  allocate(efgg_el(nrxx,3,3))
  allocate(efgr_el(nat,3,3))
  allocate(efg_io(nat,3,3))
  allocate(zion(nat))
  allocate(efg_corr_tens(3,3,nat))
  allocate(efg(3,3,nat))


!  e2 = 2.d0 ! rydberg
  e2 = 1.d0  ! hartree
  fac= fpi * e2
  aux(:)= rho(:,1)
  efgg_el(:,:,:)=(0.d0,0.d0)
  
  call cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,-1)

!
!calculation of the electic field gradient in the G-space
!
  do ig= gstart, ngm
     trace = 1.d0/3.d0 * gg(ig)
     do alpha=1,3
        efgg_el(ig,alpha,alpha)= -trace
        do beta=1,3
           efgg_el(ig,alpha,beta)=(efgg_el(ig,alpha,beta) + &
                g(alpha,ig) * g(beta,ig)) &
                * fac * (aux(nl(ig)))/gg(ig)

        enddo
     enddo
  enddo

! 
!fourier transform on the atomic position
!

efgr_el=(0.d0,0.d0)
  do alpha=1,3
     do beta=1,3
        do na=1,nat
           do ig= gstart, ngm
           arg=(tau(1,na)*g(1,ig)+tau(2,na)*g(2,ig)+tau(3,na)*g(3,ig))*tpi
           efgr_el(na,alpha,beta)=efgr_el(na,alpha,beta)+ &
                efgg_el(ig,alpha,beta) * cmplx(cos(arg),sin(arg))
           enddo
        enddo
     enddo 
  enddo

#ifdef __PARA
     call reduce (2*3*3*nat, efgr_el) !2*, efgr_el is a complex array
#endif

  write (stdout,*)
  
  do na=1,nat
     do beta=1,3

        write (stdout,1000) atm(ityp(na)),na,"efgr_el", &
             (real(efgr_el(na,alpha,beta)) , alpha =1,3 )

     enddo
     write (stdout,*)
  enddo


1000 FORMAT(1x,a,i3,2x,a,3(1x,f9.6))

!
! Ionic contribution
!

  call ewald_dipole (efg_io, zv)
  
  do na=1,nat
     do beta=1,3

        write (stdout,1000) atm(ityp(na)),na,"efg_ion", &
             (real(efg_io(1,alpha,beta)) , alpha =1,3 )

     enddo
     write (stdout,*)
  enddo

  call efg_correction(efg_corr_tens)

!symmetrize efg_tensor


  do na = 1,nat
     call trntns (efg_corr_tens(:,:,na),at, bg, -1)
  enddo
  call symz(efg_corr_tens, nsym, s, nat, irt)
  do na = 1,nat
     call trntns (efg_corr_tens(:,:,na),at, bg, 1)
  enddo

!
! print results
!

  do na=1,nat
     do beta=1,3

        write (stdout,1000) atm(ityp(na)),na,"efg_corr", &
             (2*real(efg_corr_tens(alpha,beta,1)) , alpha =1,3 )
        
     enddo
     write (stdout,*)
  enddo

  do na=1,nat
     efg(:,:,na)=real(2*efg_corr_tens(:,:,na)+efgr_el(na,:,:)+ &
          efg_io(na,:,:))
     do beta=1,3
        write (stdout,1000) atm(ityp(na)),na,"efg",&
             (efg(alpha,beta,na),alpha=1,3)
        
     enddo
     write (stdout,*)

  do alpha=1,3
     do beta=1,3
        work(beta,alpha)=cmplx(efg(alpha,beta,1),0.d0)
     enddo
  enddo

!
! diagonalise the tensor to extract the quadrupolar parameters Cq and eta
!

  call cdiagh(3,work,3,efg_eig,efg_vect)


  v(2)=efg_eig(2)
  if (abs(efg_eig(1))>abs(efg_eig(3))) then
     v(1)=efg_eig(1)
     v(3)=efg_eig(3)
  else
     v(1)=efg_eig(3)
     v(3)=efg_eig(1)
  endif

  if (abs(v(1))<1e-5) then
     eta=0.d0
  else
     eta=(v(2)-v(3))/v(1)
  endif

  Cq=v(1)*Q(ityp(na))*rytoev*2.d0*ANGSTROM_AU**2*ELECTRONVOLT_SI*1.e18/6.6262d0

  write (stdout,1200) atm(ityp(na)), na, Q(ityp(na)),Cq,eta
  write (stdout,*)
  enddo
1200 FORMAT(1x,a,1x,i3,5x,'Q= ',f6.3,5x,' Cq= ',f9.6,' MHz',5x,' eta= ',f9.6)

end subroutine do_efg
 
subroutine efg_correction(efg_corr_tens)

  use io_files ,only: nwordwfc, iunwfc
  use kinds , only: dp
  use us ,only: ap
  use parameters, only: lmaxx, ntypx
  use atom , only: r,rab,msh
  use gvect, only: g,ngm,ecutwfc
  use klist, only: nks, xk
  use cell_base, only: tpiba2
  use basis, only: nat, ntyp,ityp
  use wvfct, only:npwx, nbnd, npw, igk, g2kin
  use wavefunctions_module, only: evc
  use paw, only: paw_vkb, paw_becp, paw_nkb, aephi, psphi, paw_nh, paw_nhtol, &
       paw_nhtom, paw_indv, paw_nbeta
  use constants, only: pi

  implicit none

  integer :: j, ill, nt, ibnd, il1, il2, ik, iat, nbs1, nbs2, kkpsi
  integer :: lm,l,m,m1,m2,lm1,lm2, l1, l2,m3,m4,n,n1,nrc
  integer :: ijkb0,ijkb,ih,jh,na,np,ikb,jkb
  real(kind=dp), allocatable :: at_efg(:,:,:), work(:)
  real(kind=dp) ,intent(out):: efg_corr_tens(3,3,nat)
  complex(kind=dp) , allocatable :: efg_corr(:,:)
  real(kind=dp) :: rc

  allocate (efg_corr(lmaxx**2,nat))


  efg_corr=0.d0
  kkpsi=aephi(1,1)%kkpsi
  allocate (work(kkpsi))
 
  call init_paw_1

  allocate (at_efg(paw_nkb,paw_nkb,ntypx)) 
  allocate (paw_vkb( npwx,  paw_nkb))
  allocate (paw_becp(paw_nkb, nbnd))

rc=1.6d0
nrc=count(r(1:msh(1),1).le.rc)

  !
  ! calculate radial integration on atom site 
  ! <aephi|1/r^3|aephi>-<psphi|1/r^3|psphi>
  !

at_efg=0.d0
  do nt=1,ntyp
     do il1=1,paw_nbeta(nt)
        do il2=1,paw_nbeta(nt)
           work=0.d0
           do j = 2,nrc
           work(j)=(aephi(nt,il1)%psi(j)*aephi(nt,il2)%psi(j)-&
                psphi(nt,il1)%psi(j)*psphi(nt,il2)%psi(j))/r(j,nt)**3
           enddo
           work(1)=0.d0
           call simpson(nrc,work,rab(:,nt),at_efg(il1,il2,nt))
        enddo
     enddo
  enddo

!
!  calculation of the reconstruction part
!

  do ik = 1, nks
     call gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
     call davcio (evc, nwordwfc, iunwfc, ik, - 1)
     call init_paw_2 (npw, igk, xk (1, ik), paw_vkb)
     call ccalbec (paw_nkb, npwx, npw, nbnd, paw_becp, paw_vkb, evc)

do ibnd = 1, nbnd
   ijkb0 = 0
   do nt = 1, ntyp
      do na = 1, nat
         if (ityp (na) .eq.nt) then
            do ih = 1, paw_nh (nt)
               ikb = ijkb0 + ih
               nbs1=paw_indv(ih,nt)
               l1=paw_nhtol(ih,nt)
               m1=paw_nhtom(ih,nt)
               lm1=m1+l1**2
               do jh = 1, paw_nh (nt) 
                  jkb = ijkb0 + jh
                  nbs2=paw_indv(jh,nt)
                  l2=paw_nhtol(jh,nt)
                  m2=paw_nhtom(jh,nt)
                  lm2=m2+l2**2 
                  do lm=5,9
                     efg_corr(lm,na) =  efg_corr(lm,na) + &
                          (paw_becp(jkb,ibnd) * conjg(paw_becp(ikb,ibnd))) &
                          * at_efg(nbs1,nbs2,nt) * &
                          ap(lm,lm1,lm2)
                  enddo
               enddo
            enddo
            ijkb0 = ijkb0 + paw_nh (nt)
         endif
      enddo
   enddo
enddo
 
enddo 

!
!  transforme in cartesian coordinates
!

  efg_corr_tens(1,1,:)=real(sqrt(3.d0)*efg_corr(8,:) &
       - efg_corr(5,:))
  efg_corr_tens(2,2,:)=real(-sqrt(3.d0)*efg_corr(8,:)&
       - efg_corr(5,:))
  efg_corr_tens(3,3,:)=real(2.d0*efg_corr(5,:))
  efg_corr_tens(1,2,:)=real(-sqrt(3.d0)*efg_corr(9,:))
  efg_corr_tens(2,1,:)=efg_corr_tens(1,2,:)
  efg_corr_tens(1,3,:)=real(efg_corr(6,:)*sqrt(3.d0))
  efg_corr_tens(3,1,:)=efg_corr_tens(1,3,:)
  efg_corr_tens(2,3,:)=real(efg_corr(7,:)*sqrt(3.d0))
  efg_corr_tens(3,2,:)=efg_corr_tens(2,3,:)

  efg_corr_tens=-sqrt(4.d0*pi/5.d0)*efg_corr_tens


deallocate(work)
deallocate(efg_corr)
deallocate(at_efg)
deallocate(paw_vkb)
 
end subroutine efg_correction

