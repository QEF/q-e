program efg
  use kinds, only: DP
  use io_files, only: nd_nmbr,prefix
  use parameters, only: ntypx

  implicit none
  character (len=80) :: filerec(ntypx)
  real(kind=DP) :: Q(ntypx)
  integer :: ios

  namelist / inputpp / prefix, filerec, Q
  ! set default value

  Q=0.d0

  read (5, inputpp, err=200, iostat=ios)

200 call errore('efg.x', 'reading inputpp namelist', abs(ios))
  call start_postproc(nd_nmbr)
  call read_file
  call openfil
  call read_recon(filerec)

  call do_efg(Q) 

  call stop_pp
  stop
end program efg

subroutine read_recon(filerec)

  use read_pseudo_module
!  use wfc_utils
  use basis
!  use us
  use paw
  use atom
  use kinds, only: DP
  use parameters, only : ntypx
  USE io_global,  ONLY : stdout
  implicit none

  real(kind=DP), allocatable :: tmp(:)
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
  allocate( tmp(1:maxval(mesh(1:ntyp))) )


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
        read(14,*)  l
        call scan_end(14,'L')
        call scan_begin(14,'REC_AE',.false.)
        read(14,*) (tmp(j),j=1,kkphi)
        aephi(jtyp,i)%psi(1:kkphi)=tmp(1:kkphi)
        aephi(jtyp,i)%label%l=l
        call scan_end(14,'REC_AE')
!        call find_mt_radius(aephi(jtyp,i))
!        call save_psi(32,aephi(jtyp,i))
        allocate (psphi(jtyp,i)%psi(maxval(mesh(1:ntyp))))
        psphi(jtyp,i)%label%nt=jtyp
        psphi(jtyp,i)%label%n=i
        psphi(jtyp,i)%label%l=l
        psphi(jtyp,i)%kkpsi=kkphi
        call scan_begin(14,'REC_PS',.false.)
        read(14,*) (tmp(j),j=1,kkphi)
        psphi(jtyp,i)%psi(1:kkphi)= tmp(1:kkphi)
        call scan_end(14,'REC_PS')
!        call find_mt_radius(psphi(jtyp,i))
!        call save_psi(32,psphi(jtyp,i))
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
  use scf           !rho
  use gvect         !gvectors and parameters fot the FFT
  use brilz         !parameters of the cell
  use basis         !coordinates of the atoms
  use symme
  use pseud, only : zv !valence charge
  implicit none

  real(kind=DP) :: Q(ntypx), eta, Cq
  real(kind=DP) :: fac, trace, arg, e2
  integer :: alpha, beta, ig, na, i
!  real(kind=DP), allocatable:: aux(:,:)
  complex(kind=DP), allocatable:: aux(:)
  complex(kind=DP), allocatable:: efgg_el(:,:,:),efgr_el(:,:,:)
  complex(kind=DP), allocatable:: efg_io(:,:,:)
  real(kind=DP), allocatable:: zion(:), efg_corr_tens(:,:,:), efg(:,:,:)
  real(kind=DP):: efg_eig(3), v(3)
  complex(kind=DP) :: work(3,3), efg_vect(3,3)

!  allocate(aux(2,nrxx))
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
!  aux(2,:)=0.d0
!  aux(1,:)= rho(:,1)
  aux(:)= rho(:,1)
  efgg_el(:,:,:)=(0.d0,0.d0)
  
  call cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,-1)

!calculation of the electic field gradient in the G-space

  do ig= gstart, ngm
     trace = 1.d0/3.d0 * gg(ig)
     do alpha=1,3
        efgg_el(ig,alpha,alpha)= -trace
        do beta=1,3
           efgg_el(ig,alpha,beta)=(efgg_el(ig,alpha,beta) + &
                g(alpha,ig) * g(beta,ig)) &
                * fac * (aux(nl(ig)))/gg(ig)
!                * fac * CMPLX(aux(1,nl(ig)),aux(2,nl(ig)))/gg(ig) 
        enddo
     enddo
  enddo
!print *,trace, efgg_el(ig,1,1)
 
!fourier transform on the atomic position
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

!  zion(1)=6.0
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

  do na=1,nat
     do beta=1,3

        write (stdout,1000) atm(ityp(na)),na,"efg_corr", &
             (real(efg_corr_tens(alpha,beta,1)) , alpha =1,3 )
        
     enddo
     write (stdout,*)
  enddo

  do na=1,nat
     efg(:,:,na)=real(efg_corr_tens(:,:,na)+efgr_el(na,:,:)+ &
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
!  print *,work

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
1200 FORMAT(1x,a,1x,i3,5x,'Q= ',f6.3,5x,' Cq= ',f9.6,5x,' eta= ',f9.6)

end subroutine do_efg
 
subroutine efg_correction(efg_corr_tens)

  use io_files
  use kinds , only: dp
  use us
  use atom
!  use becmod
  use gvect
  use klist
!  use units
  use brilz
  use basis
  use wvfct
  use wavefunctions_module, only: evc
!  use wfc_utils
  use paw
  use constants, only: pi

  implicit none

  integer :: j, ill, nt, ibnd, il1, il2, ik, iat, nbs1, nbs2, kkpsi
  integer :: lm,l,m,m1,m2,lm1,lm2, l1, l2,m3,m4,n,n1,nrc
  integer :: ijkb0,ijkb,ih,jh,na,np,ikb,jkb
  real(kind=dp), allocatable :: at_efg(:,:,:), work(:)
  real(kind=dp) ,intent(out):: efg_corr_tens(3,3,nat)
  complex(kind=dp) , allocatable :: efg_corr(:,:), u(:,:)
  integer :: lmtonh(0:lmaxx-1,-lmaxx:lmaxx),mtom(lmaxx*2+1)
  real(kind=dp) :: sqrt2,rc

  allocate (efg_corr(lmaxx**2,nat))


  efg_corr=0.d0
  kkpsi=aephi(1,1)%kkpsi
  allocate (work(kkpsi))

  lm=0
  do l=0,lmaxx-1
     lm=lm+1
     lmtonh(l,0)=lm
     do m=1,l
        lm=lm+1
        lmtonh(l,m)=lm
        lm=lm+1
        lmtonh(l,-m)=lm
     enddo
  enddo
  mtom=0
  mtom(1)=0
  lm=1
  do m=1,lmaxx 
     lm=lm+1
     mtom(lm)=m
     lm=lm+1
     mtom(lm)=-m
  enddo
 
  call init_paw_1

  allocate (at_efg(paw_nkb,paw_nkb,ntypx)) 
  allocate (paw_vkb( npwx,  paw_nkb))
  allocate (paw_becp(paw_nkb, nbnd))

!        do j = 1, paw_nbeta (1)
!            do ih=1,msh(1)
!              write(59,*) r(ih,1),psphi(1,j)%psi(ih)
!           enddo
!           write(59,*)
!        enddo


rc=1.6d0
nrc=count(r(1:msh(1),1).le.rc)
!print *,'nrc',nrc, size(r(:,1)),r(nrc,1)
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
!           print *,work(1:20)
           call simpson(nrc,work,rab(:,nt),at_efg(il1,il2,nt))
!           print *,"at_efg",at_efg(il1,il2,nt), il1, il2
        enddo
     enddo
  enddo



  do ik = 1, nks
     call gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
     call davcio (evc, nwordwfc, iunwfc, ik, - 1)
     call init_paw_2 (npw, igk, xk (1, ik), paw_vkb)
     call ccalbec (paw_nkb, npwx, npw, nbnd, paw_becp, paw_vkb, evc)

     allocate(u(2*lmaxx+1, 2*lmaxx+1))
!
!  In the spin-orbit case we need the unitary matrix u which rotates the
!  real spherical harmonics and yields the complex ones.
!
     sqrt2=1.d0/dsqrt(2.d0)
     u=(0.d0,0.d0)
     l=lmaxx
     u(l+1,1)=(1.d0,0.d0)
     do n1=2,2*l+1,2
        m=n1/2
        n=l+1-m
        u(n,n1)=dcmplx((-1.d0)**m*sqrt2,0.d0)
        u(n,n1+1)=dcmplx(0.d0,-(-1.d0)**m*sqrt2)
        n=l+1+m
        u(n,n1)=dcmplx(sqrt2,0.d0)
        u(n,n1+1)=dcmplx(0.d0, sqrt2)
     enddo


!at_efg=1.d0
!paw_becp=1.d0
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
!                     do m3=-l1,l1
!                        do m4=-l2,l2
!                     efg_corr(lm,na)=efg_corr(lm,na)+ap(lm1,lm,lm2)
!                     print *,lm1,lm,lm2,ap(lm1,lm,lm2)
                     efg_corr(lm,na) =  efg_corr(lm,na) + &
                          (paw_becp(jkb,ibnd) * conjg(paw_becp(ikb,ibnd))) &
                          * at_efg(nbs1,nbs2,nt) * &
!                          (u(m3+l+1,m1)) *  conjg(u(m4+l+1,m2)) * &
                          ap(lm,lm1,lm2)
!                     enddo
!                     enddo
                  enddo
               enddo
            enddo
            ijkb0 = ijkb0 + nh (nt)
         endif
      enddo
   enddo
enddo
 
enddo 


!  print *,efg_corr(:,1)
 
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
