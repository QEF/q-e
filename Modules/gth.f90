
! Copyright (C) 2015 Sebastiano Caravati
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .!
!
! Tools for Goedecker-Teter-Hutter pseudopotentials
! Contains routine "atmlength", copyright by ABINIT group
!
module m_gth
  use kinds, only: dp
  implicit none
  !
  private
  public :: gth_parameters, readgth, vloc_gth, dvloc_gth, setlocq_gth, &
       mk_ffnl_gth, mk_dffnl_gth, deallocate_gth
  !
  type gth_parameters
     integer  :: itype, lloc, lmax
     real(dp) :: rloc, cc(4)
     integer,  pointer :: lll(:), ipr(:)
     real(dp), pointer :: rrl(:)
  end type gth_parameters
  type (gth_parameters), pointer, dimension(:), private, save :: gth_p
  !
contains
  !-----------------------------------------------------------------------
  subroutine mk_ffnl_gth(itype, ibeta, nq, qg, vq)
    !-----------------------------------------------------------------------
    !
    USE kinds,        ONLY: dp
    USE constants,    ONLY: pi, fpi, e2
    USE cell_base,    ONLY: omega
    
    implicit none
    !  
    ! I/O 
    integer,  intent(in)  :: itype, ibeta, nq
    real(dp), intent(in)  :: qg(nq)
    real(dp), intent(out) :: vq(nq)
    !     
    ! Local variables
    integer, parameter :: nprj_max(0:3)=[3, 3, 2, 1]
    integer  :: ii, my_gth, ll, iproj
    real(dp) :: rrl, qr2, fact
    !
    my_gth=0
    do ii=1,size(gth_p)
       if (gth_p(ii)%itype==itype) then
          my_gth=ii
          exit
       endif
    enddo
    if (my_gth==0) call errore('mk_ffnl_gth', 'cannot map itype in some gtp param. set', itype)
    iproj=gth_p(my_gth)%ipr(ibeta)
    ll=gth_p(my_gth)%lll(ibeta)
    rrl=gth_p(my_gth)%rrl(ll)
    if ( ll<0 .or. ll>3  ) call errore('mk_ffnl_gth', 'wrong l:', ll)
    if ( iproj>nprj_max(ll) ) call errore('mk_ffnl_gth', 'projector exceeds max. n. of projectors', iproj)
    !
    lif: if (ll==0) then     ! s channel
       !
       if(iproj==1)then
          do ii=1,nq
             qr2=(qg(ii)*rrl)**2
             vq(ii)=exp(-0.5_dp*qr2)
          end do
       else if(iproj==2)then
          do ii=1,nq
             qr2=(qg(ii)*rrl)**2
             vq(ii)=2._dp/sqrt(15._dp) * exp(-0.5_dp*qr2) * ( 3._dp-qr2 )
          end do
       else if(iproj==3)then
          do ii=1,nq
             qr2=(qg(ii)*rrl)**2
             vq(ii)=(4._dp/3._dp)/sqrt(105._dp) * exp(-0.5_dp*qr2) * &
                  &         (15._dp-10._dp*qr2 + qr2**2)
          end do
       end if
       !
    else if (ll==1) then lif ! p channel
       !
       if(iproj==1)then
          do ii=1,nq
             qr2=(qg(ii)*rrl)**2
             vq(ii)=(1._dp/sqrt(3._dp)) * exp(-0.5_dp*qr2) * qg(ii)
          end do
       else if(iproj==2)then
          do ii=1,nq
             qr2=(qg(ii)*rrl)**2
             vq(ii)=(2._dp/sqrt(105._dp)) * exp(-0.5_dp*qr2) * qg(ii)*(5._dp-qr2)
          end do
       else if(iproj==3)then
          do ii=1,nq
             qr2=(qg(ii)*rrl)**2
             vq(ii)=(4._dp/3._dp)/sqrt(1155._dp) * exp(-0.5_dp*qr2) * &
                  &         qg(ii) * (35._dp-14._dp*qr2+qr2**2)
          end do
       end if
       !
    else if (ll==2) then lif ! d channel [ ONLY 2 PROJECTORS!! ]
       !
       if(iproj==1)then
          do ii=1,nq
             qr2=(qg(ii)*rrl)**2
             vq(ii)=(1._dp/sqrt(15._dp)) * exp(-0.5_dp*qr2) * qg(ii)**2
          end do
       else if(iproj==2)then
          do ii=1,nq
             qr2=(qg(ii)*rrl)**2
             vq(ii)=(2._dp/3._dp)/sqrt(105._dp) * exp(-0.5_dp*qr2) * &
                  &         qg(ii)**2 * (7._dp-qr2)
          end do
       end if
       !
    else if (ll==3) then lif ! f channel [ ONLY 1 PROJECTOR!! ]
       !
       do ii=1,nq
          qr2=(qg(ii)*rrl)**2
          vq(ii)=qg(ii)**3 * exp(-0.5_dp*qr2)
       end do
       !
    end if lif
    !
    fact = e2 * fpi * pi**0.25_dp * sqrt( 2._dp**(ll+1) * rrl**(2*ll+3) / omega )
    vq(:)=fact*vq(:) 
    !
  end subroutine mk_ffnl_gth
  !-----------------------------------------------------------------------
  subroutine mk_dffnl_gth(itype, ibeta, nq, qg, dvq)
    !-----------------------------------------------------------------------
    !
    USE kinds,        ONLY: dp
    USE constants,    ONLY: pi, fpi, e2
    USE cell_base,    ONLY: omega, tpiba

    implicit none
    !  
    ! I/O 
    integer,  intent(in)  :: itype, ibeta, nq
    real(dp), intent(in)  :: qg(nq)
    real(dp), intent(out) :: dvq(nq)
    !     
    ! Local variables
    integer, parameter :: nprj_max(0:3)=[3, 3, 2, 1]
    integer  :: ii, my_gth, ll, iproj
    real(dp) :: rrl, rl2, qt, q1r2, q3r4, q5r6, qr2, qr4, qr6, fact, e_qr2_h
    !
    my_gth=0
    do ii=1,size(gth_p)
       if (gth_p(ii)%itype==itype) then
          my_gth=ii
          exit
       endif
    enddo
    if (my_gth==0) call errore('mk_dffnl_gth', 'cannot map itype in some gtp param. set', itype)
    iproj=gth_p(my_gth)%ipr(ibeta)
    ll=gth_p(my_gth)%lll(ibeta)
    rrl=gth_p(my_gth)%rrl(ll)
    if ( ll<0 .or. ll>3  ) call errore('mk_dffnl_gth', 'wrong l:', ll)
    if ( iproj>nprj_max(ll) ) call errore('mk_dffnl_gth', 'projector exceeds max. n. of projectors', iproj)
    !
    lif: if (ll==0) then     ! s channel
       !
       if(iproj==1)then
          do ii=1,nq
             qt = sqrt(qg(ii))*tpiba
             rl2= rrl**2
             qr2= qt*qt*rl2
             e_qr2_h=exp(-0.5_dp*qr2)
             !
             dvq(ii)=-qt*rl2*e_qr2_h
          end do
       else if(iproj==2)then
          do ii=1,nq
             qt = sqrt(qg(ii))*tpiba
             rl2= rrl**2
             q1r2=qt*rl2
             qr2= qt*q1r2
             q3r4=qr2*q1r2
             e_qr2_h=exp(-0.5_dp*qr2)
             !
             dvq(ii)=2._dp/sqrt(15._dp) * e_qr2_h * (-5._dp*q1r2+q3r4)
          end do
       else if(iproj==3)then
          do ii=1,nq
             qt = sqrt(qg(ii))*tpiba
             rl2= rrl**2
             q1r2=qt*rl2
             qr2= qt*q1r2
             q3r4=qr2*q1r2
             q5r6=qr2*q3r4
             e_qr2_h=exp(-0.5_dp*qr2)
             !
             dvq(ii)=(4._dp/3._dp)/sqrt(105._dp) * e_qr2_h * (-35._dp*q1r2 + 14._dp*q3r4 - q5r6)
          end do
       end if
       !
    else if (ll==1) then lif ! p channel
       !
       if(iproj==1)then
          do ii=1,nq
             qt = sqrt(qg(ii))*tpiba
             qr2=(qt*rrl)**2
             e_qr2_h=exp(-0.5_dp*qr2)
             !
             dvq(ii)=(1._dp/sqrt(3._dp)) * e_qr2_h * (1._dp - qr2)
          end do
       else if(iproj==2)then
          do ii=1,nq
             qt = sqrt(qg(ii))*tpiba
             qr2=(qt*rrl)**2
             qr4=qr2**2
             e_qr2_h=exp(-0.5_dp*qr2)
             !
             dvq(ii)=(2._dp/sqrt(105._dp)) * e_qr2_h * (5._dp - 8._dp*qr2 + qr4)
          end do
       else if(iproj==3)then
          do ii=1,nq
             qt = sqrt(qg(ii))*tpiba
             qr2=(qt*rrl)**2
             qr4=qr2**2
             qr6=qr4*qr2
             e_qr2_h=exp(-0.5_dp*qr2)
             !
             dvq(ii)=(4._dp/3._dp)/sqrt(1155._dp) * e_qr2_h * (35._dp - 77._dp*qr2 + 19._dp*qr4 - qr6)
          end do
       end if
       !
    else if (ll==2) then lif ! d channel [ ONLY 2 PROJECTORS!! ]
       !
       if(iproj==1)then
          do ii=1,nq
             qt = sqrt(qg(ii))*tpiba
             qr2=(qt*rrl)**2
             e_qr2_h=exp(-0.5_dp*qr2)
             !
             dvq(ii)=(1._dp/sqrt(15._dp)) * e_qr2_h * qt*(2._dp - qr2)
          end do
       else if(iproj==2)then
          do ii=1,nq
             qt = sqrt(qg(ii))*tpiba
             qr2=(qt*rrl)**2
             qr4=qr2**2
             e_qr2_h=exp(-0.5_dp*qr2)
             !
             dvq(ii)=(2._dp/3._dp)/sqrt(105._dp) * e_qr2_h * qt*(14._dp - 11._dp*qr2 + qr4)
          end do
       end if
       !
    else if (ll==3) then lif ! f channel [ ONLY 1 PROJECTOR!! ]
       !
       do ii=1,nq
          qt = qg(ii)*tpiba**2
          qr2=qt*rrl**2
          e_qr2_h=exp(-0.5_dp*qr2)
          !
          dvq(ii)=e_qr2_h * qt*(3._dp - qr2)
       end do
       !
    end if lif
    !
    fact = e2 * fpi * pi**0.25_dp * sqrt( 2._dp**(ll+1) * rrl**(2*ll+3) / omega )
    dvq(:)=fact*dvq(:) 
    !
  end subroutine mk_dffnl_gth
!-----------------------------------------------------------------------
subroutine vloc_gth(itype, zion, tpiba2, ngl, gl, omega, vloc)
  !-----------------------------------------------------------------------
  !
  USE kinds,        ONLY: dp
  USE constants,    ONLY: pi, fpi, e2, eps8

  implicit none
  !
  ! I/O
  integer,  intent(in)  :: itype, ngl
  real(dp), intent(in)  :: zion, tpiba2, omega, gl (ngl)
  real(dp), intent(out) :: vloc (ngl)
  !
  ! Local variables
  integer  :: ii, my_gth, igl, igl0
  real(dp) :: cc1, cc2, cc3, cc4, rloc, epsatm, gx, gx2, rq2, rl3, e_rq2h, fact
  !
  ! Find gtp param. set for type itype
  my_gth=0
  do ii=1,size(gth_p)
    if (gth_p(ii)%itype==itype) then
      my_gth=ii
      exit
    endif
  enddo
  if (my_gth==0) call errore('vloc_gth', 'cannot map itype in some gth param. set', itype)
  rloc=gth_p(my_gth)%rloc
  cc1=gth_p(my_gth)%cc(1)
  cc2=gth_p(my_gth)%cc(2)
  cc3=gth_p(my_gth)%cc(3)
  cc4=gth_p(my_gth)%cc(4)

  ! Compute epsatm = lim(q->0) [Vloc(q) + zion/(Pi*q^2)]
  epsatm=2._dp*pi*rloc**2*zion+(2._dp*pi)**(1.5_dp)*rloc**3*(cc1+3._dp*cc2+15._dp*cc3+105._dp*cc4)
  ! 1/\Omega * \sum_i epsatm(i) is v_loc(G=0)

  ! Compute vloc(q)
  if (gl (1) < eps8) then
     !
     ! first the G=0 term
     !
!    vloc (1) = 0._dp
     vloc (1) = epsatm
     igl0 = 2
  else
     igl0 = 1
  endif
  !
  !   here the G<>0 terms, we first compute the part of the integrand 
  !   function independent of |G| in real space
  !
  do igl = igl0, ngl
     gx     = sqrt (gl (igl) * tpiba2)
     gx2    = gx**2
     rq2    = (gx*rloc)**2
     rl3    = rloc**3
     e_rq2h = exp(-0.5_dp*rq2)
     vloc (igl) = &
         fpi * e_rq2h*(-zion/gx2 + sqrt(pi/2._dp)*rl3* &
           ( &
             cc1 + &
             cc2*(3._dp-rq2) + &
             cc3*(15._dp-10._dp*rq2+rq2**2) + &
             cc4*(105._dp-rq2*(105._dp-rq2*(21._dp-rq2))) &
           ) &
        )
  enddo
  !
  fact = e2 / omega
  vloc (:) = vloc(:) * fact
  !
end subroutine vloc_gth
!-----------------------------------------------------------------------
subroutine dvloc_gth(itype, zion, tpiba2, ngl, gl, omega, dvloc)
  !-----------------------------------------------------------------------
  !
  ! dvloc = D Vloc (g^2) / D g^2 = (1/2g) * D Vloc(g) / D g
  !
  USE kinds,        ONLY: dp
  USE constants,    ONLY: pi, tpi, e2, eps8

  implicit none
  !
  ! I/O
  integer,  intent(in)  :: itype, ngl
  real(dp), intent(in)  :: zion, tpiba2, omega, gl (ngl)
  real(dp), intent(out) :: dvloc (ngl)
  !
  ! Local variables
  integer  :: ii, my_gth, igl, igl0
  real(dp) :: cc1, cc2, cc3, cc4, rloc, &
              gx, gx2, gx3, rl2, rl3, rq2, r2q, r4g3, r6g5, e_rq2h, fact
  !
! IF ( do_comp_esm ) call errore('vloc_gth', 'ESM not implemented', itype)
  !
  ! Find gtp param. set for type itype
  my_gth=0
  do ii=1,size(gth_p)
    if (gth_p(ii)%itype==itype) then
      my_gth=ii
      exit
    endif
  enddo
  if (my_gth==0) call errore('dvloc_gth', 'cannot map itype in some gtp param. set', itype)
  rloc=gth_p(my_gth)%rloc
  cc1=gth_p(my_gth)%cc(1)
  cc2=gth_p(my_gth)%cc(2)
  cc3=gth_p(my_gth)%cc(3)
  cc4=gth_p(my_gth)%cc(4)

  ! Compute vloc(q)
  if (gl (1) < eps8) then
     !
     ! first the G=0 term
     !
     dvloc (1) = 0._dp
     igl0 = 2
  else
     igl0 = 1
  endif
  !
  !   here the G<>0 terms, we first compute the part of the integrand 
  !   function independent of |G| in real space
  !
  do igl = igl0, ngl
     gx     = sqrt (gl (igl) * tpiba2)
     gx2    = gx**2
     gx3    = gx*gx2
     rl2    = rloc**2
     rl3    = rloc*rl2
     rq2    = gx2*rl2
     r2q    = gx*rl2
     r4g3   = rl2*rl2*gx3
     r6g5   = r4g3*rl2*gx2
     e_rq2h = exp(-0.5_dp*rq2)
     dvloc (igl) = &
         e_rq2h*(zion*(rq2+2._dp)/gx3 + sqrt(pi/2._dp)*rl3* &
           ( &
             ( &
               - 2._dp*r2q* (cc2+10._dp*cc3+105._dp*cc4) &
               + 4._dp*r4g3*(cc3+21._dp*cc4) &
               - 6._dp*r6g5* cc4 &
             ) - r2q*( &
               cc1 + &
               cc2*(3._dp-rq2) + &
               cc3*(15._dp-10._dp*rq2+rq2**2) + &
               cc4*(105._dp-rq2*(105._dp-rq2*(21._dp-rq2))) &
             ) &
           ) &
        )/gx
  enddo
  !
  fact = tpi * e2 / omega
  dvloc (:) = dvloc(:) * fact
  !
end subroutine dvloc_gth
!-----------------------------------------------------------------------
subroutine setlocq_gth(itype, xq, zion, tpiba2, ngm, g, omega, vloc)
!----------------------------------------------------------------------
  !
  USE kinds,        ONLY: dp
  USE constants,    ONLY: pi, fpi, e2, eps8

  implicit none
  !  
  ! I/O
  integer,  intent(in)  :: itype, ngm
  real(dp), intent(in)  :: xq (3), zion, tpiba2, omega, g(3,ngm)
  real(dp), intent(out) :: vloc (ngm)
  !  
  ! Local variables
  integer  :: ii, ig, my_gth
  real(dp) :: cc1, cc2, cc3, cc4, rloc, g2a, gx, gx2, rq2, rl3, e_rq2h, fact
  !
  ! Find gtp param. set for type itype
  my_gth=0
  do ii=1,size(gth_p)
    if (gth_p(ii)%itype==itype) then
      my_gth=ii
      exit
    endif
  enddo
  if (my_gth==0) call errore('vloc_gth', 'cannot map itype in some gth param. set', itype)
  rloc=gth_p(my_gth)%rloc
  cc1=gth_p(my_gth)%cc(1)
  cc2=gth_p(my_gth)%cc(2)
  cc3=gth_p(my_gth)%cc(3)
  cc4=gth_p(my_gth)%cc(4)
  !
  do ig = 1, ngm
    g2a = (xq (1) + g (1, ig) ) **2 + &
          (xq (2) + g (2, ig) ) **2 + &
          (xq (3) + g (3, ig) ) **2
    if (g2a < eps8) then
      vloc (ig) = 0.d0
    else
      gx     = sqrt (g2a * tpiba2)
      gx2    = gx**2
      rq2    = (gx*rloc)**2
      rl3    = rloc**3
      e_rq2h = exp(-0.5_dp*rq2)
      vloc (ig) = &
         fpi * e_rq2h*(-zion/gx2 + sqrt(pi/2._dp)*rl3* &
           ( &
             cc1 + &
             cc2*(3._dp-rq2) + &
             cc3*(15._dp-10._dp*rq2+rq2**2) + &
             cc4*(105._dp-rq2*(105._dp-rq2*(21._dp-rq2))) &
           ) &
        )
    endif
  enddo
  !
  fact = e2 / omega
  vloc (:) = vloc(:) * fact
  !
end subroutine setlocq_gth
!-----------------------------------------------------------------------
subroutine deallocate_gth( lflag )
  !-----------------------------------------------------------------------
  !

  implicit none
  !  
  ! I/O 
  logical, intent(in) :: lflag
  !
  ! Local variables
  integer :: ii
  !
  IF ( lflag .and. ASSOCIATED( gth_p ) ) THEN
     DO ii=1, SIZE(gth_p)
        DEALLOCATE ( gth_p(ii)%lll, gth_p(ii)%ipr, gth_p(ii)%rrl )
     ENDDO
     DEALLOCATE( gth_p )
  ENDIF
  !
end subroutine deallocate_gth
!-----------------------------------------------------------------------
subroutine readgth (iunps, np, upf)
  !-----------------------------------------------------------------------
  !
  USE kinds,        ONLY: dp
  USE constants,    ONLY: e2, tpi
  USE parameters,   ONLY: lmaxx
  USE funct,        ONLY: set_dft_from_name, dft_is_hybrid
  USE pseudo_types, ONLY: pseudo_upf

  implicit none
  !
  ! I/O
  TYPE (pseudo_upf) :: upf
  integer :: iunps, np
  !
  ! Local variables
  integer  :: ios, pspdat, pspcod, pspxc, lmax, lloc, mmax, ii, jj, ll, nn, nnonloc, &
              nprl, os, ns, iv, jv
  real(dp) :: rcore, qcore, rc2, prefact, znucl, r2well, rloc, rrl, cc(4)
  character(len=256)            :: info
  character(len=  1), parameter :: ch10=char(10), spdf(0:3) = ['S','P','D','F']
  character(len=  2), external  :: atom_name
  integer,  allocatable         :: nproj(:)
  real(dp), allocatable         :: hij(:,:,:), kij(:,:,:)
  type(gth_parameters), pointer, dimension(:) :: gth_tmp_p
  !
  os=0; if (associated(gth_p)) os=size(gth_p)
  ns=os+1; allocate(gth_tmp_p(ns))
  if (os>0) then
    gth_tmp_p(1:os)=gth_p(1:os)
    deallocate(gth_p)
  end if
  nullify(gth_p)
  gth_p=>gth_tmp_p
  nullify(gth_tmp_p)
  gth_p(ns)%itype=np
  !
  upf%is_gth=.true.
  upf%generated="GTH norm-conserving PP, generated by Matthias Krack"
  upf%author   ="Goedecker/Hartwigsen/Hutter/Teter/Krack"
  upf%comment  ="GTH analytical, separable"
  upf%tvanp=.false.; upf%tpawp=.false.; upf%nlcc=.false.; upf%tcoulombp=.false.; upf%has_so=.false.
  upf%rel = 'scalar'; upf%typ = 'NC'; upf%lmax_rho = 0; upf%nwfc=0; upf%nqf = 0; upf%nqlc= 0; upf%kkbeta=-1
  upf%etotps =0._dp; upf%ecutrho=0._dp; upf%ecutwfc=0._dp
  allocate(upf%rcut(upf%nbeta), upf%rcutus(upf%nbeta), upf%lchi(upf%nwfc))
  upf%rcut(:) = 0._dp
  upf%rcutus(:) = 0._dp

  read (iunps, '(a)', end=400, err=400, iostat=ios) info
  read (iunps, *, err=400) znucl, upf%zp, pspdat
  if (upf%zp <= 0._dp .or. upf%zp > 100 ) call errore ('readgth', 'Wrong zp ', np)
  upf%psd=atom_name ( NINT(znucl) )
  call gth_grid_for_rho(upf,znucl)

  read (iunps, *, err=400) pspcod,pspxc,lmax,lloc,mmax,r2well
  IF ( pspcod /= 10 .AND. pspcod /= 12 ) &
     call errore ('readgth', 'unknown/invalid pspcod:', pspcod )
  IF ( pspcod == 12 ) THEN
     ! pseudo with NLCC
     upf%nlcc=.true.
     upf%generated="New Soft-Accurate NLCC pseudopotentials, generated by Santanu Saha"
     upf%author=upf%author//"/Saha"
  ENDIF
  IF ( lmax-1 > lmaxx ) call errore ('readgth', 'strange lmax', lmax-1)
  IF ( lmax == lloc) THEN
     upf%lmax = lmax-1
  ELSE
     upf%lmax = lmax
  ENDIF
  upf%lloc = lloc
  ! write(6, '(2f10.5,2x,i8,t47,a)' )  znucl,upf%zp,pspdat,'znucl, zion, pspdat'
  ! write(6, '(4i5,i10,f10.5,t47,a)' ) pspcod,pspxc,lmax,lloc,mmax,r2well,&
  !  'pspcod,pspxc,lmax,lloc,mmax,r2well'
  gth_p(ns)%lloc=lloc; gth_p(ns)%lmax=lmax

  IF (pspxc == 1) THEN
     upf%dft = 'PZ'
  ELSE IF (pspxc == 7) THEN
     upf%dft = 'PW'
  ELSE IF (pspxc == 11) THEN
     upf%dft = 'PBE'
  ELSE IF (pspxc == 18) THEN
     upf%dft = 'BLYP'
  ELSE IF (pspxc == -101130) THEN ! PBE from libXC
     upf%dft = 'PBE'
  ELSE
     call errore ('readgth', 'pspxc cod. cannot be understood', abs (np) )
  ENDIF
  call set_dft_from_name( upf%dft )
  !
  cc(:)=0._dp
  read (iunps, *, err=400) rloc,nn,(cc(jj),jj=1,nn)
  gth_p(ns)%rloc =rloc
  gth_p(ns)%cc(:)=cc(:)
  ! write(6, '(a,f12.7)' ) ' rloc=',rloc
  ! write(6, '(a,i1,a,4f12.7)' ) ' cc(1:',nn,')=',(cc(jj),jj=1,nn)
  read (iunps, *, err=400) nnonloc
  allocate(hij(0:lmax,3,3),kij(0:lmax,3,3),nproj(lmax+1),gth_p(ns)%rrl(0:lmax))
  hij(:,:,:)=0._dp; kij(:,:,:)=0._dp
  !
  ! Read and echo the coefficients of non-local projectors
  upf%nbeta=0
  prjloop: do ll=0,lmax
    read (iunps, *, err=400) rrl,nprl,(hij(ll,1,jj),jj=1,nprl)
    upf%nbeta = upf%nbeta + nprl
    gth_p(ns)%rrl(ll)=rrl
    do ii=2,nprl
      read (iunps, *, err=400) (hij(ll,ii,jj),jj=ii,nprl)
    end do
    nproj(ll+1)=nprl
    !   write(6, '(a,i3,a,f12.7,2a,3f12.7,2a,12x,2f12.7,2a,24x,f12.7)' )&
    !&   ' for angular momentum l =',ll,' r(l) =',rrl,ch10,&
    !&   '   h11, h12, h13 =', (hij(ll,1,jj),jj=1,3),ch10,&
    !&   '        h22, h23 =', (hij(ll,2,jj),jj=2,3),ch10,&
    !&   '             h33 =', (hij(ll,3,jj),jj=3,3)
    if (ll==0) cycle
    do ii=1,nprl
      read (iunps, *, err=400) (kij(ll,ii,jj),jj=ii,nprl)
    end do
    !   write(6, '(a,3f12.7,2a,12x,2f12.7,2a,24x,f12.7)' )&
    !&   '   k11, k12, k13 =', (kij(ll,1,jj),jj=1,3),ch10,&
    !&   '        k22, k23 =', (kij(ll,2,jj),jj=2,3),ch10,&
    !&   '             k33 =', (kij(ll,3,jj),jj=3,3)
  end do prjloop
  !
  if (upf%nlcc) then
    read (iunps, *, err=400) rcore, qcore
    ALLOCATE ( upf%rho_atc(upf%mesh) )
    rc2 = rcore**2
    prefact = qcore * (znucl-upf%zp) / (sqrt(tpi)*rcore)**3
    do ii=1,upf%mesh
      upf%rho_atc(ii) = prefact * exp(-0.5_dp * upf%r(ii)**2 / rc2)
    enddo
  end if
  !
  allocate(upf%lll(upf%nbeta), upf%els_beta(upf%nbeta), upf%dion(upf%nbeta,upf%nbeta))
  allocate(gth_p(ns)%lll(upf%nbeta), gth_p(ns)%ipr(upf%nbeta))
  iv=0
  lloop: do ll=0,lmax
    nprl=nproj(ll+1)
    iloop: do ii=1,nprl
      iv = iv+1
      gth_p(ns)%lll(iv)=ll
      gth_p(ns)%ipr(iv)=ii
      upf%lll(iv)=ll; WRITE (upf%els_beta(iv), '(I1,A1)' ) ii, spdf(ll)
      jloop: do jj=ii, nprl
        jv = iv+jj-ii
        upf%dion(iv,jv) = hij(ll,ii,jj)/e2
        if ( jj > ii ) upf%dion(jv,iv) = upf%dion(iv,jv)
      enddo jloop
    enddo iloop
  enddo lloop
  !
  deallocate(hij,kij,nproj)
  return
  !
400 call errore ('readgth', 'pseudo file is empty or wrong', abs (np) )
end subroutine readgth
!*****************************************************************************************
subroutine gth_grid_for_rho(upf,znucl)
  USE kinds,        ONLY: dp
  USE constants,    ONLY: pi, fpi
  USE pseudo_types, ONLY: pseudo_upf
  implicit none
  ! I/O
  TYPE (pseudo_upf) :: upf
  real(dp) :: znucl
  ! Local
  integer  :: i, mesh
  real(dp) :: xmin, amesh, rmax, rr, rab, length, two_l2, znorml
  !
  xmin = -7.0d0
  amesh=0.0125d0
  rmax =100.0d0
  mesh = 1 + (log(znucl*rmax)-xmin)/amesh
  mesh = (mesh/2)*2+1 ! mesh is odd (Simpson!)
  !
  call atmlength(0._dp,length,upf%zp,znucl)
  two_l2=2._dp*length**2
  znorml=fpi*upf%zp/(pi*two_l2)**1.5_dp
  ALLOCATE (upf%r(mesh), upf%rab(mesh), upf%rho_at(mesh))
  DO i=1, mesh
     rr = exp (xmin+dble(i-1)*amesh)/znucl
     rab= rr*amesh
     upf%r(i)   = rr
     upf%rab(i) = rab
     ! actually 4*pi*r**2*rho(r) !!!
     upf%rho_at(i)=znorml*exp(-rr**2/two_l2)*rr**2
  END DO
  upf%mesh  =mesh
  upf%xmin  =xmin
  upf%rmax  =rmax
  upf%zmesh =znucl
  upf%dx    =amesh
  !
end subroutine gth_grid_for_rho
!*****************************************************************************************
!{\src2tex{textfont=tt}}
!!****f* ABINIT/atmlength
!! NAME
!! atmlength
!!
!! FUNCTION
!! Return atomic decay length for one given type of atom.
!! This length is used to generate an approximate atomic gaussian density
!! in reciprocal space:   n^AT(G)=exp[-(2pi.length.G)^2]
!!
!! COPYRIGHT
!! Copyright (C) 2000-2010 ABINIT group (DCA,XG,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! densty=parameter for initialisation of the density of this atom type
!!        if densty>0, returned decay length if densty !
!! zion=charge on current type of atom (real number)
!! znucl=atomic number, for current type of atom
!!
!! OUTPUT
!! length=decay lenth
!!
!! PARENTS
!!      initro
!!
!! CHILDREN
!!
!! SOURCE

subroutine atmlength(densty,length,zion,znucl)

 use kinds,     only: dp
 use constants, only: tol10=>eps8

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: densty,zion,znucl
 real(dp),intent(out) :: length

!Local variables-------------------------------
!scalars
 integer :: nval
 real(dp) :: coreel
!arrays
 real(dp) :: data_length(16)

! *************************************************************************

!Either use the input value, or the default value, tabulated now.
 if(abs(densty)>tol10)then
   length=densty
 else

!  Count the number of core electrons.
   coreel=znucl-zion
!  Round the number of valence electrons
   nval=nint(zion)

!  For each set of core electron numbers, there are different decay lengths,
!  they start from nval=1, and proceed by group of 5, until a default is used

   if (nval==0) then
     length=0._dp

!    Bare ions : adjusted on 1h and 2he only
   else if(coreel<0.5)then
     data_length(1:4)=(/ .6_dp,.4_dp,.3_dp,.25_dp /)
     length=.2_dp
     if(nval<=4)length=data_length(nval)

!    1s2 core : adjusted on 3li, 6c, 7n, and 8o
   else if(coreel<2.5)then
     data_length(1:8)=(/ 1.8_dp,1.4_dp,1.0_dp ,.7_dp,.6_dp,&
&     .5_dp, .4_dp, .35_dp /)
     length=.3_dp
     if(nval<=8)length=data_length(nval)

!    Ne core (1s2 2s2 2p6) : adjusted on 11na, 13al, 14si and 17cl
   else if(coreel<10.5)then
     data_length(1:10)=(/ 2.0_dp,1.6_dp,1.25_dp,1.1_dp,1.0_dp,&
&     .9_dp, .8_dp, .7_dp , .7_dp, .7_dp  /)
     length=.6_dp
     if(nval<=10)length=data_length(nval)

!    Mg core (1s2 2s2 2p6 3s2) : adjusted on 19k, and on coreel==10
   else if(coreel<12.5)then
     data_length(1:10)=(/ 1.9_dp,1.5_dp,1.15_dp,1.0_dp,0.9_dp,&
&     .8_dp, .7_dp, .6_dp , .6_dp, .6_dp  /)
     length=.5_dp
     if(nval<=10)length=data_length(nval)

!    Ar core (Ne + 3s2 3p6) : adjusted on 20ca, 25mn and 30zn
   else if(coreel<18.5)then
     data_length(1:12)=(/ 2.0_dp ,1.8_dp ,1.5_dp,1.2_dp ,1.0_dp,&
&     .9_dp , .85_dp, .8_dp, .75_dp, .7_dp,&
&     .65_dp, .65_dp /)
     length=.6_dp
     if(nval<=12)length=data_length(nval)

!    Full 3rd shell core (Ar + 3d10) : adjusted on 31ga, 34se and 38sr
   else if(coreel<28.5)then
     data_length(1:14)=(/ 1.5_dp ,1.25_dp,1.15_dp,1.05_dp,1.00_dp,&
&     .95_dp, .95_dp, .9_dp , .9_dp , .85_dp,&
&     .85_dp, .80_dp, .8_dp , .75_dp         /)
     length=.7_dp
     if(nval<=14)length=data_length(nval)

!    Krypton core (Ar + 3d10 4s2 4p6) : adjusted on 39y, 42mo and 48cd
   else if(coreel<36.5)then
     data_length(1:12)=(/ 2.0_dp ,2.00_dp,1.60_dp,1.40_dp,1.25_dp,&
&     1.10_dp,1.00_dp, .95_dp, .90_dp, .85_dp,&
&     .80_dp, .75_dp /)
     length=.7_dp
     if(nval<=12)length=data_length(nval)

!    For the remaining elements, consider a function of nval only
   else
     data_length(1:12)=(/ 2.0_dp ,2.00_dp,1.55_dp,1.25_dp,1.15_dp,&
&     1.10_dp,1.05_dp,1.0_dp , .95_dp , .9_dp,&
&     .85_dp, .85_dp /)
     length=.8_dp
     if(nval<=12)length=data_length(nval)

   end if

!  End the choice between default and no-default
 end if

end subroutine atmlength
!!***
end module m_gth
