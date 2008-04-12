!--------------------------------------------------------------------
subroutine eff_pot (rho, nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, nl,   &
                     ngm, gg, gstart, nspin, alat, omega, ecutwfc,  &
                     charge, vstart, thresh_veff)
  !--------------------------------------------------------------------
  !
  !     Effective  potential (V_eff) in TF+vW scheme 
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : fpi, e2
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : prefix, nwordwfc, iunwfc
  USE klist,                ONLY : nelec
  USE gvect,                ONLY : nlm, g, qcutz, ecfixed, q2sigma
  USE wvfct,                ONLY : g2kin, wg, nbndx, et, nbnd, npwx, &
                                   igk, npw 
  USE uspp,                 ONLY : nkb
  USE scf,                  ONLY : v, vltot, vrs, rho_core
  USE gsmooth,              ONLY : nls, nlsm, nr1s, nr2s, nr3s, nrx1s,&
                                   nrx2s, nrx3s, nrxxs, doublegrid
  USE eff_v,                ONLY : rho_fft, veff
  USE mp_global,            ONLY : intra_pool_comm
  USE mp,                   ONLY : mp_sum

  USE control_flags,        ONLY : gamma_only
!  USE control_vdw,          ONLY : thresh_veff
!  USE wavefunctions_module, ONLY : evc, psic
  implicit none
  !
  !    input
  !
  integer :: nspin, nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, ngm, ngmw, &
             gstart, nl (ngm)             
  real(kind=DP) :: rho(nrxx, nspin), gg (ngm), alat, omega, ecutwfc, charge, &
                   charge_old
  !
  !    output
  !
  real(kind=DP), allocatable :: vv (:,:), rho_in(:,:), k_gamma(:) 
  !
  !    local variables
  !
  real(kind=DP) :: tpiba2, fac
  real(kind=DP), allocatable ::  aux (:,:), aux1 (:,:), psi (:,:), &
                                 psi_smooth(:,:), S(:)  
  complex(kind=DP), allocatable :: ws_psi(:,:), ws_hpsi(:,:), ws_psic(:)

  integer :: ir, is, ig
  real (kind=DP) :: avg1, avg2, eps
  integer :: nnn, nite, ite = 0
  real (kind=DP) :: vstart, thresh_veff
  character (len=10):: str_ite
  real(kind=DP), external :: erf
  !
  call start_clock('eff_pot')
  !
  allocate (aux(2,nrxx), aux1(2,ngm), psi(2,nrxx), psi_smooth(2,nrxx),S(nrxx) )
  allocate ( vv(nrxx, nspin) )
  allocate ( rho_in(nrxx, nspin) )
  !
  tpiba2 = (fpi / 2.d0 / alat) **2
  !
  rho_in = abs(rho)
  !
  !  set value for psi in real space
  !
  psi(2,:) = 0.d0
  psi(1,:) = sqrt( rho_in (:,1) )
  !
  !  bring psi to G space
  !
  call cft3 (psi, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1)
  !
  !  extract the smooth part of psi in G-space
  !
  psi_smooth(:,:) = 0.d0
  do ig = 1, ngm
    if ( (tpiba2 * gg(ig)) .lt. ( ecutwfc ) ) then 
      psi_smooth(1,nl(ig)) = psi(1,nl(ig))
      psi_smooth(2,nl(ig)) = psi(2,nl(ig))
    endif
  enddo
  !
  aux = psi_smooth
  psi = psi_smooth
  !
  !  bring psi_smooth to real space (approximation of psi)
  !
  call cft3 (psi_smooth, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
  !
  !  check the difference of the total charge density  
  !
  charge = 0.d0
  do is = 1, nspin
     do ir = 1, nrxx
        charge = charge + abs( rho_in(ir,is)-psi_smooth(is,ir)**2 )
     enddo
  enddo
  charge = charge * omega / (nr1*nr2*nr3) / nelec
#ifdef __PARA
  call mp_sum( charge, intra_pool_comm )
#endif
  WRITE( stdout, '(/,10x,"Charge difference due to FFT   ",f10.8)' ) charge
  !
  ! compute charge density using smooth wfc
  !
  do is = 1, nspin
     do ir = 1, nrxx
        rho_fft(ir,1) = psi_smooth(1,ir)**2
     enddo
  enddo
  !
  !  calculate P^2 |psi> in G-space (NB: V(G=0)=0 )
  !
  aux1(:,:) = 0.d0
  do ig = 1, ngm
     fac =  gg(ig) * tpiba2
     aux1(1,ig) = fac * aux(1,nl(ig))
     aux1(2,ig) = fac * aux(2,nl(ig))
  enddo
  !
  aux(:,:) = 0.d0
  do ig = 1, ngm
     aux(1,nl(ig)) = aux1(1,ig)
     aux(2,nl(ig)) = aux1(2,ig)
  enddo
  !
  if (gamma_only) then
     do ig = 1, ngm
        aux(1,nlm(ig)) =   aux1(1,ig)
        aux(2,nlm(ig)) = - aux1(2,ig)
     enddo
  end if
  !
  !      bring P^2 |psi>  to real space, kinetic term is kept 
  !      in aux
  !
  call cft3 (aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
  !
  !
  !  compute V_eff  potential by FT and then use it as initial
  !  potential for iteration  
  !
  avg1 = 0.d0
  avg2 = 0.d0
  nnn  = 0
  eps = 0.04d0
  is = 1
  if (.false.) then
     do ir = 1, nrxx
        vv(ir,is) = -aux(1,ir)
     end do
  else
     do ir = 1, nrxx
        if (abs(psi_smooth(1,ir)) > eps ) then
           vv(ir,is) = -aux(1,ir) / psi_smooth(1,ir)
        else
           avg1 = avg1 - aux(1,ir)
           avg2 = avg2 + psi_smooth(1,ir)
           nnn = nnn + 1
        end if
     enddo
#ifdef __PARA
     call mp_sum( avg1, intra_pool_comm )
     call mp_sum( avg2, intra_pool_comm )
     call mp_sum(nnn, intra_pool_comm) 
#endif
     if (nnn > 0 ) then
        do ir = 1, nrxx
           if (abs(psi_smooth(1,ir)) <= eps ) then
              vv(ir,is) = avg1 / avg2
           end if
        end do
     end if
  end if
  !
  ! uncomment the following loop will set local pot. as initial pot. 
  !
vstart=20000
  do ir = 1, nrxx
     vrs(ir,1) = v%of_r(ir,1) + vltot(ir)
     vv(ir,is) = erf(abs(psi_smooth(1,ir)*dble(vstart)))*vv(ir,is) + &
                (1.d0-erf( abs( psi_smooth(1,ir)*dble(vstart) ) ))*&
                vrs(ir,is)

  enddo
  !
  !  check the quality of trial potential
  !
  CALL check_v_eff(vv(1,1), charge)
  WRITE( stdout, '(/,10x,"Charge difference due to V_eff (initial)   ",f10.8,/)' ) charge
  !
  ! iterative procedure of V_eff if necessary 
  !
  nite = 0
  charge_old = 1000.d0
  do while ( charge .gt. thresh_veff )
     !
     call ite_veff_nhpsi( 1 )
     !
     !  check the quality of veff by solving equation
     !      p^2|phi> + (veff-mu)|phi> = 0
     !  whose the GS solution should be the square root 
     !  of the charge density
     !
     call check_v_eff(vv, charge)
     nite = nite + 1
!#ifdef __PARA
!     call ireduce(1, nite)
!#endif
     write( stdout, '(10x,"iter #   ", i3, "   charge diff.   ", f10.8, &
                     & "   thresh_veff   ", f10.8,/)' ) nite, charge, thresh_veff
     !
     if ( charge_old .lt. charge ) then
!        CALL io_pot( 1, TRIM( prefix )//'.veff', v, nspin )
!        CALL io_pot( 1, TRIM( prefix )//'.rho-coreff', rho_veff, nspin )
        write(stdout, '( 10x, 10("*"), "unstability happens", 10("*") )' )
!        goto 100
     endif
     !
     charge_old = charge
     !
  end do 
  !
  ! set the optmized eff. potential to veff
  !
100 continue
  veff(:,:) = vv(:,:)
  !
  deallocate ( vv, rho_in )
  deallocate (aux,aux1,psi,psi_smooth,S)
  !
  call stop_clock('eff_pot')
  !
  return
  !
  CONTAINS
     !
     !------------------------------------------------------------------------
     !
     subroutine ite_veff_nhpsi( nstep )
     !
     implicit none
     !
     integer :: nstep, nveff
     !
     real (kind=DP) :: alp, beta, s2r2, sr2, s2r, sr, r2, s2, D, Da, Db, w1
     !
     ! Compute S(r) at first step
     !
     call start_clock ('ite_veff')
     !
!write (stdout,*) ' enter ite_veff_nhpsi'
!CALL flush_unit( stdout )
     s2 = 0.d0
     do ir = 1, nrxx
        S(ir) = psi_smooth(1,ir) * ( aux(1,ir) + vv(ir,1)*psi_smooth(1,ir) )
        s2 = s2 + S(ir)**2
     enddo
#ifdef __PARA
     call mp_sum( s2, intra_pool_comm )
#endif
     !
     do nnn = 1, nstep
        !
        !
        ! Compute alpha & beta in Veff = Veff + alp*S(r) + beta
        !
        s2r2 = 0.d0
        s2r  = 0.d0
        sr2  = 0.d0
        sr   = 0.d0
        r2   = 0.d0
        !
        do ir = 1, nrxx
           r2   = r2   + psi_smooth(1,ir)**4
           s2r2 = s2r2 + ( S(ir) * psi_smooth(1,ir)**2 )**2
           sr2  = sr2  +   S(ir) * psi_smooth(1,ir)**4
           s2r  = s2r  + ( S(ir)**2) * psi_smooth(1,ir)**2
           sr   = sr   +   S(ir) * psi_smooth(1,ir)**2
        enddo
#ifdef __PARA
        call mp_sum( r2, intra_pool_comm )
        call mp_sum( s2r2, intra_pool_comm )
        call mp_sum( sr2, intra_pool_comm )
        call mp_sum( s2r, intra_pool_comm )
        call mp_sum( sr, intra_pool_comm )
#endif
        !
        D  = r2*s2r2 - sr2*sr2
        Da = sr*sr2  - s2r*r2
        Db = sr2*s2r - s2r2*sr
        !
        if (D.gt.0.d0) then
           alp = Da/D
           beta = Db/D
        else
           write(*,*) 'Det. of Hessian matrix is negative'
           stop
        endif
        !
!        if (mod(nnn,100) .eq. 0) then 
!           write(*,*)'iteration ',nnn
!           write(*,*) 's2 = ',s2
!           write(*,*) 'D = ' , D
!           write(*,*) 'Da = ', Da
!           write(*,*) 'Db = ', Db
!           write(*,*) 'alp = ',alp
!           write(*,*) 'beta = ',beta
!           write(*,*)
!        endif   
        !
        ! Update V-eff
        !  
        do ir = 1, nrxx
           vv(ir,1)= vv(ir,1) + alp*S(ir) + beta
           S(ir)   = S(ir) * (1.d0 + alp*psi_smooth(1,ir)**2) + & 
                     beta*psi_smooth(1,ir)**2
        enddo
        !
        s2 = 0.d0
        do ir = 1, nrxx
           s2 = s2 + S(ir)**2
        enddo
#ifdef __PARA
        call mp_sum( s2, intra_pool_comm )
#endif
        !
     enddo
     !
     call stop_clock ('ite_veff')
     !
     return
     !
     end subroutine ite_veff_nhpsi
     !
end subroutine eff_pot

