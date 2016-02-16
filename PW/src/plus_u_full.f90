! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Set of subroutines needed for full LDA+U calculations
! after Liechtenstein and co-workers (PRB 52, R5467 (1995)). 
! Works with two-component spinor WFs and with fully-relativistic
! pseudopotentials. 
! In the last case the WFs are projected onto: 
!   real spherical harmonics * 
!   averaged j=l+1/2, l-1/2 radial WFs *
!   up/down spinor.
! 
! A. Smogunov, C. Barreteau
!-----------------------------------------------------------------------

subroutine hubbard_matrix (lmax, L, U, J, u_matrix)
  !
  ! Build up the matrix of Coulomb integrals u_matrix(1,2,3,4)
  ! for real spherical harmonics. Implemented for s, p, d, f-shells. 
  ! Integrals with radial WFs are parametrized by U and J parameters.
  ! See Liechtenstein PRB 52, R5467 (1995), for example.  
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : rytoev, fpi 

  !
  implicit none
  !
  integer  :: lmax, L              ! max and actuel l

  real(DP), intent(in) :: U, J(3)  ! input parameters

  !  s: U
  !  p: U, J = J(1)
  !  d: U, J = J(1),  B  = J(2)
  !  f: U, J = J(1),  E2 = J(2), E3 = J(3)
 
  real(DP) :: u_matrix(2*lmax+1, 2*lmax+1, 2*lmax+1, 2*lmax+1), ak
  real(DP), allocatable :: ap(:,:,:), F(:)
  !
  integer :: n, nl, moffset, i, m1, m2, m3, m4, k, q 

!--
! number of all spher. harm.:

! from l = 0 to l = L
  nl = (L+1)**2
! from l = 0 to l = 2L
  n  = (2*L+1)**2
! up to L
  moffset = L**2 
!--

  allocate( ap(n,nl,nl) )
  allocate( F(0:6) )

!-- Set up the F_2k coefficients k = 0, 1, ... L 

  F(:) = 0.d0
  if (L.eq.0) then
   F(0) = U
  elseif (L.eq.1) then
   F(0) = U
   F(2) = 5.d0 * J(1)
  elseif (L.eq.2) then
   F(0) = U
   F(2) = 5.d0 * J(1) + 31.5d0 * J(2)
   F(4) = 9.d0 * J(1) - 31.5d0 * J(2)
  elseif (L.eq.3) then
   F(0) = U
   F(2) = 225.d0/54.d0*J(1)     + 32175.d0/42.d0*J(2)   + 2475.d0/42.d0*J(3)
   F(4) = 11.d0*J(1)            - 141570.d0/77.d0*J(2)  + 4356.d0/77.d0*J(3)
   F(6) = 7361.64d0/594.d0*J(1) + 36808.2d0/66.d0*J(2)  - 11154.d-2*J(3)
  else
   call errore( 'hubbard_matrix', &
                   & 'lda_plus_u is not implemented for L > 3 ...', 1 )
  endif
!--

  ap = 0.d0
  u_matrix = 0.d0

!-- Calculate Y_{kq} * Y_{lm} * Y_{lm'} integrals 
  call aainit_1(n, nl, ap) 
!--

  do m1 = 1, 2*l+1
    do m2 = 1, 2*l+1
      do m3 = 1, 2*l+1
        do m4 = 1, 2*l+1

         i = 0 
         do k = 0, 2*l, 2   

          ak = 0.d0
          do q = 1, 2*k + 1 
            i = i + 1
            ak = ak + ap(i,moffset+m1,moffset+m3) * ap(i,moffset+m2,moffset+m4) 
          enddo
          ak = ak * fpi / (2.d0*k+1.d0)
          u_matrix(m1,m2,m3,m4) = u_matrix(m1,m2,m3,m4) + ak*f(k) 
          i = i + 2*(k+1) + 1 
         enddo


        enddo
      enddo
    enddo
  enddo


  deallocate( ap )
  deallocate( f )

  return
end subroutine hubbard_matrix 

subroutine aainit_1(n2l, nl, ap)
    !-----------------------------------------------------------------------
    !
    ! this routine computes the expansion coefficients of 
    ! of two real spherical harmonics:
    !
    !     Y_limi(r) * Y_ljmj(r) = \sum_LM  ap(LM,limi,ljmj)  Y_LM(r)
    !
    !     using:
    !     ap(LM,limi,ljmj) = int Y_LM(r) * Y_limi(r) * Y_ljmj(r)  
    !
    !
    ! On output:
    ! ap     the expansion coefficients
    !
    ! The indices limi,ljmj and LM assume the order for real spherical
    ! harmonics given in routine ylmr2
    !
    ! The routine is similar to aainit in Modules/uspp.f90 
    !
    USE kinds,     ONLY : DP
    USE matrix_inversion

    implicit none
    !
    ! input: n2l = (2*L+1)**2,      nl = (L+1)**2    -   dimensions of 
    !                {2*L}      and      {L}             full spaces
    !  
    integer :: n2l, nl
    !
    ! local variables
    !
    integer :: li, lj, l, ir
    real(DP) , allocatable :: r(:,:), rr(:), ylm(:,:), mly(:,:)
    real(DP) :: ap(n2l, nl, nl), compute_ap_1

    allocate (r( 3, n2l ))    
    allocate (rr( n2l ))    
    allocate (ylm( n2l, n2l ))    
    allocate (mly( n2l, n2l ))    

    r(:,:)   = 0.d0
    ylm(:,:) = 0.d0
    mly(:,:) = 0.d0
    ap(:,:,:)= 0.d0

    ! - generate an array of random vectors (uniform deviate on unitary sphere)

    call gen_rndm_r_1 (n2l,r,rr)

    ! - generate the real spherical harmonics for the array: ylm(ir,lm)

    call ylmr2(n2l,n2l,r,rr,ylm)

    !-  store the inverse of ylm(ir,lm) in mly(lm,ir)

    call invmat(n2l, ylm, mly)

    !-  for each l,li,lj compute ap(l,li,lj) 
    do li = 1, nl
       do lj = 1,nl 
          do l = 1, n2l
             ap(l,li,lj) = 0.0_DP
             do ir = 1, n2l
               ap(l,li,lj) = ap(l,li,lj) + mly(l,ir)*ylm(ir,li)*ylm(ir,lj)
             end do
          end do
       end do
    end do
    
    deallocate(mly)
    deallocate(ylm)
    deallocate(rr)
    deallocate(r)
    
    return
end subroutine aainit_1


subroutine gen_rndm_r_1(llx,r,rr)
    !-----------------------------------------------------------------------
    ! - generate an array of random vectors (uniform deviate on unitary sphere)
    !
    USE kinds,     ONLY : DP
    USE constants,      ONLY: tpi
    USE random_numbers, ONLY: randy
    
    implicit none
    !
    ! first the I/O variables
    !
    integer :: llx         ! input: the dimension of r and rr
    
    real(DP) :: &
         r(3,llx),  &! output: an array of random vectors
         rr(llx)    ! output: the norm of r
    !
    ! here the local variables
    !
    integer :: ir
    real(DP) :: costheta, sintheta, phi
    
    do ir = 1, llx
       costheta = 2.0_DP * randy() - 1.0_DP
       sintheta = SQRT ( 1.0_DP - costheta*costheta)
       phi = tpi * randy()
       r (1,ir) = sintheta * cos(phi)
       r (2,ir) = sintheta * sin(phi)
       r (3,ir) = costheta
       rr(ir)   = 1.0_DP
    end do
    
    return
end subroutine gen_rndm_r_1

!-----------------------------------------------------------------------
subroutine comp_dspinldau () 
   !
   !  Initialize the spin rotation matrix d_spin_ldau for each symmetry operation. 
   !  Will be needed when symmetrizing the +U occupation matrix.  
   !

   USE kinds,     ONLY : DP
   USE ldaU,      ONLY : d_spin_ldau
   USE symm_base, ONLY : nsym, sr, t_rev, sname
   !
   implicit none

   complex(DP) :: a, b 

   integer :: isym

   d_spin_ldau = 0.d0

     do isym = 1, nsym

              call find_u(sr(1,1,isym),d_spin_ldau(1,1,isym))

!-- if time-reversal:  d_spin_ldau --> i sigma_y d_spin_ldau^* 
!      
              if (t_rev(isym)==1) then
                a = CONJG( d_spin_ldau(1,1,isym) )
                b = CONJG( d_spin_ldau(1,2,isym) )
                d_spin_ldau(1,1,isym) = CONJG( d_spin_ldau(2,1,isym) )
                d_spin_ldau(1,2,isym) = CONJG( d_spin_ldau(2,2,isym) )
                d_spin_ldau(2,1,isym) = -a
                d_spin_ldau(2,2,isym) = -b

              endif 
     enddo
!--

   return  
end subroutine comp_dspinldau   

SUBROUTINE atomic_wfc_nc_updown (ik, wfcatom)
  !-----------------------------------------------------------------------
  !
  ! For noncollinear case: builds up the superposition (for a k-point "ik") of 
  ! pure spin up or spin down atomic wavefunctions.
  ! 
  ! Based on atomic_wfc.f90

  USE kinds,      ONLY : DP
  USE constants,  ONLY : tpi, fpi, pi
  USE cell_base,  ONLY : tpiba
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp, tau
  USE basis,      ONLY : natomwfc
  USE gvect,      ONLY : mill, eigts1, eigts2, eigts3, g
  USE klist,      ONLY : xk, ngk, igk_k
  USE wvfct,      ONLY : npwx, nbnd
  USE us,         ONLY : tab_at, dq
  USE uspp_param, ONLY : upf
  USE noncollin_module, ONLY : noncolin, npol, angle1, angle2
  USE spin_orb,   ONLY : lspinorb, rot_ylm, fcoef, lmaxx, domag, &
                         starting_spin_angle
  !
  implicit none
  !
  integer, intent(in) :: ik
  complex(DP), intent(out) :: wfcatom (npwx, npol, natomwfc)
  !
  integer :: n_starting_wfc, lmax_wfc, nt, l, nb, na, m, lm, ig, iig, &
             i0, i1, i2, i3, nwfcm, npw
  real(DP), allocatable :: qg(:), ylm (:,:), chiq (:,:,:), gk (:,:)
  complex(DP), allocatable :: sk (:), aux(:)
  complex(DP) :: kphase
  real(DP) :: arg, px, ux, vx, wx

  call start_clock ('atomic_wfc')

  ! calculate max angular momentum required in wavefunctions
  lmax_wfc = 0
  do nt = 1, ntyp
     lmax_wfc = MAX ( lmax_wfc, MAXVAL (upf(nt)%lchi(1:upf(nt)%nwfc) ) )
  enddo
  !
  nwfcm = MAXVAL ( upf(1:ntyp)%nwfc )
  npw = ngk(ik)
  allocate ( ylm (npw,(lmax_wfc+1)**2), chiq(npw,nwfcm,ntyp), &
             sk(npw), gk(3,npw), qg(npw) )
  !
  do ig = 1, npw
     iig = igk_k(ig,ik)
     gk (1,ig) = xk(1, ik) + g(1,iig)
     gk (2,ig) = xk(2, ik) + g(2,iig)
     gk (3,ig) = xk(3, ik) + g(3,iig)
     qg(ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
  enddo
  !
  !  ylm = spherical harmonics
  !
  call ylmr2 ((lmax_wfc+1)**2, npw, gk, qg, ylm)
  !
  ! set now q=|k+G| in atomic units
  !
  do ig = 1, npw
     qg(ig) = sqrt(qg(ig))*tpiba
  enddo
  !
  n_starting_wfc = 0
  !
  ! chiq = radial fourier transform of atomic orbitals chi
  !
  do nt = 1, ntyp
     do nb = 1, upf(nt)%nwfc
        if ( upf(nt)%oc (nb) >= 0.d0) then
           do ig = 1, npw
              px = qg (ig) / dq - int (qg (ig) / dq)
              ux = 1.d0 - px
              vx = 2.d0 - px
              wx = 3.d0 - px
              i0 = INT( qg (ig) / dq ) + 1
              i1 = i0 + 1
              i2 = i0 + 2
              i3 = i0 + 3
              chiq (ig, nb, nt) = &
                     tab_at (i0, nb, nt) * ux * vx * wx / 6.d0 + &
                     tab_at (i1, nb, nt) * px * vx * wx / 2.d0 - &
                     tab_at (i2, nb, nt) * px * ux * wx / 2.d0 + &
                     tab_at (i3, nb, nt) * px * ux * vx / 6.d0
           enddo
        endif
     enddo
  enddo

  deallocate (qg, gk)
  allocate ( aux(npw) )
  !
  wfcatom(:,:,:) = (0.0_dp, 0.0_dp)
  !
  do na = 1, nat
     arg = (xk(1,ik)*tau(1,na) + xk(2,ik)*tau(2,na) + xk(3,ik)*tau(3,na)) * tpi
     kphase = CMPLX(cos (arg), - sin (arg) ,kind=DP)
     !
     !     sk is the structure factor
     !
     do ig = 1, npw
        iig = igk_k(ig,ik)
        sk (ig) = kphase * eigts1 (mill (1,iig), na) * &
                           eigts2 (mill (2,iig), na) * &
                           eigts3 (mill (3,iig), na)
     enddo
     !
     nt = ityp (na)
     do nb = 1, upf(nt)%nwfc
        if (upf(nt)%oc(nb) >= 0.d0) then
           l = upf(nt)%lchi(nb)
           !
           !
              IF ( upf(nt)%has_so ) THEN
                 !
                 call wfc_atom ( .true. )
                 !
              ELSE
                 !
                 call wfc_atom ( .false. )
                 !
              ENDIF
              !
        END IF
        !
     END DO
     !
  END DO

  if (n_starting_wfc /= natomwfc) call errore ('atomic_wfc_nc_updown', &
       'internal error: some wfcs were lost ', 1)

  deallocate(aux, sk, chiq, ylm)

  call stop_clock ('atomic_wfc')
  return

CONTAINS

   SUBROUTINE wfc_atom ( soc )
   !
   !
   real(DP) :: j
   real(DP), ALLOCATABLE :: chiaux(:)
   integer :: nc, ib
   logical :: soc ! .true. if the fully-relativistic pseudo
   !

!  If SOC go on only if j=l+1/2
   if (soc) j = upf(nt)%jchi(nb)
   if (soc.and.ABS(j-l+0.5_DP)<1.d-4 ) return
!

   allocate (chiaux(npw))

   if (soc) then 

!
!  Find the index for j=l-1/2
!
     if (l == 0)  then
        chiaux(:)=chiq(:,nb,nt)
     else
        do ib=1, upf(nt)%nwfc
           if ((upf(nt)%lchi(ib) == l).and. &
                        (ABS(upf(nt)%jchi(ib)-l+0.5_DP)<1.d-4)) then
              nc=ib
              exit
           endif
        enddo
!
!  Average the two radial functions 
!
        chiaux(:)=(chiq(:,nb,nt)*(l+1.0_DP)+chiq(:,nc,nt)*l)/(2.0_DP*l+1.0_DP)
     endif

   else

     chiaux(:) = chiq(:,nb,nt)

   endif


   do m = 1, 2 * l + 1
      lm = l**2 + m
      n_starting_wfc = n_starting_wfc + 1

      if (n_starting_wfc + 2*l+1 > natomwfc) call errore &
            ('atomic_wfc_nc', 'internal error: too many wfcs', 1)
      do ig=1,npw
         aux(ig) = sk(ig)*ylm(ig,lm)*chiaux(ig)
      enddo
!
      do ig=1,npw
!
         wfcatom(ig,1,n_starting_wfc) = aux(ig) 
         wfcatom(ig,2,n_starting_wfc) = 0.d0 
!
         wfcatom(ig,1,n_starting_wfc+2*l+1) = 0.d0
         wfcatom(ig,2,n_starting_wfc+2*l+1) = aux(ig) 
!
      enddo
   enddo
   n_starting_wfc = n_starting_wfc + 2*l+1

   deallocate (chiaux)
   !
   END SUBROUTINE wfc_atom
   !

END SUBROUTINE atomic_wfc_nc_updown  
