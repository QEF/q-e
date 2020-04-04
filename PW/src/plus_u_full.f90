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
!
!-----------------------------------------------------------------------
SUBROUTINE hubbard_matrix( lmax, L, U, J, u_matrix )
  !---------------------------------------------------------------------
  !! Builds up the matrix of Coulomb integrals \(\text{u_matrix}(1,2,3,4)\)
  !! for real spherical harmonics. Implemented for s, p, d, f-shells.  
  !! Integrals with radial WFs are parametrized by U and J parameters:
  !! \begin{equation}\notag
  !! \begin{split}
  !! s&: U & \\
  !! p&: U, J = J(1) & \\
  !! d&: U, J = J(1),\ B  &= J(2) \\
  !! f&: U, J = J(1),\ E2 &= J(2),\ E3 = J(3)
  !! \end{split}
  !! \end{equation}
  !! See, for example: Liechtenstein, PRB 52, R5467 (1995).
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : rytoev, fpi 
  !
  IMPLICIT NONE
  !
  INTEGER :: lmax
  !! max l
  INTEGER :: L
  !! actual l
  REAL(DP), INTENT(IN) :: U
  !! input parameter (see main comment)
  REAL(DP), INTENT(IN) :: J(3)
  !! input parameters (see main comment)
  REAL(DP) :: u_matrix(2*lmax+1, 2*lmax+1, 2*lmax+1, 2*lmax+1)
  !! matrix of Coulomb integrals
  !
  ! ... local variables
  !
  REAL(DP) :: ak
  REAL(DP), ALLOCATABLE :: ap(:,:,:), F(:)
  INTEGER :: n, nl, moffset, i, m1, m2, m3, m4, k, q 
  !
  ! ... number of all spher. harm.:
  !
  ! from l = 0 to l = L
  nl = (L+1)**2
  ! from l = 0 to l = 2L
  n = (2*L+1)**2
  ! up to L
  moffset = L**2 
  !
  ALLOCATE( ap(n,nl,nl) )
  ALLOCATE( F(0:6) )
  !
  ! ... set up the F_2k coefficients k = 0, 1, ... L 
  !
  F(:) = 0.d0
  IF (L==0) THEN
     F(0) = U
  ELSEIF (L==1) THEN
     F(0) = U
     F(2) = 5.d0 * J(1)
  ELSEIF (L==2) THEN
     F(0) = U
     F(2) = 5.d0 * J(1) + 31.5d0 * J(2)
     F(4) = 9.d0 * J(1) - 31.5d0 * J(2)
  ELSEIF (L==3) THEN
     F(0) = U
     F(2) = 225.d0/54.d0*J(1)     + 32175.d0/42.d0*J(2)   + 2475.d0/42.d0*J(3)
     F(4) = 11.d0*J(1)            - 141570.d0/77.d0*J(2)  + 4356.d0/77.d0*J(3)
     F(6) = 7361.64d0/594.d0*J(1) + 36808.2d0/66.d0*J(2)  - 11154.d-2*J(3)
  ELSE
     CALL errore( 'hubbard_matrix', &
                   & 'lda_plus_u is not implemented for L > 3 ...', 1 )
  ENDIF
  ! 
  ap = 0.d0
  u_matrix = 0.d0
  !
  ! ... Calculate Y_{kq} * Y_{lm} * Y_{lm'} integrals 
  CALL aainit_1( n, nl, ap )
  !
  DO m1 = 1, 2*l+1
    DO m2 = 1, 2*l+1
      DO m3 = 1, 2*l+1
        DO m4 = 1, 2*l+1
          !
          i = 0 
          DO k = 0, 2*l, 2   
            !
            ak = 0.d0
            DO q = 1, 2*k + 1 
              i = i + 1
              ak = ak + ap(i,moffset+m1,moffset+m3) * ap(i,moffset+m2,moffset+m4) 
            ENDDO
            ak = ak * fpi / (2.d0*k+1.d0)
            u_matrix(m1,m2,m3,m4) = u_matrix(m1,m2,m3,m4) + ak*f(k) 
            i = i + 2*(k+1) + 1 
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !
  DEALLOCATE( ap )
  DEALLOCATE( f )
  !
  RETURN
  !
END SUBROUTINE hubbard_matrix 
!
!
!----------------------------------------------------------------------------
SUBROUTINE aainit_1( n2l, nl, ap )
    !-----------------------------------------------------------------------
    !! This routine computes the expansion coefficients of two real 
    !! spherical harmonics:
    !
    !! \[ Y_{l_im_i}(r)\cdot Y_{l_jm_j}(r) = \sum_{LM} \text{ap}(LM,l_im_i,l_jm_j)\cdot Y_{LM}(r) \]
    !
    !! by using:
    !
    !! \[ \text{ap}(LM,l_im_i,l_jm_j) = \int Y_{LM}(r)\cdot Y_{l_im_i}(r)\cdot Y_{l_jm_j}(r)\ dr \]
    !
    !! The indices \(l_im_i\), \(l_jm_j\) and \(LM\) assume the order for real spherical
    !! harmonics given in routine \(\texttt{ylmr2}\).  
    !! The routine is similar to \(\texttt{aainit}\) in 'Modules/uspp.f90'.
    !
    USE kinds,           ONLY : DP
    USE matrix_inversion
    !
    IMPLICIT NONE
    !
    ! input: n2l = (2*L+1)**2,      nl = (L+1)**2 - dimensions of 
    !              {2*L}       and       {L}        full spaces
    !  
    INTEGER :: n2l
    !! \(\text{n2l} = (2L+1)^2\)
    INTEGER :: nl
    !! \(\text{nl} = (L+1)^2\)
    REAL(DP) :: ap(n2l,nl,nl)
    !! the expansion coefficients
    !
    ! ... local variables
    !
    INTEGER :: li, lj, l, ir
    REAL(DP), ALLOCATABLE :: r(:,:), rr(:), ylm(:,:), mly(:,:)
    REAL(DP) :: compute_ap_1
    !
    ALLOCATE( r(3,n2l) )
    ALLOCATE( rr(n2l)  )    
    ALLOCATE( ylm(n2l,n2l) )    
    ALLOCATE( mly(n2l,n2l) )    
    !
    r(:,:)   = 0.d0
    ylm(:,:) = 0.d0
    mly(:,:) = 0.d0
    ap(:,:,:)= 0.d0
    !
    ! - generate an array of random vectors (uniform deviate on unitary sphere)
    !
    CALL gen_rndm_r_1( n2l, r, rr )
    !
    ! - generate the real spherical harmonics for the array: ylm(ir,lm)
    !
    CALL ylmr2( n2l, n2l, r, rr, ylm )
    !
    !-  store the inverse of ylm(ir,lm) in mly(lm,ir)
    !
    CALL invmat( n2l, ylm, mly )
    !
    !-  for each l,li,lj compute ap(l,li,lj) 
    DO li = 1, nl
       DO lj = 1,nl 
          DO l = 1, n2l
             ap(l,li,lj) = 0.0_DP
             DO ir = 1, n2l
               ap(l,li,lj) = ap(l,li,lj) + mly(l,ir)*ylm(ir,li)*ylm(ir,lj)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !
    DEALLOCATE( mly )
    DEALLOCATE( ylm )
    DEALLOCATE( rr  )
    DEALLOCATE( r   )
    !
    RETURN
    !
END SUBROUTINE aainit_1
!
!----------------------------------------------------------------------------
SUBROUTINE gen_rndm_r_1( llx, r, rr )
    !-----------------------------------------------------------------------
    !! Generate an array of random vectors (uniform deviate on unitary sphere).
    !
    USE kinds,            ONLY : DP
    USE constants,        ONLY : tpi
    USE random_numbers,   ONLY : randy
    !
    IMPLICIT NONE
    !
    INTEGER :: llx
    !! input: the dimension of r and rr
    REAL(DP) :: r(3,llx)
    !! output: an array of random vectors
    REAL(DP) :: rr(llx)
    !! output: the norm of r
    !
    ! ... local variables
    !
    INTEGER :: ir
    REAL(DP) :: costheta, sintheta, phi
    !
    DO ir = 1, llx
       costheta = 2.0_DP * randy() - 1.0_DP
       sintheta = SQRT(1.0_DP - costheta*costheta)
       phi = tpi * randy()
       r(1,ir) = sintheta * COS(phi)
       r(2,ir) = sintheta * SIN(phi)
       r(3,ir) = costheta
       rr(ir)  = 1.0_DP
    ENDDO
    !
    RETURN
    !
END SUBROUTINE gen_rndm_r_1
!
!
!-----------------------------------------------------------------------
SUBROUTINE comp_dspinldau() 
   !----------------------------------------------------------------------
   !! Initializes the spin rotation matrix d_spin_ldau for each symmetry
   !! operation. It is needed when symmetrizing the +U occupation matrix.  
   !
   USE kinds,     ONLY : DP
   USE ldaU,      ONLY : d_spin_ldau
   USE symm_base, ONLY : nsym, sr, t_rev, sname
   !
   IMPLICIT NONE
   !
   ! ... local variables
   !
   COMPLEX(DP) :: a, b 
   INTEGER :: isym
   !
   d_spin_ldau = 0.d0
   !  
   DO isym = 1, nsym  
      !        
      CALL find_u( sr(1,1,isym), d_spin_ldau(1,1,isym) )  
      !        
      ! ... if time-reversal:  d_spin_ldau --> i sigma_y d_spin_ldau^*         
      !              
      IF (t_rev(isym)==1) THEN        
        a = CONJG( d_spin_ldau(1,1,isym) )        
        b = CONJG( d_spin_ldau(1,2,isym) )        
        d_spin_ldau(1,1,isym) = CONJG( d_spin_ldau(2,1,isym) )        
        d_spin_ldau(1,2,isym) = CONJG( d_spin_ldau(2,2,isym) )        
        d_spin_ldau(2,1,isym) = -a        
        d_spin_ldau(2,2,isym) = -b        
      ENDIF         
      !        
   ENDDO  
   !  
   RETURN
   !
END SUBROUTINE comp_dspinldau   
!
!
!--------------------------------------------------------------------------
SUBROUTINE atomic_wfc_nc_updown( ik, wfcatom )
  !-----------------------------------------------------------------------
  !! For noncollinear case: builds up the superposition (for a k-point 
  !! \(\text{ik}\)) of pure spin-up or spin-down atomic wavefunctions.
  ! 
  !! Based on 'atomic_wfc.f90'
  !
  USE kinds,             ONLY : DP
  USE constants,         ONLY : tpi, fpi, pi
  USE cell_base,         ONLY : tpiba
  USE ions_base,         ONLY : nat, ntyp => nsp, ityp, tau
  USE basis,             ONLY : natomwfc
  USE gvect,             ONLY : mill, eigts1, eigts2, eigts3, g
  USE klist,             ONLY : xk, ngk, igk_k
  USE wvfct,             ONLY : npwx, nbnd
  USE us,                ONLY : tab_at, dq
  USE uspp_param,        ONLY : upf
  USE noncollin_module,  ONLY : noncolin, npol, angle1, angle2
  USE spin_orb,          ONLY : lspinorb, rot_ylm, fcoef, lmaxx, domag, &
                                starting_spin_angle
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik
  !! the k-point index
  COMPLEX(DP), INTENT(OUT) :: wfcatom(npwx,npol,natomwfc)
  !! the superposition of atomic wavefunctions (up or down)
  !
  ! ... local variables
  !
  INTEGER :: n_starting_wfc, lmax_wfc, nt, l, nb, na, m, lm, ig, iig, &
             i0, i1, i2, i3, nwfcm, npw
  REAL(DP), ALLOCATABLE :: qg(:), ylm(:,:), chiq(:,:,:), gk(:,:)
  COMPLEX(DP), ALLOCATABLE :: sk(:), aux(:)
  COMPLEX(DP) :: kphase
  REAL(DP) :: arg, px, ux, vx, wx
  !
  CALL start_clock( 'atomic_wfc' )
  !
  ! ... calculate max angular momentum required in wavefunctions
  !
  lmax_wfc = 0
  do nt = 1, ntyp
     lmax_wfc = MAX( lmax_wfc, MAXVAL(upf(nt)%lchi(1:upf(nt)%nwfc)) )
  ENDDO
  !
  nwfcm = MAXVAL( upf(1:ntyp)%nwfc )
  npw = ngk(ik)
  ALLOCATE ( ylm (npw,(lmax_wfc+1)**2), chiq(npw,nwfcm,ntyp), &
             sk(npw), gk(3,npw), qg(npw) )
  !
  DO ig = 1, npw
     iig = igk_k(ig,ik)
     gk(1,ig) = xk(1,ik) + g(1,iig)
     gk(2,ig) = xk(2,ik) + g(2,iig)
     gk(3,ig) = xk(3,ik) + g(3,iig)
     qg(ig) = gk(1,ig)**2 +  gk(2,ig)**2 + gk(3,ig)**2
  ENDDO
  !
  ! ... ylm = spherical harmonics
  !
  CALL ylmr2( (lmax_wfc+1)**2, npw, gk, qg, ylm )
  !
  ! ... set now q=|k+G| in atomic units
  !
  DO ig = 1, npw
     qg(ig) = SQRT(qg(ig))*tpiba
  ENDDO
  !
  n_starting_wfc = 0
  !
  ! ... chiq = radial fourier transform of atomic orbitals chi
  !
  DO nt = 1, ntyp
     DO nb = 1, upf(nt)%nwfc
        IF ( upf(nt)%oc (nb) >= 0.d0) THEN
           DO ig = 1, npw
              px = qg (ig) / dq - INT(qg (ig) / dq)
              ux = 1.d0 - px
              vx = 2.d0 - px
              wx = 3.d0 - px
              i0 = INT( qg (ig) / dq ) + 1
              i1 = i0 + 1
              i2 = i0 + 2
              i3 = i0 + 3
              chiq(ig, nb, nt) = &
                     tab_at (i0, nb, nt) * ux * vx * wx / 6.d0 + &
                     tab_at (i1, nb, nt) * px * vx * wx / 2.d0 - &
                     tab_at (i2, nb, nt) * px * ux * wx / 2.d0 + &
                     tab_at (i3, nb, nt) * px * ux * vx / 6.d0
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  !
  DEALLOCATE( qg, gk )
  ALLOCATE( aux(npw) )
  !
  wfcatom(:,:,:) = (0.0_dp, 0.0_dp)
  !
  DO na = 1, nat
     arg = (xk(1,ik)*tau(1,na) + xk(2,ik)*tau(2,na) + xk(3,ik)*tau(3,na)) * tpi
     kphase = CMPLX(COS(arg), - SIN(arg), KIND=DP)
     !
     ! ... sk is the structure factor
     !
     DO ig = 1, npw
        iig = igk_k(ig,ik)
        sk(ig) = kphase * eigts1 (mill(1,iig), na) * &
                          eigts2 (mill(2,iig), na) * &
                          eigts3 (mill(3,iig), na)
     ENDDO
     !
     nt = ityp (na)
     DO nb = 1, upf(nt)%nwfc
        IF (upf(nt)%oc(nb) >= 0.d0) THEN
           l = upf(nt)%lchi(nb)
           !
           IF ( upf(nt)%has_so ) THEN
              !
              CALL wfc_atom( .TRUE. )
              !
           ELSE
              !
              CALL wfc_atom( .FALSE. )
              !
           ENDIF
           !
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  IF (n_starting_wfc /= natomwfc) CALL errore( 'atomic_wfc_nc_updown', &
                             'internal error: some wfcs were lost ', 1 )
  !
  DEALLOCATE( aux, sk, chiq, ylm )
  !
  CALL stop_clock ( 'atomic_wfc' )
  !
  RETURN
  !
CONTAINS
   !
   !--------------------------
   SUBROUTINE wfc_atom( soc )
      !---------------------------
      !
      LOGICAL :: soc
      !! .TRUE. if the fully-relativistic pseudo
      !
      ! ... local variables
      !
      REAL(DP) :: j
      REAL(DP), ALLOCATABLE :: chiaux(:)
      INTEGER :: nc, ib
      !
      ! ... If SOC go on only if j=l+1/2
      IF (soc) j = upf(nt)%jchi(nb)
      IF (soc .AND. ABS(j-l+0.5_DP)<1.d-4 ) RETURN
      !
      ALLOCATE( chiaux(npw) )
      !
      IF (soc) THEN 
        !
        ! ... Find the index for j=l-1/2
        !
        IF (l == 0)  THEN
           chiaux(:)=chiq(:,nb,nt)
        ELSE
           DO ib=1, upf(nt)%nwfc
              IF ((upf(nt)%lchi(ib) == l) .AND. &
                           (ABS(upf(nt)%jchi(ib)-l+0.5_DP)<1.d-4)) THEN
                 nc=ib
                 exit
              ENDIF
           ENDDO
           !
           ! ... Average the two radial functions 
           !
           chiaux(:) = (chiq(:,nb,nt)*(l+1.0_DP)+chiq(:,nc,nt)*l)/(2.0_DP*l+1.0_DP)
        ENDIF
        !
      ELSE
        !
        chiaux(:) = chiq(:,nb,nt)
        !
      ENDIF
      !
      DO m = 1, 2*l+1
         lm = l**2 + m
         n_starting_wfc = n_starting_wfc + 1

         IF (n_starting_wfc + 2*l+1 > natomwfc) CALL errore &
               ('atomic_wfc_nc', 'internal error: too many wfcs', 1)
         DO ig = 1, npw
            aux(ig) = sk(ig)*ylm(ig,lm)*chiaux(ig)
         ENDDO
         ! 
         DO ig = 1, npw
            !
            wfcatom(ig,1,n_starting_wfc) = aux(ig)
            wfcatom(ig,2,n_starting_wfc) = 0.d0
            !
            wfcatom(ig,1,n_starting_wfc+2*l+1) = 0.d0
            wfcatom(ig,2,n_starting_wfc+2*l+1) = aux(ig)
            !
         ENDDO
      ENDDO
      !
      n_starting_wfc = n_starting_wfc + 2*l+1
      !
      DEALLOCATE( chiaux )
      !
   END SUBROUTINE wfc_atom
   !
   !
END SUBROUTINE atomic_wfc_nc_updown  
