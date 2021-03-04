!
! Copyright (C) 2004-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE uspp_param
  !
  ! ... Ultrasoft and Norm-Conserving pseudopotential parameters
  !  
  USE pseudo_types, ONLY : pseudo_upf
  IMPLICIT NONE
  SAVE
  !
  INTEGER :: nsp = 0 
  TYPE (pseudo_upf),  ALLOCATABLE, TARGET :: upf(:)
  !! the upf structure contains all info on atomic pseudopotential parameters
  INTEGER, ALLOCATABLE :: nh(:)
  !! number of beta functions, with angular parts, per atomic type 
  INTEGER :: nhm
  !! max number of beta functions, including angular parts, across atoms
  INTEGER ::nbetam
  !! max number of radial beta functions
  INTEGER ::nwfcm
  !! max number of radial atomic wavefunctions across atoms
  INTEGER :: lmaxkb
  !! max angular momentum of beta functions
  INTEGER :: lmaxq
  !! max angular momentum + 1 for Q functions
  !
CONTAINS
  !
  SUBROUTINE init_uspp_dims ()
    !
    !!     calculates the number of beta functions for each atomic type
    !
    IMPLICIT NONE
    !
    INTEGER :: nt, nb
    !
    ! Check is needed, may be called more than once (but it shouldn't be!)
    ! Maybe nh should be allocated when upf is, when upf is read ?
    !
    IF ( .NOT. ALLOCATED(nh) ) ALLOCATE ( nh(nsp) )
    !
    lmaxkb = - 1
    DO nt = 1, nsp
       !
       nh (nt) = 0
       !
       ! do not add any beta projector if pseudo in 1/r fmt (AF)
       !
       IF ( upf(nt)%tcoulombp ) CYCLE 
       !
       DO nb = 1, upf(nt)%nbeta
          nh (nt) = nh (nt) + 2 * upf(nt)%lll(nb) + 1
          lmaxkb = MAX (lmaxkb, upf(nt)%lll(nb) )
       ENDDO
       !
    ENDDO
    lmaxq = 2*lmaxkb+1
    !
    ! calculate max numbers of beta functions and of atomic wavefunctions
    !
    nhm = MAXVAL (nh (1:nsp))
    nbetam = MAXVAL (upf(1:nsp)%nbeta)
    nwfcm = MAXVAL (upf(1:nsp)%nwfc)
    !
  END SUBROUTINE init_uspp_dims
  !
END MODULE uspp_param
!
!
MODULE uspp
  !
  ! Ultrasoft PPs:
  ! - Clebsch-Gordan coefficients "ap", auxiliary variables "lpx", "lpl"
  ! - beta and q functions of the solid
  !
  USE upf_kinds,   ONLY: DP
  USE upf_params,  ONLY: lmaxx, lqmax
  USE upf_spinorb, ONLY: fcoef, fcoef_d 
  IMPLICIT NONE
  PRIVATE
  SAVE
  !
  PUBLIC :: nlx, lpx, lpl, ap, aainit, indv, nhtol, nhtolm, indv_ijkb0, &
            nkb, nkbus, vkb, dvan, deeq, qq_at, qq_nt, nhtoj, ijtoh, beta, &
            becsum, ebecsum
  PUBLIC :: lpx_d, lpl_d, ap_d, indv_d, nhtol_d, nhtolm_d, indv_ijkb0_d, &
            vkb_d, dvan_d, deeq_d, qq_at_d, qq_nt_d, nhtoj_d, ijtoh_d, &
            becsum_d, ebecsum_d
  PUBLIC :: okvan, nlcc_any
  PUBLIC :: qq_so,   dvan_so,   deeq_nc,   fcoef 
  PUBLIC :: qq_so_d, dvan_so_d, deeq_nc_d, fcoef_d 
  PUBLIC :: dbeta
  !
  PUBLIC :: allocate_uspp, deallocate_uspp
  !
  ! GPU sync
  !
  PUBLIC  :: using_vkb, using_vkb_d
  LOGICAL :: vkb_d_ood = .false.    ! used to flag out of date variables
  LOGICAL :: vkb_ood   = .false.    ! used to flag out of date variables
  !
  ! Vars
  !
  INTEGER, PARAMETER :: &
       nlx  = (lmaxx+1)**2, &! maximum number of combined angular momentum
       mx   = 2*lqmax-1      ! maximum magnetic angular momentum of Q
  !
  INTEGER ::             &! for each pair of combined momenta lm(1),lm(2): 
       lpx(nlx,nlx),     &! maximum combined angular momentum LM
       lpl(nlx,nlx,mx)    ! list of combined angular momenta  LM
  REAL(DP) :: ap(lqmax*lqmax,nlx,nlx)
                          ! Clebsch-Gordan coefficients for spherical harmonics
  ! GPU vars
  INTEGER, ALLOCATABLE ::  & ! for each pair of combined momenta lm(1),lm(2): 
       lpx_d(:,:),         & ! maximum combined angular momentum LM
       lpl_d(:,:,:)          ! list of combined angular momenta  LM
  REAL(DP), ALLOCATABLE :: ap_d(:,:,:)
#if defined (__CUDA)
  attributes(DEVICE) :: lpx_d, lpl_d, ap_d
#endif
  !
  !
  INTEGER :: nkb,        &! total number of beta functions, with struct.fact.
             nkbus        ! as above, for US-PP only
  !
  INTEGER, ALLOCATABLE ::&
       indv(:,:),        &! indes linking  atomic beta's to beta's in the solid
       nhtol(:,:),       &! correspondence n <-> angular momentum l
       nhtolm(:,:),      &! correspondence n <-> combined lm index for (l,m)
       ijtoh(:,:,:),     &! correspondence beta indexes ih,jh -> composite index ijh
       indv_ijkb0(:)      ! first beta (index in the solid) for each atom 
  !
  ! GPU vars
  !
  INTEGER, ALLOCATABLE :: indv_d(:,:)
  INTEGER, ALLOCATABLE :: nhtol_d(:,:)
  INTEGER, ALLOCATABLE :: nhtolm_d(:,:)
  INTEGER, ALLOCATABLE :: ijtoh_d(:,:,:)
  INTEGER, ALLOCATABLE :: indv_ijkb0_d(:)
#if defined (__CUDA)
  attributes(DEVICE) :: indv_d, nhtol_d, nhtolm_d, ijtoh_d, indv_ijkb0_d
#endif

  LOGICAL :: &
       okvan = .FALSE.,&  ! if .TRUE. at least one pseudo is Vanderbilt
       nlcc_any=.FALSE.   ! if .TRUE. at least one pseudo has core corrections
  !
  COMPLEX(DP), ALLOCATABLE, TARGET :: &
       vkb(:,:)                ! all beta functions in reciprocal space
  REAL(DP), ALLOCATABLE :: &
       becsum(:,:,:)           ! \sum_i f(i) <psi(i)|beta_l><beta_m|psi(i)>
  REAL(DP), ALLOCATABLE :: &
       ebecsum(:,:,:)          ! \sum_i f(i) et(i) <psi(i)|beta_l><beta_m|psi(i)>
  REAL(DP), ALLOCATABLE :: &
       dvan(:,:,:),           &! the D functions of the solid
       deeq(:,:,:,:),         &! the integral of V_eff and Q_{nm} 
       qq_nt(:,:,:),          &! the integral of q functions in the solid (ONE PER NTYP) used to be the qq array
       qq_at(:,:,:),          &! the integral of q functions in the solid (ONE PER ATOM !!!!) 
       nhtoj(:,:)              ! correspondence n <-> total angular momentum
  !
  COMPLEX(DP), ALLOCATABLE :: & ! variables for spin-orbit/noncolinear case:
       qq_so(:,:,:,:),           &! Q_{nm}
       dvan_so(:,:,:,:),         &! D_{nm}
       deeq_nc(:,:,:,:)           ! \int V_{eff}(r) Q_{nm}(r) dr 
  !
  ! GPU vars
  !
  COMPLEX(DP), ALLOCATABLE :: vkb_d(:,:)
  REAL(DP),    ALLOCATABLE :: becsum_d(:,:,:)
  REAL(DP),    ALLOCATABLE :: ebecsum_d(:,:,:)
  REAL(DP),    ALLOCATABLE :: dvan_d(:,:,:)
  REAL(DP),    ALLOCATABLE :: deeq_d(:,:,:,:)
  REAL(DP),    ALLOCATABLE :: qq_nt_d(:,:,:)
  REAL(DP),    ALLOCATABLE :: qq_at_d(:,:,:)
  REAL(DP),    ALLOCATABLE :: nhtoj_d(:,:)
  COMPLEX(DP), ALLOCATABLE :: qq_so_d(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: dvan_so_d(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: deeq_nc_d(:,:,:,:)
#if defined(__CUDA)
  attributes (DEVICE) :: vkb_d, becsum_d, ebecsum_d, dvan_d, deeq_d, qq_nt_d, &
                         qq_at_d, nhtoj_d, qq_so_d, dvan_so_d, deeq_nc_d
#endif

  !
  ! spin-orbit coupling: qq and dvan are complex, qq has additional spin index
  ! noncolinear magnetism: deeq is complex (even in absence of spin-orbit)
  !
  REAL(DP), ALLOCATABLE :: &
       beta(:,:,:)           ! beta functions for CP (without struct.factor)
  REAL(DP), ALLOCATABLE :: &
       dbeta(:,:,:,:,:)      ! derivative of beta functions w.r.t. cell for CP (without struct.factor)
  !
CONTAINS
  !
  !-----------------------------------------------------------------------
  subroutine aainit(lli)
    !-----------------------------------------------------------------------
    !
    ! this routine computes the coefficients of the expansion of the product
    ! of two real spherical harmonics into real spherical harmonics.
    !
    !     Y_limi(r) * Y_ljmj(r) = \sum_LM  ap(LM,limi,ljmj)  Y_LM(r)
    !
    ! On output:
    ! ap     the expansion coefficients
    ! lpx    for each input limi,ljmj is the number of LM in the sum
    ! lpl    for each input limi,ljmj points to the allowed LM
    !
    ! The indices limi,ljmj and LM assume the order for real spherical
    ! harmonics given in routine ylmr2
    !
    USE upf_invmat
    implicit none
    !
    ! input: the maximum li considered
    !  
    integer :: lli
    !
    ! local variables
    !
    integer :: llx, l, li, lj
    real(DP) , allocatable :: r(:,:), rr(:), ylm(:,:), mly(:,:)
    ! an array of random vectors: r(3,llx)
    ! the norm of r: rr(llx)
    ! the real spherical harmonics for array r: ylm(llx,llx)
    ! the inverse of ylm considered as a matrix: mly(llx,llx)
    !
    if (lli < 0) call upf_error('aainit','lli not allowed',lli)

    if (lli*lli > nlx) call upf_error('aainit','nlx is too small ',lli*lli)

    llx = (2*lli-1)**2
    if (2*lli-1 > lqmax) &
         call upf_error('aainit','ap leading dimension is too small',llx)

    allocate (r( 3, llx ))    
    allocate (rr( llx ))    
    allocate (ylm( llx, llx ))    
    allocate (mly( llx, llx ))    

    r(:,:)   = 0.0_DP
    ylm(:,:) = 0.0_DP
    mly(:,:) = 0.0_DP
    ap(:,:,:)= 0.0_DP

    ! - generate an array of random vectors (uniform deviate on unitary sphere)

    call gen_rndm_r(llx,r,rr)

    ! - generate the real spherical harmonics for the array: ylm(ir,lm)

    call ylmr2(llx,llx,r,rr,ylm)

    !-  store the inverse of ylm(ir,lm) in mly(lm,ir)

    call invmat(llx, ylm, mly)

    !-  for each li,lj compute ap(l,li,lj) and the indices, lpx and lpl
    do li = 1, lli*lli
       do lj = 1, lli*lli
          lpx(li,lj)=0
          do l = 1, llx
             ap(l,li,lj) = compute_ap(l,li,lj,llx,ylm,mly)
             if (abs(ap(l,li,lj)) > 1.d-3) then
                lpx(li,lj) = lpx(li,lj) + 1
                if (lpx(li,lj) > mx) &
                     call upf_error('aainit','mx dimension too small', lpx(li,lj))
                lpl(li,lj,lpx(li,lj)) = l
             end if
          end do
       end do
    end do
    ! 
    deallocate(mly)
    deallocate(ylm)
    deallocate(rr)
    deallocate(r)
    !
#if defined (__CUDA)
    IF (ALLOCATED(ap_d)) DEALLOCATE(ap_d)
    ALLOCATE(ap_d, SOURCE=ap)
    IF (ALLOCATED(lpx_d)) DEALLOCATE(lpx_d)
    ALLOCATE(lpx_d, SOURCE=lpx)
    IF (ALLOCATED(lpl_d)) DEALLOCATE(lpl_d)
    ALLOCATE(lpl_d, SOURCE=lpl)
#endif
    return
  end subroutine aainit
  !
  !-----------------------------------------------------------------------
  subroutine gen_rndm_r(llx,r,rr)
    !-----------------------------------------------------------------------
    ! - generate an array of random vectors (uniform deviate on unitary sphere)
    !
    USE upf_const,  ONLY: tpi
    implicit none
    !
    ! first the I/O variables
    !
    integer :: llx   ! input: the dimension of r and rr
    real(DP) :: &
         r(3,llx),  &! output: an array of random vectors
         rr(llx)     ! output: the norm of r
    !
    ! here the local variables
    !
    integer  :: ir
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
  end subroutine gen_rndm_r
!
!------------------------------------------------------------------------
   FUNCTION randy ( irand )
      !------------------------------------------------------------------------
      !   
      ! x=randy(n): reseed with initial seed idum=n ( 0 <= n <= ic, see below)
      !             if randy is not explicitly initialized, it will be
      !             initialized with seed idum=0 the first time it is called
      ! x=randy() : generate uniform real(DP) numbers x in [0,1]
      !   
      use upf_kinds, only : DP
      implicit none
      !
      REAL(DP) :: randy
      INTEGER, optional    :: irand
      !   
      INTEGER , PARAMETER  :: m    = 714025, &
                              ia   = 1366, &
                              ic   = 150889, &
                              ntab = 97
      REAL(DP), PARAMETER  :: rm = 1.0_DP / m 
      INTEGER              :: j
      INTEGER, SAVE        :: ir(ntab), iy, idum=0
      LOGICAL, SAVE        :: first=.true.
      !   
      IF ( present(irand) ) THEN
         idum = MIN( ABS(irand), ic) 
         first=.true.
      END IF

      IF ( first ) THEN
         !   
         first = .false.
         idum = MOD( ic - idum, m ) 
         !   
         DO j=1,ntab
            idum=mod(ia*idum+ic,m)
            ir(j)=idum
         END DO
         idum=mod(ia*idum+ic,m)
         iy=idum
      END IF
      j=1+(ntab*iy)/m
      IF( j > ntab .OR. j <  1 ) call upf_error('randy','j out of range',ABS(j)+1)
      iy=ir(j)
      randy=iy*rm
      idum=mod(ia*idum+ic,m)
      ir(j)=idum
      !   
      RETURN
      !   
   END FUNCTION randy
  !
  !-----------------------------------------------------------------------
  function compute_ap(l,li,lj,llx,ylm,mly)
    !-----------------------------------------------------------------------
    !-  given an l and a li,lj pair compute ap(l,li,lj)
    implicit none
    !
    ! first the I/O variables
    !
    integer :: &
         llx,         &! the dimension of ylm and mly
         l,li,lj       ! the arguments of the array ap
    
    real(DP) :: &
         compute_ap,  &! this function
         ylm(llx,llx),&! the real spherical harmonics for array r
         mly(llx,llx)  ! the inverse of ylm considered as a matrix
    !
    ! here the local variables
    !
    integer :: ir
    
    compute_ap = 0.0_DP
    do ir = 1,llx
       compute_ap = compute_ap + mly(l,ir)*ylm(ir,li)*ylm(ir,lj)
    end do
    
    return
  end function compute_ap
  !
  !-----------------------------------------------------------------------
  subroutine allocate_uspp(use_gpu,noncolin,lspinorb,tqr,nhm,nsp,nat,nspin)
    !-----------------------------------------------------------------------
    implicit none
    logical, intent(in) :: use_gpu
    logical, intent(in) :: noncolin,lspinorb,tqr
    integer, intent(in) :: nhm,nsp,nat,nspin
    !
    !if (nhm_/=nhm) call errore("allocate_uspp","invalid nhm",1)
    !if (nsp_/=nsp) call errore("allocate_uspp","invalid nsp",1)
    !
    allocate( indv(nhm,nsp)   )
    allocate( nhtol(nhm,nsp)  )
    allocate( nhtolm(nhm,nsp) )
    allocate( nhtoj(nhm,nsp)  )
    allocate( ijtoh(nhm,nhm,nsp) )
    allocate( deeq(nhm,nhm,nat,nspin) )
    if ( noncolin ) then
       allocate( deeq_nc(nhm,nhm,nat,nspin) )
    endif
    allocate( qq_at(nhm,nhm,nat) )
    allocate( qq_nt(nhm,nhm,nsp) )
    if ( lspinorb ) then
       allocate( qq_so(nhm,nhm,4,nsp) )
       allocate( dvan_so(nhm,nhm,nspin,nsp) )
       allocate( fcoef(nhm,nhm,2,2,nsp) )
    else
       allocate( dvan(nhm,nhm,nsp) )
    endif
    allocate(becsum( nhm*(nhm+1)/2, nat, nspin))
    if (tqr) then
       allocate(ebecsum( nhm*(nhm+1)/2, nat, nspin))
    endif
    allocate( indv_ijkb0(nat) )
    !
    ! GPU-vars (protecting zero-size allocations)
    !
    if (use_gpu) then
      !
      if (nhm>0) then
        allocate( indv_d(nhm,nsp)   )
        allocate( nhtol_d(nhm,nsp)  )
        allocate( nhtolm_d(nhm,nsp) )
        allocate( nhtoj_d(nhm,nsp)  )
        allocate( ijtoh_d(nhm,nhm,nsp) )
        allocate( deeq_d(nhm,nhm,nat,nspin) )
        if ( noncolin ) then
           allocate( deeq_nc_d(nhm,nhm,nat,nspin) )
        endif
        allocate( qq_at_d(nhm,nhm,nat) )
        allocate( qq_nt_d(nhm,nhm,nsp) )
        if ( lspinorb ) then
           allocate( qq_so_d(nhm,nhm,4,nsp) )
           allocate( dvan_so_d(nhm,nhm,nspin,nsp) )
           allocate( fcoef_d(nhm,nhm,2,2,nsp) )
        else
           allocate( dvan_d(nhm,nhm,nsp) )
        endif
        allocate(becsum_d( nhm*(nhm+1)/2, nat, nspin))
        if (tqr) then
           allocate(ebecsum_d( nhm*(nhm+1)/2, nat, nspin))
        endif
        !
      endif
      allocate( indv_ijkb0_d(nat) )
      !
    endif
    !
  end subroutine allocate_uspp
  !
  !-----------------------------------------------------------------------
  SUBROUTINE deallocate_uspp()
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    IF( ALLOCATED( nhtol ) )      DEALLOCATE( nhtol )
    IF( ALLOCATED( indv ) )       DEALLOCATE( indv )
    IF( ALLOCATED( nhtolm ) )     DEALLOCATE( nhtolm )
    IF( ALLOCATED( nhtoj ) )      DEALLOCATE( nhtoj )
    IF( ALLOCATED( indv_ijkb0 ) ) DEALLOCATE( indv_ijkb0 )
    IF( ALLOCATED( ijtoh ) )      DEALLOCATE( ijtoh )
    IF( ALLOCATED( vkb ) )        DEALLOCATE( vkb )
    IF( ALLOCATED( becsum ) )     DEALLOCATE( becsum )
    IF( ALLOCATED( ebecsum ) )    DEALLOCATE( ebecsum )
    IF( ALLOCATED( qq_at ) )      DEALLOCATE( qq_at )
    IF( ALLOCATED( qq_nt ) )      DEALLOCATE( qq_nt )
    IF( ALLOCATED( dvan ) )       DEALLOCATE( dvan )
    IF( ALLOCATED( deeq ) )       DEALLOCATE( deeq )
    IF( ALLOCATED( qq_so ) )      DEALLOCATE( qq_so )
    IF( ALLOCATED( dvan_so ) )    DEALLOCATE( dvan_so )
    IF( ALLOCATED( deeq_nc ) )    DEALLOCATE( deeq_nc )
    IF( ALLOCATED( fcoef ) )      DEALLOCATE( fcoef )
    IF( ALLOCATED( beta ) )       DEALLOCATE( beta )
    IF( ALLOCATED( dbeta ) )      DEALLOCATE( dbeta )
    !
    ! GPU variables
    IF( ALLOCATED( ap_d ) )       DEALLOCATE( ap_d )
    IF( ALLOCATED( lpx_d ) )      DEALLOCATE( lpx_d )
    IF( ALLOCATED( lpl_d ) )      DEALLOCATE( lpl_d )
    IF( ALLOCATED( indv_d ) )     DEALLOCATE( indv_d )
    IF( ALLOCATED( nhtol_d ) )    DEALLOCATE( nhtol_d )
    IF( ALLOCATED( nhtolm_d ) )   DEALLOCATE( nhtolm_d )
    IF( ALLOCATED( ijtoh_d ) )    DEALLOCATE( ijtoh_d )
    IF( ALLOCATED( indv_ijkb0_d)) DEALLOCATE( indv_ijkb0_d )
    IF( ALLOCATED( vkb_d ) )      DEALLOCATE( vkb_d )
    IF( ALLOCATED( becsum_d ) )   DEALLOCATE( becsum_d )
    IF( ALLOCATED( ebecsum_d ) )  DEALLOCATE( ebecsum_d )
    IF( ALLOCATED( dvan_d ) )     DEALLOCATE( dvan_d )
    IF( ALLOCATED( deeq_d ) )     DEALLOCATE( deeq_d )
    IF( ALLOCATED( qq_nt_d ) )    DEALLOCATE( qq_nt_d )
    IF( ALLOCATED( qq_at_d ) )    DEALLOCATE( qq_at_d )
    IF( ALLOCATED( nhtoj_d ) )    DEALLOCATE( nhtoj_d )
    IF( ALLOCATED( qq_so_d ) )    DEALLOCATE( qq_so_d )
    IF( ALLOCATED( dvan_so_d ) )  DEALLOCATE( dvan_so_d )
    IF( ALLOCATED( deeq_nc_d ) )  DEALLOCATE( deeq_nc_d )
    IF( ALLOCATED( fcoef_d ) )    DEALLOCATE( fcoef_d )
    !
  END SUBROUTINE deallocate_uspp
  !
  ! In the following, allocations are not checked and assume
  ! to be dimensioned correctly.
  !
  ! intento is used to specify what the variable will  be used for :
  !  0 -> in , the variable needs to be synchronized but won't be changed
  !  1 -> inout , the variable needs to be synchronized AND will be changed
  !  2 -> out , NO NEED to synchronize the variable, everything will be overwritten
  ! 
  SUBROUTINE using_vkb(intento)
      !
      !USE uspp, ONLY : vkb, vkb_d
      implicit none
      INTEGER, INTENT(IN) :: intento
      IF (size(vkb)==0) return
#if defined(__CUDA)  || defined(__CUDA_GNU)
      IF (vkb_ood) THEN
          IF (intento < 2) vkb = vkb_d
          vkb_ood = .false.
      ENDIF
      IF (intento > 0)    vkb_d_ood = .true.
#endif
  END SUBROUTINE using_vkb
  !
  SUBROUTINE using_vkb_d(intento)
      !
      !USE uspp, ONLY : vkb, vkb_d
      implicit none
      INTEGER, INTENT(IN) :: intento
      IF (size(vkb)==0) return
#if defined(__CUDA) || defined(__CUDA_GNU)
      !
      IF (vkb_d_ood) THEN
          IF (intento < 2) vkb_d = vkb
          vkb_d_ood = .false.
      ENDIF
      IF (intento > 0)    vkb_ood = .true.
#else
      CALL errore('using_vkb_d', 'no GPU support', 1)
#endif
  END SUBROUTINE using_vkb_d
  !   
END MODULE uspp

