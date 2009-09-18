! OBM :
! 050608 Modified for calbec interface in v4.0
!        gamma_only correction
!        tvanp --> upf%tvanp
! 150608 A lot of changes to variables. WARNING see Changelog-4.0 5646, kkbeta may cause problems
!           -a large quantitiy of parameters moved to upf still in uspp_param not sure about kkbeta
!           -atom r,rab moved to rgrid still inside atom
!           -nbrx is depreciated now, only gipaw_module uses it, no need to compile all the magnetic
!                 resonance code for just this parameter, moving it to lr_variables
! 160608 compute_qdipol degeneracy with pw/compute_qdipol prefixing the native one with lr
!        reduce-->mp_sum
! 160709 K point correction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE sm1_psi( recalc, ik, lda, n, m, psi, spsi)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !----------------------------------------------------------------------------
  !
  !    This routine applies the S^{-1} matrix to m wavefunctions psi
  !    and puts the results in spsi.
  !    Requires the products of psi with all beta functions
  !    in array becp(nkb,m) (calculated in h_psi or by ccalbec)
  ! input:
  !     recalc decides if the overlap of beta functions is recalculated or not.
  !            this is needed e.g. if ions are moved and the overlap changes accordingly
  !     lda   leading dimension of arrays psi, spsi
  !     n     true dimension of psi, spsi
  !     m     number of states psi
  !     psi
  ! output:
  !     spsi  S^{-1}*psi
  !

  USE kinds,      ONLY : DP
  USE control_flags,      ONLY : gamma_only
  USE uspp,       ONLY : okvan, vkb, nkb, qq
  USE uspp_param, ONLY : nh, upf
  USE wvfct,      ONLY : igk, g2kin
  USE gsmooth,    ONLY : nls, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nrxxs
  USE ldaU,       ONLY : lda_plus_u
  USE ions_base,  ONLY : ityp,nat,ntyp=>nsp
  use mp,         only : mp_sum
  USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm
  USE lr_variables,   ONLY : lr_verbosity
  USE io_global,      ONLY : stdout
  !
  IMPLICIT NONE
  !
  ! ... First the dummy variables
  !
  LOGICAL          :: recalc
  INTEGER          :: lda, n, m, ik
  COMPLEX(KIND=DP) :: psi(lda,m), spsi(lda,m)
  !
  If (lr_verbosity > 5) THEN
    WRITE(stdout,'("<lr_sm1_psi>")')
  endif
  !
  CALL start_clock( 'lr_sm1_psi' )  
  !
  IF ( gamma_only ) THEN
     CALL sm1_psi_gamma()
  ELSE
     !
     CALL sm1_psi_k()
     !
  END IF
  !
  CALL stop_clock( 'lr_sm1_psi' )
  !
  RETURN
  !
CONTAINS
  !-----------------------------------------------------------------------
  SUBROUTINE sm1_psi_gamma()
    !-----------------------------------------------------------------------
    !
    USE becmod,               ONLY : bec_type,becp,calbec
    !use real_beta,            only : ccalbecr_gamma,s_psir,fft_orbital_gamma
    !use lr_variables,         only : real_space
    USE realus,              ONLY : real_space, fft_orbital_gamma, initialisation_level, &
                                    bfft_orbital_gamma, calbec_rs_gamma, add_vuspsir_gamma, &
                                    v_loc_psir, s_psir_gamma,check_fft_orbital_gamma, real_space_debug

    !
    IMPLICIT NONE
    !
    ! ... local variables
    !
    INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0, ibnd, ii
    ! counters
    real(KIND=DP), ALLOCATABLE :: ps(:,:)
    real(kind=dp), allocatable, save :: BB_(:,:)
    logical, save :: first_entry = .true.
    if(first_entry) then
      if(allocated(BB_)) deallocate(BB_)
      first_entry = .false.
    endif


    ! the product vkb and psi
    !
    ! ... initialize  spsi
    !
    CALL ZCOPY( lda * m, psi, 1, spsi, 1 )
    !
    ! ... The product with the beta functions
    !
    IF ( nkb == 0 .OR. .NOT. okvan ) RETURN
    !
    !BB_ = sum
    !if (allocated(BB_)) then 
    !  print *, "BB is allocated, ", BB_(1,1)
    !else
    !  print *, "BB is not allocated"
    !endif 
    if(recalc .and. .not.allocated(BB_)) then
       allocate(BB_(nkb,nkb))
       BB_=0.d0
       call errore('sm1_psi','recalculating BB_ matrix',-1)
       !print *, "did you see the recalculating message?"
       if (lr_verbosity > 1) then
          WRITE(stdout,'(5X,"Calculating S^-1")')
       endif
       !call pw_gemm('Y',nkb,nkb,n,vkb,lda,vkb,lda,BB_,nkb) 
       call calbec (n,vkb,vkb,BB_,nkb)
       ALLOCATE( ps( nkb, nkb ) )    
       ps(:,:) = (0.d0)
       ijkb0 = 0
       do nt=1,ntyp
          if (upf(nt)%tvanp) then
             do na=1,nat
                if(ityp(na).eq.nt) then
                   do ii=1,nkb
                      do jh=1,nh(nt)
                         jkb=ijkb0 + jh
                         do ih=1,nh(nt)
                            ikb = ijkb0 + ih
                            ps(ikb,ii) = ps(ikb,ii) + qq(ih,jh,nt)*BB_(jkb,ii)
                         enddo
                      enddo
                   enddo
                   ijkb0 = ijkb0+nh(nt)
                endif
             enddo
          else
             DO na = 1, nat
                IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
             END DO
          endif
       enddo
       
       do ii=1,nkb
          ps(ii,ii) = ps(ii,ii) + 1.d0
       enddo
       
       call dinv_matrix(ps,nkb)
       BB_(:,:) = 0.d0
       ijkb0 = 0
       do nt=1,ntyp
          if (upf(nt)%tvanp) then
             do na=1,nat
                if(ityp(na).eq.nt) then
                   do ii=1,nkb
                      do jh=1,nh(nt)
                         jkb=ijkb0 + jh
                         do ih=1,nh(nt)
                            ikb = ijkb0 + ih
                            BB_(ii,jkb) = BB_(ii,jkb) - ps(ii,ikb)*qq(ih,jh,nt)
                         enddo
                      enddo
                   enddo
                   ijkb0 = ijkb0+nh(nt)
                endif
             enddo
          else
             DO na = 1, nat
                IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
             END DO
          endif
       enddo
       
       deallocate(ps)
    endif
    !print *, "BB is now, ", BB_(1,1)
   
    if (real_space_debug>3) then !was 3
      do ibnd=1,m,2
       call fft_orbital_gamma(psi,ibnd,m)
       call calbec_rs_gamma(ibnd,m,becp%r)
      enddo
    else
     call calbec(n,vkb,psi,becp,m)
    !call pw_gemm('Y',nkb,m,n,vkb,lda,psi,lda,rbecp,nkb) 
    endif
    !
    ALLOCATE( ps( nkb, m ) )    
!    ps(:,:) = 0.D0
    !
!    do ibnd=1,m
!       do jkb=1,nkb
!          do ii=1,nkb
!             ps(jkb,ibnd) = ps(jkb,ibnd)+BB_(jkb,ii)*rbecp(ii,ibnd)
!          enddo
!       enddo
!    enddo
    !
    call DGEMM( 'N','N',nkb,m,nkb,1.d0,BB_,nkb,becp%r,nkb,0.d0,ps,nkb)


!   do ibnd=1,m
!      do ii=1,nkb
!          call ZAXPY(n,cmplx(ps(ii,ibnd),0.0d0,dp),vkb(1,ii),1,spsi(1,ibnd),1) 
!       enddo
!    enddo 
    call DGEMM('N','N',2*n,m,nkb,1.d0,vkb,2*lda,ps,nkb,1.d0,spsi,2*lda)

    !
    DEALLOCATE( ps )
    !
    RETURN
    !
  END SUBROUTINE sm1_psi_gamma
  !
  !-----------------------------------------------------------------------
  SUBROUTINE sm1_psi_k()
    !-----------------------------------------------------------------------
    !
    ! ... k-points version
    !
    USE becmod,        ONLY : bec_type,becp,calbec
    !USE lr_variables,    ONLY: igk_k, npw_k
    use realus,        only : igk_k,npw_k
    USE klist,         only : nks, xk
    !
    IMPLICIT NONE
    !
    ! ... local variables
    !
    INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0, ibnd, ii, ik1
    ! counters
    complex(KIND=DP), ALLOCATABLE :: ps(:,:)
    complex(kind=dp), allocatable, save :: BB_(:,:,:)

    ! the product vkb and psi
    !
    ! ... initialize  spsi
    !
    CALL ZCOPY( lda * m, psi, 1, spsi, 1 )
    !
    ! ... The product with the beta functions
    !
    IF ( nkb == 0 .OR. .NOT. okvan ) RETURN
    !
    if(recalc .and. .not. allocated(BB_)) then
       allocate(BB_(nkb,nkb,nks))
       BB_=(0.d0,0.d0)
       call errore('sm1_psi','recalculating BB_ matrix',-1)
       
       ALLOCATE( ps( nkb, nkb ) )    
       
       do ik1 = 1,nks
          call init_us_2(npw_k(ik1),igk_k(:,ik1),xk(1,ik1),vkb)
          call zgemm('C','N',nkb,nkb,npw_k(ik1),(1.d0,0.d0),vkb,lda,vkb,lda,(0.d0,0.d0),BB_(1,1,ik1),nkb)
#ifdef __PARA
          !CALL reduce( 2 * nkb * nkb, BB_(:,:,ik1) )
          call mp_sum(BB_(:,:,ik1), intra_pool_comm)
#endif
       
          ps(:,:) = (0.d0,0.d0)
          ijkb0 = 0
          do nt=1,ntyp
             if (upf(nt)%tvanp) then
                do na=1,nat
                   if(ityp(na).eq.nt) then
                      do ii=1,nkb
                         do jh=1,nh(nt)
                            jkb=ijkb0 + jh
                            do ih=1,nh(nt)
                               ikb = ijkb0 + ih
                               ps(ikb,ii) = ps(ikb,ii) + BB_(jkb,ii, ik1)*qq(ih,jh,nt)
                            enddo
                         enddo
                      enddo
                      ijkb0 = ijkb0+nh(nt)
                   endif
                enddo
             else
                DO na = 1, nat
                   IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
                END DO
             endif
          enddo
       
          do ii=1,nkb
             ps(ii,ii) = ps(ii,ii) + (1.d0,0.d0)
          enddo
       
          call zinv_matrix(ps,nkb)
          BB_(:,:,ik1) = (0.d0,0.d0)
          ijkb0 = 0
          do nt=1,ntyp
             if (upf(nt)%tvanp) then
                do na=1,nat
                   if(ityp(na).eq.nt) then
                      do ii=1,nkb
                         do jh=1,nh(nt)
                            jkb=ijkb0 + jh
                            do ih=1,nh(nt)
                               ikb = ijkb0 + ih
                               BB_(ii,jkb,ik1) = BB_(ii,jkb,ik1) - ps(ii,ikb)*qq(ih,jh,nt)
                            enddo
                         enddo
                      enddo
                      ijkb0 = ijkb0+nh(nt)
                   endif
                enddo
             else
                DO na = 1, nat
                   IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
                END DO
             endif
          enddo
       enddo
       deallocate(ps)
    endif
    
    call init_us_2(npw_k(ik),igk_k(:,ik),xk(1,ik),vkb)
    !call ccalbec( nkb, lda, n, m, becp, vkb, psi )
    call calbec(n,vkb,psi,becp,m)

    !
    ALLOCATE( ps( nkb, m ) )    
    ps(:,:) = (0.d0,0.d0)
    !
    do ibnd=1,m
       do jkb=1,nkb
          do ii=1,nkb
             ps(jkb,ibnd) = ps(jkb,ibnd)+BB_(jkb,ii,ik)*becp%k(ii,ibnd)
          enddo
       enddo
    enddo
    !
    !
    CALL ZGEMM( 'N', 'N', n, m, nkb, (1.D0, 0.D0), vkb, &
         lda, ps, nkb, (1.D0, 0.D0), spsi, lda )


    DEALLOCATE( ps )
    !
    RETURN
    !
  END SUBROUTINE sm1_psi_k
  !
END subroutine sm1_psi
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dinv_matrix(M,N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  USE kinds,      ONLY : DP

  implicit none

  integer :: N                              !  matrix dimension
  real(kind=dp), dimension(0:N-1,0:N-1) :: M ! MAtrix to be inverted

  real(kind=dp), dimension(:), allocatable :: work
  integer, dimension(:), allocatable :: ipiv

  integer :: i,lwork,info
  integer, save :: lworkfact

  data lworkfact /64/

  lwork = lworkfact*N

  allocate(ipiv(0:N-1))
  allocate(work(1:lwork))

! Factorize Matrix M

  call dgetrf( N, N, M, N, ipiv, info )
  if (info.ne.0) then
     call errore('dinv_matrix','error in dgetrf',info)
  endif

! Invert Matrix

  call dgetri( N, M, N, ipiv, work, lwork, info )
  if (info.ne.0) then
     call errore('dinv_matrix','error in dgetri',info)
  else
     lworkfact = int(work(1)/N)
  endif
     
  deallocate(work)
  deallocate(ipiv)

end subroutine dinv_matrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zinv_matrix(M,N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  USE kinds,      ONLY : DP

  implicit none

  integer :: N                            !  matrix dimension
  complex(kind=dp), dimension(0:N-1,0:N-1) :: M ! MAtrix to be inverted

  complex(kind=dp), dimension(:), allocatable :: work
  integer, dimension(:), allocatable :: ipiv

  integer :: i,lwork,info
  integer, save :: lworkfact

  data lworkfact /64/

  lwork = lworkfact*N

  allocate(ipiv(0:N-1))
  allocate(work(1:lwork))

! Factorize Matrix M

  call zgetrf( N, N, M, N, ipiv, info )
  if (info.ne.0) then
     call errore('zinv_matrix','error in zgetrf',info)
  endif

! Invert Matrix

  call zgetri( N, M, N, ipiv, work, lwork, info )
  if (info.ne.0) then
     call errore('zinv_matrix','error in zgetri',info)
  else
     lworkfact = int(work(1)/N)
  endif
     
  deallocate(work)
  deallocate(ipiv)

end subroutine zinv_matrix
!---------------------------------------------------------------------- 
subroutine lr_adddvepsi_us_gamma(becp1,becp2,ipol,kpoint,dvpsi) 
  !
  ! This subdoutine adds to dvpsi the terms which depend on the augmentation 
  ! charge. It assume that the variable dpqq, has been set. 
  ! 
#include "f_defs.h" 
 
use ions_base,             only : ityp,nat,ntyp=>nsp
use cell_base,             only : at
use uspp
use uspp_param
use wvfct, only : npw, npwx, nbnd
use control_flags, only : gamma_only
USE kinds, only : DP
 
implicit none 
 
integer, intent(in) :: ipol, kpoint 
real(kind=dp), intent(in) :: becp1(nkb,nbnd),becp2(nkb,nbnd) 
complex(kind=dp) :: dvpsi(npwx,nbnd) 
 
real(kind=dp), allocatable :: dpqq(:,:,:,:) 
real(kind=dp) :: fact 
real(kind=dp), allocatable :: ps(:) 
integer:: ijkb0, nt, na, ih, jh, ikb, jkb, ibnd 
 
allocate (dpqq( nhm, nhm, 3, ntyp)) 
 
allocate (ps(nbnd))     
call lr_compute_qdipol(dpqq) 
 
ijkb0 = 0 
do nt = 1, ntyp 
   do na = 1, nat 
      if (ityp(na).eq.nt) then 
         do ih = 1, nh (nt) 
            ikb = ijkb0 + ih 
            ps = 0.0d0
            do jh = 1, nh (nt) 
               jkb = ijkb0 + jh 
               fact=at(1,ipol)*dpqq(ih,jh,1,nt)+  & 
                    at(2,ipol)*dpqq(ih,jh,2,nt)+  & 
                    at(3,ipol)*dpqq(ih,jh,3,nt) 
               do ibnd=1, nbnd 
                  ps(ibnd) = ps(ibnd)                             & 
                     + becp2(jkb,ibnd)*qq(ih,jh,nt)+  & 
                       becp1(jkb,ibnd)*fact 
               enddo 
            enddo 
            do ibnd = 1, nbnd
               call ZAXPY(npw,cmplx(ps(ibnd),0.d0,DP),vkb(1,ikb),1,dvpsi(1,ibnd),1) 
            enddo 
         enddo 
         ijkb0=ijkb0+nh(nt) 
      endif 
   enddo 
enddo 
if (jkb.ne.nkb) call errore ('lr_adddvepsi_us', 'unexpected error', 1) 
 
deallocate(ps) 
deallocate(dpqq) 

return 
end subroutine lr_adddvepsi_us_gamma
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine lr_adddvepsi_us_k(becp1,becp2,ipol,kpoint,dvpsi) 
  !
  ! This subdoutine adds to dvpsi the terms which depend on the augmentation 
  ! charge. It assume that the variable dpqq, has been set. 
  ! 
#include "f_defs.h" 
 
use ions_base,             only : ityp,nat,ntyp=>nsp
use cell_base,             only : at
use uspp
use uspp_param
use wvfct, only : npw, npwx, nbnd
use control_flags, only : gamma_only
USE kinds, only : DP
 
implicit none 
 
integer, intent(in) :: ipol, kpoint 
complex(kind=dp), intent(in) :: becp1(nkb,nbnd),becp2(nkb,nbnd) 
complex(kind=dp) :: dvpsi(npwx,nbnd) 
 
real(kind=dp), allocatable :: dpqq(:,:,:,:) 
real(kind=dp) :: fact 
complex(kind=dp), allocatable :: ps(:) 
integer:: ijkb0, nt, na, ih, jh, ikb, jkb, ibnd 
 
allocate (dpqq( nhm, nhm, 3, ntyp)) 
 
allocate (ps(nbnd))     
call lr_compute_qdipol(dpqq) 
 
ijkb0 = 0 
do nt = 1, ntyp 
   do na = 1, nat 
      if (ityp(na).eq.nt) then 
         do ih = 1, nh (nt) 
            ikb = ijkb0 + ih 
            ps = (0.d0,0.d0) 
            do jh = 1, nh (nt) 
               jkb = ijkb0 + jh 
               fact=at(1,ipol)*dpqq(ih,jh,1,nt)+  & 
                    at(2,ipol)*dpqq(ih,jh,2,nt)+  & 
                    at(3,ipol)*dpqq(ih,jh,3,nt) 
               do ibnd=1, nbnd 
                  ps(ibnd) = ps(ibnd)                             & 
                     + becp2(jkb,ibnd)*qq(ih,jh,nt)+  & 
                       becp1(jkb,ibnd)*fact 
               enddo 
            enddo 
            do ibnd = 1, nbnd 
               call ZAXPY(npw,ps(ibnd),vkb(1,ikb),1,dvpsi(1,ibnd),1) 
            enddo 
         enddo 
         ijkb0=ijkb0+nh(nt) 
      endif 
   enddo 
enddo 
if (jkb.ne.nkb) call errore ('lr_adddvepsi_us', 'unexpected error', 1) 
 
deallocate(ps) 
deallocate(dpqq) 
 
return 
end subroutine lr_adddvepsi_us_k
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine lr_compute_qdipol(dpqq)
  !
  ! This routine computes the term dpqq, i.e. the dipole moment of the
  ! augmentation charge
  !
  USE kinds, only: DP
  USE constants, ONLY: fpi
  USE atom, ONLY: rgrid!r, rab
  USE ions_base, ONLY: ntyp => nsp
  USE uspp, only: nhtol, nhtolm, indv, nlx, ap
  USE uspp_param, only: upf, nh, nhm !nbeta, lll, kkbeta, qfunc, rinner, &
  !     qfcoef, nqf, upf, nh, nhm, nbrx
  USE lr_variables, only: nbrx ! look for a way to remove dependency to nbrx, how is it done in pwscf? 
  implicit none

  real(kind=dp) :: dpqq(nhm, nhm, 3, ntyp)
  real(DP), allocatable :: qrad2(:,:,:), qtot(:,:,:), aux(:)
  real(DP) :: fact
  integer :: nt, l, ir, nb, mb, ijv, ilast, ipol, ih, ivl, jh, jvl, lp, ndm

  call start_clock('cmpt_qdipol')
  ndm = MAXVAL (upf(1:ntyp)%kkbeta) !MAXVAL (kkbeta(1:ntyp))
  allocate (qrad2( nbrx , nbrx, ntyp))    
  allocate (aux( ndm))    
  allocate (qtot( ndm, nbrx, nbrx))    

  qrad2(:,:,:)=0.d0
  dpqq=0.d0

  do nt = 1, ntyp
     if (upf(nt)%tvanp ) then
        l=1
!
!   Only l=1 terms enter in the dipole of Q
!
        do nb = 1, upf(nt)%nbeta
           do mb = nb, upf(nt)%nbeta
              ijv = mb * (mb-1) /2 + nb
              if ((l.ge.abs(upf(nt)%lll(nb)-upf(nt)%lll(mb))) .and. &
                   (l.le.upf(nt)%lll(nb)+upf(nt)%lll(mb))      .and. &
                   (mod (l+upf(nt)%lll(nb)+upf(nt)%lll(mb),2) .eq.0) ) then
                 do ir = 1, upf(nt)%kkbeta
                    if (rgrid(nt)%r(ir).ge.upf(nt)%rinner(l+1)) then
                       qtot(ir, nb, mb)=upf(nt)%qfunc(ir,ijv)
                    else
                       ilast = ir
                    endif
                 enddo
                 if (upf(nt)%rinner(l+1).gt.0.d0) &
                      call setqf(upf(nt)%qfcoef (1, l+1, nb, mb), &
                      qtot(1,nb,mb), rgrid(nt)%r(1), upf(nt)%nqf,l,ilast)
              endif
           enddo
        enddo
        do nb=1, upf(nt)%nbeta
           !
           !    the Q are symmetric with respect to indices
           !
           do mb=nb, upf(nt)%nbeta
              if ( (l.ge.abs(upf(nt)%lll(nb)-upf(nt)%lll(mb) ) )    .and.  &
                   (l.le.upf(nt)%lll(nb) + upf(nt)%lll(mb) )        .and.  &
                   (mod(l+upf(nt)%lll(nb)+upf(nt)%lll(mb), 2).eq.0) ) then
                 do ir = 1, upf(nt)%kkbeta
                    aux(ir)=rgrid(nt)%r(ir)*qtot(ir, nb, mb)
                 enddo
                 call simpson (upf(nt)%kkbeta,aux,rgrid(nt)%rab(1),qrad2(nb,mb,nt))
              endif
           enddo
        enddo
     endif
     ! ntyp
  enddo
  
  do ipol = 1,3
     fact=-sqrt(fpi/3.d0)
     if (ipol.eq.1) lp=3
     if (ipol.eq.2) lp=4
     if (ipol.eq.3) then
        lp=2
        fact=-fact
     endif
     do nt = 1,ntyp
        if (upf(nt)%tvanp) then
           do ih = 1, nh(nt)
              ivl = nhtolm(ih, nt)
              mb = indv(ih, nt)
              do jh = ih, nh (nt)
                 jvl = nhtolm(jh, nt)
                 nb=indv(jh,nt)
                 if (ivl > nlx) call errore('lr_compute_qdipol',' ivl > nlx', ivl)
                 if (jvl > nlx) call errore('lr_compute_qdipol',' jvl > nlx', jvl)
                 if (nb > nbrx) call errore('lr_compute_qdipol',' nb > nbrx', nb)
                 if (mb > nbrx) call errore('lr_compute_qdipol',' mb > nbrx', mb)
                 if (mb > nb) call errore('lr_compute_qdipol',' mb > nb', 1)
                 dpqq(ih,jh,ipol,nt)=fact*ap(lp,ivl,jvl)*qrad2(mb,nb,nt)
                 dpqq(jh,ih,ipol,nt)=dpqq(ih,jh,ipol,nt)
                 ! WRITE( stdout,'(3i5,2f15.9)') ih,jh,ipol,dpqq(ih,jh,ipol,nt)
              enddo
           enddo
        endif
     enddo
  enddo
  deallocate(qtot)
  deallocate(aux)
  deallocate(qrad2)
  call stop_clock('cmpt_qdipol')

  return
end subroutine lr_compute_qdipol
