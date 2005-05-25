!
! Copyright (C) 2002-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#include "f_defs.h"
!
#define ZERO ( 0.D0, 0.D0 )
!
#undef DEBUG
!-----------------------------------------------------------------------
SUBROUTINE mix_rho_nc( rhout, rhoin, nsout, nsin, alphamix, &
                       dr2, tr2_min, iter, n_iter, filename, conv )
  !-----------------------------------------------------------------------
  !
  ! Modified Broyden's method for charge density mixing
  !             d.d. johnson prb 38, 12807 (1988)
  ! On output: the mixed density is in rhoin, rhout is UNCHANGED
  !
  USE ions_base,  ONLY : nat, ityp
  USE kinds, only : DP
  USE wavefunctions_module, ONLY : psic => psic_nc
  use pwcom
  USE control_flags,        ONLY : imix, ngm0, tr2
  USE parser,               ONLY : find_free_unit
  implicit none
  !
  !   First the I/O variable
  !
  character (len=256) ::  &
                filename     !  (in) I/O filename for mixing history
                             !  if absent everything is kept in memory
  integer ::    &
                iter,       &!  (in)  counter of the number of iterations
                n_iter       !  (in)  numb. of iterations used in mixing

  real (kind=DP) :: &
                rhout(nrxx,nspin), &! (in) the "out" density
                rhoin(nrxx,nspin), &! (in) the "in" density; (out) the new dens.
                nsout(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat), &!
                nsin(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat),  &!
                alphamix,          &! (in) mixing factor
                dr2                 ! (out) the estimated errr on the energy
  REAL (KIND=DP) :: &
    tr2_min       ! estimated error from diagonalization. If the estimated scf 
                  ! error is smaller than this, exit: a more accurate 
                  ! diagonalization is needed
  logical ::    &
                conv        ! (out) if true the convergence has been reached
  !
  integer, parameter :: &
                maxmix =25  ! max number of iterations for charge mixing

  !
  !   Here the local variables
  !
  integer ::    &
                iunmix,    &! I/O unit number of charge density file
                iunmix2,   &! I/O unit number of ns file
                iunit,     &! counter on I/O unit numbers
                iter_used, &! actual number of iterations used
                ipos,      &! index of the present iteration
                inext,     &! index of the next iteration
                i, j,      &! counters on number of iterations
                is,        &! counter on spin component
                ig,        &! counter on G-vectors
                iwork(maxmix),&! dummy array used as output by libr. routines
                info        ! flag saying if the exec. of libr. routines was ok

  complex (kind=DP), allocatable :: rhocin(:,:), rhocout(:,:), &
                rhoinsave(:,:), rhoutsave(:,:), &
                nsinsave(:,:,:,:),  nsoutsave(:,:,:,:)
  complex (kind=DP), allocatable, save :: df(:,:,:), dv(:,:,:), &
                                      df_ns(:,:,:,:,:), dv_ns(:,:,:,:,:)
                ! rhocin(ngm0,nspin)
                ! rhocout(ngm0,nspin)
                ! rhoinsave(ngm0*nspin): work space
                ! rhoutsave(ngm0*nspin): work space
                ! df(ngm0,nspin,n_iter): information from preceding iterations
                ! dv(ngm0,nspin,n_iter):    "  "       "     "        "  "
                ! df_ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat,n_iter):idem
                ! dv_ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat,n_iter):idem

  integer :: ldim

  real (kind=DP) :: betamix(maxmix,maxmix), gamma0, work(maxmix)

  logical ::    &
                saveonfile, &! save intermediate steps on file "filename"
                opnd,       &! if true the file is already opened
                exst         ! if true the file exists

  real (kind=DP), external :: rho_dot_product_nc, ns_dot_product_nc

  call start_clock('mix_rho')

  if (iter < 1) call errore('mix_rho_nc','iter is wrong',1)
  if (n_iter > maxmix) call errore('mix_rho_nc','n_iter too big',1)
  if (lda_plus_u) ldim = 2 * Hubbard_lmax + 1

  saveonfile = ( filename /= ' ' )
  !
  allocate(rhocin(ngm0,nspin), rhocout(ngm0,nspin))
  !
  ! psic is used as work space - must be already allocated !
  !
  DO is=1,nspin
     !
     psic(:,1) = DCMPLX (rhoin(:,is), 0.d0)
     !
     call cft3(psic(1,1),nr1,nr2,nr3,nrx1,nrx2,nrx3,-1)
     !
     rhocin(:,is) = psic(nl(:),1)
     !
     psic(:,1) = DCMPLX (rhout(:,is), 0.d0)
     !
     call cft3( psic(1,1), nr1, nr2, nr3, nrx1, nrx2, nrx3, -1 )
     !
     rhocout(:,is) = psic(nl(:),1) - rhocin(:,is)
     !
  END DO
  !
  IF (lda_plus_u) nsout(:,:,:,:) = nsout(:,:,:,:) - nsin(:,:,:,:)
  !
  dr2 = rho_dot_product_nc(rhocout,rhocout) + ns_dot_product_nc(nsout,nsout) 
  !
  ! ... if the self-consistency error (dr2) is smaller than the estimated 
  ! ... error due to diagonalization (tr2_min), exit and leave rhoin and 
  ! ... rhout unchanged
  !
  IF ( dr2 < tr2_min ) THEN
     !
     DEALLOCATE( rhocin, rhocout )
     !
     CALL stop_clock( 'mix_rho' )
     !
     RETURN
     !
  END IF
  !
  conv = (dr2 < tr2)
  !
  IF (saveonfile) THEN
     !
     iunmix = find_free_unit()
     !
     IF (lda_plus_u) iunmix2 = find_free_unit()
     !
     IF (conv) THEN
        !
        CALL diropn (iunmix, filename, ( 2 * ngm0 * nspin ), exst)
        CLOSE( UNIT = iunmix, STATUS = 'DELETE')
        !
        IF ( lda_plus_u ) THEN
           !
           CALL diropn( iunmix2, TRIM( filename ) // '.ns', &
                        ( ldim * ldim * nspin * nat ), exst )
           CLOSE( UNIT = iunmix2, STATUS = 'DELETE' )
           !
        END IF
        !
        DEALLOCATE (rhocin, rhocout)
        !
        CALL stop_clock('mix_rho')
        !
        RETURN
        !
     END IF
     !
     CALL diropn( iunmix, filename, ( 2 * ngm0 * nspin ), exst )
     !
     IF ( lda_plus_u ) &
        CALL diropn( iunmix2, TRIM( filename ) // '.ns', &
                     ( ldim * ldim * nspin * nat ), exst )
     !
     IF ( iter > 1 .AND. .NOT. exst ) THEN
        !
        call errore('mix_rho_nc','file not found, restarting',-1)
        iter=1
        !
     END IF
     !
     ALLOCATE (df(ngm0,nspin,n_iter), dv(ngm0,nspin,n_iter))
     IF (lda_plus_u) &
        ALLOCATE (df_ns(ldim,ldim,nspin,nat,n_iter), &
                  dv_ns(ldim,ldim,nspin,nat,n_iter))
  ELSE
     !
     IF (iter == 1) THEN
        !
        ALLOCATE (df(ngm0,nspin,n_iter), dv(ngm0,nspin,n_iter))
        !
        IF (lda_plus_u) &
           ALLOCATE (df_ns(ldim,ldim,nspin,nat,n_iter),&
                     dv_ns(ldim,ldim,nspin,nat,n_iter))
     END If
     !
     IF (conv) THEN
        !
        IF (lda_plus_u) DEALLOCATE(df_ns, dv_ns)
        !
        DEALLOCATE (df, dv)
        !
        DEALLOCATE (rhocin, rhocout)
        !
        CALL stop_clock('mix_rho')
        !
        RETURN
        !
     END If
     !
     ALLOCATE (rhoinsave(ngm0,nspin), rhoutsave(ngm0,nspin))
     !
     IF (lda_plus_u) &
        ALLOCATE(nsinsave (ldim,ldim,nspin,nat), &
                 nsoutsave(ldim,ldim,nspin,nat))
  END IF
  !
  ! copy only the high frequency Fourier component into rhoin
  !
  DO is=1,nspin
     !
     psic(:,1) = ZERO
     !
     psic(nl(:),1) = rhocin(:,is) + rhocout(:,is)
     !
     CALL cft3(psic(1,1),nr1,nr2,nr3,nrx1,nrx2,nrx3,+1)
     !
     rhoin(:,is) = rhout(:,is) - psic(:,1)
     !
  END DO
  !
  ! iter_used = iter-1  if iter <= n_iter
  ! iter_used = n_iter  if iter >  n_iter
  !
  iter_used = MIN( ( iter - 1), n_iter )
  !
  ! ipos is the position in which results from the present iteration
  ! are stored. ipos=iter-1 until ipos=n_iter, then back to 1,2,...
  !
  ipos = iter - 1 - ( ( iter - 2 ) / n_iter ) * n_iter
  !
  IF ( iter > 1 ) THEN
     !
     IF (saveonfile) THEN
        !
        CALL davcio( df(1,1,ipos), 2*ngm0*nspin, iunmix, 1, -1 )
        CALL davcio( dv(1,1,ipos), 2*ngm0*nspin, iunmix, 2, -1 )
        !
        IF (lda_plus_u) THEN
           !
           CALL davcio(df_ns(1,1,1,1,ipos),ldim*ldim*nspin*nat,iunmix2,1,-1)
           CALL davcio(dv_ns(1,1,1,1,ipos),ldim*ldim*nspin*nat,iunmix2,2,-1)
           !
        END IF
        !
     END If
     !
     df(:,:,ipos) = df(:,:,ipos) - rhocout(:,:)
     dv(:,:,ipos) = dv(:,:,ipos) - rhocin (:,:)
     !
     IF (lda_plus_u) THEN
        !
        df_ns(:,:,:,:,ipos) = df_ns(:,:,:,:,ipos) - nsout
        dv_ns(:,:,:,:,ipos) = dv_ns(:,:,:,:,ipos) - nsin
        !
     END IF
     !
  END If
  !
  IF ( saveonfile ) THEN
     !
     DO i = 1, iter_used
        !
        IF ( i /= ipos) THEN
           !
           CALL davcio( df(1,1,i), 2*ngm0*nspin, iunmix, 2*i+1, -1 )
           CALL davcio( dv(1,1,i), 2*ngm0*nspin, iunmix, 2*i+2, -1 )
           !
           IF (lda_plus_u) THEN
              !
              CALL davcio(df_ns(1,1,1,1,i),ldim*ldim*nspin*nat,iunmix2,2*i+1,-1)
              CALL davcio(dv_ns(1,1,1,1,i),ldim*ldim*nspin*nat,iunmix2,2*i+2,-1)
              !
           END IF
           !
        END IF
        !
     END DO
     !
     CALL davcio( rhocout, 2*ngm0*nspin, iunmix, 1, 1 )
     CALL davcio( rhocin , 2*ngm0*nspin, iunmix, 2, 1 )
     !
     IF ( iter > 1) THEN
        !
        CALL davcio( df(1,1,ipos), 2*ngm0*nspin, iunmix, 2*ipos+1, 1 )
        CALL davcio( dv(1,1,ipos), 2*ngm0*nspin, iunmix, 2*ipos+2, 1 )
        !
     END IF
     !
     IF ( lda_plus_u ) THEN
        !
        CALL davcio( nsout, ldim*ldim*nspin*nat, iunmix2, 1, 1 )
        CALL davcio( nsin , ldim*ldim*nspin*nat, iunmix2, 2, 1 )
        !
        IF ( iter > 1) THEN
           !
           CALL davcio( df_ns(1,1,1,1,ipos), ldim*ldim*nspin*nat, &
                        iunmix2, 2*ipos+1, 1 )
           CALL davcio( dv_ns(1,1,1,1,ipos), ldim*ldim*nspin*nat, &
                        iunmix2, 2*ipos+2, 1 )
        END IF
        !
     END IF
     !
  ELSE
     !
     rhoinsave = rhocin
     rhoutsave = rhocout
     !
     IF (lda_plus_u) THEN
        !
        nsinsave  = nsin
        nsoutsave = nsout
        !
     END IF
     !
  END IF
  !
  DO i = 1, iter_used
     !
     DO j = i, iter_used
        !
        betamix(i,j) = rho_dot_product_nc( df(1,1,j), df(1,1,i) )  
        !
        IF ( lda_plus_u ) &
           betamix(i,j) = betamix(i,j) + &
                     ns_dot_product_nc( df_ns(1,1,1,1,j), df_ns(1,1,1,1,i) )
        !
     END DO
     !
  END DO
  !
  call DSYTRF ('U',iter_used,betamix,maxmix,iwork,work,maxmix,info)
  call errore('mix_rho_nc','factorization',info)
  !
  call DSYTRI ('U',iter_used,betamix,maxmix,iwork,work,info)
  call errore('mix_rho_nc','DSYTRI',info)
  !
  do i=1,iter_used
     do j=i+1,iter_used
        betamix(j,i)=betamix(i,j)
     end do
  end do
  !
  DO i = 1, iter_used
     !
     work(i) = rho_dot_product_nc( df(1,1,i), rhocout )
     !
     IF ( lda_plus_u ) &
        work(i) = work(i) + ns_dot_product_nc( df_ns(1,1,1,1,i), nsout )
     !
  END DO
  !
  do i=1,iter_used
     !
     gamma0 = SUM( betamix(1:iter_used,i) * work(1:iter_used) )
     !
     rhocin (:,:) = rhocin (:,:) - gamma0 * dv(:,:,i)
     rhocout(:,:) = rhocout(:,:) - gamma0 * df(:,:,i)
     !
     IF (lda_plus_u) THEN
        !
        nsin  = nsin  - gamma0 * dv_ns(:,:,:,:,i)
        nsout = nsout - gamma0 * df_ns(:,:,:,:,i)
        !
     END IF
     !
  END DO
  !
  ! - auxiliary vectors dv and df not needed anymore
  !
  IF (saveonfile) THEN
     !
     IF (lda_plus_u) THEN
        !
        CLOSE( iunmix2, STATUS = 'KEEP' )
        !
        DEALLOCATE (df_ns, dv_ns)
        !
     END If
     !
     CLOSE( iunmix, STATUS = 'KEEP' )
     !
     DEALLOCATE (df, dv)
     !
  ELSE
     !
     inext = iter - ( ( iter - 1 ) / n_iter ) * n_iter
     !
     IF (lda_plus_u) THEN
        !
        df_ns(:,:,:,:,inext) = nsoutsave
        dv_ns(:,:,:,:,inext) = nsinsave
        !
        DEALLOCATE( nsinsave, nsoutsave )
        !
     END IF
     !
     df(:,:,inext) = rhoutsave(:,:)
     dv(:,:,inext) = rhoinsave(:,:)
     !
     DEALLOCATE( rhoinsave, rhoutsave )

  END IF

  ! - preconditioning the new search direction (if imix.gt.0)

  IF ( imix == 1 ) THEN
     !
     CALL approx_screening_nc(rhocout)
     !
  ELSE IF ( imix == 2 ) THEN
     !
     CALL approx_screening2_nc(rhocout,rhocin)
     !
  END IF
  !
  ! - set new trial density
  !
  rhocin = rhocin + alphamix * rhocout
  !
  IF ( lda_plus_u ) nsin = nsin + alphamix * nsout
  !
  ! ... back ro real space
  !
  do is=1,nspin
     !
     psic(:,1) = ZERO
     !
     psic(nl(:),1) = rhocin(:,is)
     !
     CALL cft3( psic(1,1), nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
     !
     rhoin(:,is) = rhoin(:,is) + DBLE( psic(:,1) )

  END DO
  !
  ! - clean up
  !
  DEALLOCATE( rhocout )
  DEALLOCATE( rhocin )
  CALL stop_clock( 'mix_rho' )

  RETURN
END SUBROUTINE mix_rho_nc

!
!--------------------------------------------------------------------
FUNCTION rho_dot_product_nc (rho1,rho2)
  !--------------------------------------------------------------------
  ! this function evaluates the dot product between two input densities
  !
  USE kinds, only : DP
  use pwcom
  USE control_flags,      ONLY : ngm0
  implicit none
  !
  ! I/O variables
  !
  real (kind=DP) :: rho_dot_product_nc ! (out) the function value

  complex (kind=DP), intent(in) :: rho1(ngm0,nspin), rho2(ngm0,nspin) ! (in) the two densities

  !
  ! and the local variables
  !
  real (kind=DP) :: fac   ! a multiplicative factors

  integer  :: is, ig

  if (nspin.ne.4)  call errore &
      ('rho_dot_product_nc','it shoud be used only in the noncolinear case ',4)

  rho_dot_product_nc = 0.d0

  do ig = gstart,ngm0
     fac = e2*fpi / (tpiba2*gg(ig))
     rho_dot_product_nc = rho_dot_product_nc +  fac * &
                          DREAL( conjg( rho1(ig,1))*(rho2(ig,1) ))
  end do

  fac = e2*fpi / (tpi**2)  ! lambda=1 a.u.
  do ig = 1,ngm0
     rho_dot_product_nc = rho_dot_product_nc +  fac * &
                        ( DREAL( conjg( rho1(ig,2))*(rho2(ig,2) )) + &
                          DREAL( conjg( rho1(ig,3))*(rho2(ig,3) )) + &
                          DREAL( conjg( rho1(ig,4))*(rho2(ig,4) )) )
  end do

  rho_dot_product_nc = rho_dot_product_nc * omega / 2.d0
  !
#ifdef __PARA
  call reduce(1,rho_dot_product_nc)
#endif

  return
end function rho_dot_product_nc

!
!--------------------------------------------------------------------
function ns_dot_product_nc (ns1,ns2)
  !--------------------------------------------------------------------
  ! this function evaluates the dot product between two input densities
  !
  USE kinds, only : DP
  use pwcom
  USE ions_base,  ONLY : nat, ityp
  !
  implicit none
  !
  ! I/O variables
  !
  real (kind=DP) :: ns_dot_product_nc ! (out) the function value

  real (kind=DP) :: ns1(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat), &
                    ns2(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat) 
                    ! (in) the two ns 
  !
  ! and the local variables
  !
  real (kind=DP) :: sum
  integer  :: na, nt, is, m1, m2

  ns_dot_product_nc = 0.d0
  if (.not. lda_plus_u ) return

  call errore('ns_dot_product_nc','not tested',1)
  do na = 1, nat
     nt = ityp (na)
     if (Hubbard_U(nt).ne.0.d0 .or. Hubbard_alpha(nt).ne.0.d0) then
        sum =0.d0
        do is = 1,nspin
           do m1 = 1, 2 * Hubbard_l(nt) + 1
              do m2 = m1, 2 * Hubbard_l(nt) + 1
                 sum = sum + ns1(m1,m2,is,na)*ns2(m2,m1,is,na)
              enddo
           enddo
        end do
        ns_dot_product_nc = ns_dot_product_nc + 0.5d0*Hubbard_U(nt) * sum
     endif
  end do
  if (nspin.eq.1) ns_dot_product_nc = 2.d0 * ns_dot_product_nc

  return
end function ns_dot_product_nc

!--------------------------------------------------------------------
subroutine approx_screening_nc (drho)
  !--------------------------------------------------------------------
  ! apply an average TF preconditioning to drho
  !
  USE kinds, only : DP
  use pwcom
  USE control_flags,      ONLY : ngm0
  implicit none
  !
  ! I/O
  !
  complex (kind=DP) drho(ngm0,nspin) ! (in/out)
  !
  ! and the local variables
  !
  real (kind=DP) :: rrho, rmag, rs, agg0

  integer :: is, ig

  rs = (3.d0*omega/fpi/nelec)**(1.d0/3.d0)
  agg0 = (12.d0/pi)**(2.d0/3.d0)/tpiba2/rs

#ifdef DEBUG
  WRITE( stdout,'(a,f12.6,a,f12.6)') ' avg rs  =', rs, ' avg rho =', nelec/omega
#endif

  if (nspin == 1) then
     is = 1
     do ig = 1,ngm0
        drho(ig,is) =  drho(ig,is) * gg(ig)/(gg(ig)+agg0)
     end do
  else
     do ig = 1,ngm0
        rrho = (drho(ig,1) + drho(ig,2)) * gg(ig)/(gg(ig)+agg0)
        rmag = (drho(ig,1) - drho(ig,2))
        drho(ig,1) =  0.5d0 * (rrho + rmag)
        drho(ig,2) =  0.5d0 * (rrho - rmag)
     end do
  end if

  return
end subroutine approx_screening_nc

!
!--------------------------------------------------------------------
  subroutine approx_screening2_nc (drho,rhobest)
  !--------------------------------------------------------------------
  ! apply a local-density dependent TF preconditioning to drho
  !
  USE kinds, only : DP
  USE wavefunctions_module, ONLY : psic => psic_nc
  use pwcom
  USE control_flags,      ONLY : ngm0
  !
  ! I/O
  !
  implicit none
  complex (kind=DP) ::  drho(ngm0,nspin), rhobest(ngm0,nspin)
  !
  !    and the local variables
  !
  integer, parameter :: mmx = 12
  integer :: iwork(mmx),i,j,m,info, nspin_save
  real (kind=DP) :: rs, min_rs, max_rs, avg_rsm1, target, &
                    dr2_best, ccc, cbest, l2smooth
  real (kind=DP) :: aa(mmx,mmx), invaa(mmx,mmx), bb(mmx), work(mmx), &
                    vec(mmx),agg0
  complex (kind=DP) :: rrho, rmag

  complex (kind=DP), allocatable :: v(:,:), w(:,:), dv(:), &
                                vbest(:), wbest(:)
  ! v(ngm0,mmx), w(ngm0,mmx), dv(ngm0), vbest(ngm0), wbest(ngm0)
  real (kind=DP), allocatable :: alpha(:)
  ! alpha(nrxx)

  integer :: is, ir, ig

  real (kind=DP), external :: rho_dot_product_nc

  if (nspin == 2) then
     do ig=1,ngm0
        rrho       = drho(ig,1) + drho(ig,2)
        rmag       = drho(ig,1) - drho(ig,2)
        drho(ig,1) = rrho
        drho(ig,2) = rmag
     end do
  end if

  nspin_save = nspin
  nspin = 1
  is = 1
  target = 0.d0

!  WRITE( stdout,*) ' eccoci qua '

  if (gg(1) < 1.d-8) drho(1,is) = (0.d0,0.d0)

  allocate (alpha(nrxx), v(ngm0,mmx), w(ngm0,mmx), &
            dv(ngm0), vbest(ngm0), wbest(ngm0))

  v(:,:) = (0.d0,0.d0)
  w(:,:) = (0.d0,0.d0)
  dv(:) = (0.d0,0.d0)
  vbest(:)= (0.d0,0.d0)
  wbest(:)= (0.d0,0.d0)

  !
  ! - calculate alpha from density smoothed with a lambda=0 a.u.
  !
  l2smooth = 0.d0
  psic(:,1) = (0.d0,0.d0)
  if (nspin == 1) then
     do ig=1,ngm0
        psic(nl(ig),1) = rhobest(ig,1) * exp(-0.5*l2smooth*tpiba2*gg(ig))
     end do
  else
     do ig=1,ngm0
        psic(nl(ig),1) =(rhobest(ig,1) + rhobest(ig,2)) &
                                    * exp(-0.5*l2smooth*tpiba2*gg(ig))
     end do
  end if
  call cft3(psic(1,1),nr1,nr2,nr3,nrx1,nrx2,nrx3,+1)
  alpha(:) = real(psic(:,1))

  min_rs = (3.d0*omega/fpi/nelec)**(1.d0/3.d0)
  max_rs = min_rs
  avg_rsm1 = 0.d0

  do ir=1,nrxx

     alpha(ir)=abs(alpha(ir))
     rs = (3.d0/fpi/alpha(ir))**(1.d0/3.d0)

     min_rs = min(min_rs,rs)
     avg_rsm1 =avg_rsm1 + 1.d0/rs
     max_rs = max(max_rs,rs)

     alpha(ir) = rs

  end do

#ifdef __PARA
  call reduce  (1, avg_rsm1)
  call extreme (min_rs, -1)
  call extreme (max_rs, +1)
#endif

  call DSCAL(nrxx, 3.d0*(tpi/3.d0)**(5.d0/3.d0), alpha, 1)

  avg_rsm1 = (nr1*nr2*nr3)/avg_rsm1
  rs = (3.d0*omega/fpi/nelec)**(1.d0/3.d0)
  agg0 = (12.d0/pi)**(2.d0/3.d0)/tpiba2/avg_rsm1
#ifdef DEBUG
  WRITE( stdout,'(a,5f12.6)') ' min/avgm1/max rs  =', min_rs,avg_rsm1,max_rs,rs
#endif

  !
  ! - calculate deltaV and the first correction vector
  !
  psic(:,1) = (0.d0,0.d0)
  do ig=1,ngm0
     psic(nl(ig),1) = drho(ig,is)
  end do
  call cft3(psic(1,1),nr1,nr2,nr3,nrx1,nrx2,nrx3,+1)
  do ir=1,nrxx
    psic(ir,1) = psic(ir,1) * alpha(ir)
  end do
  call cft3(psic(1,1),nr1,nr2,nr3,nrx1,nrx2,nrx3,-1)
  do ig=1,ngm0
     dv(ig) = psic(nl(ig),1)*gg(ig)*tpiba2
     v(ig,1)= psic(nl(ig),1)*gg(ig)/(gg(ig)+agg0)
  end do
  m=1
  ccc = rho_dot_product_nc(dv,dv)
  aa(:,:) = 0.d0
  bb(:) = 0.d0

3 continue
  !
  ! - generate the vector w
  !
  do ig=1,ngm0
     w(ig,m) = gg(ig)*tpiba2*v(ig,m)
  end do
  psic(:,1) = (0.d0,0.d0)
  do ig=1,ngm0
     psic(nl(ig),1) = v(ig,m)
  end do
  call cft3(psic(1,1),nr1,nr2,nr3,nrx1,nrx2,nrx3,+1)
  do ir=1,nrxx
     psic(ir,1) = psic(ir,1)*fpi*e2/alpha(ir)
  end do
  call cft3(psic(1,1),nr1,nr2,nr3,nrx1,nrx2,nrx3,-1)
  do ig=1,ngm0
     w(ig,m) = w(ig,m) + psic(nl(ig),1)
  end do

  !
  ! - build the linear system
  !
  do i=1,m
     aa(i,m) = rho_dot_product_nc(w(1,i),w(1,m))
     aa(m,i) = aa(i,m)
  end do
  bb(m) = rho_dot_product_nc(w(1,m),dv)

  !
  ! - solve it -> vec
  !
  call DCOPY (mmx*mmx,aa,1,invaa,1)
  call DSYTRF ('U',m,invaa,mmx,iwork,work,mmx,info)
  call errore('approx_screening2_nc','factorization',info)
  call DSYTRI ('U',m,invaa,mmx,iwork,work,info)
  call errore('approx_screening2_nc','DSYTRI',info)
  !
  do i=1,m
     do j=i+1,m
        invaa(j,i)=invaa(i,j)
     end do
  end do
  do i=1,m
     vec(i) = 0.d0
     do j=1,m
        vec(i) = vec(i) + invaa(i,j)*bb(j)
     end do
  end do
  ! -
  vbest(:) = (0.d0,0.d0)
  wbest(:) = dv(:)
  do i=1,m
     call DAXPY(2*ngm0, vec(i), v(1,i),1, vbest,1)
     call DAXPY(2*ngm0,-vec(i), w(1,i),1, wbest,1)
  end do

  cbest = ccc
  do i=1,m
     cbest = cbest - bb(i)*vec(i)
  end do

  dr2_best= rho_dot_product_nc(wbest,wbest)
  if (target == 0.d0) target = 1.d-6 * dr2_best
!  WRITE( stdout,*) m, dr2_best, cbest

  if (dr2_best < target) then
!     WRITE( stdout,*) ' last', dr2_best/target * 1.d-6
     psic(:,1) = (0.d0,0.d0)
     do ig=1,ngm0
        psic(nl(ig),1) = vbest(ig)
     end do
     call cft3(psic(1,1),nr1,nr2,nr3,nrx1,nrx2,nrx3,+1)
     do ir=1,nrxx
        psic(ir,1) = psic(ir,1)/alpha(ir)
     end do
     call cft3(psic(1,1),nr1,nr2,nr3,nrx1,nrx2,nrx3,-1)
     do ig=1,ngm0
        drho(ig,is) = psic(nl(ig),1)
     end do
     nspin = nspin_save
     if (nspin == 2) then
        do ig=1,ngm0
           rrho = drho(ig,1)
           rmag = drho(ig,2)
           drho(ig,1) = 0.5d0 * ( rrho + rmag )
           drho(ig,2) = 0.5d0 * ( rrho - rmag )
        end do
     end if
     deallocate (alpha, v, w, dv, vbest, wbest)
     return
  else if (m >= mmx) then
!     WRITE( stdout,*) m, dr2_best, cbest
     m=1
     do ig=1,ngm0
        v(ig,m)=vbest(ig)
     end do
     aa(:,:) = 0.d0
     bb(:) = 0.d0
     go to 3
  end if

  m = m + 1
  do ig=1,ngm0
     v(ig,m)=wbest(ig)/(gg(ig)+agg0)
  end do

  go to 3

end subroutine approx_screening2_nc
