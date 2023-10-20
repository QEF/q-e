!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#if defined(__CUDA)
#define PINMEM ,PINNED 
#else
#define PINMEM
#endif

!=----------------------------------------------------------------------------=!
   MODULE wavefunctions
!=----------------------------------------------------------------------------=!
     !! Wavefunction arrays.
     !
     USE kinds, ONLY :  DP
#if defined (__CUDA)
     USE cudafor
#endif

     IMPLICIT NONE
     SAVE

     !
     COMPLEX(DP), ALLOCATABLE, TARGET :: &
       evc(:,:)
       !! wavefunctions in the PW basis set.  
       !! noncolinear case: first index is a combined PW + spin index
       !
#if defined(__CUDA)
       attributes(PINNED) :: evc
#endif
     !
     COMPLEX(DP) , ALLOCATABLE, TARGET :: psic(:)
     !! additional memory for FFT
     COMPLEX(DP) , ALLOCATABLE, TARGET :: psic_nc(:,:)
     !! additional memory for FFT for the noncolinear case
     !
     !
     ! electronic wave functions, CPV code
     ! distributed over gvector and bands
     !
!dir$ attributes align: 4096 :: c0_bgrp, cm_bgrp, phi
     COMPLEX(DP), ALLOCATABLE PINMEM :: c0_bgrp(:,:)  ! wave functions at time t
     COMPLEX(DP), ALLOCATABLE PINMEM :: cm_bgrp(:,:)  ! wave functions at time t-delta t
     COMPLEX(DP), ALLOCATABLE :: phi(:,:) ! |phi> = s'|c0> = |c0> + sum q_ij |i><j|c0>
     COMPLEX(DP), ALLOCATABLE :: c0_d(:,:)  ! wave functions at time t
     COMPLEX(DP), ALLOCATABLE :: cm_d(:,:)  ! wave functions at time t-delta t
#if defined (__CUDA)
     ATTRIBUTES(DEVICE) :: c0_d, cm_d, phi
#endif
     ! for hybrid functionals in CP with Wannier functions
     COMPLEX(DP), ALLOCATABLE :: cv0(:,:) ! Lingzhu Kong

   CONTAINS

      SUBROUTINE deallocate_wavefunctions
       IF( ALLOCATED( cv0) ) DEALLOCATE( cv0)   ! Lingzhu Kong
       IF( ALLOCATED( c0_bgrp ) ) DEALLOCATE( c0_bgrp )
       IF( ALLOCATED( cm_bgrp ) ) DEALLOCATE( cm_bgrp )
       IF( ALLOCATED( phi ) ) DEALLOCATE( phi )
       IF( ALLOCATED( psic_nc ) ) DEALLOCATE( psic_nc )
       IF( ALLOCATED( psic ) ) DEALLOCATE( psic )
       IF( ALLOCATED( evc ) ) DEALLOCATE( evc )
#if defined (__CUDA)
       IF( ALLOCATED( c0_d ) ) DEALLOCATE( c0_d )
       IF( ALLOCATED( cm_d ) ) DEALLOCATE( cm_d )
#endif
     END SUBROUTINE deallocate_wavefunctions

     SUBROUTINE allocate_cp_wavefunctions( ngw, nbspx, vnbsp, lwfpbe0nscf )
       INTEGER, INTENT(IN) :: ngw, nbspx, vnbsp
       LOGICAL, INTENT(IN) :: lwfpbe0nscf
       INTEGER :: ierr
       ALLOCATE( c0_bgrp( ngw, nbspx ), STAT=ierr )
       IF( ierr /= 0 ) &
         CALL errore( ' allocate_cp_wavefunctions ', ' allocating on CPU ', ABS( ierr ) )
       c0_bgrp = (0_DP,0_DP)
       ALLOCATE( cm_bgrp( ngw, nbspx ), STAT=ierr )
       IF( ierr /= 0 ) &
         CALL errore( ' allocate_cp_wavefunctions ', ' allocating on CPU ', ABS( ierr ) )
       cm_bgrp = (0_DP,0_DP)
       ALLOCATE( phi( ngw, nbspx ), STAT=ierr )
       IF( ierr /= 0 ) &
         CALL errore( ' allocate_cp_wavefunctions ', ' allocating on CPU ', ABS( ierr ) )
       phi = (0_DP,0_DP)
       IF(lwfpbe0nscf) THEN
         ALLOCATE(cv0( ngw, vnbsp ), STAT=ierr )   ! Lingzhu Kong
         IF( ierr /= 0 ) &
           CALL errore( ' allocate_cp_wavefunctions ', ' allocating on CPU ', ABS( ierr ) )
         cv0 = (0_DP,0_DP)
       END IF
#if defined (__CUDA)
       ALLOCATE( c0_d( ngw, nbspx ), STAT=ierr )
       IF( ierr /= 0 ) &
         CALL errore( ' allocate_cp_wavefunctions ', ' allocating on GPU ', ABS( ierr ) )
       ALLOCATE( cm_d( ngw, nbspx ), STAT=ierr )
       IF( ierr /= 0 ) &
         CALL errore( ' allocate_cp_wavefunctions ', ' allocating on GPU ', ABS( ierr ) )
#endif
     END SUBROUTINE

!=----------------------------------------------------------------------------=!
   END MODULE wavefunctions
!=----------------------------------------------------------------------------=!

   module atomic_wfc_init

      implicit none
   
   contains
   
      SUBROUTINE atomic_wfc_nog(ik, wfcatom, omega, tpiba, nat, ntyp, ityp, tau, natomwfc, &
         mill, eigts1_, eigts2_, eigts3_, g, xk, igk_k, ngk, npwx, &
         upf, nwfcm, noncolin, domag, npol, angle1, angle2, starting_spin_angle, &
         rot_ylm)
   
   
         USE kinds,            ONLY : DP
         USE pseudo_types,     ONLY : pseudo_upf ! type
         USE constants,        ONLY : tpi, fpi, pi ! parameter
         USE upf_spinorb,      ONLY :  lmaxx ! parameter
   
   
         implicit none
         ! old global variables
         REAL(DP), INTENT(IN) :: omega, tpiba, tau(:,:), g(:,:), xk(:,:), angle1(:), angle2(:)
         integer, intent(IN) :: nat, ntyp, ityp(:), natomwfc, mill(:,:), igk_k(:,:), ngk(:), &
            npwx, nwfcm, npol
         complex(DP), intent(in),target, contiguous :: eigts1_(:,:), eigts2_(:,:), eigts3_(:,:)
         complex(dp), intent(in) :: rot_ylm(:,:)
         complex(dp), pointer :: eigts1(:,:), eigts2(:,:), eigts3(:,:)
         type(pseudo_upf), intent(in) :: upf(:)
         logical, intent(in) :: noncolin, domag, starting_spin_angle
   
         INTEGER, INTENT(IN) :: ik
         !! k-point index
         COMPLEX(DP), INTENT(INOUT) :: wfcatom( npwx, npol, natomwfc )
         !! Superposition of atomic wavefunctions
         !
         ! ... local variables
         !
         INTEGER :: eigts_shape(2)
         INTEGER :: n_starting_wfc, lmax_wfc, nt, l, nb, na, m, lm, ig, iig, &
            i0, i1, i2, i3, npw
         REAL(DP),    ALLOCATABLE :: qg(:), ylm (:,:), chiq (:,:,:), gk (:,:)
         COMPLEX(DP), ALLOCATABLE :: sk (:), aux(:)
         COMPLEX(DP) :: kphase, lphase
         REAL(DP)    :: arg
   
   
         CALL start_clock( 'atomic_wfc' )
   
         !Fortran has problems to transmit the information about array's bounds. Here we use a pointer to walk around the issue
         eigts_shape = shape(eigts1_)
         eigts1(-eigts_shape(1)/2:eigts_shape(1)/2,1:eigts_shape(2) )=>eigts1_
         eigts_shape = shape(eigts2_)
         eigts2(-eigts_shape(1)/2:eigts_shape(1)/2,1:eigts_shape(2) )=>eigts2_
         eigts_shape = shape(eigts3_)
         eigts3(-eigts_shape(1)/2:eigts_shape(1)/2,1:eigts_shape(2) )=>eigts3_
   
         ! calculate max angular momentum required in wavefunctions
         lmax_wfc = 0
         DO nt = 1, ntyp
            lmax_wfc = MAX( lmax_wfc, MAXVAL( upf(nt)%lchi(1:upf(nt)%nwfc) ) )
         END DO
         !
         npw = ngk(ik)
         !
         ALLOCATE( ylm (npw,(lmax_wfc+1)**2), chiq(npw,nwfcm,ntyp), &
            gk(3,npw), qg(npw) )
         !
         DO ig = 1, npw
            iig = igk_k (ig,ik)
            gk (1,ig) = xk(1, ik) + g(1,iig)
            gk (2,ig) = xk(2, ik) + g(2,iig)
            gk (3,ig) = xk(3, ik) + g(3,iig)
            qg(ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
         END DO
         !
         !  ylm = spherical harmonics
         !
         CALL ylmr2( (lmax_wfc+1)**2, npw, gk, qg, ylm )
         !
         ! set now q=|k+G| in atomic units
         !
         DO ig = 1, npw
            qg(ig) = SQRT( qg(ig) )*tpiba
         END DO
         !
         CALL interp_atwfc ( npw, qg, nwfcm, chiq )
         !
         DEALLOCATE( qg, gk )
         ALLOCATE( aux(npw), sk(npw) )
         !
         wfcatom(:,:,:) = (0.0_dp, 0.0_dp)
         n_starting_wfc = 0
         !
         DO na = 1, nat
            arg = (xk(1,ik)*tau(1,na) + xk(2,ik)*tau(2,na) + xk(3,ik)*tau(3,na)) * tpi
            kphase = CMPLX( COS(arg), - SIN(arg) ,KIND=DP)
            !
            !     sk is the structure factor
            !
            DO ig = 1, npw
               iig = igk_k (ig,ik)
               sk (ig) = kphase * eigts1 (mill (1,iig), na) * &
                  eigts2 (mill (2,iig), na) * &
                  eigts3 (mill (3,iig), na)
            END DO
            !
            nt = ityp (na)
            DO nb = 1, upf(nt)%nwfc
               IF ( upf(nt)%oc(nb) >= 0.d0 ) THEN
                  l = upf(nt)%lchi(nb)
                  lphase = (0.d0,1.d0)**l
                  !
                  !  the factor i^l MUST BE PRESENT in order to produce
                  !  wavefunctions for k=0 that are real in real space
                  !
                  IF ( noncolin ) THEN
                     !
                     IF ( upf(nt)%has_so ) THEN
                        !
                        IF (starting_spin_angle.OR..NOT.domag) THEN
                           CALL atomic_wfc_so( )
                        ELSE
                           CALL atomic_wfc_so_mag( )
                        END IF
                        !
                     ELSE
                        !
                        CALL atomic_wfc_nc( )
                        !
                     END IF
                     !
                  ELSE
                     !
                     CALL atomic_wfc___( )
                     !
                  END IF
                  !
               END IF
               !
            END DO
            !
         END DO
   
         IF ( n_starting_wfc /= natomwfc) call errore ('atomic_wfc', &
            'internal error: some wfcs were lost ', 1 )
   
         DEALLOCATE( aux, sk, chiq, ylm )
   
         CALL stop_clock( 'atomic_wfc' )
   
         RETURN
   
      CONTAINS
         !----------------------------------------------------------------
         SUBROUTINE atomic_wfc_so( )
            !------------------------------------------------------------
            !! Spin-orbit case.
            !
            implicit none
            REAL(DP) :: fact(2), j
            REAL(DP), EXTERNAL :: spinor
            INTEGER :: ind, ind1, n1, is, sph_ind
            !
            j = upf(nt)%jchi(nb)
            DO m = -l-1, l
               fact(1) = spinor(l,j,m,1)
               fact(2) = spinor(l,j,m,2)
               IF ( ABS(fact(1)) > 1.d-8 .OR. ABS(fact(2)) > 1.d-8 ) THEN
                  n_starting_wfc = n_starting_wfc + 1
                  IF (n_starting_wfc > natomwfc) CALL errore &
                     ('atomic_wfc_so', 'internal error: too many wfcs', 1)
                  DO is=1,2
                     IF (abs(fact(is)) > 1.d-8) THEN
                        ind=lmaxx+1+sph_ind(l,j,m,is)
                        aux=(0.d0,0.d0)
                        DO n1=1,2*l+1
                           ind1=l**2+n1
                           if (abs(rot_ylm(ind,n1)) > 1.d-8) &
                              aux(:)=aux(:)+rot_ylm(ind,n1)*ylm(:,ind1)
                        ENDDO
                        do ig = 1, npw
                           wfcatom(ig,is,n_starting_wfc) = lphase*fact(is)*&
                              sk(ig)*aux(ig)*chiq (ig, nb, nt)
                        END DO
                     ELSE
                        wfcatom(:,is,n_starting_wfc) = (0.d0,0.d0)
                     END IF
                  END DO
               END IF
            END DO
            !
         END SUBROUTINE atomic_wfc_so
         !
         SUBROUTINE atomic_wfc_so_mag( )
            !
            !! Spin-orbit case, magnetization along "angle1" and "angle2"
            !! In the magnetic case we always assume that magnetism is much larger
            !! than spin-orbit and average the wavefunctions at l+1/2 and l-1/2
            !! filling then the up and down spinors with the average wavefunctions,
            !! according to the direction of the magnetization, following what is
            !! done in the noncollinear case.
            !
            implicit none
            REAL(DP) :: alpha, gamman, j
            COMPLEX(DP) :: fup, fdown
            REAL(DP), ALLOCATABLE :: chiaux(:)
            INTEGER :: nc, ib
            !
            j = upf(nt)%jchi(nb)
            !
            !  This routine creates two functions only in the case j=l+1/2 or exit in the
            !  other case
            !
            IF (ABS(j-l+0.5_DP)<1.d-4) RETURN
   
            ALLOCATE(chiaux(npw))
            !
            !  Find the functions j=l-1/2
            !
            IF (l == 0)  THEN
               chiaux(:)=chiq(:,nb,nt)
            ELSE
               DO ib=1, upf(nt)%nwfc
                  IF ((upf(nt)%lchi(ib) == l).AND. &
                     (ABS(upf(nt)%jchi(ib)-l+0.5_DP)<1.d-4)) THEN
                     nc=ib
                     EXIT
                  ENDIF
               ENDDO
               !
               !  Average the two functions
               !
               chiaux(:)=(chiq(:,nb,nt)*(l+1.0_DP)+chiq(:,nc,nt)*l)/(2.0_DP*l+1.0_DP)
               !
            ENDIF
            !
            !  and construct the starting wavefunctions as in the noncollinear case.
            !
            alpha = angle1(nt)
            gamman = - angle2(nt) + 0.5d0*pi
            !
            DO m = 1, 2 * l + 1
               lm = l**2 + m
               n_starting_wfc = n_starting_wfc + 1
               IF ( n_starting_wfc + 2*l+1 > natomwfc ) CALL errore &
                  ('atomic_wfc_so_mag', 'internal error: too many wfcs', 1)
               DO ig = 1, npw
                  aux(ig) = sk(ig)*ylm(ig,lm)*chiaux(ig)
               END DO
               !
               ! now, rotate wfc as needed
               ! first : rotation with angle alpha around (OX)
               !
               DO ig = 1, npw
                  fup = cos(0.5d0*alpha)*aux(ig)
                  fdown = (0.d0,1.d0)*sin(0.5d0*alpha)*aux(ig)
                  !
                  ! Now, build the orthogonal wfc
                  ! first rotation with angle (alpha+pi) around (OX)
                  !
                  wfcatom(ig,1,n_starting_wfc) = (cos(0.5d0*gamman) &
                     +(0.d0,1.d0)*sin(0.5d0*gamman))*fup
                  wfcatom(ig,2,n_starting_wfc) = (cos(0.5d0*gamman) &
                     -(0.d0,1.d0)*sin(0.5d0*gamman))*fdown
                  !
                  ! second: rotation with angle gamma around (OZ)
                  !
                  ! Now, build the orthogonal wfc
                  ! first rotation with angle (alpha+pi) around (OX)
                  !
                  fup = cos(0.5d0*(alpha+pi))*aux(ig)
                  fdown = (0.d0,1.d0)*sin(0.5d0*(alpha+pi))*aux(ig)
                  !
                  ! second, rotation with angle gamma around (OZ)
                  !
                  wfcatom(ig,1,n_starting_wfc+2*l+1) = (cos(0.5d0*gamman) &
                     +(0.d0,1.d0)*sin(0.5d0 *gamman))*fup
                  wfcatom(ig,2,n_starting_wfc+2*l+1) = (cos(0.5d0*gamman) &
                     -(0.d0,1.d0)*sin(0.5d0*gamman))*fdown
               END DO
            END DO
            n_starting_wfc = n_starting_wfc + 2*l+1
            DEALLOCATE( chiaux )
            !
         END SUBROUTINE atomic_wfc_so_mag
         !
         SUBROUTINE atomic_wfc_nc( )
            !
            !! noncolinear case, magnetization along "angle1" and "angle2"
            !
            implicit none
            REAL(DP) :: alpha, gamman
            COMPLEX(DP) :: fup, fdown
            !
            alpha = angle1(nt)
            gamman = - angle2(nt) + 0.5d0*pi
            !
            DO m = 1, 2 * l + 1
               lm = l**2 + m
               n_starting_wfc = n_starting_wfc + 1
               IF ( n_starting_wfc + 2*l+1 > natomwfc) CALL errore &
                  ('atomic_wfc_nc', 'internal error: too many wfcs', 1)
               DO ig = 1, npw
                  aux(ig) = sk(ig)*ylm(ig,lm)*chiq(ig,nb,nt)
               END DO
               !
               ! now, rotate wfc as needed
               ! first : rotation with angle alpha around (OX)
               !
               DO ig = 1, npw
                  fup = cos(0.5d0*alpha)*aux(ig)
                  fdown = (0.d0,1.d0)*sin(0.5d0*alpha)*aux(ig)
                  !
                  ! Now, build the orthogonal wfc
                  ! first rotation with angle (alpha+pi) around (OX)
                  !
                  wfcatom(ig,1,n_starting_wfc) = (cos(0.5d0*gamman) &
                     +(0.d0,1.d0)*sin(0.5d0*gamman))*fup
                  wfcatom(ig,2,n_starting_wfc) = (cos(0.5d0*gamman) &
                     -(0.d0,1.d0)*sin(0.5d0*gamman))*fdown
                  !
                  ! second: rotation with angle gamma around (OZ)
                  !
                  ! Now, build the orthogonal wfc
                  ! first rotation with angle (alpha+pi) around (OX)
                  !
                  fup = cos(0.5d0*(alpha+pi))*aux(ig)
                  fdown = (0.d0,1.d0)*sin(0.5d0*(alpha+pi))*aux(ig)
                  !
                  ! second, rotation with angle gamma around (OZ)
                  !
                  wfcatom(ig,1,n_starting_wfc+2*l+1) = (cos(0.5d0*gamman) &
                     +(0.d0,1.d0)*sin(0.5d0 *gamman))*fup
                  wfcatom(ig,2,n_starting_wfc+2*l+1) = (cos(0.5d0*gamman) &
                     -(0.d0,1.d0)*sin(0.5d0*gamman))*fdown
               END DO
            END DO
            n_starting_wfc = n_starting_wfc + 2*l+1
            !
         END SUBROUTINE atomic_wfc_nc
   
         SUBROUTINE atomic_wfc___( )
            !
            ! ... LSDA or nonmagnetic case
            !
            implicit none
            DO m = 1, 2 * l + 1
               lm = l**2 + m
               n_starting_wfc = n_starting_wfc + 1
               IF ( n_starting_wfc > natomwfc) CALL errore &
                  ('atomic_wfc___', 'internal error: too many wfcs', 1)
               !
               DO ig = 1, npw
                  wfcatom (ig, 1, n_starting_wfc) = lphase * &
                     sk (ig) * ylm (ig, lm) * chiq (ig, nb, nt)
               ENDDO
               !
            END DO
            !
         END SUBROUTINE atomic_wfc___
         !
      END SUBROUTINE atomic_wfc_nog
   
      SUBROUTINE atomic_wfc_cp(evc, omega, tpiba, nat, ntyp, ityp, tau, natomwfc, &
         mill, eigts1, eigts2, eigts3, g, iupdwn, &
         npwx, &
         upf, nwfcm, npol )
   
   
         USE kinds,        ONLY : DP
         USE pseudo_types, ONLY : pseudo_upf
         USE mp_global,    ONLY : intra_bgrp_comm
         use uspp_data,    ONLY : nqx, tab_at, dq
         use uspp_param,   ONLY : nsp
         USE gvecw,        ONLY : ecutwfc
   
         implicit none
         ! old global variables
         REAL(DP), INTENT(IN) :: omega, tpiba, tau(:,:), g(:,:)
         integer, intent(IN) :: nat, ntyp, ityp(:), natomwfc, mill(:,:), &
            npwx, nwfcm, npol, iupdwn(2)
         complex(DP), intent(in) :: eigts1(:,:), eigts2(:,:), eigts3(:,:)
         type(pseudo_upf), intent(in) :: upf(:)
   
         COMPLEX(DP), intent(inout) :: evc (:,:)
         COMPLEX(DP)  :: wfcatom( npwx, npol, natomwfc )
         !! Superposition of atomic wavefunctions
   
         !cp specific settings (gamma only)
         ! xk is 0,0,0, igk_k is the identical permutation, and ngk is 1
         ! wfcatom has 2 dimensions (npwx, nbnd*nspin)
         ! looks like only nspin=1,2 are implemented. The layout of the wfc is different:
         ! in cp it is the equivalent of (npwx, nbnd, nspin ), while in pw is (npwx, nspin, nbnd)
         real(dp) :: xk(3,1)
         integer :: ngk(1), igk_k(npwx,1), i, ik, ipol, sh(2)
   
         complex(dp), allocatable :: rot_ylm(:,:)
         REAL(DP), allocatable :: angle1(:), angle2(:)
         !! dummy, unused and not allocated variables
   
         ! identity permutation
         do i=1,npwx
            igk_k(i,1)=i
         enddo
         ! gamma point only
         xk=0.d0
         ngk(1) = npwx
         ik=1
         nqx = INT( (SQRT(ecutwfc) / dq + 4) * 1.d0 )
         allocate(tab_at(nqx,nwfcm,nsp))
         call init_tab_atwfc(omega, intra_bgrp_comm)
   
         ! only nospin / LSDA in CP
         call atomic_wfc_nog(ik, wfcatom, omega, tpiba, nat, ntyp, ityp, tau, natomwfc, &
            mill, eigts1, eigts2, eigts3, g, xk, igk_k, ngk, npwx, &
            upf, nwfcm, .false., .false., npol, angle1, angle2, .false., &
            rot_ylm)
   
         sh = shape(evc)
   
         !write the result in the correct order in evc
         do i=1,min(natomwfc,sh(2)/npol)
            do ipol = 1, npol
               evc(:,i + iupdwn(ipol)-1) = wfcatom(:,ipol,i)
            enddo
         enddo
   
         deallocate (tab_at)
         
      end subroutine atomic_wfc_cp
   
   end module atomic_wfc_init