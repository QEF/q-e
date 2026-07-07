!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------
SUBROUTINE PAW_dsymmetrize(dbecsum)
   !---------------------------------------------------------------------------
   !! This routine, similar to PAW_symmetrize, symmetrizes the change of
   !! dbecsum in a linear-response calculation. The symmetry representation
   !! of the perturbations should be already defined in the variables
   !! lr_npert, upert, upert_mq in module lr_symm_base. See for example
   !! the routine PHonon/PH/ph_set_upert.f90.
   !---------------------------------------------------------------------------
   !
   USE kinds,             ONLY : DP
   USE mp,                ONLY : mp_sum
   USE mp_images,         ONLY : nproc_image, me_image, intra_image_comm
   USE noncollin_module,  ONLY : nspin_mag, nspin_lsda, domag
   USE lsda_mod,          ONLY : nspin
   USE uspp_param,        ONLY : nhm
   USE ions_base,         ONLY : nat, ityp
   USE cell_base,         ONLY : at, bg
   USE symm_base,         ONLY : irt, d1, d2, d3, t_rev, sname, s, invs, inverse_s
   USE constants,         ONLY : tpi
   USE uspp,              ONLY : nhtolm, nhtol, ijtoh
   USE uspp_param,        ONLY : nh, upf
   USE io_global,         ONLY : stdout, ionode
   USE qpoint,            ONLY : xq
   USE lr_symm_base,      ONLY : nsymq, lr_npert, upert, rtau, minus_q
   !
   IMPLICIT NONE
   !
   COMPLEX(DP), INTENT(INOUT) :: dbecsum(nhm*(nhm+1)/2, nat, nspin_mag, lr_npert)
   !! cross band occupations
   !
   ! ... local variables
   !
   COMPLEX(DP) :: becsym(nhm*(nhm+1)/2, nat, nspin_mag, lr_npert)
   !! symmetrized becsum
   REAL(DP) :: pref, usym
   !
   INTEGER :: ia, mykey,ia_s,ia_e   ! atoms counters and indexes
   INTEGER :: is, nt       ! counters on spin, atom-type
   INTEGER :: ma           ! atom symmetric to na
   INTEGER :: ih,jh, ijh   ! counters for augmentation channels
   INTEGER :: lm_i, lm_j, &! angular momentums of non-symmetrized becsum
      l_i, l_j, m_i, m_j
   INTEGER :: m_o, m_u     ! counters for sums on m
   INTEGER :: oh, uh, ouh  ! auxiliary indexes corresponding to m_o and m_u
   INTEGER :: isym
   !! counter for symmetry operation
   INTEGER :: ipert, jpert
   !! counters for perturbations
   INTEGER :: ipol, jpol, kpol
   COMPLEX(DP) :: fase(48,nat), mb(3)
   REAL(DP) :: arg, segno
   !
   ! The following mess is necessary because the symmetrization operation
   ! in DFT+Hubbard code is simpler than in PAW, so the required quantities are
   ! represented in a simple but not general way.
   ! I will fix this when everything works.
   REAL(DP), TARGET :: d0(1,1,48)
   TYPE symmetrization_tensor
      REAL(DP),POINTER :: d(:,:,:)
   END TYPE symmetrization_tensor
   TYPE(symmetrization_tensor) :: D(0:3)
   !
   ! Symmetrize with time-reversal symmetry if present
   !
   IF (minus_q) CALL PAW_dmqsymmetrize(dbecsum)
   !
   IF( nsymq==1 ) RETURN
   d0(1,1,:) = 1._dp
   D(0)%d => d0 ! d0(1,1,48)
   D(1)%d => d1 ! d1(3,3,48)
   D(2)%d => d2 ! d2(5,5,48)
   D(3)%d => d3 ! d3(7,7,48)
   !
! => lm = l**2 + m
! => ih = lm + (l+proj)**2  <-- if the projector index starts from zero!
!       = lm + proj**2 + 2*l*proj
!       = m + l**2 + proj**2 + 2*l*proj
!        ^^^
! Known ih and m_i I can compute the index oh of a different m = m_o but
! the same augmentation channel (l_i = l_o, proj_i = proj_o):
!  oh = ih - m_i + m_o
! this expression should be general inside pwscf.
!
!#define __DEBUG_PAW_SYM
   !
   CALL start_clock( 'PAW_dsymm' )
   !
   becsym(:, :, :, :) = (0.0_DP,0.0_DP)
   usym = 1._dp / DBLE(nsymq)
   !
   DO ia = 1, nat
      DO isym = 1, nsymq
         arg = 0.0_DP
         DO ipol = 1, 3
            arg = arg + xq(ipol) * rtau(ipol, isym, ia)
         ENDDO
         arg = arg * tpi
         fase(isym, ia) = CMPLX(COS(arg), SIN(arg), KIND=DP)
      ENDDO
   ENDDO
   !
   ! Parallel: divide among processors for the same image
   CALL block_distribute(nat, me_image, nproc_image, ia_s, ia_e, mykey)
   !
   DO is = 1, nspin_lsda
      !
      atoms: DO ia = ia_s, ia_e
         nt = ityp(ia)
         ! No need to symmetrize non-PAW atoms
         IF ( .NOT. upf(nt)%tpawp ) CYCLE
         !
         DO ih = 1, nh(nt)
            DO jh = ih, nh(nt) ! note: jh >= ih
               !ijh = nh(nt)*(ih-1) - ih*(ih-1)/2 + jh
               ijh = ijtoh(ih,jh,nt)
               !
               lm_i = nhtolm(ih,nt)
               lm_j = nhtolm(jh,nt)
               !
               l_i  = nhtol(ih,nt)
               l_j  = nhtol(jh,nt)
               !
               m_i  = lm_i - l_i**2
               m_j  = lm_j - l_j**2
               !
               DO isym = 1, nsymq
                  ma = irt(isym,ia)
                  DO m_o = 1, 2*l_i +1
                     DO m_u = 1, 2*l_j +1
                        oh = ih - m_i + m_o
                        uh = jh - m_j + m_u
                        ouh = ijtoh(oh,uh,nt)
                        ! In becsum off-diagonal terms are multiplied by 2, I have
                        ! to neutralize this factor and restore it later
                        IF ( oh == uh ) THEN
                           pref = 2._dp * usym
                        ELSE
                           pref = usym
                        ENDIF
                        !
                        DO ipol = 1, lr_npert
                           DO jpol = 1, lr_npert
                              becsym(ijh, ia, is, ipol) = becsym(ijh, ia, is,ipol) &
                                 + D(l_i)%d(m_o,m_i, isym) * D(l_j)%d(m_u,m_j, isym) &
                                 * pref * dbecsum(ouh, ma, is, jpol) * &
                                 upert(jpol, ipol, isym) * fase(isym,ia)
                           ENDDO
                        ENDDO
                     ENDDO ! m_o
                  ENDDO ! m_u
               ENDDO ! isym
               !
               ! Put the prefactor back in:
               IF ( ih == jh ) becsym(ijh,ia,is,:) = .5_dp * becsym(ijh,ia,is,:)
            ENDDO ! ih
         ENDDO ! jh
         !
      ENDDO atoms ! nat
      !
   ENDDO ! nspin
   !
   IF (nspin==4 .AND. domag) THEN
      !
      CALL inverse_s( )
      !
      becsym(:, :, 2:4, :) = 0._dp
      DO ia = 1, nat
         nt = ityp(ia)
         ! No need to symmetrize non-PAW atoms
         IF ( .NOT. upf(nt)%tpawp ) CYCLE
         !
         !  Bring the magnetization in the basis of the crystal
         !
         DO ijh = 1, (nh(nt)*(nh(nt)+1))/2
            DO ipert = 1, lr_npert
               DO jpol = 1, 3
                  mb(jpol) = dbecsum(ijh,ia,jpol+1,ipert)
               ENDDO
               DO jpol = 1, 3
                  dbecsum(ijh,ia,jpol+1,ipert) &
                     = bg(1,jpol)*mb(1) + bg(2,jpol)*mb(2) + bg(3,jpol)*mb(3)
               ENDDO
            ENDDO
         ENDDO
         !
      ENDDO
      !
      DO ia = ia_s, ia_e
         nt = ityp(ia)
         ! No need to symmetrize non-PAW atoms
         IF ( .NOT. upf(nt)%tpawp ) CYCLE
         !
         DO ih = 1, nh(nt)
            DO jh = ih, nh(nt) ! note: jh >= ih
               !ijh = nh(nt)*(ih-1) - ih*(ih-1)/2 + jh
               ijh = ijtoh(ih,jh,nt)
               !
               lm_i  = nhtolm(ih,nt)
               lm_j  = nhtolm(jh,nt)
               !
               l_i   = nhtol(ih,nt)
               l_j   = nhtol(jh,nt)
               !
               m_i   = lm_i - l_i**2
               m_j   = lm_j - l_j**2
               !
               DO isym = 1,nsymq
                  ma = irt(isym,ia)
                  DO m_o = 1, 2*l_i +1
                     DO m_u = 1, 2*l_j +1
                        oh = ih - m_i + m_o
                        uh = jh - m_j + m_u
                        ouh = ijtoh(oh,uh,nt)
                        ! In becsum off-diagonal terms are multiplied by 2, I have
                        ! to neutralize this factor and restore it later
                        IF ( oh == uh ) THEN
                           pref = 2._dp * usym
                        ELSE
                           pref = usym
                        ENDIF
                        !
                        segno = 1.0_DP
                        IF (sname(isym)(1:3)=='inv') segno = -segno
                        IF (t_rev(isym)==1)  segno = -segno
                        !
                        DO ipert = 1, lr_npert
                           DO jpert = 1, lr_npert
                              DO is = 1, 3
                                 DO kpol = 1, 3
                                    becsym(ijh,ia,is+1,ipert)=becsym(ijh,ia,is+1,ipert) &
                                       + D(l_i)%d(m_o,m_i,isym)*D(l_j)%d(m_u,m_j,isym)*  &
                                       pref*dbecsum(ouh,ma,kpol+1,jpert)*               &
                                       upert(jpert, ipert, isym)*fase(isym,ia)*            &
                                       segno*s(kpol,is,invs(isym))
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                        !
                     ENDDO ! m_o
                  ENDDO ! m_u
               ENDDO ! isym
               ! Put the prefactor back in:
               IF ( ih == jh ) becsym(ijh,ia,is,:) = 0.5_dp * becsym(ijh,ia,is,:)
               !
            ENDDO ! ih
         ENDDO ! jh
         !
      ENDDO ! nat
      !
      !
      DO ia = ia_s, ia_e
         nt = ityp(ia)
         ! No need to symmetrize non-PAW atoms
         IF ( .NOT. upf(nt)%tpawp ) CYCLE
         !
         !  Bring the magnetization in cartesian basis
         !
         DO ijh = 1, (nh(nt)*(nh(nt)+1))/2
            DO ipert = 1, lr_npert
               DO jpol = 1, 3
                  mb(jpol) = becsym(ijh,ia,jpol+1,ipert)
               ENDDO
               DO jpol = 1, 3
                  becsym(ijh,ia,jpol+1,ipert) = at(jpol,1)*mb(1) + at(jpol,2)*mb(2) + &
                     at(jpol,3)*mb(3)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !
   ENDIF ! (nspin=4 and domag)
   !
   !
#if defined(__MPI)
   IF( mykey /= 0 ) becsym = 0.0_dp
   CALL mp_sum( becsym, intra_image_comm )
#endif
   !
#if defined(__DEBUG_PAW_SYM)
   WRITE(stdout,*) "------------"
   IF (ionode) THEN
      ia = 1
      nt = ityp(ia)
      DO is = 1, nspin_lsda
         WRITE(*,*) is
         DO ih = 1, nh(nt)
            DO jh = 1, nh(nt)
               ijh = ijtoh(ih,jh,nt)
               DO ipert = 1, lr_npert
                  WRITE(stdout,"(1f10.3)", ADVANCE='no') becsym(ijh,ia,is,ipert)
               ENDDO
            ENDDO
            WRITE(stdout,*)
         ENDDO
         WRITE(stdout,*)
      ENDDO
   endif
   WRITE(stdout,*) "------------"
#endif
   !
   ! Apply symmetrization:
   dbecsum(:,:,:,:) = becsym(:,:,:,:)
   !
   CALL stop_clock( 'PAW_dsymm' )
   !
!------------------------------------------------------------------------------
END SUBROUTINE PAW_dsymmetrize
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
SUBROUTINE PAW_dmqsymmetrize(dbecsum)
   !---------------------------------------------------------------------------
   !! This routine, similar to PAW_symmetrize, symmetrizes the change of
   !! dbecsum in a linear-response calculation for time-reversal symmetry
   !! (mapping q to -q). The symmetry index irotmq and the representation of the
   !! perturbations upert_mq should be already defined in module lr_symm_base.
   !! See for example the routine PHonon/PH/ph_set_upert.f90.
   !---------------------------------------------------------------------------
   !
   USE kinds,             ONLY : DP
   USE mp,                ONLY : mp_sum
   USE mp_images,         ONLY : nproc_image, me_image, intra_image_comm
   USE noncollin_module,  ONLY : nspin_lsda, nspin_mag
   USE uspp_param,        ONLY : nhm
   USE ions_base,         ONLY : nat, ityp
   USE constants,         ONLY : tpi
   USE symm_base,         ONLY : irt, d1, d2, d3
   USE uspp,              ONLY : nhtolm, nhtol, ijtoh
   USE uspp_param,        ONLY : nh, upf
   USE io_global,         ONLY : stdout, ionode
   USE qpoint,            ONLY : xq
   USE lr_symm_base,      ONLY : irotmq, lr_npert, upert_mq, rtau
   !
   !
   COMPLEX(DP), INTENT(INOUT) :: dbecsum(nhm*(nhm+1)/2, nat, nspin_mag, lr_npert)
   !! cross band occupations
   !
   ! ... local variables
   !
   COMPLEX(DP) :: becsym(nhm*(nhm+1)/2, nat, nspin_mag, lr_npert)! symmetrized becsum
   REAL(DP) :: pref
   !
   INTEGER :: ia, mykey,ia_s,ia_e   ! atoms counters and indexes
   INTEGER :: is, nt       ! counters on spin, atom-type
   INTEGER :: ma           ! atom symmetric to na
   INTEGER :: ih,jh, ijh   ! counters for augmentation channels
   INTEGER :: lm_i, lm_j, &! angular momentums of non-symmetrized becsum
      l_i, l_j, m_i, m_j
   INTEGER :: m_o, m_u     ! counters for sums on m
   INTEGER :: oh, uh, ouh  ! auxiliary indexes corresponding to m_o and m_u
   INTEGER :: ipol
   INTEGER :: ipert, jpert
   !! counters for perturbations
   REAL(DP) :: arg
   COMPLEX(DP) :: fase(nat)
   !
   ! The following mess is necessary because the symmetrization operation
   ! in DFT+Hubbard code is simpler than in PAW, so the required quantities are
   ! represented in a simple but not general way.
   ! I will fix this when everything works.
   REAL(DP), TARGET :: d0(1,1,48)
   TYPE symmetrization_tensor
      REAL(DP),POINTER :: d(:,:,:)
   END TYPE symmetrization_tensor
   TYPE(symmetrization_tensor) :: D(0:3)
   !
   IF (nspin_mag==4) CALL errore( 'PAW_dmqsymmetrize', &
   & 'This should not happen (magnetic calculations do not have time-reversal symmetry)', 1 )
   !
   CALL start_clock( 'PAW_dmqsymm' )
   !
   d0(1,1,:) = 1._dp
   D(0)%d => d0 ! d0(1,1,48)
   D(1)%d => d1 ! d1(3,3,48)
   D(2)%d => d2 ! d2(5,5,48)
   D(3)%d => d3 ! d3(7,7,48)
   !
! => lm = l**2 + m
! => ih = lm + (l+proj)**2  <-- if the projector index starts from zero!
!       = lm + proj**2 + 2*l*proj
!       = m + l**2 + proj**2 + 2*l*proj
!        ^^^
! Known ih and m_i I can compute the index oh of a different m = m_o but
! the same augmentation channel (l_i = l_o, proj_i = proj_o):
!  oh = ih - m_i + m_o
! this expression should be general inside pwscf.
!
!#define __DEBUG_PAW_SYM
   !
   becsym(:,:,:,:) = (0.0_DP,0.0_DP)
   DO ia = 1, nat
      arg = 0.0_DP
      DO ipol = 1, 3
         arg = arg + xq (ipol) *  rtau(ipol, irotmq, ia)
      ENDDO
      arg = arg * tpi
      fase(ia) = CMPLX( COS(arg), SIN(arg), KIND=DP)
   ENDDO
   !
   ! Parallel: divide among processors for the same image
   CALL block_distribute( nat, me_image, nproc_image, ia_s, ia_e, mykey )
   !
   DO is = 1, nspin_lsda
      !
      atoms: DO ia = ia_s, ia_e
         nt = ityp(ia)
         ! No need to symmetrize non-PAW atoms
         IF ( .NOT. upf(nt)%tpawp ) CYCLE
         !
         DO ih = 1, nh(nt)
            DO jh = ih, nh(nt) ! note: jh >= ih
               !ijh = nh(nt)*(ih-1) - ih*(ih-1)/2 + jh
               ijh = ijtoh(ih,jh,nt)
               !
               lm_i  = nhtolm(ih,nt)
               lm_j  = nhtolm(jh,nt)
               !
               l_i   = nhtol(ih,nt)
               l_j   = nhtol(jh,nt)
               !
               m_i   = lm_i - l_i**2
               m_j   = lm_j - l_j**2
               !
               ma = irt(irotmq,ia)
               !
               DO m_o = 1, 2*l_i +1
                  DO m_u = 1, 2*l_j +1
                     oh = ih - m_i + m_o
                     uh = jh - m_j + m_u
                     ouh = ijtoh(oh,uh,nt)
                     ! In becsum off-diagonal terms are multiplied by 2, I have
                     ! to neutralize this factor and restore it later
                     IF ( oh == uh ) THEN
                        pref = 2._dp
                     ELSE
                        pref = 1._DP
                     ENDIF
                     !
                     DO ipert = 1, lr_npert
                        DO jpert = 1, lr_npert
                           becsym(ijh, ia, is, ipert) = becsym(ijh, ia, is, ipert)     &
                              + D(l_i)%d(m_o,m_i, irotmq) * D(l_j)%d(m_u,m_j, irotmq) &
                              * pref * dbecsum(ouh, ma, is, jpert) *               &
                              upert_mq(jpert, ipert) * fase(ia)
                        ENDDO
                     ENDDO
                     !
                  ENDDO ! m_o
               ENDDO ! m_u
               !
               ! Put the prefactor back in:
               IF ( ih == jh ) becsym(ijh,ia,is,:) = .5_dp * becsym(ijh,ia,is,:)
               becsym(ijh, ia, is, :) = ( CONJG(becsym(ijh, ia, is, :)) &
                                        + dbecsum(ijh, ia, is, :) ) * 0.5_DP
            ENDDO ! ih
         ENDDO ! jh
         !
      ENDDO atoms ! nat
      !
   ENDDO ! nspin
   !
   !
#if defined(__MPI)
   IF( mykey /= 0 ) becsym = 0.0_dp
   CALL mp_sum( becsym, intra_image_comm )
#endif
   !
#if defined(__DEBUG_PAW_SYM)
   WRITE(stdout,*) "------------"
   IF (ionode) THEN
      ia = 1
      nt = ityp(ia)
      DO is = 1, nspin_mag
         WRITE(*,*) is
         DO ih = 1, nh(nt)
            DO jh = 1, nh(nt)
               ijh = ijtoh(ih,jh,nt)
               DO ipert = 1, lr_npert
                  WRITE(stdout,"(1f10.3)", ADVANCE='no') becsym(ijh,ia,is,ipert)
               ENDDO
            ENDDO
            WRITE(stdout,*)
         ENDDO
         WRITE(stdout,*)
      ENDDO
   ENDIF
   WRITE(stdout,*) "------------"
#endif
   !
   ! Apply symmetrization:
   dbecsum(:,:,:,:) = becsym(:,:,:,:)
   !
   CALL stop_clock( 'PAW_dmqsymm' )
   !
!------------------------------------------------------------------------------
END SUBROUTINE PAW_dmqsymmetrize
!------------------------------------------------------------------------------
