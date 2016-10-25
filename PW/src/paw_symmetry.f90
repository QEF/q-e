!
! Copyright (C) 2007-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE paw_symmetry
    !
    USE kinds,          ONLY : DP
    USE mp_images,      ONLY : nproc_image, me_image, intra_image_comm
    USE mp,             ONLY : mp_sum
    !
    IMPLICIT NONE

    ! entry points:
    PUBLIC :: PAW_symmetrize ! symmetrize becsums
    PUBLIC :: PAW_symmetrize_ddd ! symmetrize the D coefficients
    PUBLIC :: PAW_desymmetrize! symmetrize dbecsums for electric field
    PUBLIC :: PAW_dusymmetrize! symmetrize dbecsums for phonon modes
    PUBLIC :: PAW_dumqsymmetrize! symmetrize dbecsums for phonon modes
                             ! with respect to minus_q
    !
    PRIVATE

 CONTAINS

SUBROUTINE PAW_symmetrize(becsum)
    USE lsda_mod,          ONLY : nspin
    USE cell_base,         ONLY : at, bg
    USE noncollin_module,  ONLY : nspin_mag, nspin_lsda
    USE spin_orb,          ONLY : domag
    USE uspp_param,        ONLY : nhm
    USE ions_base,         ONLY : nat, ityp
    USE symm_base,         ONLY : nsym, irt, d1, d2, d3, t_rev, sname, s, &
                                  invs, inverse_s
    USE uspp,              ONLY : nhtolm,nhtol,ijtoh
    USE uspp_param,        ONLY : nh, upf
    USE io_global,         ONLY : stdout, ionode

    REAL(DP), INTENT(INOUT) :: becsum(nhm*(nhm+1)/2,nat,nspin)! cross band occupations

    REAL(DP)                :: becsym(nhm*(nhm+1)/2,nat,nspin)! symmetrized becsum
    REAL(DP) :: pref, usym, segno
    REAL(DP) :: mb(3)

    INTEGER :: ia,mykey,ia_s,ia_e 
                            ! atoms counters and indexes
    INTEGER :: is, nt       ! counters on spin, atom-type
    INTEGER :: ma           ! atom symmetric to na
    INTEGER :: ih,jh, ijh   ! counters for augmentation channels
    INTEGER :: lm_i, lm_j, &! angular momentums of non-symmetrized becsum
               l_i, l_j, m_i, m_j
    INTEGER :: m_o, m_u     ! counters for sums on m
    INTEGER :: oh, uh, ouh  ! auxiliary indexes corresponding to m_o and m_u
    INTEGER :: isym         ! counter for symmetry operation
    INTEGER :: ipol, kpol
    INTEGER :: table(48, 48)

    ! The following mess is necessary because the symmetrization operation
    ! in LDA+U code is simpler than in PAW, so the required quantities are
    ! represented in a simple but not general way.
    ! I will fix this when everything works.
    REAL(DP), TARGET :: d0(1,1,48)
    TYPE symmetrization_tensor
        REAL(DP),POINTER :: d(:,:,:)
    END TYPE symmetrization_tensor
    TYPE(symmetrization_tensor) :: D(0:3)


    IF( nsym==1 ) RETURN
    d0(1,1,:) = 1._dp
    D(0)%d => d0 ! d0(1,1,48)
    D(1)%d => d1 ! d1(3,3,48)
    D(2)%d => d2 ! d2(5,5,48)
    D(3)%d => d3 ! d3(7,7,48)

! => lm = l**2 + m
! => ih = lm + (l+proj)**2  <-- if the projector index starts from zero!
!       = lm + proj**2 + 2*l*proj
!       = m + l**2 + proj**2 + 2*l*proj
!        ^^^
! Known ih and m_i I can compute the index oh of a different m = m_o but
! the same augmentation channel (l_i = l_o, proj_i = proj_o):
!  oh = ih - m_i + m_o
! this expression should be general inside pwscf.

!#define __DEBUG_PAW_SYM


    CALL start_clock('PAW_symme')

    becsym(:,:,:) = 0._dp
    usym = 1._dp / DBLE(nsym)

    ! Parallel: divide among processors for the same image
    CALL block_distribute( nat, me_image, nproc_image, ia_s, ia_e, mykey )

    DO is = 1, nspin_lsda
    !
    atoms: DO ia = ia_s, ia_e
        nt = ityp(ia)
        ! No need to symmetrize non-PAW atoms
        IF ( .not. upf(nt)%tpawp ) CYCLE
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
            DO isym = 1,nsym
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
                    becsym(ijh, ia, is) = becsym(ijh, ia, is) &
                        + D(l_i)%d(m_o,m_i, isym) * D(l_j)%d(m_u,m_j, isym) &
                          * pref * becsum(ouh, ma, is)
                ENDDO ! m_o
                ENDDO ! m_u
            ENDDO ! isym
            !
            ! Put the prefactor back in:
            IF ( ih == jh ) becsym(ijh,ia,is) = .5_dp * becsym(ijh,ia,is)
        ENDDO ! ih
        ENDDO ! jh
    ENDDO atoms ! nat
    ENDDO ! nspin

    IF (nspin==4.and.domag) THEN
       !
       call inverse_s( )

       becsym(:,:,2:4) = 0._dp
       DO ia = 1, nat
          nt = ityp(ia)
          ! No need to symmetrize non-PAW atoms
          IF ( .not. upf(nt)%tpawp ) CYCLE
          !
          !  Bring the magnetization in the basis of the crystal
          !        
          DO ijh=1,(nh(nt)*(nh(nt)+1))/2
             DO ipol=1,3
                mb(ipol)=becsum(ijh,ia,ipol+1)
             ENDDO
             DO ipol=1,3
                becsum(ijh,ia,ipol+1)=bg(1,ipol)*mb(1)+bg(2,ipol)*mb(2) + &
                                      bg(3,ipol)*mb(3) 
             END DO
          END DO
       END DO
       atoms_1: DO ia = ia_s, ia_e
          nt = ityp(ia)
          ! No need to symmetrize non-PAW atoms
          IF ( .not. upf(nt)%tpawp ) CYCLE
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
               DO isym = 1,nsym
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
                      segno=1.0_DP
                      IF (sname(isym)(1:3)=='inv') segno=-segno
                      IF (t_rev(isym)==1)  segno=-segno
           
                      DO is=1,3
                      DO kpol=1,3
                         becsym(ijh, ia, is+1) = becsym(ijh, ia, is+1) &
                       + D(l_i)%d(m_o,m_i, isym) * D(l_j)%d(m_u,m_j, isym) &
                          * pref * becsum(ouh, ma, kpol+1)*&
                                   s(kpol,is,invs(isym))* &
                            segno
                      ENDDO
                      ENDDO
                  ENDDO ! m_o
                  ENDDO ! m_u
               ENDDO ! isym
            !
            ! Put the prefactor back in:
            IF ( ih == jh ) becsym(ijh,ia,2:4) = .5_dp * becsym(ijh,ia,2:4)
        ENDDO ! ih
        ENDDO ! jh
      ENDDO atoms_1 ! nat

      DO ia = ia_s, ia_e
         nt = ityp(ia)
         ! No need to symmetrize non-PAW atoms
         IF ( .not. upf(nt)%tpawp ) CYCLE
         !
         !  Bring the magnetization in cartesian basis
         !        
         DO ijh=1,(nh(nt)*(nh(nt)+1))/2
            DO ipol=1,3
               mb(ipol)=becsym(ijh,ia,ipol+1)
            ENDDO
            DO ipol=1,3
               becsym(ijh,ia,ipol+1)=at(ipol,1)*mb(1)+at(ipol,2)*mb(2) + &
                                     at(ipol,3)*mb(3)
            END DO
         END DO
      END DO
   END IF

#if defined(__MPI)
    IF( mykey /= 0 ) becsym = 0.0_dp
    CALL mp_sum(becsym, intra_image_comm)
#endif

#if defined(__DEBUG_PAW_SYM)
   write(stdout,*) "------------"
    if(ionode) then
        ia = 1
        nt = ityp(ia)
        DO is = 1, nspin
            write(*,*) is
        DO ih = 1, nh(nt)
        DO jh = 1, nh(nt)
            ijh = ijtoh(ih,jh,nt)
            write(stdout,"(1f10.3)", advance='no') becsym(ijh,ia,is)
        ENDDO
            write(stdout,*)
        ENDDO
            write(stdout,*)
        ENDDO
    endif
   write(stdout,*) "------------"
#endif

    ! Apply symmetrization:
    becsum(:,:,:) = becsym(:,:,:)

    CALL stop_clock('PAW_symme')

END SUBROUTINE PAW_symmetrize

SUBROUTINE PAW_symmetrize_ddd(ddd)
    USE lsda_mod,          ONLY : nspin
    USE cell_base,         ONLY : at, bg
    USE noncollin_module,  ONLY : nspin_mag, nspin_lsda
    USE spin_orb,          ONLY : domag
    USE uspp_param,        ONLY : nhm
    USE ions_base,         ONLY : nat, ityp
    USE symm_base,         ONLY : nsym, irt, d1, d2, d3, t_rev, sname, s, &
                                  invs, inverse_s
    USE uspp,              ONLY : nhtolm,nhtol,ijtoh
    USE uspp_param,        ONLY : nh, upf
    USE io_global,         ONLY : stdout, ionode

    REAL(DP), INTENT(INOUT) :: ddd(nhm*(nhm+1)/2,nat,nspin)! cross band occupations

    REAL(DP)                :: dddsym(nhm*(nhm+1)/2,nat,nspin)! symmetrized becsum
    REAL(DP) :: usym, segno
    REAL(DP) :: mb(3)

    INTEGER :: ia,mykey,ia_s,ia_e 
                            ! atoms counters and indexes
    INTEGER :: is, nt       ! counters on spin, atom-type
    INTEGER :: ma           ! atom symmetric to na
    INTEGER :: ih,jh, ijh   ! counters for augmentation channels
    INTEGER :: lm_i, lm_j, &! angular momentums of non-symmetrized becsum
               l_i, l_j, m_i, m_j
    INTEGER :: m_o, m_u     ! counters for sums on m
    INTEGER :: oh, uh, ouh  ! auxiliary indexes corresponding to m_o and m_u
    INTEGER :: isym         ! counter for symmetry operation
    INTEGER :: ipol, kpol
    INTEGER :: table(48, 48)

    ! The following mess is necessary because the symmetrization operation
    ! in LDA+U code is simpler than in PAW, so the required quantities are
    ! represented in a simple but not general way.
    ! I will fix this when everything works.
    REAL(DP), TARGET :: d0(1,1,48)
    TYPE symmetrization_tensor
        REAL(DP),POINTER :: d(:,:,:)
    END TYPE symmetrization_tensor
    TYPE(symmetrization_tensor) :: D(0:3)

    IF( nsym==1 ) RETURN
    d0(1,1,:) = 1._dp
    D(0)%d => d0 ! d0(1,1,48)
    D(1)%d => d1 ! d1(3,3,48)
    D(2)%d => d2 ! d2(5,5,48)
    D(3)%d => d3 ! d3(7,7,48)

! => lm = l**2 + m
! => ih = lm + (l+proj)**2  <-- if the projector index starts from zero!
!       = lm + proj**2 + 2*l*proj
!       = m + l**2 + proj**2 + 2*l*proj
!        ^^^
! Known ih and m_i I can compute the index oh of a different m = m_o but
! the same augmentation channel (l_i = l_o, proj_i = proj_o):
!  oh = ih - m_i + m_o
! this expression should be general inside pwscf.

!#define __DEBUG_PAW_SYM


    CALL start_clock('PAW_symme')
    
    dddsym(:,:,:) = 0._dp
    usym = 1._dp / DBLE(nsym)

    ! Parallel: divide among processors for the same image
    CALL block_distribute( nat, me_image, nproc_image, ia_s, ia_e, mykey )

    DO is = 1, nspin_lsda
    !
    atoms: DO ia = ia_s, ia_e
        nt = ityp(ia)
        ! No need to symmetrize non-PAW atoms
        IF ( .not. upf(nt)%tpawp ) CYCLE
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
            DO isym = 1,nsym
                ma = irt(isym,ia)
                DO m_o = 1, 2*l_i +1
                DO m_u = 1, 2*l_j +1
                    oh = ih - m_i + m_o
                    uh = jh - m_j + m_u
                    ouh = ijtoh(oh,uh,nt)
                    !
                    dddsym(ijh, ia, is) = dddsym(ijh, ia, is) &
                        + D(l_i)%d(m_o,m_i, isym) * D(l_j)%d(m_u,m_j, isym) &
                          * usym * ddd(ouh, ma, is)
                ENDDO ! m_o
                ENDDO ! m_u
            ENDDO ! isym
            !
        ENDDO ! ih
        ENDDO ! jh
    ENDDO atoms ! nat
    ENDDO ! nspin

    IF (nspin==4.and.domag) THEN
       !
       call inverse_s( )

       dddsym(:,:,2:4) = 0._dp
       DO ia = 1, nat
          nt = ityp(ia)
          ! No need to symmetrize non-PAW atoms
          IF ( .not. upf(nt)%tpawp ) CYCLE
          !
          !  Bring the magnetization in the basis of the crystal
          !        
          DO ijh=1,(nh(nt)*(nh(nt)+1))/2
             DO ipol=1,3
                mb(ipol)=ddd(ijh,ia,ipol+1)
             ENDDO
             DO ipol=1,3
                ddd(ijh,ia,ipol+1)=bg(1,ipol)*mb(1)+bg(2,ipol)*mb(2) + &
                                   bg(3,ipol)*mb(3) 
             END DO
          END DO
       END DO
       atoms_1: DO ia = ia_s, ia_e
          nt = ityp(ia)
          ! No need to symmetrize non-PAW atoms
          IF ( .not. upf(nt)%tpawp ) CYCLE
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
               DO isym = 1,nsym
                  ma = irt(isym,ia)
                  segno=1.0_DP
                  IF (sname(invs(isym))(1:3)=='inv') segno=-segno
                  IF (t_rev(invs(isym))==1)  segno=-segno
                  DO m_o = 1, 2*l_i +1
                  DO m_u = 1, 2*l_j +1
                      oh = ih - m_i + m_o
                      uh = jh - m_j + m_u
                      ouh = ijtoh(oh,uh,nt)
                      !
                      DO is=1,3
                      DO kpol=1,3
                         dddsym(ijh, ia, is+1) = dddsym(ijh, ia, is+1) &
                       + D(l_i)%d(m_o,m_i, isym) * D(l_j)%d(m_u,m_j, isym) &
                          * usym * ddd(ouh, ma, kpol+1)*&
                                   s(kpol,is,invs(isym))*segno
                      ENDDO
                      ENDDO
                  ENDDO ! m_o
                  ENDDO ! m_u
               ENDDO ! isym
            !
        ENDDO ! ih
        ENDDO ! jh
      ENDDO atoms_1 ! nat

      DO ia = ia_s, ia_e
         nt = ityp(ia)
         ! No need to symmetrize non-PAW atoms
         IF ( .not. upf(nt)%tpawp ) CYCLE
         !
         !  Bring the magnetization in cartesian basis
         !        
         DO ijh=1,(nh(nt)*(nh(nt)+1))/2
            DO ipol=1,3
               mb(ipol)=dddsym(ijh,ia,ipol+1)
            ENDDO
            DO ipol=1,3
               dddsym(ijh,ia,ipol+1)=at(ipol,1)*mb(1)+at(ipol,2)*mb(2) + &
                                     at(ipol,3)*mb(3)
            END DO
         END DO
      END DO
   END IF

#if defined(__MPI)
    IF( mykey /= 0 ) dddsym = 0.0_dp
    CALL mp_sum(dddsym, intra_image_comm)
#endif

#if defined(__DEBUG_PAW_SYM)
   write(stdout,*) "------------"
    if(ionode) then
        ia = 1
        nt = ityp(ia)
        DO is = 1, nspin
            write(*,*) is
        DO ih = 1, nh(nt)
        DO jh = 1, nh(nt)
            ijh = ijtoh(ih,jh,nt)
            write(stdout,"(1f10.3)", advance='no') dddsym(ijh,ia,is)
        ENDDO
            write(stdout,*)
        ENDDO
            write(stdout,*)
        ENDDO
    endif
   write(stdout,*) "------------"
#endif

    ! Apply symmetrization:
    ddd(:,:,:) = dddsym(:,:,:)

    CALL stop_clock('PAW_symme')

END SUBROUTINE PAW_symmetrize_ddd

SUBROUTINE PAW_desymmetrize(dbecsum)
!
! This routine similar to PAW_symmetrize, symmetrize the change of 
! dbecsum due to an electric field perturbation. 
!
    USE lsda_mod,          ONLY : nspin
    USE uspp_param,        ONLY : nhm
    USE ions_base,         ONLY : nat, ityp
    USE noncollin_module,  ONLY : nspin_lsda, nspin_mag
    USE cell_base,         ONLY : at, bg
    USE spin_orb,          ONLY : domag
    USE symm_base,         ONLY : nsym, irt, d1, d2, d3, s, t_rev, sname, &
                                  invs, inverse_s
    USE uspp,              ONLY : nhtolm,nhtol,ijtoh
    USE uspp_param,        ONLY : nh, upf
    USE io_global,         ONLY : stdout, ionode

    COMPLEX(DP), INTENT(INOUT) :: dbecsum(nhm*(nhm+1)/2,nat,nspin_mag,3)! cross band occupations

    COMPLEX(DP)                :: becsym(nhm*(nhm+1)/2,nat,nspin_mag,3)! symmetrized becsum
    COMPLEX(DP) :: mb(3)
    REAL(DP) :: pref, usym, segno

    INTEGER :: ia, mykey,ia_s,ia_e   ! atoms counters and indexes
    INTEGER :: is, nt       ! counters on spin, atom-type
    INTEGER :: ma           ! atom symmetric to na
    INTEGER :: ih,jh, ijh   ! counters for augmentation channels
    INTEGER :: lm_i, lm_j, &! angular momentums of non-symmetrized becsum
               l_i, l_j, m_i, m_j
    INTEGER :: m_o, m_u     ! counters for sums on m
    INTEGER :: oh, uh, ouh  ! auxiliary indexes corresponding to m_o and m_u
    INTEGER :: isym         ! counter for symmetry operation
    INTEGER :: ipol, jpol, kpol
    INTEGER :: table(48, 48)

    ! The following mess is necessary because the symmetrization operation
    ! in LDA+U code is simpler than in PAW, so the required quantities are
    ! represented in a simple but not general way.
    ! I will fix this when everything works.
    REAL(DP), TARGET :: d0(1,1,48)
    TYPE symmetrization_tensor
        REAL(DP),POINTER :: d(:,:,:)
    END TYPE symmetrization_tensor
    TYPE(symmetrization_tensor) :: D(0:3)

    IF( nsym == 1 ) RETURN
    d0(1,1,:) = 1._dp
    D(0)%d => d0 ! d0(1,1,48)
    D(1)%d => d1 ! d1(3,3,48)
    D(2)%d => d2 ! d2(5,5,48)
    D(3)%d => d3 ! d3(7,7,48)

! => lm = l**2 + m
! => ih = lm + (l+proj)**2  <-- if the projector index starts from zero!
!       = lm + proj**2 + 2*l*proj
!       = m + l**2 + proj**2 + 2*l*proj
!        ^^^
! Known ih and m_i I can compute the index oh of a different m = m_o but
! the same augmentation channel (l_i = l_o, proj_i = proj_o):
!  oh = ih - m_i + m_o
! this expression should be general inside pwscf.

!#define __DEBUG_PAW_SYM


    CALL start_clock('PAW_dsymme')

    becsym(:,:,:,:) = (0.0_DP,0.0_DP)
    usym = 1._dp / DBLE(nsym)

    ! Parallel: divide among processors for the same image
    CALL block_distribute( nat, me_image, nproc_image, ia_s, ia_e, mykey )

    DO is = 1, nspin_lsda
    !
    atoms: DO ia = ia_s, ia_e
        nt = ityp(ia)
        ! No need to symmetrize non-PAW atoms
        IF ( .not. upf(nt)%tpawp ) CYCLE
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
            DO isym = 1,nsym
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
                    DO ipol=1,3
                       DO jpol=1,3
                          becsym(ijh, ia, is, ipol) = becsym(ijh, ia, is,ipol) &
                        + D(l_i)%d(m_o,m_i, isym) * D(l_j)%d(m_u,m_j, isym) &
                          * pref * dbecsum(ouh, ma, is, jpol) * s(ipol,jpol,isym)
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
    ENDDO atoms ! nat
    ENDDO ! nspin

    IF (nspin==4.and.domag) THEN
       !
       !
       call inverse_s ( )

       becsym(:,:,2:4,1:3) = 0._dp
       DO ia = 1, nat
          nt = ityp(ia)
          ! No need to symmetrize non-PAW atoms
          IF ( .not. upf(nt)%tpawp ) CYCLE
          !
          !  Bring the magnetization in the basis of the crystal
          !        
          DO ijh=1,(nh(nt)*(nh(nt)+1))/2
             DO ipol=1,3
                DO jpol=1,3
                   mb(jpol)=dbecsum(ijh,ia,jpol+1,ipol)
                ENDDO
                DO jpol=1,3
                   dbecsum(ijh,ia,jpol+1,ipol)=bg(1,jpol)*mb(1) +  &
                                  bg(2,jpol)*mb(2) + bg(3,jpol)*mb(3)
                ENDDO
             ENDDO
          ENDDO
       ENDDO

       DO ia = ia_s, ia_e
          nt = ityp(ia)
          ! No need to symmetrize non-PAW atoms
          IF ( .not. upf(nt)%tpawp ) CYCLE
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
             DO isym = 1,nsym
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
                   segno=1.0_DP
                   IF (sname(isym)(1:3)=='inv') segno=-segno
                   IF (t_rev(isym)==1)  segno=-segno
                   !
                   DO ipol=1,3
                      DO jpol=1,3
                         DO is=1,3
                          DO kpol=1,3
                           becsym(ijh,ia,is+1,ipol)=becsym(ijh,ia,is+1,ipol) &
                         + D(l_i)%d(m_o,m_i,isym)*D(l_j)%d(m_u,m_j,isym)* &
                           pref*dbecsum(ouh,ma,kpol+1,jpol)*s(ipol,jpol,isym)*&
                           segno*s(kpol,is,invs(isym))
                          END DO
                         END DO   
                      END DO
                   END DO
                END DO ! m_o
                END DO ! m_u
            END DO ! isym
            !
            ! Put the prefactor back in:
            IF ( ih == jh ) becsym(ijh,ia,is,:) = .5_dp * becsym(ijh,ia,is,:)
        ENDDO ! ih
        ENDDO ! jh
     ENDDO  ! nat
!
      DO ia = ia_s, ia_e
         nt = ityp(ia)
         ! No need to symmetrize non-PAW atoms
         IF ( .not. upf(nt)%tpawp ) CYCLE
         !
         !  Bring the magnetization in cartesian basis
         !        
         DO ijh=1,(nh(nt)*(nh(nt)+1))/2
            DO ipol=1,3
               DO jpol=1,3
                  mb(jpol)=becsym(ijh,ia,jpol+1,ipol)
               ENDDO
               DO jpol=1,3
                  becsym(ijh,ia,jpol+1,ipol)=at(jpol,1)*mb(1)+at(jpol,2)*mb(2)+&
                                             at(jpol,3)*mb(3)
               END DO
            END DO
         END DO
      END DO
   ENDIF

#if defined(__MPI)
    IF( mykey /= 0 ) becsym = 0.0_dp
    CALL mp_sum(becsym, intra_image_comm)
#endif

#if defined(__DEBUG_PAW_SYM)
   write(stdout,*) "------------"
    if(ionode) then
        ia = 1
        nt = ityp(ia)
        DO is = 1, nspin_mag
            write(*,*) is
        DO ih = 1, nh(nt)
        DO jh = 1, nh(nt)
            ijh = ijtoh(ih,jh,nt)
            DO ipol=1,3
               write(stdout,"(1f10.3)", advance='no') becsym(ijh,ia,is,ipol)
            ENDDO
        ENDDO
            write(stdout,*)
        ENDDO
            write(stdout,*)
        ENDDO
    endif
   write(stdout,*) "------------"
#endif

    ! Apply symmetrization:
    dbecsum(:,:,:,:) = becsym(:,:,:,:)

    CALL stop_clock('PAW_dsymme')

END SUBROUTINE PAW_desymmetrize

SUBROUTINE PAW_dusymmetrize(dbecsum,npe,irr,npertx,nsymq,rtau,xq,t)
!
! This routine similar to PAW_symmetrize, symmetrize the change of 
! dbecsum due to an electric field perturbation. 
!
    USE noncollin_module,  ONLY : nspin_mag, nspin_lsda
    USE lsda_mod,          ONLY : nspin
    USE uspp_param,        ONLY : nhm
    USE ions_base,         ONLY : nat, ityp
    USE cell_base,         ONLY : at, bg
    USE symm_base,         ONLY : irt, d1, d2, d3, t_rev, sname, s, nsym, &
                                  invs, inverse_s
    USE spin_orb,          ONLY : domag
    USE constants,         ONLY : tpi
    USE uspp,              ONLY : nhtolm,nhtol,ijtoh
    USE uspp_param,        ONLY : nh, upf
    USE io_global,         ONLY : stdout, ionode

    INTEGER, INTENT(IN) :: npe, irr, npertx, nsymq
    REAL(DP), INTENT(IN) :: rtau(3,48,nat), xq(3)
    COMPLEX(DP), INTENT(IN) :: t(npertx, npertx, 48, 3*nat)
    COMPLEX(DP), INTENT(INOUT) :: dbecsum(nhm*(nhm+1)/2,nat,nspin_mag,npe)! cross band occupations

    COMPLEX(DP)                :: becsym(nhm*(nhm+1)/2,nat,nspin_mag,npe)! symmetrized becsum
    REAL(DP) :: pref, usym

    INTEGER :: ia, mykey,ia_s,ia_e   ! atoms counters and indexes
    INTEGER :: is, nt       ! counters on spin, atom-type
    INTEGER :: ma           ! atom symmetric to na
    INTEGER :: ih,jh, ijh   ! counters for augmentation channels
    INTEGER :: lm_i, lm_j, &! angular momentums of non-symmetrized becsum
               l_i, l_j, m_i, m_j
    INTEGER :: m_o, m_u     ! counters for sums on m
    INTEGER :: oh, uh, ouh  ! auxiliary indexes corresponding to m_o and m_u
    INTEGER :: isym, irot   ! counter for symmetry operation
    INTEGER :: ipol, jpol
    COMPLEX(DP) :: fase(48,nat), mb(3)
    REAL(DP) :: arg, ft(3), segno
    INTEGER :: kpol
    INTEGER :: table(48, 48)

    ! The following mess is necessary because the symmetrization operation
    ! in LDA+U code is simpler than in PAW, so the required quantities are
    ! represented in a simple but not general way.
    ! I will fix this when everything works.
    REAL(DP), TARGET :: d0(1,1,48)
    TYPE symmetrization_tensor
        REAL(DP),POINTER :: d(:,:,:)
    END TYPE symmetrization_tensor
    TYPE(symmetrization_tensor) :: D(0:3)

    IF( nsymq==1 ) RETURN
    d0(1,1,:) = 1._dp
    D(0)%d => d0 ! d0(1,1,48)
    D(1)%d => d1 ! d1(3,3,48)
    D(2)%d => d2 ! d2(5,5,48)
    D(3)%d => d3 ! d3(7,7,48)

! => lm = l**2 + m
! => ih = lm + (l+proj)**2  <-- if the projector index starts from zero!
!       = lm + proj**2 + 2*l*proj
!       = m + l**2 + proj**2 + 2*l*proj
!        ^^^
! Known ih and m_i I can compute the index oh of a different m = m_o but
! the same augmentation channel (l_i = l_o, proj_i = proj_o):
!  oh = ih - m_i + m_o
! this expression should be general inside pwscf.

!#define __DEBUG_PAW_SYM

    CALL start_clock('PAW_dusymm')

    becsym(:,:,:,:) = (0.0_DP,0.0_DP)
    usym = 1._dp / DBLE(nsymq)

    do ia=1,nat
       do isym=1,nsymq
          irot = isym
          arg = 0.0_DP
          do ipol = 1, 3
             arg = arg + xq (ipol) *  rtau(ipol,irot,ia)
          enddo
          arg = arg * tpi
          fase(irot,ia) = CMPLX(cos (arg),  sin (arg) ,kind=DP)
       enddo
    enddo

    ! Parallel: divide among processors for the same image
    CALL block_distribute( nat, me_image, nproc_image, ia_s, ia_e, mykey )

    DO is = 1, nspin_lsda
    !
    atoms: DO ia = ia_s, ia_e
        nt = ityp(ia)
        ! No need to symmetrize non-PAW atoms
        IF ( .not. upf(nt)%tpawp ) CYCLE
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
                irot = isym
                ma = irt(irot,ia)
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
                    DO ipol=1,npe
                       DO jpol=1,npe
                          becsym(ijh, ia, is, ipol) = becsym(ijh, ia, is,ipol) &
                        + D(l_i)%d(m_o,m_i, irot) * D(l_j)%d(m_u,m_j, irot) &
                          * pref * dbecsum(ouh, ma, is, jpol) * &
                          t(jpol,ipol,irot,irr) * fase(irot,ia)
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
    ENDDO atoms ! nat
    ENDDO ! nspin

    IF (nspin==4.and.domag) THEN
       !
       call inverse_s ()
       !
       becsym(:,:,2:4,1:npe) = 0._dp
       DO ia = 1, nat
          nt = ityp(ia)
          ! No need to symmetrize non-PAW atoms
          IF ( .not. upf(nt)%tpawp ) CYCLE
          !
          !  Bring the magnetization in the basis of the crystal
          !        
          DO ijh=1,(nh(nt)*(nh(nt)+1))/2
             DO ipol=1,npe
                DO jpol=1,3
                   mb(jpol)=dbecsum(ijh,ia,jpol+1,ipol)
                END DO
                DO jpol=1,3
                   dbecsum(ijh,ia,jpol+1,ipol)=bg(1,jpol)*mb(1) +  &
                                  bg(2,jpol)*mb(2) + bg(3,jpol)*mb(3)
                END DO
             END DO
          END DO
       END DO

       DO ia = ia_s, ia_e
          nt = ityp(ia)
          ! No need to symmetrize non-PAW atoms
          IF ( .not. upf(nt)%tpawp ) CYCLE
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
                irot = isym
                ma = irt(irot,ia)
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
                    segno=1.0_DP
                    IF (sname(isym)(1:3)=='inv') segno=-segno
                    IF (t_rev(isym)==1)  segno=-segno

                    DO ipol=1,npe
                       DO jpol=1,npe
                          DO is=1, 3
                             DO kpol=1,3
                          becsym(ijh,ia,is+1,ipol)=becsym(ijh,ia,is+1,ipol) &
                        + D(l_i)%d(m_o,m_i,irot)*D(l_j)%d(m_u,m_j,irot)*    &
                          pref*dbecsum(ouh,ma,kpol+1,jpol)*                 &
                          t(jpol,ipol,irot,irr)*fase(irot,ia)*              &
                          segno*s(kpol,is,invs(isym))
                             ENDDO
                          ENDDO
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
      ENDDO  ! nat

      DO ia = ia_s, ia_e
         nt = ityp(ia)
         ! No need to symmetrize non-PAW atoms
         IF ( .not. upf(nt)%tpawp ) CYCLE
         !
         !  Bring the magnetization in cartesian basis
         !        
         DO ijh=1,(nh(nt)*(nh(nt)+1))/2
            DO ipol=1,npe
               DO jpol=1,3
                  mb(jpol)=becsym(ijh,ia,jpol+1,ipol)
               ENDDO
               DO jpol=1,3
                  becsym(ijh,ia,jpol+1,ipol)=at(jpol,1)*mb(1)+at(jpol,2)*mb(2)+&
                                             at(jpol,3)*mb(3)
               END DO
            END DO
         END DO
      END DO
   END IF



#if defined(__MPI)
    IF( mykey /= 0 ) becsym = 0.0_dp
    CALL mp_sum(becsym, intra_image_comm)
#endif

#if defined(__DEBUG_PAW_SYM)
   write(stdout,*) "------------"
    if(ionode) then
        ia = 1
        nt = ityp(ia)
        DO is = 1, nspin_lsda
            write(*,*) is
        DO ih = 1, nh(nt)
        DO jh = 1, nh(nt)
            ijh = ijtoh(ih,jh,nt)
            DO ipol=1,npe
               write(stdout,"(1f10.3)", advance='no') becsym(ijh,ia,is,ipol)
            ENDDO
        ENDDO
            write(stdout,*)
        ENDDO
            write(stdout,*)
        ENDDO
    endif
   write(stdout,*) "------------"
#endif

    ! Apply symmetrization:
    dbecsum(:,:,:,:) = becsym(:,:,:,:)

    CALL stop_clock('PAW_dusymm')

END SUBROUTINE PAW_dusymmetrize

SUBROUTINE PAW_dumqsymmetrize(dbecsum,npe,irr,npertx,isymq,rtau,xq,tmq)
!
! This routine similar to PAW_symmetrize, symmetrize the change of 
! dbecsum due to an electric field perturbation. 
!
    USE noncollin_module,  ONLY : nspin_lsda, nspin_mag
    USE lsda_mod,          ONLY : nspin
    USE uspp_param,        ONLY : nhm
    USE ions_base,         ONLY : nat, ityp
    USE constants,         ONLY : tpi
    USE symm_base,         ONLY : nsym, irt, d1, d2, d3
    USE uspp,              ONLY : nhtolm,nhtol,ijtoh
    USE uspp_param,        ONLY : nh, upf
    USE io_global,         ONLY : stdout, ionode

    INTEGER, INTENT(IN) :: npe, irr, npertx
    INTEGER, INTENT(IN) :: isymq         ! counter for symmetry operation
    COMPLEX(DP), INTENT(IN) :: tmq(npertx, npertx, 3*nat)
    COMPLEX(DP), INTENT(INOUT) :: dbecsum(nhm*(nhm+1)/2,nat,nspin_mag,npe)! cross band occupations
    REAL(DP), INTENT(IN) :: rtau(3,48,nat), xq(3)

    COMPLEX(DP)                :: becsym(nhm*(nhm+1)/2,nat,nspin_mag,npe)! symmetrized becsum
    REAL(DP) :: pref

    INTEGER :: ia, mykey,ia_s,ia_e   ! atoms counters and indexes
    INTEGER :: is, nt       ! counters on spin, atom-type
    INTEGER :: ma           ! atom symmetric to na
    INTEGER :: ih,jh, ijh   ! counters for augmentation channels
    INTEGER :: lm_i, lm_j, &! angular momentums of non-symmetrized becsum
               l_i, l_j, m_i, m_j
    INTEGER :: m_o, m_u     ! counters for sums on m
    INTEGER :: oh, uh, ouh  ! auxiliary indexes corresponding to m_o and m_u
    INTEGER :: ipol, jpol
    REAL(DP) :: arg
    COMPLEX(DP) :: fase(nat)

    ! The following mess is necessary because the symmetrization operation
    ! in LDA+U code is simpler than in PAW, so the required quantities are
    ! represented in a simple but not general way.
    ! I will fix this when everything works.
    REAL(DP), TARGET :: d0(1,1,48)
    TYPE symmetrization_tensor
        REAL(DP),POINTER :: d(:,:,:)
    END TYPE symmetrization_tensor
    TYPE(symmetrization_tensor) :: D(0:3)

    IF (nspin_mag==4) call errore('PAW_dumqsymmetrize',&
                     & 'This should not happen',1)

    CALL start_clock('PAW_dumqsym')

    d0(1,1,:) = 1._dp
    D(0)%d => d0 ! d0(1,1,48)
    D(1)%d => d1 ! d1(3,3,48)
    D(2)%d => d2 ! d2(5,5,48)
    D(3)%d => d3 ! d3(7,7,48)

! => lm = l**2 + m
! => ih = lm + (l+proj)**2  <-- if the projector index starts from zero!
!       = lm + proj**2 + 2*l*proj
!       = m + l**2 + proj**2 + 2*l*proj
!        ^^^
! Known ih and m_i I can compute the index oh of a different m = m_o but
! the same augmentation channel (l_i = l_o, proj_i = proj_o):
!  oh = ih - m_i + m_o
! this expression should be general inside pwscf.

!#define __DEBUG_PAW_SYM


    becsym(:,:,:,:) = (0.0_DP,0.0_DP)
    do ia=1,nat
       arg = 0.0_DP
       do ipol = 1, 3
          arg = arg + xq (ipol) *  rtau(ipol,isymq,ia)
       enddo
       arg = arg * tpi
       fase(ia) = CMPLX(cos (arg),  sin (arg) ,kind=DP)
    enddo


    ! Parallel: divide among processors for the same image
    CALL block_distribute( nat, me_image, nproc_image, ia_s, ia_e, mykey )

    DO is = 1, nspin_lsda
    !
    atoms: DO ia = ia_s, ia_e
        nt = ityp(ia)
        ! No need to symmetrize non-PAW atoms
        IF ( .not. upf(nt)%tpawp ) CYCLE
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
            ma = irt(isymq,ia)
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
               DO ipol=1,npe
                  DO jpol=1,npe
                     becsym(ijh, ia, is, ipol) = becsym(ijh, ia, is,ipol) &
                        + D(l_i)%d(m_o,m_i, isymq) * D(l_j)%d(m_u,m_j, isymq) &
                          * pref * dbecsum(ouh, ma, is, jpol) * &
                          tmq(jpol,ipol,irr)*fase(ia)
                  ENDDO
               ENDDO
            ENDDO ! m_o
            ENDDO ! m_u
            !
            ! Put the prefactor back in:
            IF ( ih == jh ) becsym(ijh,ia,is,:) = .5_dp * becsym(ijh,ia,is,:)
            becsym(ijh, ia, is,:)=(CONJG(becsym(ijh, ia, is, :))+ &
                                            dbecsum(ijh, ia, is, :))*0.5_DP
        ENDDO ! ih
        ENDDO ! jh
    ENDDO atoms ! nat
    ENDDO ! nspin

#if defined(__MPI)
    IF( mykey /= 0 ) becsym = 0.0_dp
    CALL mp_sum(becsym, intra_image_comm)
#endif

#if defined(__DEBUG_PAW_SYM)
   write(stdout,*) "------------"
    if(ionode) then
        ia = 1
        nt = ityp(ia)
        DO is = 1, nspin_mag
            write(*,*) is
        DO ih = 1, nh(nt)
        DO jh = 1, nh(nt)
            ijh = ijtoh(ih,jh,nt)
            DO ipol=1,npe
               write(stdout,"(1f10.3)", advance='no') becsym(ijh,ia,is,ipol)
            ENDDO
        ENDDO
            write(stdout,*)
        ENDDO
            write(stdout,*)
        ENDDO
    endif
   write(stdout,*) "------------"
#endif

    ! Apply symmetrization:
    dbecsum(:,:,:,:) = becsym(:,:,:,:)

    CALL stop_clock('PAW_dumqsym')

END SUBROUTINE PAW_dumqsymmetrize

END MODULE paw_symmetry
