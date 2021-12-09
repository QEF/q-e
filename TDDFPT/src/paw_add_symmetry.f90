!
! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!  This modules adds functionalities not available in paw_symmetry
!  but used by thermo_pw.
!
!
MODULE paw_add_symmetry
    !
    USE kinds,          ONLY : DP
    USE mp_images,      ONLY : nproc_image, me_image, intra_image_comm
    USE mp,             ONLY : mp_sum
    !
    IMPLICIT NONE

    PUBLIC :: paw_deqsymmetrize! symmetrize dbecsums for electric field
    !
    PRIVATE

 CONTAINS

SUBROUTINE paw_deqsymmetrize(dbecsum)
!
! This routine similar to PAW_symmetrize, symmetrize the change of 
! dbecsum due to an electric field perturbation. 
!
    USE lsda_mod,          ONLY : nspin
    USE uspp_param,        ONLY : nhm
    USE ions_base,         ONLY : nat, ityp
    USE noncollin_module,  ONLY : domag, nspin_lsda, nspin_mag
    USE cell_base,         ONLY : at, bg
    USE symm_base,         ONLY : nsym, irt, d1, d2, d3, s, t_rev, sname, &
                                  invs, inverse_s
    USE uspp,              ONLY : nhtolm,nhtol,ijtoh
    USE uspp_param,        ONLY : nh, upf
    USE io_global,         ONLY : stdout, ionode

    COMPLEX(DP), INTENT(INOUT) :: dbecsum(nhm*(nhm+1)/2,nat,nspin_mag)! cross band occupations

    COMPLEX(DP)                :: becsym(nhm*(nhm+1)/2,nat,nspin_mag)! symmetrized becsum
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
    INTEGER :: jpol, kpol
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

    becsym(:,:,:) = (0.0_DP,0.0_DP)
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
                          * pref * dbecsum(ouh, ma, is) 
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
       !
       call inverse_s ( )

       becsym(:,:,2:4) = (0._DP,0.0_DP)
       DO ia = 1, nat
          nt = ityp(ia)
          ! No need to symmetrize non-PAW atoms
          IF ( .not. upf(nt)%tpawp ) CYCLE
          !
          !  Bring the magnetization in the basis of the crystal
          !        
          DO ijh=1,(nh(nt)*(nh(nt)+1))/2
             DO jpol=1,3
                mb(jpol)=dbecsum(ijh,ia,jpol+1)
             ENDDO
             DO jpol=1,3
                dbecsum(ijh,ia,jpol+1)=bg(1,jpol)*mb(1) +  &
                                  bg(2,jpol)*mb(2) + bg(3,jpol)*mb(3)
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
                   DO is=1,3
                      DO kpol=1,3
                         becsym(ijh,ia,is+1)=becsym(ijh,ia,is+1) &
                         + D(l_i)%d(m_o,m_i,isym)*D(l_j)%d(m_u,m_j,isym)* &
                           pref*dbecsum(ouh,ma,kpol+1)*&
                           segno*s(kpol,is,invs(isym))
                      END DO
                   END DO   
                END DO ! m_o
                END DO ! m_u
            END DO ! isym
            !
            ! Put the prefactor back in:
            IF ( ih == jh ) becsym(ijh,ia,is) = .5_dp * becsym(ijh,ia,is)
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
            DO jpol=1,3
               mb(jpol)=becsym(ijh,ia,jpol+1)
            ENDDO
            DO jpol=1,3
               becsym(ijh,ia,jpol+1)=at(jpol,1)*mb(1)+at(jpol,2)*mb(2)+&
                                             at(jpol,3)*mb(3)
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
    dbecsum(:,:,:) = becsym(:,:,:)

    CALL stop_clock('PAW_dsymme')

END SUBROUTINE paw_deqsymmetrize

END MODULE paw_add_symmetry
