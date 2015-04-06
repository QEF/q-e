!
! Copyright (C) 2015 GWW group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
    !
  SUBROUTINE adduspos_gamma_r(iw,jw,r_ij,ik,becp_iw,becp_jw)
  !----------------------------------------------------------------------
  !
  !  This routine adds the US term < Psi_iw|r><r|Psi_jw>
  !  to the array r_ij
  !  this is a GAMMA only routine (i.e. r_ij is real)
  !
  USE kinds,                ONLY : DP
  USE realus,               ONLY : tabp
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE fft_base,             ONLY : dfftp
  USE gvect,                ONLY : ngm, nl, nlm, gg, g
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE uspp,                 ONLY : okvan, becsum, nkb, ijtoh, indv_ijkb0
  USE uspp_param,           ONLY : upf, lmaxq, nh
  USE wvfct,                ONLY : wg
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : psic
  USE io_global,            ONLY : stdout
  USE cell_base,            ONLY : omega
  !
  USE mp_bands,         ONLY : intra_bgrp_comm
  USE mp,               ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  !
  INTEGER, INTENT(in) :: iw,jw!the states indices
  REAL(kind=DP), INTENT(inout) :: r_ij(dfftp%nnr)!where to add the us term
  INTEGER, INTENT(in) :: ik!spin index for spin polarized calculations NOT IMPLEMENTED YET
  REAL(kind=DP), INTENT(in) ::  becp_iw( nkb)!overlap of wfcs with us  projectors
  REAL(kind=DP), INTENT(in) ::  becp_jw( nkb)!overlap of wfcs with us  projectors

  !     here the local variables
  !

  INTEGER :: na, nt, nhnt, ir, ih, jh, is , ia, mbia, irb, iqs, sizeqsave
  INTEGER :: ikb, jkb, np
  ! counters

  ! work space for rho(G,nspin)
  ! Fourier transform of q

  IF (.not.okvan) RETURN
  IF( .not.gamma_only) THEN
     WRITE(stdout,*) ' adduspos_gamma_r is a gamma ONLY routine'
     STOP
  ENDIF

  DO is=1,nspin
     !
     DO np = 1, ntyp
        !
        iqs = 0
        !
        IF ( upf(np)%tvanp ) THEN
           !
           DO ia = 1, nat
              !
              mbia = tabp(ia)%maxbox
              nt = ityp(ia)
              nhnt = nh(nt)
              !
              IF ( ityp(ia) /= np ) iqs=iqs+(nhnt+1)*nhnt*mbia/2
              IF ( ityp(ia) /= np ) CYCLE
              !
              DO ih = 1, nhnt
                 !
                 ikb = indv_ijkb0(ia) + ih
                 !
                 DO jh = ih, nhnt
                    !
                    jkb = indv_ijkb0(ia) + jh
                    !
                    DO ir = 1, mbia
                       !
                       irb = tabp(ia)%box(ir)
                       iqs = iqs + 1
                       !
                       r_ij(irb) = r_ij(irb) + tabp(ia)%qr(ir,ijtoh(ih,jh,np))&
                                              *becp_iw(ikb)*becp_jw(jkb)*omega
                       !
                       IF ( ih /= jh ) THEN
                          r_ij(irb) = r_ij(irb) + tabp(ia)%qr(ir,ijtoh(ih,jh,np))&
                                                 *becp_iw(jkb)*becp_jw(ikb)*omega
                       ENDIF
                    ENDDO
                 ENDDO
              ENDDO
              !
           ENDDO
           !
        ENDIF
     ENDDO
  ENDDO
  !
  RETURN
  !
  END SUBROUTINE adduspos_gamma_r
  !
  SUBROUTINE adduspos_r(r_ij,becp_iw,becp_jw)
  !----------------------------------------------------------------------
  !
  !  This routine adds the US term < Psi_iw|r><r|Psi_jw>
  !  to the array r_ij
  USE kinds,                ONLY : DP
  USE realus,               ONLY : tabp
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE fft_base,             ONLY : dfftp
  USE gvect,                ONLY : ngm, nl, nlm, gg, g
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE uspp,                 ONLY : okvan, becsum, nkb, ijtoh, indv_ijkb0
  USE uspp_param,           ONLY : upf, lmaxq, nh
  USE wvfct,                ONLY : wg
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : psic
  USE cell_base,            ONLY : omega
  !
  IMPLICIT NONE
  !
  COMPLEX(kind=DP), INTENT(inout) :: r_ij(dfftp%nnr)!where to add the us term
  COMPLEX(kind=DP), INTENT(in) ::  becp_iw( nkb)!overlap of wfcs with us  projectors
  COMPLEX(kind=DP), INTENT(in) ::  becp_jw( nkb)!overlap of wfcs with us  projectors
  !     here the local variables
  !
  INTEGER :: na, ia, nt, nhnt, ir, ih, jh, is, mbia, irb, iqs
  INTEGER :: ikb, jkb, np
  ! counters

  ! work space for rho(G,nspin)
  ! Fourier transform of q

  IF (.not.okvan) RETURN

  DO is=1,nspin
     !
     DO np = 1, ntyp
        !
        iqs = 0
        !
        IF ( upf(np)%tvanp ) THEN
           !
           DO ia = 1, nat
              !
              mbia = tabp(ia)%maxbox
              nt = ityp(ia)
              nhnt = nh(nt)
              !
              IF ( ityp(ia) /= np ) iqs=iqs+(nhnt+1)*nhnt*mbia/2
              IF ( ityp(ia) /= np ) CYCLE
              !
              DO ih = 1, nhnt
                 !
                 ikb = indv_ijkb0(ia) + ih
                 DO jh = ih, nhnt
                    !
                    jkb = indv_ijkb0(ia) + jh
                    !
                    DO ir = 1, mbia
                       !
                       irb = tabp(ia)%box(ir)
                       iqs = iqs + 1
                       !
                       r_ij(irb) = r_ij(irb) + tabp(ia)%qr(ir,ijtoh(ih,jh,np))&
                                              *conjg(becp_iw(ikb))*becp_jw(jkb)*omega
                       !
                       IF ( ih /= jh ) THEN
                          r_ij(irb) = r_ij(irb) + tabp(ia)%qr(ir,ijtoh(ih,jh,np))&
                                                 *conjg(becp_iw(jkb))*becp_jw(ikb)*omega
                       ENDIF
                    ENDDO
                 ENDDO
              ENDDO
              !
           ENDDO
           !
        ENDIF
     ENDDO
  ENDDO
  !
  RETURN
  END SUBROUTINE adduspos_r
  !
  SUBROUTINE adduspos_real(sca,qq_op,becp_iw,becp_jw)
  !----------------------------------------------------------------------
  !
  !  This routine adds the US term < Psi_iw|r><r|Psi_jw>
  !  to the array r_ij
  USE kinds,                ONLY : DP
  USE realus,               ONLY : tabp
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE gvect,                ONLY : ngm, nl, nlm, gg, g
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE uspp,                 ONLY : okvan, becsum, nkb, qq, indv_ijkb0
  USE uspp_param,           ONLY : upf, lmaxq, nh, nhm
  USE wvfct,                ONLY : wg
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : psic
  USE cell_base,            ONLY : omega
  !
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  REAL(kind=DP), INTENT(inout) :: sca!where to add the us term
  REAL(kind=DP), INTENT(in) ::  becp_iw( nkb)!overlap of wfcs with us  projectors
  REAL(kind=DP), INTENT(in) ::  becp_jw( nkb)!overlap of wfcs with us  projectors
  REAL(kind=DP), INTENT(in) ::  qq_op(nhm, nhm,nat)!US augmentation charges

  !     here the local variables
  !

  INTEGER :: na, ia, nhnt, nt, ih, jh, is, mbia
  INTEGER :: ikb, jkb, np
  ! counters

  ! work space for rho(G,nspin)
  ! Fourier transform of q

  IF (.not.okvan) RETURN

  DO is=1,nspin
     !
     DO np = 1, ntyp
        !
        IF ( upf(np)%tvanp ) THEN
           !
           DO ia = 1, nat
              !
              IF ( ityp(ia) /= np ) CYCLE
              !
              mbia = tabp(ia)%maxbox
              nt = ityp(ia)
              nhnt = nh(nt)
              !
              DO ih = 1, nhnt
                 !
                 ikb = indv_ijkb0(ia) + ih
                 DO jh = ih, nhnt
                    !
                    jkb = indv_ijkb0(ia) + jh
                    !
                    sca = sca + qq_op(ih,jh,ia) * becp_iw(ikb)*becp_jw(jkb)
                    !
                    IF ( ih /= jh ) THEN
                       sca = sca + qq_op(jh,ih,ia) * becp_iw(ikb)*becp_jw(jkb)
                    ENDIF
                    !
                 ENDDO
              ENDDO
              !
           ENDDO
           !
        ENDIF
     ENDDO
  ENDDO
  !
  RETURN
  !
  END SUBROUTINE adduspos_real
  !
