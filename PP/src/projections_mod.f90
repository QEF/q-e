!
! Copyright (C) 2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE projections
  USE kinds, ONLY : DP
  
  TYPE wfc_label
     INTEGER na, n, l, m, ind
     REAL (DP) jj
     CHARACTER(LEN=2) els
  END TYPE wfc_label
  TYPE(wfc_label), ALLOCATABLE :: nlmchi(:)
  
  REAL (DP),    ALLOCATABLE :: proj (:,:,:)
  COMPLEX (DP), ALLOCATABLE :: proj_aux (:,:,:)
  COMPLEX (DP), ALLOCATABLE :: ovps_aux (:,:,:)
  
  CONTAINS
    !
    SUBROUTINE fill_nlmchi ( natomwfc, lmax_wfc )
      !
      USE ions_base,  ONLY : ityp, nat
      USE upf_ions,   ONLY : n_atom_wfc
      USE upf_utils,  ONLY : l_to_spdf
      USE uspp_param, ONLY: upf
      USE noncollin_module, ONLY: noncolin, lspinorb
      !
      IMPLICIT NONE
      INTEGER, INTENT (OUT) :: natomwfc, lmax_wfc
      !
      INTEGER :: nwfc, na, nt, n, n1, n2, l, m, ind
      REAL(dp) :: jj, fact(2)
      REAL(dp), EXTERNAL :: spinor
      CHARACTER(LEN=2) :: label
      INTEGER :: nn(0:3)
      !
      natomwfc = n_atom_wfc (nat, ityp, noncolin)
      !
      ALLOCATE (nlmchi(natomwfc))
      nwfc=0
      lmax_wfc = 0
      DO na = 1, nat
         nt = ityp (na)
         n2 = 0
         nn = [1,2,3,4]
         DO n = 1, upf(nt)%nwfc
            IF (upf(nt)%oc (n) >= 0.d0) THEN
               label = upf(nt)%els(n)
               l = upf(nt)%lchi (n)
               ! the following lines guess the label if absent
               IF ( label =='Xn' ) THEN
                  WRITE(label,'(I1,A1)') nn(l), l_to_spdf(l)
                  nn(l) = nn(l) + 1
               END IF
               lmax_wfc = max (lmax_wfc, l )
               IF (lspinorb) THEN
                  IF (upf(nt)%has_so) THEN
                     jj = upf(nt)%jchi (n)
                     ind = 0
                     DO m = -l-1, l
                        fact(1) = spinor(l,jj,m,1)
                        fact(2) = spinor(l,jj,m,2)
                        IF (abs(fact(1)) > 1.d-8 .or. abs(fact(2)) > 1.d-8) THEN
                           nwfc = nwfc + 1
                           ind = ind + 1
                           nlmchi(nwfc)%na = na
                           nlmchi(nwfc)%n  =  n
                           nlmchi(nwfc)%l  =  l
                           nlmchi(nwfc)%m  =  m
                           nlmchi(nwfc)%ind = ind
                           nlmchi(nwfc)%jj  = jj
                           nlmchi(nwfc)%els = label
                        ENDIF
                     ENDDO
                  ELSE
                     DO n1 = l, l+1
                        jj= dble(n1) - 0.5d0
                        ind = 0
                        IF (jj>0.d0)  THEN
                           n2 = n2 + 1
                           DO m = -l-1, l
                              fact(1) = spinor(l,jj,m,1)
                              fact(2) = spinor(l,jj,m,2)
                              IF (abs(fact(1)) > 1.d-8 .or. abs(fact(2)) > 1.d-8) THEN
                                 nwfc = nwfc + 1
                                 ind = ind + 1
                                 nlmchi(nwfc)%na = na
                                 nlmchi(nwfc)%n  =  n2
                                 nlmchi(nwfc)%l  =  l
                                 nlmchi(nwfc)%m  =  m
                                 nlmchi(nwfc)%ind = ind
                                 nlmchi(nwfc)%jj  = jj
                                 nlmchi(nwfc)%els = label
                              ENDIF
                           ENDDO
                        ENDIF
                     ENDDO
                  ENDIF
               ELSE
                  DO m = 1, 2 * l + 1
                     nwfc=nwfc+1
                     nlmchi(nwfc)%na = na
                     nlmchi(nwfc)%n  =  n
                     nlmchi(nwfc)%l  =  l
                     nlmchi(nwfc)%m  =  m
                     nlmchi(nwfc)%ind=  m
                     nlmchi(nwfc)%jj =  0.d0
                     nlmchi(nwfc)%els=  label
                  ENDDO
                  IF ( noncolin) THEN
                     DO m = 1, 2 * l + 1
                        nwfc=nwfc+1
                        nlmchi(nwfc)%na = na
                        nlmchi(nwfc)%n  =  n
                        nlmchi(nwfc)%l  =  l
                        nlmchi(nwfc)%m  =  m
                        nlmchi(nwfc)%ind=  m+2*l+1
                        nlmchi(nwfc)%jj =  0.d0
                        nlmchi(nwfc)%els = label
                     END DO
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      !
      IF (lmax_wfc > 3) CALL errore ('fill_nlmchi', 'l > 3 not yet implemented',1)
      IF (nwfc /= natomwfc) CALL errore ('fill_nlmchi','wrong # of atomic wfcs',1)
      
    END SUBROUTINE fill_nlmchi
    !
    SUBROUTINE fill_nlmbeta( nkb, nwfc )
      !
      USE ions_base, ONLY : ityp, nat, ntyp => nsp
      USE uspp_param, ONLY: upf
      !
      IMPLICIT NONE
      INTEGER, INTENT(in)  :: nkb
      INTEGER, INTENT(out) :: nwfc
      !
      INTEGER :: na, nt, nb, l, m
      !
      ALLOCATE (nlmchi(nkb))
      nwfc=0
      DO nt=1,ntyp
         DO na=1,nat
            IF (ityp(na)==nt) THEN
               DO nb=1,upf(nt)%nbeta
                  l=upf(nt)%lll(nb)
                  DO m = 1, 2 * l + 1
                     nwfc=nwfc+1
                     nlmchi(nwfc)%na = na
                     nlmchi(nwfc)%n =  nb
                     nlmchi(nwfc)%l = l
                     nlmchi(nwfc)%m = m
                     nlmchi(nwfc)%ind = m
                     nlmchi(nwfc)%jj  =  0.d0
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      !
    END SUBROUTINE fill_nlmbeta
    !
    !-----------------------------------------------------------------------
    SUBROUTINE sym_proj_g (rproj0, proj_out)
      !-----------------------------------------------------------------------
      !
      USE kinds,      ONLY : DP
      USE basis,      ONLY : natomwfc
      USE wvfct,      ONLY : nbnd
      USE symm_base,  ONLY : nsym, irt, t_rev, d1, d2, d3
      !
      IMPLICIT NONE
      REAL(DP),    INTENT(IN) ::rproj0 (natomwfc, nbnd)
      REAL   (DP), INTENT(OUT):: proj_out(natomwfc, nbnd)
      !
      INTEGER :: na, nb, n, l, m, m1, isym, nwfc, nwfc1, ibnd
      REAL   (DP), ALLOCATABLE :: rwork1(:)
      !
      ! initialize D_Sl for l=1, l=2 and l=3, for l=0 D_S0 is 1
      !
      CALL d_matrix (d1, d2, d3)
      proj_out(:,:) = 0.d0
      !
      ALLOCATE(rwork1 (nbnd) )
      !
      DO nwfc = 1, natomwfc
         !
         !  atomic wavefunction nwfc is on atom na
         !
         na= nlmchi(nwfc)%na
         n = nlmchi(nwfc)%n
         l = nlmchi(nwfc)%l
         m = nlmchi(nwfc)%m
         !
         DO isym = 1, nsym
            nb = irt (isym, na)
            DO nwfc1 =1, natomwfc
               IF (nlmchi(nwfc1)%na == nb             .and. &
                    nlmchi(nwfc1)%n == nlmchi(nwfc)%n .and. &
                    nlmchi(nwfc1)%l == nlmchi(nwfc)%l .and. &
                    nlmchi(nwfc1)%m == 1 ) GOTO 10
            ENDDO
            CALL errore('sym_proj_g','cannot symmetrize',1)
10          nwfc1=nwfc1-1
            !
            !  nwfc1 is the first rotated atomic wfc corresponding to nwfc
            !
            IF (l == 0) THEN
               rwork1(:) = rproj0 (nwfc1 + 1,:)
            ELSEIF (l == 1) THEN
               rwork1(:) = 0.d0
               DO m1 = 1, 3
                  rwork1(:)=rwork1(:)+d1(m1,m,isym)*rproj0(nwfc1+m1,:)
               ENDDO
            ELSEIF (l == 2) THEN
               rwork1(:) = 0.d0
               DO m1 = 1, 5
                  rwork1(:)=rwork1(:)+d2(m1,m,isym)*rproj0(nwfc1+m1,:)
               ENDDO
            ELSEIF (l == 3) THEN
               rwork1(:) = 0.d0
               DO m1 = 1, 7
                  rwork1(:)=rwork1(:)+d3(m1,m,isym)*rproj0(nwfc1+m1,:)
               ENDDO
            ENDIF
            DO ibnd = 1, nbnd
               proj_out (nwfc, ibnd ) = proj_out (nwfc, ibnd ) + &
                    rwork1(ibnd) * rwork1(ibnd) / nsym
            ENDDO
         ENDDO
      ENDDO
      !
      DEALLOCATE(rwork1 )
      !  
    END SUBROUTINE sym_proj_g
    !
    !-----------------------------------------------------------------------
    SUBROUTINE sym_proj_k (proj0, proj_out)
      !-----------------------------------------------------------------------
      !
      USE kinds,      ONLY : DP
      USE basis,      ONLY : natomwfc
      USE wvfct,      ONLY : nbnd
      USE symm_base,  ONLY : nsym, irt, t_rev, d1, d2, d3
      !
      IMPLICIT NONE
      COMPLEX(DP), INTENT(IN) :: proj0 (natomwfc, nbnd)
      REAL   (DP), INTENT(OUT):: proj_out(natomwfc, nbnd)
      !
      INTEGER :: na, nb, n, l, m, m1, isym, nwfc, nwfc1, ibnd
      COMPLEX(DP), ALLOCATABLE ::  work1(:)
      !
      ! initialize D_Sl for l=1, l=2 and l=3, for l=0 D_S0 is 1
      !
      CALL d_matrix (d1, d2, d3)
      proj_out(:,:) = 0.d0
      !
      ALLOCATE(work1 (nbnd) )
      !
      DO nwfc = 1, natomwfc
         !
         !  atomic wavefunction nwfc is on atom na
         !
         na= nlmchi(nwfc)%na
         n = nlmchi(nwfc)%n
         l = nlmchi(nwfc)%l
         m = nlmchi(nwfc)%m
         !
         DO isym = 1, nsym
            nb = irt (isym, na)
            DO nwfc1 =1, natomwfc
               IF (nlmchi(nwfc1)%na == nb             .and. &
                    nlmchi(nwfc1)%n == nlmchi(nwfc)%n .and. &
                    nlmchi(nwfc1)%l == nlmchi(nwfc)%l .and. &
                    nlmchi(nwfc1)%m == 1 ) GOTO 10
            ENDDO
            CALL errore('sym_proj_k','cannot symmetrize',1)
10          nwfc1=nwfc1-1
            !
            !  nwfc1 is the first rotated atomic wfc corresponding to nwfc
            !
            IF (l == 0) THEN
               work1(:) = proj0 (nwfc1 + 1,:)
            ELSEIF (l == 1) THEN
               work1(:) = 0.d0
               DO m1 = 1, 3
                  work1(:)=work1(:)+d1(m1,m,isym)*proj0(nwfc1+m1,:)
               ENDDO
            ELSEIF (l == 2) THEN
               work1(:) = 0.d0
               DO m1 = 1, 5
                  work1(:)=work1(:)+d2(m1,m,isym)*proj0(nwfc1+m1,:)
               ENDDO
            ELSEIF (l == 3) THEN
               work1(:) = 0.d0
               DO m1 = 1, 7
                  work1(:)=work1(:)+d3(m1,m,isym)*proj0(nwfc1+m1,:)
               ENDDO
            ENDIF
            DO ibnd = 1, nbnd
               proj_out (nwfc, ibnd) = proj_out (nwfc, ibnd) + &
                    work1(ibnd) * conjg (work1(ibnd)) / nsym
            ENDDO
         ENDDO
      ENDDO
      !
      DEALLOCATE(work1 )
      !  
    END SUBROUTINE sym_proj_k
    !
    !-----------------------------------------------------------------------
    SUBROUTINE sym_proj_so (domag, proj0, proj_out  )
      !-----------------------------------------------------------------------
      !
      USE kinds,      ONLY : DP
      USE basis,      ONLY : natomwfc
      USE wvfct,      ONLY : nbnd
      USE symm_base,  ONLY : nsym, irt, t_rev
      !
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: domag
      COMPLEX(DP), INTENT(IN) :: proj0 (natomwfc, nbnd)
      REAL   (DP), INTENT(OUT):: proj_out(natomwfc, nbnd)
      !
      INTEGER :: na, nb, n, l, m, m1, ind, ind0, isym, nwfc, nwfc1, ibnd
      REAL(DP) :: jj
      COMPLEX(DP), ALLOCATABLE ::  work1(:)
      COMPLEX(DP) :: d12(2, 2, 48), d32(4, 4, 48), d52(6, 6, 48), d72(8, 8, 48)
      !
      ! initialize D_Sj for j=1/2, j=3/2, j=5/2 and j=7/2
      !
      CALL d_matrix_so (d12, d32, d52, d72)
      !
      proj_out(:,:) = 0.d0
      !
      ALLOCATE(work1 (nbnd) )
      !
      DO nwfc = 1, natomwfc
         !
         !  atomic wavefunction nwfc is on atom na
         !
         na= nlmchi(nwfc)%na
         n = nlmchi(nwfc)%n
         l = nlmchi(nwfc)%l
         m = nlmchi(nwfc)%m
         ind0 = nlmchi(nwfc)%ind
         jj = nlmchi(nwfc)%jj
         !
         DO isym = 1, nsym
            !-- check for the time reversal
            IF (t_rev(isym) == 1) THEN
               ind = 2*jj + 2 - ind0
            ELSE
               ind = ind0
            ENDIF
            !--
            nb = irt (isym, na)
            DO nwfc1 =1, natomwfc
               IF (nlmchi(nwfc1)%na == nb             .and. &
                    nlmchi(nwfc1)%n == nlmchi(nwfc)%n .and. &
                    nlmchi(nwfc1)%l == nlmchi(nwfc)%l .and. &
                    nlmchi(nwfc1)%jj == nlmchi(nwfc)%jj .and. &
                    nlmchi(nwfc1)%ind == 1 ) GOTO 10
            ENDDO
            CALL errore('sym_proj_so','cannot symmetrize',1)
10          nwfc1=nwfc1-1
            !
            !  nwfc1 is the first rotated atomic wfc corresponding to nwfc
            !
            IF (abs(jj-0.5d0)<1.d-8) THEN
               work1(:) = 0.d0
               DO m1 = 1, 2
                  work1(:)=work1(:)+d12(m1,ind,isym)*proj0(nwfc1+m1,:)
               ENDDO
            ELSEIF (abs(jj-1.5d0)<1.d-8) THEN
               work1(:) = 0.d0
               DO m1 = 1, 4
                  work1(:)=work1(:)+d32(m1,ind,isym)*proj0(nwfc1 + m1,:)
               ENDDO
            ELSEIF (abs(jj-2.5d0)<1.d-8) THEN
               work1(:) = 0.d0
               DO m1 = 1, 6
                  work1(:)=work1(:)+d52(m1,ind,isym)*proj0(nwfc1+m1,:)
               ENDDO
            ELSEIF (abs(jj-3.5d0)<1.d-8) THEN
               work1(:) = 0.d0
               DO m1 = 1, 8
                  work1(:)=work1(:)+d72(m1,ind,isym)*proj0(nwfc1+m1,:)
               ENDDO
            ENDIF
            DO ibnd = 1, nbnd
               proj_out (nwfc, ibnd) = proj_out (nwfc, ibnd) + &
                    work1(ibnd) * conjg (work1(ibnd)) / nsym
            ENDDO
            ! on symmetries
            !--  in a nonmagnetic case - another loop with the time reversal
            IF ( .not.domag .and. ind==ind0 ) THEN
               ind = 2*jj + 2 - ind0
               nwfc1 = nwfc1 + 1
               GOTO 10
            ENDIF
         ENDDO
         !--  in a nonmagnetic case - rescale
         IF (.not.domag) THEN
            DO ibnd = 1, nbnd
               proj_out(nwfc,ibnd) = 0.5_dp*proj_out(nwfc,ibnd)
            ENDDO
         ENDIF
         ! on atomic wavefunctions
      END DO
      !
      DEALLOCATE (work1)
      !
    END SUBROUTINE sym_proj_so
    !-----------------------------------------------------------------------
    SUBROUTINE sym_proj_nc ( proj0, proj_out  )
      !
      USE kinds,      ONLY : DP
      USE basis,      ONLY : natomwfc
      USE wvfct,      ONLY : nbnd
      USE symm_base,  ONLY : nsym, irt, t_rev
      !
      IMPLICIT NONE
      COMPLEX(DP), INTENT(IN) :: proj0 (natomwfc, nbnd)
      REAL   (DP), INTENT(OUT):: proj_out(natomwfc, nbnd)
      !
      INTEGER :: na, nb, n, l, m, m1, ind, ind0, jj, isym, nwfc, nwfc1, ibnd
      COMPLEX(DP), ALLOCATABLE ::  work1(:)
      COMPLEX(DP) :: d012(2, 2, 48), d112(6, 6, 48), d212(10, 10, 48), &
           d312(14, 14, 48)
      !
      ! initialize D_Sl for l=0, l=1, l=2 and l=3
      !
      CALL d_matrix_nc (d012, d112, d212, d312)
      !
      proj_out(:,:) = 0.d0
      !
      ALLOCATE(work1 (nbnd) )
      !
      DO nwfc = 1, natomwfc
         na= nlmchi(nwfc)%na
         n = nlmchi(nwfc)%n
         l = nlmchi(nwfc)%l
         m = nlmchi(nwfc)%m
         ind0 = nlmchi(nwfc)%ind
         !
         DO isym = 1, nsym
            !-- check for the time reversal
            IF (t_rev(isym) == 1) THEN
               ind = 2*m - ind0 + 2*l + 1
            ELSE
               ind = ind0
            ENDIF
            !--
            nb = irt (isym, na)
            DO nwfc1 =1, natomwfc
               IF (nlmchi(nwfc1)%na == nb             .and. &
                    nlmchi(nwfc1)%n == nlmchi(nwfc)%n .and. &
                    nlmchi(nwfc1)%l == nlmchi(nwfc)%l .and. &
                    nlmchi(nwfc1)%m == 1 .and. &
                    nlmchi(nwfc1)%ind == 1) GOTO 15
            ENDDO
            CALL errore('sym_proj_nc','cannot symmetrize',1)
15          nwfc1=nwfc1-1
            IF (l == 0) THEN
               work1(:) = 0.d0
               DO m1 = 1, 2
                  work1(:) = work1(:) + d012 (m1, ind, isym) * &
                       proj0 (nwfc1 + m1,:)
               ENDDO
            ELSEIF (l == 1) THEN
               work1(:) = 0.d0
               DO m1 = 1, 6
                  work1(:) = work1(:) + d112 (m1, ind, isym) * &
                       proj0 (nwfc1 + m1,:)
               ENDDO
            ELSEIF (l == 2) THEN
               work1(:) = 0.d0
               DO m1 = 1, 10
                  work1(:) = work1(:) + d212 (m1, ind, isym) * &
                       proj0 (nwfc1 + m1,:)
               ENDDO
            ELSEIF (l == 3) THEN
               work1(:) = 0.d0
               DO m1 = 1, 14
                  work1(:) = work1(:) + d312 (m1, ind, isym) * &
                       proj0 (nwfc1 + m1,:)
               ENDDO
            ENDIF
            DO ibnd = 1, nbnd
               proj_out (nwfc, ibnd) = proj_out (nwfc, ibnd) + &
                    work1(ibnd) * conjg (work1(ibnd)) / nsym
            ENDDO
            ! on symmetries
         ENDDO
         ! on atomic wavefunctions
      END DO
      !
      DEALLOCATE (work1)
      !
    END SUBROUTINE sym_proj_nc
    !
    !-----------------------------------------------------------------------
    FUNCTION compute_mj(j,l,m)
      !-----------------------------------------------------------------------
      USE kinds, ONLY: DP
      IMPLICIT NONE
      !
      REAL(DP) :: compute_mj, j
      INTEGER  :: l, m

      IF (abs(j-l-0.5d0)<1.d-4) THEN
          compute_mj=m+0.5d0
      ELSEIF (abs(j-l+0.5d0)<1.d-4) THEN
        compute_mj=m-0.5d0
      ELSE
        CALL errore('compute_mj','l and j not compatible',1)
      ENDIF

      RETURN
    END FUNCTION compute_mj
    !
    SUBROUTINE compute_zdistmat( npw, n, nx, v, w, dm, &
      idesc, rank_ip, idesc_ip )
      !
      !  This subroutine compute <vi|wj> and store the
      !  result in distributed matrix dm
      !
      USE mp, ONLY : mp_root_sum
      USE mp_pools, ONLY: intra_pool_comm
      !
      IMPLICIT NONE
      !
      include 'laxlib.fh'
      !
      INTEGER, INTENT(in) :: idesc(LAX_DESC_SIZE)
      INTEGER, INTENT(in) :: rank_ip( :, : )
      INTEGER, INTENT(in) :: idesc_ip( :, :, : )
      !
      INTEGER :: ipc, ipr
      INTEGER :: nr, nc, ir, ic, root, ldv, ldw
      INTEGER, INTENT(in) :: npw ! local number of plane wave
      INTEGER, INTENT(in) :: n   ! global dimension of matrix dm
      INTEGER, INTENT(in) :: nx  ! local leading dimension of matrix dm
                                 ! WARNING: nx is the same on all proc, SIZE( dm, 1 ) NO!
      COMPLEX(DP), INTENT(out) :: dm( :, : )
      COMPLEX(DP) :: v(:,:), w(:,:)
      COMPLEX(DP), ALLOCATABLE :: work( :, : )
      !
      ALLOCATE( work( nx, nx ) )
      !
      work = (0.0_dp, 0.0_dp)
      !
      ldv = size( v, 1 )
      ldw = size( w, 1 )
      !
      DO ipc = 1, idesc(LAX_DESC_NPC) !  loop on column procs
         !
         nc = idesc_ip( LAX_DESC_NC, 1, ipc )
         ic = idesc_ip( LAX_DESC_IC, 1, ipc )
         !
         DO ipr = 1, ipc ! desc( la_npr_ ) ! ipc ! use symmetry for the loop on row procs
            !
            nr = idesc_ip( LAX_DESC_NR, ipr, ipc )
            ir = idesc_ip( LAX_DESC_IR, ipr, ipc )
            !
            !  rank of the processor for which this block (ipr,ipc) is destinated
            !
            root = rank_ip( ipr, ipc )
 
            ! use blas subs. on the matrix block
 
            CALL ZGEMM( 'C', 'N', nr, nc, npw, (1.0_dp,0.0_dp) , &
                        v(1,ir), ldv, w(1,ic), ldw, (0.0_dp,0.0_dp), work, nx )
 
            ! accumulate result on dm of root proc.
 
            CALL mp_root_sum( work, dm, root, intra_pool_comm )
 
         ENDDO
         !
      ENDDO
      !
      CALL laxlib_zsqmher( n, dm, nx, idesc )
      !
      DEALLOCATE( work )
      !
      RETURN
    END SUBROUTINE compute_zdistmat
    !
    SUBROUTINE compute_ddistmat( npw, n, nx, v, w, dm, &
      idesc, rank_ip, idesc_ip )
      !
      !  This subroutine compute <vi|wj> and store the
      !  result in distributed matrix dm
      !
      USE mp, ONLY : mp_root_sum
      USE mp_pools, ONLY: intra_pool_comm
      USE gvect, ONLY : gstart
      USE wvfct, ONLY : npwx
      !
      IMPLICIT NONE
      !
      include 'laxlib.fh'
      !
      INTEGER, INTENT(in) :: idesc(LAX_DESC_SIZE)
      INTEGER, INTENT(in) :: rank_ip( :, : )
      INTEGER, INTENT(in) :: idesc_ip( :, :, : )
      !
      INTEGER :: ipc, ipr
      INTEGER :: nr, nc, ir, ic, root, ldv, ldw, npw2, npwx2
      INTEGER, INTENT(in) :: npw ! local number of plane wave
      INTEGER, INTENT(in) :: n   ! global dimension of matrix dm
      INTEGER, INTENT(in) :: nx  ! local leading dimension of matrix dm
                                ! WARNING: nx is the same on all proc, SIZE( dm, 1 ) NO!
      REAL(DP), INTENT(out) :: dm( :, : )
      COMPLEX(DP) :: v(:,:), w(:,:)
      REAL(DP), ALLOCATABLE :: work( :, : )
      !
      ALLOCATE( work( nx, nx ) )
      !
      npw2  = 2*npw
      npwx2 = 2*npwx
      !
      work = (0.0_dp, 0.0_dp)
      !
      ldv = size( v, 1 )
      ldw = size( w, 1 )
      !
      DO ipc = 1, idesc(LAX_DESC_NPC) !  loop on column procs
        !
        nc = idesc_ip( LAX_DESC_NC, 1, ipc )
        ic = idesc_ip( LAX_DESC_IC, 1, ipc )
        !
        DO ipr = 1, ipc ! desc( la_npr_ ) ! ipc ! use symmetry for the loop on row procs
            !
            nr = idesc_ip( LAX_DESC_NR, ipr, ipc )
            ir = idesc_ip( LAX_DESC_IR, ipr, ipc )
            !
            !  rank of the processor for which this block (ipr,ipc) is destinated
            !
            root = rank_ip( ipr, ipc )

            ! use blas subs. on the matrix block

            ! use blas subs. on the matrix block

            CALL DGEMM( 'T', 'N', nr, nc, npw2, 2.D0 , &
                        v(1,ir), npwx2, w(1,ic), npwx2, 0.D0, work, nx )

            IF ( gstart == 2 ) &
              CALL DGER( nr, nc, -1.D0, v(1,ir), npwx2, w(1,ic), npwx2, work, nx )

            ! accumulate result on dm of root proc.

            CALL mp_root_sum( work, dm, root, intra_pool_comm )

        ENDDO
        !
      ENDDO
      !
      CALL laxlib_dsqmsym( n, dm, nx, idesc )
      !
      DEALLOCATE( work )
      !
      RETURN
    END SUBROUTINE compute_ddistmat
    !
    SUBROUTINE wf_times_overlap( nx, npw, swfc, ovr, wfc, &
      idesc, rank_ip, idesc_ip, la_proc )
      ! Compute |wfc> = |swfc> ovr
      ! on input swfc is a total matrix, ovr is a distributed matrix
      ! on output wfc is a total matrix
      !
      USE mp_pools, ONLY: intra_pool_comm
      USE mp, ONLY: mp_bcast
      !
      IMPLICIT NONE
      !
      include 'laxlib.fh'
      !
      INTEGER, INTENT(in) :: idesc(LAX_DESC_SIZE)
      INTEGER, INTENT(in) :: rank_ip( :, : )
      INTEGER, INTENT(in) :: idesc_ip( :, :, : )
      LOGICAL, INTENT(in) :: la_proc
      !
      INTEGER, INTENT(in) :: nx, npw
      COMPLEX(DP) :: swfc( :, : ), ovr( :, : ), wfc( :, : )
      !
      INTEGER :: npwx
      INTEGER :: ipc, ipr
      INTEGER :: nr, nc, ir, ic, root
      COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )
      COMPLEX(DP) :: beta

      ALLOCATE( vtmp( nx, nx ) )
      !
      npwx = SIZE(swfc,1)
      DO ipc = 1, idesc(LAX_DESC_NPC) !  loop on column procs
        !
        nc = idesc_ip( LAX_DESC_NC, 1, ipc )
        ic = idesc_ip( LAX_DESC_IC, 1, ipc )
        !
        beta = (0.0_dp, 0.0_dp)

        DO ipr = 1, idesc(LAX_DESC_NPR)
            !
            nr = idesc_ip( LAX_DESC_NR, ipr, ipc )
            ir = idesc_ip( LAX_DESC_IR, ipr, ipc )
            !
            root = rank_ip( ipr, ipc )

            IF( ipr-1 == idesc(LAX_DESC_MYR) .and. ipc-1 == idesc(LAX_DESC_MYC) .and. la_proc ) THEN
              !
              !  this proc sends his block
              !
              CALL mp_bcast( ovr, root, intra_pool_comm )
              CALL ZGEMM( 'N', 'N', npw, nc, nr, (1.0_dp,0.0_dp), &
                          swfc(1,ir), npwx, ovr, nx, beta, wfc(1,ic), npwx )
            ELSE
              !
              !  all other procs receive
              !
              CALL mp_bcast( vtmp, root, intra_pool_comm )
              CALL ZGEMM( 'N', 'N', npw, nc, nr, (1.0_dp,0.0_dp), &
                        swfc(1,ir), npwx, vtmp, nx, beta, wfc(1,ic), npwx )
            ENDIF
            !
            beta = (1.0_dp,0.0_dp)

        ENDDO
        !
      ENDDO
      !
      DEALLOCATE( vtmp )

      RETURN
  
    END SUBROUTINE wf_times_overlap
  
    !
    SUBROUTINE wf_times_roverlap( nx, npw, swfc, ovr, wfc, &
      idesc, rank_ip, idesc_ip, la_proc )
      !
      USE gvect, ONLY : gstart
      USE mp_pools, ONLY : intra_pool_comm
      USE mp, ONLY : mp_bcast
      !
      IMPLICIT NONE
      !
      include 'laxlib.fh'
      !
      INTEGER, INTENT(in) :: idesc(LAX_DESC_SIZE)
      INTEGER, INTENT(in) :: rank_ip( :, : )
      INTEGER, INTENT(in) :: idesc_ip( :, :, : )
      LOGICAL, INTENT(in) :: la_proc
      !
      INTEGER, INTENT(in) :: nx, npw
      COMPLEX(DP) :: swfc( :, : ), wfc( :, : )
      REAL(DP)    :: ovr( :, : )
      !
      INTEGER :: ipc, ipr, npw2, npwx2
      INTEGER :: nr, nc, ir, ic, root
      REAL(DP), ALLOCATABLE :: vtmp( :, : )
      REAL(DP) :: beta

      npw2  = 2*npw
      npwx2 = 2*SIZE(swfc,1)

      ALLOCATE( vtmp( nx, nx ) )
      !
      DO ipc = 1, idesc(LAX_DESC_NPC) !  loop on column procs
        !
        nc = idesc_ip( LAX_DESC_NC, 1, ipc )
        ic = idesc_ip( LAX_DESC_IC, 1, ipc )
        !
        beta = 0.0d0

        DO ipr = 1, idesc(LAX_DESC_NPR)
            !
            nr = idesc_ip( LAX_DESC_NR, ipr, ipc )
            ir = idesc_ip( LAX_DESC_IR, ipr, ipc )
            !
            root = rank_ip( ipr, ipc )

            IF( ipr-1 == idesc(LAX_DESC_MYR) .and. ipc-1 == idesc(LAX_DESC_MYC) .and. la_proc ) THEN
              !
              !  this proc sends his block
              !
              CALL mp_bcast( ovr, root, intra_pool_comm )
              CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0, &
                          swfc(1,ir), npwx2, ovr, nx, beta, wfc(1,ic), npwx2 )
              !
            ELSE
              !
              !  all other procs receive
              !
              CALL mp_bcast( vtmp, root, intra_pool_comm )
              CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0, &
                          swfc(1,ir), npwx2, vtmp, nx, beta, wfc(1,ic), npwx2 )
              !
            ENDIF
            !
            beta = 1.0d0

        ENDDO
        !
      ENDDO
      !
      DEALLOCATE( vtmp )

      RETURN
  
    END SUBROUTINE wf_times_roverlap

END MODULE projections
