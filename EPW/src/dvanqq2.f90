  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  !
  ! Copyright (C) 2001 PWSCF group
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  ! Adapted from PH/dvanqq (QE)
  ! 
  !----------------------------------------------------------------------
  SUBROUTINE dvanqq2
  !----------------------------------------------------------------------
  !!
  !! New
  !! This routine calculates two integrals of the Q functions and
  !! its derivatives with c V_loc and V_eff which are used
  !! to compute term dV_bare/dtau * psi  in addusdvqpsi.
  !! The result is stored in int1,int2. The routine is called
  !! for each q in nqc. 
  !! 
  !! RM - Nov/Dec 2014 
  !! Imported the noncolinear case implemented by xlzhang
  !!
  !
  USE ions_base,        ONLY : nat, ityp, ntyp => nsp
  USE io_files,         ONLY : prefix, tmp_dir
  USE io_global,        ONLY : stdout
  USE io_epw,           ONLY : iurecover
  USE pwcom,            ONLY : lspinorb, gg, ngm, tpiba2, nl, eigts3, &
                               eigts2, eigts1, omega, tpiba, g, nspin, mill, & 
                               domag
  USE scf,              ONLY : v, vltot
  USE noncollin_module, ONLY : noncolin
  USE kinds,            ONLY : DP
  USE phcom,            ONLY : int1, int2, int4, int5, recover,  &
                               int5_so, vlocq, int4_nc
  USE qpoint,           ONLY : xq, eigqts
  USE uspp_param,       ONLY : upf, lmaxq, nh
  USE mp_global,        ONLY : my_pool_id, npool, intra_pool_comm
  USE mp,               ONLY : mp_sum
  USE uspp,             ONLY : okvan
  USE fft_base,         ONLY : dfftp
  USE fft_interfaces,   ONLY : fwfft

  implicit none
  !
  !   And the local variables
  !
  INTEGER :: nt, na, nb, ig, nta, ntb, ir, ih, jh, ijh, ipol, jpol, is, nspin0
  !
  real(kind=DP), ALLOCATABLE :: qmod (:), qmodg (:), qpg (:,:), &
       ylmkq (:,:), ylmk0 (:,:)
  ! the modulus of q+G
  ! the modulus of G
  ! the  q+G vectors
  ! the spherical harmonics

  COMPLEX(kind=DP) :: fact, fact1, ZDOTC
  COMPLEX(kind=DP), ALLOCATABLE :: aux1 (:), aux2 (:),&
       aux3 (:), aux5 (:), veff (:,:), sk(:)
  ! work space
  complex(kind=DP), ALLOCATABLE, TARGET :: qgm(:)
  ! the augmentation function at G
  complex(kind=DP), POINTER :: qgmq (:)
  ! the augmentation function at q+G
  character (len=256) :: tempfile
  character (len=3) :: filelab
  logical :: exst

  IF (.not.okvan) RETURN

!  if (recover.and..not.ldisp) return

  nspin0=nspin
  IF (nspin==4.and..not.domag) nspin0=1

  CALL start_clock ('dvanqq2')
  int1(:,:,:,:,:) = (0.d0, 0.d0)
  int2(:,:,:,:,:) = (0.d0, 0.d0)
  int4(:,:,:,:,:) = (0.d0, 0.d0)
  int5(:,:,:,:,:) = (0.d0, 0.d0)
  ALLOCATE (sk  (  ngm))    
  ALLOCATE (aux1(  ngm))    
  ALLOCATE (aux2(  ngm))    
  ALLOCATE (aux3(  ngm))    
  ALLOCATE (aux5(  ngm))    
  ALLOCATE (qmodg( ngm))    
  ALLOCATE (ylmk0( ngm , lmaxq * lmaxq))    
  ALLOCATE (qgm  ( ngm))    
  ALLOCATE (ylmkq(ngm , lmaxq * lmaxq))    
  ALLOCATE (qmod( ngm))    
  ALLOCATE (qgmq( ngm))    
  !
  !     compute spherical harmonics
  !
  CALL ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
  DO ig = 1, ngm
     qmodg (ig) = sqrt (gg (ig) )
  ENDDO
  ALLOCATE (qpg (3, ngm))    
  CALL setqmod (ngm, xq, g, qmod, qpg)
  CALL ylmr2 (lmaxq * lmaxq, ngm, qpg, qmod, ylmkq)
  DEALLOCATE (qpg)
  DO ig = 1, ngm
     qmod (ig) = sqrt (qmod (ig) )
  ENDDO
  !   we start by computing the FT of the effective potential
  !
  ALLOCATE (veff ( dfftp%nnr , nspin))    
  DO is = 1, nspin
     IF (nspin.ne.4.or.is==1) THEN
        DO ir = 1, dfftp%nnr
           veff (ir, is) = CMPLX (vltot (ir) + v%of_r (ir, is), 0.d0, kind=DP)
        ENDDO
     ELSE
        DO ir = 1, dfftp%nnr
           veff (ir, is) = CMPLX (v%of_r (ir, is), 0.d0, kind=DP)
        ENDDO
     ENDIF
     CALL fwfft ('Dense', veff(:,is), dfftp)
  ENDDO
  !
  !     We compute here two of the three integrals needed in the phonon
  !
  fact1 = CMPLX (0.d0, - tpiba * omega, kind=DP)
  !
  tempfile = trim(tmp_dir) // trim(prefix) // '.recover' 
#if defined(__MPI)
  CALL set_ndnmbr (0,my_pool_id+1,1,npool,filelab)
  tempfile = trim(tmp_dir) // trim(prefix) // '.recover' // filelab
#endif
  !
  IF (recover) THEN
     WRITE (stdout,*) "    Using recover mode for USPP"
     inquire(file = tempfile, exist=exst)
     IF (.not. exst ) CALL errore( 'dvanqq2', 'recover files not found ', 1)
     OPEN (iurecover, file = tempfile, form = 'unformatted')
     IF (noncolin) THEN
        READ (iurecover) int1, int2, int4, int5
     ELSE
        READ (iurecover) int1, int2
     ENDIF
     CLOSE(iurecover)
  ELSE
  DO ntb = 1, ntyp
     IF (upf(ntb)%tvanp ) THEN
        ijh = 0
        DO ih = 1, nh (ntb)
           DO jh = ih, nh (ntb)
              ijh = ijh + 1
              !
              !    compute the augmentation function
              !
              CALL qvan2 (ngm, ih, jh, ntb, qmodg, qgm, ylmk0)
              CALL qvan2 (ngm, ih, jh, ntb, qmod, qgmq, ylmkq)
              !
              !     NB: for this integral the moving atom and the atom of Q
              !     do not necessarily coincide
              !
              DO nb = 1, nat
                 IF (ityp (nb) == ntb) THEN
                    DO ig = 1, ngm
                       aux1 (ig) = qgmq (ig) * eigts1 (mill(1,ig), nb) &
                                             * eigts2 (mill(2,ig), nb) &
                                             * eigts3 (mill(3,ig), nb)
                    ENDDO
                    DO na = 1, nat
                       fact = eigqts (na) * CONJG(eigqts (nb) )
                       !
                       !    nb is the atom of the augmentation function
                       !
                       nta = ityp (na)
                       DO ig=1, ngm
                          sk(ig)=vlocq(ig,nta) * eigts1(mill(1,ig), na) &
                                               * eigts2(mill(2,ig), na) &
                                               * eigts3(mill(3,ig), na) 
                       ENDDO
                       DO ipol = 1, 3
                          DO ig=1, ngm
                            aux5(ig)= sk(ig) * (g (ipol, ig) + xq (ipol) )
                          ENDDO
                          int2 (ih, jh, ipol, na, nb) = fact * fact1 * &
                                ZDOTC (ngm, aux1, 1, aux5, 1)
                          DO jpol = 1, 3
                             IF (jpol >= ipol) THEN
                                DO ig = 1, ngm
                                   aux3 (ig) = aux5 (ig) * &
                                               (g (jpol, ig) + xq (jpol) )
                                ENDDO
                                int5 (ijh, ipol, jpol, na, nb) = &
                                     CONJG(fact) * tpiba2 * omega * &
                                     ZDOTC (ngm, aux3, 1, aux1, 1)
                             ELSE
                                int5 (ijh, ipol, jpol, na, nb) = &
                                     int5 (ijh, jpol, ipol, na, nb)
                             ENDIF
                          ENDDO
                       ENDDO
                    ENDDO
                    DO ig = 1, ngm
                       aux1 (ig) = qgm (ig) * eigts1 (mill(1,ig), nb) &
                                            * eigts2 (mill(2,ig), nb) &
                                            * eigts3 (mill(3,ig), nb)
                    ENDDO
                    DO is = 1, nspin0
                       DO ipol = 1, 3
                          DO ig = 1, ngm
                             aux2 (ig) = veff (nl (ig), is) * g (ipol, ig)
                          ENDDO
                          int1 (ih, jh, ipol, nb, is) = - fact1 * &
                               ZDOTC (ngm, aux1, 1, aux2, 1)
                          DO jpol = 1, 3
                             IF (jpol >= ipol) THEN
                                DO ig = 1, ngm
                                   aux3 (ig) = aux2 (ig) * g (jpol, ig)
                                ENDDO
                                int4 (ijh, ipol, jpol, nb, is) = - tpiba2 * &
                                     omega * ZDOTC (ngm, aux3, 1, aux1, 1)
                             ELSE
                                int4 (ijh, ipol, jpol, nb, is) = &
                                     int4 (ijh, jpol, ipol, nb, is)
                             ENDIF
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
        DO ih = 1, nh (ntb)
           DO jh = ih + 1, nh (ntb)
              !
              !    We use the symmetry properties of the integral factor
              !
              DO nb = 1, nat
                 IF (ityp (nb) == ntb) THEN
                    DO ipol = 1, 3
                       DO is = 1, nspin
                          int1(jh,ih,ipol,nb,is) = int1(ih,jh,ipol,nb,is)
                       ENDDO
                       DO na = 1, nat
                          int2(jh,ih,ipol,na,nb) = int2(ih,jh,ipol,na,nb)
                       ENDDO
                    ENDDO
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDIF
  ENDDO
  CALL mp_sum(  int1, intra_pool_comm )
  CALL mp_sum(  int2, intra_pool_comm )
  call mp_sum(  int4, intra_pool_comm )
  call mp_sum(  int5, intra_pool_comm )
  OPEN (iurecover, file = tempfile, form = 'unformatted')
  IF (noncolin) THEN 
     WRITE (iurecover) int1, int2, int4, int5
  ELSE
     WRITE (iurecover) int1, int2
  ENDIF
  CLOSE(iurecover)
  ENDIF
  IF (noncolin) THEN
     CALL set_int12_nc(0)
     int4_nc = (0.d0, 0.d0)
     IF (lspinorb) int5_so = (0.d0, 0.d0)
     DO nt = 1, ntyp
        IF ( upf(nt)%tvanp ) THEN
           DO na = 1, nat
              IF (ityp(na)==nt) THEN
                 IF (upf(nt)%has_so)  THEN
                    CALL transform_int4_so(int4,na)
                    CALL transform_int5_so(int5,na)
                 ELSE
                    CALL transform_int4_nc(int4,na)
                    IF (lspinorb) CALL transform_int5_nc(int5,na)
                 ENDIF
              ENDIF
           ENDDO
        ENDIF
     ENDDO
  ENDIF
  !
  DEALLOCATE (veff)
  DEALLOCATE (qgmq)
  DEALLOCATE (qmod)
  DEALLOCATE (ylmkq)
  DEALLOCATE (qgm)
  DEALLOCATE (ylmk0)
  DEALLOCATE (qmodg)
  DEALLOCATE (aux5)
  DEALLOCATE (aux3)
  DEALLOCATE (aux2)
  DEALLOCATE (aux1)
  DEALLOCATE (sk)

  CALL stop_clock ('dvanqq2')
  CALL print_clock('dvanqq2')
  RETURN
END SUBROUTINE dvanqq2
