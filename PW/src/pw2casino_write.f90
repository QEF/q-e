!
! Copyright (C) 2004-2009 Dario Alfe' and Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE write_casino_wfn(gather,blip,multiplicity,binwrite,single_precision_blips,n_points_for_test,postfix)

   USE kinds, ONLY: DP,sgl
   USE ions_base, ONLY : nat, ntyp => nsp, ityp, tau, zv, atm
   USE cell_base, ONLY: omega, alat, tpiba2, at, bg
   USE run_info,  ONLY: title    ! title of the run
   USE constants, ONLY: tpi, e2, eps6
   USE ener, ONLY: ewld, ehart, etxc, vtxc, etot, etxcc, demet, ef
   USE fft_base,  ONLY: dfftp
   USE fft_interfaces, ONLY : fwfft
   USE gvect, ONLY: ngm, gstart, g, gg, gcutm, nl, nlm, igtongl
   USE klist , ONLY: nks, nelec, xk, wk, degauss, ngauss, igk_k, ngk
   USE lsda_mod, ONLY: lsda, nspin
   USE scf, ONLY: rho, rho_core, rhog_core, v
   USE ldaU, ONLY : eth
   USE vlocal, ONLY: vloc, strf
   USE wvfct, ONLY: npwx, nbnd, wg, et
   USE gvecw, ONLY: ecutwfc
   USE control_flags, ONLY : gamma_only
   USE uspp, ONLY: nkb, vkb, dvan
   USE uspp_param, ONLY: nh
   USE io_global, ONLY: stdout, ionode, ionode_id
   USE io_files, ONLY: nd_nmbr, nwordwfc, iunwfc, prefix, tmp_dir, seqopn
   USE wavefunctions_module, ONLY : evc
   USE funct, ONLY : dft_is_meta
   USE mp_pools, ONLY: inter_pool_comm, intra_pool_comm, nproc_pool, me_pool
   USE mp_bands, ONLY: intra_bgrp_comm
   USE mp, ONLY: mp_sum, mp_gather, mp_bcast, mp_get
   USE buffers,              ONLY : get_buffer

   USE pw2blip

   IMPLICIT NONE
   LOGICAL, INTENT(in) :: gather,blip,binwrite,single_precision_blips
   REAL(dp), INTENT(in) :: multiplicity
   INTEGER, INTENT(in) :: n_points_for_test
   CHARACTER(*), INTENT(in) :: postfix

   INTEGER, PARAMETER :: n_overlap_tests = 12
   REAL(dp), PARAMETER :: eps = 1.d-10
   INTEGER, PARAMETER :: io = 77, iob = 78
   INTEGER :: npw, ig, ibnd, ik, ispin, nbndup, nbnddown, &
              nk, ig7, ikk, id, ip, iorb, iorb_node, inode, ierr, norb
   INTEGER :: jk(nproc_pool), jspin(nproc_pool), jbnd(nproc_pool)
   INTEGER :: jk2(nproc_pool), jspin2(nproc_pool), jbnd2(nproc_pool)
   INTEGER, ALLOCATABLE :: idx(:), igtog(:), gtoig(:)
   LOGICAL :: exst,dowrite
   REAL(DP) :: ek, eloc, enl, etot_
   INTEGER, EXTERNAL :: atomic_number
   REAL (DP), EXTERNAL :: ewald, w1gauss

   ! number of g vectors (union of all k points)
   INTEGER ngtot_l ! on this processor
   INTEGER, ALLOCATABLE :: ngtot_d(:), ngtot_cumsum(:), indx(:)
   INTEGER ngtot_g ! sum over processors
   REAL(DP), ALLOCATABLE :: g_l(:,:), g_g(:,:), g2(:)
   COMPLEX(DP), ALLOCATABLE :: evc_l(:), evc_g(:), evc_g2(:), avc_tmp(:,:,:), cavc_tmp(:,:,:)
   LOGICAL dotransform

   REAL(dp) :: av_overlap(5,2),avsq_overlap(5,2)

!----------------------------------------------------------------------------!
! Random number generator, using the method suggested by D.E. Knuth in       !
! Seminumerical Algorithms (vol 2 of The Art of Computer Programming).       !
! The method is based on lagged Fibonacci sequences with subtraction.        !
!----------------------------------------------------------------------------!
   INTEGER,PARAMETER :: KK=100,LL=37 ! Leave these.
   REAL(DP) :: ranstate(kk)  ! Determines output of gen_ran_array.

   INTEGER,PARAMETER :: default_seed=310952  ! Random seed, betw. 0 & 2^30-3.
   INTEGER,PARAMETER :: Nran=1009,Nkeep=100 ! See comment on p. 188 of Knuth.
   INTEGER,SAVE :: ran_array_idx=-1
   REAL(DP),SAVE :: ran_array(Nran)

   dowrite=ionode.or..not.(gather.or.blip)

   ALLOCATE (idx (ngm) )
   ALLOCATE (igtog (ngm) )
   ALLOCATE (gtoig (ngm) )
   idx(:) = 0
   igtog(:) = 0

   IF( lsda )THEN
      nbndup = nbnd
      nbnddown = nbnd
      nk = nks/2
      !     nspin = 2
   ELSE
      nbndup = nbnd
      nbnddown = 0
      nk = nks
      !     nspin = 1
   ENDIF

   CALL calc_energies

   DO ispin = 1, nspin
      DO ik = 1, nk
         ikk = ik + nk*(ispin-1)
         idx( igk_k(1:ngk(ikk),ikk) ) = 1
      ENDDO
   ENDDO

   ngtot_l = 0
   DO ig = 1, ngm
      IF( idx(ig) >= 1 )THEN
         ngtot_l = ngtot_l + 1
         igtog(ngtot_l) = ig
         gtoig(ig) = ngtot_l
      ENDIF
   ENDDO

   DEALLOCATE (idx)

   IF(dowrite)THEN
      IF(blip)THEN
         IF(binwrite)THEN
            WRITE (6,'(a)')'Writing file '//trim(prefix)//'.bwfn.data.b1'//trim(postfix)//' for program CASINO.'
            OPEN( iob, file=trim(tmp_dir)//trim(prefix)//'.bwfn.data.b1'//trim(postfix), &
                  form='unformatted', action='write', access='sequential')
         ELSE
            WRITE (6,'(a)')'Writing file '//trim(prefix)//'.bwfn.data'//trim(postfix)//' for program CASINO.'
            OPEN( io, file=trim(tmp_dir)//trim(prefix)//'.bwfn.data'//trim(postfix), &
                  form='formatted', action='write', access='sequential')
         ENDIF
      ELSE
         IF(gather)THEN
            WRITE (6,'(a)')'Writing file '//trim(prefix)//'.pwfn.data'//trim(postfix)//' for program CASINO.'
            OPEN( io, file=trim(tmp_dir)//trim(prefix)//'.pwfn.data'//trim(postfix), & 
                  form='formatted', action='write', access='sequential')
         ELSE
            WRITE (6,'(a)')'Writing one file per node '//trim(prefix)//'.pwfn.data'//trim(postfix)//'.XX for program CASINO'
            CALL seqopn( io, 'pwfn.data'//trim(postfix), 'formatted',exst)
         ENDIF
      ENDIF
      WRITE (6,'(a)')
   ENDIF

   ALLOCATE ( g_l(3,ngtot_l), evc_l(ngtot_l) )
   DO ig = 1, ngtot_l
      g_l(:,ig) = g(:,igtog(ig))
   ENDDO

   IF(gather.or.blip)THEN
      ALLOCATE ( ngtot_d(nproc_pool), ngtot_cumsum(nproc_pool) )
      CALL mp_gather( ngtot_l, ngtot_d, ionode_id, intra_pool_comm )
      CALL mp_bcast( ngtot_d, ionode_id, intra_pool_comm )
      id = 0
      DO ip = 1,nproc_pool
         ngtot_cumsum(ip) = id
         id = id + ngtot_d(ip)
      ENDDO
      ngtot_g = id

      ALLOCATE ( g_g(3,ngtot_g), g2(ngtot_g), evc_g(ngtot_g) )
      IF(blip.and.gamma_only)THEN
         ALLOCATE( evc_g2(ngtot_g) )
      ENDIF
      CALL mp_gather( g_l, g_g, ngtot_d, ngtot_cumsum, ionode_id, intra_pool_comm)

      IF(blip)THEN
         CALL mp_bcast( g_g, ionode_id, intra_pool_comm )
         g2(:) = sum(g_g(:,:)**2,dim=1)
         CALL pw2blip_init(ngtot_g,g_g,multiplicity)
         IF(dowrite)THEN
            WRITE (6,'(a)')'Blip grid: '//trim(i2s(blipgrid(1)))//'x'//trim(i2s(blipgrid(2)))//'x'//trim(i2s(blipgrid(3)))
            WRITE (6,'(a)')
         ENDIF
      ELSEIF(dowrite)THEN
         ALLOCATE ( indx(ngtot_g) )
         CALL create_index2(g_g,indx)
      ENDIF
   ELSEIF(dowrite)THEN
      ALLOCATE ( indx(ngtot_l) )
      CALL create_index2(g_l,indx)
   ENDIF

   IF(dowrite)THEN
      CALL write_header
      IF(blip)THEN
         CALL write_gvecs_blip
      ELSEIF(gather)THEN
         CALL write_gvecs(g_g,indx)
      ELSE
         CALL write_gvecs(g_l,indx)
      ENDIF
      CALL write_wfn_head
   ENDIF

   IF(dowrite.and.blip.and.binwrite)THEN
      IF(gamma_only)THEN
         ALLOCATE(avc_tmp(blipgrid(1),blipgrid(2),blipgrid(3)))
      ELSE
         ALLOCATE(cavc_tmp(blipgrid(1),blipgrid(2),blipgrid(3)))
      ENDIF
   ENDIF

   ! making some assumptions about the parallel layout:
   IF(ionode_id/=0)CALL errore('write_casino_wfn','ionode_id/=0: ',ionode_id)

   iorb = 0
   norb = nk*nspin*nbnd

   DO ik = 1, nk
      DO ispin = 1, nspin
         ikk = ik + nk*(ispin-1)
         npw = ngk(ikk)
         IF( nks > 1 ) CALL get_buffer(evc,nwordwfc,iunwfc,ikk)
         DO ibnd = 1, nbnd
            evc_l(:) = (0.d0, 0d0)
            evc_l(gtoig(igk_k(1:npw,ikk))) = evc(1:npw,ibnd)
            IF(blip)THEN
               iorb = iorb + 1
               IF(gamma_only)THEN
                  iorb_node = mod((iorb-1)/2,nproc_pool) ! the node that should compute this orbital
                  IF(mod(iorb,2)==0)THEN
                     jk2(iorb_node+1) = ik
                     jspin2(iorb_node+1) = ispin
                     jbnd2(iorb_node+1) = ibnd
                     dotransform = (iorb_node==nproc_pool-1)
                  ELSE
                     jk(iorb_node+1) = ik
                     jspin(iorb_node+1) = ispin
                     jbnd(iorb_node+1) = ibnd
                     dotransform = .false.
                  ENDIF
               ELSE
                  iorb_node = mod(iorb-1,nproc_pool) ! the node that should compute this orbital
                  jk(iorb_node+1) = ik
                  jspin(iorb_node+1) = ispin
                  jbnd(iorb_node+1) = ibnd
                  dotransform=(iorb_node==nproc_pool-1)
               ENDIF
               DO inode=0,nproc_pool-1
                  IF(gamma_only.and.mod(iorb,2)==0)THEN
                     CALL mp_get(&
                        evc_g2(ngtot_cumsum(inode+1)+1:ngtot_cumsum(inode+1)+ngtot_d(inode+1)),&
                        evc_l(:),me_pool,iorb_node,inode,1234,intra_pool_comm)
                  ELSE
                     CALL mp_get(&
                        evc_g(ngtot_cumsum(inode+1)+1:ngtot_cumsum(inode+1)+ngtot_d(inode+1)),&
                        evc_l(:),me_pool,iorb_node,inode,1234,intra_pool_comm)
                  ENDIF
               ENDDO
               IF(dotransform .or. iorb == norb)THEN
                  IF(me_pool <= iorb_node)THEN
                     IF(gamma_only.and.(me_pool/=iorb_node.or.iorb/=norb.or.mod(norb,2)==0))THEN
                        CALL pw2blip_transform2(evc_g(:),evc_g2(:))
                     ELSE
                        CALL pw2blip_transform(evc_g(:))
                     ENDIF
                  ENDIF
                  IF(me_pool <= iorb_node) CALL test_overlap
                  DO inode=0,iorb_node
                     CALL pw2blip_get(inode)
                     IF(gamma_only)THEN
                        IF(ionode)WRITE(6,*)"Transformed real orbital k="//trim(i2s(jk(inode+1)))//&
                           &", spin="//trim(i2s(jspin(inode+1)))//&
                           &", band="//trim(i2s(jbnd(inode+1)))//" on node "//trim(i2s(inode))
                        CALL print_overlap(inode,1)
                        IF(blipreal==2)THEN
                           IF(ionode)WRITE(6,*)"Transformed real orbital k="//trim(i2s(jk2(inode+1)))//&
                              &", spin="//trim(i2s(jspin2(inode+1)))//&
                              &", band="//trim(i2s(jbnd2(inode+1)))//" on node "//trim(i2s(inode))
                        ENDIF
                        CALL print_overlap(inode,2)
                     ELSE
                        IF(ionode)WRITE(6,*)"Transformed complex orbital k="//trim(i2s(jk(inode+1)))//&
                           &", spin="//trim(i2s(jspin(inode+1)))//&
                           &", band="//trim(i2s(jbnd(inode+1)))//" on node "//trim(i2s(inode))
                        CALL print_overlap(inode,1)
                     ENDIF
                     IF(gamma_only)THEN
                        IF(ionode)CALL write_bwfn_data_gamma(1,jk(inode+1),jspin(inode+1),jbnd(inode+1))
                        IF(blipreal==2)THEN
                           IF(ionode)CALL write_bwfn_data_gamma(2,jk2(inode+1),jspin2(inode+1),jbnd2(inode+1))
                        ENDIF
                     ELSE
                        IF(ionode)CALL write_bwfn_data(jk(inode+1),jspin(inode+1),jbnd(inode+1))
                     ENDIF
                  ENDDO
               ENDIF

            ELSEIF(gather)THEN
               CALL mp_gather( evc_l, evc_g, ngtot_d, ngtot_cumsum, ionode_id, intra_pool_comm)
               IF(dowrite)CALL write_pwfn_data(ik,ispin,ibnd,evc_g,indx)
            ELSE
               CALL write_pwfn_data(ik,ispin,ibnd,evc_l,indx)
            ENDIF
         ENDDO
      ENDDO
   ENDDO
   IF(dowrite)THEN
      IF(binwrite)THEN
         CLOSE(iob)
      ELSE
         CLOSE(io)
      ENDIF
   ENDIF
   IF(dowrite.and.blip.and.binwrite)THEN
      IF(gamma_only)THEN
         DEALLOCATE(avc_tmp)
      ELSE
         DEALLOCATE(cavc_tmp)
      ENDIF
   ENDIF

   IF(blip)CALL pw2blip_cleanup
   DEALLOCATE (igtog, g_l, evc_l )
   IF(blip.or.gather) DEALLOCATE ( ngtot_d, ngtot_cumsum, g_g, evc_g )
   IF(dowrite.and..not.blip) DEALLOCATE (indx)

CONTAINS

   SUBROUTINE calc_energies
      USE becmod, ONLY: becp, calbec, allocate_bec_type, deallocate_bec_type
      USE exx,    ONLY : exxenergy2, fock2
      USE funct,  ONLY : dft_is_hybrid

      COMPLEX(DP), ALLOCATABLE :: aux(:)
      INTEGER :: npw, ibnd, j, ig, ik,ikk, ispin, na, nt, ijkb0, ikb,jkb, ih,jh
      REAL(dp), ALLOCATABLE :: g2kin(:)
      REAL(DP) :: charge, etotefield, elocg

      ALLOCATE (aux(dfftp%nnr))
      CALL allocate_bec_type ( nkb, nbnd, becp )

      ek  = 0.d0
      eloc= 0.d0
      enl = 0.d0
      demet=0.d0
      fock2=0.d0
      !
      ALLOCATE ( g2kin(npwx) )
      DO ispin = 1, nspin
         !
         !     calculate the local contribution to the total energy
         !
         !      bring rho to G-space
         !
         aux(:) = cmplx( rho%of_r(:,ispin), 0.d0,kind=DP)
         CALL fwfft ('Dense', aux, dfftp)
         !
         DO nt=1,ntyp
            DO ig = 1, ngm
               elocg = vloc(igtongl(ig),nt) * &
                       dble ( strf(ig,nt) * conjg(aux(nl(ig))) )
               eloc = eloc + elocg
               IF( gamma_only .and. ig>=gstart) eloc = eloc + elocg
            ENDDO
         ENDDO

         DO ik = 1, nk
            ikk = ik + nk*(ispin-1)
            npw = ngk(ikk)
            IF( nks > 1 ) CALL get_buffer (evc, nwordwfc, iunwfc, ikk )
            CALL init_us_2 (npw, igk_k(1,ikk), xk (1, ikk), vkb)
            CALL calbec ( npw, vkb, evc, becp )
            !
            ! -TS term for metals (if any)
            !
            IF( degauss > 0.0_dp)THEN
               DO ibnd = 1, nbnd
                  demet = demet + wk (ik) * &
                     degauss * w1gauss ( (ef-et(ibnd,ik)) / degauss, ngauss)
               ENDDO
            ENDIF
            !
            ! calculate the kinetic energy
            !
            g2kin(1:npw) = ( ( xk(1,ikk) + g(1,igk_k(1:npw,ikk)) )**2 + &
                             ( xk(2,ikk) + g(2,igk_k(1:npw,ikk)) )**2 + &
                             ( xk(3,ikk) + g(3,igk_k(1:npw,ikk)) )**2 ) * tpiba2
            DO ibnd = 1, nbnd
               DO j = 1, npw
                  IF(gamma_only)THEN !.and.j>1)then
                     ek = ek +  2*conjg(evc(j,ibnd)) * evc(j,ibnd) * &
                                    g2kin(j) * wg(ibnd,ikk)
                  ELSE
                     ek = ek +  conjg(evc(j,ibnd)) * evc(j,ibnd) * &
                                    g2kin(j) * wg(ibnd,ikk)
                  ENDIF
               ENDDO

               !
               ! Calculate Non-local energy
               !
               ijkb0 = 0
               DO nt = 1, ntyp
                  DO na = 1, nat
                     IF(ityp (na) == nt)THEN
                        DO ih = 1, nh (nt)
                           ikb = ijkb0 + ih
                           IF(gamma_only)THEN
                              enl=enl+becp%r(ikb,ibnd)*becp%r(ikb,ibnd) &
                                 *wg(ibnd,ikk)* dvan(ih,ih,nt)
                           ELSE
                              enl=enl+conjg(becp%k(ikb,ibnd))*becp%k(ikb,ibnd) &
                                 *wg(ibnd,ikk)* dvan(ih,ih,nt)
                           ENDIF
                           DO jh = ( ih + 1 ), nh(nt)
                              jkb = ijkb0 + jh
                              IF(gamma_only)THEN
                                 enl=enl + &
                                    (becp%r(ikb,ibnd)*becp%r(jkb,ibnd)+&
                                       becp%r(jkb,ibnd)*becp%r(ikb,ibnd))&
                                    * wg(ibnd,ikk) * dvan(ih,jh,nt)
                              ELSE
                                 enl=enl + &
                                    (conjg(becp%k(ikb,ibnd))*becp%k(jkb,ibnd)+&
                                       conjg(becp%k(jkb,ibnd))*becp%k(ikb,ibnd))&
                                    * wg(ibnd,ikk) * dvan(ih,jh,nt)
                              ENDIF

                           ENDDO

                        ENDDO
                        ijkb0 = ijkb0 + nh (nt)
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      DEALLOCATE ( g2kin )
#if defined(__MPI)
      CALL mp_sum( eloc,  intra_bgrp_comm )
      CALL mp_sum( ek,    intra_bgrp_comm )
      CALL mp_sum( ek,    inter_pool_comm )
      CALL mp_sum( enl,   inter_pool_comm )
      CALL mp_sum( demet, inter_pool_comm )
#endif
      eloc = eloc * omega
      !
      ! compute ewald contribution
      !
      ewld = ewald( alat, nat, ntyp, ityp, zv, at, bg, tau, omega, &
           g, gg, ngm, gcutm, gstart, gamma_only, strf )
      !
      ! compute hartree and xc contribution
      !
      CALL v_of_rho( rho, rho_core, rhog_core, &
                     ehart, etxc, vtxc, eth, etotefield, charge, v )
      !
      ! compute exact exchange contribution (if present)
      !
      IF(dft_is_hybrid()) fock2 = 0.5_DP * exxenergy2()
      !
      etot_=(ek + (etxc-etxcc)+ehart+eloc+enl+ewld)+demet+fock2
      !
      IF ( ABS(etot-etot_) > ABS(eps6*etot) ) THEN
         WRITE (stdout,'(5X,"Etot: ",f15.8," Ry from PWscf vs ", &
                & f15.8," Ry from pw2casino!")') etot, etot_
         CALL errore("pw2casino","Mismatch in computed energy",1)
      ELSE
         etot = etot_
      END IF
      !
      CALL deallocate_bec_type (becp)
      DEALLOCATE (aux)

      WRITE (stdout,*)
      WRITE (stdout,*) 'Energies determined by pw2casino tool'
      WRITE (stdout,*) '-------------------------------------'
      WRITE (stdout,*) 'Kinetic energy   ', ek/e2, ' au  =  ', ek, ' Ry'
      WRITE (stdout,*) 'Local energy     ', eloc/e2, ' au  =  ', eloc, ' Ry'
      WRITE (stdout,*) 'Non-Local energy ', enl/e2, ' au  =  ', enl, ' Ry'
      WRITE (stdout,*) 'Ewald energy     ', ewld/e2, ' au  =  ', ewld, ' Ry'
      WRITE (stdout,*) 'xc contribution  ',(etxc-etxcc)/e2, ' au  =  ', etxc-etxcc, ' Ry'
      WRITE (stdout,*) 'hartree energy   ', ehart/e2, ' au  =  ', ehart, ' Ry'
      IF(dft_is_hybrid()) & 
           WRITE (stdout,*) 'EXX energy       ', fock2/e2, ' au  =  ', fock2, ' Ry' 
      IF( degauss > 0.0_dp ) &
         WRITE (stdout,*) 'Smearing (-TS)   ', demet/e2, ' au  =  ', demet, ' Ry'
      WRITE (stdout,*) 'Total energy     ', etot/e2, ' au  =  ', etot, ' Ry'
      WRITE (stdout,*)

   END SUBROUTINE calc_energies


   SUBROUTINE test_overlap
! Carry out the overlap test described in the CASINO manual.
! Repeat the whole test n_overlap_tests times, to compute error bars.
      INTEGER i,j,k
      REAL(dp) r(3)
      COMPLEX(dp) xb(5),xp(5) ! 1->val, 2:4->grad, 5->lap
      REAL(dp) xbb(5,2),xpp(5,2)
      COMPLEX(dp) xbp(5,2)
      REAL(dp) overlap(5,2),sum_overlap(5,2),sumsq_overlap(5,2)

      IF(n_points_for_test<=0)RETURN
      IF(n_overlap_tests<=0)RETURN

      CALL init_rng(12345678)

      sum_overlap(:,:)=0.d0 ; sumsq_overlap(:,:)=0.d0
      DO j=1,n_overlap_tests
         xbb(:,:)=0.d0 ; xpp(:,:)=0.d0 ; xbp(:,:)=0.d0

         DO i=1,n_points_for_test
            r(1)=ranx() ; r(2)=ranx() ; r(3)=ranx()
            CALL blipeval(r,xb(1),xb(2:4),xb(5))
            CALL pweval(r,xp(1),xp(2:4),xp(5))

            IF(gamma_only)THEN
               xbb(:,1)=xbb(:,1)+dble(xb(:))**2
               xbp(:,1)=xbp(:,1)+dble(xb(:))*dble(xp(:))
               xpp(:,1)=xpp(:,1)+dble(xp(:))**2
               IF(blipreal==2)THEN
                  ! two orbitals - use complex and imaginary part independently
                  xbb(:,2)=xbb(:,2)+aimag(xb(:))**2
                  xbp(:,2)=xbp(:,2)+aimag(xb(:))*aimag(xp(:))
                  xpp(:,2)=xpp(:,2)+aimag(xp(:))**2
               ENDIF
            ELSE
               xbb(:,1)=xbb(:,1)+dble(xb(:))**2+aimag(xb(:))**2
               xbp(:,1)=xbp(:,1)+xb(:)*conjg(xp(:))
               xpp(:,1)=xpp(:,1)+dble(xp(:))**2+aimag(xp(:))**2
            ENDIF
         ENDDO ! i
         overlap(:,:)=0.d0
         DO k=1,5
            IF(xbb(k,1)/=0.d0.and.xpp(k,1)/=0.d0)THEN
               overlap(k,1)=(dble(xbp(k,1))**2+aimag(xbp(k,1))**2)/(xbb(k,1)*xpp(k,1))
            ENDIF ! xb & xd nonzero
         ENDDO ! k

         IF(blipreal==2)THEN
            DO k=1,5
               IF(xbb(k,2)/=0.d0.and.xpp(k,2)/=0.d0)THEN
                  overlap(k,2)=(dble(xbp(k,2))**2+aimag(xbp(k,2))**2)/(xbb(k,2)*xpp(k,2))
               ENDIF ! xb & xd nonzero
            ENDDO ! k
         ELSE
         ENDIF
         sum_overlap(:,:)=sum_overlap(:,:)+overlap(:,:)
         sumsq_overlap(:,:)=sumsq_overlap(:,:)+overlap(:,:)**2
      ENDDO ! j
      av_overlap(:,:)=sum_overlap(:,:)/dble(n_overlap_tests)
      avsq_overlap(:,:)=sumsq_overlap(:,:)/dble(n_overlap_tests)

   END SUBROUTINE test_overlap


   SUBROUTINE pweval(r,val,grad,lap)
      DOUBLE PRECISION,INTENT(in) :: r(3)
      COMPLEX(dp),INTENT(out) :: val,grad(3),lap

      INTEGER ig
      REAL(dp) dot_prod
      COMPLEX(dp) eigr,eigr2

      REAL(dp),PARAMETER :: pi=3.141592653589793238462643d0
      COMPLEX(dp),PARAMETER :: iunity=(0.d0,1.d0)

      val=0.d0 ; grad(:)=0.d0 ; lap=0.d0
      DO ig=1,ngtot_g
         dot_prod=tpi*sum(dble(g_int(:,ig))*r(:))
         eigr=evc_g(ig)*cmplx(cos(dot_prod),sin(dot_prod),KIND=dp)
         IF(.not.gamma_only)THEN
            val=val+eigr
            grad(:)=grad(:)+(eigr*iunity)*dble(g_int(:,ig))
            lap=lap-eigr*g2(ig)
         ELSEIF(blipreal==1)THEN
            IF(all(g_int(:,ig)==0))eigr=eigr*0.5d0
            val=val+dble(eigr)
            grad(:)=grad(:)-aimag(eigr)*dble(g_int(:,ig))
            lap=lap-dble(eigr)*g2(ig)
         ELSEIF(blipreal==2)THEN
            eigr2=evc_g2(ig)*cmplx(cos(dot_prod),sin(dot_prod),KIND=dp)
            IF(all(g_int(:,ig)==0))THEN
               eigr=eigr*0.5d0
               eigr2=eigr2*0.5d0
            ENDIF
            val=val+cmplx(dble(eigr),dble(eigr2),KIND=dp)
            grad(:)=grad(:)+cmplx(-aimag(eigr),-aimag(eigr2),KIND=dp)*dble(g_int(:,ig))
            lap=lap-cmplx(dble(eigr),dble(eigr2),KIND=dp)*g2(ig)
         ENDIF
      ENDDO ! ig
      IF(gamma_only)THEN
         val = val*2.d0
         grad(:) = grad(:)*2.d0
         lap = lap*2.d0
      ENDIF
      grad(:)=matmul(bg(:,:),grad(:))*(tpi/alat)
      lap=lap*(tpi/alat)**2
   END SUBROUTINE pweval


   SUBROUTINE print_overlap(inode,whichband)
!-------------------------------------------------------------------------!
! Write out the overlaps of the value, gradient and Laplacian of the blip !
! orbitals.  Give error bars where possible.                              !
!-------------------------------------------------------------------------!
      INTEGER,INTENT(in) :: inode
      INTEGER,INTENT(in) :: whichband ! 1 or 2, indexing within a pair of real orbitals
      REAL(dp) :: av(5),avsq(5),err(5)
      INTEGER k
      CHARACTER(12) char12_arr(5)

      IF(n_points_for_test<=0)RETURN
      IF(n_overlap_tests<=0)RETURN

      CALL mp_get(av(:),av_overlap(:,whichband),me_pool,ionode_id,inode,6434,intra_pool_comm)
      CALL mp_get(avsq(:),avsq_overlap(:,whichband),me_pool,ionode_id,inode,6434,intra_pool_comm)

      IF(.not.ionode)RETURN
      IF(blipreal==1.and.whichband==2)RETURN

      IF(n_overlap_tests<2)THEN
         WRITE(stdout,*)'Error: need at least two overlap tests, to estimate error bars.'
         STOP
      ENDIF ! Too few overlap tests
      err(:)=sqrt(max(avsq(:)-av(:)**2,0.d0)/dble(n_overlap_tests-1))
      DO k=1,5
         char12_arr(k)=trim(write_mean(av(k),err(k)))
      ! Not room to quote error bar.  Just quote mean.
         IF(index(char12_arr(k),')')==0)WRITE(char12_arr(k),'(f12.9)')av(k)
      ENDDO ! k
      WRITE(stdout,'(2(1x,a),2x,3(1x,a))')char12_arr(1:5)
   END SUBROUTINE print_overlap


   FUNCTION to_c80(c)
      CHARACTER(*),INTENT(in) :: c
      CHARACTER(80) :: to_c80
      to_c80=c
   END FUNCTION to_c80

   SUBROUTINE write_header
      INTEGER j, na, nt, at_num
      REAL(dp) :: kvec(3,nk),ksq(nk),kprod(6,nk)

      IF(blip.and.binwrite)THEN
         WRITE(iob)&
            to_c80(title)    ,&
            to_c80("PWSCF")  ,&
            to_c80("DFT")    ,&
            to_c80("unknown"),&
            to_c80("unknown"),&
            dble(ecutwfc/2)  ,&
            lsda             ,&
            dble(etot/e2)    ,&
            dble(ek/e2)      ,&
            dble(eloc/e2)    ,&
            dble(enl/e2)     ,&
            dble(ehart/e2)   ,&
            dble(ewld/e2)    ,&
            nint(nelec)      ,&
            nat              ,&
            ngtot_g          ,&
            nk               ,&
            blipgrid(1:3)    ,&
            nbnd             ,&
            gamma_only       ,&
            .true.           ,&
            (/0,0/)          ,&
            alat*at(1:3,1)   ,&
            alat*at(1:3,2)   ,&
            alat*at(1:3,3)   ,&
            2                ,&
            nbnd

!     some old PGI compiler seems to choke on this commented version....
!             to_c80(title)    ,& ! title
!             to_c80("PWSCF")  ,& ! code
!             to_c80("DFT")    ,& ! method
!             to_c80("unknown"),& ! functional
!             to_c80("unknown"),& ! pseudo_type
!             dble(ecutwfc/2)  ,&  ! plane_wave_cutoff
!             lsda             ,&  ! spin_polarized,
!             dble(etot/e2)    ,&  ! total_energy
!             dble(ek/e2)      ,&  ! kinetic_energy
!             dble(eloc/e2)    ,&  ! local_potential_energy
!             dble(enl/e2)     ,&  ! non_local_potential_energy
!             dble(ehart/e2)   ,&  ! electron_electron_energy
!             dble(ewld/e2)    ,&  ! eionion
!             nint(nelec)      ,&  ! num_electrons
!             nat              ,&  ! nbasis
!             ngtot_g          ,&  ! nwvec
!             nk               ,&  ! nkvec
!             blipgrid(1:3)    ,&  ! nr
!             nbnd             ,&  ! maxband
!             gamma_only       ,&  ! gamma_only
!             .true.           ,&  ! ext_orbs_present
!             (/0,0/)          ,&  ! no_loc_orbs
!             alat*at(1:3,1)   ,&  ! pa1
!             alat*at(1:3,2)   ,&  ! pa2
!             alat*at(1:3,3)   ,&  ! pa3
!             2                ,&  ! nspin_check
!             nbnd                 ! num_nonloc_max

         kvec(:,:) = tpi/alat*xk(1:3,1:nk)
         kprod(1,:)=kvec(1,:)*kvec(1,:)
         kprod(2,:)=kvec(2,:)*kvec(2,:)
         kprod(3,:)=kvec(3,:)*kvec(3,:)
         kprod(4,:)=kvec(1,:)*kvec(2,:)
         kprod(5,:)=kvec(1,:)*kvec(3,:)
         kprod(6,:)=kvec(2,:)*kvec(3,:)
         ksq(:)=kprod(1,:)+kprod(2,:)+kprod(3,:)

         WRITE(iob)&
            kvec                                          ,&
            ksq                                           ,&
            kprod                                         ,&
            (atomic_number(trim(atm(ityp(na)))),na=1,nat) ,&
            (alat*tau(1:3,na),na=1,nat)                   ,&
            (nbnd,j=1,nk*2)                               ,&
            et(1:nbnd,1:nk*nspin)/e2                      ,&
            (.true.,j=1,nbnd*nk*nspin)                    ,&
            (/nbnd,nbnd/)

!             kvec                                          ,& ! kvec
!             ksq                                           ,& ! ksq
!             kprod                                         ,& ! kprod
!             (atomic_number(trim(atm(ityp(na)))),na=1,nat) ,& ! atno   -- atomic numbers
!             (alat*tau(1:3,na),na=1,nat)                   ,& ! basis  -- atom positions
!             (nbnd,j=1,nk*2)                               ,& ! nband
!             et(1:nbnd,1:nk*nspin)/e2                      ,& ! eigenvalue
!             (.true.,j=1,nbnd*nk*nspin)                    ,& ! on_this_cpu
!             (/nbnd,nbnd/)                                    ! num_nonloc

         WRITE(iob)single_precision_blips                    ! single_precision_blips

         ! IF(no_loc_orbs>0)THEN
         !    ...
         ! ENDIF

         WRITE(iob)&
          (0,j=1,nbnd*nk*2) ,&
          (0,j=1,nbnd*nk*2) ,&
          (0,j=1,nbnd*nk*2) ,&
          (0,j=1,nbnd*nk*2)

!           (0,j=1,nbnd*nk*2) ,& ! orb_map_band
!           (0,j=1,nbnd*nk*2) ,& ! orb_map_ik
!           (0,j=1,nbnd*nk*2) ,& ! orb_map_iorb
!           (0,j=1,nbnd*nk*2)    ! occupied
         RETURN
      ENDIF

      WRITE(io,'(a)') title
      WRITE(io,'(a)')
      WRITE(io,'(a)') ' BASIC INFO'
      WRITE(io,'(a)') ' ----------'
      WRITE(io,'(a)') ' Generated by:'
      WRITE(io,'(a)') '  PWSCF'
      WRITE(io,'(a)') ' Method:'
      WRITE(io,'(a)') '  DFT'
      WRITE(io,'(a)') ' DFT Functional:'
      WRITE(io,'(a)') '  unknown'
      WRITE(io,'(a)') ' Pseudopotential'
      WRITE(io,'(a)') '  unknown'
      WRITE(io,'(a)') ' Plane wave cutoff (au)'
      WRITE(io,*) ecutwfc/2
      WRITE(io,'(a)') ' Spin polarized:'
      WRITE(io,*)lsda
      IF( degauss > 0.0_dp )THEN
         WRITE(io,'(a)') ' Total energy (au per primitive cell; includes -TS term)'
         WRITE(io,*)etot/e2, demet/e2
      ELSE
         WRITE(io,'(a)') ' Total energy (au per primitive cell)'
         WRITE(io,*)etot/e2
      ENDIF
      WRITE(io,'(a)') ' Kinetic energy (au per primitive cell)'
      WRITE(io,*)ek/e2
      WRITE(io,'(a)') ' Local potential energy (au per primitive cell)'
      WRITE(io,*)eloc/e2
      WRITE(io,'(a)') ' Non local potential energy(au per primitive cell)'
      WRITE(io,*)enl/e2
      WRITE(io,'(a)') ' Electron electron energy (au per primitive cell)'
      WRITE(io,*)ehart/e2
      WRITE(io,'(a)') ' Ion-ion energy (au per primitive cell)'
      WRITE(io,*)ewld/e2
      WRITE(io,'(a)') ' Number of electrons per primitive cell'
      WRITE(io,*)nint(nelec)
      ! uncomment the following ifyou want the Fermi energy - KN 2/4/09
      !  WRITE(io,'(a)') ' Fermi energy (au)'
      !  WRITE(io,*) ef/e2
      WRITE(io,'(a)') ' '
      WRITE(io,'(a)') ' GEOMETRY'
      WRITE(io,'(a)') ' -------- '
      WRITE(io,'(a)') ' Number of atoms per primitive cell '
      WRITE(io,*) nat
      WRITE(io,'(a)')' Atomic number and position of the atoms(au) '
      DO na = 1, nat
         nt = ityp(na)
         at_num = atomic_number(trim(atm(nt)))
         WRITE(io,'(i6,3f20.14)') at_num, (alat*tau(j,na),j=1,3)
      ENDDO
      WRITE(io,'(a)') ' Primitive lattice vectors (au) '
      WRITE(io,100) alat*at(1,1), alat*at(2,1), alat*at(3,1)
      WRITE(io,100) alat*at(1,2), alat*at(2,2), alat*at(3,2)
      WRITE(io,100) alat*at(1,3), alat*at(2,3), alat*at(3,3)
      WRITE(io,'(a)') ' '

  100 FORMAT (3(1x,f20.15))

   END SUBROUTINE write_header


   SUBROUTINE write_gvecs(g,indx)
      REAL(DP),INTENT(in) :: g(:,:)
      INTEGER,INTENT(in) :: indx(:)
      INTEGER ig

      IF(binwrite)RETURN

      WRITE(io,'(a)') ' G VECTORS'
      WRITE(io,'(a)') ' ---------'
      WRITE(io,'(a)') ' Number of G-vectors'
      WRITE(io,*) size(g,2)
      WRITE(io,'(a)') ' Gx Gy Gz (au)'
      DO ig = 1, size(g,2)
         WRITE(io,'(3(1x,f20.15))') &
         &tpi/alat*g(1,indx(ig)),tpi/alat*g(2,indx(ig)),tpi/alat*g(3,indx(ig))
      ENDDO

      WRITE(io,'(a)') ' '
   END SUBROUTINE write_gvecs


   SUBROUTINE write_gvecs_blip
      IF(binwrite)RETURN

      WRITE(io,'(a)') ' G VECTORS'
      WRITE(io,'(a)') ' ---------'
      WRITE(io,'(a)') ' Number of G-vectors'
      WRITE(io,*) 0
      WRITE(io,'(a)') ' Gx Gy Gz (au)'
      WRITE(io,'(a)') ' Blip grid'
      WRITE(io,'(3(1x,3i4))') blipgrid

      WRITE(io,'(a)') ' '
   END SUBROUTINE write_gvecs_blip


   SUBROUTINE write_wfn_head
      IF(binwrite)RETURN

      WRITE(io,'(a)') ' WAVE FUNCTION'
      WRITE(io,'(a)') ' -------------'
      WRITE(io,'(a)') ' Number of k-points'
      WRITE(io,*) nk
   END SUBROUTINE write_wfn_head


   SUBROUTINE write_pwfn_data(ik,ispin,ibnd,evc,indx)
      INTEGER,INTENT(in) :: ik,ispin,ibnd
      COMPLEX(DP),INTENT(in) :: evc(:)
      INTEGER,INTENT(in) :: indx(:)
      INTEGER ig,j,ikk

      IF(binwrite)RETURN

      ikk = ik + nk*(ispin-1)
      IF(ispin==1.and.ibnd==1)THEN
         WRITE(io,'(a)') ' k-point # ; # of bands (up spin/down spin); &
               &           k-point coords (au)'
         WRITE(io,'(3i4,3f20.16)') ik, nbndup, nbnddown, &
               (tpi/alat*xk(j,ik),j=1,3)
      ENDIF
      IF(binwrite)RETURN

      ! KN: if you want to print occupancies, replace these two lines ...
      WRITE(io,'(a)') ' Band, spin, eigenvalue (au)'
      WRITE(io,*) ibnd, ispin, et(ibnd,ikk)/e2
      ! ...with the following two - KN 2/4/09
      ! WRITE(io,'(a)') ' Band, spin, eigenvalue (au), occupation number'
      ! WRITE(io,*) ibnd, ispin, et(ibnd,ikk)/e2, wg(ibnd,ikk)/wk(ikk)
      WRITE(io,'(a)') ' Eigenvectors coefficients'
      DO ig=1, size(indx,1)
         WRITE(io,*)evc(indx(ig))
      ENDDO
   END SUBROUTINE write_pwfn_data


   SUBROUTINE write_bwfn_data(ik,ispin,ibnd)
      INTEGER,INTENT(in) :: ik,ispin,ibnd
      INTEGER lx,ly,lz,ikk,j,l1,l2,l3

      IF(binwrite)THEN
         DO l3=1,blipgrid(3)
            DO l2=1,blipgrid(2)
               DO l1=1,blipgrid(1)
                  cavc_tmp(l1,l2,l3) = cavc(l1-1,l2-1,l3-1)
               ENDDO
            ENDDO
         ENDDO

         IF(single_precision_blips)THEN
            WRITE(iob)cmplx(cavc_tmp(:,:,:),kind=sgl)
         ELSE
            WRITE(iob)cmplx(cavc_tmp(:,:,:),kind=DP)
         ENDIF
         RETURN
      ENDIF

      ikk = ik + nk*(ispin-1)
      IF(ispin==1.and.ibnd==1)THEN
         WRITE(io,'(a)') ' k-point # ; # of bands (up spin/down spin); &
               &           k-point coords (au)'
         WRITE(io,'(3i4,3f20.16)') ik, nbndup, nbnddown, &
               (tpi/alat*xk(j,ik),j=1,3)
      ENDIF
      ! KN: if you want to print occupancies, replace these two lines ...
      WRITE(io,'(a)') ' Band, spin, eigenvalue (au), localized'
      WRITE(io,*) ibnd, ispin, et(ibnd,ikk)/e2,'F'
      ! ...with the following two - KN 2/4/09
      ! WRITE(io,'(a)') ' Band, spin, eigenvalue (au), occupation number'
      ! WRITE(io,*) ibnd, ispin, et(ibnd,ikk)/e2, wg(ibnd,ikk)/wk(ikk)
      WRITE(io,*)'Complex blip coefficients for extended orbital'
      DO lx=0,blipgrid(1)-1
         DO ly=0,blipgrid(2)-1
            DO lz=0,blipgrid(3)-1
               WRITE(io,*)cavc(lx,ly,lz)
            ENDDO ! lz
         ENDDO ! ly
      ENDDO ! lx
   END SUBROUTINE write_bwfn_data


   SUBROUTINE write_bwfn_data_gamma(re_im,ik,ispin,ibnd)
      INTEGER,INTENT(in) :: ik,ispin,ibnd,re_im
      INTEGER lx,ly,lz,ikk,j,l1,l2,l3

      IF(binwrite)THEN
         IF(re_im==1)THEN
            DO l3=1,blipgrid(3)
               DO l2=1,blipgrid(2)
                  DO l1=1,blipgrid(1)
                     avc_tmp(l1,l2,l3) = avc1(l1-1,l2-1,l3-1)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            DO l3=1,blipgrid(3)
               DO l2=1,blipgrid(2)
                  DO l1=1,blipgrid(1)
                     avc_tmp(l1,l2,l3) = avc2(l1-1,l2-1,l3-1)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

         IF(single_precision_blips)THEN
            WRITE(iob)real(avc_tmp(:,:,:),kind=sgl)
         ELSE
            WRITE(iob)real(avc_tmp(:,:,:),kind=DP)
         ENDIF
         RETURN
      ENDIF

      ikk = ik + nk*(ispin-1)
      IF(ispin==1.and.ibnd==1)THEN
         WRITE(io,'(a)') ' k-point # ; # of bands (up spin/down spin); &
               &           k-point coords (au)'
         WRITE(io,'(3i4,3f20.16)') ik, nbndup, nbnddown, &
               (tpi/alat*xk(j,ik),j=1,3)
      ENDIF
      ! KN: if you want to print occupancies, replace these two lines ...
      WRITE(io,'(a)') ' Band, spin, eigenvalue (au), localized'
      WRITE(io,*) ibnd, ispin, et(ibnd,ikk)/e2,'F'
      ! ...with the following two - KN 2/4/09
      ! WRITE(io,'(a)') ' Band, spin, eigenvalue (au), occupation number'
      ! WRITE(io,*) ibnd, ispin, et(ibnd,ikk)/e2, wg(ibnd,ikk)/wk(ikk)
      WRITE(io,*)'Real blip coefficients for extended orbital'
      DO lx=0,blipgrid(1)-1
         DO ly=0,blipgrid(2)-1
            DO lz=0,blipgrid(3)-1
               IF(re_im==1)THEN
                  WRITE(io,*)avc1(lx,ly,lz)
               ELSE
                  WRITE(io,*)avc2(lx,ly,lz)
               ENDIF
            ENDDO ! lz
         ENDDO ! ly
      ENDDO ! lx
   END SUBROUTINE write_bwfn_data_gamma


   SUBROUTINE create_index2(y,x_index)
      DOUBLE PRECISION,INTENT(in) :: y(:,:)
      INTEGER,INTENT(out) :: x_index(size(y,2))
      DOUBLE PRECISION y2(size(y,2))
      INTEGER i
      DO i = 1,size(y,2)
         y2(i) = sum(y(:,i)**2)
      ENDDO
      CALL create_index(y2,x_index)
   END SUBROUTINE create_index2


   SUBROUTINE create_index(y,x_index)
 !-----------------------------------------------------------------------------!
 ! This subroutine creates an index array x_index for the n items of data in   !
 ! the array y.  Adapted from Numerical Recipes.                               !
 ! Copied from merge_pwfn.f90, included with CASINO distribution               !
 !-----------------------------------------------------------------------------!
      IMPLICIT NONE
      DOUBLE PRECISION,INTENT(in) :: y(:)
      INTEGER,INTENT(out) :: x_index(:)
      INTEGER,PARAMETER :: ins_sort_thresh=7,stacksize=80
      INTEGER n,i,x_indexj,ir,itemp,j,jstack,k,l,lp1,istack(stacksize)
      DOUBLE PRECISION yj
      n=size(x_index)
      DO j=1,n
         x_index(j)=j
      ENDDO ! j
      IF(n<=1)RETURN
      jstack=0
      l=1
      ir=n
      DO
         IF(ir-l<ins_sort_thresh)THEN
    jloop : DO j=l+1,ir
               x_indexj=x_index(j) ; yj=y(x_indexj)
               DO i=j-1,l,-1
                  IF(y(x_index(i))<=yj)THEN
                     x_index(i+1)=x_indexj
                     CYCLE jloop
                  ENDIF! y(x_index(i))<=yj
                  x_index(i+1)=x_index(i)
               ENDDO ! i
               x_index(l)=x_indexj
            ENDDO jloop ! j
            IF(jstack==0)RETURN
            ir=istack(jstack)
            l=istack(jstack-1)
            jstack=jstack-2
         ELSE
            k=(l+ir)/2
            lp1=l+1
            itemp=x_index(k)    ; x_index(k)=x_index(lp1)  ; x_index(lp1)=itemp
            IF(y(x_index(l))>y(x_index(ir)))THEN
               itemp=x_index(l)   ; x_index(l)=x_index(ir)   ; x_index(ir)=itemp
            ENDIF
            IF(y(x_index(lp1))>y(x_index(ir)))THEN
               itemp=x_index(lp1) ; x_index(lp1)=x_index(ir) ; x_index(ir)=itemp
            ENDIF
            IF(y(x_index(l))>y(x_index(lp1)))THEN
               itemp=x_index(l)   ; x_index(l)=x_index(lp1)  ; x_index(lp1)=itemp
            ENDIF
            i=lp1
            j=ir
            x_indexj=x_index(lp1)
            yj=y(x_indexj)
            DO
               DO
                  i=i+1
                  IF(y(x_index(i))>=yj)exit
               ENDDO ! i
               DO
                  j=j-1
                  IF(y(x_index(j))<=yj)exit
               ENDDO ! j
               IF(j<i)exit
               itemp=x_index(i) ; x_index(i)=x_index(j) ; x_index(j)=itemp
            ENDDO
            x_index(lp1)=x_index(j)
            x_index(j)=x_indexj
            jstack=jstack+2
            IF(jstack>stacksize)THEN
               WRITE(6,*)'stacksize is too small.'
               STOP
            ENDIF! jstack>stacksize
            IF(ir-i+1>=j-l)THEN
               istack(jstack)=ir
               istack(jstack-1)=i
               ir=j-1
            ELSE
               istack(jstack)=j-1
               istack(jstack-1)=l
               l=i
            ENDIF! ir-i+1>=j-l
         ENDIF! ir-l<ins_sort_thresh
      ENDDO
   END SUBROUTINE create_index


   CHARACTER(20) FUNCTION i2s(n)
      INTEGER,INTENT(in) :: n
      INTEGER m,j

      m = abs(n)
      DO j=len(i2s),2,-1
         i2s(j:j)=achar(ichar('0')+mod(m,10))
         m=m/10
         IF(m==0)exit
      ENDDO

      IF(n<0)THEN
         j = j-1
         i2s(j:j)='-'
      ENDIF

      i2s=i2s(j:len(i2s))
   END FUNCTION i2s


   CHARACTER(72) FUNCTION write_mean(av,std_err_in_mean,err_prec_in)
!-----------------------------------------------------------------------------!
! Write out a mean value with the standard error in the mean in the form      !
! av(std_err_in_mean), e.g. 0.123546(7).  err_prec_in specifies the number of !
! digits of precision to which the error should be quoted (by default 1).     !
!-----------------------------------------------------------------------------!
      DOUBLE PRECISION,INTENT(in) :: av,std_err_in_mean
      INTEGER,INTENT(in),OPTIONAL :: err_prec_in
      INTEGER lowest_digit_to_quote,err_quote,err_prec,int_part,dec_part,i
      INTEGER,PARAMETER :: err_prec_default=1
      DOUBLE PRECISION av_quote
      CHARACTER(1) sgn
      CHARACTER(72) zero_pad

      IF(std_err_in_mean<=0.d0)THEN
         write_mean='ERROR: NON-POSITIVE ERROR BAR!!!'
         RETURN
      ENDIF ! Error is negative

      IF(present(err_prec_in))THEN
         IF(err_prec_in>=1)THEN
            err_prec=err_prec_in
         ELSE
            write_mean='ERROR: NON-POSITIVE PRECISION!!!'
            RETURN
         ENDIF ! err_prec_in sensible.
      ELSE
         err_prec=err_prec_default
      ENDIF ! Accuracy of error supplied.

! Work out lowest digit of precision that should be retained in the
! mean (i.e. the digit in terms of which the error is specified).
! Calculate the error in terms of this digit and round.
      lowest_digit_to_quote=floor(log(std_err_in_mean)/log(10.d0))+1-err_prec
      err_quote=nint(std_err_in_mean*10.d0**dble(-lowest_digit_to_quote))
      IF(err_quote==10**err_prec)THEN
         lowest_digit_to_quote=lowest_digit_to_quote+1
         err_quote=err_quote/10
      ENDIF ! err_quote rounds up to next figure.

      IF(err_quote>=10**err_prec.or.err_quote<10**(err_prec-1))THEN
         write_mean='ERROR: BUG IN WRITE_MEAN!!!'
         RETURN
      ENDIF ! Check error is in range.

! Truncate the mean to the relevant precision.  Establish its sign,
! then take the absolute value and work out the integer part.
      av_quote=anint(av*10.d0**dble(-lowest_digit_to_quote)) &
         &*10.d0**dble(lowest_digit_to_quote)
      IF(av_quote<0.d0)THEN
         sgn='-'
         av_quote=-av_quote
      ELSE
         sgn=' '
      ENDIF ! Sign
      IF(aint(av_quote)>dble(huge(1)))THEN
         write_mean='ERROR: NUMBERS ARE TOO LARGE IN WRITE_MEAN!'
         RETURN
      ENDIF ! Vast number
      int_part=floor(av_quote)

      IF(lowest_digit_to_quote<0)THEN
! If the error is in a decimal place then construct string using
! integer part and decimal part, noting that the latter may need to
! be padded with zeros, e.g. if we want "0001" rather than "1".
         IF(anint((av_quote-dble(int_part)) &
            &*10.d0**dble(-lowest_digit_to_quote))>dble(huge(1)))THEN
            write_mean='ERROR: NUMBERS ARE TOO LARGE IN WRITE_MEAN!'
            RETURN
         ENDIF ! Vast number
         dec_part=nint((av_quote-dble(int_part))*10.d0**dble(-lowest_digit_to_quote))
         zero_pad=' '
         IF(dec_part<0)THEN
            write_mean='ERROR: BUG IN WRITE_MEAN! (2)'
            RETURN
         ENDIF ! dec
         DO i=1,-lowest_digit_to_quote-no_digits_int(dec_part)
            zero_pad(i:i)='0'
         ENDDO ! i
         write_mean=sgn//trim(i2s(int_part))//'.'//trim(zero_pad) &
         &//trim(i2s(dec_part))//'('//trim(i2s(err_quote))//')'
      ELSE
! If the error is in a figure above the decimal point then, of
! course, we don't have to worry about a decimal part.
         write_mean=sgn//trim(i2s(int_part))//'(' &
            &//trim(i2s(err_quote*10**lowest_digit_to_quote))//')'
      ENDIF ! lowest_digit_to_quote<0

   END FUNCTION write_mean


   INTEGER FUNCTION no_digits_int(i)
   !----------------------------------------------------------------------!
   ! Calculate the number of digits in integer i.  For i>0 this should be !
   ! floor(log(i)/log(10))+1, but sometimes rounding errors cause this    !
   ! expression to give the wrong result.                                 !
   !----------------------------------------------------------------------!
      INTEGER,INTENT(in) :: i
      INTEGER j,k
      j=i ; k=1
      DO
         j=j/10
         IF(j==0)exit
         k=k+1
      ENDDO
      no_digits_int=k
   END FUNCTION no_digits_int



   SUBROUTINE init_rng(seed)
!--------------------------------------------!
! Initialize the RNG: see Knuth's ran_start. !
!--------------------------------------------!
      INTEGER,INTENT(in) :: seed
      INTEGER j,s,t,sseed
      INTEGER,PARAMETER :: MM=2**30,TT=70
      REAL(DP) ss,x(KK+KK-1)
      REAL(DP),PARAMETER :: ULP=1.d0/2.d0**52,ULP2=2.d0*ULP
      IF(seed<0)THEN
         sseed=MM-1-mod(-1-seed,MM)
      ELSE
         sseed=mod(seed,MM)
      ENDIF ! seed<0
      ss=ULP2*dble(sseed+2)
      DO j=1,KK
         x(j)=ss
         ss=ss+ss
         IF(ss>=1.d0)ss=ss-1.d0+ULP2
      ENDDO ! j
      x(2)=x(2)+ULP
      s=sseed
      t=TT-1
      DO
         DO j=KK,2,-1
            x(j+j-1)=x(j)
            x(j+j-2)=0.d0
         ENDDO ! j
         DO j=KK+KK-1,KK+1,-1
            x(j-(KK-LL))=mod(x(j-(KK-LL))+x(j),1.d0)
            x(j-KK)=mod(x(j-KK)+x(j),1.d0)
         ENDDO ! j
         IF(mod(s,2)==1)THEN
            DO j=KK,1,-1
               x(j+1)=x(j)
            ENDDO ! j
            x(1)=x(KK+1)
            x(LL+1)=mod(x(LL+1)+x(KK+1),1.d0)
         ENDIF ! s odd
         IF(s/=0)THEN
            s=s/2
         ELSE
            t=t-1
         ENDIF ! s/=0
         IF(t<=0)exit
      ENDDO
      ranstate(1+KK-LL:KK)=x(1:LL)
      ranstate(1:KK-LL)=x(LL+1:KK)
      DO j=1,10
         CALL gen_ran_array(x,KK+KK-1)
      ENDDO ! j
      ran_array_idx=Nkeep
   END SUBROUTINE init_rng


   REAL(dp) FUNCTION ranx()
!------------------------------------------------------------------------------!
! Return a random number uniformly distributed in [0,1).                       !
! Uses M. Luescher's suggestion: generate 1009 random numbers at a time using  !
! Knuth's algorithm, but only use the first 100.                               !
!------------------------------------------------------------------------------!
      IF(ran_array_idx==-1)THEN
         CALL init_rng(default_seed) ! Initialize the RNG.
      ENDIF ! First call.
      IF(ran_array_idx==Nkeep)THEN
         CALL gen_ran_array(ran_array,Nran) ! Generate a new array of random nos.
         ran_array_idx=0
      ENDIF ! i=Nkeep
      ran_array_idx=ran_array_idx+1
      ranx=ran_array(ran_array_idx)
   END FUNCTION ranx


   SUBROUTINE gen_ran_array(ran_array,N)
!---------------------------------------------------------------!
! Generate an array of N random numbers: see Knuth's ran_array. !
!---------------------------------------------------------------!
      INTEGER,INTENT(in) :: N
      REAL(DP),INTENT(out) :: ran_array(N)
      INTEGER j
      ran_array(1:KK)=ranstate(1:KK)
      DO j=KK+1,N
         ran_array(j)=mod(ran_array(j-KK)+ran_array(j-LL),1.d0)
      ENDDO ! j
      DO j=1,LL
         ranstate(j)=mod(ran_array(N+j-KK)+ran_array(N+j-LL),1.d0)
      ENDDO ! j
      DO j=LL+1,KK
         ranstate(j)=mod(ran_array(N+j-KK)+ranstate(j-LL),1.d0)
      ENDDO ! j
   END SUBROUTINE gen_ran_array


END SUBROUTINE write_casino_wfn
