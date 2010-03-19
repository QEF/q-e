!
! Copyright (C) 2004-2009 Dario Alfe' and Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE write_casino_wfn(gather,blip,multiplicity,binwrite,single_precision_blips)

   USE kinds, ONLY: DP
   USE ions_base, ONLY : nat, ntyp => nsp, ityp, tau, zv, atm
   USE cell_base, ONLY: omega, alat, tpiba2, at, bg
   USE printout_base, ONLY: title    ! title of the run
   USE constants, ONLY: tpi, e2
   USE ener, ONLY: ewld, ehart, etxc, vtxc, etot, etxcc, demet, ef
   USE gvect, ONLY: ngm, gstart, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                    nrxx, g, gg, ecutwfc, gcutm, nl, nlm, igtongl
   USE klist , ONLY: nks, nelec, xk, wk, degauss, ngauss
   USE lsda_mod, ONLY: lsda, nspin
   USE scf, ONLY: rho, rho_core, rhog_core, vnew
   USE ldaU, ONLY : eth
   USE vlocal, ONLY: vloc, strf
   USE wvfct, ONLY: npw, npwx, nbnd, igk, g2kin, wg, et
   USE control_flags, ONLY : gamma_only
   USE uspp, ONLY: nkb, vkb, dvan
   USE uspp_param, ONLY: nh
   USE io_global, ONLY: stdout, ionode, ionode_id
   USE io_files, ONLY: nd_nmbr, nwordwfc, iunwfc
   USE wavefunctions_module, ONLY : evc
   USE funct, ONLY : dft_is_meta
   USE mp_global, ONLY: inter_pool_comm, intra_pool_comm, nproc_pool, me_pool
   USE mp, ONLY: mp_sum, mp_gather, mp_bcast, mp_get
   USE dfunct, ONLY : newd

   USE pw2blip

   IMPLICIT NONE
   LOGICAL, INTENT(in) :: gather,blip,binwrite,single_precision_blips
   REAL(dp), INTENT(in) :: multiplicity

   REAL(dp), PARAMETER :: eps = 1.d-10
   INTEGER, PARAMETER :: io = 77, iob = 78
   INTEGER :: ig, ibnd, ik, ispin, nbndup, nbnddown, &
              nk, ig7, ikk, id, ip, iorb, iorb_node, inode, ierr, norb
   INTEGER :: jk(nproc_pool), jspin(nproc_pool), jbnd(nproc_pool)
   INTEGER :: jk2(nproc_pool), jspin2(nproc_pool), jbnd2(nproc_pool)
   INTEGER, ALLOCATABLE :: idx(:), igtog(:)
   LOGICAL :: exst,dowrite
   REAL(DP) :: ek, eloc, enl
   INTEGER, EXTERNAL :: atomic_number
   REAL (DP), EXTERNAL :: ewald, w1gauss

   ! number of g vectors (union of all k points)
   INTEGER ngtot_l ! on this processor
   INTEGER, ALLOCATABLE :: ngtot_d(:), ngtot_cumsum(:), indx(:)
   INTEGER ngtot_g ! sum over processors
   REAL(DP), ALLOCATABLE :: g_l(:,:), g_g(:,:), g2(:)
   COMPLEX(DP), ALLOCATABLE :: evc_l(:), evc_g(:), evc_g2(:), avc_tmp(:,:,:), cavc_tmp(:,:,:)
   LOGICAL dotransform

   CALL init_us_1
   CALL newd

   dowrite=ionode.or..not.(gather.or.blip)

   ALLOCATE (idx (ngm) )
   ALLOCATE (igtog (ngm) )
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
         CALL gk_sort (xk (1:3, ikk), ngm, g(1:3,1:ngm), ecutwfc / tpiba2, & ! input
                      &npw, igk, g2kin)                                      ! output
         idx( igk(1:npw) ) = 1
      ENDDO
   ENDDO

   ngtot_l = 0
   DO ig = 1, ngm
      IF( idx(ig) >= 1 )THEN
         ngtot_l = ngtot_l + 1
         igtog(ngtot_l) = ig
      ENDIF
   ENDDO

   DEALLOCATE (idx)

   if(gamma_only)then
      blipreal=-2
   elseif(nk==1.and.all(abs(xk(1:3,1))<eps))
      blipreal=2
   else
      blipreal=0
   endif

   IF(dowrite)THEN
      IF(blip)THEN
         IF(binwrite)THEN
            WRITE (6,'(/,5x,''Writing file bwfn.data.b1 for program CASINO'')')
            CALL seqopn( iob, 'bwfn.data.b1', 'unformatted',exst)
         ELSE
            WRITE (6,'(/,5x,''Writing file bwfn.data for program CASINO'')')
            CALL seqopn( io, 'bwfn.data', 'formatted',exst)
         ENDIF
      ELSE
         WRITE (6,'(/,5x,''Writing file pwfn.data for program CASINO'')')
         CALL seqopn( io, 'pwfn.data', 'formatted',exst)
      ENDIF
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

      ALLOCATE ( g_g(3,ngtot_g), evc_g(ngtot_g) )
      IF(blip.and.blipreal/=0)THEN
         ALLOCATE( evc_g2(ngtot_g) )
      ENDIF
      CALL mp_gather( g_l, g_g, ngtot_d, ngtot_cumsum, ionode_id, intra_pool_comm)

      IF(blip)THEN
         CALL mp_bcast( g_g, ionode_id, intra_pool_comm )
         CALL pw2blip_init(ngtot_g,g_g,multiplicity)
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
      if(blipreal/=0)then
         ALLOCATE(avc_tmp(blipgrid(1),blipgrid(2),blipgrid(3)))
      else
         ALLOCATE(cavc_tmp(blipgrid(1),blipgrid(2),blipgrid(3)))
      endif
   ENDIF

   ! making some assumptions about the parallel layout:
   IF(ionode_id/=0)CALL errore('write_casino_wfn','ionode_id/=0: ',ionode_id)

   iorb = 0
   norb = nk*nspin*nbnd

   DO ik = 1, nk
      DO ispin = 1, nspin
         ikk = ik + nk*(ispin-1)
         IF( nks > 1 )THEN
            CALL gk_sort (xk (1:3, ikk), ngm, g(1:3,1:ngm), ecutwfc / tpiba2, & ! input
                         &npw, igk, g2kin)                                      ! output
            CALL davcio(evc,nwordwfc,iunwfc,ikk,-1)
         ENDIF
         DO ibnd = 1, nbnd
            evc_l(:) = (0.d0, 0d0)
            DO ig=1, ngtot_l
               ! now for all G vectors find the PW coefficient for this k-point
               find_ig: DO ig7 = 1, npw
                  IF( igk(ig7) == igtog(ig) )THEN
                     evc_l(ig) = evc(ig7,ibnd)
                     exit find_ig
                  ENDIF
               ENDDO find_ig
            ENDDO
            IF(blip)THEN
               iorb = iorb + 1
               if(blipreal/=0)then
                  iorb_node = mod((iorb-1)/2,nproc_pool) ! the node that should compute this orbital
                  if(mod(iorb,2)==1)then
                     jk(iorb_node+1) = ik
                     jspin(iorb_node+1) = ispin
                     jbnd(iorb_node+1) = ibnd
                     dotransform = .false.
                  else
                     jk2(iorb_node+1) = ik
                     jspin2(iorb_node+1) = ispin
                     jbnd2(iorb_node+1) = ibnd
                     dotransform=(iorb_node==nproc_pool-1)
                  endif
               else
                  iorb_node = mod(iorb-1,nproc_pool) ! the node that should compute this orbital
                  jk(iorb_node+1) = ik
                  jspin(iorb_node+1) = ispin
                  jbnd(iorb_node+1) = ibnd
                  dotransform=(iorb_node==nproc_pool-1)
               endif
               DO inode=0,nproc_pool-1
                  if(blipreal/=0.and.mod(iorb,2)==0)then
                     CALL mp_get(&
                        evc_g2(ngtot_cumsum(inode+1)+1:ngtot_cumsum(inode+1)+ngtot_d(inode+1)),&
                        evc_l(:),me_pool,iorb_node,inode,1234,intra_pool_comm)
                  else
                     CALL mp_get(&
                        evc_g(ngtot_cumsum(inode+1)+1:ngtot_cumsum(inode+1)+ngtot_d(inode+1)),&
                        evc_l(:),me_pool,iorb_node,inode,1234,intra_pool_comm)
                  endif
               ENDDO
               IF(dotransform .or. iorb == norb)THEN
                  IF(me_pool <= iorb_node)THEN
                     IF(blipreal.and.(me_pool/=iorb_node.or.iorb/=norb.or.mod(norb,2)==0)then
                        CALL pw2blip_transform2(evc_g(:),evc_g2(:))
                     ELSE
                        CALL pw2blip_transform(evc_g(:))
                     ENDIF
                  ENDIF
                  DO inode=0,iorb_node
                     CALL pw2blip_get(inode)
                     if(ionode)then
                        if(blipreal/=0)then
                           write(6,*)"Transformed real orbital k=",jk(inode+1),",spin=",jspin(inode+1),&
                              &",band=",jbnd(inode+1))," on node ",inode
                           CALL pw2blip_stat(inode)
                           CALL write_bwfn_data_gamma(1,jk(inode+1),jspin(inode+1),jbnd(inode+1))
                           if(modulo(blipreal,2)==0)then
                              write(6,*)"Transformed real orbital k=",jk(inode+1),",spin=",jspin(inode+1),&
                                 &",band=",jbnd(inode+1))," on node ",inode
                              CALL pw2blip_stat2(inode)
                              CALL write_bwfn_data_gamma(2,jk2(inode+1),jspin2(inode+1),jbnd2(inode+1))
                           endif
                        else
                           write(6,*)"Transformed complex orbital k=",jk(inode+1),",spin=",jspin(inode+1),&
                              &",band=",jbnd(inode+1))," on node ",inode
                           CALL pw2blip_stat(inode)
                           CALL write_bwfn_data(jk(inode+1),jspin(inode+1),jbnd(inode+1))
                        endif
                     endif
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
      if(blipreal)then
         DEALLOCATE(avc_tmp)
      else
         DEALLOCATE(cavc_tmp)
      endif
   ENDIF

   IF(blip)CALL pw2blip_cleanup
   DEALLOCATE (igtog, g_l, evc_l )
   IF(blip.or.gather) DEALLOCATE ( ngtot_d, ngtot_cumsum, g_g, evc_g )
   IF(dowrite.and..not.blip) DEALLOCATE (indx)

CONTAINS

   SUBROUTINE calc_energies
      USE becmod, ONLY: becp, calbec, allocate_bec_type, deallocate_bec_type

      COMPLEX(DP), ALLOCATABLE :: aux(:)
      INTEGER :: ibnd, j, ig, ik, ikk, ispin, na, nt, ijkb0, ikb, ih, jh, jkb

      REAL(DP) :: charge, etotefield

      ALLOCATE (aux(nrxx))
      CALL allocate_bec_type ( nkb, nbnd, becp )

      ek  = 0.d0
      eloc= 0.d0
      enl = 0.d0
      demet=0.d0
      !
      DO ispin = 1, nspin
         !
         !     calculate the local contribution to the total energy
         !
         !      bring rho to G-space
         !
         aux(:) = cmplx( rho%of_r(:,ispin), 0.d0,kind=DP)
         CALL cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,-1)
         !
         DO nt=1,ntyp
            DO ig = 1, ngm
               if(gamma_only)then !.and.ig>1)then
                  eloc = eloc + vloc(igtongl(ig),nt) * strf(ig,nt) &
                     * conjg(aux(nl(ig)))
                  eloc = eloc + vloc(igtongl(ig),nt) * strf(ig,nt) &
                     * conjg(aux(nlm(ig)))
               else
                  eloc = eloc + vloc(igtongl(ig),nt) * strf(ig,nt) &
                     * conjg(aux(nl(ig)))
               endif
            ENDDO
         ENDDO

         DO ik = 1, nk
            ikk = ik + nk*(ispin-1)
            CALL gk_sort (xk (1, ikk), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
            CALL davcio (evc, nwordwfc, iunwfc, ikk, - 1)
            CALL init_us_2 (npw, igk, xk (1, ikk), vkb)
            CALL calbec ( npw, vkb, evc, becp )
            !
            ! -TS term for metals (ifany)
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
            DO ibnd = 1, nbnd
               DO j = 1, npw
                  if(gamma_only)then !.and.j>1)then
                     ek = ek +  2*conjg(evc(j,ibnd)) * evc(j,ibnd) * &
                                    g2kin(j) * wg(ibnd,ikk)
                  else
                     ek = ek +  conjg(evc(j,ibnd)) * evc(j,ibnd) * &
                                    g2kin(j) * wg(ibnd,ikk)
                  endif
               ENDDO

               !
               ! Calculate Non-local energy
               !
               ijkb0 = 0
               DO nt = 1, ntyp
                  DO na = 1, nat
                     IF(ityp (na) .eq. nt)THEN
                        DO ih = 1, nh (nt)
                           ikb = ijkb0 + ih
                           if(gamma_only)then
                              enl=enl+becp%r(ikb,ibnd)*becp%r(ikb,ibnd) &
                                 *wg(ibnd,ikk)* dvan(ih,ih,nt)
                           else
                              enl=enl+conjg(becp%k(ikb,ibnd))*becp%k(ikb,ibnd) &
                                 *wg(ibnd,ikk)* dvan(ih,ih,nt)
                           endif
                           DO jh = ( ih + 1 ), nh(nt)
                              jkb = ijkb0 + jh
                              if(gamma_only)then
                                 enl=enl + &
                                    (becp%r(ikb,ibnd)*becp%r(jkb,ibnd)+&
                                       becp%r(jkb,ibnd)*becp%r(ikb,ibnd))&
                                    * wg(ibnd,ikk) * dvan(ih,jh,nt)
                              else
                                 enl=enl + &
                                    (conjg(becp%k(ikb,ibnd))*becp%k(jkb,ibnd)+&
                                       conjg(becp%k(jkb,ibnd))*becp%k(ikb,ibnd))&
                                    * wg(ibnd,ikk) * dvan(ih,jh,nt)
                              endif

                           ENDDO

                        ENDDO
                        ijkb0 = ijkb0 + nh (nt)
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

#ifdef __PARA
      CALL mp_sum( eloc,  intra_pool_comm )
      CALL mp_sum( ek,    intra_pool_comm )
      CALL mp_sum( ek,    inter_pool_comm )
      CALL mp_sum( enl,   inter_pool_comm )
      CALL mp_sum( demet, inter_pool_comm )
#endif
      eloc = eloc * omega
      ek = ek * tpiba2

      !
      ! compute ewald contribution
      !
      ewld = ewald( alat, nat, ntyp, ityp, zv, at, bg, tau, omega, &
           g, gg, ngm, gcutm, gstart, gamma_only, strf )
      !
      ! compute hartree and xc contribution
      !
      CALL v_of_rho( rho, rho_core, rhog_core, &
                     ehart, etxc, vtxc, eth, etotefield, charge, vnew )
      !
      etot=(ek + (etxc-etxcc)+ehart+eloc+enl+ewld)+demet

      CALL deallocate_bec_type (becp)
      DEALLOCATE (aux)

      WRITE (stdout,*) 'Kinetic energy   ', ek/e2
      WRITE (stdout,*) 'Local energy     ', eloc/e2
      WRITE (stdout,*) 'Non-Local energy ', enl/e2
      WRITE (stdout,*) 'Ewald energy     ', ewld/e2
      WRITE (stdout,*) 'xc contribution  ',(etxc-etxcc)/e2
      WRITE (stdout,*) 'hartree energy   ', ehart/e2
      IF( degauss > 0.0_dp ) &
         WRITE (stdout,*) 'Smearing (-TS)   ', demet/e2
      WRITE (stdout,*) 'Total energy     ', etot/e2


   END SUBROUTINE calc_energies

   FUNCTION to_c80(c)
      CHARACTER(*),INTENT(in) :: c
      CHARACTER(80) :: to_c80
      to_c80=c
   END FUNCTION to_c80

   SUBROUTINE write_header
      INTEGER j, na, nt, at_num
      REAL(dp) :: kvec(3,nk),ksq(nk),kprod(6,nk)

      IF(binwrite)THEN
         WRITE(iob)&
            to_c80(title)    ,& ! title
            to_c80("PWSCF")  ,& ! code
            to_c80("DFT")    ,& ! method
            to_c80("unknown"),& ! functional
            to_c80("unknown"),& ! pseudo_type
            dble(ecutwfc/2)  ,&  ! plane_wave_cutoff
            lsda             ,&  ! spin_polarized,
            dble(etot/e2)    ,&  ! total_energy
            dble(ek/e2)      ,&  ! kinetic_energy
            dble(eloc/e2)    ,&  ! local_potential_energy
            dble(enl/e2)     ,&  ! non_local_potential_energy
            dble(ehart/e2)   ,&  ! electron_electron_energy
            dble(ewld/e2)    ,&  ! eionion
            nint(nelec)      ,&  ! num_electrons
            nat              ,&  ! nbasis
            ngtot_g          ,&  ! nwvec
            nk               ,&  ! nkvec
            blipgrid(1:3)    ,&  ! nr
            nbnd             ,&  ! maxband
            blipreal         ,&  ! gamma_only
            .true.           ,&  ! ext_orbs_present
            (/0,0/)          ,&  ! no_loc_orbs
            alat*at(1:3,1)   ,&  ! pa1
            alat*at(1:3,2)   ,&  ! pa2
            alat*at(1:3,3)   ,&  ! pa3
            2                ,&  ! nspin_check
            nbnd                 ! num_nonloc_max

         kvec(:,:) = tpi/alat*xk(1:3,1:nk)
         kprod(1,:)=kvec(1,:)*kvec(1,:)
         kprod(2,:)=kvec(2,:)*kvec(2,:)
         kprod(3,:)=kvec(3,:)*kvec(3,:)
         kprod(4,:)=kvec(1,:)*kvec(2,:)
         kprod(5,:)=kvec(1,:)*kvec(3,:)
         kprod(6,:)=kvec(2,:)*kvec(3,:)
         ksq(:)=kprod(1,:)+kprod(2,:)+kprod(3,:)

         WRITE(iob)&
            kvec                                          ,& ! kvec
            ksq                                           ,& ! ksq
            kprod                                         ,& ! kprod
            (atomic_number(trim(atm(ityp(na)))),na=1,nat) ,& ! atno   -- atomic numbers
            (alat*tau(1:3,na),na=1,nat)                   ,& ! basis  -- atom positions
            (nbnd,j=1,nk*2)                               ,& ! nband
            et(1:nbnd,1:nk*nspin)/e2                      ,& ! eigenvalue
            (.true.,j=1,nbnd*nk*nspin)                    ,& ! on_this_cpu
            (/nbnd,nbnd/)                                    ! num_nonloc
         WRITE(iob)single_precision_blips                    ! single_precision_blips

         ! IF(no_loc_orbs>0)THEN
         !    ...
         ! ENDIF

         WRITE(iob)&
          (0,j=1,nbnd*nk*2) ,& ! orb_map_band
          (0,j=1,nbnd*nk*2) ,& ! orb_map_ik
          (0,j=1,nbnd*nk*2) ,& ! orb_map_iorb
          (0,j=1,nbnd*nk*2)    ! occupied
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
                  cavc_tmp(l1,l2,l3) = cavc(l1,l2,l3)
               ENDDO
            ENDDO
         ENDDO

         IF(single_precision_blips)THEN
            WRITE(iob)cmplx(cavc_tmp(:,:,:),kind=kind(1.))
         ELSE
            WRITE(iob)cmplx(cavc_tmp(:,:,:),kind=kind(1.d0))
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
      DO lx=1,blipgrid(1)
         DO ly=1,blipgrid(2)
            DO lz=1,blipgrid(3)
               WRITE(io,*)cavc(lx,ly,lz)
            ENDDO ! lz
         ENDDO ! ly
      ENDDO ! lx
   END SUBROUTINE write_bwfn_data


   SUBROUTINE write_bwfn_data_gamma(re_im,ik,ispin,ibnd)
      INTEGER,INTENT(in) :: ik,ispin,ibnd,re_im
      INTEGER lx,ly,lz,ikk,j,l1,l2,l3

      IF(binwrite)THEN
         if(re_im==1)then
            DO l3=1,blipgrid(3)
               DO l2=1,blipgrid(2)
                  DO l1=1,blipgrid(1)
                     avc_tmp(l1,l2,l3) = avc1(l1,l2,l3)
                  ENDDO
               ENDDO
            ENDDO
         else
            DO l3=1,blipgrid(3)
               DO l2=1,blipgrid(2)
                  DO l1=1,blipgrid(1)
                     avc_tmp(l1,l2,l3) = avc2(l1,l2,l3)
                  ENDDO
               ENDDO
            ENDDO
         endif

         IF(single_precision_blips)THEN
            WRITE(iob)real(avc_tmp(:,:,:),kind=kind(1.))
         ELSE
            WRITE(iob)real(avc_tmp(:,:,:),kind=kind(1.d0))
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
      DO lx=1,blipgrid(1)
         DO ly=1,blipgrid(2)
            DO lz=1,blipgrid(3)
               if(re_im==1)then
                  WRITE(io,*)avc1(lx,ly,lz)
               else
                  WRITE(io,*)avc2(lx,ly,lz)
               endif
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

END SUBROUTINE write_casino_wfn
