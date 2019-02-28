  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! 
  !---------------------------------------------------------------------
  SUBROUTINE elphel2_shuffle( npe, imode0, dvscfins, gmapsym, eigv, isym, xq0, timerev )
  !---------------------------------------------------------------------
  !!
  !!      Calculation of the electron-phonon matrix elements el_ph_mat
  !!      <\psi(k+q)|dV_{SCF}/du^q_{i a}|\psi(k)>
  !!
  !!      Written by Feliciano Giustino based on the routine PH/elphon.f90/elphel. 
  !!      Main difference w.r.t. to original routine is gauge fixing, 
  !!      shuffle (umklapp) mode and all-q implementation.
  !!
  !!      Shuffle mode implemented on may 7 2006
  !!
  !!      Nota Bene: this subroutine is intended only for one proc per pool, 
  !!      i.e. with no G-vector parallelization (some work on the igkq is 
  !!      required for that in the g-mapping)
  !!
  !!      In order to allow a pool reading the wfc file of another
  !!      pool, I had to modify the bound npwx in PW/n_plane_waves.f90
  !!      which is now the max across all pools. In this way lrwfc is
  !!      the same for all pools.
  !!
  !!      RM - Nov/Dec 2014
  !!      Imported the noncolinear case implemented by xlzhang
  !!
  !!      SP - Nov 2015
  !!      We want g(k,Sq) = < k+S(q)(r) | V_S(q)(r) | k(r) >
  !!                      = < k+S(q)(r) | V_q({S|v}^-1 r) | k(r) > 
  !!                      = < k+S(q)({S|v}r) | V_q (r) | k({S|v}r) > 
  !!
  !!      It is important to note that the KB projectors that are applied to the V need
  !!      to be computed at (r) and not ({S|v}r). Therefore, for the KB proj (computed in 
  !!      init_us_2, we need to provide the < Sk+q (r)| and |Sk (r)>.
  !!      See Eq. 11.40 and 11.41 of the R. Martin Electronic Structure book.
  !! 
  !!      Note that in QE Sq is defined as S^-1(q)                
  !! 
  !!      In case of time-reversal
  !!      ------------------------
  !!       g(k,-Sq) = < k-S(q)({S|v}r) | V^loc_-q (r) + (V^nloc_q)* | k({S|v}r) > 
  !!       where V^loc_{-q} is obtained with setlocq and V^nloc_q = CONGJ(u_pattern)*dvscfins*u_pattern.
  !!       We have to do this splitting because we do not have V^nloc_-q and
  !!       V^loc has to be computed at -q to be mappable with the vkb of the wavefunctions
  !!       computed in init_us_2.
  !! 
  !!      Roxana Margine - Jan 2019: Updated based on QE 6.3 for US 
  !!
  !---------------------------------------------------------------------
  !
  USE kinds,         ONLY : DP 
  USE mp_global,     ONLY : my_pool_id, nproc_pool, intra_pool_comm, &
                            inter_pool_comm, inter_image_comm, world_comm
  USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
  USE io_global,     ONLY : stdout
  USE wavefunctions, ONLY : evc
  USE io_files,      ONLY : diropn, seqopn
  USE wvfct,         ONLY : npwx
  USE pwcom,         ONLY : current_spin, isk, lsda, nbnd, xk, nks
  USE cell_base,     ONLY : tpiba
  USE gvect,         ONLY : ngm, g
  USE uspp,          ONLY : vkb
  USE symm_base,     ONLY : s
  USE modes,         ONLY : u  
  USE qpoint,        ONLY : xq, npwq
  USE eqv,           ONLY : dvpsi, evq
  USE units_lr,      ONLY : lrwfc, iuwfc
  USE phus,          ONLY : alphap
  USE lrus,          ONLY : becp1
  USE becmod,        ONLY : calbec 
  USE elph2,         ONLY : shift, gmap, el_ph_mat, umat, umatq, igk_k_all, &
                            umat_all, xk_all, et_all, xkq, etq, igkq, igk, &
                            ngk_all, lower_band, upper_band
  USE fft_base,      ONLY : dffts
  USE constants_epw, ONLY : czero, cone, ci, zero
  USE control_flags, ONLY : iverbosity
  USE klist,         ONLY : nkstot
  USE division,      ONLY : kpointdivision, fkbounds, fkbounds_bnd
  USE noncollin_module, ONLY : noncolin, npol, nspin_mag
  ! 
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: npe
  !! Number of perturbations for this irr representation
  INTEGER, INTENT(in) :: imode0
  !! Current mode number
  INTEGER, INTENT(in) :: gmapsym(ngm,48)
  !! Correspondence  G->S(G)  
  INTEGER, INTENT(in) :: isym
  !! The symmetry which generates the current q in the star  
  REAL(kind=DP), INTENT(in) :: xq0(3)
  !! The first q-point in the star (cartesian coords.)
  COMPLEX(kind=DP), INTENT(in) :: dvscfins(dffts%nnr, nspin_mag, npe)
  !! Delta scf potential
  COMPLEX(kind=DP), INTENT(in) :: eigv (ngm,48)
  !! $e^{iGv}$ for 1...nsym (v the fractional translation)
  LOGICAL, INTENT(in) :: timerev
  !!  true if we are using time reversal
  !
  ! Local variables
  !
  LOGICAL :: exst
  !! logical variable to check file exists
  !
  INTEGER :: ik
  !! Counter on k-points in the pool
  INTEGER :: ik0
  !! Index of the first k-point block in this pool - 1
  INTEGER :: ibnd
  !! Counter on bands
  INTEGER :: jbnd
  !! Counter on bands
  INTEGER :: ipert
  !! Counter on change of Vscf due to perturbations
  INTEGER :: mode
  !! Counter on modes plus pertubations
  INTEGER :: ig
  !! Counter on G-vectors
  INTEGER :: ipooltmp
  !! Index of pool for k
  INTEGER :: ipool
  !! Index of pool for k+q
  INTEGER :: igkq_tmp(npwx)
  !! Correspondence k+q+G <-> G 
  INTEGER :: imap
  !! Index in gmap for G-sphere igkq translation with G_0
  INTEGER :: ipol
  !! Counter on polarizations
  INTEGER :: npw
  !! Number of k+G-vectors inside 'ecut sphere' 
  INTEGER :: ng0vec
  !! Number of G_0 vectors 
  INTEGER :: ngxx
  !! Maximum number of G-vectors over all pools
  INTEGER :: lower_bnd
  !! Lower bounds index after k paral
  INTEGER :: upper_bnd
  !! Upper bounds index after k paral
  INTEGER :: nkk
  !! Index of k-point in the pool
  INTEGER :: nkk_abs
  !! Absolute index of k-point
  INTEGER :: nkq
  !! Index of k+q-point in the pool
  INTEGER :: nkq_abs
  !! Absolute index of k+q-point
  !
  ! Local variables for rotating the wavefunctions (in order to use q in the irr wedge)
  REAL(kind=DP) :: xkqtmp(3)
  !! Temporary k+q vector for KB projectors 
  REAL(kind=DP) :: sxk(3)
  !! Rotated k-point xk
  REAL(kind=DP) :: g0vec_all_r(3,125)
  !! G_0 vectors needed to fold the k+q grid into the k grid, cartesian coord.
  REAL(kind=DP) :: zero_vect(3)
  !! Temporary zero vector 
  ! 
  COMPLEX(kind=DP), ALLOCATABLE :: aux1(:,:), aux2(:,:), aux3(:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: eptmp(:,:), elphmat(:,:,:)
  !! arrays for e-ph matrix elements 
  !
!DBSP - NAG complains ...
  COMPLEX(kind=DP), EXTERNAL :: zdotc
!DBSP
!  REAL(kind=DP) :: b, c, d
!END
  !
  IF ( .not. ALLOCATED(elphmat) ) ALLOCATE( elphmat(nbnd, nbnd, npe) ) 
  IF ( .not. ALLOCATED(eptmp) )   ALLOCATE( eptmp(nbnd, nbnd) ) 
  IF ( .not. ALLOCATED(aux1) )    ALLOCATE( aux1(dffts%nnr, npol) )
  IF ( .not. ALLOCATED(aux2) )    ALLOCATE( aux2(npwx*npol, nbnd) )
  elphmat(:,:,:) = czero
  eptmp(:,:) = czero
  aux1(:,:) = czero
  aux2(:,:) = czero
  zero_vect = zero
  !
  IF (ALLOCATED(xkq) ) DEALLOCATE(xkq)                  
  IF (.not. ALLOCATED(xkq) ) ALLOCATE( xkq(3,nkstot) ) 
  xkq(:,:) = zero
  !
  IF ( nproc_pool>1 ) CALL errore &
    ('elphel2_shuffle', 'only one proc per pool in shuffle mode', 1)
  !
  ! find the bounds of k-dependent arrays in the parallel case in each pool
  CALL fkbounds( nkstot, lower_bnd, upper_bnd )
  !
  ! SP: Bound for band parallelism
  CALL fkbounds_bnd( nbnd, lower_band, upper_band )
  !
  IF ( .not. ALLOCATED(aux3) )  ALLOCATE( aux3(npwx*npol, lower_band:upper_band) )
  IF ( .not. ALLOCATED(dvpsi) ) ALLOCATE( dvpsi(npwx*npol, lower_band:upper_band) )
  aux3(:,:) = czero
  dvpsi(:,:) = czero
  !
  ! setup for k+q folding
  !
  CALL kpointdivision( ik0 )
  CALL readgmap( nkstot, ngxx, ng0vec, g0vec_all_r, lower_bnd )
  !
  IF (imode0.eq.0 .AND. iverbosity.eq.1) WRITE(stdout,5) ngxx
5 FORMAT (5x,'Estimated size of gmap: ngxx =',i5)
  !
  ! close all sequential files in order to re-open them as direct access
  ! close all .wfc files in order to prepare shuffled read
  !
  CLOSE(iuwfc, status = 'keep')
  ! never remove this barrier
  CALL mp_barrier(inter_pool_comm)
  !
  DO ik = 1, nks
     !
     IF (lsda) current_spin = isk(ik)
     elphmat(:,:,:) = czero
!DBSP
!     b = zero
!     c = zero
!     d = zero
!END
     !
     ! find index, and possibly pool, of k+q 
     ! the index nkq (nkq_abs) takes into account the even/odd ordering 
     ! of the nscf calc
     ! we also redefine the ikq points and the corresponding energies
     ! (we need to make sure that xk(:,ikq) is really k+q for the KB projectors
     ! below and also that the eigenvalues are taken correctly in ephwann)
     !
     CALL ktokpmq( xk(:,ik), xq, +1, ipool, nkq, nkq_abs )
     !
     !   we define xkq(:,ik) and etq(:,ik) for the current xq
     !
     IF (ALLOCATED(etq)) DEALLOCATE(etq)
     IF (.not. ALLOCATED(etq) ) ALLOCATE( etq(nbnd,nks) )
     etq(:,:) = zero
     !
     xkq(:,ik) = xk_all(:,nkq_abs)
     etq(:,ik) = et_all(:,nkq_abs) 
     !
     ipooltmp = my_pool_id + 1
     !
     ! in serial execution ipool is not used in the called subroutines, 
     ! in parallel ipooltmp is for k and ipool is for k+q
     !
     ! read unperturbed wavefunctions psi(k) and psi(k+q)
     !
     CALL readwfc( ipooltmp, ik, evc )
     CALL readwfc( ipool, nkq, evq )
     !
     ! Now we define the igk and igkq from the global igk_k_all
     ! 
     npw  = ngk_all(ik+lower_bnd-1)
     npwq = ngk_all(nkq_abs)
     ! 
     IF (ALLOCATED(igk)) DEALLOCATE(igk)
     IF (ALLOCATED(igkq)) DEALLOCATE(igkq)
     ALLOCATE( igk(npw)  )
     ALLOCATE( igkq(npwq) )
     ! 
     igk = igk_k_all(1:npw,ik+lower_bnd-1)
     igkq = igk_k_all(1:npwq,nkq_abs)
     !
     IF ( nks.gt.1 .AND. maxval(igkq(1:npwq)).gt.ngxx ) &
       CALL errore('elphel2_shuffle', 'ngxx too small', 1 )
     !
     ! ----------------------------------------------------------------
     ! Set the gauge for the eigenstates: unitary transform and phases
     ! ----------------------------------------------------------------
     !
     ! With this option, different compilers and different machines
     ! should always give the same wavefunctions.
     !
     CALL ktokpmq( xk(:,ik),  zero_vect, +1, ipool, nkk, nkk_abs )
     CALL ktokpmq( xkq(:,ik), zero_vect, +1, ipool, nkk, nkq_abs )
     !
     IF ( .not. ALLOCATED(umat) )  ALLOCATE( umat(nbnd,nbnd,nks) )
     IF ( .not. ALLOCATED(umatq) ) ALLOCATE( umatq(nbnd,nbnd,nks) )
     umat(:,:,ik)  = umat_all(:,:,nkk_abs)
     umatq(:,:,ik) = umat_all(:,:,nkq_abs)
     !
     ! the k-vector needed for the KB projectors
     xkqtmp = xkq(:,ik)
     !
     ! --------------------------------------------------
     !   Fourier translation of the G-sphere igkq
     ! --------------------------------------------------
     !
     !  Translate by G_0 the G-sphere where evq is defined, 
     !  none of the G-points are lost.
     !
     DO ig = 1, npwq
        imap = ng0vec * ( igkq(ig) - 1 ) + shift(ik+ik0)
        igkq_tmp(ig) = gmap(imap)
        !  the old matrix version... 
        !  igkq_tmp(ig) = gmap( igkq(ig), shift(ik+ik0) )
     ENDDO
     igkq = igkq_tmp
     !
     !  find k+q from k+q+G_0
     !  (this is needed in the calculation of the KB terms
     !  for nonlocal pseudos)
     !
     xkqtmp = xkq(:,ik) - g0vec_all_r(:,shift(ik+ik0))
     !
     ! ---------------------------------------------------------------------
     ! phase factor arising from fractional traslations
     ! ---------------------------------------------------------------------
     !
     !  u_{k+q+G_0} carries an additional factor e^{i G_0 v}
     !
     CALL fractrasl( npw,  igk,  evc, eigv(:,isym), cone )
     CALL fractrasl( npwq, igkq, evq, eigv(:,isym), cone )
     !
     ! ---------------------------------------------------------------------
     ! wave function rotation to generate matrix elements for the star of q
     ! ---------------------------------------------------------------------
     !
     ! ps. don't use npwx instead of npw, npwq since the unused elements
     ! may be large and blow up gmapsym (personal experience)
     !
     igk (1:npw ) = gmapsym( igk (1:npw ), isym )
     igkq(1:npwq) = gmapsym( igkq(1:npwq), isym )
     !
     ! In dvqpsi_us_only3 we need becp1 and alphap for the rotated wfs. 
     ! The other quantities (deeq and qq) do not depend on the wfs, in
     ! particular in the KB case (not ultrasoft), the deeq's are the
     ! unscreened coefficients, and the qq's are zero.
     !
     ! For the KB part, remember dV_NL[q_0] ~ |S^-1(k)+q_0> <S^-1(k)|
     ! the total momentum transfer must be q_0 and the rotation 
     ! tranforms k+Sq_0 into S^-1(k)+q_0, k into S^-1(k)
     ! [see Eqs. (A9),(A14) Baroni et al. RMP]
     ! 
     ! Since in QE a normal rotation s is defined as S^-1 we have here
     ! sxk = S(k).  
     !
     CALL rotate_cart( xk(:,ik), s(:,:,isym), sxk )
     !
     ! here we generate vkb on the igk() set and for k ...
     CALL init_us_2( npw, igk, sxk, vkb )
     !
     ! ... and we recompute the becp terms with the wfs (rotated through igk)
     !
     CALL calbec( npw, vkb, evc, becp1(ik) )
     !
     ! we also recompute the derivative of the becp terms with the (rotated) wfs
     !
     DO ipol = 1, 3
        aux2 = czero
        DO ibnd = 1, nbnd
           DO ig = 1, npw
              aux2(ig,ibnd) = evc(ig,ibnd) * tpiba * ci * & 
                             ( sxk(ipol) + g(ipol,igk(ig)) )
           END DO
           IF (noncolin) THEN
              DO ig = 1, npw
                 aux2(ig+npwx,ibnd) = evc(ig+npwx,ibnd) * tpiba * ci * &
                             ( sxk(ipol) + g(ipol,igk(ig)) )
              ENDDO
           ENDIF
        ENDDO
        CALL calbec( npw, vkb, aux2, alphap(ipol,ik) )
     ENDDO
     !
     ! now we generate vkb on the igkq() set because dvpsi is needed on that set
     ! we need S(k)+q_0 in the KB projector: total momentum transfer must be q_0
     !
     xkqtmp = sxk + xq0
     CALL init_us_2( npwq, igkq, xkqtmp, vkb )
     !
     ! --------------------------------------------------
     !   Calculation of the matrix element
     ! --------------------------------------------------
     !
     DO ipert = 1, npe
        !
        !  recalculate dvbare_q*psi_k 
        !  the call to dvqpsi_us3 differs from the old one to dvqpsi_us 
        !  only the xkqtmp passed. 
        !
        !  we have to use the first q in the star in the dvqpsi_us3 call below (xq0)
        !  
        mode = imode0 + ipert
        IF (timerev) THEN
          CALL dvqpsi_us3( ik, conjg(u(:,mode)), .false., xkqtmp, xq0 )
        ELSE
          CALL dvqpsi_us3( ik, u(:,mode), .false., xkqtmp, xq0 )
        ENDIF
!DBSP 
!        b = b+SUM((REAL(REAL(dvpsi(:,:))))**2)+SUM((REAL(AIMAG(dvpsi(:,:))))**2)
!END
        !
        !  calculate dvscf_q*psi_k
        !
        CALL start_clock('dvscf_q*psi_k')
        ! 
        aux3 = czero
        DO ibnd = lower_band, upper_band
          CALL invfft_wave(npw, igk, evc(:,ibnd), aux1)
          IF (timerev) THEN
            CALL apply_dpot(dffts%nnr, aux1, conjg(dvscfins(:,:,ipert)), current_spin)
          ELSE
            CALL apply_dpot(dffts%nnr, aux1, dvscfins(:,:,ipert), current_spin)
          ENDIF
          CALL fwfft_wave(npwq, igkq, aux3(:,ibnd), aux1)
        ENDDO
        dvpsi = dvpsi + aux3
        !
!DBSP
!        c = c+SUM((REAL(REAL(dvpsi(:,:))))**2)+SUM((REAL(AIMAG(dvpsi(:,:))))**2)
!END
        !
        CALL adddvscf2( ipert, ik )
!DBRM
!        d = c+SUM((REAL(REAL(dvpsi(:,:))))**2)+SUM((REAL(AIMAG(dvpsi(:,:))))**2)
!END
        !
        ! calculate elphmat(j,i)=<psi_{k+q,j}|dvscf_q*psi_{k,i}> for this pertur
        !
        ! 
        DO ibnd =lower_band, upper_band
           DO jbnd = 1, nbnd
              elphmat(jbnd,ibnd,ipert) = &
                    zdotc( npwq, evq(1,jbnd), 1, dvpsi(1,ibnd), 1 )
              IF (noncolin) &
                 elphmat(jbnd,ibnd,ipert) = elphmat(jbnd,ibnd,ipert) + &
                    zdotc( npwq, evq(npwx+1,jbnd), 1, dvpsi(npwx+1,ibnd), 1 )
           ENDDO
        ENDDO
     ENDDO
     !
     CALL mp_sum(elphmat, intra_pool_comm)
     CALL mp_sum(elphmat, inter_image_comm)
     !
!DBSP
!     IF (ik==2) THEN
!       write(*,*)'SUM dvpsi b ', b
!       write(*,*)'SUM dvpsi c ', c
!       write(*,*)'SUM dvpsi d ', d
!       write(*,*)'elphmat(:,:,:)**2', SUM((REAL(REAL(elphmat(:,:,:))))**2)+SUM((REAL(AIMAG(elphmat(:,:,:))))**2)
!     ENDIF
!END
     !
     !  Rotate elphmat with the gauge matrices (this should be equivalent 
     !  to calculate elphmat with the truely rotated eigenstates)
     ! 
     DO ipert = 1, npe
        !
        ! the two zgemm call perform the following ops:
        !  elphmat = umat(k+q)^\dagger * [ elphmat * umat(k) ]
        !
        CALL zgemm( 'n', 'n', nbnd, nbnd, nbnd, cone, elphmat(:,:,ipert), & 
                    nbnd, umat(:,:,ik), nbnd, czero, eptmp, nbnd )
        CALL zgemm( 'c', 'n', nbnd, nbnd, nbnd, cone, umatq(:,:,ik), & 
                    nbnd, eptmp, nbnd, czero, elphmat(:,:,ipert), nbnd )
        !
     ENDDO
     !
     !  save eph matrix elements into el_ph_mat
     !
     DO ipert = 1, npe
        DO jbnd = 1, nbnd
           DO ibnd = 1, nbnd
              el_ph_mat(ibnd,jbnd,ik,ipert+imode0) = elphmat(ibnd,jbnd,ipert)
           ENDDO
        ENDDO
     ENDDO
     !
  ENDDO
  !
  !  restore original configuration of files
  !
  CALL diropn(iuwfc, 'wfc', lrwfc, exst) 
  ! never remove this barrier - > insures that wfcs are restored to each pool before moving on
  CALL mp_barrier(world_comm)
  !
  DEALLOCATE(elphmat)
  DEALLOCATE(eptmp)
  DEALLOCATE(aux1) 
  DEALLOCATE(aux2)
  DEALLOCATE(aux3)
  DEALLOCATE(gmap)
  DEALLOCATE(shift)
  !
  END SUBROUTINE elphel2_shuffle
  !
  !------------------------------------------------------------
  SUBROUTINE fractrasl( npw, igk, evc, eigv1, eig0v )
  !------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE wvfct, ONLY : nbnd, npwx
  USE gvect, ONLY : ngm
  USE noncollin_module, ONLY : noncolin, npol
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: npw, igk(npwx)
  COMPLEX(kind=DP), INTENT(inout) :: evc(npwx*npol, nbnd)
  COMPLEX(kind=DP), INTENT(in) :: eigv1(ngm), eig0v
  !
  INTEGER :: ig
  !! Counter on G-vectors
  INTEGER :: ibnd
  !! Counter on bands
  ! 
  DO ibnd = 1, nbnd
     DO ig = 1, npw
        evc(ig,ibnd) = evc(ig,ibnd) * eigv1( igk(ig) ) * eig0v
        IF (noncolin) THEN
           evc(ig+npwx,ibnd) = evc(ig+npwx,ibnd) * eigv1( igk(ig) ) * eig0v
        ENDIF
     ENDDO
  ENDDO
  !
  END SUBROUTINE fractrasl
  !
  !------------------------------------------------------------
  SUBROUTINE rotate_cart( x, s, sx )
  !------------------------------------------------------------
  !
  ! a simple symmetry operation in cartesian coordinates 
  ! ( s is integer and in crystal coord!)
  !
  USE kinds, ONLY : DP
  USE cell_base, ONLY : at, bg
  !
  IMPLICIT NONE
  !
  REAL(kind=DP), INTENT(in) :: x(3)
  !! Input x
  INTEGER, INTENT(in) :: s(3,3)
  !! Symmetry matrix
  REAL(kind=DP), INTENT(out) :: sx(3)
  !! Output rotated x
  !
  REAL(DP) :: xcrys(3)
  !! x in cartesian coords
  INTEGER :: i
  !
  xcrys = x
  CALL cryst_to_cart(1, xcrys, at, -1)
  DO i = 1, 3
     sx(i) = dble(s(i,1)) * xcrys(1) &
           + dble(s(i,2)) * xcrys(2) &
           + dble(s(i,3)) * xcrys(3)
  ENDDO
  CALL cryst_to_cart(1, sx, bg, +1)
  !
  END SUBROUTINE rotate_cart
  !------------------------------------------------------------
