!
! Copyright (C) 2013 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Nov 2013, D. Ceresoli: Z2 calculation according Soluyanov and Vanderbilt [1]
!
! This routine is largerly inspired from bp_c_phase.f90, and it includes
! routines from the Z2PACK code by A. Soluyanov [2]
!
!##############################################################################!
!#                                                                            #!
!#                                                                            #!
!#   This is the main one of a set of Fortran 90 files designed to compute    #!
!#   the Z2 invariant without inversion symmetry.                             #!
!#                                                                            #!
!#                                                                            #!
!#   BRIEF SUMMARY OF THE METHODOLOGY                                         #!
!#   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                         #!
!#   The standard procedure would be for the user to first perform a          #!
!#   self-consistent (sc) calculation to obtain a converged charge density.   #!
!#   With well-converged sc charge density, the user would then run one       #!
!#   or more non-self consistent (or "band structure") calculations,          #!
!#   using the same main code, but with a flag to ask for the Z2 calculation. #!
!#   Each such run would calculate the Z2 invariant of TRS face of the        #!
!#   Brillouin zone.                                                          #!
!#                                                                            #!
!#   Accurate calculation of the Z2 invariant requires overlaps between       #!
!#   wavefunctions along fairly dense lines (or "strings") in k-space in the  #!
!#   direction of the primitive G-vector parallel to the BZ surface.          #!
!#                                                                            #!
!#                                                                            #!
!#   FUNCTIONALITY/COMPATIBILITY                                              #!
!#   ~~~~~~~~~~~~~~~~~~~~~~~~~~~                                              #!
!#   * Spin-polarized systems supported.                                      #!
!#                                                                            #!
!#   * (Ultrasoft and) norm-conserving pseudopotentials supported             #!
!#                                                                            #!
!#   * Calculation must be non collinear and including spin orbit             #!
!#                                                                            #!
!#                                                                            #!
!#   NEW INPUT PARAMETERS                                                     #!
!#   ~~~~~~~~~~~~~~~~~~~~                                                     #!
!#   * lcalc_z2 (.TRUE. or .FALSE.)                                           #!
!#     Tells PWSCF that a Z2 invariant calcultion is desired.                 #!
!#                                                                            #!
!#   * gdir (1, 2, or 3)                                                      #!
!#     Specifies the direction of the k-point strings in reciprocal space.    #!
!#     '1' refers to the first reciprocal lattice vector, '2' to the          #!
!#     second, and '3' to the third.                                          #!
!#                                                                            #!
!#   * nppstr (integer)                                                       #!
!#     Specifies the number of k-points to be calculated along each           #!
!#     symmetry-reduced string.                                               #!
!#                                                                            #!
!#   * z2_m_threshodl (float, default 0.8)                                    #!
!#     Threshold for SVD decomposition. If SV < z2_m_threshold, the number    #!
!#     k-points in the string should be increased                             #!
!#                                                                            #!
!#   * z2_z_threshold (float, default 0.05)                                   #!
!#     Threshold for detecting Wannier centers jumps. In case of warning,     #!
!#     you must increase the number of strings                                #!
!#                                                                            #!
!#                                                                            #!
!#   EXPLANATION OF K-POINT MESH                                              #!
!#   ~~~~~~~~~~~~~~~~~~~~~~~~~~~                                              #!
!#   See directory PW/examples/example12 for a thourough explanation.         #!
!#                                                                            #!
!#                                                                            #!
!#   BIBLIOGRAPHY                                                             #!
!#   ~~~~~~~~~~~~                                                             #!
!#   [1] A. A. Soluyanov and D. Vanderbilt, "Computing topological            #!
!#       invariants without inversion symmetry",                              #!
!#       Phys Rev B 83, 235401 (2011).                                        #!
!#                                                                            #!
!#   [2] http://www.physics.rutgers.edu/z2pack/                               #!
!#                                                                            #!
!#                                                                            #!
!##############################################################################!


!======================================================================
SUBROUTINE c_phase_z2
!======================================================================
   !
   !   Geometric phase calculation along a strip of nppstr k-points
   !   averaged over a 1D grid of nkort k-points ortogonal to nppstr 
   !  --- Make use of the module with common information ---
   !
   USE kinds,                ONLY : DP
   USE io_global,            ONLY : stdout
   USE io_files,             ONLY : iunwfc, nwordwfc
   USE buffers,              ONLY : get_buffer
   USE cell_base,            ONLY : at, tpiba, omega, tpiba2
   USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau
   USE constants,            ONLY : pi, tpi
   USE gvect,                ONLY : ngm, g, gcutm, ngm_g, ig_l2g
   USE fft_base,             ONLY : dfftp
   USE uspp,                 ONLY : nkb, vkb, okvan
   USE uspp_param,           ONLY : upf, lmaxq, nbetam, nh, nhm
   USE klist,                ONLY : nks, xk, wk
   USE wvfct,                ONLY : npwx, nbnd, ecutwfc
   USE wavefunctions_module, ONLY : evc
   USE bp,                   ONLY : gdir, nppstr, mapgm_global, &
                                    z2_m_threshold, z2_z_threshold
   USE becmod,               ONLY : calbec, bec_type, allocate_bec_type, &
                                    deallocate_bec_type
   USE noncollin_module,     ONLY : noncolin, npol, nspin_lsda
   USE spin_orb,             ONLY : lspinorb
   USE mp_bands,             ONLY : intra_bgrp_comm, nproc_bgrp
   USE mp,                   ONLY : mp_sum
   USE control_flags,        ONLY : iverbosity

!  --- Avoid implicit definitions ---
   IMPLICIT NONE

!  --- Internal definitions ---
   INTEGER :: i
   INTEGER :: igk1(npwx)
   INTEGER :: igk0(npwx)
   INTEGER :: ig
   INTEGER :: is
   INTEGER :: istring
   INTEGER :: iv
   INTEGER :: j
   INTEGER :: jkb
   INTEGER :: jkb_bp
   INTEGER :: jkb1
   INTEGER :: jv
   INTEGER :: kindex
   INTEGER :: kort
   INTEGER :: kpar
   INTEGER :: kpoint
   INTEGER :: kstart
   INTEGER :: mb
   INTEGER :: mk1
   INTEGER :: mk2
   INTEGER :: mk3
   INTEGER , ALLOCATABLE :: ln(:,:,:)
   INTEGER :: n1
   INTEGER :: n2
   INTEGER :: n3
   INTEGER :: na
   INTEGER :: nb
   INTEGER :: ng
   INTEGER :: nhjkb
   INTEGER :: nhjkbm
   INTEGER :: nkbtona(nkb)
   INTEGER :: nkbtonh(nkb)
   INTEGER :: nkort
   INTEGER :: np
   INTEGER :: npw1
   INTEGER :: npw0
   INTEGER :: nstring
   INTEGER :: nbnd_occ
   INTEGER :: nt
   INTEGER, ALLOCATABLE :: map_g(:)
   LOGICAL :: l_para
   LOGICAL, ALLOCATABLE :: l_cal(:) ! flag for occupied/empty states
   REAL(DP) :: dk(3)
   REAL(DP) :: dkmod
   REAL(DP), parameter :: eps = 1d-6
   REAL(DP) :: fac
   REAL(DP) :: g2kin_bp(npwx)
   REAL(DP) :: gpar(3)
   REAL(DP) :: gtr(3)
   REAL(DP) :: gvec
   REAL(DP) :: qrad_dk(nbetam,nbetam,lmaxq,ntyp)
   REAL(DP) :: ylm_dk(lmaxq*lmaxq)
   COMPLEX(DP), ALLOCATABLE :: aux(:)
   COMPLEX(DP), ALLOCATABLE :: aux_g(:)
   COMPLEX(DP), ALLOCATABLE :: aux0(:)
   TYPE (bec_type) :: becp0
   TYPE (bec_type) :: becp_bp
   COMPLEX(DP) :: mat(nbnd,nbnd)
   COMPLEX(DP) :: pref
   COMPLEX(DP), ALLOCATABLE :: psi(:,:)
   COMPLEX(DP), ALLOCATABLE :: q_dk_so(:,:,:,:)
   COMPLEX(DP) :: q_dk(nhm,nhm,ntyp)
   COMPLEX(DP) :: struc(nat)
   COMPLEX(DP), external :: zdotc

!  -------------------------------------------------------------------------   !
!                               Z2 variables
!  -------------------------------------------------------------------------   !
   complex(dp), allocatable :: UU(:,:), VT(:,:), lambda(:,:), eig(:)
   real(dp), allocatable :: SV(:), numbers(:,:)
   real(dp) :: counter

!  -------------------------------------------------------------------------   !
!                               INITIALIZATIONS
!  -------------------------------------------------------------------------   !
   ALLOCATE (psi(npwx*npol,nbnd))
   ALLOCATE (aux(ngm*npol))
   ALLOCATE (aux0(ngm*npol))
   IF (okvan) THEN
      CALL allocate_bec_type ( nkb, nbnd, becp0 )
      CALL allocate_bec_type ( nkb, nbnd, becp_bp )
      IF (lspinorb) ALLOCATE(q_dk_so(nhm,nhm,4,ntyp))
   END IF

   l_para= (nproc_bgrp > 1 .AND. gdir /= 3)
   IF (l_para) THEN
      ALLOCATE ( aux_g(ngm_g*npol) )
   ELSE
      ALLOCATE ( map_g(ngm) )
   ENDIF

!  --- Write header ---
   WRITE( stdout,"(/,/,/,15X,50('='))")
   WRITE( stdout,"(27X,'Z2 INVARIANT CALCULATION')")
   WRITE( stdout,"(25X,'!!! NOT THOROUGHLY TESTED !!!')")
   WRITE( stdout,"(15X,50('='),/)")
   if (.not. lspinorb) call errore('z2_c_phase', 'Z2 needs spin orbit', 1)
   if (nspin_lsda /= 1) call errore('z2_c_phase', 'internal error: nspin_lsda=', nspin_lsda)

!  --- Recalculate FFT correspondence (see ggen.f90) ---
   ALLOCATE (ln (-dfftp%nr1:dfftp%nr1, -dfftp%nr2:dfftp%nr2, -dfftp%nr3:dfftp%nr3) )
   DO ng=1,ngm
      mk1=nint(g(1,ng)*at(1,1)+g(2,ng)*at(2,1)+g(3,ng)*at(3,1))
      mk2=nint(g(1,ng)*at(1,2)+g(2,ng)*at(2,2)+g(3,ng)*at(3,2))
      mk3=nint(g(1,ng)*at(1,3)+g(2,ng)*at(2,3)+g(3,ng)*at(3,3))
      ln(mk1,mk2,mk3) = ng
   END DO

   if(okvan) then
!  --- Initialize arrays ---
      jkb_bp=0
      DO nt=1,ntyp
         DO na=1,nat
            IF (ityp(na).eq.nt) THEN
               DO i=1, nh(nt)
                  jkb_bp=jkb_bp+1
                  nkbtona(jkb_bp) = na
                  nkbtonh(jkb_bp) = i
               END DO
            END IF
         END DO
      END DO
   endif
!  --- Get the number of strings ---
   nstring=nks/nppstr
   nkort=nstring   !/nspin_lsda

!  -------------------------------------------------------------------------   !
!           electronic polarization: set values for k-points strings           !
!  -------------------------------------------------------------------------   !

!  --- Find vector along strings ---
   gpar(1)=xk(1,nppstr)-xk(1,1)
   gpar(2)=xk(2,nppstr)-xk(2,1)
   gpar(3)=xk(3,nppstr)-xk(3,1)
   gvec=dsqrt(gpar(1)**2+gpar(2)**2+gpar(3)**2)*tpiba

!  --- Find vector between consecutive points in strings ---
   dk(1)=xk(1,2)-xk(1,1)
   dk(2)=xk(2,2)-xk(2,1) 
   dk(3)=xk(3,2)-xk(3,1)
   dkmod=SQRT(dk(1)**2+dk(2)**2+dk(3)**2)*tpiba
   IF (ABS(dkmod-gvec/(nppstr-1)) > eps) & 
     CALL errore('c_phase_z2','Wrong k-strings?',1)

!  --- Check that k-points form strings ---
   DO i=1,nspin_lsda*nkort
      DO j=2,nppstr
         kindex=j+(i-1)*nppstr
         IF (ABS(xk(1,kindex)-xk(1,kindex-1)-dk(1)) > eps) &
            CALL errore('c_phase_z2','Wrong k-strings?',1)
         IF (ABS(xk(2,kindex)-xk(2,kindex-1)-dk(2)) > eps) &
            CALL errore('c_phase_z2','Wrong k-strings?',1)
         IF (ABS(xk(3,kindex)-xk(3,kindex-1)-dk(3)) > eps) &
            CALL errore('c_phase_z2','Wrong k-strings?',1)
         IF (ABS(wk(kindex)-wk(kindex-1)) > eps) &
            CALL errore('c_phase_z2','Wrong k-strings weights?',1)
      END DO
   END DO

   ! TODO: check k-points are on the edge of BZ


!  -------------------------------------------------------------------------   !
!                  electronic polarization: structure factor                   !
!  -------------------------------------------------------------------------   !

!  --- Calculate structure factor e^{-i dk*R} ---
   DO na=1,nat
      fac=(dk(1)*tau(1,na)+dk(2)*tau(2,na)+dk(3)*tau(3,na))*tpi 
      struc(na)=CMPLX(cos(fac),-sin(fac),kind=DP)
   END DO

!  -------------------------------------------------------------------------   !
!                     electronic polarization: form factor                     !
!  -------------------------------------------------------------------------   !
   if(okvan) then
!  --- Calculate Bessel transform of Q_ij(|r|) at dk [Q_ij^L(|r|)] ---
      CALL calc_btq(dkmod,qrad_dk,0)

!  --- Calculate the q-space real spherical harmonics at dk [Y_LM] --- 
      dkmod=dk(1)**2+dk(2)**2+dk(3)**2
      CALL ylmr2(lmaxq*lmaxq, 1, dk, dkmod, ylm_dk)
!  --- Form factor: 4 pi sum_LM c_ij^LM Y_LM(Omega) Q_ij^L(|r|) ---
      q_dk = (0.d0, 0.d0)
      DO np =1, ntyp
         if( upf(np)%tvanp ) then
            DO iv = 1, nh(np)
               DO jv = iv, nh(np)
                  call qvan3(iv,jv,np,pref,ylm_dk,qrad_dk)
                  q_dk(iv,jv,np) = omega*pref
                  q_dk(jv,iv,np) = omega*pref
               ENDDO
            ENDDO
         endif
      ENDDO
      IF (lspinorb) CALL transform_qq_so(q_dk,q_dk_so)
   endif

!  -------------------------------------------------------------------------   !
!                   electronic polarization: strings phases                    !
!  -------------------------------------------------------------------------   !
   kpoint=0
   allocate (SV(nbnd), UU(nbnd,nbnd), VT(nbnd,nbnd), lambda(nbnd,nbnd), eig(nbnd))
   allocate (numbers(nkort, nbnd+2))
   numbers = 0.d0

   allocate (l_cal(nbnd))  ! l_cal(n) = .true./.false. if n-th state is occupied/empty
   nbnd_occ = nbnd
   l_cal(1:nbnd) = .true.
   !call weights()

!  --- Start loop over spin ---
   DO is=1,nspin_lsda

!     --- Start loop over orthogonal k-points ---
      DO kort=1,nkort
         write(stdout,*)
         write(stdout,*)
         write(stdout,'(5X,''========== kort='',I4,'' =========='')') kort

!        --- Index for this string ---
         istring=kort+(is-1)*nkort

!        --- Initialize Lambda matrix ---
         lambda = (0.d0, 0.d0)
         do nb = 1, nbnd
            lambda(nb,nb) = (1.d0, 0.d0)
         enddo

!        --- Start loop over parallel k-points ---
         DO kpar = 1,nppstr

!           --- Set index of k-point ---
            kpoint = kpoint + 1

!           --- Calculate dot products between wavefunctions and betas ---
            IF (kpar /= 1) THEN

!              --- Dot wavefunctions and betas for PREVIOUS k-point ---
               CALL gk_sort(xk(1,kpoint-1),ngm,g,ecutwfc/tpiba2, &
                            npw0,igk0,g2kin_bp) 
               CALL get_buffer (psi,nwordwfc,iunwfc,kpoint-1)
               if (okvan) then
                  CALL init_us_2 (npw0,igk0,xk(1,kpoint-1),vkb)
                  CALL calbec (npw0, vkb, psi, becp0)
               endif
!              --- Dot wavefunctions and betas for CURRENT k-point ---
               IF (kpar /= nppstr) THEN
                  CALL gk_sort(xk(1,kpoint),ngm,g,ecutwfc/tpiba2, &
                               npw1,igk1,g2kin_bp)        
                  CALL get_buffer(evc,nwordwfc,iunwfc,kpoint)
                  if (okvan) then
                     CALL init_us_2 (npw1,igk1,xk(1,kpoint),vkb)
                     CALL calbec (npw1, vkb, evc, becp_bp)
                  endif
               ELSE
                  kstart = kpoint-nppstr+1
                  CALL gk_sort(xk(1,kstart),ngm,g,ecutwfc/tpiba2, &
                               npw1,igk1,g2kin_bp)  
                  CALL get_buffer(evc,nwordwfc,iunwfc,kstart)
                  if (okvan) then
                     CALL init_us_2 (npw1,igk1,xk(1,kstart),vkb)
                     CALL calbec(npw1, vkb, evc, becp_bp)
                  endif
               ENDIF

               IF (kpar == nppstr .AND. .NOT. l_para) THEN
                  map_g(:) = 0
                  DO ig=1,npw1
!                          --- If k'=k+G_o, the relation psi_k+G_o (G-G_o) ---
!                          --- = psi_k(G) is used, gpar=G_o, gtr = G-G_o ---

                     gtr(1)=g(1,igk1(ig)) - gpar(1)
                     gtr(2)=g(2,igk1(ig)) - gpar(2)
                     gtr(3)=g(3,igk1(ig)) - gpar(3)
!                          --- Find crystal coordinates of gtr, n1,n2,n3 ---
!                          --- and the position ng in the ngm array ---

                     IF (gtr(1)**2+gtr(2)**2+gtr(3)**2 <= gcutm) THEN
                        n1=NINT(gtr(1)*at(1,1)+gtr(2)*at(2,1) &
                             +gtr(3)*at(3,1))
                        n2=NINT(gtr(1)*at(1,2)+gtr(2)*at(2,2) &
                             +gtr(3)*at(3,2))
                        n3=NINT(gtr(1)*at(1,3)+gtr(2)*at(2,3) &
                             +gtr(3)*at(3,3))
                        ng=ln(n1,n2,n3)

                        IF ( (ABS(g(1,ng)-gtr(1)) > eps) .OR. &
                             (ABS(g(2,ng)-gtr(2)) > eps) .OR. &
                             (ABS(g(3,ng)-gtr(3)) > eps) ) THEN
                           WRITE(6,*) ' error: translated G=', &
                                gtr(1),gtr(2),gtr(3), &
                                &     ' with crystal coordinates',n1,n2,n3, &
                                &     ' corresponds to ng=',ng,' but G(ng)=', &
                                &     g(1,ng),g(2,ng),g(3,ng)
                           WRITE(6,*) ' probably because G_par is NOT', &
                                &    ' a reciprocal lattice vector '
                           WRITE(6,*) ' Possible choices as smallest ', &
                                ' G_par:'
                           DO i=1,50
                              WRITE(6,*) ' i=',i,'   G=', &
                                   g(1,i),g(2,i),g(3,i)
                           ENDDO
                           CALL errore('c_phase_z2','wrong g',1)
                        ENDIF
                     ELSE
                        WRITE(6,*) ' |gtr| > gcutm  for gtr=', &
                             gtr(1),gtr(2),gtr(3)
                        CALL errore('c_phase_z2','wrong gtr',1)
                     END IF
                     map_g(ig)=ng
                  END DO
               END IF

!              --- Matrix elements calculation ---

               mat(:,:) = (0.d0, 0.d0)
               DO mb=1,nbnd
                  IF ( .NOT. l_cal(mb) ) THEN
                      mat(mb,mb)=(1.d0, 0.d0)
                  ELSE
                     aux(:) = (0.d0, 0.d0)
                     IF (kpar /= nppstr) THEN
                        DO ig=1,npw1
                           aux(igk1(ig))=evc(ig,mb)
                        ENDDO
                        IF (noncolin) THEN
                           DO ig=1,npw1
                              aux(igk1(ig)+ngm)=evc(ig+npwx,mb)
                           ENDDO
                        ENDIF
                     ELSEIF (.NOT. l_para) THEN
                        DO ig=1,npw1
                           aux(map_g(ig))=evc(ig,mb)
                        ENDDO
                        IF (noncolin) THEN
                           DO ig=1,npw1
                              aux(map_g(ig)+ngm)=evc(ig+npwx,mb)
                           ENDDO
                        ENDIF
                     ELSE
!
!   In this case this processor might not have the G-G_0
!
                        aux_g=(0.d0,0.d0)
                        DO ig=1,npw1
                           aux_g(mapgm_global(ig_l2g(igk1(ig)),gdir)) &
                                                =evc(ig,mb)
                        ENDDO
                        IF (noncolin) THEN
                           DO ig=1,npw1
                              aux_g(mapgm_global(ig_l2g(igk1(ig)),gdir) &
                                                + ngm_g) =evc(ig+npwx,mb)
                           ENDDO
                        ENDIF
                        CALL mp_sum(aux_g(:), intra_bgrp_comm )
                        DO ig=1,ngm
                           aux(ig) = aux_g(ig_l2g(ig))
                        ENDDO
                        IF (noncolin) THEN
                           DO ig=1,ngm
                              aux(ig+ngm) = aux_g(ig_l2g(ig)+ngm_g)
                           ENDDO
                        ENDIF
                     ENDIF
!
                     DO nb=1,nbnd
                        IF ( l_cal(nb) ) THEN
                           aux0(:)= (0.d0, 0.d0)
                           DO ig=1,npw0
                              aux0(igk0(ig))=psi(ig,nb)
                           END DO
                           IF (noncolin) THEN
                              DO ig=1,npw0
                                aux0(igk0(ig)+ngm)=psi(ig+npwx,nb)
                              END DO
                           ENDIF
                           mat(nb,mb) = zdotc (ngm*npol,aux0,1,aux,1)
                        END IF
                     END DO
                  END IF
               END DO
               !
               call mp_sum( mat, intra_bgrp_comm )
               !
               DO nb=1,nbnd
                  DO mb=1,nbnd
!                    --- Calculate the augmented part: ij=KB projectors, ---
!                    --- R=atom index: SUM_{ijR} q(ijR) <u_nk|beta_iR>   ---
!                    --- <beta_jR|u_mk'> e^i(k-k')*R =                   ---
!                    --- also <u_nk|beta_iR>=<psi_nk|beta_iR> = becp^*   ---
                     IF ( l_cal(nb) .AND. l_cal(mb) ) THEN
                        if (okvan) then
                           pref = (0.d0,0.d0)
                           DO jkb=1,nkb
                              nhjkb = nkbtonh(jkb)
                              na = nkbtona(jkb)
                              np = ityp(na)
                              nhjkbm = nh(np)
                              jkb1 = jkb - nhjkb
                              DO j = 1,nhjkbm
                                 IF (noncolin) THEN
                                    IF (lspinorb) THEN
                                       pref = pref+(CONJG(becp0%nc(jkb,1,nb))* &
                                                  becp_bp%nc(jkb1+j,1,mb)  &
                                            *q_dk_so(nhjkb,j,1,np)   &
                                            +CONJG(becp0%nc(jkb,1,nb))* &
                                                   becp_bp%nc(jkb1+j,2,mb)  &
                                            *q_dk_so(nhjkb,j,2,np) &
                                            +CONJG(becp0%nc(jkb,2,nb))* &
                                                   becp_bp%nc(jkb1+j,1,mb)  &
                                            *q_dk_so(nhjkb,j,3,np) &
                                            +CONJG(becp0%nc(jkb,2,nb))* &
                                                   becp_bp%nc(jkb1+j,2,mb)   &
                                            *q_dk_so(nhjkb,j,4,np))*struc(na)
                                    ELSE
                                       pref = pref+(CONJG(becp0%nc(jkb,1,nb))* &
                                                   becp_bp%nc(jkb1+j,1,mb) + &
                                             CONJG(becp0%nc(jkb,2,nb))* &
                                                   becp_bp%nc(jkb1+j,2,mb))  &
                                             *q_dk(nhjkb,j,np)*struc(na)
                                    END IF
                                 ELSE
                                    pref = pref+CONJG(becp0%k(jkb,nb))* &
                                           becp_bp%k(jkb1+j,mb) &
                                      *q_dk(nhjkb,j,np)*struc(na)
                                 END IF
                              ENDDO
                           ENDDO
                           mat(nb,mb) = mat(nb,mb) + pref
                        endif
                     endif
                  ENDDO
               ENDDO

               CALL Z2PACK_Z2_MAIN1
!           --- End of dot products between wavefunctions and betas ---
            ENDIF

!        --- End loop over parallel k-points ---
         END DO  ! kpar

         CALL Z2PACK_Z2_MAIN2

!     --- End loop over orthogonal k-points ---
      END DO  ! kort

!  --- End loop over spin ---
   END DO

   CALL Z2PACK_Z2_FINAL(counter)

!  -------------------------------------------------------------------------   !
!                           write output information                           !
!  -------------------------------------------------------------------------   !

!  --- Information about the k-points string used ---
   WRITE( stdout,"(/,21X,'K-POINTS STRINGS USED IN CALCULATIONS')")
   WRITE( stdout,"(21X,37('~'),/)")
   WRITE( stdout,"(7X,'G-vector along string (2 pi/a):',3F9.5)") &
           gpar(1),gpar(2),gpar(3)
   WRITE( stdout,"(7X,'Modulus of the vector (1/bohr):',F9.5)") &
           gvec
   WRITE( stdout,"(7X,'Number of k-points per string:',I4)") nppstr
   WRITE( stdout,"(7X,'Number of different strings  :',I4)") nkort
   WRITE( stdout,"(7X,'Threshold for SVD            :',F9.5)") z2_m_threshold
   WRITE( stdout,"(7X,'Threshold for z2_main        :',F9.5)") z2_z_threshold

   WRITE( stdout,* )
   WRITE( stdout,"(7X,'The parity for the current TRS surface is:',F10.3)") counter
   if (counter > 0.99d0) then
      WRITE( stdout,"(7X,'The Z2 invariant for the current TRS surface is:  0')")
   elseif (counter < -0.99d0) then
      WRITE( stdout,"(7X,'The Z2 invariant for the current TRS surface is: -1')")
   else
      WRITE(stdout,"(7X,'Something is wrong with parity')")
   endif
   WRITE(stdout,*)
   WRITE(stdout,"(7X,'[A. A. Soluyanov and D. Vanderbilt, Phys Rev B 83, 235401 (2011)]')")         

!  --- End of information relative to polarization calculation ---
   WRITE( stdout,"(/,/,15X,50('=')/,/)")

!  -------------------------------------------------------------------------   !
!                                  finalization                                !
!  -------------------------------------------------------------------------   !

!  --- Free memory ---
   DEALLOCATE(ln)
   DEALLOCATE(aux)
   DEALLOCATE(aux0)
   DEALLOCATE(psi)
   IF (l_para) THEN
      DEALLOCATE ( aux_g )
   ELSE
      DEALLOCATE ( map_g )
   ENDIF
   IF (okvan) THEN
      CALL deallocate_bec_type ( becp0 )
      CALL deallocate_bec_type ( becp_bp )
      IF (lspinorb) DEALLOCATE(q_dk_so)
   END IF


CONTAINS
!  -------------------------------------------------------------------------   !
!  THE FOLLOWING CODE IS TAKEN FROM z2_main.f90 BY ALEXEY SOLUYANOV
!  -------------------------------------------------------------------------   !
   SUBROUTINE Z2PACK_Z2_MAIN1
   implicit none
   complex(dp), allocatable :: work(:)
   real(dp), allocatable :: rwork(:)
   integer :: lwork, info

!  --- Calculate SVD of M matrix: M = UU * SV * VV^* ---
   allocate (work(1))
   lwork = -1
   call ZGESVD('A','A',nbnd,nbnd,aux,nbnd,SV,UU,nbnd,VT,nbnd,work,lwork,rwork,info)
   lwork = int(dble(work(1)))
   deallocate (work)

   allocate (work(lwork), rwork(5*nbnd))
   call ZGESVD('A','A',nbnd,nbnd,mat,nbnd,SV,UU,nbnd,VT,nbnd,work,lwork,rwork,info)
   deallocate (work, rwork)
   if (info > 0) call errore('c_phase_z2','error in SVD factorization, info > 0', info)
   if (info < 0) call errore('c_phase_z2','error in SVD factorization, info < 0', -info)

   if (iverbosity > 1) then
      write(stdout,'(5X,''kort,kpar='',2I4,10X,''k,kp='',2I4)') kort, kpar, kpoint-1, kpoint
      write(stdout,'(5X,''sigmas:'')')
      write(stdout,'(''  '',8F9.4)') (SV(nb), nb=1,nbnd)
      write(stdout,*)
   endif

   ! test for Mmn threshold
   do nb = 1, nbnd
      if (SV(nb) < z2_m_threshold) &
        write(stdout,'(5X,''warning: k-point string is too coarse at kpar='',I3)') kpar
   enddo

   mat = matmul(UU, VT)
   lambda = matmul(lambda, mat)

   END SUBROUTINE Z2PACK_Z2_MAIN1



!  -------------------------------------------------------------------------   !
!  THE FOLLOWING CODE IS TAKEN FROM z2_main.f90 BY ALEXEY SOLUYANOV
!  -------------------------------------------------------------------------   !
   SUBROUTINE Z2PACK_Z2_MAIN2
   implicit none
   complex(dp), allocatable :: work(:)
   real(dp), allocatable :: rwork(:), zz(:), gaps(:)
   integer, allocatable :: ind(:)
   integer :: lwork, info, i
   real(dp) :: point, max_gap
   integer :: kk

   ! diagonalize lamdba matrix
   lwork = -1
   allocate(work(1))
   call ZGEEV('N','N',nbnd,aux,nbnd,eig,VT,nbnd,UU,nbnd,work,lwork,rwork,info)
   lwork = int(dble(work(1)))
   deallocate (work)

   allocate (work(lwork), rwork(2*nbnd))
   call ZGEEV('N','N',nbnd,lambda,nbnd,eig,VT,nbnd,UU,nbnd,work,lwork,rwork,info)
   deallocate (work, rwork)
   if (info > 0) call errore('c_phase_z2','error in diagonalization, info > 0', info)
   if (info < 0) call errore('c_phase_z2','error in diagonalization, info < 0', -info)

   if (iverbosity > 1) then
      write(stdout,'(5X,''eigenvalues (real part):'')')
      write(stdout,'(''  '',8F9.4)') (dble(eig(nb)), nb=1,nbnd)
      write(stdout,'(5X,''eigenvalues (imag part):'')')
      write(stdout,'(''  '',8F9.4)') (dimag(eig(nb)), nb=1,nbnd)
   endif

   allocate(zz(nbnd), gaps(nbnd), ind(nbnd))
   do nb = 1, nbnd
       if (dreal(cdlog(eig(nb))) > 1d-8) &
          write(stdout,'(5X,''warning: eigenvalue'',I4,'' has real part'')') nb
       zz(nb) = dimag(cdlog(eig(nb))) / (2.d0*PI)
   enddo

   ! putting eigenvalues within the [-0.5,0.5) window
   do nb = 1, nbnd
      i = int(zz(nb))
      zz(nb) = zz(nb) - i
      if (zz(nb) <= -0.5d0) zz(nb) = zz(nb) + 1.d0
      if (zz(nb) > 0.5d0) zz(nb) = zz(nb) - 1.d0
   enddo
 
   call hpsort(nbnd, zz, ind)
   write(stdout,'(5X,''Wannier centers:'')')
   write(stdout,'(''  '',8F9.4)') (zz(nb), nb=1,nbnd)

   ! gaps
   do nb = 1, nbnd
      if (nb < nbnd) then
         gaps(nb) = zz(nb+1)-zz(nb)
      else
         gaps(nb) = zz(1)+1.d0-zz(nb)
      endif
   enddo

   ! find largest gap
   max_gap = maxval(gaps)
   do nb = 1, nbnd
      if (gaps(nb) == max_gap) kk = nb
   enddo
   if (kk /= nbnd) then
       point = (zz(kk+1) + zz(kk))/2.d0
   else
       point = (zz(1) + 1.d0 + zz(kk))/2.d0
   endif

   write(stdout,'(5X,''largest gap midpoint:'',F9.4)') point
   numbers(kort,1:nbnd) = zz(1:nbnd)
   numbers(kort,nbnd+1) = zz(1) + 1.d0
   numbers(kort,nbnd+2) = point

   deallocate(zz, ind, gaps)

   END SUBROUTINE Z2PACK_Z2_MAIN2


!  -------------------------------------------------------------------------   !
!  THE FOLLOWING CODE IS TAKEN FROM z2_final.f90 BY ALEXEY SOLUYANOV
!  -------------------------------------------------------------------------   !
   SUBROUTINE Z2PACK_Z2_FINAL(counter)
   implicit none
   real(dp), intent(out) :: counter
   real(dp) :: mid1, mid2, theta1, theta2, p, den, center

   counter = 0
   do kort = 2, nkort
      mid1 = numbers(kort-1, nbnd+2)
      mid2 = numbers(kort, nbnd+2)
      theta1 = min(mid1, mid2)
      theta2 = max(mid1, mid2)
      do nb = 1, nbnd
         if (dabs(numbers(kort,nb)-mid1) < z2_z_threshold) &
            write(stdout,'(5X,''warning: string spacing is too coarse'')')
         if (numbers(kort,nb) > theta1 .and. numbers(kort,nb) < theta2) &
            counter = counter + 1
      enddo          
   enddo

   counter = 1
   do kort = 2, nkort
      mid1 = 2.d0*pi*numbers(kort-1, nbnd+2)
      mid2 = 2.d0*pi*numbers(kort, nbnd+2)
      do nb = 1, nbnd
         center = 2.d0*pi*numbers(kort,nb)
         den = dsin(mid2-mid1) + dsin(center-mid2) + dsin(mid1-center)
         if (dabs(den) < eps) then
            p = 1.d0
         else
            p = den / dabs(den)
         endif
         counter = counter * p
       enddo
    enddo
   END SUBROUTINE Z2PACK_Z2_FINAL


END SUBROUTINE c_phase_z2

