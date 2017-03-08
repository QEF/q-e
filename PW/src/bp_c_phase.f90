!
! Copyright (C) 2004 Vanderbilt's group at Rutgers University, NJ
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! April 2012, A. Dal Corso: parallelization for gdir /= 3 imported
!                           from c_phase_field.f90
! May 2012, A. Dal Corso: Noncollinear/spin-orbit case allowed (experimental).
!
!##############################################################################!
!#                                                                            #!
!#                                                                            #!
!#   This is the main one of a set of Fortran 90 files designed to compute    #!
!#   the electrical polarization in a crystaline solid.                       #!
!#                                                                            #!
!#                                                                            #!
!#   AUTHORS                                                                  #!
!#   ~~~~~~~                                                                  #!
!#   This set of subprograms is based on code written in an early Fortran     #!
!#   77 version of PWSCF by Alessio Filippetti. These were later ported       #!
!#   into another version by Lixin He. Oswaldo Dieguez, in collaboration      #!
!#   with Lixin He and Jeff Neaton, ported these routines into Fortran 90     #!
!#   version 1.2.1 of PWSCF. He, Dieguez, and Neaton were working at the      #!
!#   time in David Vanderbilt's group at Rutgers, The State University of     #!
!#   New Jersey, USA.                                                         #!
!#                                                                            #!
!#                                                                            #!
!#   LIST OF FILES                                                            #!
!#   ~~~~~~~~~~~~~                                                            #!
!#   The complete list of files added to the PWSCF distribution is:           #!
!#   * ../PW/bp_calc_btq.f90                                                  #!
!#   * ../PW/bp_c_phase.f90                                                   #!
!#   * ../PW/bp_qvan3.f90                                                     #!
!#   * ../PW/bp_strings.f90                                                   #!
!#                                                                            #!
!#   The PWSCF files that needed (minor) modifications were:                  #!
!#   * ../PW/electrons.f90                                                    #!
!#   * ../PW/input.f90                                                        #!
!#   * ../PW/pwcom.f90                                                        #!
!#   * ../PW/setup.f90                                                        #!
!#                                                                            #!
!#   Present in the original version and later removed:                       #!
!#   * bp_ylm_q.f bp_dbess.f bp_radin.f bp_bess.f                             #!
!#                                                                            #!
!#   BRIEF SUMMARY OF THE METHODOLOGY                                         #!
!#   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                         #!
!#   The spontaneous polarization has two contibutions, electronic            #!
!#   and ionic. With these additional routines, PWSCF will output both.       #!
!#                                                                            #!
!#   The ionic contribution is relatively trivial to compute, requiring       #!
!#   knowledge only of the atomic positions and core charges. The new         #!
!#   subroutines focus mainly on evaluating the electronic contribution,      #!
!#   computed as a Berry phase, i.e., a global phase property that can        #!
!#   be computed from inner products of Bloch states at neighboring           #!
!#   points in k-space.                                                       #!
!#                                                                            #!
!#   The standard procedure would be for the user to first perform a          #!
!#   self-consistent (sc) calculation to obtain a converged charge density.   #!
!#   With well-converged sc charge density, the user would then run one       #!
!#   or more non-self consistent (or "band structure") calculations,          #!
!#   using the same main code, but with a flag to ask for the polarization.   #!
!#   Each such run would calculate the projection of the polarization         #!
!#   onto one of the three primitive reciprocal lattice vectors. In           #!
!#   cases of high symmetry (e.g. a tetragonal ferroelectric phase), one      #!
!#   such run would suffice. In the general case of low symmetry, the         #!
!#   user would have to submit up to three jobs to compute the three          #!
!#   components of polarization, and would have to obtain the total           #!
!#   polarization "by hand" by summing these contributions.                   #!
!#                                                                            #!
!#   Accurate calculation of the electronic or "Berry-phase" polarization     #!
!#   requires overlaps between wavefunctions along fairly dense lines (or     #!
!#   "strings") in k-space in the direction of the primitive G-vector for     #!
!#   which one is calculating the projection of the polarization. The         #!
!#   code would use a higher-density k-mesh in this direction, and a          #!
!#   standard-density mesh in the two other directions. See below for         #!
!#   details.                                                                 #!
!#                                                                            #!
!#                                                                            #!
!#   FUNCTIONALITY/COMPATIBILITY                                              #!
!#   ~~~~~~~~~~~~~~~~~~~~~~~~~~~                                              #!
!#   * Berry phases for a given G-vector.                                     #!
!#                                                                            #!
!#   * Contribution to the polarization (in relevant units) for a given       #!
!#     G-vector.                                                              #!
!#                                                                            #!
!#   * Spin-polarized systems supported.                                      #!
!#                                                                            #!
!#   * Ultrasoft and norm-conserving pseudopotentials supported.              #!
!#                                                                            #!
!#   * The value of the "polarization quantum" and the ionic contribution     #!
!#     to the polarization are reported.                                      #!
!#                                                                            #!
!#                                                                            #!
!#   NEW INPUT PARAMETERS                                                     #!
!#   ~~~~~~~~~~~~~~~~~~~~                                                     #!
!#   * lberry (.TRUE. or .FALSE.)                                             #!
!#     Tells PWSCF that a Berry phase calcultion is desired.                  #!
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
!#                                                                            #!
!#   EXPLANATION OF K-POINT MESH                                              #!
!#   ~~~~~~~~~~~~~~~~~~~~~~~~~~~                                              #!
!#   If gdir=1, the program takes the standard input specification of the     #!
!#   k-point mesh (nk1 x nk2 x nk3) and stops if the k-points in dimension    #!
!#   1 are not equally spaced or if its number is not equal to nppstr,        #!
!#   working with a mesh of dimensions (nppstr x nk2 x nk3).  That is, for    #!
!#   each point of the (nk2 x nk3) two-dimensional mesh, it works with a      #!
!#   string of nppstr k-points extending in the third direction.  Symmetry    #!
!#   will be used to reduce the number of strings (and assign them weights)   #!
!#   if possible.  Of course, if gdir=2 or 3, the variables nk2 or nk3 will   #!
!#   be overridden instead, and the strings constructed in those              #!
!#   directions, respectively.                                                #!
!#                                                                            #!
!#                                                                            #!
!#   BIBLIOGRAPHY                                                             #!
!#   ~~~~~~~~~~~~                                                             #!
!#   The theory behind this implementation is described in:                   #!
!#   [1] R D King-Smith and D Vanderbilt, "Theory of polarization of          #!
!#       crystaline solids", Phys Rev B 47, 1651 (1993).                      #!
!#                                                                            #!
!#   Other relevant sources of information are:                               #!
!#   [2] D Vanderbilt and R D King-Smith, "Electronic polarization in the     #!
!#       ultrasoft pseudopotential formalism", internal report (1998),        #!
!#   [3] D Vanderbilt, "Berry phase theory of proper piezoelectric            #!
!#       response", J Phys Chem Solids 61, 147 (2000).                        #!
!#                                                                            #!
!#                                                                            #!
!#                                              dieguez@physics.rutgers.edu   #!
!#                                                             09 June 2003   #!
!#                                                                            #!
!#                                                                            #!
!##############################################################################!


!======================================================================!

SUBROUTINE c_phase

!----------------------------------------------------------------------!

!   Geometric phase calculation along a strip of nppstr k-points
!   averaged over a 2D grid of nkort k-points ortogonal to nppstr 

!  --- Make use of the module with common information ---
   USE kinds,                ONLY : DP
   USE io_global,            ONLY : stdout
   USE io_files,             ONLY : iunwfc, nwordwfc
   USE buffers,              ONLY : get_buffer
   USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau, zv, atm
   USE cell_base,            ONLY : at, alat, tpiba, omega
   USE constants,            ONLY : pi, tpi
   USE gvect,                ONLY : ngm, g, gcutm, ngm_g, ig_l2g
   USE fft_base,             ONLY : dfftp
   USE uspp,                 ONLY : nkb, vkb, okvan
   USE uspp_param,           ONLY : upf, lmaxq, nbetam, nh, nhm
   USE lsda_mod,             ONLY : nspin
   USE klist,                ONLY : nelec, degauss, nks, xk, wk, igk_k, ngk
   USE wvfct,                ONLY : npwx, nbnd, wg
   USE wavefunctions_module, ONLY : evc
   USE bp,                   ONLY : gdir, nppstr, mapgm_global, pdl_tot
   USE becmod,               ONLY : calbec, bec_type, allocate_bec_type, &
                                    deallocate_bec_type
   USE noncollin_module,     ONLY : noncolin, npol, nspin_lsda
   USE spin_orb,             ONLY : lspinorb
   USE mp_bands,             ONLY : intra_bgrp_comm, nproc_bgrp
   USE mp,                   ONLY : mp_sum
   USE qexsd_module,         ONLY : qexsd_init_berryPhaseOutput, qexsd_bp_obj
!  --- Avoid implicit definitions ---
   IMPLICIT NONE

!  --- Internal definitions ---
   INTEGER :: i
   INTEGER :: igk1(npwx)
   INTEGER :: igk0(npwx)
   INTEGER :: ig
   INTEGER :: ind1
   INTEGER :: info
   INTEGER :: is
   INTEGER :: istring
   INTEGER :: iv
   INTEGER :: ivpt(nbnd)
   INTEGER :: j
   INTEGER :: jkb
   INTEGER :: jkb_bp
   INTEGER :: jkb1
   INTEGER :: job
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
   INTEGER , ALLOCATABLE :: mod_elec(:)
   INTEGER , ALLOCATABLE :: ln(:,:,:)
   INTEGER :: mod_elec_dw
   INTEGER :: mod_elec_tot
   INTEGER :: mod_elec_up
   INTEGER :: mod_ion(nat)
   INTEGER :: mod_ion_tot
   INTEGER :: mod_tot
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
   LOGICAL :: lodd
   LOGICAL :: l_para
   LOGICAL, ALLOCATABLE :: l_cal(:) ! flag for occupied/empty states
   REAL(DP) :: dk(3)
   REAL(DP) :: dkmod
   REAL(DP) :: el_loc
   REAL(DP) :: eps
   REAL(DP) :: fac
   REAL(DP) :: gpar(3)
   REAL(DP) :: gtr(3)
   REAL(DP) :: gvec
   REAL(DP), ALLOCATABLE :: loc_k(:)
   REAL(DP), ALLOCATABLE :: pdl_elec(:)
   REAL(DP), ALLOCATABLE :: phik(:)
   REAL(DP) :: phik_ave
   REAL(DP) :: qrad_dk(nbetam,nbetam,lmaxq,ntyp)
   REAL(DP) :: weight
   REAL(DP) :: upol(3)
   REAL(DP) :: pdl_elec_dw
   REAL(DP) :: pdl_elec_tot
   REAL(DP) :: pdl_elec_up
   REAL(DP) :: pdl_ion(nat)
   REAL(DP) :: pdl_ion_tot
   REAL(DP) :: phidw
   REAL(DP) :: phiup
   REAL(DP) :: rmod
   REAL(DP), ALLOCATABLE :: wstring(:)
   REAL(DP) :: ylm_dk(lmaxq*lmaxq)
   REAL(DP) :: zeta_mod
   COMPLEX(DP), ALLOCATABLE :: aux(:)
   COMPLEX(DP), ALLOCATABLE :: aux_g(:)
   COMPLEX(DP), ALLOCATABLE :: aux0(:)
   TYPE (bec_type) :: becp0
   TYPE (bec_type) :: becp_bp
   COMPLEX(DP) :: cave
   COMPLEX(DP) , ALLOCATABLE :: cphik(:)
   COMPLEX(DP) :: det
   COMPLEX(DP) :: dtheta
   COMPLEX(DP) :: mat(nbnd,nbnd)
   COMPLEX(DP) :: pref
   COMPLEX(DP), ALLOCATABLE :: psi(:,:)
   COMPLEX(DP), ALLOCATABLE :: q_dk_so(:,:,:,:)
   COMPLEX(DP) :: q_dk(nhm,nhm,ntyp)
   COMPLEX(DP) :: struc(nat)
   COMPLEX(DP) :: theta0
   COMPLEX(DP) :: zdotc
   COMPLEX(DP) :: zeta

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
   WRITE( stdout,"(28X,'POLARIZATION CALCULATION')")
   WRITE( stdout,"(25X,'!!! NOT THOROUGHLY TESTED !!!')")
   WRITE( stdout,"(15X,50('-'),/)")

!  --- Check that we are working with an insulator with no empty bands ---
   IF ( degauss > 0.0_dp ) CALL errore('c_phase', &
                'Polarization only for insulators',1)

!  --- Define a small number ---
   eps=1.0E-6_dp

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
   nkort=nstring/nspin_lsda

!  --- Allocate memory for arrays ---
   ALLOCATE(phik(nstring))
   ALLOCATE(loc_k(nstring))
   ALLOCATE(cphik(nstring))
   ALLOCATE(wstring(nstring))
   ALLOCATE(pdl_elec(nstring))
   ALLOCATE(mod_elec(nstring))

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
     CALL errore('c_phase','Wrong k-strings?',1)

!  --- Check that k-points form strings ---
   DO i=1,nspin_lsda*nkort
      DO j=2,nppstr
         kindex=j+(i-1)*nppstr
         IF (ABS(xk(1,kindex)-xk(1,kindex-1)-dk(1)) > eps) &
            CALL errore('c_phase','Wrong k-strings?',1)
         IF (ABS(xk(2,kindex)-xk(2,kindex-1)-dk(2)) > eps) &
            CALL errore('c_phase','Wrong k-strings?',1)
         IF (ABS(xk(3,kindex)-xk(3,kindex-1)-dk(3)) > eps) &
            CALL errore('c_phase','Wrong k-strings?',1)
         IF (ABS(wk(kindex)-wk(kindex-1)) > eps) &
            CALL errore('c_phase','Wrong k-strings weights?',1)
      END DO
   END DO

!  -------------------------------------------------------------------------   !
!                   electronic polarization: weight strings                    !
!  -------------------------------------------------------------------------   !

!  --- Calculate string weights, normalizing to 1 (no spin or noncollinear)
!       or 1+1 (spin) ---
   DO is=1,nspin_lsda
      weight=0.0_dp
      DO kort=1,nkort
         istring=kort+(is-1)*nkort
         wstring(istring)=wk(nppstr*istring)
         weight=weight+wstring(istring)
      END DO
      DO kort=1,nkort
         istring=kort+(is-1)*nkort
         wstring(istring)=wstring(istring)/weight
      END DO
   END DO

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

   el_loc=0.d0
   kpoint=0
   ALLOCATE ( l_cal(nbnd) ) 
   CALL weights()

!  --- Start loop over spin ---
   DO is=1,nspin_lsda

      ! l_cal(n) = .true./.false. if n-th state is occupied/empty
      nbnd_occ=0
      DO nb = 1, nbnd
         l_cal(nb) = (wg(nb,1+nks*(is-1)/2) > eps)
         IF (l_cal(nb)) nbnd_occ = nbnd_occ + 1
      END DO

!     --- Start loop over orthogonal k-points ---
      DO kort=1,nkort

!        --- Index for this string ---
         istring=kort+(is-1)*nkort

!        --- Initialize expectation value of the phase operator ---
         zeta=(1.d0,0.d0)
         zeta_mod = 1.d0

!        --- Start loop over parallel k-points ---
         DO kpar = 1,nppstr

!           --- Set index of k-point ---
            kpoint = kpoint + 1

!           --- Calculate dot products between wavefunctions and betas ---
            IF (kpar /= 1) THEN

!              --- Dot wavefunctions and betas for PREVIOUS k-point ---
               npw0 = ngk(kpoint-1)
               igk0(:) = igk_k(:,kpoint-1)
               CALL get_buffer (psi,nwordwfc,iunwfc,kpoint-1)
               if (okvan) then
                  CALL init_us_2 (npw0,igk0,xk(1,kpoint-1),vkb)
                  CALL calbec (npw0, vkb, psi, becp0)
               endif
!              --- Dot wavefunctions and betas for CURRENT k-point ---
               IF (kpar /= nppstr) THEN
                  npw1 = ngk(kpoint)
                  igk1(:) = igk_k(:,kpoint)
                  CALL get_buffer(evc,nwordwfc,iunwfc,kpoint)
                  if (okvan) then
                     CALL init_us_2 (npw1,igk1,xk(1,kpoint),vkb)
                     CALL calbec (npw1, vkb, evc, becp_bp)
                  endif
               ELSE
                  kstart = kpoint-nppstr+1
                  npw1 = ngk(kstart)
                  igk1(:) = igk_k(:,kstart)
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
                           CALL errore('c_phase','wrong g',1)
                        ENDIF
                     ELSE
                        WRITE(6,*) ' |gtr| > gcutm  for gtr=', &
                             gtr(1),gtr(2),gtr(3)
                        CALL errore('c_phase','wrong gtr',1)
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

!              --- Calculate matrix determinant ---
               CALL ZGETRF (nbnd,nbnd,mat,nbnd,ivpt,info)
               CALL errore('c_phase','error in factorization',abs(info))
               det=(1.d0,0.d0)
               do nb=1,nbnd
                  det = det*mat(nb,nb)
                  if(nb.ne.ivpt(nb)) det=-det
               enddo
!              --- Multiply by the already calculated determinants ---
               zeta=zeta*det

!           --- End of dot products between wavefunctions and betas ---
            ENDIF

!        --- End loop over parallel k-points ---
         END DO 

!        --- Calculate the phase for this string ---
         phik(istring)=AIMAG(LOG(zeta))
         cphik(istring)=COS(phik(istring))*(1.0_dp,0.0_dp) &
                     +SIN(phik(istring))*(0.0_dp,1.0_dp)

!        --- Calculate the localization for current kort ---
         zeta_mod= DBLE(CONJG(zeta)*zeta)

         loc_k(istring)= - (nppstr-1) / gvec**2 / nbnd_occ *log(zeta_mod)

!     --- End loop over orthogonal k-points ---
      END DO

!  --- End loop over spin ---
   END DO
   DEALLOCATE ( l_cal ) 

!  -------------------------------------------------------------------------   !
!                    electronic polarization: phase average                    !
!  -------------------------------------------------------------------------   !

!  --- Start loop over spins ---
   DO is=1,nspin_lsda

!  --- Initialize average of phases as complex numbers ---
      cave=(0.0_dp,0.0_dp)
      phik_ave=(0.0_dp,0.0_dp)

!     --- Start loop over strings with same spin ---
      DO kort=1,nkort

!        --- Calculate string index ---
         istring=kort+(is-1)*nkort

!        --- Average phases as complex numbers ---
         cave=cave+wstring(istring)*cphik(istring)

!     --- End loop over strings with same spin ---
      END DO

!     --- Get the angle corresponding to the complex numbers average ---
      theta0=atan2(AIMAG(cave), DBLE(cave))
!     --- Put the phases in an around theta0 ---
      DO kort=1,nkort
        istring=kort+(is-1)*nkort
        cphik(istring)=cphik(istring)/cave
        dtheta=atan2(AIMAG(cphik(istring)), DBLE(cphik(istring)))
        phik(istring)=theta0+dtheta
        phik_ave=phik_ave+wstring(istring)*phik(istring)
      END DO

!     --- Assign this angle to the corresponding spin phase average ---
      IF (nspin == 1) THEN
         phiup=phik_ave !theta0+dtheta
         phidw=phik_ave !theta0+dtheta
      ELSE IF (nspin == 2) THEN
         IF (is == 1) THEN
            phiup=phik_ave !theta0+dtheta
         ELSE IF (is == 2) THEN
            phidw=phik_ave !theta0+dtheta
         END IF
      ELSE IF (nspin==4 ) THEN
         phiup=phik_ave
         phidw=0.0_DP
      END IF

!  --- End loop over spins
   END DO

!  -------------------------------------------------------------------------   !
!                     electronic polarization: remap phases                    !
!  -------------------------------------------------------------------------   !

!  --- Remap string phases to interval [-0.5,0.5) ---
   pdl_elec=phik/(2.0_dp*pi)
   mod_elec=1

!  --- Remap spin average phases to interval [-0.5,0.5) ---
   pdl_elec_up=phiup/(2.0_dp*pi)
   mod_elec_up=1
   pdl_elec_dw=phidw/(2.0_dp*pi)
   mod_elec_dw=1

!  --- Depending on nspin, remap total phase to [-1,1) or [-0.5,0.5) ---
   pdl_elec_tot=pdl_elec_up+pdl_elec_dw
   IF (nspin == 1) THEN
      pdl_elec_tot=pdl_elec_tot-2.0_dp*NINT(pdl_elec_tot/2.0_dp)
      mod_elec_tot=2
   ELSE IF (nspin == 2 .OR. nspin == 4) THEN
      pdl_elec_tot=pdl_elec_tot-1.0_dp*NINT(pdl_elec_tot/1.0_dp)
      mod_elec_tot=1
   END IF

!  -------------------------------------------------------------------------   !
!                              ionic polarization                              !
!  -------------------------------------------------------------------------   !

!  --- Look for ions with odd number of charges ---
   mod_ion=2
   lodd=.FALSE.
   DO na=1,nat
      IF (MOD(NINT(zv(ityp(na))),2) == 1) THEN
         mod_ion(na)=1
         lodd=.TRUE.
      END IF
   END DO

!  --- Calculate ionic polarization phase for every ion ---
   pdl_ion=0.0_dp
   DO na=1,nat
      DO i=1,3
         pdl_ion(na)=pdl_ion(na)+zv(ityp(na))*tau(i,na)*gpar(i)
      ENDDO
      IF (mod_ion(na) == 1) THEN
         pdl_ion(na)=pdl_ion(na)-1.0_dp*nint(pdl_ion(na)/1.0_dp)
      ELSE IF (mod_ion(na) == 2) THEN
         pdl_ion(na)=pdl_ion(na)-2.0_dp*nint(pdl_ion(na)/2.0_dp)
      END IF
   ENDDO

!  --- Add up the phases modulo 2 iff the ionic charges are even numbers ---
   pdl_ion_tot=SUM(pdl_ion(1:nat))
   IF (lodd) THEN
      pdl_ion_tot=pdl_ion_tot-1.d0*nint(pdl_ion_tot/1.d0)
      mod_ion_tot=1
   ELSE
      pdl_ion_tot=pdl_ion_tot-2.d0*nint(pdl_ion_tot/2.d0)
      mod_ion_tot=2
   END IF

!  -------------------------------------------------------------------------   !
!                              total polarization                              !
!  -------------------------------------------------------------------------   !

!  --- Add electronic and ionic contributions to total phase ---
   pdl_tot=pdl_elec_tot+pdl_ion_tot
   IF ((.NOT.lodd).AND.(nspin == 1)) THEN
      mod_tot=2
   ELSE
      mod_tot=1
   END IF

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

!  --- Information about ionic polarization phases ---
   WRITE( stdout,"(2/,31X,'IONIC POLARIZATION')")
   WRITE( stdout,"(31X,18('~'),/)")
   WRITE( stdout,"(8X,'Note: (mod 1) means that the phases (angles ranging from' &
           & /,8X,'-pi to pi) have been mapped to the interval [-1/2,+1/2) by',&
           & /,8X,'dividing by 2*pi; (mod 2) refers to the interval [-1,+1)',&
           & /)")
   WRITE( stdout,"(2X,76('='))")
   WRITE( stdout,"(4X,'Ion',4X,'Species',4X,'Charge',14X, &
           & 'Position',16X,'Phase')")
   WRITE( stdout,"(2X,76('-'))")
   DO na=1,nat
      WRITE( stdout,"(3X,I3,8X,A2,F12.3,5X,3F8.4,F12.5,' (mod ',I1,')')") &
           & na,atm(ityp(na)),zv(ityp(na)), &
           & tau(1,na),tau(2,na),tau(3,na),pdl_ion(na),mod_ion(na)
   END DO
   WRITE( stdout,"(2X,76('-'))")
   WRITE( stdout,"(47X,'IONIC PHASE: ',F9.5,' (mod ',I1,')')") pdl_ion_tot,mod_ion_tot
   WRITE( stdout,"(2X,76('='))")

!  --- Information about electronic polarization phases ---
   WRITE( stdout,"(2/,28X,'ELECTRONIC POLARIZATION')")
   WRITE( stdout,"(28X,23('~'),/)")
   WRITE( stdout,"(8X,'Note: (mod 1) means that the phases (angles ranging from' &
           & /,8X,'-pi to pi) have been mapped to the interval [-1/2,+1/2) by',&
           & /,8X,'dividing by 2*pi; (mod 2) refers to the interval [-1,+1)',&
           & /)")
   WRITE( stdout,"(2X,76('='))")
   WRITE( stdout,"(3X,'Spin',4X,'String',5X,'Weight',6X, &
            &  'First k-point in string',9X,'Phase')")
   WRITE( stdout,"(2X,76('-'))")
   DO istring=1,nstring/nspin_lsda
      ind1=1+(istring-1)*nppstr
      WRITE( stdout,"(3X,' up ',3X,I5,F14.6,4X,3(F8.4),F12.5,' (mod ',I1,')')") &
          &  istring,wstring(istring), &
          &  xk(1,ind1),xk(2,ind1),xk(3,ind1),pdl_elec(istring),mod_elec(istring)
   END DO
   WRITE( stdout,"(2X,76('-'))")
!  --- Treat unpolarized/polarized spin cases ---
   IF (nspin_lsda == 1) THEN
!     --- In unpolarized spin, just copy again the same data ---
      DO istring=1,nstring
         ind1=1+(istring-1)*nppstr
         WRITE( stdout,"(3X,'down',3X,I5,F14.6,4X,3(F8.4),F12.5,' (mod ',I1,')')") &
              istring,wstring(istring), xk(1,ind1),xk(2,ind1),xk(3,ind1), &
              pdl_elec(istring),mod_elec(istring)
      END DO
   ELSE IF (nspin_lsda == 2) THEN
!     --- If there is spin polarization, write information for new strings ---
      DO istring=nstring/2+1,nstring
         ind1=1+(istring-1)*nppstr
         WRITE( stdout,"(3X,'down',3X,I4,F15.6,4X,3(F8.4),F12.5,' (mod ',I1,')')") &
           &    istring,wstring(istring), xk(1,ind1),xk(2,ind1),xk(3,ind1), &
           &    pdl_elec(istring),mod_elec(istring)
      END DO
   END IF
   WRITE( stdout,"(2X,76('-'))")
   IF (noncolin) THEN
      WRITE( stdout,"(42X,'Average phase   : ',F9.5,' (mod ',I1,')')") & 
        pdl_elec_up,mod_elec_up
   ELSE
      WRITE( stdout,"(40X,'Average phase (up): ',F9.5,' (mod ',I1,')')") & 
        pdl_elec_up,mod_elec_up
      WRITE( stdout,"(38X,'Average phase (down): ',F9.5,' (mod ',I1,')')")& 
        pdl_elec_dw,mod_elec_dw
      WRITE( stdout,"(42X,'ELECTRONIC PHASE: ',F9.5,' (mod ',I1,')')") & 
        pdl_elec_tot,mod_elec_tot
   ENDIF
   WRITE( stdout,"(2X,76('='))")

!  --- Information about total phase ---
   WRITE( stdout,"(2/,31X,'SUMMARY OF PHASES')")
   WRITE( stdout,"(31X,17('~'),/)")
   WRITE( stdout,"(26X,'Ionic Phase:',F9.5,' (mod ',I1,')')") &
        pdl_ion_tot,mod_ion_tot
   WRITE( stdout,"(21X,'Electronic Phase:',F9.5,' (mod ',I1,')')") &
        pdl_elec_tot,mod_elec_tot
   WRITE( stdout,"(26X,'TOTAL PHASE:',F9.5,' (mod ',I1,')')") &
        pdl_tot,mod_tot

!  --- Information about the value of polarization ---
   WRITE( stdout,"(2/,29X,'VALUES OF POLARIZATION')")
   WRITE( stdout,"(29X,22('~'),/)")
   WRITE( stdout,"( &
      &   8X,'The calculation of phases done along the direction of vector ',I1, &
      &   /,8X,'of the reciprocal lattice gives the following contribution to', &
      &   /,8X,'the polarization vector (in different units, and being Omega', &
      &   /,8X,'the volume of the unit cell):')") &
          gdir
!  --- Calculate direction of polarization and modulus of lattice vector ---
   rmod=SQRT(at(1,gdir)*at(1,gdir)+at(2,gdir)*at(2,gdir) &
            +at(3,gdir)*at(3,gdir))
   upol(:)=at(:,gdir)/rmod
   rmod=alat*rmod
!  --- Give polarization in units of (e/Omega).bohr ---
   fac=rmod
   WRITE( stdout,"(/,11X,'P = ',F11.7,'  (mod ',F11.7,')  (e/Omega).bohr')") &
        fac*pdl_tot,fac*DBLE(mod_tot)
!  --- Give polarization in units of e.bohr ---
   fac=rmod/omega
   WRITE( stdout,"(/,11X,'P = ',F11.7,'  (mod ',F11.7,')  e/bohr^2')") &
        fac*pdl_tot,fac*DBLE(mod_tot)
!  --- Give polarization in SI units (C/m^2) ---
   fac=(rmod/omega)*(1.60097E-19_dp/5.29177E-11_dp**2)
   WRITE( stdout,"(/,11X,'P = ',F11.7,'  (mod ',F11.7,')  C/m^2')") &
        fac*pdl_tot,fac*DBLE(mod_tot)
!  --- Write polarization direction ---
   WRITE( stdout,"(/,8X,'The polarization direction is:  (', &
       &  F8.5,' ,',F8.5,' ,',F8.5,' )')") upol(1),upol(2),upol(3)

!  --- End of information relative to polarization calculation ---
   WRITE( stdout,"(/,/,15X,50('=')/,/)")
!------------------------------------------------------------------------------
!                            INITIALIZE  QEXSD OUTPUT ELEMENT
! Here we write all output information in a berry_phase_type variable to print
! them in the XML output  P.D. april 2016
!------------------------------------------------------------------------------
  CALL qexsd_init_berryPhaseOutput(qexsd_bp_obj, gpar, gvec, nppstr, nkort, xk, pdl_ion, mod_ion,  &
                                  pdl_ion_tot, mod_ion_tot, nstring, pdl_elec , mod_elec, wstring, &
                                  pdl_elec_up, mod_elec_up, pdl_elec_dw, mod_elec_dw, pdl_elec_tot,&
                                  mod_elec_tot, pdl_tot, mod_tot, upol, rmod)
!  -------------------------------------------------------------------------   !
!                                  finalization                                !
!  -------------------------------------------------------------------------   !

!  --- Free memory ---
   DEALLOCATE(mod_elec)
   DEALLOCATE(pdl_elec)
   DEALLOCATE(wstring)
   DEALLOCATE(cphik)
   DEALLOCATE(loc_k)
   DEALLOCATE(phik)
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

!------------------------------------------------------------------------------!

END SUBROUTINE c_phase

!==============================================================================!
