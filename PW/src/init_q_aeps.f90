!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE init_q_aeps ( )
   !------------------------------------------------------------------------
   !
   ! Initialization of the pseudopotential-dependent quantities needed for
   ! LDA+U method with projections computed through <beta|psi>:
   !  q_ae = integral of the AE wfc up to r_core
   !  q_ps = integral of the PS wfc up to r_core (not used at the moment)
   !
   USE kinds,      ONLY : DP
   USE ions_base,  ONLY : ntyp => nsp, ityp, nat
   USE atom,       ONLY : rgrid, msh
   USE lsda_mod,   ONLY : nspin
   USE ldaU,       ONLY : q_ae, q_ps, Hubbard_l, &
                          U_projection, is_hubbard, nwfcU, offsetU
   USE uspp_param, ONLY : nbetam, nh, nhm, upf
   USE uspp,       ONLY : indv, nhtol, nhtolm, nkb
   USE control_flags, ONLY : iverbosity
   USE io_global, ONLY : ionode
   !
   IMPLICIT NONE
   ! LOCAL
   INTEGER :: l, m, mb, nb, ndm, cnt, kk, iwfc, jwfc
   INTEGER :: nt, nt_, na, ih, jh, ib, jb, lH, nchiH, nbH
   REAL(DP), ALLOCATABLE :: aux (:), qq_ae(:,:,:), qq_ps(:,:,:)
   REAL(DP) :: psint, aeint, wsgn
   !
   !
   ALLOCATE ( q_ae(nwfcU,nhm,nat), q_ps(nwfcU,nhm,nat) )
   !
   ndm = MAXVAL (msh(1:ntyp))
   ALLOCATE ( aux(ndm), qq_ae(nbetam,nbetam,ntyp), qq_ps(nbetam,nbetam,ntyp) )
   !
   qq_ae(:,:,:) = 0.0_DP
   qq_ps(:,:,:) = 0.0_DP
   q_ae(:,:,:)  = 0.0_DP
   q_ps(:,:,:)  = 0.0_DP
   !
   !
   ! Compute the integrals of the AE and PS wavefunctions up to core radii
   ! (only for atomic types entering in the Hubbard Hamiltonian)
   !
   DO nt = 1, ntyp
      !
      IF ( .NOT. is_hubbard(nt) ) CYCLE
      !
      IF ( .NOT.upf(nt)%has_wfc ) CALL errore('init_q_aeps', &
         "All-electron atomic-wavefunctions needed for pseudo U_projection",1)
      !
      DO nb = 1, upf(nt)%nbeta
         !
         DO mb = nb, upf(nt)%nbeta
            !
            IF ( upf(nt)%lll(mb) == upf(nt)%lll(nb) ) then
               !
               kk = MAX(upf(nt)%kbeta(mb),upf(nt)%kbeta(nb)) ! needed ???
               aux(1:msh(nt)) = upf(nt)%aewfc(1:msh(nt),mb)*upf(nt)%aewfc(1:msh(nt),nb)
               CALL simpson (upf(nt)%kbeta(nb),aux,rgrid(nt)%rab,aeint)
               qq_ae(nb,mb,nt) = aeint
               qq_ae(mb,nb,nt) = aeint
               aux(1:msh(nt)) = upf(nt)%pswfc(1:msh(nt),mb)*upf(nt)%pswfc(1:msh(nt),nb)
               CALL simpson (upf(nt)%kbeta(nb),aux,rgrid(nt)%rab,psint)
               qq_ps(nb,mb,nt) = psint
               qq_ps(mb,nb,nt) = psint
               !
            ENDIF
         ENDDO
      ENDDO

      !!! WARNING: when generated with lsave_wfc, the PP file contains the 
      !!! AE and PS wfcs for every beta projector (in principle more than one
      !!! for each l). We identify which beta corresponds to the bound state 
      !!! by checking the norm of the difference |pswfc(r) - chi(r)|
      lH = Hubbard_l(nt)
      nbH = -1
      !
      ! select chi corresponding to bound states (the same used to build initial
      ! wfcs) AND with l = Hubbard_l (only for species with Hubbard_l defined)
      !
      IF ( lH .GE. 0 ) THEN
         !
         !!! NOTE: one might run into troubles when using a PP with semicore 
         !!! states with same l as valence states (also otherwhere for LDA+U)
         DO nb = 1, upf(nt)%nwfc
            IF (upf(nt)%lchi(nb) == lH .AND. upf(nt)%oc(nb) >= 0.d0)   nchiH = nb
         ENDDO
         !
         DO nb = 1, upf(nt)%nbeta
            !
            IF (upf(nt)%lll(nb) == lH) THEN
               ! check if chi and pswfc have the same sign or not
               aux(1:msh(nt)) = upf(nt)%pswfc(:,nb)*upf(nt)%chi(:,nchiH)
               CALL simpson(msh(nt),aux,rgrid(nt)%rab,psint)
               wsgn = sign(1.0_DP,psint)
               ! compute norm of the difference [pswfc(r) - chi(r)]
               aux(1:msh(nt)) = (upf(nt)%pswfc(:,nb) - wsgn*upf(nt)%chi(:,nchiH))**2
               CALL simpson(msh(nt),aux,rgrid(nt)%rab,psint)
               IF ( abs(psint) .LE. 1.d-9 ) nbH = nb
            ENDIF
            !
         ENDDO
         !
         !!! DEBUG
         if ( ionode .AND. iverbosity == 1 ) then
            write(*,*) '> QQ_AE matrix:'
            do nb = 1, upf(nt)%nbeta
               write(*,'(99F9.6)') qq_ae(nb,1:upf(nt)%nbeta,nt)
            enddo
            write(*,*) "nbH=",nbH,", lH",lH
         endif
         !!!
         !
         IF ( nbH .EQ. -1 ) CALL errore("init_q_aeps", "could not set nbH", 1)
         !
      ENDIF

      cnt = 0
      !
      ! initialize q_ae and q_ps for U projectors on beta functions (in the solid)
      DO na = 1, nat ! on atoms
         !
         nt_ = ityp(na)
         ! offset for atomic wavefunctions (initialized in offset_atom_wfc)
         iwfc = offsetU(na)

         IF ( nt_ == nt .AND. lH .GE. 0 ) THEN
            !
            !!! DEBUG
            if ( ionode .AND. iverbosity == 1 ) then
               write(*,*) "na, ityp, lH=",na,ityp(na),lH
               write(*,*) "nbH,lH,offset",nbH,lH,iwfc
            endif
            !!!

            DO jh = 1, nh(nt)
               !
               IF (nhtol(jh,nt) .NE. lH) CYCLE
               jb = indv(jh,nt)
               !
               DO ih = 1, nh(nt)
                  !
                  ib = indv(ih,nt)
                  IF (nhtol(ih,nt) .NE. lH) CYCLE
                  IF (ib .NE. nbH) CYCLE
                  IF ( nhtolm(ih,nt)==nhtolm(jh,nt) ) THEN
                     !
                     m=nhtolm(ih,nt)-lH*lH
                     !!! DEBUG
                     if ( ionode .AND. iverbosity == 1 ) write(*,'(A,6I3,F9.6)') &
                        "jh,ih,nhtolm,lH,m,iwfc+m,qq",jh,ih,nhtolm(ih,nt),lH,m,iwfc+m,qq_ae(jb,ib,nt)
                     !!!
                     !
                     q_ae(iwfc+m,jh,na) = qq_ae(jb,ib,nt)
                     q_ps(iwfc+m,jh,na) = qq_ps(jb,ib,nt)
                     !
                     !!! DEBUG
                     if ( ionode .AND. iverbosity == 1 ) THEN
                        write(*,'(A,3I3,2F9.6)') "iwfc,jh,na,q_ae,qq_ae", &
                           iwfc+m,jh,na,q_ae(iwfc+m,jh,na),qq_ae(jb,ib,nt)
                        write(*,'(A,3I3,2F9.6)') "iwfc,jh,na,q_ps,qq_ps", &
                           iwfc+m,jh,na,q_ps(iwfc+m,jh,na),qq_ps(jb,ib,nt)
                     endif
                     !!!
                     !
                  ENDIF
                  !
               ENDDO
               ! ih
            ENDDO
            ! jh
         ENDIF
         ! ityp
      ENDDO
      ! on atoms
      !
   ENDDO 
   ! on atomic types

   !!! DEBUG
   if ( ionode .AND. iverbosity == 1 ) then
      iwfc=0
      do na = 1,nat
         nt = ityp(na)
         write(*,*) ">>> atom ",na,", type ",nt
         jwfc=iwfc
         write(*,*) "    q_ae matrix"
         do nb = 1, upf(nt)%nwfc
            if (upf(nt)%oc(nb) >= 0.d0) then
               l = upf(nt)%lchi(nb)
               do m = 1,2*l+1
                  jwfc=jwfc+1
                  write(*,'(2I1,99F6.3)') l,m,q_ae(jwfc,:,na)
                  !
               enddo
            endif
         enddo
         !
         jwfc=iwfc
         write(*,*) "    q_ps matrix"
         do nb = 1, upf(nt)%nwfc
            if (upf(nt)%oc(nb) >= 0.d0) then
               l = upf(nt)%lchi(nb)
               do m = 1,2*l+1
                  jwfc=jwfc+1
                  write(*,'(2I1,99F6.3)') l,m,q_ps(jwfc,:,na)
                  !
               enddo
            endif
         enddo
         !
         iwfc=jwfc
         !
      enddo
   endif
   !!!
   !
   !
   deallocate( aux, qq_ae, qq_ps )
   !
   RETURN

END SUBROUTINE init_q_aeps
!
