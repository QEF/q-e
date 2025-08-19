MODULE oscdft_context
#if defined (__OSCDFT)
   USE kinds,     ONLY : DP
   USE io_global, ONLY : ionode, ionode_id, stdout
   USE mp,        ONLY : mp_bcast
   USE mp_images, ONLY : intra_image_comm

   USE oscdft_input,        ONLY : oscdft_input_type
   USE oscdft_indices,      ONLY : oscdft_indices_type,&
                                   oscdft_orbital_indices_type,&
                                   oscdft_constr_indices_type
   USE oscdft_wavefunction, ONLY : oscdft_wavefunction_type
   USE oscdft_forces,       ONLY : oscdft_forces_type
   USE oscdft_enums

   PRIVATE
   PUBLIC oscdft_context_type,&
          oscdft_init_context,&
          oscdft_ns_type,&
          oscdft_alloc_nst,&
          oscdft_init_indices

   TYPE oscdft_ns_type
      LOGICAL               :: eigvects_set
      REAL(DP), ALLOCATABLE :: ns(:,:,:),& ! ioscdft
                               numbers(:),& ! ioscdft
                               eigvals(:,:),& ! ioscdft
                               eigvects(:,:,:),& ! ioscdft
                               occup_numbers(:),& ! iconstr
                               occup_eigvals(:,:),& ! iconstr
                               occup_eigvects(:,:,:),& ! iconstr
                               gradient(:) ! iconstr
   END TYPE oscdft_ns_type

   TYPE oscdft_context_type
      LOGICAL                        :: initialized    = .false.,&
                                        warming_up     = .true.,&
                                        recalculate_ns = .false.,&
                                        wfc_allocated = .false.,&
                                        has_debug,&
                                        preserve_iter_when_restart_done = .false.,&
                                        conv_thr_not_minimum
      INTEGER                        :: warm_up_iter,&
                                        iter,&
                                        multiplier_iter,&
                                        global_start_index
      REAL(DP)                       :: conv_thr,&
                                        energy,&
                                        old_total_energy,&
                                        total_energy
      REAL(DP), ALLOCATABLE          :: multipliers(:),& ! iconstr
                                        old_multipliers(:)
      TYPE(oscdft_input_type)        :: inp
      TYPE(oscdft_indices_type)      :: idx
      TYPE(oscdft_ns_type)           :: nst
      TYPE(oscdft_forces_type)       :: forces
      TYPE(oscdft_wavefunction_type) :: wfcO, wfcS

      REAL(DP), ALLOCATABLE          :: force_oscdft(:,:)

      ! oscdft_type=2
      LOGICAL, ALLOCATABLE           :: constraining(:)
      LOGICAL                        :: conv, &
                                        is_constraint
      REAL(DP), ALLOCATABLE          :: constraint(:,:,:,:) 

   END TYPE oscdft_context_type

   CONTAINS
      SUBROUTINE parse_orbital_desc(idx, inp)
         IMPLICIT NONE

         TYPE(oscdft_indices_type), TARGET, INTENT(INOUT) :: idx
         TYPE(oscdft_input_type),           INTENT(INOUT) :: inp
         TYPE(oscdft_orbital_indices_type), POINTER       :: orbs

         INTEGER :: ioscdft, ierr, testvar, istr, iorb

         orbs => idx%orbs

         orbs%norbs = 0
         DO ioscdft=1,inp%noscdft
            inp%orbital_desc(ioscdft) = TRIM(inp%orbital_desc(ioscdft))
            orbs%norbs = orbs%norbs + count_norbitals(inp%orbital_desc(ioscdft), ioscdft)
         END DO
         ALLOCATE(orbs%iat(orbs%norbs),&
                  orbs%n(orbs%norbs),&
                  orbs%l(orbs%norbs),&
                  orbs%iorb_start(inp%noscdft),&
                  orbs%iorb_end(inp%noscdft))
         iorb = 1
         DO ioscdft=1,inp%noscdft
            orbs%iorb_start(ioscdft) = iorb
            CALL parse_norbitals(orbs, inp%orbital_desc(ioscdft), ioscdft, iorb)
            orbs%iorb_end(ioscdft) = iorb - 1
            IF (iorb == orbs%iorb_start(ioscdft)) THEN
               orbs%iorb_start(ioscdft) = 0
               orbs%iorb_end(ioscdft) = 0
            END IF
         END DO

         CALL mp_bcast(orbs%norbs,      ionode_id, intra_image_comm)
         CALL mp_bcast(orbs%n,          ionode_id, intra_image_comm)
         CALL mp_bcast(orbs%l,          ionode_id, intra_image_comm)
         CALL mp_bcast(orbs%iorb_start, ionode_id, intra_image_comm)
         CALL mp_bcast(orbs%iorb_end,   ionode_id, intra_image_comm)

         CONTAINS
            SUBROUTINE parse_norbitals(orbs, desc, ioscdft, iorb)
               USE ions_base, ONLY : nat

               IMPLICIT NONE

               TYPE(oscdft_orbital_indices_type), INTENT(INOUT) :: orbs
               CHARACTER(LEN=*), INTENT(IN)    :: desc
               INTEGER,          INTENT(IN)    :: ioscdft
               INTEGER,          INTENT(INOUT) :: iorb
               INTEGER                         :: at ! 0: parse atom, 1: parse principal, 2: parse azimuthal
               INTEGER                         :: istr, jstr, testvar, ierr
               INTEGER                         :: iat, n, l
               LOGICAL                         :: is_paren

               istr = 1
               iat  = -1
               n    = -1
               l    = -1
               at   = 0
               DO WHILE (istr <= LEN_TRIM(desc))
                  SELECT CASE (at)
                     CASE (0)
                        jstr = istr
                        DO WHILE (istr <= LEN_TRIM(desc))
                           is_paren = (desc(istr:istr) == '(')
                           istr = istr + 1
                           IF (is_paren) EXIT
                        END DO
                        READ(desc(jstr:istr-2), *, IOSTAT=ierr) testvar
                        IF (ierr /= 0) THEN
                           CALL errore("oscdft_parse_orbital_desc",&
                                       "expect integer for atom number: " // desc,&
                                       ioscdft)
                        END IF
                        at = 1
                        iat = testvar
                        IF (iat <= 0 .OR. iat > nat) THEN
                           WRITE(stdout, 100) iat
                           CALL errore("oscdft_parse_orbital_desc",&
                                       "malformed orbital descriptor: " // desc,&
                                       ioscdft)
                        END IF
                     CASE (1)
                        READ(desc(istr:istr), *, IOSTAT=ierr) testvar
                        IF (ierr /= 0) THEN
                           CALL errore("oscdft_parse_orbital_desc",&
                                       "expect integer for principal quantum number: " // desc,&
                                       ioscdft)
                        END IF
                        istr = istr + 1
                        at = 2
                        n = testvar
                     CASE (2)
                        SELECT CASE (desc(istr:istr))
                           CASE ('S')
                              l = 0
                           CASE ('P')
                              l = 1
                           CASE ('D')
                              l = 2
                           CASE DEFAULT
                              CALL errore("oscdft_parse_orbital_desc",&
                                          "expect S,P,D: " // desc,&
                                          ioscdft)
                        END SELECT
                        IF (iat == -1 .OR. n == -1 .OR. l == -1) THEN
                              CALL errore("oscdft_parse_orbital_desc",&
                                          "malformed orbital descriptor: " // desc,&
                                          ioscdft)
                        END IF
                        orbs%iat(iorb) = iat
                        orbs%n(iorb)   = n
                        orbs%l(iorb)   = l

                        n = -1
                        l = -1
                        istr = istr + 1
                        iorb = iorb + 1
                        IF (istr > LEN_TRIM(desc)) THEN
                           CALL errore("oscdft_parse_orbital_desc",&
                                       "unexpected end of descriptor: " // desc,&
                                       ioscdft)
                        END IF
                        IF (desc(istr:istr) == ')') THEN
                           at   = 0
                           istr = istr + 1
                           iat  = -1
                        ELSE
                           at = 1
                        END IF
                  END SELECT
               END DO
               100 FORMAT("OSCDFT ERROR: invalid atom index: ", I0)
            END SUBROUTINE parse_norbitals

            FUNCTION count_norbitals(desc, ioscdft) result(norbs)
               IMPLICIT NONE

               INTEGER                      :: norbs
               CHARACTER(LEN=*), INTENT(IN) :: desc
               INTEGER,          INTENT(IN) :: ioscdft
               INTEGER                      :: at ! 0: parse atom, 1: parse principal, 2: parse azimuthal
               INTEGER                      :: istr, jstr, testvar, ierr
               LOGICAL                      :: is_paren

               norbs = 0
               istr  = 1
               at    = 0
               DO WHILE (istr <= LEN_TRIM(desc))
                  SELECT CASE (at)
                     CASE (0)
                        jstr = istr
                        DO WHILE (istr <= LEN_TRIM(desc))
                           is_paren = (desc(istr:istr) == '(')
                           istr = istr + 1
                           IF (is_paren) EXIT
                        END DO
                        READ(desc(jstr:istr-2), *, IOSTAT=ierr) testvar
                        IF (ierr /= 0) THEN
                           CALL errore("oscdft_parse_orbital_desc",&
                                       "expect integer for atom number: " // desc,&
                                       ioscdft)
                        END IF
                        at = 1
                     CASE (1)
                        READ(desc(istr:istr), *, IOSTAT=ierr) testvar
                        IF (ierr /= 0) THEN
                           CALL errore("oscdft_parse_orbital_desc",&
                                       "expect integer for principal quantum number: " // desc,&
                                       ioscdft)
                        END IF
                        istr = istr + 1
                        at = 2
                     CASE (2)
                        SELECT CASE (desc(istr:istr))
                           CASE ('S')
                           CASE ('P')
                           CASE ('D')
                           CASE DEFAULT
                              CALL errore("oscdft_parse_orbital_desc",&
                                          "expect S,P,D: " // desc,&
                                          ioscdft)
                        END SELECT
                        istr = istr + 1
                        norbs = norbs + 1
                        IF (istr > LEN_TRIM(desc)) THEN
                           CALL errore("oscdft_parse_orbital_desc",&
                                       "unexpected end of descriptor: " // desc,&
                                       ioscdft)
                        END IF
                        IF (desc(istr:istr) == ')') THEN
                           at = 0
                           istr = istr + 1
                        ELSE
                           at = 1
                        END IF
                  END SELECT
               END DO
            END FUNCTION count_norbitals

      END SUBROUTINE parse_orbital_desc

      SUBROUTINE oscdft_init_indices(idx, inp)
         USE uspp_param,    ONLY : upf, nwfcm
         USE ions_base,     ONLY : ityp, nat, nsp
         USE symm_base,     ONLY : nsym, irt
         USE klist,         ONLY : nkstot, nks
         USE control_flags, ONLY : gamma_only
         IMPLICIT NONE

         TYPE(oscdft_indices_type), TARGET, INTENT(INOUT) :: idx
         TYPE(oscdft_input_type),           INTENT(INOUT) :: inp
         INTEGER, EXTERNAL                                :: set_Hubbard_l
         INTEGER                                          :: iconstr, ioscdft, iorb,&
                                                             na, nt, nwfc, n, l, i, start,&
                                                             maxl, maxn, itype, iwfc,&
                                                             isum, osum, nwfc_max, isym,&
                                                             iconstr_orbital, ioccup, m,&
                                                             idx1, idx2
         LOGICAL                                          :: largest_n, found, symmetry_error
         LOGICAL, ALLOCATABLE                             :: symmetry_error_list(:)
         TYPE(oscdft_constr_indices_type),  POINTER       :: constr
         TYPE(oscdft_orbital_indices_type), POINTER       :: orbs
         CHARACTER (len=1)  :: l_label(0:3)=(/'s','p','d','f'/)
         CHARACTER (len=7)  :: lm_label(1:7,0:3)=reshape( (/ &
            '       ','       ','       ','       ','       ','       ','       ', &
            'z      ','x      ','y      ','       ','       ','       ','       ', &
            'z2     ','xz     ','yz     ','x2-y2  ','xy     ','       ','       ', &
            'z3     ','xz2    ','yz2    ','zx2-zy2','xyz    ','x3-3xy2','3yx2-y3' /), (/7,4/) )
         CHARACTER (len=40) :: orb_label(0:3)=(/ &
            '                                        ', &
            '(z,x,y)                                 ', &
            '(z2,xz,yz,x2-y2,xy)                     ', &
            '(z3,xz2,yz2,zx2-zy2,xyz,x3-3xy2,3yx2-y3)' /)

         constr => idx%constr
         orbs   => idx%orbs

         ALLOCATE(idx%nchi(nwfcm,nsp))
         idx%nchi = 0
         ! parsing the principal quantum number from upf(nt)%els(iwfc)
         ! save to idx%nchi(iwfc,nt)
         IF (inp%debug_print) THEN
            WRITE(stdout, *) ""
            WRITE(stdout, 200)
         END IF
         DO nt=1,nsp
            IF (inp%debug_print) WRITE(stdout, 202) nt, upf(nt)%psd
            DO iwfc=1,upf(nt)%nwfc
               IF (upf(nt)%oc(iwfc) >= 0.D0) THEN
                  READ(upf(nt)%els(iwfc)(1:1), '(I1)') n
                  idx%nchi(iwfc,nt) = n
                  IF (inp%debug_print) THEN 
                     WRITE(stdout, 203) iwfc,&
                                        upf(nt)%els(iwfc),&
                                        idx%nchi(iwfc,nt),&
                                        upf(nt)%lchi(iwfc)
                  END IF
               END IF
            END DO
         END DO
         IF (inp%debug_print) WRITE(stdout, 201) 
         CALL mp_bcast(idx%nchi, ionode_id, intra_image_comm)

         ! parse inp%orbital_desc to the idx%orbs data structure
         ! example:
         ! TARGET_OCCUPATION_NUMBER
         ! T 1 4(3d) 5 0.99 0.0
         ! F 2 4(3d)
         ! F 1 5(2s2p)
         ! F 2 5(2s2p)
         !
         ! | iorb | ioscdft | orbs%iat | spin | orbital | n | l | 
         ! | 1    | 1       | 4        | UP   | 3d      | 3 | 2 |
         ! | 2    | 2       | 4        | DOWN | 3d      | 3 | 2 |
         ! | 3    | 3       | 5        | UP   | 2s      | 2 | 0 |
         ! | 4    | 3       | 5        | UP   | 2p      | 2 | 1 |
         ! | 5    | 4       | 5        | DOWN | 2s      | 2 | 0 |
         ! | 6    | 4       | 5        | DOWN | 2p      | 2 | 1 |
         ! orbs%iorb_start(1:inp%noscdft) = (1, 2, 3, 5)
         ! orbs%iorb_end  (1:inp%noscdft) = (1, 2, 4, 6)
         CALL parse_orbital_desc(idx, inp)

         ALLOCATE(orbs%iorb2ioscdft(orbs%norbs))
         DO ioscdft=1,inp%noscdft
            IF (orbs%iorb_start(ioscdft) == 0) THEN
               CALL errore("oscdft_init_indices", "orbital not specified", ioscdft)
            END IF
            DO iorb=orbs%iorb_start(ioscdft),orbs%iorb_end(ioscdft)
               na = orbs%iat(iorb)
               nt = ityp(na)

               orbs%iorb2ioscdft(iorb) = ioscdft
               ! IF (inp%debug_print) THEN
               !    WRITE(stdout, 301) na, nt, iorb, orbs%n(iorb), orbs%l(iorb)
               ! END IF
               found = .false.
               DO iwfc=1,upf(nt)%nwfc
                  IF (upf(nt)%oc(iwfc) < 0.D0) CYCLE

                  n = idx%nchi(iwfc,nt)
                  l = upf(nt)%lchi(iwfc)
                  IF (n == orbs%n(iorb) .AND. l == orbs%l(iorb)) THEN
                     found = .true.
                  END IF
               END DO
               IF (.NOT.found) THEN
                  WRITE(stdout, 301) na, nt, ioscdft, iorb, orbs%n(iorb), orbs%l(iorb)
                  CALL errore("oscdft_init_indices", "no suitable n and l pair found in upf", iorb)
               END IF
            END DO
         END DO
         CALL mp_bcast(orbs%iorb2ioscdft, ionode_id, intra_image_comm)

         idx%nconstr = COUNT(inp%constraint_applied/=CONSTR_FALSE)
         CALL mp_bcast(idx%nconstr, ionode_id, intra_image_comm)

         ALLOCATE(idx%iconstr2ioscdft(idx%nconstr))
         ALLOCATE(idx%ioscdft2iconstr(inp%noscdft))

         iconstr = 0
         DO ioscdft=1,inp%noscdft
            IF (inp%constraint_applied(ioscdft).NE.CONSTR_FALSE) THEN
               iconstr = iconstr + 1
               idx%iconstr2ioscdft(iconstr) = ioscdft
            ENDIF
         ENDDO

         idx%ioscdft2iconstr = 0
         DO iconstr=1,idx%nconstr
            ioscdft = idx%iconstr2ioscdft(iconstr)
            idx%ioscdft2iconstr(ioscdft) = iconstr
         ENDDO
         CALL mp_bcast(idx%iconstr2ioscdft, ionode_id, intra_image_comm)
         CALL mp_bcast(idx%ioscdft2iconstr, ionode_id, intra_image_comm)

         ALLOCATE(idx%ns_dim(inp%noscdft))
         DO ioscdft=1,inp%noscdft
            idx%ns_dim(ioscdft) = 0
            DO iorb=orbs%iorb_start(ioscdft),orbs%iorb_end(ioscdft)
               idx%ns_dim(ioscdft) = idx%ns_dim(ioscdft) + 2 * orbs%l(iorb) + 1
            END DO
         END DO
         idx%max_ns_dim = MAXVAL(idx%ns_dim)
         CALL mp_bcast(idx%ns_dim,     ionode_id, intra_image_comm)
         CALL mp_bcast(idx%max_ns_dim, ionode_id, intra_image_comm)

         DO ioscdft=1,inp%noscdft
            IF (inp%occup_index(ioscdft) > idx%ns_dim(ioscdft)) THEN
               CALL errore("oscdft_init_indices", "invalid inp%occup_index", ioscdft)
            END IF
            IF (inp%constraint_applied(ioscdft) /= CONSTR_FALSE) THEN
               IF (inp%occup_index(ioscdft) == 0) THEN
                  CALL errore("oscdft_init_indices", "invalid inp%occup_index", ioscdft)
               END IF
            END IF
            IF (inp%occup_index(ioscdft) == OCCUP_SUM) THEN
               IF (inp%occup_index_sum(2,ioscdft) > idx%ns_dim(ioscdft)) THEN
                  CALL errore("oscdft_init_indices", "invalid inp%occup_index", ioscdft)
               END IF
               DO isum=2,inp%occup_index_sum(1,ioscdft)
                  osum = inp%occup_index_sum(isum+1,ioscdft)
                  IF ((osum <= 0) .OR. (osum > inp%noscdft)) THEN
                     CALL errore("oscdft_init_indices", "invalid inp%occup_index_sum", ioscdft)
                  ENDIF
                  IF (inp%occup_index(osum) == 0) THEN
                     CALL errore("oscdft_init_indices", "invalid inp%occup_index", ioscdft)
                  END IF
               END DO
            END IF
         END DO

         ALLOCATE(orbs%iat_sym(nsym,orbs%norbs))
         ALLOCATE(symmetry_error_list(nat))
         symmetry_error_list = .false.
         symmetry_error = .false.
         DO iorb=1,orbs%norbs
            ioscdft = orbs%iorb2ioscdft(iorb)
            na = orbs%iat(iorb)
            IF (irt(1,na) /= na) CALL errore("oscdft_base", "internal error", 1)
            DO isym=1,nsym
               orbs%iat_sym(isym,iorb) = irt(isym,na)
            END DO
            IF (inp%constraint_applied(ioscdft) == CONSTR_FALSE) CYCLE
            IF (nsym <= 1) CYCLE
            DO isym=2,nsym
               na = orbs%iat(iorb)
               IF (orbs%iat_sym(1,iorb) /= orbs%iat_sym(isym,iorb)) THEN
                  symmetry_error = .true.
                  IF (.NOT.symmetry_error_list(na)) THEN
                     WRITE(stdout, 401) na
                     symmetry_error_list(na) = .true.
                  END IF
               END IF
            END DO
         END DO
         IF (symmetry_error) THEN
            CALL errore("oscdft_init_indices", "symmetry_error", 1)
         END IF
         DEALLOCATE(symmetry_error_list)

         constr%norbs = 0
         DO iconstr=1,idx%nconstr
            ioscdft = idx%iconstr2ioscdft(iconstr)
            IF (orbs%iorb_start(ioscdft) == 0) CYCLE
            constr%norbs = constr%norbs +&
               orbs%iorb_end(ioscdft) + 1 - orbs%iorb_start(ioscdft)
         END DO
         CALL mp_bcast(constr%norbs, ionode_id, intra_image_comm)
         ALLOCATE(constr%icorb_start(idx%nconstr),&
                  constr%icorb_end(idx%nconstr),&
                  constr%icorb2iorb(constr%norbs),&
                  constr%iorb2icorb(orbs%norbs))
         iconstr_orbital = 1
         DO iconstr=1,idx%nconstr
            constr%icorb_start(iconstr) = iconstr_orbital
            iconstr_orbital = iconstr_orbital + &
               orbs%iorb_end(ioscdft) + 1 - orbs%iorb_start(ioscdft)
            constr%icorb_end(iconstr) = iconstr_orbital - 1
         END DO
         DO iconstr=1,idx%nconstr
            ioscdft = idx%iconstr2ioscdft(iconstr)
            iorb    = orbs%iorb_start(ioscdft)
            DO iconstr_orbital=constr%icorb_start(iconstr),constr%icorb_end(iconstr)
               constr%icorb2iorb(iconstr_orbital) = iorb
               iorb = iorb + 1
            END DO
         END DO
         constr%iorb2icorb = 0
         DO iconstr_orbital=1,constr%norbs
            constr%iorb2icorb(constr%icorb2iorb(iconstr_orbital)) = iconstr_orbital
         END DO
         CALL mp_bcast(constr%icorb_start, ionode_id, intra_image_comm)
         CALL mp_bcast(constr%icorb_end,   ionode_id, intra_image_comm)
         CALL mp_bcast(constr%icorb2iorb,  ionode_id, intra_image_comm)
         CALL mp_bcast(constr%iorb2icorb,  ionode_id, intra_image_comm)

         ALLOCATE(orbs%ins2iorb(idx%max_ns_dim,inp%noscdft),&
                  orbs%ins2ioff(idx%max_ns_dim,inp%noscdft),&
                  constr%ins2iorb(idx%max_ns_dim,idx%nconstr),&
                  constr%ins2ioff(idx%max_ns_dim,idx%nconstr))
         DO ioscdft=1,inp%noscdft
            ioccup = 1
            start = orbs%iorb_start(ioscdft)
            DO iorb=start,orbs%iorb_end(ioscdft)
               na = orbs%iat(iorb)
               m = 2 * orbs%l(iorb) + 1
               DO i=1,m
                  orbs%ins2iorb(ioccup, ioscdft) = iorb
                  orbs%ins2ioff(ioccup, ioscdft) = i
                  ioccup = ioccup + 1
               END DO
            END DO
         END DO
         DO iconstr=1,idx%nconstr
            ioscdft = idx%iconstr2ioscdft(iconstr)
            DO i=1,idx%ns_dim(ioscdft)
               constr%ins2iorb(i,iconstr) = constr%iorb2icorb(orbs%ins2iorb(i,ioscdft))
               constr%ins2ioff(i,iconstr) = orbs%ins2ioff(i,ioscdft)
            END DO
         END DO
         CALL mp_bcast(orbs%ins2iorb,   ionode_id, intra_image_comm)
         CALL mp_bcast(orbs%ins2ioff,   ionode_id, intra_image_comm)
         CALL mp_bcast(constr%ins2iorb, ionode_id, intra_image_comm)
         CALL mp_bcast(constr%ins2ioff, ionode_id, intra_image_comm)

         WRITE(stdout, *) ""
         WRITE(stdout, 700)
         ! OSCDFT: ioscdft 1
         ! OSCDFT: |- 1: atom #  5; 3d(z2,xz,yz,x2-y2,xy) orbital
         ! OSCDFT: |- 1: atom #  5; 3d(z2,xz,yz,x2-y2,xy) orbital
         DO ioscdft=1,inp%noscdft
            idx1 = orbs%iorb_start(ioscdft)
            idx2 = orbs%iorb_end(ioscdft)
            WRITE(stdout, 701) ioscdft, MERGE("UP", "DW", inp%spin_index(ioscdft) == 1)
            DO iorb=idx1,idx2
               WRITE(stdout, 702) iorb, orbs%iat(iorb),&
                                  orbs%n(iorb),&
                                  l_label(orbs%l(iorb)),&
                                  TRIM(orb_label(orbs%l(iorb)))
            END DO
         END DO
         WRITE(stdout, *) ""
         ! DO ioscdft=1,inp%noscdft
         !    idx1 = orbs%iorb_start(ioscdft)
         !    idx2 = orbs%iorb_end(ioscdft)
         !    WRITE(stdout, 701) ioscdft, ((orbs%iat(iorb),&
         !                                  orbs%n(iorb),&
         !                                  l_label(orbs%l(iorb)),&
         !                                  TRIM(lm_label(m,orbs%l(iorb))),&
         !                                  m=1,2*orbs%l(iorb)+1),&
         !                                  iorb=idx1,idx2)
         ! END DO

         IF (inp%debug_print) THEN
            WRITE(stdout, 102) "iorb2isocdft", orbs%iorb2ioscdft
            WRITE(stdout, 100) "nconstr", idx%nconstr
            WRITE(stdout, 102) "iconstr2ioscdft", idx%iconstr2ioscdft
            WRITE(stdout, 102) "ioscdft2iconstr", idx%ioscdft2iconstr
            WRITE(stdout, 102) "ns_dim",     idx%ns_dim
            WRITE(stdout, 102) "max_ns_dim", idx%max_ns_dim
            WRITE(stdout, 102) "occupation index", inp%occup_index
            DO isym=1,nsym
               WRITE(stdout, 400) isym, orbs%iat_sym(isym,:)
            END DO
            WRITE(stdout, 102) "constr%norb", constr%norbs
            WRITE(stdout, 102) "constr%icorb2iorb", constr%icorb2iorb
            WRITE(stdout, 102) "constr%iorb2icorb", constr%iorb2icorb
            WRITE(stdout, 102) "constr%icorb_start", constr%icorb_start
            WRITE(stdout, 102) "constr%icorb_end", constr%icorb_end
            DO ioscdft=1,inp%noscdft
               m = idx%ns_dim(ioscdft)
               WRITE(stdout, 500) "ns_dim  ", ioscdft, idx%ns_dim(ioscdft)
               WRITE(stdout, 500) "ins2iorb", ioscdft, orbs%ins2iorb(:m,ioscdft)
               WRITE(stdout, 500) "ins2ioff", ioscdft, orbs%ins2ioff(:m,ioscdft)
            END DO
            DO iconstr=1,idx%nconstr
               ioscdft = idx%iconstr2ioscdft(iconstr)
               m = idx%ns_dim(ioscdft)
               WRITE(stdout, 500) "constr: ins2iorb", iconstr, constr%ins2iorb(:m,iconstr)
               WRITE(stdout, 500) "constr: ins2ioff", iconstr, constr%ins2ioff(:m,iconstr)
            END DO
         END IF

         IF (inp%orthogonalize_ns) THEN
            IF (gamma_only) THEN
               ALLOCATE(idx%overlap_gam(idx%max_ns_dim,idx%max_ns_dim,nsym,inp%noscdft,nkstot))
            ELSE
               ALLOCATE(idx%overlap_k  (idx%max_ns_dim,idx%max_ns_dim,nsym,inp%noscdft,nkstot))
            END IF
            ALLOCATE(idx%coeffs(idx%max_ns_dim,idx%max_ns_dim,nsym,inp%noscdft,nks))
            IF (idx%nconstr > 0) THEN
               CALL errore("oscdft_init_indices", "orthogonalize_ns with constraint not implemented", 1)
            END IF
         END IF


         100 FORMAT("OSCDFT DEBUG: ", A, ": ", I3)
         101 FORMAT("OSCDFT DEBUG: ", A, ": ", ES14.7)
         102 FORMAT("OSCDFT DEBUG: ", A, ": ", *(I3, :, " "))
         103 FORMAT("OSCDFT DEBUG: ", A, ": ", *(ES14.7, :, " "))

         200 FORMAT("OSCDFT DEBUG: ===UPF pseudopotentials===")
         201 FORMAT("OSCDFT DEBUG: ==========================")
         202 FORMAT("OSCDFT DEBUG: nt:", I3, "; psd: ", A)
         203 FORMAT("OSCDFT DEBUG: |-", I3, ": ", A2, " (n =", I3, "; l =", I3, ")")

         301 FORMAT("OSCDFT ERORR: na: ", I0, "; nt: ", I0, "; ioscdft: ", I0, "; iorb: ", I0, "; n:", I0, "; l: ", I0)

         400 FORMAT("OSCDFT DEBUG: isym: ", I3, "; atom_sym: ", *(I3, :, " "))
         401 FORMAT("OSCDFT ERROR: na: ", I5, "; must be a separate type to break cell symmetry")

         500 FORMAT("OSCDFT DEBUG: ", A, "(:,", I5, "): ", *(I5, :, " "))

         700 FORMAT("OSCDFT: Orbital Information")
         701 FORMAT("OSCDFT: ioscdft: ", I5, " spin: ", A)
         702 FORMAT("OSCDFT: |-", I3, ": atom #", I5, "; orbital ", I1, A1, A)
         ! 701 FORMAT("OSCDFT: ioscdft: ", I5, " atom+orbitals: ", *(I0, ":", I1, A1, A, :, ", "))
      END SUBROUTINE oscdft_init_indices

      SUBROUTINE oscdft_debug_nks
         USE mp_pools,  ONLY : me_pool, npool, inter_pool_comm
         USE mp_images, ONLY : intra_image_comm
         USE mp_world,  ONLY : mpime
         USE klist,     ONLY : nks
         USE mp,        ONLY : mp_barrier, mp_gather
         IMPLICIT NONE
         INTEGER :: ipool
         INTEGER, ALLOCATABLE :: nks_gather(:)

         IF (me_pool == 0) THEN
            ALLOCATE(nks_gather(npool))
            CALL mp_gather(nks, nks_gather, 0, inter_pool_comm)
            IF (ionode) THEN
               WRITE(stdout, 100) npool, nks_gather
            END If
         END IF
         100 FORMAT("OSCDFT DEBUG: k-point npool: ", I3, "; distribution: ", *(I5, :, ", "))
      END SUBROUTINE oscdft_debug_nks

      SUBROUTINE oscdft_init_context(ctx)
         USE mp_pools,         ONLY : npool, nproc_pool
         USE mp_bands,         ONLY : nbgrp, nproc_bgrp
         USE symm_base,        ONLY : nsym, d1, d2, d3
         USE control_flags,    ONLY : max_cg_iter, restart
         USE wvfct,            ONLY : nbnd
         USE klist,            ONLY : nelup, neldw
         USE noncollin_module, ONLY : npol, noncolin
         IMPLICIT NONE

         CLASS(oscdft_context_type), INTENT(INOUT), TARGET :: ctx

         INTEGER :: iconstr, na, ioscdft,row,&
                    nwordwfcO, ierr, i, curr_dim, oidx,&
                    iorb, isym, m
         TYPE(oscdft_input_type),   POINTER :: inp
         TYPE(oscdft_indices_type), POINTER :: idx

         IF (ctx%initialized) RETURN
         CALL start_clock("oscdft_init")

         inp => ctx%inp
         idx => ctx%idx

         WRITE(stdout, 701) inp%oscdft_type

         IF (.NOT.(inp%oscdft_type==1)) RETURN

         IF (npol .NE. 1) CALL errore("oscdft_init", "current value of npol not implemented", 1)
         IF (noncolin) CALL errore("oscdft_init", "noncolin not supported", 1)

         IF (inp%iteration_type == ITER_MULTIPLIERS_RHO) THEN
            WRITE(stdout, 700) inp%iteration_type, "micro-iteration of OS-CDFT multipliers inside iteration of rho"
         ELSE
            WRITE(stdout, 700) inp%iteration_type, "micro-iteration of rho inside iteration of OS-CDFT multipliers"
         END IF

         CALL oscdft_debug_nks
         CALL oscdft_init_indices(idx, inp)

         IF (ANY(inp%constraint_applied == CONSTR_LE3) .OR. ANY(inp%constraint_applied == CONSTR_GE3)) THEN
            IF (inp%convergence_type /= CONV_GRADIENT) THEN
               CALL errore("oscdft_init", "LE3 and GE3 must use gradient as convergence_type", 1)
            END IF
         END IF

         ALLOCATE(ctx%multipliers(idx%nconstr))
         ALLOCATE(ctx%old_multipliers(idx%nconstr))
         ctx%old_multipliers = 0
         ctx%multipliers = inp%initial_multipliers(idx%iconstr2ioscdft)

         CALL d_matrix(d1, d2, d3)

         ! IF (inp%debug_print) WRITE(stdout, 500) npool, nproc_pool, nbgrp, nproc_bgrp
         IF (inp%debug_print) WRITE(stdout, 200) nelup, neldw, nbnd

         ctx%conv_thr = inp%max_conv_thr
         CALL oscdft_alloc_nst(ctx%nst, idx%max_ns_dim, idx%nconstr, inp%noscdft)
         ctx%warming_up = .true.

         ctx%has_debug = ANY(inp%constraint_applied >= CONSTR_D1)

         ! IF (inp%debug_print) WRITE(stdout, 100) "iteration_type", inp%iteration_type

         IF (inp%swapping_technique.NE.OSCDFT_NONE) THEN
            DO ioscdft=1,inp%noscdft
               IF (inp%constraint_applied(ioscdft).NE.CONSTR_TRUE) CYCLE
               oidx = inp%occup_index(ioscdft)
               IF ((oidx <= 0).AND.&
                   (oidx /= OCCUP_TRACE).AND.&
                   (oidx /= OCCUP_SUM)) THEN
                   CALL errore("oscdft_init",&
                               "occup_index not implemented yet for this swapping technique",&
                                ioscdft)
               ENDIF
            ENDDO
         ENDIF

         ctx%global_start_index = 1
         max_cg_iter = 100
         IF (inp%debug_print) WRITE(stdout, 100) "max_cg_iter", max_cg_iter

         IF (restart) THEN
            CALL read_oscdft_save(ctx)
         END IF

         ctx%initialized = .true.
         CALL stop_clock("oscdft_init")

         100 FORMAT("OSCDFT DEBUG: ", A, ": ", I3)
         101 FORMAT("OSCDFT DEBUG: ", A, ": ", ES14.7)
         102 FORMAT("OSCDFT DEBUG: ", A, ": ", *(I3, " "))
         103 FORMAT("OSCDFT DEBUG: ", A, ": ", *(ES14.7, " "))
         104 FORMAT("OSCDFT SYMMETRY ", I1, ": idx: ", I5, "; vals: ", *(I1, " "))
         200 FORMAT("OSCDFT DEBUG: NELUP: ", ES14.7, "; NELDW: ", ES14.7, "; NBND: ", I5)
         500 FORMAT("OSCDFT DEBUG: npool:", I5, "; nproc_pool: ", I5,&
                    "; nbgrp: ", I5, "; nproc_bgrp:", I5)
         600 FORMAT("OSCDFT DEBUG: ", A, "(", I5, "): ", *(F8.5, " "))
         601 FORMAT("=======================================================================================")
         700 FORMAT("OSCDFT: iteration type ", I1, ": ", A)
         701 FORMAT(/5x,"OSCDFT: oscdft_type ", I1)
      END SUBROUTINE oscdft_init_context

      SUBROUTINE oscdft_alloc_nst(nst, max_ns_dim, nconstr, noscdft)
         IMPLICIT NONE

         TYPE(oscdft_ns_type), INTENT(INOUT) :: nst
         INTEGER,              INTENT(IN)    :: max_ns_dim, nconstr, noscdft

         ALLOCATE(nst%ns(max_ns_dim,max_ns_dim,noscdft),&
                  nst%numbers(noscdft),&
                  nst%eigvals(max_ns_dim,noscdft),&
                  nst%eigvects(max_ns_dim,max_ns_dim,noscdft),&
                  nst%occup_numbers(nconstr),&
                  nst%occup_eigvals(max_ns_dim,nconstr),&
                  nst%occup_eigvects(max_ns_dim,max_ns_dim,nconstr),&
                  nst%gradient(nconstr))

         nst%ns(:,:,:)             = 0.D0
         nst%numbers(:)            = 0.D0
         nst%eigvals(:,:)          = 0.D0
         nst%eigvects(:,:,:)       = 0.D0
         nst%occup_numbers(:)      = 0.D0
         nst%occup_eigvals(:,:)    = 0.D0
         nst%occup_eigvects(:,:,:) = 0.D0
         nst%gradient(:)           = 0.D0
      END SUBROUTINE oscdft_alloc_nst

      SUBROUTINE read_oscdft_save(ctx)
         USE klist,     ONLY : nelup, neldw
         USE mp,        ONLY : mp_bcast
         USE io_global, ONLY : ionode_id
         USE io_files,  ONLY : restart_dir
         IMPLICIT NONE

         TYPE (oscdft_context_type), INTENT(INOUT) :: ctx
         INTEGER, EXTERNAL                         :: find_free_unit
         INTEGER                                   :: u, dummy, curr_dim, row, ioscdft
         REAL(DP)                                  :: nelupDP, neldwDP
         REAL(DP)                                  :: energy
         LOGICAL                                   :: exst

         IF (ionode) THEN
            u = find_free_unit()
            INQUIRE(file=TRIM(restart_dir())//'/oscdft_save', exist=exst)
            IF (.NOT.exst) THEN
               CALL errore("oscdft_read_oscdft_save",&
                  "missing "//TRIM(restart_dir())//"/oscdft_save file", 1)
            END IF
            OPEN(unit=u, FILE=TRIM(restart_dir())//'/oscdft_save', status="old")

            READ(u, 101) nelupDP, neldwDP
            READ(u, 101) energy
            IF (ctx%idx%nconstr > 0) THEN
               READ(u, 101) ctx%multipliers(1:ctx%idx%nconstr)
            ELSE
               READ(u, 102) dummy
            END IF
            READ(u, 101) ctx%conv_thr
            READ(u, 102) ctx%warm_up_iter, ctx%global_start_index, ctx%multiplier_iter
            READ(u, 103) ctx%inp%get_ground_state_first

            ctx%nst%ns(:,:,:) = 0.D0
            DO ioscdft=1,ctx%inp%noscdft
               curr_dim = ctx%idx%ns_dim(ioscdft)
               DO row=1,curr_dim
                  READ(u, 101) ctx%nst%ns(row,1:curr_dim,ioscdft)
               END DO
            END DO

            CLOSE(u)
         END IF
         CALL mp_bcast(ctx%energy,            ionode_id,intra_image_comm)
         IF (ctx%idx%nconstr > 0) THEN
            CALL mp_bcast(ctx%multipliers,    ionode_id,intra_image_comm)
         END IF
         CALL mp_bcast(ctx%conv_thr,          ionode_id,intra_image_comm)
         CALL mp_bcast(ctx%warm_up_iter,      ionode_id,intra_image_comm)
         CALL mp_bcast(ctx%global_start_index,ionode_id,intra_image_comm)
         CALL mp_bcast(ctx%inp%get_ground_state_first,ionode_id,intra_image_comm)
         CALL mp_bcast(ctx%nst%ns,            ionode_id,intra_image_comm)

         IF (ionode) THEN
            WRITE(stdout, 200)
            WRITE(stdout, 201) "nel",                    nelupDP, neldwDP
            WRITE(stdout, 201) "total_energy",           energy
            IF (ctx%idx%nconstr > 0) THEN
               WRITE(stdout, 201) "multipliers",         ctx%multipliers(1:ctx%idx%nconstr)
            END IF
            WRITE(stdout, 201) "conv_thr",               ctx%conv_thr
            WRITE(stdout, 202) "warm_up_iter",           ctx%warm_up_iter
            WRITE(stdout, 202) "global_start_index",     ctx%global_start_index
            WRITE(stdout, 202) "multiplier_iter",        ctx%multiplier_iter
            WRITE(stdout, 203) "get_ground_state_first", ctx%inp%get_ground_state_first
            DO ioscdft=1,ctx%inp%noscdft
               WRITE(stdout, 202) "OCCUPATION MATRIX", ioscdft
               curr_dim = ctx%idx%ns_dim(ioscdft)
               DO row=1,curr_dim
                  WRITE(stdout, 201) "OCCUPATION MATRIX", ctx%nst%ns(row,1:curr_dim,ioscdft)
               END DO
            END DO
         END IF

         101 FORMAT(*(ES14.7))
         102 FORMAT(*(I5))
         103 FORMAT(L)

         200 FORMAT("OSCDFT: reading oscdft_save from file")
         201 FORMAT("OSCDFT: ", A, ": ", *(ES14.7, " "))
         202 FORMAT("OSCDFT: ", A, ": ", *(I5, " "))
         203 FORMAT("OSCDFT: ", A, ": ", *(L, " "))
      END SUBROUTINE read_oscdft_save
#endif
END MODULE oscdft_context
