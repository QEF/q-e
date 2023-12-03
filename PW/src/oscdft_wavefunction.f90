MODULE oscdft_wavefunction
#if defined (__OSCDFT)
   USE kinds, ONLY : DP

   PRIVATE
   PUBLIC oscdft_wavefunction_type, check_bec_type_unallocated,&
          check_bec_type_unallocated_gpu

   TYPE oscdft_wavefunction_type
      INTEGER                  :: iun = 0, n, nword
      INTEGER,     ALLOCATABLE :: offset(:,:)
      INTEGER,     ALLOCATABLE :: source(:,:)
      COMPLEX(DP), ALLOCATABLE :: wfc(:,:)
      LOGICAL                  :: constr_index
      LOGICAL                  :: use_sym
      LOGICAL                  :: have_S ! if .true. then wfc has the S matrix (beta projectors) applied
      LOGICAL                  :: initialized = .false.

      CONTAINS
         PROCEDURE :: get_offset_ioscdft => get_offset_ioscdft
         PROCEDURE :: get_offset_iconstr => get_offset_iconstr
         GENERIC :: get_offset => get_offset_ioscdft, get_offset_iconstr
   END TYPE oscdft_wavefunction_type

   CONTAINS
      FUNCTION get_offset_ioscdft(this, orbs, m, ioscdft, isym) RESULT(offset)
         USE oscdft_indices, ONLY : oscdft_orbital_indices_type
         IMPLICIT NONE
         INTEGER :: offset
         CLASS(oscdft_wavefunction_type),   INTENT(IN) :: this
         TYPE(oscdft_orbital_indices_type), INTENT(IN) :: orbs
         INTEGER,                           INTENT(IN) :: m, ioscdft
         INTEGER,                           INTENT(IN) :: isym

         INTEGER :: iorb, ioff, isym_
         LOGICAL :: isym_present

         iorb = orbs%ins2iorb(m,ioscdft)
         ioff = orbs%ins2ioff(m,ioscdft)

         isym_present = isym /= -1

         IF (this%use_sym .AND. .NOT.isym_present) THEN
            CALL errore("oscdft_get_offset", "internal error need to use sym", 1)
         ELSE IF (.NOT.this%use_sym .AND. isym_present) THEN
            CALL errore("oscdft_get_offset", "internal error no sym information", 1)
         END IF
         IF (this%constr_index) CALL errore("oscdft_get_offset", "use ioscdft index", 1)

         isym_ = MERGE(isym, 1, isym_present)
         offset = this%offset(isym_,iorb) + ioff
         RETURN
      END FUNCTION get_offset_ioscdft

      FUNCTION get_offset_iconstr(this, constr, m, iconstr, isym) RESULT(offset)
         USE oscdft_indices, ONLY : oscdft_constr_indices_type
         IMPLICIT NONE
         INTEGER :: offset
         CLASS(oscdft_wavefunction_type),  INTENT(IN) :: this
         TYPE(oscdft_constr_indices_type), INTENT(IN) :: constr
         INTEGER,                          INTENT(IN) :: m, iconstr
         INTEGER,                          INTENT(IN) :: isym

         INTEGER :: iorb, ioff, isym_
         LOGICAL :: isym_present

         iorb = constr%ins2iorb(m,iconstr)
         ioff = constr%ins2ioff(m,iconstr)

         isym_present = isym /= -1

         IF (this%use_sym .AND. .NOT.isym_present) THEN
            CALL errore("oscdft_get_offset", "internal error need to use sym", 1)
         ELSE IF (.NOT.this%use_sym .AND. isym_present) THEN
            CALL errore("oscdft_get_offset", "internal error no sym information", 1)
         END IF
         IF (.NOT.this%constr_index) CALL errore("oscdft_get_offset", "use iconstr index", 1)

         isym_ = MERGE(isym, 1, isym_present)
         offset = this%offset(isym_,iorb) + ioff
         RETURN
      END FUNCTION get_offset_iconstr

      SUBROUTINE check_bec_type_unallocated(bec)
         USE becmod, ONLY : bec_type
         IMPLICIT NONE

         TYPE(bec_type), INTENT(INOUT) :: bec

         IF (ALLOCATED(bec%r)) THEN
            CALL errore("check_bec_type_unallocated", "bec%r allocated", 1)
         END IF
         IF (ALLOCATED(bec%k)) THEN
            CALL errore("check_bec_type_unallocated", "bec%k allocated", 1)
         END IF
      END SUBROUTINE check_bec_type_unallocated
#endif
END MODULE oscdft_wavefunction
