!------------------------------------------------------------------------------!
    MODULE kinds
!------------------------------------------------------------------------------!

      USE io_global,  ONLY : stdout
      
      IMPLICIT NONE
      SAVE
! ... kind definitions
      INTEGER, PARAMETER :: dbl = selected_real_kind(14,200)
      INTEGER, PARAMETER :: sgl = selected_real_kind(6,30)
      INTEGER, PARAMETER :: i4b = selected_int_kind(9)
      INTEGER, PARAMETER :: DP = kind(0.0d0) ! compatibility with PWSCF
      PRIVATE
      PUBLIC :: i4b, dbl, sgl, DP, print_kind_info
!
!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!
!
!!   Print information about the used data types.
      SUBROUTINE print_kind_info
!
!------------------------------------------------------------------------------!
!
        IMPLICIT NONE
!
        WRITE( stdout,'(/,T2,A)') 'DATA TYPE INFORMATION:'
!
        WRITE( stdout,'(/,T2,A,T78,A,2(/,T2,A,T75,I6),3(/,T2,A,T67,E14.8))') &
          'REAL: Data type name:', 'dbl', '      Kind value:', kind(0.0_dbl), &
          '      Precision:', precision(0.0_dbl), &
          '      Smallest nonnegligible quantity relative to 1:', &
          epsilon(0.0_dbl), '      Smallest positive number:', tiny(0.0_dbl), &
          '      Largest representable number:', huge(0.0_dbl)
        WRITE( stdout,'(/,T2,A,T78,A,2(/,T2,A,T75,I6),3(/,T2,A,T67,E14.8))') &
          '      Data type name:', 'sgl', '      Kind value:', kind(0.0_sgl), &
          '      Precision:', precision(0.0_sgl), &
          '      Smallest nonnegligible quantity relative to 1:', &
          epsilon(0.0_sgl), '      Smallest positive number:', tiny(0.0_sgl), &
          '      Largest representable number:', huge(0.0_sgl)
        WRITE( stdout,'(/,T2,A,T72,A,4(/,T2,A,T61,I20))') &
          'INTEGER: Data type name:', '(default)', '         Kind value:', &
          kind(0), '         Bit size:', bit_size(0), &
          '         Largest representable number:', huge(0)
        WRITE( stdout,'(/,T2,A,T72,A,/,T2,A,T75,I6,/)') 'LOGICAL: Data type name:', &
          '(default)', '         Kind value:', kind(.TRUE.)
        WRITE( stdout,'(/,T2,A,T72,A,/,T2,A,T75,I6,/)') &
          'CHARACEER: Data type name:', '(default)', '           Kind value:', &
          kind('C')
!
      END SUBROUTINE print_kind_info
!
!------------------------------------------------------------------------------!
    END MODULE kinds
!------------------------------------------------------------------------------!
