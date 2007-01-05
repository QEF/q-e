!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE stre
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE 
  SAVE
  !
  REAL(DP) :: stress(3,3)
  !
END MODULE stre
!
!----------------------------------------------------------------------------
MODULE dqrad_mod
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE 
  SAVE
  !
  REAL(DP), ALLOCATABLE :: dqrad(:,:,:,:,:,:)
  !
  CONTAINS
  !
  SUBROUTINE deallocate_dqrad_mod()
    !
    IF ( ALLOCATED( dqrad ) ) DEALLOCATE( dqrad )
    !
  END SUBROUTINE deallocate_dqrad_mod
  !
END MODULE dqrad_mod
!
!----------------------------------------------------------------------------
module betax
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE 
  SAVE
  !
  INTEGER              :: mmx = 5000
  REAL(DP)             :: refg
  REAL(DP),ALLOCATABLE :: betagx(:,:,:), dbetagx(:,:,:), &
                          qradx(:,:,:,:), dqradx(:,:,:,:)
  !
  CONTAINS
  !
  SUBROUTINE deallocate_betax()
    !
    IF ( ALLOCATED( betagx ) )  DEALLOCATE( betagx )
    IF ( ALLOCATED( dbetagx ) ) DEALLOCATE( dbetagx )
    IF ( ALLOCATED( qradx ) )   DEALLOCATE( qradx )
    IF ( ALLOCATED( dqradx ) )  DEALLOCATE( dqradx )
    !
  END SUBROUTINE deallocate_betax
  !
END MODULE betax
!
!----------------------------------------------------------------------------
MODULE cpr_subroutines
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  SAVE
  !
  CONTAINS
  !
  subroutine compute_stress( stress, detot, h, omega )
    real(8) :: stress(3,3), detot(3,3), h(3,3), omega
    integer :: i, j
         do i=1,3
            do j=1,3
               stress(i,j)=-1.d0/omega*(detot(i,1)*h(j,1)+              &
     &                      detot(i,2)*h(j,2)+detot(i,3)*h(j,3))
            enddo
         enddo
    return
  end subroutine compute_stress

  subroutine print_atomic_var( var, na, nsp, head, iunit )
    use io_global, only: stdout
    real(8) :: var(:,:)
    integer :: na(:), nsp
    integer, optional :: iunit
    character(len=*), optional :: head
    integer :: i, ia, is, iu, isa
    if( present( iunit ) ) then
      iu = iunit
    else
      iu = stdout
    end if
    if( present( head ) ) then 
      WRITE( iu,*) head
    end if
    isa = 0
    DO is = 1, nsp
      DO ia = 1, na(is)
        isa = isa + 1
        WRITE( iu,'(3f14.8)') ( var(i,isa), i=1, 3 )
      END DO
    END DO
    return
  end subroutine print_atomic_var

  subroutine print_cell_var( var, head, iunit )
    use io_global, only: stdout
    real(8) :: var(3,3)
    integer, optional :: iunit
    character(len=*), optional :: head
    integer :: i, j, iu
    if( present( iunit ) ) then
      iu = iunit
    else
      iu = stdout
    end if
    if( present( head ) ) then 
      WRITE( iu,*)
      WRITE( iu,*) head
      WRITE( iu, 5555 ) ((var(i,j),j=1,3),i=1,3)
 5555    format(1x,f12.5,1x,f12.5,1x,f12.5/                             &
     &          1x,f12.5,1x,f12.5,1x,f12.5/                             &
     &          1x,f12.5,1x,f12.5,1x,f12.5//)
    else
      write(iu,3340) ((var(i,j),i=1,3),j=1,3)
 3340     format(9(1x,f9.5))
    end if
    return
  end subroutine print_cell_var
  !
  !--------------------------------------------------------------------------
  SUBROUTINE ions_cofmsub( tausp, iforce, na, nsp, cdm, cdm0 )
    !--------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(INOUT) :: tausp(:,:)
    INTEGER,        INTENT(IN)    :: iforce(:,:)
    INTEGER,        INTENT(IN)    :: na(:), nsp
    REAL(DP), INTENT(IN)    :: cdm(:), cdm0(:)
    !
    INTEGER :: i, ia, is, isa
    !
    !
    isa = 0
    !
    DO is = 1, nsp
       !
       DO ia = 1, na(is)
          !
          isa = isa + 1
          !
          DO i = 1, 3
             !
             tausp(i,isa) = tausp(i,isa) + &
                            DBLE( iforce(i,isa) ) * ( cdm0(i) - cdm(i) )
             !
          END DO
          !
       END DO
       !
    END DO
    !
    RETURN
    !
  END SUBROUTINE ions_cofmsub
  !

 
  SUBROUTINE print_lambda( lambda, n, nshow, ccc, iunit )
    USE io_global,         ONLY: stdout, ionode
    USE cp_main_variables, ONLY: collect_lambda, descla
    USE electrons_base,    ONLY: nudx
    IMPLICIT NONE
    real(8), intent(in) :: lambda(:,:,:), ccc
    integer, intent(in) :: n, nshow
    integer, intent(in), optional :: iunit
    integer :: nnn, j, un, i, is
    real(8), allocatable :: lambda_repl(:,:)
    if( present( iunit ) ) then
      un = iunit
    else
      un = stdout
    end if
    nnn = min( nudx, nshow )
    ALLOCATE( lambda_repl( nudx, nudx ) )
    IF( ionode ) WRITE( un,*)
    DO is = 1, SIZE( lambda, 3 )
       CALL collect_lambda( lambda_repl, lambda(:,:,is), descla(:,is) )
       IF( ionode ) THEN
          WRITE( un,3370) '    lambda   nudx, spin = ', nudx, is
          IF( nnn < n ) WRITE( un,3370) '    print only first ', nnn
          DO i=1,nnn
             WRITE( un,3380) (lambda_repl(i,j)*ccc,j=1,nnn)
          END DO
       END IF
    END DO
    DEALLOCATE( lambda_repl )
3370   FORMAT(26x,a,i4)
3380   FORMAT(9f8.4)
    RETURN
  END SUBROUTINE print_lambda

   subroutine add_thermal_stress( stress, pmass, omega, h, vels, nsp, na )
     real(8) :: stress(3,3)
     real(8), intent(in) :: pmass(:), omega, h(3,3), vels(:,:)
     integer, intent(in) :: nsp, na(:)
     integer :: i, j, is, ia, isa
     isa = 0
     do is=1,nsp
       do ia=1,na(is)
       isa = isa + 1
         do i=1,3
           do j=1,3
             stress(i,j)=stress(i,j)+pmass(is)/omega*           &
      &        ((h(i,1)*vels(1,isa)+h(i,2)*vels(2,isa)+    &
      &          h(i,3)*vels(3,isa))*(h(j,1)*vels(1,isa)+  &
      &          h(j,2)*vels(2,isa)+h(j,3)*vels(3,isa)))
           enddo
         enddo
       enddo
     enddo
     return
   end subroutine add_thermal_stress

END MODULE cpr_subroutines
