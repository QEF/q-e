!
! Copyright (C) 2020-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------!
module upf_auxtools
!------------------------------------------------------------------------------!
  ! 
  use upf_kinds
  USE pseudo_types, ONLY : pseudo_upf
  implicit none
  private

  public :: upf_get_pp_format
  public :: upf_check_atwfc_norm

contains

!
!-----------------------------------------------------------------------
FUNCTION upf_get_pp_format (psfile) result(pp_format)
  !-----------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER :: pp_format
  CHARACTER (LEN=*) :: psfile
  INTEGER :: l
  !
  l = LEN_TRIM (psfile)
  pp_format = 5
  IF (l > 3) THEN
     IF (psfile (l-3:l) =='.xml' .OR. psfile (l-3:l) =='.XML') THEN
        pp_format = 0
     ELSE IF (psfile (l-3:l) =='.upf' .OR. psfile (l-3:l) =='.UPF') THEN
        pp_format = 1
     ELSE IF (psfile (l-3:l) =='.vdb' .OR. psfile (l-3:l) =='.van') THEN
        pp_format = 2
     ELSE IF (psfile (l-3:l) =='.gth') THEN
        pp_format = 3
     ELSE IF (psfile (l-3:l) =='.cpi' .OR. psfile (l-3:l) =='.fhi') THEN
        pp_format = 6
     ELSE IF (psfile (l-4:l) =='.cpmd') THEN
        pp_format = 7
     ELSE IF (l > 5) THEN
        If (psfile (l-5:l) =='.RRKJ3') pp_format = 4
     END IF
  END IF
  !
END FUNCTION upf_get_pp_format
!
!---------------------------------------------------------------
SUBROUTINE upf_check_atwfc_norm(upf,psfile)
  !---------------------------------------------------------------
  !  check for the presence of zero wavefunctions first
  !  check the normalization of the atomic wfc (only those with non-negative
  !  occupations) and renormalize them if the calculated norm is incorrect 
  !  by more than eps6 (10^{-6})
  !
  USE upf_kinds,    ONLY : dp
  USE upf_const,    ONLY : eps6, eps8
  USE upf_io,       ONLY : stdout
  implicit none
  type(pseudo_upf), intent(inout) :: upf
  character(LEN=*), optional, intent(in) :: psfile
  !
  integer ::             &
     mesh, kkbeta,       & ! auxiliary indices of integration limits
     l,                  & ! orbital angular momentum 
     iwfc, ir,           & ! counter on atomic wfcs and on radial mesh
     ibeta, ibeta1, ibeta2 ! counters on betas
  logical :: &
     match                 ! a logical variable 
  real(DP) :: &
     norm,               & ! the norm
     j                     ! total (spin+orbital) angular momentum
  real(DP), allocatable :: &
     work(:), gi(:)        ! auxiliary variable for becp
  character (len=80) :: renorm
  !
  allocate (work(upf%nbeta), gi(upf%mesh) )

  ! define indices for integration limits
  mesh = upf%mesh
  kkbeta = upf%kkbeta
  !
  renorm = ' '
  DO iwfc = 1, upf%nwfc
     l = upf%lchi(iwfc)
     if ( upf%has_so ) j = upf%jchi(iwfc)
     !
     ! the smooth part first ..
     gi(1:mesh) = upf%chi(1:mesh,iwfc) * upf%chi(1:mesh,iwfc)
     call simpson (mesh, gi, upf%rab, norm)
     !
     IF ( norm < eps8 ) then
        WRITE( stdout,'(5X,"WARNING: atomic wfc # ",i2, &
             & " for atom type",a," has zero norm")') iwfc, trim(upf%psd)
       !
       ! set occupancy to a small negative number so that this wfc
       ! is not going to be used for starting wavefunctions
       !
       upf%oc (iwfc) = -eps8
     END IF
     !
     IF ( upf%oc(iwfc) < 0.d0) CYCLE ! only occupied states are normalized
     !
     if (  upf%tvanp ) then
        !
        ! the US part if needed
        do ibeta = 1, upf%nbeta
           match = l.eq.upf%lll(ibeta)
           if (upf%has_so) match=match.and.abs(j-upf%jjj(ibeta)) < eps6
           if (match) then
              gi(1:kkbeta)= upf%beta(1:kkbeta,ibeta) * &
                            upf%chi (1:kkbeta,iwfc) 
              call simpson (kkbeta, gi, upf%rab, work(ibeta))
           else
              work(ibeta)=0.0_dp
           endif
        enddo
        do ibeta1=1,upf%nbeta
           do ibeta2=1,upf%nbeta
              norm=norm+upf%qqq(ibeta1,ibeta2)*work(ibeta1)*work(ibeta2)  
           enddo
        enddo
     end if
     norm=sqrt(norm)
     if (abs(norm-1.0_dp) > eps6 ) then
        renorm = TRIM(renorm) // ' ' // upf%els(iwfc)
        upf%chi(1:mesh,iwfc)=upf%chi(1:mesh,iwfc)/norm
     end if
  end do
  deallocate (work, gi )
  if ( len_trim(renorm) > 0 ) then
     if (present(psfile)) then
        write(stdout, '(5x,"file ",a,": wavefunction(s) ",a," renormalized")') &
              trim(psfile),trim(renorm)
     else
        write(stdout, '(5x,"specie ",a,": wavefunction(s) ",a," renormalized")') &
              trim(upf%psd),trim(renorm)
     endif
  endif
  return
  !
end subroutine upf_check_atwfc_norm
  !
end module
