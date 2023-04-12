!
! Copyright (C) 2023 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file 'License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
PROGRAM upftest
  !---------------------------------------------------------------------
  !
  !     Pseudopotential testing utility
  !
  USE pseudo_types, ONLY : pseudo_upf, deallocate_pseudo_upf
  USE casino_pp,    ONLY : conv_upf2casino, write_casino_tab
  USE write_upf_new,ONLY : write_upf
  !
  IMPLICIT NONE
  TYPE(pseudo_upf) :: upf
  INTEGER :: nargs, prefix_len
  CHARACTER(LEN=256) :: filein, fileout
  !
  nargs = command_argument_count()
  if ( nargs < 1 ) then
     write(6,"('pseudopotential file > ')",ADVANCE="no") 
     read (5,"(A)") filein
  else  if ( nargs == 1 ) then
     CALL get_command_argument(1, filein)
  else  if ( nargs > 1 ) then
     write(6,*) 'Usage: upftest "pseudopotential file"'
     stop
  end if
  !
  CALL read_ps ( filein, upf )
  !
  IF ( INDEX(TRIM(filein),'.UPF' ) > 0) THEN 
     prefix_len = INDEX(TRIM(filein),'.UPF') - 1
  ELSE IF (INDEX(TRIM(filein),'.upf') > 0 ) THEN
     prefix_len = INDEX(TRIM(filein),'.upf') - 1
  ELSE 
     prefix_len = LEN_TRIM(filein)
  ENDIF
  fileout = filein(1:prefix_len)
  !
  CALL test_ps_loc ( fileout, upf )
  CALL test_ps_beta ( fileout, upf )
  !
END PROGRAM upftest
!
SUBROUTINE test_ps_loc ( fileout, upf )
  !
  USE pseudo_types, ONLY : pseudo_upf
  !
  USE upf_kinds,  ONLY : dp
  USE upf_const,  ONLY : fpi, e2
  !
  IMPLICIT NONE
  !
  TYPE(pseudo_upf), intent(in) :: upf
  CHARACTER(LEN=256) :: fileout
  !
  REAL(DP) :: dq = 0.01_dp
  INTEGER, PARAMETER :: nq = 4000
  REAL(DP) :: q(nq), vq(nq), vq2(nq), aux(upf%mesh), aux1(upf%mesh)
  REAL(DP) :: vq0, vq1, vlq
  INTEGER :: iq, ir, msh
  !
  open( unit = 100, file=trim(fileout) // '-vofr.dat', status='unknown', &
       form='formatted' )
  write(100,'("#   r,  vloc(r), vloc+Ze^2/r")') 
  do ir=1,upf%mesh
     write(100,'(f10.4,2f25.12)') upf%r(ir), upf%vloc(ir), &
          upf%vloc(ir) + upf%zp*e2/max(0.1,upf%r(ir))
  end do
  close(unit=100)

  do iq = 1, nq
     q(iq) = (iq-1)*dq
  end do
  !
  ! Compute vloc(q=0) for the local potential without the divergence
  !
  do ir = 1, upf%mesh
     aux (ir) = upf%r(ir) * (upf%r(ir)*upf%vloc(ir) + upf%zp * e2)
     if ( upf%r(ir) < 10 ) msh = ir
  enddo
  !
  print '("msh, mesh = ",2i5)',msh,upf%mesh
  call simpson (upf%mesh, aux, upf%rab, vq0)
  call simpson (msh, aux, upf%rab, vq1)
  !
  ! Compute q^2*vloc(q). We first store the part of the integrand 
  !   function independent of q in real space
  !
  do ir = 1, upf%mesh
     aux1 (ir) = upf%r(ir) * upf%vloc(ir) + upf%zp * e2 * erf (upf%r(ir) )
  enddo
  !
  !    and here we perform the integral, after multiplication
  !    with the q-dependent part
  !
  do iq = 1, nq
     do ir = 1, upf%mesh
        aux(ir) = aux1(ir) * sin (q(iq)*upf%r(ir)) * q(iq)
     enddo
     call simpson (upf%mesh, aux, upf%rab, vlq)
     !
     !   here we re-add the analytic fourier transform of the erf function
     !
     vq(iq) = vlq  - e2 * upf%zp * exp ( - q(iq)**2/4.0_dp) !!! / q(iq)**2
     !
     call simpson (msh, aux, upf%rab, vlq)
     vq2(iq)= vlq  - e2 * upf%zp * exp ( - q(iq)**2/4.0_dp) !!! / q(iq)**2
  enddo
  !!! vq(:) = vq(:) * fpi / omega
  !
  open( unit = 100, file=trim(fileout) // '-vofq.dat', status='unknown', &
       form='formatted' )
  write(100,'("#   q,  q^2*v(q) ((old)  a^2*v(q) (cp90)  v(0) = ",2f15.8)') &
       vq0,vq1
  do iq=1,nq
     write(100,'(f10.4,2f25.12)') q(iq), vq(iq), vq2(iq)
  end do
  close(unit=100)
  !  
END SUBROUTINE test_ps_loc
!
SUBROUTINE test_ps_beta ( fileout, upf )
  !
  USE pseudo_types, ONLY : pseudo_upf
  !
  USE upf_kinds,  ONLY : dp
  USE upf_const,  ONLY : fpi, e2
  !
  IMPLICIT NONE
  !
  TYPE(pseudo_upf), intent(in) :: upf
  CHARACTER(LEN=256) :: fileout
  !
  REAL(DP) :: dq = 0.01_dp
  INTEGER, PARAMETER :: nq = 4000
  REAL(DP) :: q(nq), aux(upf%mesh), aux1(upf%mesh)
  REAL(DP) :: vq0, vq1
  INTEGER :: iq, ir, nb, iun
  character(len=256) :: filename
  character(len=2) :: number
  !
  print '("kkbeta = ",1i5)',upf%kkbeta
  do nb=1,upf%nbeta
     write(filename,'(A,"-betar-",i1,".dat")') trim(fileout),nb
     iun = 199 + nb
     open( unit=iun, file=filename, status='unknown', form='formatted' )
     write(iun,'("#   r,  vbeta(r,lm)")') 
     do ir=1,upf%kkbeta
        write(iun,'(f10.4,f25.12)') upf%r(ir), upf%beta(ir,nb)
     end do
     close(unit=iun)
  end do
  
  do iq = 1, nq
     q(iq) = (iq-1)*dq
  end do
  !
  ! Compute q^2*beta(q)
  !
  do nb=1,upf%nbeta
     write(filename,'(A,"-betaq-",i1,".dat")') trim(fileout),nb
     iun = 299 + nb
     open( unit=iun, file=filename, status='unknown', form='formatted' )
     write(iun,'("#   q,  q^2*v(q) ((old)  a^2*v(q) (cp90)")')
     do iq = 1, nq
        do ir = 1, upf%kkbeta
           aux(ir) = upf%beta(ir,nb) * sin (q(iq)*upf%r(ir)) * q(iq)
        enddo
        call simpson (upf%kkbeta, aux, upf%rab, vq0)
        call simpson_cp90 (upf%kkbeta, aux, upf%rab, vq1)
        write(iun,'(f10.4,2f25.12)') q(iq), vq0, vq1
     end do
     close(unit=iun)
  enddo
  !
END SUBROUTINE test_ps_beta
