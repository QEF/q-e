!
! Copyright (C) 2001-2014 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
!
! This routine reads lambda*.dat and compute a^2F, phonon DOS, lambda,
! & omega_ln 
!
!---------------------------------------------------------------------
MODULE alpha2f_vals
  !-------------------------------------------------------------------
  !
  ! This MODULE contains global variables for alpha2f.x 
  ! 
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER :: nfreq
  !
  REAL(DP),ALLOCATABLE,SAVE :: &
  & omg(:,:), & ! (nmodes,nqs) Phonon frequencies on irreducible q
  & lam(:,:), & ! (nmodes,nqs) El-Ph coupling  on irreducible q
  & pol(:,:,:)  ! (nmodes/3,nmodes,nqs) 
  !
END MODULE alpha2f_vals
!
!----------------------------------------------------------------
MODULE alpha2f_routines
  !--------------------------------------------------------------
  !
  ! This module contains SUBROUTINEs for alpha2f.x
  !
  IMPLICIT NONE
  !
CONTAINS
!
!-------------------------------------------------------------------------
SUBROUTINE read_polarization()
  !-----------------------------------------------------------------------
  !
  ! This routine read the polarization vectors 
  ! from [prefix].dyn* & lambda*.dat
  !
  USE kinds,      ONLY : DP
  USE disp,       ONLY : nqs
  USE modes,      ONLY : nmodes
  USE output,     ONLY : fildyn
  USE ions_base,  ONLY : nat, amass, ityp
  USE io_dyn_mat, ONLY : read_dyn_mat_param, read_dyn_mat_tail
  USE control_ph, ONLY : xmldyn
  USE constants,  ONLY : AMU_RY
  USE symm_base,  ONLY : irt, nsym
  !  
  USE alpha2f_vals, ONLY : pol
  !
  IMPLICIT NONE
  !
  INTEGER :: iq, imode, fi, iat, idummy1, idummy2, isym
  REAL(DP) :: norm, rdummy(nmodes), pol_sym(nmodes,nqs)
  COMPLEX(DP) :: cpol(3,nat, nmodes)
  CHARACTER :: star
  CHARACTER(10) :: ctmp
  CHARACTER(256) :: filin
  !
  CHARACTER(LEN=6) :: int_to_char
  INTEGER, EXTERNAL :: find_free_unit
  !
  ALLOCATE(pol(nat,nmodes,nqs))
  !
  IF(.NOT. xmldyn) fi = find_free_unit()
  !
  DO iq = 1, nqs
     !
     IF(xmldyn) THEN
        filin = TRIM(fildyn) // TRIM( int_to_char( iq ) )
        CALL read_dyn_mat_param(filin,idummy1,idummy2)
        CALL read_dyn_mat_tail(nat,rdummy,cpol)
     ELSE
        !
        OPEN(fi, file = TRIM(fildyn) // TRIM(int_to_char(iq)))
        !
        DO imode = 1, 100000
           !
           READ(fi,*) star
           !
           IF (star == "*") then
              exit
           END IF
           !
        END DO
        !
        DO imode = 1, nmodes
           !
           READ(fi,*) star
           !
           DO iat = 1, nat
              READ(fi,'(a2,6f10.6)') ctmp, cpol(1:3, iat, imode)
           END DO
           !
        END DO
        !
        CLOSE(fi)
        !
     END IF !(.NOT. xmldyn)
     !
     DO imode = 1, nmodes
        !
        DO iat = 1, nat
           pol(iat,imode,iq) = REAL(DOT_PRODUCT(cpol(1:3,iat,imode), &
           &                                    cpol(1:3,iat,imode)), DP) &
           &                 * amass(ityp(iat)) * AMU_RY
        END DO
        !
        ! .. Normalize Pol
        !
        norm = SUM(pol(1:nat,imode,iq))
        pol(1:nat,imode,iq) = pol(1:nat,imode,iq) / norm
        !
     END DO
     !
  END DO ! iq
  !
  ! .. Symmetrize polaization vector
  !
  DO isym = 1, nsym
     !
     DO iat = 1, nat
        pol_sym(1:nmodes,1:nqs) = ( pol(         iat ,1:nmodes,1:nqs) &
        &                         + pol(irt(isym,iat),1:nmodes,1:nqs) ) * 0.5_dp
        pol(         iat ,1:nmodes,1:nqs) = pol_sym(1:nmodes,1:nqs)
        pol(irt(isym,iat),1:nmodes,1:nqs) = pol_sym(1:nmodes,1:nqs)
     END DO
     !
  END DO
  !
END SUBROUTINE read_polarization
!
!--------------------------------------------------------------------
SUBROUTINE read_lam()
  !------------------------------------------------------------------
  !  
  ! This routine reads lambad_{q nu} & omega_{q nu} from lambda*.dat
  !
  USE kinds,     ONLY : DP
  USE disp,      ONLY : nqs, x_q
  USE cell_base, ONLY : at
  USE modes,     ONLY : nmodes
  USE output,    ONLY : fildyn
  !
  USE alpha2f_vals, ONLY : omg, lam
  !
  IMPLICIT NONE
  !
  INTEGER :: iq, fi = 10, im
  REAL(DP) :: x_q0(3)
  character(100) :: ctmp
  !
  CHARACTER(LEN=6) :: int_to_char
  INTEGER, EXTERNAL :: find_free_unit
  !
  ALLOCATE(omg(nmodes,nqs), lam(nmodes,nqs))
  !
  fi = find_free_unit()
  !
  DO iq = 1, nqs
     !
     OPEN(fi, file = TRIM(fildyn) // TRIM(int_to_char(iq)) // ".elph." // TRIM(int_to_char(iq)))
     !
     READ(fi,*) x_q0(1:3), nmodes
     !
     IF (ALL(ABS(x_q0(1:3) - x_q(1:3,iq)) > 1d-8)) THEN
        WRITE(*,*) "iq : ", iq
        WRITE(*,*) "q(file) : ", x_q0(1:3)
        WRITE(*,*) "q(grid) : ", x_q(1:3,iq)
        STOP "read_lam : qv "
     END IF
     !
     READ(fi,*) omg(1:nmodes,iq)
     !
     READ(fi,*) ctmp
     READ(fi,*) ctmp
     !
     DO im = 1, nmodes
        READ(fi,'(a19,f8.4)') ctmp, lam(im,iq)
        IF (omg(im,iq) > 0.0_dp) THEN
           omg(im,iq) = SQRT(omg(im,iq))
        ELSE
           omg(im,iq) = 1.0_dp
           lam(im,iq) = 0.0_dp
        END IF
     END DO
     !
     CLOSE(fi)
     !
  END DO
  !
END SUBROUTINE read_lam
!
!-------------------------------------------------------------------
SUBROUTINE compute_a2F()
  !-----------------------------------------------------------------
  !
  ! This routine writes a2F and phonon DOS to file (a2F.dat).
  !
  USE kinds,     ONLY : DP
  USE ions_base, ONLY : nat
  USE modes,     ONLY : nmodes
  USE ktetra,    ONLY : ntetra, tetra, opt_tetra_init, opt_tetra_partialdos
  USE wvfct,     ONLY : et, nbnd
  USE io_global, ONLY : stdout
  USE cell_base, ONLY : at, bg
  USE klist,     ONLY : nkstot
  USE disp,      ONLY : nq1, nq2, nq3, nqs, x_q
  USE symm_base, ONLY : nsym, s, time_reversal, t_rev
  USE lsda_mod,  ONLY : nspin
  USE constants, ONLY : rytoev
  USE io_files,  ONLY : prefix
  !
  USE alpha2f_vals, ONLY : omg, pol, lam, nfreq
  !
  IMPLICIT NONE
  !
  INTEGER :: ie, fo
  REAL(DP) :: DeltaE, lambda, omglog, freq, proj(1+nat,nmodes, nqs)
  LOGICAL :: kresolveddos
  REAL(DP),ALLOCATABLE :: dos(:), pdos(:,:)
  !
  INTEGER, EXTERNAL :: find_free_unit
  !
  WRITE(stdout,*) ""
  WRITE(stdout,'(a)') "   Calculation of alpha2F"
  WRITE(stdout,*) ""
  !
  DeltaE = MAXVAL(omg(1:nmodes,1:nqs)) / REAL(nfreq, DP)
  WRITE(stdout,*) "    Number of Frequencies : ", nfreq
  WRITE(stdout,*) "    Frequency Step [Ry] : ", DeltaE
  !
  ALLOCATE(dos(0:nfreq), pdos(0:nfreq,nat+1))
  dos(0:nfreq) = 0.0_dp
  pdos(0:nfreq,1:nat+1) = 0.0_dp
  !
  ntetra = nq1 * nq2 * nq3 * 6
  CALL opt_tetra_init(nsym, s, time_reversal, t_rev, at, bg, nqs, &
  &                1, 1, 1, nq1, nq2, nq3, nqs, x_q, 1)
  !
  IF(ALLOCATED(et)) DEALLOCATE(et)
  ALLOCATE(et(nmodes,nqs))
  proj(1:nat, 1:nmodes, 1:nqs) = pol(1:nat, 1:nmodes, 1:nqs)
  proj(1+nat, 1:nmodes, 1:nqs) = lam(       1:nmodes, 1:nqs) &
  &                            * omg(       1:nmodes, 1:nqs) * 0.5_dp 
  et(         1:nmodes, 1:nqs) = omg(       1:nmodes, 1:nqs)
  !
  kresolveddos = .FALSE.
  nspin = 1
  nbnd = nmodes
  nkstot = nqs
  CALL opt_tetra_partialdos(1, kresolveddos, nfreq, nat + 1, 1, &
  &                         0.0_dp, DeltaE, proj, pdos, dos, 1)
  !
  dos( 0:nfreq         ) = dos( 0:nfreq         ) * rytoev
  pdos(0:nfreq, 1:nat+1) = pdos(0:nfreq, 1:nat+1) * rytoev
  !
  ! .. Write a2F, dos, pdos two a file
  !
  WRITE(stdout,*)
  WRITE(stdout,'(a,a)') "     Writing a2F to a file ", TRIM(prefix) // ".a2F.dat"
  WRITE(stdout,*)
  !
  fo = find_free_unit()
  OPEN(fo, file = TRIM(prefix) // ".a2F.dat")
  !
  WRITE(fo,*) "# Frequency[Ry], a2F, DOS[Ry], Partial DOS(Atom 1), PDOS(Atom 2), .."
  !
  lambda = 0.0_dp
  omglog = 0.0_dp
  DO ie = 1, nfreq
     freq = REAL(ie, dp) * DeltaE
     WRITE(fo,'(200e25.15)') freq, pdos(ie, 1+nat), dos(ie), pdos(ie, 1:nat)
     lambda = lambda + DeltaE * 2.0_dp * pdos(ie, 1+nat) / freq
     omglog = omglog + DeltaE * 2.0_dp * pdos(ie, 1+nat) / freq * LOG(freq)
  END DO
  !
  CLOSE(fo)
  !
  omglog = omglog / lambda
  omglog = EXP(omglog)
  !
  WRITE(stdout,*) "    Compute lambda and omega_ln from a2F to verify it."
  WRITE(stdout,*) "             lambda : ", lambda
  WRITE(stdout,*) "      omega_ln [Ry] : ", omglog
  !
END SUBROUTINE compute_a2F
!
!-----------------------------------------------------------------
SUBROUTINE compute_lambda()
  !---------------------------------------------------------------
  !
  ! This routine computes omega_ln & lambda
  !
  USE kinds,     ONLY : DP
  USE modes,     ONLY : nmodes  
  USE constants, ONLY : RY_TO_THZ, K_BOLTZMANN_RY
  USE ktetra,    ONLY : ntetra, wlsm, tetra
  USE disp,      ONLY : nqs
  USE io_global, ONLY : stdout
  USE io_files,  ONLY : prefix
  !
  USE alpha2f_vals, ONLY : omg, lam
  !
  IMPLICIT NONE
  !
  INTEGER :: ii, itetra, iq, fo
  REAL(DP) :: omglog, lambda, wq(nqs)
  !
  INTEGER, EXTERNAL :: find_free_unit
  !
  wq(1:nqs) = 0.0_dp
  DO itetra = 1, ntetra
     DO ii = 1, 20
        wq(tetra(ii, itetra)) = wq(tetra(ii, itetra)) + SUM(0.25_dp * wlsm(1:4,ii))
     END DO
  END DO
  wq(1:nqs) = wq(1:nqs) / REAL(ntetra, DP)
  !
  lambda = 0.0_dp
  omglog = 0.0_dp
  DO iq = 1, nqs
     lambda = lambda + wq(iq) * SUM(lam(1:nmodes, iq))
     omglog = omglog + wq(iq) * SUM(LOG(omg(1:nmodes,iq)) * lam(1:nmodes,iq))
  END DO
  omglog = omglog / lambda
  omglog = EXP(omglog)
  !
  WRITE(stdout,*) ""
  WRITE(stdout,'(a)') "   Compute lambda and omega_ln directly from omega_q and lambda_q"
  WRITE(stdout,*) ""
  WRITE(stdout,*) "            lambda : ", lambda
  WRITE(stdout,*) "    omega_ln [Ryd] : ", omglog
  WRITE(stdout,*) "    omega_ln [THz] : ", omglog * RY_TO_THZ
  WRITE(stdout,*) "      omega_ln [K] : ", omglog / K_BOLTZMANN_RY
  !
  ! Tc from McMillan's formula
  !
  fo = find_free_unit()
  OPEN(fo, file = TRIM(prefix) // ".McMillan.gp")
  WRITE(stdout,*)
  WRITE(stdout,'(a,a)') "   For plotting T_c from the McMillan formula, please type"
  WRITE(stdout,*) "    $ gnuplot ", TRIM(prefix) // ".McMillan.gp"
  !
  WRITE(fo,'(a)') 'unset key'
  WRITE(fo,'(a)') 'set xlabel "mu*"'
  WRITE(fo,'(a)') 'set ylabel "T_c [K]"'
  WRITE(fo,'(a,e15.5,a)') "set xrange [0.0:", lambda / (1.0_dp + 0.62_dp * lambda), "]"
  WRITE(fo,'(a,e15.5,a,e15.5,a,e15.5,a,e15.5,a)') &
  & "plot ", omglog / K_BOLTZMANN_RY / 1.2_dp, &
  & "*exp(", -1.04_dp *(1.0_dp + lambda), "/(", lambda, "-x*", 1.0_dp + 0.62_dp * lambda, "))"
  WRITE(fo,'(a)') 'pause -1'
  !tc = omglog * 1.5788737d5 / 1.2_dp * &
  !exp(-1.04_dp *(1.0_dp + lambda)/(lambda - mustar * (1.0_dp + 0.62_dp * lambda)) )
  !
  CLOSE(fo)
  !
END SUBROUTINE compute_lambda
!
END MODULE alpha2f_routines
!
!--------------------------------------------------------------------------------
PROGRAM alpha2f
  !------------------------------------------------------------------------------
  !
  ! This routine reads lambda*.dat and compute a^2F, phonon DOS, lambda,
  ! & omega_ln 
  !
  USE mp_global,      ONLY : mp_startup, mp_global_end
  USE environment,    ONLY : environment_start, environment_end
  USE elph_tetra_mod, ONLY : in_alpha2f
  USE io_global,      ONLY : ionode
  !
  USE alpha2f_vals,     ONLY : nfreq
  USE alpha2f_routines, ONLY : read_lam, compute_a2f, compute_lambda, read_polarization
  !
  implicit none
  !
  CHARACTER (LEN=256) :: auxdyn
  !
  NAMELIST /INPUTA2F/ nfreq
  !
#if defined(__MPI)
  CALL mp_startup()
#endif
  CALL environment_start('ALPHA2F')
  in_alpha2f = .TRUE.
  !
  CALL phq_readin()
  !
  IF(ionode) READ( 5, INPUTA2F)
  !
  CALL check_initial_status(auxdyn)
  !
  IF(ionode) THEN
     !
     CALL read_polarization()
     CALL read_lam()
     !
     CALL compute_a2f()
     !
     CALL compute_lambda()
     !
  END IF
  !
  CALL environment_end('ALPHA2F')
#if defined(__MPI)
  CALL mp_global_end()
#endif
  !
END PROGRAM alpha2f
