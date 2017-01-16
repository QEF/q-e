!
! Copyright (C) 2002-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

!=----------------------------------------------------------------------------=!
        SUBROUTINE ecutoffs_setup( ecutwfc_, ecutrho_, ecfixed_, qcutz_, &
                                   q2sigma_, refg_ )
!------------------------------------------------------------------------------! 
          USE kinds,           ONLY: DP
          USE constants,       ONLY: eps8
          USE gvecw,           ONLY: ecutwfc
          USE gvecw,           ONLY: ecfixed, qcutz, q2sigma
          USE gvect,           ONLY: ecutrho
          USE gvecs,           ONLY: ecuts, dual, doublegrid
          USE pseudopotential, only: tpstab
          USE io_global,       only: stdout, ionode
          USE uspp,            only: okvan
          use betax,           only: mmx, refg

          IMPLICIT NONE
          REAL(DP), INTENT(IN) ::  ecutwfc_, ecutrho_, ecfixed_, qcutz_, &
                                   q2sigma_, refg_

          ecutwfc = ecutwfc_

          IF ( ecutrho_ <= 0.D0 ) THEN
             !
             dual = 4.D0
             !
          ELSE
             !
             dual = ecutrho_ / ecutwfc
             !
             IF ( dual <= 1.D0 ) &
                CALL errore( ' ecutoffs_setup ', ' invalid dual? ', 1 )
             !
          END IF

          doublegrid = ( dual > 4.D0 )
          IF ( doublegrid .AND. .NOT. okvan ) &
             CALL errore( 'setup', 'No USPP: set ecutrho=4*ecutwfc', 1 )
          ecutrho = dual * ecutwfc
          !
          IF ( doublegrid ) THEN
             !
             ecuts = 4.D0 * ecutwfc
             !
          ELSE
             !
             ecuts = ecutrho
             !
          END IF
          !
          ecfixed = ecfixed_
          qcutz   = qcutz_
          q2sigma = q2sigma_

          IF( refg_ < 0.0001d0 ) THEN
             tpstab = .FALSE.
             refg   = 0.05d0
          ELSE
             refg   = refg_
          END IF

          CALL set_interpolation_table_size( mmx, refg, ecutrho )

          RETURN
        END SUBROUTINE ecutoffs_setup


        SUBROUTINE set_interpolation_table_size( mmx, refg, gmax )
          USE control_flags,   only: thdyn
          USE kinds,           only: DP
          IMPLICIT NONE
          INTEGER, INTENT(OUT) :: mmx
          REAL(DP), INTENT(IN) :: refg
          REAL(DP), INTENT(IN) :: gmax
          IF( thdyn ) THEN
             !  ... a larger table is used when cell is moving to allow 
             !  ... large volume fluctuation
             mmx  = NINT( 2.0d0 * gmax / refg )
          ELSE
             mmx  = NINT( 1.2d0 * gmax / refg )
          END IF
          RETURN
        END SUBROUTINE set_interpolation_table_size


        SUBROUTINE gcutoffs_setup( alat, tk_inp, nk_inp, kpoints_inp )

!  (describe briefly what this routine does...)
!  ----------------------------------------------

          USE kinds, ONLY: DP
          USE gvecw, ONLY: ecutwfc,  gcutw
          USE gvect, ONLY: ecutrho,  gcutm
          USE gvecs, ONLY: ecuts, gcutms
          USE gvecw, ONLY: ekcut, gkcut
          USE constants, ONLY: eps8, pi

          IMPLICIT NONE

! ...     declare subroutine arguments
          REAL(DP), INTENT(IN) :: alat
          LOGICAL, INTENT(IN) :: tk_inp
          INTEGER, INTENT(IN) :: nk_inp
          REAL(DP), INTENT(IN) :: kpoints_inp(3,*)

! ...     declare other variables
          INTEGER   :: i
          REAL(DP) :: kcut, ksq
          REAL(DP) :: tpiba

!  end of declarations
!  ----------------------------------------------

! ...   Set Values for the cutoff


          IF( alat < eps8 ) THEN
            CALL errore(' cut-off setup ', ' alat too small ', 0)
          END IF

          tpiba = 2.0d0 * pi / alat 

          ! ...  Constant cutoff simulation parameters

          gcutw = ecutwfc / tpiba**2  ! wave function cut-off
          gcutm = ecutrho / tpiba**2  ! potential cut-off
          gcutms= ecuts   / tpiba**2  ! smooth mesh cut-off

          kcut = 0.0_DP
          IF ( tk_inp ) THEN
! ...       augment plane wave cutoff to include all k+G's
            DO i = 1, nk_inp
! ...         calculate modulus
              ksq = kpoints_inp( 1, i ) ** 2 + kpoints_inp( 2, i ) ** 2 + kpoints_inp( 3, i ) ** 2
              IF ( ksq > kcut ) kcut = ksq
            END DO
          END IF

          gkcut = ( sqrt( kcut ) + sqrt( gcutw ) ) ** 2

          ekcut = gkcut * tpiba ** 2

          RETURN
        END SUBROUTINE gcutoffs_setup

!  ----------------------------------------------

      SUBROUTINE cutoffs_print_info()

        !  Print out information about different cut-offs

        USE gvecw, ONLY: ecutwfc,  gcutw
        USE gvect, ONLY: ecutrho,  gcutm
        USE gvecw, ONLY: ecfixed, qcutz, q2sigma
        USE gvecw, ONLY: ekcut, gkcut
        USE gvecs, ONLY: ecuts, gcutms
        use betax, only: mmx, refg
        USE io_global, ONLY: stdout
        USE input_parameters, ONLY: ref_cell, ref_alat

        WRITE( stdout, 100 ) ecutwfc, ecutrho, ecuts, sqrt(gcutw), &
                             sqrt(gcutm), sqrt(gcutms)

        IF(ref_cell) WRITE( stdout,'(3X,"Reference Cell alat is",F14.8,1X,"A.U is used to Compute Gcutoffs:")') ref_alat ! BS : debug

        IF( qcutz > 0.0d0 ) THEN
          WRITE( stdout, 150 ) qcutz, q2sigma, ecfixed
        END IF

        WRITE( stdout,200) refg, mmx

100     FORMAT(/,3X,'Energy Cut-offs',/ &
                ,3X,'---------------',/ &
                ,3X,'Ecutwfc = ',F6.1,' Ry,   ', 3X,'Ecutrho = ',F6.1,' Ry,   ', 3X,'Ecuts = ',F6.1,' Ry',/ &
                ,3X,'Gcutwfc = ',F6.1,'     , ', 3X,'Gcutrho = ',F6.1,'       ', 3X,'Gcuts = ',F6.1)
150     FORMAT(  3X,'modified kinetic energy functional, with parameters:',/,   &
                 3X,'ecutz = ',f8.4,'  ecsig = ', f7.4,'  ecfix = ',f6.2)
200     FORMAT(  3X,'NOTA BENE: refg, mmx = ', f10.6,I6 )

        RETURN
      END SUBROUTINE cutoffs_print_info

!  ----------------------------------------------

      SUBROUTINE orthogonalize_info( )
        USE control_flags, ONLY: ortho_eps, ortho_max
        USE io_global, ONLY: stdout
        IMPLICIT NONE
           WRITE(stdout, 585)
           WRITE(stdout, 511) ortho_eps, ortho_max
  511   FORMAT(   3X,'Orthog. with lagrange multipliers : eps = ',E10.2, ',  max = ',I3)
  585   FORMAT(   3X,'Eigenvalues calculated without the kinetic term contribution')
        RETURN
      END SUBROUTINE orthogonalize_info


!  ----------------------------------------------


      SUBROUTINE electrons_print_info( )

          USE kinds, ONLY: DP
          USE electrons_base, ONLY: nbnd, nspin, nel, nelt, nupdwn, iupdwn, &
                                    f, qbac
          USE io_global, ONLY: stdout
          USE ions_base, ONLY: zv, nsp, na

          IMPLICIT NONE
          INTEGER :: i,is

          IF( nspin == 1) THEN
            WRITE(stdout,6) nelt, nbnd
            WRITE(stdout,7) ( f( i ), i = 1, nbnd )
          ELSE
            WRITE(stdout,8) nelt
            WRITE(stdout,9) nel(1)
            WRITE(stdout,7) ( f( i ), i = 1, nupdwn(1))
            WRITE(stdout,10) nel(2)
            WRITE(stdout,7) ( f( i ), i = iupdwn(2), ( iupdwn(2) + nupdwn(2) - 1 ) )
          END IF

         qbac=0.
         do is=1,nsp
           qbac=qbac+na(is)*zv(is)
         end do
         qbac=qbac-nelt
         if(qbac.ne.0) write(stdout,11) qbac


6         FORMAT(/,3X,'Electronic states',/  &
                  ,3X,'-----------------',/  &
                  ,3X,'Number of Electrons= ',I5,', of States = ',I5,/ &
                  ,3X,'Occupation numbers :')
7         FORMAT(2X,10F5.2)
8         FORMAT(/,3X,'Electronic states',/  &
                  ,3X,'-----------------',/  &
                  ,3X,'Local Spin Density calculation',/ &
                  ,3X,'Number of Electrons= ',I5)
9         FORMAT(  3X,'Spins up   = ', I5, ', occupations: ')
10        FORMAT(  3X,'Spins down = ', I5, ', occupations: ')
11        FORMAT(/,3X,'WARNING: system charge = ',F12.6)
          RETURN
      END SUBROUTINE electrons_print_info


!  ----------------------------------------------


      SUBROUTINE exch_corr_print_info()

        USE funct, ONLY: write_dft_name
        USE io_global, ONLY: stdout

        IMPLICIT NONE

        WRITE(stdout,800)
        call write_dft_name ( )
800 FORMAT(//,3X,'Exchange and correlations functionals',/ &
             ,3X,'-------------------------------------')

        RETURN
      END SUBROUTINE exch_corr_print_info



!  ----------------------------------------------



       SUBROUTINE ions_print_info( )
            
         !  Print info about input parameter for ion dynamic

         USE io_global,     ONLY: ionode, stdout
         USE control_flags, ONLY: tranp, amprp, tnosep, tolp, tfor, tsdp, &
                                  tzerop, tv0rd, taurdr, nbeg, tcp, tcap
         USE ions_base,     ONLY: tau_srt, if_pos, ind_srt, nsp, na, &
                                  amass, nat, fricp, greasp, rcmax
         USE ions_nose,     ONLY: tempw, ndega
         USE constants,     ONLY: amu_au

         IMPLICIT NONE
              
         integer is, ia, k, ic, isa
         LOGICAL :: ismb( 3 ) 
                
         WRITE( stdout, 50 ) 

         IF( .NOT. tfor ) THEN
           WRITE( stdout, 518 )
         ELSE
           WRITE( stdout, 520 )
           IF( tsdp ) THEN
             WRITE( stdout, 521 )
           ELSE
             WRITE( stdout, 522 )
           END IF
           WRITE( stdout, 523 ) ndega
           WRITE( stdout, 524 ) fricp, greasp
           IF( tv0rd ) THEN
              WRITE( stdout, 850 ) 
           ELSE IF ( tzerop ) THEN
               WRITE( stdout, 635 )
           ENDIF
         END IF 
              
         DO is = 1, nsp
           IF( tranp(is) ) THEN
             WRITE( stdout,510)
             WRITE( stdout,512) is, amprp(is)
           END IF
         END DO

         WRITE(stdout,660) 
         isa = 0
         DO IS = 1, nsp
           WRITE(stdout,1000) is, na(is), amass(is)*amu_au, amass(is), rcmax(is)
           DO IA = 1, na(is)
             isa = isa + 1
             WRITE(stdout,1010) ( tau_srt(k,isa), K = 1,3 )
           END DO
         END DO    

         IF ( ( nbeg > -1 ) .AND. ( .NOT. taurdr ) ) THEN
            WRITE(stdout,661)
         ELSE
            WRITE(stdout,662)
         ENDIF

         IF( tfor ) THEN

            IF( ANY( ( if_pos( 1:3, 1:nat ) == 0 )  ) ) THEN

              WRITE(stdout,1020)
              WRITE(stdout,1022)

              DO isa = 1, nat
                ia = ind_srt( isa )
                ismb( 1 ) = ( if_pos(1,ia) /= 0 )
                ismb( 2 ) = ( if_pos(2,ia) /= 0 )
                ismb( 3 ) = ( if_pos(3,ia) /= 0 )
                IF( .NOT. ALL( ismb ) ) THEN
                  WRITE( stdout, 1023 ) isa, ( ismb(k), K = 1, 3 )
                END IF
              END DO

            ELSE

              WRITE(stdout,1021)

            END IF
         END IF

         IF( tfor ) THEN
           if( ( tcp .or. tcap .or. tnosep ) .and. tsdp ) then
             call errore(' ions_print_info', &
               ' Temperature control not allowed with steepest descent',1)
           endif
           IF(.not. tcp .and. .not. tcap .and. .not. tnosep ) THEN
              WRITE( stdout,550)
           ELSE IF( tcp .and. tcap ) then
             call errore(' ions_print_info', ' Velocity rescaling not' &
                         //' compatible with random velocity initialization',1)
           ELSE IF( tcp .and. tnosep ) then
             call errore(' ions_print_info', ' Velocity rescaling and' &
                         //' Nose thermostat are incompatible',1)
           ELSE IF(tcap .and. tnosep ) then
             call errore(' ions_print_info', ' Nose thermostat not' &
                         //' compatible with random velocity initialization',1)
           ELSE IF(tcp) THEN
             WRITE( stdout,555) tempw,tolp
           ELSE IF(tcap) THEN
             WRITE( stdout,560) tempw,tolp
           ELSE IF(tnosep) THEN
             WRITE( stdout,595)
           ELSE
             WRITE( stdout,550)
           END IF
         END IF

   50 FORMAT(//,3X,'Ions Simulation Parameters',/ &
               ,3X,'--------------------------')

  510 FORMAT(   3X,'Initial random displacement of ionic coordinates',/, & 
                3X,' specie  amplitude')
  512 FORMAT(   3X,I7,2X,F9.6)

  518 FORMAT(   3X,'Ions are not allowed to move')
  520 FORMAT(   3X,'Ions are allowed to move')
  521 FORMAT(   3X,'Ions dynamics with steepest descent')
  522 FORMAT(   3X,'Ions dynamics with newton equations')
  523 format(   3X,'the temperature is computed for ',i5,' degrees of freedom')
  524 format(   3X,'ion dynamics with fricp = ',f7.4,' and greasp = ',f7.4)
  550 FORMAT(   3X,'Ionic temperature is not controlled')
  555 FORMAT(   3X,'Ionic temperature control via ', &
                   'rescaling of velocities :',/ &
               ,3X,'temperature required = ',F10.5,'K, ', &
                   'tolerance = ',F10.5,'K')
  560 FORMAT(   3X,'Ionic temperature control via ', &
                   'canonical velocities rescaling :',/ &
               ,3X,'temperature required = ',F10.5,'K, ', &
                   'tolerance = ',F10.5,'K')
  595 FORMAT(   3X,'Ionic temperature control via nose thermostat')
  635 FORMAT(   3X,'Zero initial momentum for ions')

  660 FORMAT(   3X,'Ionic position (from input)', /, &
                3X,'sorted by specie, and converted to real a.u. coordinates')
  661 FORMAT(   3X,'Ionic position will be re-read from restart file')
  662 FORMAT(   3X,'Ionic position read from input file')

  850 FORMAT(   3X,'Initial ion velocities read from input')

 1000 FORMAT(3X,'Species ',I3,' atoms = ',I4,' mass = ',F12.2, ' (a.u.), ', &
               & F12.2, ' (amu)', ' rcmax = ', F6.2, ' (a.u.)' )
 1010 FORMAT(3X,3(1X,F12.6))
 1020 FORMAT(/,3X,'NOT all atoms are allowed to move ')
 1021 FORMAT(/,3X,'All atoms are allowed to move')
 1022 FORMAT(  3X,' indx  ..x.. ..y.. ..z..')
 1023 FORMAT(  3X,I4,3(1X,L5))



         RETURN
       END SUBROUTINE ions_print_info


!  ----------------------------------------------

        subroutine cell_print_info( )

          USE constants, ONLY: au_gpa
          USE control_flags, ONLY: thdyn, tsdc, tzeroc, tbeg, nbeg, tpre
          USE control_flags, ONLY: tnoseh
          USE io_global, ONLY: stdout
          USE cell_base, ONLY: press, frich, greash, wmass

          IMPLICIT NONE

          WRITE(stdout,545 )
          IF ( tpre ) WRITE( stdout, 600 )
          IF ( tbeg ) THEN
            WRITE(stdout,546)
          ELSE
            WRITE(stdout,547)
            IF( nbeg > -1 ) WRITE( stdout, 548 )
          END IF

          IF( .NOT. thdyn ) THEN
            WRITE( stdout,525)
            WRITE( stdout,606)
          ELSE
            IF( tsdc ) THEN
              WRITE( stdout,526)
            ELSE
              IF( frich /= 0.0d0 ) THEN
                WRITE( stdout,602) frich, greash
              ELSE
                WRITE( stdout,527)
              END IF
              IF( tnoseh ) then
                WRITE( stdout,604) 
              ELSE
                WRITE( stdout,565)
              END IF
              IF( tzeroc ) THEN
                WRITE( stdout,563)
              ENDIF
            END IF
            WRITE( stdout,530) press * au_gpa, wmass
          END IF


 545     FORMAT(//,3X,'Cell Dynamics Parameters (from STDIN)',/ &
                  ,3X,'-------------------------------------')
 546     FORMAT(   3X,'Simulation cell read from STDIN')
 547     FORMAT(   3X,'Starting cell generated from CELLDM')
 548     FORMAT(   3X,'Cell parameters will be re-read from restart file')
 525     FORMAT(   3X,'Constant VOLUME Molecular dynamics')
 606     format(   3X,'cell parameters are not allowed to move')
 526     FORMAT(   3X,'Volume dynamics with steepest descent')
 527     FORMAT(   3X,'Volume dynamics with newton equations')
 530     FORMAT(   3X,'Constant PRESSURE Molecular dynamics:',/ &
                  ,3X,'External pressure (GPa) = ',F11.2,/ &
                  ,3X,'Volume mass             = ',F11.2)
 563     FORMAT(   3X,'Zero initial momentum for cell variables')
 565     FORMAT(   3X,'Volume dynamics: the temperature is not controlled')
 604     format(   3X,'cell parameters dynamics with nose` temp. control' )

 600  format( 3X, 'internal stress tensor calculated')
 602  format( 3X, 'cell parameters dynamics with frich = ',f7.4,            &
     &        3X, 'and greash = ',f7.4 )

        return
      end subroutine cell_print_info


!----------------------------------------------
SUBROUTINE gmeshinfo( )
!----------------------------------------------
   !
   !   Print out the number of g vectors for the different mesh
   !
   USE kinds,     ONLY: DP
   USE mp_global, ONLY: nproc_bgrp, intra_bgrp_comm
   USE io_global, ONLY: ionode, ionode_id, stdout
   USE mp,        ONLY: mp_max, mp_gather
   use smallbox_gvec,     only: ngb
   USE gvecw,     only: ngw_g, ngw, ngwx
   USE gvecs,     only: ngms_g, ngms, ngsx
   USE gvect,     only: ngm, ngm_g, ngmx

   IMPLICIT NONE

   INTEGER :: ip, ng_snd(3), ng_rcv( 3, nproc_bgrp )
   INTEGER :: ierr, min_val, max_val, i
   REAL(DP) :: avg_val

   IF(ionode) THEN
      WRITE( stdout,*)
      WRITE( stdout,*) '  Reciprocal Space Mesh'
      WRITE( stdout,*) '  ---------------------'
   END IF

   ng_snd(1) = ngm_g
   ng_snd(2) = ngm
   ng_snd(3) = ngmx
   CALL mp_gather(ng_snd, ng_rcv, ionode_id, intra_bgrp_comm)
   !
   IF(ionode) THEN
      min_val = MINVAL( ng_rcv(2,:) )
      max_val = MAXVAL( ng_rcv(2,:) )
      avg_val = REAL(SUM( ng_rcv(2,:) ))/nproc_bgrp
      WRITE( stdout,1000)
      WRITE( stdout,1011) ng_snd(1), min_val, max_val, avg_val
   END IF
   !
   ng_snd(1) = ngms_g
   ng_snd(2) = ngms
   ng_snd(3) = ngsx
   CALL mp_gather(ng_snd, ng_rcv, ionode_id, intra_bgrp_comm)
   !
   ierr = 0
   !
   IF(ionode) THEN
      WRITE( stdout,1001)
      min_val = MINVAL( ng_rcv(2,:) )
      max_val = MAXVAL( ng_rcv(2,:) )
      avg_val = REAL(SUM( ng_rcv(2,:) ))/nproc_bgrp
      WRITE( stdout,1011) ng_snd(1), min_val, max_val, avg_val
      IF( min_val < 1 ) ierr = ip
   END IF
   !
   CALL mp_max( ierr, intra_bgrp_comm )
   !
   IF( ierr > 0 ) &
      CALL errore( " gmeshinfo ", " Wow! some processors have no G-vectors ", ierr )
   !
   ng_snd(1) = ngw_g
   ng_snd(2) = ngw
   ng_snd(3) = ngwx
   CALL mp_gather(ng_snd, ng_rcv, ionode_id, intra_bgrp_comm)
   !
   IF(ionode) THEN
      WRITE( stdout,1002)
      min_val = MINVAL( ng_rcv(2,:) )
      max_val = MAXVAL( ng_rcv(2,:) )
      avg_val = REAL(SUM( ng_rcv(2,:) ))/nproc_bgrp
      WRITE( stdout,1011) ng_snd(1), min_val, max_val, avg_val
      IF( min_val < 1 ) ierr = ip
   END IF
   !
   CALL mp_max( ierr, intra_bgrp_comm )
   !
   IF( ierr > 0 ) &
      CALL errore( " gmeshinfo ", " Wow! some processors have no G-vectors ", ierr )
   !
   IF(ionode .AND. ngb > 0 ) THEN
      WRITE( stdout,1050)
      WRITE( stdout,1060) ngb
   END IF

   1000    FORMAT(3X,'Large Mesh',/, &
           '     Global(ngm_g)    MinLocal       MaxLocal      Average') 
   1001    FORMAT(3X,'Smooth Mesh',/, &
           '     Global(ngms_g)   MinLocal       MaxLocal      Average') 
   1002    FORMAT(3X,'Wave function Mesh',/, &
           '     Global(ngw_g)    MinLocal       MaxLocal      Average') 
   1011    FORMAT(  3I15, F15.2 )
   1050    FORMAT(/,3X,'Small box Mesh')
   1060    FORMAT( 3X, 'ngb = ', I12, ' not distributed to processors' )

   RETURN

END SUBROUTINE gmeshinfo

!----------------------------------------------
SUBROUTINE constraint_info()
!----------------------------------------------
   USE kinds,              ONLY: DP
   USE constraints_module, ONLY: nconstr, constr_tol, &
                                 constr_type, constr, constr_target
   USE io_global,          ONLY: ionode, stdout
   USE control_flags,      ONLY: lconstrain
   !
   IMPLICIT NONE
   !
   INTEGER :: ic
   !
   IF( lconstrain .AND. ionode ) THEN
      !
      WRITE( stdout, 10 ) 
      WRITE( stdout, 20 ) nconstr, constr_tol
      !
      DO ic = 1, nconstr
         !
         IF( constr_type( ic ) == 3 ) THEN
            !
            ! distance
            !
            WRITE( stdout, 30 ) ic
            WRITE( stdout, 40 ) NINT( constr(1,ic) ), &
                                NINT( constr(2,ic) ), constr_target(ic)
            !
         END IF
         !
      END DO
      !
   END IF
   !
10 FORMAT( 3X, "Using constrained dynamics")
20 FORMAT( 3X, "number of constrain and tolerance: ", I5, D10.2)
30 FORMAT( 3X, "constrain ", I5, " type distance ")
40 FORMAT( 3X, "  atoms ", I5, I5, " target dist ", F10.5)
   !
END SUBROUTINE constraint_info


SUBROUTINE new_atomind_constraints()
   !
   USE kinds,              ONLY: DP
   USE constraints_module, ONLY: constr
   USE ions_base,          ONLY: ind_bck
   !
   IMPLICIT NONE
   !
   INTEGER  :: ic, ia
   INTEGER  :: iaa
   REAL(DP) :: aa
   !
   !  Substitute the atom index given in the input file
   !  with the new atom index, after the sort in the
   !  atomic coordinates.
   !
   DO ic = 1, SIZE( constr, 2 )
      DO ia = 1, SIZE( constr, 1 )
         IF( constr( ia, ic ) > 0.0d0 ) THEN
            iaa = NINT( constr( ia, ic ) )
            aa  = DBLE( ind_bck( iaa ) )
            constr( ia, ic ) = aa
         END IF
      END DO
   END DO
   !
   RETURN
   !
END SUBROUTINE new_atomind_constraints


SUBROUTINE compute_stress_x( stress, detot, h, omega )
   USE kinds, ONLY : DP
   IMPLICIT NONE
   REAL(DP), INTENT(OUT) :: stress(3,3)
   REAL(DP), INTENT(IN)  :: detot(3,3), h(3,3), omega
   integer :: i, j
   do i=1,3 
      do j=1,3
         stress(i,j)=-1.d0/omega*(detot(i,1)*h(j,1)+              &
     &                      detot(i,2)*h(j,2)+detot(i,3)*h(j,3))
      enddo
   enddo
   return
END SUBROUTINE compute_stress_x
!-----------------------------------------------------------------------
subroutine formf( tfirst, eself )
  !-----------------------------------------------------------------------

  !computes (a) the self-energy eself of the ionic pseudocharges;
  !         (b) the form factors of: (i) pseudopotential (vps),
  !             (ii) ionic pseudocharge (rhops)
  !         also calculated the derivative of vps with respect to
  !         g^2 (dvps)
  ! 
  USE kinds,           ONLY : DP
  use mp,              ONLY : mp_sum
  use control_flags,   ONLY : iprint, tpre, iverbosity
  use io_global,       ONLY : stdout
  use mp_global,       ONLY : intra_bgrp_comm
  use gvecs,           ONLY : ngms
  use cell_base,       ONLY : omega, tpiba2, tpiba
  use ions_base,       ONLY : rcmax, zv, nsp, na
  use local_pseudo,    ONLY : vps, vps0, rhops, dvps, drhops
  use atom,            ONLY : rgrid
  use uspp_param,      ONLY : upf, oldvan
  use pseudo_base,     ONLY : compute_rhops, formfn, formfa, compute_eself
  use pseudopotential, ONLY : tpstab, vps_sp, dvps_sp
  use splines,         ONLY : spline
  use gvect,           ONLY : gstart, gg
  use constants,       ONLY : autoev
  !
  implicit none
  logical      :: tfirst
  real(DP)    :: eself, DeltaV0
  !
  real(DP)    :: vpsum, rhopsum
  integer      :: is, ig
  REAL(DP)    :: cost1, xg

  call start_clock( 'formf' )
  !
  IF( .NOT. ALLOCATED( rgrid ) ) &
     CALL errore( ' formf ', ' rgrid not allocated ', 1 )
  IF( .NOT. ALLOCATED( upf ) ) &
     CALL errore( ' formf ', ' upf not allocated ', 1 )
  !
  ! calculation of gaussian selfinteraction
  !
  eself = compute_eself( na, zv, rcmax, nsp )

  if( tfirst .or. ( iverbosity > 2 ) )then
     WRITE( stdout, 1200 ) eself
  endif
  !
  1200 format(/,3x,'formf: eself=',f12.5)
  !
  do is = 1, nsp

     IF( tpstab ) THEN
        !
        !  Use interpolation table, with cubic spline
        !
        cost1 = 1.0d0/omega
        !
        IF( gstart == 2 ) THEN
           vps (1,is) =  vps_sp(is)%y(1) * cost1
           dvps(1,is) = dvps_sp(is)%y(1) * cost1
        END IF
        !
        DO ig = gstart, ngms
           xg = SQRT( gg(ig) ) * tpiba
           vps (ig,is) = spline(  vps_sp(is), xg ) * cost1
           dvps(ig,is) = spline( dvps_sp(is), xg ) * cost1
        END DO
        !
     ELSE

        call formfn( rgrid(is)%r, rgrid(is)%rab, &
                     upf(is)%vloc(1:rgrid(is)%mesh), zv(is), rcmax(is), gg, &
                     omega, tpiba2, rgrid(is)%mesh, ngms, oldvan(is), tpre, &
                     vps(:,is), vps0(is), dvps(:,is) )

! obsolete BHS form
! call formfa( vps(:,is), dvps(:,is), rc1(is), rc2(is), wrc1(is), wrc2(is), &
!              rcl(:,is,lloc(is)), al(:,is,lloc(is)), bl(:,is,lloc(is)),    &
!              zv(is), rcmax(is), g, omega, tpiba2, ngms, gstart, tpre )

     END IF
     !
     !     fourier transform of local pp and gaussian nuclear charge
     !
     call compute_rhops( rhops(:,is), drhops(:,is), zv(is), rcmax(is), gg, &
                         omega, tpiba2, ngms, tpre )

     if( tfirst .or. ( iverbosity > 2 ) )then
        vpsum = SUM( vps( 1:ngms, is ) )
        rhopsum = SUM( rhops( 1:ngms, is ) )
        call mp_sum( vpsum, intra_bgrp_comm )
        call mp_sum( rhopsum, intra_bgrp_comm )
        WRITE( stdout,1250) vps(1,is),rhops(1,is)
        WRITE( stdout,1300) vpsum,rhopsum
     endif
     !
  end do
  ! 
  ! ... DeltaV0 is the shift to be applied to eigenvalues
  ! ... in order to align them to other plane wave codes
  !
  DeltaV0 = 0.0_dp
  DO is = 1, nsp
     !
     ! ...  na(is)/omega is the structure factor at G=0
     !
     DeltaV0 = DeltaV0 + na(is) / omega * vps0(is)
  END DO
  !
  IF ( tfirst .or. ( iverbosity > 2 ) ) THEN
      write(6,'("   Delta V(G=0): ",f10.6,"Ry, ",f11.6,"eV")') &
         deltaV0, deltaV0*autoev
  END IF
  !
  call stop_clock( 'formf' )
  !
  1250 format(3x,'formf:     vps(g=0)=',f12.7,'     rhops(g=0)=',f12.7)
  1300 format(3x,'formf: sum_g vps(g)=',f12.7,' sum_g rhops(g)=',f12.7)
  !
  return
end subroutine formf
!
!-----------------------------------------------------------------------
SUBROUTINE newnlinit()
  !-----------------------------------------------------------------------
  !
  ! ... this routine calculates arrays beta, qq, qgb, rhocb
  ! ... and derivatives w.r.t. cell parameters dbeta
  ! ... See also comments in nlinit
  !
  use control_flags,    ONLY : tpre
  use pseudopotential,  ONLY : tpstab
  use cp_interfaces,    ONLY : interpolate_beta, interpolate_qradb, compute_qradx, compute_betagx, &
                               exact_beta, check_tables, exact_qradb, build_pstab, build_cctab
  use betax,            only : mmx, refg
  use kinds,            only : dp
  use io_global,        only : ionode, stdout
  !
  IMPLICIT NONE
  !
  LOGICAL  :: recompute_table
  REAL(DP) :: gmax
  ! 
  ! ... initialization for vanderbilt species
  !
  CALL start_clock( 'newnlinit' )

  IF( tpstab ) THEN

     recompute_table = tpre .AND. check_tables( gmax )
     !
     IF ( recompute_table ) THEN

        IF( ionode ) &
           WRITE( stdout, * ) "newnliinit: recomputing the pseudopotentials tables" 
        !"!

        CALL set_interpolation_table_size( mmx, refg, gmax )

        CALL compute_qradx( tpre )

        call compute_betagx( tpre )

        call build_pstab()
        !
        call build_cctab()

     END IF
     !
     !     initialization that is common to all species
     !
     CALL interpolate_beta( tpre )
     !
     CALL interpolate_qradb( tpre )
     !
  ELSE
     !
     ! ... this is mainly for testing
     !
     CALL exact_beta( tpre )
     !
     CALL exact_qradb( tpre )
     !
  END IF
  !
  ! ... non-linear core-correction   ( rhocb(ig,is) )
  !
  CALL core_charge_ftr( tpre )

  CALL stop_clock( 'newnlinit' )
  !
  RETURN
  !
END SUBROUTINE newnlinit
!
!-----------------------------------------------------------------------
subroutine nlfh_x( stress, bec_bgrp, dbec, lambda, descla )
  !-----------------------------------------------------------------------
  !
  !     contribution to the internal stress tensor due to the constraints
  !
  USE kinds,             ONLY : DP
  use uspp,              ONLY : nkb, qq
  use uspp_param,        ONLY : nh, nhm, nvb, ish
  use ions_base,         ONLY : na
  use electrons_base,    ONLY : nbspx, nbsp, nudx, nspin, nupdwn, iupdwn, ibgrp_g2l
  use cell_base,         ONLY : omega, h
  use constants,         ONLY : pi, fpi, au_gpa
  use io_global,         ONLY : stdout
  use control_flags,     ONLY : iverbosity
  USE descriptors,       ONLY : la_descriptor
  USE mp,                ONLY : mp_sum
  USE mp_global,         ONLY : intra_bgrp_comm, inter_bgrp_comm

!
  implicit none

  TYPE(la_descriptor), INTENT(IN) :: descla(:)
  REAL(DP), INTENT(INOUT) :: stress(3,3) 
  REAL(DP), INTENT(IN)    :: bec_bgrp( :, : ), dbec( :, :, :, : )
  REAL(DP), INTENT(IN)    :: lambda( :, :, : )
!
  INTEGER  :: i, j, ii, jj, inl, iv, jv, ia, is, iss, nss, istart
  INTEGER  :: jnl, ir, ic, nr, nc, ibgrp_i, nrcx
  REAL(DP) :: fpre(3,3), TT, T1, T2
  !
  REAL(DP), ALLOCATABLE :: tmpbec(:,:), tmpdh(:,:), temp(:,:), bec(:,:,:)
  !
  nrcx = MAXVAL( descla( : )%nrcx )
  !
  ALLOCATE( bec( nkb, nrcx, nspin ) )
  !
  DO iss = 1, nspin
     IF( descla( iss )%active_node > 0 ) THEN
        nss = nupdwn( iss )
        istart = iupdwn( iss )
        ic = descla( iss )%ic
        nc = descla( iss )%nc
        DO i=1,nc
           ibgrp_i = ibgrp_g2l( i+istart-1+ic-1 )
           IF( ibgrp_i > 0 ) THEN
              bec( :, i, iss ) = bec_bgrp( :, ibgrp_i )
           ELSE
              bec( :, i, iss ) = 0.0d0
           END IF
        END DO
     ELSE
        bec(:,:,iss)   = 0.0d0
     END IF
  END DO

  CALL mp_sum( bec, inter_bgrp_comm )
  !
  IF (nspin == 1) THEN
     IF( ( descla( 1 )%active_node > 0 ) ) THEN
        ALLOCATE ( tmpbec(nhm,nrcx), tmpdh(nrcx,nhm), temp(nrcx,nrcx) )
     ENDIF
  ELSEIF (nspin == 2) THEN
     IF( ( descla( 1 )%active_node > 0 ) .OR. ( descla( 2 )%active_node > 0 ) ) THEN
        ALLOCATE ( tmpbec(nhm,nrcx), tmpdh(nrcx,nhm), temp(nrcx,nrcx) )
     END IF
  ENDIF
  !
  fpre = 0.d0
  !
  do ii=1,3

     do jj=1,3

        do is=1,nvb

           do ia=1,na(is)

              do iss = 1, nspin
                 !
                 istart = iupdwn( iss )
                 nss    = nupdwn( iss )
                 !
                 IF( descla( iss )%active_node > 0 ) THEN

                    nr = descla( iss )%nr
                    nc = descla( iss )%nc
                    ir = descla( iss )%ir
                    ic = descla( iss )%ic

                    tmpbec = 0.d0
                    tmpdh  = 0.d0
!
                    do iv=1,nh(is)
                       do jv=1,nh(is)
                          inl=ish(is)+(jv-1)*na(is)+ia
                          if(abs(qq(iv,jv,is)).gt.1.e-5) then
                             do i = 1, nc
                                tmpbec(iv,i) = tmpbec(iv,i) +  qq(iv,jv,is) * bec( inl, i, iss  )
                             end do
                          endif
                       end do
                    end do

                    do iv=1,nh(is)
                       inl=ish(is)+(iv-1)*na(is)+ia
                       do i = 1, nr
                          tmpdh(i,iv) = dbec( inl, i + (iss-1)*nrcx, ii, jj )
                       end do
                    end do

                    if(nh(is).gt.0)then

                       CALL dgemm &
                       ( 'N', 'N', nr, nc, nh(is), 1.0d0, tmpdh, nrcx, tmpbec, nhm, 0.0d0, temp, nrcx )

                       do j = 1, nc
                          do i = 1, nr
                             fpre(ii,jj) = fpre(ii,jj) + 2D0 * temp( i, j ) * lambda(i,j,iss)
                          end do
                       end do
                    endif

                 END IF
                 !
              end do
              !
           end do
           !
        end do
        !
     end do
     !
  end do

  CALL mp_sum( fpre, intra_bgrp_comm )

  do i=1,3
     do j=1,3
        stress(i,j)=stress(i,j)+ &
                    (fpre(i,1)*h(j,1)+fpre(i,2)*h(j,2)+fpre(i,3)*h(j,3))/omega
     enddo
  enddo

  IF (allocated(tmpbec)) THEN
     DEALLOCATE ( tmpbec, tmpdh, temp )
  END IF

  DEALLOCATE( bec )


  IF( iverbosity > 1 ) THEN
     WRITE( stdout,*) 
     WRITE( stdout,*) "constraints contribution to stress"
     WRITE( stdout,5555) ((-fpre(i,j),j=1,3),i=1,3)
     fpre = MATMUL( fpre, TRANSPOSE( h ) ) / omega * au_gpa * 10.0d0
     WRITE( stdout,5555) ((fpre(i,j),j=1,3),i=1,3)
     WRITE( stdout,*) 
  END IF
!

5555  FORMAT(1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5//)

  return
end subroutine nlfh_x


!-----------------------------------------------------------------------
subroutine nlinit
  !-----------------------------------------------------------------------
  !
  !     this routine allocates and initializes arrays beta, qq, qgb,
  !     rhocb, and derivatives w.r.t. cell parameters dbeta
  !
  !       beta(ig,l,is) = 4pi/sqrt(omega) y^r(l,q^)
  !                               int_0^inf dr r^2 j_l(qr) betar(l,is,r)
  !
  !       Note that beta(g)_lm,is = (-i)^l*beta(ig,l,is) (?)
  !
  !       qq_ij=int_0^r q_ij(r)=omega*qg(g=0)
  !
  !     beta and qradb are first calculated on a fixed linear grid in |G|
  !     (betax, qradx) then calculated on the box grid by interpolation
  !     (this is done in routine newnlinit)
  !
      use kinds,           ONLY : dp
      use control_flags,   ONLY : iprint, tpre
      use io_global,       ONLY : stdout, ionode
      use gvecw,           ONLY : ngw
      use core,            ONLY : rhocb, allocate_core
      use constants,       ONLY : pi, fpi
      use ions_base,       ONLY : na, nsp
      use uspp,            ONLY : aainit, beta, qq, dvan, nhtol, nhtolm, indv,&
                                  dbeta
      use uspp_param,      ONLY : upf, lmaxq, nbetam, lmaxkb, nhm, nh, ish, nvb
      use atom,            ONLY : rgrid
      use qgb_mod,         ONLY : qgb, dqgb
      use smallbox_gvec,   ONLY : ngb
      use gvect,           ONLY : ngm
      use cp_interfaces,   ONLY : pseudopotential_indexes, compute_dvan, &
                                  compute_betagx, compute_qradx, build_pstab, build_cctab
      USE fft_base,        ONLY : dfftp
      use pseudopotential, ONLY : tpstab

!
      implicit none
!
      integer  is, il, l, ir, iv, jv, lm, ind, ltmp, i0
      real(dp), allocatable:: fint(:), jl(:),  jltmp(:), djl(:),    &
     &              dfint(:)
      real(dp) xg, xrg, fac

      CALL start_clock( 'nlinit' )

      IF( ionode ) THEN
        WRITE( stdout, 100 )
 100    FORMAT( //, &
                3X,'Pseudopotentials initialization',/, &
                3X,'-------------------------------' )
      END IF

      IF( .NOT. ALLOCATED( rgrid ) ) &
         CALL errore( ' nlinit ', ' rgrid not allocated ', 1 )
      IF( .NOT. ALLOCATED( upf ) ) &
         CALL errore( ' nlinit ', ' upf not allocated ', 1 )
      !
      !   initialize indexes
      !
      CALL pseudopotential_indexes( )
      !
      !   initialize array ap
      !
      call aainit( lmaxkb + 1 )
      !
      CALL allocate_core( dfftp%nnr, ngm, ngb, nsp )
      !
      !
      allocate( beta( ngw, nhm, nsp ) )
      allocate( qgb( ngb, nhm*(nhm+1)/2, nsp ) )
      allocate( qq( nhm, nhm, nsp ) )
      qq  (:,:,:) =0.d0
      IF (tpre) THEN
         allocate( dqgb( ngb, nhm*(nhm+1)/2, nsp, 3, 3 ) )
         allocate( dbeta( ngw, nhm, nsp, 3, 3 ) )
      END IF
      !
      !     initialization for vanderbilt species
      !
      CALL compute_qradx( tpre )
      !    
      !     initialization that is common to all species
      !   
      WRITE( stdout, fmt="(//,3X,'Common initialization' )" )

      do is = 1, nsp
         WRITE( stdout, fmt="(/,3X,'Specie: ',I5)" ) is
         !     fac converts ry to hartree
         fac=0.5d0
         do iv = 1, nh(is)
            WRITE( stdout,901) iv, indv(iv,is), nhtol(iv,is)
         end do
 901     format(2x,i2,'  indv= ',i2,'   ang. mom= ',i2)
         !
         WRITE( stdout,*)
         WRITE( stdout,'(20x,a)') '    dion '
         do iv = 1, upf(is)%nbeta
            WRITE( stdout,'(8f9.4)') ( fac*upf(is)%dion(iv,jv), jv = 1, upf(is)%nbeta )
         end do
         !
      end do
      !
      !   calculation of array  betagx(ig,iv,is)
      !
      call compute_betagx( tpre )
      !
      !   calculate array  dvan(iv,jv,is)
      !
      call compute_dvan()
      !
      IF( tpstab ) THEN

         call build_pstab()
         !
         call build_cctab()
         !
      END IF
      !
      ! newnlinit stores qgb and qq, calculates arrays  beta  rhocb
      ! and derivatives wrt cell dbeta
      !
      call newnlinit()

      CALL stop_clock( 'nlinit' )

      return
end subroutine nlinit

!-------------------------------------------------------------------------
subroutine qvan2b(ngy,iv,jv,is,ylm,qg,qradb)
  !--------------------------------------------------------------------------
  !
  !     q(g,l,k) = sum_lm (-i)^l ap(lm,l,k) yr_lm(g^) qrad(g,l,l,k)
  !
  USE kinds,         ONLY : DP
  use control_flags, ONLY : iprint, tpre
  use uspp,          ONLY : nlx, lpx, lpl, ap, indv, nhtolm
  use smallbox_gvec,         ONLY : ngb
  use uspp_param,    ONLY : lmaxq, nbetam
  use ions_base,     ONLY : nsp
! 
  implicit none
  !
  integer,     intent(in)  :: ngy, iv, jv, is
  real(DP),    intent(in)  :: ylm( ngb, lmaxq*lmaxq )
  real(DP),    intent(in)  :: qradb( ngb, nbetam*(nbetam+1)/2, lmaxq, nsp )
  complex(DP), intent(out) :: qg( ngb )
!
  integer      :: ivs, jvs, ijvs, ivl, jvl, i, ii, ij, l, lp, ig
  complex(DP) :: sig
  ! 
  !       iv  = 1..8     s_1 p_x1 p_z1 p_y1 s_2 p_x2 p_z2 p_y2
  !       ivs = 1..4     s_1 s_2 p_1 p_2
  !       ivl = 1..4     s p_x p_z p_y
  ! 
  ivs=indv(iv,is)
  jvs=indv(jv,is)
  if (ivs >= jvs) then
     ijvs = ivs*(ivs-1)/2 + jvs
  else
     ijvs = jvs*(jvs-1)/2 + ivs
  end if
  ! ijvs is the packed index for (ivs,jvs)
  ivl=nhtolm(iv,is)
  jvl=nhtolm(jv,is)
  if (ivl > nlx .OR. jvl > nlx) &
       call errore (' qvan2b ', ' wrong dimensions', MAX(ivl,jvl))
  !
  qg(:) = (0.d0, 0.d0)
  !
  !     lpx = max number of allowed y_lm
  !     lp  = composite lm to indentify them
  !
  do i=1,lpx(ivl,jvl)
     lp=lpl(ivl,jvl,i)
     if (lp > lmaxq*lmaxq) call errore(' qvan2b ',' lp out of bounds ',lp)
     !
     !     extraction of angular momentum l from lp:  
     !     l = int ( sqrt( DBLE(l-1) + epsilon) ) + 1
     !
     if (lp == 1) then
        l=1         
     else if ((lp >= 2) .and. (lp <= 4)) then
        l=2
     else if ((lp >= 5) .and. (lp <= 9)) then
        l=3
     else if ((lp >= 10).and.(lp <= 16)) then
        l=4
     else if ((lp >= 17).and.(lp <= 25)) then
        l=5
     else if ((lp >= 26).and.(lp <= 36)) then 
        l=6
     else if ((lp >= 37).and.(lp <= 49)) then 
        l=7
     else
        call errore(' qvan2b ',' not implemented ',lp)
     endif
     !     
     !       sig= (-i)^l
     !
     sig=(0.d0,-1.d0)**(l-1)
     sig=sig*ap(lp,ivl,jvl)
     do ig=1,ngy
        qg(ig)=qg(ig)+sig*ylm(ig,lp)*qradb(ig,ijvs,l,is)
     end do
  end do

  return
end subroutine qvan2b

!-------------------------------------------------------------------------
subroutine dqvan2b(ngy,iv,jv,is,ylm,dylm,dqg,dqrad,qradb)
  !--------------------------------------------------------------------------
  !
  !     dq(i,j) derivatives wrt to h(i,j) of q(g,l,k) calculated in qvan2b
  !
  USE kinds,         ONLY : DP
  use control_flags, ONLY : iprint, tpre
  use uspp,          ONLY : nlx, lpx, lpl, ap, indv, nhtolm
  use smallbox_gvec,         ONLY : ngb
  use uspp_param,    ONLY : lmaxq, nbetam
  use ions_base,     ONLY : nsp

  implicit none

  integer,     intent(in)  :: ngy, iv, jv, is
  REAL(DP),    INTENT(IN)  :: ylm( ngb, lmaxq*lmaxq ), dylm( ngb, lmaxq*lmaxq, 3, 3 )
  complex(DP), intent(out) :: dqg( ngb, 3, 3 )
  REAL(DP),    INTENT(IN)  :: dqrad( ngb, nbetam*(nbetam+1)/2, lmaxq, nsp, 3, 3 )
  real(DP),    intent(in)  :: qradb( ngb, nbetam*(nbetam+1)/2, lmaxq, nsp )

  integer      :: ivs, jvs, ijvs, ivl, jvl, i, ii, ij, l, lp, ig
  complex(DP) :: sig, z1, z2, zfac
  !
  ! 
  !       iv  = 1..8     s_1 p_x1 p_z1 p_y1 s_2 p_x2 p_z2 p_y2
  !       ivs = 1..4     s_1 s_2 p_1 p_2
  !       ivl = 1..4     s p_x p_z p_y
  ! 
  ivs=indv(iv,is)
  jvs=indv(jv,is)
  !
  if (ivs >= jvs) then
     ijvs = ivs*(ivs-1)/2 + jvs
  else
     ijvs = jvs*(jvs-1)/2 + ivs
  end if
  !
  ! ijvs is the packed index for (ivs,jvs)
  !
  ivl=nhtolm(iv,is)
  jvl=nhtolm(jv,is)
  !
  if (ivl > nlx .OR. jvl > nlx) &
       call errore (' qvan2 ', ' wrong dimensions (2)', MAX(ivl,jvl))
  !
  dqg(:,:,:) = (0.d0, 0.d0)

  !  lpx = max number of allowed y_lm
  !  lp  = composite lm to indentify them

  z1 = 0.0d0
  z2 = 0.0d0
  do i=1,lpx(ivl,jvl)
     lp=lpl(ivl,jvl,i)
     if (lp > lmaxq*lmaxq) call errore(' dqvan2b ',' lp out of bounds ',lp)

     !  extraction of angular momentum l from lp:  
     !  l = int ( sqrt( DBLE(l-1) + epsilon) ) + 1
     !
     if (lp == 1) then
        l=1         
     else if ((lp >= 2) .and. (lp <= 4)) then
        l=2
     else if ((lp >= 5) .and. (lp <= 9)) then
        l=3
     else if ((lp >= 10).and.(lp <= 16)) then
        l=4
     else if ((lp >= 17).and.(lp <= 25)) then
        l=5
     else if ((lp >= 26).and.(lp <= 36)) then 
        l=6
     else if ((lp >= 37).and.(lp <= 49)) then 
        l=7
     else
        call errore(' qvan2b ',' not implemented ',lp)
     endif
     !     
     !       sig= (-i)^l
     !
     sig = (0.0d0,-1.0d0)**(l-1)
     !
     sig = sig * ap( lp, ivl, jvl ) 
     !
     do ij=1,3
        do ii=1,3
           do ig=1,ngy
              zfac = ylm(ig,lp) * dqrad(ig,ijvs,l,is,ii,ij)
              zfac = zfac - dylm(ig,lp,ii,ij) * qradb(ig,ijvs,l,is)
              dqg(ig,ii,ij) = dqg(ig,ii,ij) +  sig * zfac
           end do
        end do
     end do
  end do
  !
  ! WRITE(6,*) 'DEBUG dqvan2b: ', z1, z2
  !
  return
end subroutine dqvan2b

!-----------------------------------------------------------------------
subroutine dylmr2_( nylm, ngy, g, gg, ainv, dylm )
  !-----------------------------------------------------------------------
  !
  ! temporary CP interface for PW routine dylmr2
  ! dylmr2  calculates d Y_{lm} /d G_ipol
  ! dylmr2_ calculates G_ipol \sum_k h^(-1)(jpol,k) (dY_{lm} /dG_k)
  !
  USE kinds, ONLY: DP

  implicit none
  !
  integer,   intent(IN)  :: nylm, ngy
  real(DP), intent(IN)  :: g (3, ngy), gg (ngy), ainv(3,3)
  real(DP), intent(OUT) :: dylm (ngy, nylm, 3, 3)
  !
  integer :: ipol, jpol, lm, ig
  real(DP), allocatable :: dylmaux (:,:,:)
  !
  allocate ( dylmaux(ngy,nylm,3) )
  !
  dylmaux(:,:,:) = 0.d0
  !
  do ipol =1,3
     call dylmr2 (nylm, ngy, g, gg, dylmaux(1,1,ipol), ipol)
  enddo
  !
  do ipol =1,3
     do jpol =1,3
        do lm=1,nylm
           do ig = 1, ngy
              dylm (ig,lm,ipol,jpol) = (dylmaux(ig,lm,1) * ainv(jpol,1) + & 
                                        dylmaux(ig,lm,2) * ainv(jpol,2) + &
                                        dylmaux(ig,lm,3) * ainv(jpol,3) ) &
                                       * g(ipol,ig)
           end do
        end do
     end do
  end do
  !
  deallocate ( dylmaux )
  !
  return
  !
end subroutine dylmr2_


SUBROUTINE print_lambda_x( lambda, descla, n, nshow, ccc, iunit )
    USE kinds, ONLY : DP
    USE descriptors,       ONLY: la_descriptor
    USE io_global,         ONLY: stdout, ionode
    USE cp_interfaces,     ONLY: collect_lambda
    USE electrons_base,    ONLY: nudx
    IMPLICIT NONE
    real(DP), intent(in) :: lambda(:,:,:), ccc
    TYPE(la_descriptor), INTENT(IN) :: descla(:)
    integer, intent(in) :: n, nshow
    integer, intent(in), optional :: iunit
    !
    integer :: nnn, j, un, i, is
    real(DP), allocatable :: lambda_repl(:,:)
    if( present( iunit ) ) then
      un = iunit
    else
      un = stdout
    end if
    nnn = min( nudx, nshow )
    ALLOCATE( lambda_repl( nudx, nudx ) )
    IF( ionode ) WRITE( un,*)
    DO is = 1, SIZE( lambda, 3 )
       CALL collect_lambda( lambda_repl, lambda(:,:,is), descla(is) )
       IF( ionode ) THEN
          WRITE( un,3370) '    lambda   nudx, spin = ', nudx, is
          IF( nnn < n ) WRITE( un,3370) '    print only first ', nnn
          DO i=1,nnn
             WRITE( un,3380) (lambda_repl(i,j)*ccc,j=1,nnn)
          END DO
       END IF
    END DO
    DEALLOCATE( lambda_repl )
3370   FORMAT(26x,a,2i4)
3380   FORMAT(9f8.4)
    RETURN
END SUBROUTINE print_lambda_x
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
   SUBROUTINE denlcc_x( nnr, nspin, vxcr, sfac, drhocg, dcc )
!-----------------------------------------------------------------------
!
! derivative of non linear core correction exchange energy wrt cell 
! parameters h 
! Output in dcc
!
      USE kinds,              ONLY: DP
      USE ions_base,          ONLY: nsp
      USE gvect, ONLY: gstart, g, gg
      USE gvecs,              ONLY: ngms
      USE gvect,              ONLY: ngm, nl
      USE cell_base,          ONLY: omega, ainv, tpiba2
      USE mp,                 ONLY: mp_sum
      USE mp_global,          ONLY: intra_bgrp_comm
      USE uspp_param,         ONLY: upf
      USE fft_interfaces,     ONLY: fwfft
      USE fft_base,           ONLY: dfftp

      IMPLICIT NONE

      ! input

      INTEGER,     INTENT(IN) :: nnr, nspin
      REAL(DP),    INTENT(IN) :: vxcr( :, : )
      COMPLEX(DP), INTENT(IN) :: sfac( :, : )
      REAL(DP),    INTENT(IN) :: drhocg( :, : )

      ! output

      REAL(DP), INTENT(OUT) ::  dcc( :, : )

      ! local

      INTEGER     :: i, j, ig, is
      COMPLEX(DP) :: srhoc
      REAL(DP)    :: vxcc
      !
      COMPLEX(DP), ALLOCATABLE :: vxc( : )
!
      dcc = 0.0d0
      !
      ALLOCATE( vxc( nnr ) )
      !
      vxc(:) = vxcr(:,1)
      !
      IF( nspin > 1 ) vxc(:) = vxc(:) + vxcr(:,2)
      !
      CALL fwfft( 'Dense', vxc, dfftp )
      !
      DO i=1,3
         DO j=1,3
            DO ig = gstart, ngms
               srhoc = 0.0d0
               DO is = 1, nsp
                 IF( upf(is)%nlcc ) srhoc = srhoc + sfac( ig, is ) * drhocg( ig, is )
               ENDDO
               vxcc = DBLE( CONJG( vxc( nl( ig ) ) ) * srhoc ) / SQRT( gg(ig) * tpiba2 )
               dcc(i,j) = dcc(i,j) + vxcc * &
     &                      2.d0 * tpiba2 * g(i,ig) *                  &
     &                    (g(1,ig)*ainv(j,1) +                         &
     &                     g(2,ig)*ainv(j,2) +                         &
     &                     g(3,ig)*ainv(j,3) )
            ENDDO
         ENDDO
      ENDDO

      DEALLOCATE( vxc )

      dcc = dcc * omega

      CALL mp_sum( dcc( 1:3, 1:3 ), intra_bgrp_comm )

      RETURN
   END SUBROUTINE denlcc_x



!-----------------------------------------------------------------------
      SUBROUTINE dotcsc_x( eigr, cp, ngw, n )
!-----------------------------------------------------------------------
!
      USE kinds,              ONLY: DP
      USE ions_base,          ONLY: na, nsp, nat
      USE io_global,          ONLY: stdout
      USE gvect, ONLY: gstart
      USE uspp,               ONLY: nkb, qq
      USE uspp_param,         ONLY: nh, ish, nvb
      USE mp,                 ONLY: mp_sum
      USE mp_global,          ONLY: intra_bgrp_comm, nbgrp, inter_bgrp_comm
      USE cp_interfaces,      ONLY: nlsm1
      USE electrons_base,     ONLY: ispin, ispin_bgrp, nbspx_bgrp, ibgrp_g2l, iupdwn, nupdwn, nbspx
!
      IMPLICIT NONE
!
      INTEGER,     INTENT(IN) :: ngw, n
      COMPLEX(DP), INTENT(IN) :: eigr(:,:), cp(:,:)
! local variables
      REAL(DP) rsum, csc(n) ! automatic array
      COMPLEX(DP) temp(ngw) ! automatic array
 
      REAL(DP), ALLOCATABLE::  becp(:,:), cp_tmp(:), becp_tmp(:)
      INTEGER i,kmax,nnn,k,ig,is,ia,iv,jv,inl,jnl
      INTEGER :: ibgrp_i, ibgrp_k
!
      ALLOCATE( becp( nkb, nbspx_bgrp ) )
      ALLOCATE( cp_tmp( SIZE( cp, 1 ) ) )
      ALLOCATE( becp_tmp( nkb ) )
!
!     < beta | phi > is real. only the i lowest:
!

      CALL nlsm1( nbspx_bgrp, 1, nvb, eigr, cp, becp )

      nnn = MIN( 12, n )

      DO i = nnn, 1, -1

         csc = 0.0d0

         ibgrp_i = ibgrp_g2l( i )
         IF( ibgrp_i > 0 ) THEN
            cp_tmp = cp( :, ibgrp_i )
         ELSE 
            cp_tmp = 0.0d0
         END IF

         CALL mp_sum( cp_tmp, inter_bgrp_comm )

         kmax = i
!
         DO k=1,kmax
            ibgrp_k = ibgrp_g2l( k )
            IF( ibgrp_k > 0 ) THEN
               DO ig=1,ngw
                  temp(ig)=CONJG(cp(ig,ibgrp_k))*cp_tmp(ig)
               END DO
               csc(k)=2.d0*DBLE(SUM(temp))
               IF (gstart == 2) csc(k)=csc(k)-DBLE(temp(1))
            END IF
         END DO

         CALL mp_sum( csc( 1:kmax ), intra_bgrp_comm )

         IF( ibgrp_i > 0 ) THEN
            becp_tmp = becp( :, ibgrp_i )
         ELSE 
            becp_tmp = 0.0d0
         END IF

         CALL mp_sum( becp_tmp, inter_bgrp_comm )

         DO k=1,kmax
            rsum=0.d0
            ibgrp_k = ibgrp_g2l( k )
            IF( ibgrp_k > 0 ) THEN
               DO is=1,nvb
                  DO iv=1,nh(is)
                     DO jv=1,nh(is)
                        DO ia=1,na(is)
                           inl=ish(is)+(iv-1)*na(is)+ia
                           jnl=ish(is)+(jv-1)*na(is)+ia
                           rsum = rsum + qq(iv,jv,is)*becp_tmp(inl)*becp(jnl,ibgrp_k)
                        END DO
                     END DO
                  END DO
               END DO
            END IF
            csc(k)=csc(k)+rsum
         END DO
!
         CALL mp_sum( csc( 1:kmax ), inter_bgrp_comm )
         WRITE( stdout,'("dotcsc =",12f18.15)') (csc(k),k=1,i)
!
      END DO
      WRITE( stdout,*)
!
      DEALLOCATE(becp)
      DEALLOCATE(cp_tmp)
      DEALLOCATE(becp_tmp)
!
      RETURN
      END SUBROUTINE dotcsc_x


!
!-----------------------------------------------------------------------
   FUNCTION enkin_x( c, f, n )
!-----------------------------------------------------------------------
      !
      ! calculation of kinetic energy term
      !
      USE kinds,              ONLY: DP
      USE constants,          ONLY: pi, fpi
      USE gvecw,              ONLY: ngw
      USE gvect,              ONLY: gstart
      USE gvecw,              ONLY: g2kin
      USE mp,                 ONLY: mp_sum
      USE mp_global,          ONLY: intra_bgrp_comm
      USE cell_base,          ONLY: tpiba2

      IMPLICIT NONE

      REAL(DP)                :: enkin_x

      ! input

      INTEGER,     INTENT(IN) :: n
      COMPLEX(DP), INTENT(IN) :: c( :, : )
      REAL(DP),    INTENT(IN) :: f( : )
      !
      ! local

      INTEGER  :: ig, i
      REAL(DP) :: sk(n)  ! automatic array
      !
      DO i=1,n
         sk(i)=0.0d0
         DO ig=gstart,ngw
            sk(i)=sk(i)+DBLE(CONJG(c(ig,i))*c(ig,i))*g2kin(ig)
         END DO
      END DO

      CALL mp_sum( sk(1:n), intra_bgrp_comm )

      enkin_x=0.0d0
      DO i=1,n
         enkin_x=enkin_x+f(i)*sk(i)
      END DO

      ! ... reciprocal-space vectors are in units of alat/(2 pi) so a
      ! ... multiplicative factor (2 pi/alat)**2 is required

      enkin_x = enkin_x * tpiba2
!
      RETURN
   END FUNCTION enkin_x

!-------------------------------------------------------------------------
      SUBROUTINE nlfl_bgrp_x( bec_bgrp, becdr_bgrp, lambda, descla, fion )
!-----------------------------------------------------------------------
!     contribution to fion due to the orthonormality constraint
! 
!
      USE kinds,             ONLY: DP
      USE io_global,         ONLY: stdout
      USE ions_base,         ONLY: na, nsp, nat
      USE uspp,              ONLY: nhsa=>nkb, qq
      USE uspp_param,        ONLY: nhm, nh, ish, nvb
      USE electrons_base,    ONLY: nspin, iupdwn, nupdwn, nbspx_bgrp, ibgrp_g2l, i2gupdwn_bgrp, nbspx, &
                                   iupdwn_bgrp, nupdwn_bgrp
      USE constants,         ONLY: pi, fpi
      USE descriptors,       ONLY: la_descriptor
      USE mp,                ONLY: mp_sum
      USE mp_global,         ONLY: intra_bgrp_comm, inter_bgrp_comm
!
      IMPLICIT NONE
      REAL(DP) :: bec_bgrp(:,:), becdr_bgrp(:,:,:)
      REAL(DP), INTENT(IN) :: lambda(:,:,:)
      TYPE(la_descriptor), INTENT(IN) :: descla(:)
      REAL(DP), INTENT(INOUT) :: fion(:,:)

!
      INTEGER :: k, is, ia, iv, jv, i, j, inl, isa, iss, nss, istart, ir, ic, nr, nc, ibgrp_i
      INTEGER :: n1, n2, m1, m2, nrcx
      REAL(DP), ALLOCATABLE :: temp(:,:), tmpbec(:,:),tmpdr(:,:) 
      REAL(DP), ALLOCATABLE :: fion_tmp(:,:)
      REAL(DP), ALLOCATABLE :: bec(:,:,:)
      REAL(DP), ALLOCATABLE :: becdr(:,:,:,:)
      REAL(DP), ALLOCATABLE :: bec_g(:,:)
      REAL(DP), ALLOCATABLE :: becdr_g(:,:,:)
      !
      CALL start_clock( 'nlfl' )
      !
      ALLOCATE( fion_tmp( 3, nat ) )
      !
      fion_tmp = 0.0d0
      !
      nrcx = MAXVAL( descla( : )%nrcx )
      !
      ALLOCATE( temp( nrcx, nrcx ), tmpbec( nhm, nrcx ), tmpdr( nrcx, nhm ) )
      ALLOCATE( bec( nhsa, nrcx, nspin ), becdr( nhsa, nrcx, nspin, 3 ) )

      ! redistribute bec, becdr according to the ortho subgroup
      ! this is required because they are combined with "lambda" matrixes
      
      DO iss = 1, nspin
         IF( descla( iss )%active_node > 0 ) THEN
            nss = nupdwn( iss )
            istart = iupdwn( iss )
            ic = descla( iss )%ic
            nc = descla( iss )%nc
            DO i=1,nc
               ibgrp_i = ibgrp_g2l( i+istart-1+ic-1 )
               IF( ibgrp_i > 0 ) THEN
                  bec( :, i, iss ) = bec_bgrp( :, ibgrp_i )
               ELSE
                  bec( :, i, iss ) = 0.0d0
               END IF
            END DO
            ir = descla( iss )%ir
            nr = descla( iss )%nr
            DO i=1,nr
               ibgrp_i = ibgrp_g2l( i+istart-1+ir-1 )
               IF( ibgrp_i > 0 ) THEN
                  becdr(:,i,iss,1) = becdr_bgrp( :, ibgrp_i, 1 )
                  becdr(:,i,iss,2) = becdr_bgrp( :, ibgrp_i, 2 )
                  becdr(:,i,iss,3) = becdr_bgrp( :, ibgrp_i, 3 )
               ELSE
                  becdr(:,i,iss,1) = 0.0d0
                  becdr(:,i,iss,2) = 0.0d0
                  becdr(:,i,iss,3) = 0.0d0
               END IF
            END DO
         ELSE
            bec(:,:,iss)   = 0.0d0
            becdr(:,:,iss,1) = 0.0d0
            becdr(:,:,iss,2) = 0.0d0
            becdr(:,:,iss,3) = 0.0d0
         END IF
      END DO

      CALL mp_sum( bec, inter_bgrp_comm )
      CALL mp_sum( becdr, inter_bgrp_comm )
      !
      DO k=1,3
         isa = 0
         DO is=1,nvb
            DO ia=1,na(is)
               isa = isa + 1
               !
               DO iss = 1, nspin
                  !
                  nss = nupdwn( iss )
                  istart = iupdwn( iss )
                  !
                  tmpbec = 0.d0
                  tmpdr  = 0.d0
                  !
                  IF( descla( iss )%active_node > 0 ) THEN
                     ! tmpbec distributed by columns
                     ic = descla( iss )%ic
                     nc = descla( iss )%nc
                     DO iv=1,nh(is)
                        DO jv=1,nh(is)
                           inl=ish(is)+(jv-1)*na(is)+ia
                           IF(ABS(qq(iv,jv,is)).GT.1.e-5) THEN
                              DO i=1,nc
                                 tmpbec(iv,i)=tmpbec(iv,i) + qq(iv,jv,is)*bec(inl,i,iss)
                              END DO
                           ENDIF
                        END DO
                     END DO
                     ! tmpdr distributed by rows
                     ir = descla( iss )%ir
                     nr = descla( iss )%nr
                     DO iv=1,nh(is)
                        inl=ish(is)+(iv-1)*na(is)+ia
                        DO i=1,nr
                           tmpdr(i,iv) = becdr( inl, i, iss, k )
                        END DO
                     END DO
                  END IF
                  !
                  IF(nh(is).GT.0)THEN
                     !
                     IF( descla( iss )%active_node > 0 ) THEN
                        ir = descla( iss )%ir
                        ic = descla( iss )%ic
                        nr = descla( iss )%nr
                        nc = descla( iss )%nc
                        CALL dgemm( 'N', 'N', nr, nc, nh(is), 1.0d0, tmpdr, nrcx, tmpbec, nhm, 0.0d0, temp, nrcx )
                        DO j = 1, nc
                           DO i = 1, nr
                              fion_tmp(k,isa) = fion_tmp(k,isa) + 2D0 * temp( i, j ) * lambda( i, j, iss )
                           END DO
                        END DO
                     END IF
!
                  ENDIF

               END DO
!
            END DO
         END DO
      END DO
      !
      DEALLOCATE( bec, becdr )
      DEALLOCATE( temp, tmpbec, tmpdr )
      !
      CALL mp_sum( fion_tmp, intra_bgrp_comm )
      !
      fion = fion + fion_tmp
      !
      DEALLOCATE( fion_tmp )
      !
      CALL stop_clock( 'nlfl' )
      !
      RETURN

      END SUBROUTINE nlfl_bgrp_x


!
!-----------------------------------------------------------------------
      SUBROUTINE pbc(rin,a1,a2,a3,ainv,rout)
!-----------------------------------------------------------------------
!
!     brings atoms inside the unit cell
!
      USE kinds,  ONLY: DP

      IMPLICIT NONE
! input
      REAL(DP) rin(3), a1(3),a2(3),a3(3), ainv(3,3)
! output
      REAL(DP) rout(3)
! local
      REAL(DP) x,y,z
!
! bring atomic positions to crystal axis
!
      x = ainv(1,1)*rin(1)+ainv(1,2)*rin(2)+ainv(1,3)*rin(3)
      y = ainv(2,1)*rin(1)+ainv(2,2)*rin(2)+ainv(2,3)*rin(3)
      z = ainv(3,1)*rin(1)+ainv(3,2)*rin(2)+ainv(3,3)*rin(3)
!
! bring x,y,z in the range between -0.5 and 0.5
!
      x = x - NINT(x)
      y = y - NINT(y)
      z = z - NINT(z)
!
! bring atomic positions back in cartesian axis
!
      rout(1) = x*a1(1)+y*a2(1)+z*a3(1)
      rout(2) = x*a1(2)+y*a2(2)+z*a3(2)
      rout(3) = x*a1(3)+y*a2(3)+z*a3(3)
!
      RETURN
      END SUBROUTINE pbc

!
!-------------------------------------------------------------------------
      SUBROUTINE prefor_x(eigr,betae)
!-----------------------------------------------------------------------
!
!     input :        eigr =  e^-ig.r_i
!     output:        betae_i,i(g) = (-i)**l beta_i,i(g) e^-ig.r_i 
!
      USE kinds,      ONLY : DP
      USE ions_base,  ONLY : nsp, na
      USE gvecw,      ONLY : ngw
      USE uspp,       ONLY : beta, nhtol
      USE uspp_param, ONLY : nh, ish
!
      IMPLICIT NONE
      COMPLEX(DP), INTENT(IN) :: eigr( :, : )
      COMPLEX(DP), INTENT(OUT) :: betae( :, : )
!
      INTEGER     :: is, iv, ia, inl, ig, isa
      COMPLEX(DP) :: ci
!
      CALL start_clock( 'prefor' )
      isa = 0
      DO is=1,nsp
         DO iv=1,nh(is)
            ci=(0.0d0,-1.0d0)**nhtol(iv,is)
            DO ia=1,na(is)
               inl=ish(is)+(iv-1)*na(is)+ia
               DO ig=1,ngw
                  betae(ig,inl)=ci*beta(ig,iv,is)*eigr(ig,ia+isa)
               END DO
            END DO
         END DO
         isa = isa + na(is)
      END DO
      CALL stop_clock( 'prefor' )
!
      RETURN
      END SUBROUTINE prefor_x

!------------------------------------------------------------------------
    SUBROUTINE collect_bec_x( bec_repl, bec_dist, desc, nspin )
!------------------------------------------------------------------------
       USE kinds,       ONLY : DP
       USE mp_global,   ONLY : intra_bgrp_comm
       USE mp,          ONLY : mp_sum
       USE descriptors, ONLY : la_descriptor
       USE io_global,   ONLY : stdout
       REAL(DP), INTENT(OUT) :: bec_repl(:,:)
       REAL(DP), INTENT(IN)  :: bec_dist(:,:)
       TYPE(la_descriptor), INTENT(IN)  :: desc(:)
       INTEGER,  INTENT(IN)  :: nspin
       INTEGER :: i, ir, n, nrcx, iss
       !
       bec_repl = 0.0d0
       !
       !  bec is distributed across row processor, the first column is enough
       !
       IF( desc( 1 )%active_node > 0 .AND. ( desc( 1 )%myc == 0 ) ) THEN
          ir = desc( 1 )%ir
          DO i = 1, desc( 1 )%nr
             bec_repl( :, i + ir - 1 ) = bec_dist( :, i )
          END DO
          IF( nspin == 2 ) THEN
             n  = desc( 1 )%n   ! number of states with spin==1 ( nupdw(1) )
             nrcx = desc( 1 )%nrcx ! array elements reserved for each spin ( bec(:,2*nrcx) )
             ir = desc( 2 )%ir
             DO i = 1, desc( 2 )%nr
                bec_repl( :, i + ir - 1 + n ) = bec_dist( :, i + nrcx )
             END DO
          END IF
       END IF
       !
       CALL mp_sum( bec_repl, intra_bgrp_comm )
       !
       RETURN
    END SUBROUTINE collect_bec_x

!------------------------------------------------------------------------
    SUBROUTINE distribute_lambda_x( lambda_repl, lambda_dist, desc )
!------------------------------------------------------------------------
       USE kinds,       ONLY : DP
       USE descriptors
       REAL(DP), INTENT(IN)  :: lambda_repl(:,:)
       REAL(DP), INTENT(OUT) :: lambda_dist(:,:)
       TYPE(la_descriptor), INTENT(IN)  :: desc
       INTEGER :: i, j, ic, ir
       IF( desc%active_node > 0 ) THEN
          ir = desc%ir
          ic = desc%ic
          DO j = 1, desc%nc
             DO i = 1, desc%nr
                lambda_dist( i, j ) = lambda_repl( i + ir - 1, j + ic - 1 )
             END DO
          END DO
       END IF
       RETURN
    END SUBROUTINE distribute_lambda_x
    !
!------------------------------------------------------------------------
    SUBROUTINE distribute_bec_x( bec_repl, bec_dist, desc, nspin )
!------------------------------------------------------------------------
       USE kinds,       ONLY : DP
       USE descriptors
       REAL(DP), INTENT(IN)  :: bec_repl(:,:)
       REAL(DP), INTENT(OUT) :: bec_dist(:,:)
       TYPE(la_descriptor), INTENT(IN)  :: desc(:)
       INTEGER,  INTENT(IN)  :: nspin
       INTEGER :: i, ir, n, nrcx
       !
       IF( desc( 1 )%active_node > 0 ) THEN
          !
          bec_dist = 0.0d0
          !
          ir = desc( 1 )%ir
          DO i = 1, desc( 1 )%nr
             bec_dist( :, i ) = bec_repl( :, i + ir - 1 )
          END DO
          !
          IF( nspin == 2 ) THEN
             n     = desc( 1 )%n  !  number of states with spin 1 ( nupdw(1) )
             nrcx  = desc( 1 )%nrcx   !  array elements reserved for each spin ( bec(:,2*nrcx) )
             ir = desc( 2 )%ir
             DO i = 1, desc( 2 )%nr
                bec_dist( :, i + nrcx ) = bec_repl( :, i + ir - 1 + n )
             END DO
          END IF
          !
       END IF
       RETURN
    END SUBROUTINE distribute_bec_x
    !
!------------------------------------------------------------------------
    SUBROUTINE distribute_zmat_x( zmat_repl, zmat_dist, desc )
!------------------------------------------------------------------------
       USE kinds,       ONLY : DP
       USE descriptors
       REAL(DP), INTENT(IN)  :: zmat_repl(:,:)
       REAL(DP), INTENT(OUT) :: zmat_dist(:,:)
       TYPE(la_descriptor), INTENT(IN)  :: desc
       INTEGER :: i, ii, j, me, np
       me = desc%mype
       np = desc%npc * desc%npr
       IF( desc%active_node > 0 ) THEN
          DO j = 1, desc%n
             ii = me + 1
             DO i = 1, desc%nrl
                zmat_dist( i, j ) = zmat_repl( ii, j )
                ii = ii + np
             END DO
          END DO
       END IF
       RETURN
    END SUBROUTINE distribute_zmat_x
    !
!------------------------------------------------------------------------
    SUBROUTINE collect_lambda_x( lambda_repl, lambda_dist, desc )
!------------------------------------------------------------------------
       USE kinds,       ONLY : DP
       USE mp_global,   ONLY: intra_bgrp_comm
       USE mp,          ONLY: mp_sum
       USE descriptors
       REAL(DP), INTENT(OUT) :: lambda_repl(:,:)
       REAL(DP), INTENT(IN)  :: lambda_dist(:,:)
       TYPE(la_descriptor), INTENT(IN)  :: desc
       INTEGER :: i, j, ic, ir
       lambda_repl = 0.0d0
       IF( desc%active_node > 0 ) THEN
          ir = desc%ir
          ic = desc%ic
          DO j = 1, desc%nc
             DO i = 1, desc%nr
                lambda_repl( i + ir - 1, j + ic - 1 ) = lambda_dist( i, j )
             END DO
          END DO
       END IF
       CALL mp_sum( lambda_repl, intra_bgrp_comm )
       RETURN
    END SUBROUTINE collect_lambda_x
    !
!------------------------------------------------------------------------
    SUBROUTINE collect_zmat_x( zmat_repl, zmat_dist, desc )
!------------------------------------------------------------------------
       USE kinds,       ONLY : DP
       USE mp_global,   ONLY: intra_bgrp_comm
       USE mp,          ONLY: mp_sum
       USE descriptors
       REAL(DP), INTENT(OUT) :: zmat_repl(:,:)
       REAL(DP), INTENT(IN)  :: zmat_dist(:,:)
       TYPE(la_descriptor), INTENT(IN)  :: desc
       INTEGER :: i, ii, j, me, np, nrl
       zmat_repl = 0.0d0
       me = desc%mype
       np = desc%npc * desc%npr
       nrl = desc%nrl
       IF( desc%active_node > 0 ) THEN
          DO j = 1, desc%n
             ii = me + 1
             DO i = 1, nrl
                zmat_repl( ii, j ) = zmat_dist( i, j )
                ii = ii + np
             END DO
          END DO
       END IF
       CALL mp_sum( zmat_repl, intra_bgrp_comm )
       RETURN
    END SUBROUTINE collect_zmat_x
    !
!------------------------------------------------------------------------
    SUBROUTINE setval_lambda_x( lambda_dist, i, j, val, desc )
!------------------------------------------------------------------------
       USE kinds,       ONLY : DP
       USE descriptors
       REAL(DP), INTENT(OUT) :: lambda_dist(:,:)
       INTEGER,  INTENT(IN)  :: i, j
       REAL(DP), INTENT(IN)  :: val
       TYPE(la_descriptor), INTENT(IN)  :: desc
       IF( desc%active_node > 0 ) THEN
          IF( ( i >= desc%ir ) .AND. ( i - desc%ir + 1 <= desc%nr ) ) THEN
             IF( ( j >= desc%ic ) .AND. ( j - desc%ic + 1 <= desc%nc ) ) THEN
                lambda_dist( i - desc%ir + 1, j - desc%ic + 1 ) = val
             END IF
          END IF
       END IF
       RETURN
    END SUBROUTINE setval_lambda_x


!------------------------------------------------------------------------

