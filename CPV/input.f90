!
! Copyright (C) 2002 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE input_cp

   IMPLICIT NONE
   SAVE

   PRIVATE

   PUBLIC :: iosys

CONTAINS

!-----------------------------------------------------------------------
      subroutine iosys( nbeg_ , ndr_ , ndw_ , nomore_ , iprint_                       &
     & , delt_ , emass_ , emaec_  , tsde_ , frice_ , grease_ , twall_                        &
     & , tortho_ , eps_ , max_ , trane_ , ampre_ , tranp_ , amprp_                           &
     & , tfor_ , tsdp_ , fricp_ , greasp_ , tcp_ , tcap_ , tolp_ , trhor_ , trhow_ , tvlocw_ &
     & , tnosep_ , qnp_ , tempw_ , tnosee_ , qne_ , ekincw_                                &
     & , tpre_ , thdyn_ , thdiag_ , twmass_ , wmass_ , frich_ , greash_ , press_           &
     & , tnoseh_ , qnh_ , temph_ , celldm_ , ibrav_ , tau0_ , ecutw_ , ecut_ , iforce_ &
     & , nat_ , nsp_ , na_ , pmass_ , rcmax_ , ipp_ , f_ , nel_ , nspin_ , nupdwn_  &
     & , iupdwn_ , n_ , nx_, nr1_ , nr2_ , nr3_ , omega_ , alat_ , a1_ , a2_ , a3_  & 
     & , nr1b_ , nr2b_ , nr3b_ , nr1s_ , nr2s_ , nr3s_ , agg_ , sgg_ , e0gg_ &
     & , psfile_ , pseudo_dir_, iprsta_, ispin_ )

!-----------------------------------------------------------------------
!   this subroutine reads control variables from standard input (unit 5)
!     ------------------------------------------------------------------

      use input_parameters, only: &
           nr1, nr2, nr3, greash, press, nr2s, nr3s, nr1s, tolp, temph, grease, &
           tempw, fnoseh, amprp, greasp, twall, tranp, atomic_positions, nelec, &
           if_pos, rd_ht, nelup, neldw, occupations, f_inp, rd_pos, nr3b, pseudo_dir, &
           nr1b, nr2b, atom_ptyp, sp_pos, atom_mass, atom_pfile, iprint, orthogonalization, &
           electron_velocities, startingwfc, ndr, ndw, ion_dynamics, ion_damping, &
           cell_velocities, electron_dynamics, electron_damping, ion_velocities, &
           celldm, nbnd, nspin, calculation, ntyp, ibrav, restart_mode, ion_positions, &
           nstep, ecutwfc, ecutrho, ampre, ortho_eps, ortho_max, wmass, qcutz, q2sigma, &
           ecfixed, ekincw, fnosep, nat, tstress, disk_io, fnosee, ion_temperature, &
           cell_temperature, cell_dofree, cell_dynamics, cell_damping, electron_temperature, &
           dt, emass, emass_cutoff, ion_radius, isave, verbosity

      use read_namelists_module, only: read_namelists
      use read_cards_module, only: read_cards

      use constants, only: pi, scmass, factem, eps8
      use parameters, only: nsx, natx, nbndxx
      use io_global, only: ionode
      use mp, only: mp_bcast

      !
      implicit none
      !
      !
      real(kind=8) :: ampre_ , delt_ , ekincw_ , emass_ , emaec_ , eps_ ,       &
     &       frice_ , fricp_ , frich_ , grease_ , greasp_ , greash_ ,        &
     &       press_ , qnp_ , qne_ , qnh_ , tempw_ , temph_ , tolp_ , wmass_ ,    &
             amprp_ ( nsx ), celldm_ ( 6 ), tau0_ ( 3, natx, nsx ), ecut_ , ecutw_

      integer :: nbeg_ , ndr_ , ndw_ , nomore_ , iprint_ , max_ , iforce_( 3, natx, nsx )

      logical :: trane_ , tsde_ , twall_ , tortho_ , tnosee_ , tfor_ , tsdp_ , tcp_ , &
           tcap_ , tnosep_ , trhor_ , trhow_ , tvlocw_ , tpre_ , thdyn_ , thdiag_ ,   &
           twmass_ , tnoseh_ , tranp_ ( nsx )

      integer :: nat_ , nsp_ , na_ ( nsx ), ipp_ ( nsx ), nel_ ( 2 ), nspin_ , &
     &     nupdwn_ ( 2 ), iupdwn_ ( 2 ), n_ , nx_ , nr1_ , nr2_ , nr3_ , &
     &     nr1b_ , nr2b_ , nr3b_ , nr1s_ , nr2s_ , nr3s_ , ibrav_, iprsta_

      real(kind=8) :: pmass_ ( nsx ), rcmax_ ( nsx ), f_ ( nbndxx ), ispin_ ( nbndxx ), &
     &     omega_ , alat_ , a1_ ( 3 ), a2_ ( 3 ), a3_ ( 3 ), agg_ , sgg_ , e0gg_


      character(len=80) :: psfile_ ( nsx ) , pseudo_dir_

      !
      ! local variables
      !

      real(kind=8), parameter:: terahertz = 2.418D-5
      real(kind=8) :: taus( 3, natx, nsx ), ocp, fsum
      integer :: unit = 5, ionode_id = 0, i, ia, ios, is, iss, in


      CALL read_namelists( 'CP' )

      CALL read_cards( 'CP' )

      IF( TRIM( calculation ) == 'nscf' ) trhor_ = .true.
     
!
! translate from input to internals of CP, various checks

      ! ...   Set the number of species

      nsp_ = ntyp

      ! ...   IBRAV and CELLDM

      ibrav_  = ibrav
      celldm_ = celldm

      ! ...   Set Values for bands and spin

      n_     = nbnd * nspin
      nspin_ = nspin

      ! ...   Set Values for the cutoff

      ecutw_ = ecutwfc
      if ( ecutrho <= 0.d0 ) ecutrho = 4.d0 * ecutwfc
      ecut_ = ecutrho

      ! ...   nbeg

      ampre_ = ampre
      SELECT CASE ( restart_mode ) 
         CASE ('from_scratch')
            nbeg_ = -2
            if ( ion_positions == 'from_input' ) nbeg_ = -1
            nomore_ = nstep
            trane_  = ( startingwfc == 'random' )
            if ( ampre_ == 0.d0 ) ampre_ = 0.02
         CASE ('reset_counters')
            nbeg_ = 0
            nomore_ = nstep
         CASE ('restart')
            nbeg_ = 1
            nomore_ = nstep
         CASE DEFAULT
            CALL errore(' iosys ',' unknown restart_mode '//trim(restart_mode), 1 )
      END SELECT

      ndr_ = ndr
      ndw_ = ndw
      iprint_ = iprint


      ! ...   TORTHO

      SELECT CASE ( orthogonalization ) 
      CASE ('Gram-Schmidt')
         tortho_ = .FALSE.
      CASE ('ortho')
         tortho_ = .TRUE.
      CASE DEFAULT
         CALL errore(' iosys ',' unknown orthogonalization '//&
              trim(orthogonalization), 1 )
      END SELECT


      SELECT CASE ( electron_velocities ) 
      CASE ('default')
         continue
      CASE ('zero')
         print '("Warning: electron_velocities keyword has no effect")'
      CASE DEFAULT
         CALL errore(' iosys ',' electron_velocities='// &
              trim(electron_velocities)//' not implemented', 1 )
      END SELECT

      ! ...   TSDE

      SELECT CASE ( electron_dynamics ) 
      CASE ('sd')
         tsde_  = .TRUE.
         frice_ = 0.d0
      CASE ('verlet')
         tsde_  = .FALSE.
         frice_ = 0.d0
      CASE ('damp')
         tsde_  = .FALSE.
         frice_ = electron_damping
      CASE ('none')
         tsde_  = .FALSE.
         frice_ = 0.d0
      CASE DEFAULT
         CALL errore(' iosys ',' unknown electron_dynamics '//&
              trim(electron_dynamics),1)
      END SELECT

      ! Ion velocities

      SELECT CASE ( ion_velocities ) 
      CASE ('default')
         tcap_ = .false.
      CASE ('random')
         tcap_ = .true.
      CASE ('zero')
         print '("Warning: ion_velocities = zero not yet implemented")'
      CASE DEFAULT
         CALL errore(' iosys ',' unknown ion_velocities '//trim(ion_velocities),1)
      END SELECT

      ! ...   TFOR TSDP

      SELECT CASE ( ion_dynamics ) 
      CASE ('sd')
         tsdp_ = .TRUE.
         tfor_ = .TRUE.
         fricp_= 0.d0
      CASE ('verlet')
         tsdp_ = .FALSE.
         tfor_ = .TRUE.
         fricp_= 0.d0
      CASE ('damp')
         tsdp_ = .FALSE.
         tfor_ = .TRUE.
         fricp_= ion_damping
      CASE ('none')
         tsdp_ = .FALSE.
         tfor_ = .FALSE.
         fricp_= 0.d0
      CASE DEFAULT
         CALL errore(' iosys ',' unknown ion_dynamics '//trim(ion_dynamics), 1 )
      END SELECT

      !

      SELECT CASE ( cell_velocities ) 
      CASE ('default')
         continue
      CASE ('zero')
         print '("Warning: cell_velocities = zero not yet implemented")'
      CASE DEFAULT
         CALL errore(' iosys ',' unknown cell_velocities '//trim(cell_velocities),1)
      END SELECT

      !
      
      SELECT CASE ( cell_dynamics ) 
      CASE ('sd')
         tpre_ = .TRUE.
         thdyn_= .TRUE.
         frich_= 0.d0
      CASE ('pr')
         tpre_ = .TRUE.
         thdyn_= .TRUE.
         frich_= 0.d0
      CASE ('damp-pr')
         tpre_  = .TRUE.
         thdyn_= .TRUE.
         frich_ = cell_damping
      CASE ('none')
         tpre_ = .FALSE.
         thdyn_= .FALSE.
         frich_= 0.d0
      CASE DEFAULT
         CALL errore(' iosys ',' unknown cell_dynamics '//trim(cell_dynamics), 1 )
      END SELECT

      !

      SELECT CASE ( electron_temperature ) 
         !         temperature control of electrons via Nose' thermostat
         !         EKINW (REAL(DBL))  average kinetic energy (in atomic units)
         !         FNOSEE (REAL(DBL))  frequency (in terahertz)
      CASE ('nose')
         tnosee_ = .TRUE.
      CASE ('not_controlled')
         tnosee_ = .FALSE.
      CASE DEFAULT
         CALL errore(' iosys ',' unknown electron_temperature '//&
              trim(electron_temperature), 1 )
      END SELECT

      !

      SELECT CASE ( ion_temperature ) 
         !         temperature control of ions via Nose' thermostat
         !         TEMPW (REAL(DBL))  frequency (in which units?)
         !         FNOSEP (REAL(DBL))  temperature (in which units?)
      CASE ('nose')
         tnosep_ = .TRUE.
         tcp_ = .false.
      CASE ('not_controlled')
         tnosep_ = .FALSE.
         tcp_ = .false.
      CASE ('rescaling' )
         tnosep_ = .FALSE.
         tcp_ = .true.
      CASE DEFAULT
         CALL errore(' iosys ',' unknown ion_temperature '//&
              trim(ion_temperature), 1 )
      END SELECT


      SELECT CASE ( cell_temperature ) 
         !         cell temperature control of ions via Nose' thermostat
         !         FNOSEH (REAL(DBL))  frequency (in which units?)
         !         TEMPH (REAL(DBL))  temperature (in which units?)
      CASE ('nose')
         tnoseh_ = .TRUE.
      CASE ('not_controlled')
         tnoseh_ = .FALSE.
      CASE DEFAULT
         CALL errore(' iosys ',' unknown cell_temperature '//&
              trim(cell_temperature), 1 )
      END SELECT


      SELECT CASE ( cell_dofree )
      CASE ('all')
         thdiag_ =.false.
      CASE ('xyz')
         thdiag_ =.true.
      CASE DEFAULT
         CALL errore(' iosys ',' unknown cell_dofree '//trim(cell_dofree), 1 )
      END SELECT


      ! ...  radii, masses

      DO is = 1, nsp_
         rcmax_ (is) = ion_radius(is)
         IF( ion_radius(is) <= 0.d0 ) THEN
            CALL errore(' iosys ',' invalid  ion_radius ', is) 
         END IF
      END DO

      !
      ! compatibility between FPMD and CP90
      !
      iprint_ = isave 
      if ( trim( verbosity ) == 'high' ) then
         iprsta_ = 3
      else
         iprsta_ = 1
      end if
      delt_   = dt
      emass_ = emass
      emaec_  = emass_cutoff
      agg_ = qcutz
      sgg_ = q2sigma
      e0gg_= ecfixed
      eps_ = ortho_eps
      max_ = ortho_max
      ! wmass is calculated in "init"
      wmass_ = wmass
      twmass_ = ( wmass == 0.d0 )
      if ( tstress ) tpre_ = .true.
      trhow_ = ( trim( disk_io ) == 'high' )
      tvlocw_ = .false. ! temporaneo
      !
      qne_ = 0.0d0
      qnp_ = 0.0d0
      qnh_ = 0.0d0
      if( fnosee > 0.0d0 ) qne_ = 4.d0*ekincw/(fnosee*(2.d0*pi)*terahertz)**2
      if( fnosep > 0.0d0 ) qnp_ = 2.d0*(3*nat)*tempw/factem/(fnosep*(2.d0*pi)*terahertz)**2
      if( fnoseh > 0.0d0 ) qnh_ = 2.d0*(3*3  )*temph/factem/(fnoseh*(2.d0*pi)*terahertz)**2
      tempw_ = tempw
      temph_ = temph
      ekincw_ = ekincw

      grease_ = grease
      twall_ = twall
      tranp_ ( 1 : nsp_ ) =  tranp ( 1 : nsp_ )
      amprp_ ( 1 : nsp_ ) =  amprp ( 1 : nsp_ )
 
      greasp_ = greasp
      tolp_ = tolp
      greash_ = greash
      press_ = press

      nr1_ = nr1
      nr2_ = nr2
      nr3_ = nr3
      nr1s_ = nr1s
      nr2s_ = nr2s
      nr3s_ = nr3s
      nr1b_ = nr1b
      nr2b_ = nr2b
      nr3b_ = nr3b

      nat_ = nat
      pseudo_dir_ = pseudo_dir
      
      ! cards parameters

      tau0_  = 0.0
      iforce_= 0
      ipp_   = 0
      psfile_= ' '
      a1_    = 0.0
      a2_    = 0.0
      a3_    = 0.0

      pmass_ ( 1:nsp_ ) = atom_mass( 1:nsp_ )
      psfile_ ( 1:nsp_ ) = atom_pfile( 1:nsp_ )
      ipp_ ( 1:nsp_ ) = atom_ptyp( 1:nsp_ )

      na_ = 0
      do is = 1, nsp_
          do ia = 1, nat_
             if ( sp_pos(ia) == is) then
                na_(is) = na_(is) + 1
                if( na_(is) > natx ) call errore(' cards',' na > natx', na_ (is) )
                do i = 1, 3
                   tau0_ (i, na_ (is), is ) = rd_pos(i, ia)
                   iforce_ (i, na_ (is), is ) = if_pos(i, ia)
                end do
             end if
          end do
       end do

      !
      ! set up atomic positions and crystal lattice
      !
      if ( ibrav_ == 0 ) then
         a1_ = rd_ht( 1, 1:3 )
         a2_ = rd_ht( 2, 1:3 )
         a3_ = rd_ht( 3, 1:3 )
         if ( celldm_ (1) == 0.d0 ) then
            celldm_ (1) = sqrt( a1_ (1) ** 2 + a1_ (2) ** 2 + a1_ (3) ** 2 )
            a1_(:) = a1_(:) / celldm_(1)
            a2_(:) = a2_(:) / celldm_(1)
            a3_(:) = a3_(:) / celldm_(1)
         end if
      else
         call latgen( ibrav_ , celldm_ , a1_ , a2_ , a3_ , omega_ )
      end if
      alat_ = celldm_ (1)

      !
      SELECT CASE ( atomic_positions ) 
         !
         !  convert input atomic positions to internally used format:
         !  tau0 in atomic units
         !
         CASE ('alat')
            !
            !  input atomic positions are divided by a0
            !
            tau0_ = tau0_ * alat_
         CASE ('bohr')
            !
            !  input atomic positions are in a.u.: do nothing
            !
            continue
         CASE ('crystal')
            !
            !  input atomic positions are in crystal axis ("scaled"):
            !
            taus = tau0_
            do is = 1, nsp_
               do ia = 1, na_(is)
                  do i = 1, 3
                     tau0_ ( i, ia, is ) = a1_ (i) * taus( 1, ia, is) &
                                         + a2_ (i) * taus( 2, ia, is) &
                                         + a3_ (i) * taus( 3, ia, is)
                  end do
               end do
            end do
         CASE ('angstrom')
            !
            !  atomic positions in A
            !
            tau0_ = tau0_ / 0.529177
         CASE DEFAULT
            CALL errore(' iosys ',' atomic_positions='//trim(atomic_positions)// &
                 ' not implemented ', 1 )
      END SELECT

      !
      !  set occupancies
      !
      
      IF( nelec < 1 ) THEN
         CALL errore(' iosys ',' nelec less than 1 ', int(nelec) )
      END IF
      IF( nint(nelec) - nelec > eps8 ) THEN
         CALL errore(' iosys ',' nelec must be integer', int(nelec) )
      END IF

      if( mod( n_ , 2 ) .ne. 0 ) then
         nx_ = n_ + 1
      else
         nx_= n_
      end if

      iupdwn_ ( 1 ) = 1
      nel_ = 0

      SELECT CASE ( TRIM(occupations) ) 
      CASE ('bogus')
         !
         ! empty-states calculation: occupancies have a (bogus) finite value
         !
         ! bogus to ensure \sum_i f_i = Nelec  (nelec is integer)
         !
         f_ ( : ) = nelec / n_         
         nel_ (1) = nint(nelec)
         nupdwn_ (1) = n_
         if ( nspin_ == 2 ) then
            !
            ! bogus to ensure Nelec = Nup + Ndw
            !
            nel_ (1) = ( nint(nelec) + 1 ) / 2
            nel_ (2) =   nint(nelec)       / 2
            nupdwn_ (1)=nbnd
            nupdwn_ (2)=nbnd
            iupdwn_ (2)=nbnd+1
         end if
      CASE ('from_input')
         !
         ! occupancies have been read from input
         !
         f_ ( 1:nbnd ) = f_inp( 1:nbnd, 1 )
         if( nspin_ == 2 ) f_ ( nbnd+1 : 2*nbnd ) = f_inp( 1:nbnd, 2 ) 
         if( nelec == 0.d0 ) nelec = SUM ( f_ ( 1:n_ ) )
         if( nspin_ == 2 .and. nelup == 0) nelup = SUM ( f_ ( 1:nbnd ) )
         if( nspin_ == 2 .and. neldw == 0) neldw = SUM ( f_ ( nbnd+1 : 2*nbnd ) )

         if( nspin_ == 1 ) then 
           nel_ (1) = nint(nelec)
           nupdwn_ (1) = n_
         else
           IF ( ABS (nelup + neldw - nelec) > eps8 ) THEN
              CALL errore(' iosys ',' wrong # of up and down spin', 1 )
           END IF
           nel_ (1) = nint(nelup)
           nel_ (2) = nint(neldw)
           nupdwn_ (1)=nbnd
           nupdwn_ (2)=nbnd
           iupdwn_ (2)=nbnd+1
         end if

      CASE ('fixed')

         if( nspin_ == 1 ) then
            nel_ (1) = nint(nelec)
            nupdwn_ (1) = n_
         else
            IF ( nelup + neldw /= nelec  ) THEN
               CALL errore(' iosys ',' wrong # of up and down spin', 1 )
            END IF
            nel_ (1) = nint(nelup)
            nel_ (2) = nint(neldw)
            nupdwn_ (1)=nbnd
            nupdwn_ (2)=nbnd
            iupdwn_ (2)=nbnd+1
         end if

         ! ocp = 2 for spinless systems, ocp = 1 for spin-polarized systems
         ocp = 2.d0 / nspin_
         ! default filling: attribute ocp electrons to each states
         !                  until the good number of electrons is reached
         do iss = 1, nspin_
            fsum = 0.0d0
            do in = iupdwn_ ( iss ), iupdwn_ ( iss ) - 1 + nupdwn_ ( iss )
               if ( fsum + ocp < nel_ ( iss ) + 0.0001 ) then
                  f_ (in) = ocp
               else
                  f_ (in) = max( nel_ ( iss ) - fsum, 0.d0 )
               end if
                fsum=fsum + f_(in)
            end do
         end do

      CASE DEFAULT
         CALL errore(' iosys ',' occupation method not implemented', 1 )
      END SELECT

      do iss = 1, nspin_
         do in = iupdwn_(iss), iupdwn_(iss) - 1 + nupdwn_(iss)
            ispin_(in) = iss
         end do
      end do

!
!     --------------------------------------------------------
!     print out heading
!
      write(6,500) nbeg_ , nomore_ , iprint_ , ndr_ , ndw_
      write(6,505) delt_
      write(6,510) emass_ , emaec_
!
      if( tortho_ ) then
         write(6,511) eps_ , max_
      else
         write(6,512)
      endif
!
      if( tsde_ ) then
         write(6,513)
      else
         if ( tnosee_ ) frice_ = 0.
         write(6,509)
         write(6,514) frice_ , grease_
      endif
!
      if ( trhor_ ) then
         write(6,720)
      endif
!
      if( .not. trhor_ .and. trhow_ )then
         write(6,721)
      endif
!
      if( tvlocw_ )then
         write(6,722)
      endif
!
      if( trane_ ) then
         write(6,515) ampre_
      endif
      write(6,516)
      do is =1, nsp_
         if(tranp_(is)) write(6,517) is, amprp_(is)
      end do
!
      if(tfor_) then
         if(tnosep_) fricp_ = 0.
         write(6,520)
         if(tsdp_)then
            write(6,521)
         else
            write(6,522) fricp_ , greasp_
         endif
      else
         write(6,518)
      endif
!
      if( tfor_ ) then
         if(( tcp_ .or. tcap_ .or. tnosep_ ) .and. tsdp_ ) then
            call errore(' main',' t contr. for ions when tsdp=.t.',0)
         endif
         if(.not. tcp_ .and. .not. tcap_ .and. .not. tnosep_ ) then
            write(6,550)
         else if(tcp_ .and. tcap_ ) then
            call errore(' main',' tcp and tcap both true',0)
         else if(tcp_ .and. tnosep_ ) then
            call errore(' main',' tcp and tnosep both true',0)
         else if(tcap_ .and. tnosep_ ) then
            call errore(' main',' tcap and tnosep both true',0)
         else if(tcp_ ) then
            write(6,555) tempw_ , tolp_
         else if(tcap_) then
            write(6,560) tempw_ , tolp_
         else if(tnosep_ ) then
            write(6,562) tempw_ , qnp_
         end if
         if(tnosee_) then
            write(6,566) ekincw_ , qne_
         end if
      end if
!
      if(tpre_) then
         write(6,600)
         if(thdyn_) then
            if(thdiag_) write(6,608)
            if(tnoseh_) then
               frich_=0.
               write(6,604) temph_,qnh_,press_
            else
               write(6,602) frich_,greash_,press_
            endif
         else
            write(6,606)
         endif
      endif
      if ( agg_ .ne. 0.d0) then
            write(6,650) agg_, sgg_, e0gg_
      end if
      write(6,700) iprsta_

!     
 500  format(//                                                         &
     &       ' nbeg=',i3,' nomore=',i7,3x,' iprint=',i4,/               &
     &       ' reads from',i3,' writes on',i3)
 505  format(' time step = ',f9.4/)
 510  format(' parameters for electron dynamics:'/                      &
     &       ' emass= ',f10.2,2x,'emaec= ',f10.2,'ry')
 511  format(' orthog. with lagrange multipliers: eps=',e10.2,          &
     &         ' max=',i3)
 512  format(' orthog. with gram-schmidt')
 513  format(' electron dynamics with steepest descent')
 509  format(' verlet algorithm for electron dynamics')
 514  format(' with friction frice = ',f7.4,' , grease = ',f7.4)
 720  format(' charge density is read from unit 47',/)
 721  format(' charge density is written in unit 47',/)
 722  format(' local potential is written in unit 46',/)
 515  format(' initial random displacement of el. coordinates with ',   &
     &       ' amplitude=',f10.6,/                                      &
     &       ' trane not to be used with mass preconditioning')
 516  format(/)
 517  format(' initial random displacement of ionic coord. for species ',&
     &       i4,' : amplitude=',f10.6)
 518  format(' ions are not allowed to move'/)
 520  format(' ions are allowed to move')
 521  format(' ion dynamics with steepest descent')
 522  format(' ion dynamics with fricp = ',f7.4,' and greasp = ',f7.4)
 550  format(' ion dynamics: the temperature is not controlled'//)
 555  format(' ion dynamics with rescaling of velocities:'/             &
     &       ' temperature required=',f10.5,'(kelvin)',' tolerance=',   &
     &       f10.5//)
 560  format(' ion dynamics with canonical temp. control:'/             &
     &       ' temperature required=',f10.5,'(kelvin)',' tolerance=',   &
     &       f10.5//)
 562  format(' ion dynamics with nose` temp. control:'/                 &
     &       ' temperature required=',f10.5,'(kelvin)',' nose` mass = ',&
     &       f10.3//)
 566  format(' electronic dynamics with nose` temp. control:'/          &
     &       ' elec. kin. en. required=',f10.5,'(hartree)',             &
     &       ' nose` mass = ',f10.3//)
 600  format(' internal stress tensor calculated')
 602  format(' cell parameters dynamics with frich = ',f7.4,            &
     &       ' and greash = ',f7.4,/                                    &
     &       ' external pressure = ',f11.7,'(gpa)'//)
 604  format(' cell parameters dynamics with nose` temp. control:'/     &
     &       ' cell temperature required = ',f10.5,'(kelvin)',          &
     &       ' nose` mass = ',f10.3,/                                   &
     &       ' external pressure = ',f11.7,'(gpa)'//)
 606  format(' cell parameters are not allowed to move'//)
 608  format(' frozen off-diagonal cell parameters'//)
 650  format(' modified kinetic energy functional, with parameters:'/   &
           & ' agg = ',f8.4,'  sgg = ', f7.4,'  e0gg = ',f6.2)
 700  format(' iprsta = ',i2/)
      return
      end subroutine
!
END MODULE input_cp
