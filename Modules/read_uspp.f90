!
! Copyright (C) 2006-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
MODULE read_uspp_module
  !---------------------------------------------------------------------
  !
  !  routines reading ultrasoft pseudopotentials in older formats:
  !  Vanderbilt's code and Andrea's RRKJ3 format
  !     
  USE kinds, ONLY: DP
  USE parameters, ONLY: lmaxx, lqmax
  USE io_global, ONLY: stdout
  USE funct, ONLY: set_dft_from_name, dft_is_hybrid, dft_is_meta, &
       set_dft_from_indices
  USE matrix_inversion
  !
  ! Variables above are not modified, variables below are
  !
  USE uspp_param, ONLY: oldvan
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  PUBLIC :: readvan, readrrkj
  !
CONTAINS
  !---------------------------------------------------------------------
  subroutine readvan( iunps, is, upf )
    !---------------------------------------------------------------------
    !
    !     Read Vanderbilt pseudopotential from unit "iunps"
    !     for species "is" into the structure "upf"
    !     info on DFT level in module "funct"
    !
    !     ------------------------------------------------------
    !     Important:
    !     ------------------------------------------------------
    !     The order of all l-dependent objects is always s,p,d
    !     ------------------------------------------------------
    !     potentials, e.g. vloc_at, are really r*v(r)
    !     wave funcs, e.g. chi, are really proportional to r*psi(r)
    !     and are normalized so int (chi**2) dr = 1
    !     thus psi(r-vec)=(1/r)*chi(r)*y_lm(theta,phi)
    !     conventions carry over to beta, etc
    !     charge dens, e.g. rho_atc, really 4*pi*r**2*rho
    !
    !     ------------------------------------------------------
    !     Notes on qfunc and qfcoef:
    !     ------------------------------------------------------
    !     Since Q_ij(r) is the product of two orbitals like
    !     psi_{l1,m1}^star * psi_{l2,m2}, it can be decomposed by
    !     total angular momentum L, where L runs over | l1-l2 | ,
    !     | l1-l2 | +2 , ... , l1+l2.  (L=0 is the only component
    !     needed by the atomic program, which assumes spherical
    !     charge symmetry.)
    !
    !     Recall  qfunc(r) = y1(r) * y2(r)  where y1 and y2 are the
    !     radial parts of the wave functions defined according to
    !
    !       psi(r-vec) = (1/r) * y(r) * Y_lm(r-hat)  .
    !
    !     For each total angular momentum L, we pseudize qfunc(r)
    !     inside rc as:
    !
    !       qfunc(r) = r**(L+2) * [ a_1 + a_2*r**2 + a_3*r**4 ]
    !
    !     in such a way as to match qfunc and its 1'st derivative at
    !     rc, and to preserve
    !
    !       integral dr r**L * qfunc(r)   ,
    !
    !     i.e., to preserve the L'th moment of the charge.  The array
    !     qfunc has been set inside rc to correspond to this pseudized
    !     version using the minimal L, namely L = | l1-l2 | (e.g., L=0
    !     for diagonal elements).  The coefficients a_i (i=1,2,3)
    !     are stored in the array qfcoef(i,L+1,j,k) for each L so that
    !     the correctly pseudized versions of qfunc can be reconstructed
    !     for each L.  (Note that for given l1 and l2, only the values
    !     L = | l1-l2 | , | l1-l2 | +2 , ... , l1+l2 are ever used.)
    !     ------------------------------------------------------
    !
    USE constants, ONLY : fpi
    USE pseudo_types
    !
    implicit none
    !
    !    First the arguments passed to the subroutine
    !  
    TYPE (pseudo_upf) :: upf
    integer                                                           &
         &      is,        &! The number of the pseudopotential
         &      iunps       ! The unit of the pseudo file
    !
    !   Local variables

    real(DP)                                                     &
         &       exfact,        &! index of the exchange and correlation used 
         &       etotpseu,      &! total pseudopotential energy
         &       eloc,          &! energy of the local potential
         &       dummy,         &! dummy real variable
         &       rinner1,       &! rinner if only one is present
         &       rcloc           ! the cut-off radius of the local potential 
    real(DP), allocatable::  &
         &       ee(:),         &! the energy of the valence states
         &       rc(:),         &! the cut-off radii of the pseudopotential
         &       eee(:),        &! energies of the beta function
         &       ddd(:,:)        ! the screened D_{\mu,\nu} parameters
    integer, allocatable ::  &
         &       nnlz(:),       &! The nlm values of the valence states
         &       iptype(:)       ! more recent parameters 
    integer                                                           &
         &       iver(3),       &! contains the version of generating code
         &       idmy(3),       &! contains the date of creation of the pseudo
         &       ifpcor,        &! for core correction, 0 otherwise
         &       ios,           &! integer variable for I/O control
         &       i,             &! dummy counter 
         &       keyps,         &! the type of pseudopotential. Only US allowed
         &       irel,          &! says if the pseudopotential is relativistic
         &       ifqopt,        &! level of Q optimization
         &       npf,           &! as above
         &       nang,          &! number of angular momenta in pseudopotentials
         &       lloc,          &! angular momentum of the local part of PPs
         &       lp,            &! counter on Q angular momenta
         &       l,             &! counter on angular momenta
         &       iv, jv, ijv,   &! beta function counter
         &       ir              ! mesh points counter
    !
    character(len=20) :: title
    character(len=60) fmt
    !
    !     We first check the input variables
    !
    if (is <= 0) &
         call errore('readvan','routine called with wrong 1st argument', 1)
    if (iunps <= 0 .or. iunps >= 100000) &
         call errore('readvan','routine called with wrong 2nd argument', 1)
    !
    read(iunps, *, err=100, iostat=ios ) &
         (iver(i),i=1,3), (idmy(i),i=1,3)
    write(upf%generated, &
         "('Generated by Vanderbilt code, v. ',i1,'.',i1,'.',i1)") iver
    !
    if ( iver(1) > 7 .or. iver(1) < 1 .or. &
         iver(2) > 9 .or. iver(2) < 0 .or. &
         iver(3) > 9 .or. iver(3) < 0 ) & 
         call errore('readvan','wrong file version read',1)
    !
    read( iunps, '(a20,3f15.9)', err=100, iostat=ios ) &
         title, upf%zmesh, upf%zp, exfact 
    !
    upf%psd = title(1:2)
    !
    if ( upf%zmesh < 1 .or. upf%zmesh > 100.0_DP) &
         call errore( 'readvan','wrong zmesh read', is )
    if ( upf%zp <= 0.0_DP .or. upf%zp > 100.0_DP) &
         call errore('readvan','wrong atomic charge read', is )
    if ( exfact < -6 .or. exfact > 6) &
         &     call errore('readvan','Wrong xc in pseudopotential',1)
    ! convert from "our" conventions to Vanderbilt conventions
    call dftname_cp (nint(exfact), upf%dft)
    call set_dft_from_name( upf%dft )
    IF ( dft_is_meta() ) &
         CALL errore( 'readvan ', 'META-GGA not implemented', 1 )
    !
    read( iunps, '(2i5,1pe19.11)', err=100, iostat=ios ) &
         upf%nwfc, upf%mesh, etotpseu
    if ( upf%nwfc < 0 ) &
         call errore( 'readvan', 'wrong nchi read', upf%nwfc )
    if ( upf%mesh < 0 ) &
         call errore( 'readvan','wrong mesh', is )
    !
    !     info on pseudo eigenstates - energies are not used
    !
    ALLOCATE ( upf%oc(upf%nwfc), upf%lchi(upf%nwfc) ) 
    ALLOCATE ( nnlz(upf%nwfc), ee(upf%nwfc) )
    read( iunps, '(i5,2f15.9)', err=100, iostat=ios ) &
         ( nnlz(iv), upf%oc(iv), ee(iv), iv=1,upf%nwfc )
    do iv = 1, upf%nwfc
       i = nnlz(iv) / 100
       upf%lchi(iv) = nnlz(iv)/10 - i * 10
    enddo
    read( iunps, '(2i5,f15.9)', err=100, iostat=ios ) &
         keyps, ifpcor, rinner1
    upf%nlcc = (ifpcor == 1)
    !
    !     keyps= 0 --> standard hsc pseudopotential with exponent 4.0
    !            1 --> standard hsc pseudopotential with exponent 3.5
    !            2 --> vanderbilt modifications using defaults
    !            3 --> new generalized eigenvalue pseudopotentials
    !            4 --> frozen core all-electron case
    if ( keyps < 0 .or. keyps > 4 ) then
       call errore('readvan','wrong keyps',keyps)
    else if (keyps == 4) then
       call errore('readvan','keyps not implemented',keyps)
    end if
    upf%tvanp = (keyps == 3)
    upf%tpawp = .false.
    !
    !     Read information on the angular momenta, and on Q pseudization
    !     (version > 3.0)
    !
    if (iver(1) >= 3) then
       read( iunps, '(2i5,f9.5,2i5,f9.5)', err=100, iostat=ios )  &
            nang, lloc, eloc, ifqopt, upf%nqf, dummy
!!! PWSCF: lmax(is)=nang, lloc(is)=lloc
       !
       !    NB: In the Vanderbilt atomic code the angular momentum goes 
       !        from 1 to nang
       !
       if ( nang < 0 ) &
            call errore(' readvan', 'Wrong nang read', nang)
       if ( lloc == -1 ) lloc = nang+1
       if ( lloc > nang+1 .or. lloc < 0 ) &
            call errore( 'readvan', 'wrong lloc read', is )
       if ( upf%nqf < 0 ) &
            call errore(' readvan', 'Wrong nqf read', upf%nqf)
       if ( ifqopt < 0 ) &
            call errore( 'readvan', 'wrong ifqopt read', is )
    else
       ! old format: no distinction between nang and nchi
       nang = upf%nwfc
    end if
    !
    !     Read and test the values of rinner (version > 5.1)
    !     rinner = radius at which to cut off partial core or q_ij
    !
    ALLOCATE ( upf%rinner(2*nang-1) ) 
    if (10*iver(1)+iver(2) >= 51) then
       !
       read( iunps, *, err=100, iostat=ios ) &
            (upf%rinner(lp), lp=1,2*nang-1 )
       !
       do lp = 1, 2*nang-1
          if (upf%rinner(lp) < 0.0_DP) &
               call errore('readvan','Wrong rinner read', is )
       enddo
    else if (iver(1) > 3) then
       do lp = 2, 2*nang-1
          upf%rinner(lp)=rinner1
       end do
    end if
    !
    if (iver(1) >= 4) &
         read( iunps, '(i5)',err=100, iostat=ios ) irel
    !     
    !       set the number of angular momentum terms in q_ij to read in
    !
    if (iver(1) == 1) then
       oldvan(is) = .TRUE.
       ! old format: no optimization of q_ij => 3-term taylor series
       upf%nqf=3
       upf%nqlc=5
    else if (iver(1) == 2) then
       upf%nqf=3
       upf%nqlc = 2*nang - 1
    else
       upf%nqlc = 2*nang - 1
    end if
    !
    if ( upf%nqlc > lqmax .or. upf%nqlc < 0 ) &
         call errore(' readvan', 'Wrong  nqlc read', upf%nqlc )
    !
    ALLOCATE ( rc(nang) )
    read( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
         ( rc(l), l=1,nang )
    !
    !     reads the number of beta functions 
    !
    read( iunps, '(2i5)', err=100, iostat=ios ) &
         upf%nbeta, upf%kkbeta
    !
    ALLOCATE ( upf%kbeta(upf%nbeta) )
    upf%kbeta(:) = upf%kkbeta
    !
    if( upf%nbeta < 0 ) &
         call errore( 'readvan','nbeta wrong', is )
    if( upf%kkbeta > upf%mesh .or. upf%kkbeta < 0 ) &
         call errore( 'readvan','kkbeta wrong or too large', is )
    !
    !    Now reads the main Vanderbilt parameters
    !
    ALLOCATE ( upf%lll(upf%nbeta) )
    ALLOCATE ( upf%beta(upf%mesh,upf%nbeta) )
    ALLOCATE ( upf%dion(upf%nbeta,upf%nbeta), upf%qqq(upf%nbeta,upf%nbeta) )
    ALLOCATE ( upf%qfunc(upf%mesh,upf%nbeta*(upf%nbeta+1)/2) )
    ALLOCATE ( upf%qfcoef(upf%nqf, upf%nqlc, upf%nbeta, upf%nbeta) )
    ALLOCATE ( eee(upf%nbeta), ddd(upf%nbeta,upf%nbeta) )
    do iv=1,upf%nbeta
       read( iunps, '(i5)',err=100, iostat=ios ) upf%lll(iv)
       read( iunps, '(1p4e19.11)',err=100, iostat=ios ) &
            eee(iv), ( upf%beta(ir,iv), ir=1,upf%kkbeta )
       do ir=upf%kkbeta+1,upf%mesh
          upf%beta(ir,iv)=0.0_DP
       enddo
       if ( upf%lll(iv) > lmaxx .or. upf%lll(iv) < 0 ) &
            call errore( 'readvan', 'lll wrong or too large ', is )
       do jv=iv,upf%nbeta
          !
          !  the symmetric matric Q_{nb,mb} is stored in packed form
          !  Q(iv,jv) => qfunc(ijv) as defined below (for jv >= iv)
          !
          ijv = jv * (jv-1) / 2 + iv
          read( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
               upf%dion(iv,jv), ddd(iv,jv), upf%qqq(iv,jv), &
               (upf%qfunc(ir,ijv),ir=1,upf%kkbeta),         &
               ((upf%qfcoef(i,lp,iv,jv),i=1,upf%nqf),lp=1,upf%nqlc)
          do ir=upf%kkbeta+1,upf%mesh
            upf%qfunc(ir,ijv)=0.0_DP
          enddo
          !
          !     Use the symmetry of the coefficients
          !
          if ( iv /= jv ) then
             upf%dion(jv,iv)=upf%dion(iv,jv)
             upf%qqq(jv,iv) =upf%qqq(iv,jv)
             upf%qfcoef(:,:,jv,iv)=upf%qfcoef(:,:,iv,jv)
          end if
       enddo
    enddo
    !
    ! Set additional, not present, variables to dummy values
    ALLOCATE(upf%els(upf%nwfc))
    upf%els(:) = 'nX'
    ALLOCATE(upf%els_beta(upf%nbeta))
    upf%els_beta(:) = 'nX'
    ALLOCATE(upf%rcut(upf%nbeta), upf%rcutus(upf%nbeta))
    upf%rcut(:) = 0._dp
    upf%rcutus(:) = 0._dp

    DEALLOCATE (ddd)
    !
    !    for versions later than 7.2
    !
    if (10*iver(1)+iver(2) >= 72) then
       ALLOCATE (iptype(upf%nbeta))
       read( iunps, '(6i5)',err=100, iostat=ios ) &
            (iptype(iv), iv=1,upf%nbeta)
       read( iunps, '(i5,f15.9)',err=100, iostat=ios ) &
            npf, dummy
       DEALLOCATE (iptype)
    end if
    !
    !   read the local potential
    !
    ALLOCATE ( upf%vloc(upf%mesh) )
    read( iunps, '(1p4e19.11)',err=100, iostat=ios ) &
         rcloc, ( upf%vloc(ir), ir=1,upf%mesh )
    !
    !   If present reads the core charge rho_atc(r)=4*pi*r**2*rho_core(r)
    !
    if ( upf%nlcc ) then 
       ALLOCATE ( upf%rho_atc(upf%mesh) )
       if (iver(1) >= 7) &
            read( iunps, '(1p4e19.11)', err=100, iostat=ios ) dummy
       read( iunps, '(1p4e19.11)', err=100, iostat=ios )  &
            ( upf%rho_atc(ir), ir=1,upf%mesh )
    endif
    !
    !     Read the screened local potential (not used)
    !
    ALLOCATE ( upf%rho_at(upf%mesh) )
    read( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
         (upf%rho_at(ir), ir=1,upf%mesh)
    !
    !     Read the valence atomic charge
    !
    read( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
         (upf%rho_at(ir), ir=1,upf%mesh)
    !
    !     Read the logarithmic mesh (if version > 1)
    !
    ALLOCATE ( upf%r(upf%mesh), upf%rab(upf%mesh) ) 
    if (iver(1) >1) then
       read( iunps, '(1p4e19.11)',err=100, iostat=ios ) &
            (upf%r(ir),ir=1,upf%mesh)
       read( iunps, '(1p4e19.11)',err=100, iostat=ios ) &
            (upf%rab(ir),ir=1,upf%mesh)
    else
       !
       !     generate herman-skillman mesh (if version = 1)
       !
       call herman_skillman_grid &
          ( upf%mesh, upf%zmesh, upf%r, upf%rab )
    end if
    !
    !     convert vloc to the conventions used in the rest of the code
    !     (as read from Vanderbilt's format it is r*v_loc(r))
    !
    do ir = 2, upf%mesh
       upf%vloc (ir) = upf%vloc (ir) / upf%r(ir)
    enddo
    upf%vloc (1) = upf%vloc (2)
    !
    !     set rho_atc(r)=rho_core(r)  (without 4*pi*r^2 factor,
    !     for compatibility with rho_atc in the non-US case)
    !
    if (upf%nlcc) then
       upf%rho_atc(1) = 0.0_DP
       do ir=2,upf%mesh
          upf%rho_atc(ir) = upf%rho_atc(ir)/fpi/upf%r(ir)**2
       enddo
    end if
    !
    !    Read the wavefunctions of the atom
    !
    if (iver(1) >= 7) then
       read( iunps, *, err=100, iostat=ios ) i
       if (i /= upf%nwfc) &
            call errore('readvan','unexpected or unimplemented case',1)
    end if
    !
    ALLOCATE ( upf%chi(upf%mesh, upf%nwfc) )
    if (iver(1) >= 6) &
         read( iunps, *, err=100, iostat=ios ) &
         ( (upf%chi(ir,iv), ir=1,upf%mesh), iv=1,upf%nwfc )
    !
    if (iver(1) == 1) then
       !
       !   old version: read the q_l(r) and fit them with the Vanderbilt's form
       ! 
       call fit_qrl ( )
       !
    end if
    !
    !    Here we write on output information on the pseudopotential 
    !
    WRITE( stdout,200) is
200 format (/4x,60('=')/4x,'|  pseudopotential report',               &
         &        ' for atomic species:',i3,11x,'|')
    WRITE( stdout,300) 'pseudo potential version', &
         iver(1), iver(2), iver(3)
300 format (4x,'|  ',1a30,3i4,13x,' |' /4x,60('-'))
    WRITE( stdout,400) title, upf%dft
400 format (4x,'|  ',2a20,' exchange-corr  |')
    WRITE( stdout,500) upf%zmesh, is, upf%zp, exfact
500 format (4x,'|  z =',f5.0,4x,'zv(',i2,') =',f5.0,4x,'exfact =',    &
         &     f10.5, 9x,'|')
    WRITE( stdout,600) ifpcor, etotpseu
600 format (4x,'|  ifpcor = ',i2,10x,' atomic energy =',f10.5,        &
         &     ' Ry',6x,'|')
    WRITE( stdout,700)
700 format(4x,'|  index    orbital      occupation    energy',14x,'|')
    WRITE( stdout,800) ( iv, nnlz(iv), upf%oc(iv), ee(iv), iv=1,upf%nwfc )
    DEALLOCATE (ee, nnlz)
800 format(4x,'|',i5,i11,5x,f10.2,f12.2,15x,'|')
    if (iver(1) >= 3 .and. nang > 0) then
       IF (nang < 4) THEN
          write(fmt,900) 2*nang-1, 40-8*(2*nang-2)
       ELSE
          write(fmt,900) 2*nang-1, 1
       ENDIF
900    format('(4x,"|  rinner =",',i1,'f8.4,',i2,'x,"|")')
       WRITE( stdout,fmt)  (upf%rinner(lp),lp=1,2*nang-1)
    end if
    WRITE( stdout,1000)
1000 format(4x,'|    new generation scheme:',32x,'|')
    WRITE( stdout,1100) upf%nbeta, upf%kkbeta, rcloc
1100 format(4x,'|    nbeta = ',i2,5x,'kkbeta =',i5,5x,'rcloc =',f10.4,4x,&
         &     '|'/4x,'|    ibeta    l     epsilon   rcut',25x,'|')
    do iv = 1, upf%nbeta
       lp=upf%lll(iv)+1
       WRITE( stdout,1200) iv,upf%lll(iv),eee(iv),rc(lp)
1200   format(4x,'|',5x,i2,6x,i2,4x,2f7.2,25x,'|')
    enddo
    WRITE( stdout,1300)
1300 format (4x,60('='))
    !
    DEALLOCATE (eee, rc)
    return
100 call errore('readvan','error reading pseudo file', abs(ios) )
  !
  CONTAINS
  !-----------------------------------------------------------------------
  subroutine fit_qrl ( )
    !-----------------------------------------------------------------------
    !
    ! find coefficients qfcoef that fit the pseudized qrl in US PP
    ! these coefficients are written to file in newer versions of the 
    ! Vanderbilt PP generation code but not in some ancient versions
    !
    implicit none
    !
    real (kind=DP), allocatable :: qrl(:,:), a(:,:), ainv(:,:), b(:), x(:)
    integer :: iv, jv, ijv, lmin, lmax, l, ir, irinner, i,j
    !
    !
    allocate ( a(upf%nqf,upf%nqf), ainv(upf%nqf,upf%nqf) )
    allocate ( b(upf%nqf), x(upf%nqf) )
    ALLOCATE ( qrl(upf%kkbeta, upf%nqlc) )
    !
    do iv=1,upf%nbeta
       do jv=iv,upf%nbeta
          !
          ! original version, assuming lll(jv) >= lll(iv) 
          !   lmin=lll(jv,is)-lll(iv,is)+1
          !   lmax=lmin+2*lll(iv,is)
          ! note that indices run from 1 to Lmax+1, not from 0 to Lmax
          !
          lmin = ABS( upf%lll(jv) - upf%lll(iv) ) + 1
          lmax =      upf%lll(jv) + upf%lll(iv)   + 1
          IF ( lmin < 1 .OR. lmax >  SIZE(qrl,2)) &
               CALL errore ('fit_qrl', 'bad 2rd dimension for array qrl', 1)
          !
          !  read q_l(r) for all l
          !
          read(iunps,*, err=100) &
                  ( (qrl(ir,l),ir=1,upf%kkbeta), l=lmin,lmax)
          !
          ijv = jv * (jv-1) / 2 + iv
          !
          do l=lmin,lmax
             !
             ! reconstruct rinner
             !
             do ir=upf%kkbeta,1,-1
                if ( abs(qrl(ir,l)-upf%qfunc(ir,ijv)) > 1.0d-6) go to 10
             end do
10           irinner = ir+1
             upf%rinner(l) = upf%r(irinner)
             !
             ! least square minimization: find
             ! qrl = sum_i c_i r^{l+1}r^{2i-2} for r < rinner
             !
             a(:,:) = 0.0_DP
             b(:)   = 0.0_DP
             do i = 1, upf%nqf
                do ir=1,irinner
                   b(i) = b(i) + upf%r(ir)**(2*i-2+l+1) * qrl(ir,l)
                end do
                do j = i, upf%nqf
                   do ir=1,irinner
                      a(i,j) = a(i,j) + upf%r(ir)**(2*i-2+l+1) * &
                                        upf%r(ir)**(2*j-2+l+1) 
                   end do
                   if (j > i) a(j,i) = a(i,j) 
                end do
             end do
             !
             call invmat (upf%nqf, a, ainv)
             !
             do i = 1, upf%nqf
                upf%qfcoef(i,l,iv,jv) = dot_product(ainv(i,:),b(:))
                if (iv /= jv) upf%qfcoef(i,l,jv,iv) = upf%qfcoef(i,l,iv,jv)
             end do
          end do
       end do
    end do
    !
    deallocate ( qrl, x, b , ainv, a )
    return
    !
100 call errore('readvan','error reading Q_L(r)', 1 )
  end subroutine fit_qrl
  !
  end subroutine readvan
  !-----------------------------------------------------------------------
  SUBROUTINE herman_skillman_grid (mesh,z,r,rab)
    !-----------------------------------------------------------------------
    !
    !     generate Herman-Skillman radial grid (obsolescent)
    !     c    - 0.88534138/z**(1/3)
    !
    IMPLICIT NONE
    !
    INTEGER mesh
    REAL(DP) :: z, r(mesh), rab(mesh)
    !
    REAL(DP) :: deltax,pi
    INTEGER :: nblock,i,j,k
    !
    pi=4.0_DP*ATAN(1.0_DP)
    nblock = mesh/40
    i=1
    r(i)=0.0_DP
    deltax=0.0025_DP*0.5_DP*(3.0_DP*pi/4.0_DP)**(2.0_DP/3.0_DP)/z**(1.0_DP/3.0_DP)
    DO j=1,nblock
       DO k=1,40
          i=i+1
          r(i)=r(i-1)+deltax
          rab(i)=deltax
       END DO
       deltax=deltax+deltax
    END DO
    !
    RETURN
  END SUBROUTINE herman_skillman_grid
  !
  !---------------------------------------------------------------------
  subroutine readrrkj( iunps, is, upf )
    !---------------------------------------------------------------------
    !
    !     This routine reads Vanderbilt pseudopotentials produced by the
    !     code of Andrea Dal Corso. Hard PPs are first generated
    !     according to the Rabe Rappe Kaxiras Johannopoulos recipe.
    !     Ultrasoft PP's are subsequently generated from the hard PP's.
    !
    !     Output parameters in module "uspp_param"
    !     info on DFT level in module "dft"
    !
    USE constants, ONLY : fpi
    USE pseudo_types
    !
    implicit none
    !
    !    First the arguments passed to the subroutine
    !
    TYPE (pseudo_upf) :: upf
    integer :: &
         is,        &! The index of the pseudopotential
         iunps       ! the unit from with pseudopotential is read
    !
    !    Local variables
    !
    integer::  iexch, icorr, igcx, igcc

    integer:: &
         nb,mb, ijv,&! counters on beta functions
         n,         &! counter on mesh points
         ir,        &! counters on mesh points
         pseudotype,&! the type of pseudopotential
         ios,       &! I/O control
         ndum,      &! dummy integer variable
         l           ! counter on angular momentum
    real(DP):: &
         x,         &! auxiliary variable
         etotps,    &! total energy of the pseudoatom
         rdum        ! dummy real variable
    !
    logical :: & 
         rel         ! if true the atomic calculation is relativistic
    !
    character(len=75) :: &
         titleps     ! the title of the pseudo
    !
    integer :: &
         lmax       ! max angular momentum
    character(len=2) :: &
         adum       ! dummy character variable
    !
    !     We first check the input variables
    !
    if (is <= 0) &
         call errore('readrrkj','routine called with wrong 1st argument', 1)
    if (iunps <= 0 .or. iunps >= 100000) &
         call errore('readrrkj','routine called with wrong 2nd argument', 1)
    !
    read( iunps, '(a75)', err=100, iostat=ios ) &
         titleps
    upf%psd = titleps(7:8)
    !
    read( iunps, '(i5)',err=100, iostat=ios ) &
         pseudotype
    upf%tvanp = (pseudotype == 3)
    upf%tpawp = .false.

    if ( upf%tvanp ) then
       upf%generated = &
         "RRKJ3 Ultrasoft PP, generated by Andrea Dal Corso code"
    else
       upf%generated = &
         "RRKJ3 norm-conserving PP, generated by Andrea Dal Corso code"
    endif

    read( iunps, '(2l5)',err=100, iostat=ios ) &
         rel, upf%nlcc
    read( iunps, '(4i5)',err=100, iostat=ios ) &
         iexch, icorr, igcx,  igcc
    !
    ! workaround to keep track of which dft was read
    ! See also upf2internals
    !
    write( upf%dft, "('INDEX:',4i1)") iexch,icorr,igcx,igcc
    call set_dft_from_indices(iexch,icorr,igcx,igcc, 0) ! Cannot read nonlocal in this format

    read( iunps, '(2e17.11,i5)') &
         upf%zp, etotps, lmax
    if ( upf%zp < 1 .or. upf%zp > 100 ) &
         call errore('readrrkj','wrong potential read',is)
    !
    read( iunps, '(4e17.11,i5)',err=100, iostat=ios ) &
         upf%xmin, rdum, upf%zmesh, upf%dx, upf%mesh
    !
    if ( upf%mesh < 0) &
         call errore('readrrkj', 'wrong mesh',is)
    !
    read( iunps, '(2i5)', err=100, iostat=ios ) &
         upf%nwfc, upf%nbeta
    !
    if ( upf%nbeta < 0) &
         call errore('readrrkj', 'wrong nbeta', is)
    if ( upf%nwfc < 0 ) &
         call errore('readrrkj', 'wrong nchi', is)
    !
    read( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
         ( rdum, nb=1,upf%nwfc )
    read( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
         ( rdum, nb=1,upf%nwfc )
    !
    ALLOCATE ( upf%oc(upf%nwfc), upf%lchi(upf%nwfc), upf%lll(upf%nwfc) ) 
    !
    do nb=1,upf%nwfc
       read(iunps,'(a2,2i3,f6.2)',err=100,iostat=ios) &
            adum, ndum, upf%lchi(nb), upf%oc(nb)
       upf%lll(nb)=upf%lchi(nb)
       !
       ! oc < 0 distinguishes between bound states from unbound states
       !
       if ( upf%oc(nb) <= 0.0_DP) upf%oc(nb) = -1.0_DP
    enddo
    !
    ALLOCATE ( upf%kbeta(upf%nbeta) )
    ALLOCATE ( upf%dion(upf%nbeta,upf%nbeta), upf%qqq(upf%nbeta,upf%nbeta) )
    ALLOCATE ( upf%beta(upf%mesh,upf%nbeta) )
    ALLOCATE ( upf%qfunc(upf%mesh,upf%nbeta*(upf%nbeta+1)/2) )
    upf%kkbeta = 0
    do nb=1,upf%nbeta
       read ( iunps, '(i6)',err=100, iostat=ios ) upf%kbeta(nb)
       upf%kkbeta = MAX ( upf%kkbeta, upf%kbeta(nb) )
       read ( iunps, '(1p4e19.11)',err=100, iostat=ios ) &
            ( upf%beta(ir,nb), ir=1,upf%kbeta(nb))
       do ir=upf%kbeta(nb)+1,upf%mesh
          upf%beta(ir,nb)=0.0_DP
       enddo
       do mb=1,nb
          ! 
          ! the symmetric matric Q_{nb,mb} is stored in packed form
          ! Q(nb,mb) => qfunc(ijv) as defined below (for mb <= nb)
          ! 
          ijv = nb * (nb - 1) / 2 + mb
          read( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
               upf%dion(nb,mb)
          if (pseudotype == 3) then
             read(iunps,'(1p4e19.11)',err=100,iostat=ios) &
                  upf%qqq(nb,mb)
             read(iunps,'(1p4e19.11)',err=100,iostat=ios) &
                  (upf%qfunc(n,ijv),n=1,upf%mesh)
          else
             upf%qqq(nb,mb)=0.0_DP
             upf%qfunc(:,ijv)=0.0_DP
          endif
          if ( mb /= nb ) then
             upf%dion(mb,nb)=upf%dion(nb,mb)
             upf%qqq(mb,nb)=upf%qqq(nb,mb)
          end if
       enddo
    enddo
    !
    !   reads the local potential 
    !
    ALLOCATE ( upf%vloc(upf%mesh) )
    read( iunps, '(1p4e19.11)',err=100, iostat=ios ) &
         rdum, ( upf%vloc(ir), ir=1,upf%mesh )
    !
    !   reads the atomic charge
    !
    ALLOCATE ( upf%rho_at(upf%mesh) )
    read( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
         ( upf%rho_at(ir), ir=1,upf%mesh )
    !
    !   if present reads the core charge
    !
    if ( upf%nlcc ) then 
       ALLOCATE ( upf%rho_atc(upf%mesh) )
       read( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
            ( upf%rho_atc(ir), ir=1,upf%mesh )
    endif
    !
    !   read the pseudo wavefunctions of the atom
    !  
    ALLOCATE ( upf%chi(upf%mesh, upf%nwfc) )
    read( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
         ((upf%chi(ir,nb),ir=1,upf%mesh),nb=1,upf%nwfc)
    !
    !    set several variables for compatibility with the rest of the code
    !
    upf%nqf=0
    upf%nqlc=2*lmax+1
    if ( upf%nqlc > lqmax .or. upf%nqlc < 0 ) &
         call errore(' readrrkj', 'Wrong  nqlc', upf%nqlc )
    ALLOCATE ( upf%rinner(upf%nqlc) )
    do l=1,upf%nqlc
       upf%rinner(l)=0.0_DP
    enddo
    !
    !    compute the radial mesh
    !
    ALLOCATE ( upf%r(upf%mesh), upf%rab(upf%mesh) )
    do ir = 1, upf%mesh
       x = upf%xmin + DBLE(ir-1) * upf%dx
       upf%r(ir) = EXP(x) / upf%zmesh
       upf%rab(ir) = upf%dx * upf%r(ir)
    end do
    !
    !     set rho_atc(r)=rho_core(r)  (without 4*pi*r^2 factor)
    !
    if ( upf%nlcc ) then
       do ir=1,upf%mesh
          upf%rho_atc(ir) = upf%rho_atc(ir)/fpi/upf%r(ir)**2
       enddo
    end if
    !
    ! Set additional, not present, variables to dummy values
    allocate(upf%els(upf%nwfc))
    upf%els(:) = 'nX'
    allocate(upf%els_beta(upf%nbeta))
    upf%els_beta(:) = 'nX'
    allocate(upf%rcut(upf%nbeta), upf%rcutus(upf%nbeta))
    upf%rcut(:) = 0._dp
    upf%rcutus(:) = 0._dp
    !
    return
100 call errore('readrrkj','Reading pseudo file',abs(ios))
  end subroutine readrrkj
  !
end module read_uspp_module
