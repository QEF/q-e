
! Copyright (C) 2006 Quantum-Espresso group
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
  USE parameters, ONLY: nchix, lmaxx, nbrx, ndmx, nsx, lqmax, nqfx
  USE io_global, ONLY: stdout, meta_ionode
  USE funct, ONLY: set_dft_from_name, dft_is_hybrid, dft_is_meta, &
       set_dft_from_indices
  !
  ! Variables above are not modified, variables below are
  !
  USE uspp_param, ONLY: qfunc, qfcoef, qqq, betar, dion, vloc_at, &
       rinner, kkbeta, lll, nbeta, nqf, nqlc, tvanp, oldvan, iver, psd
  USE atom, ONLY: zmesh, nchi, chi, lchi, r, rab, mesh, nlcc, oc, &
        rho_at, rho_atc, xmin, dx
  USE ions_base, ONLY: zv
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  PUBLIC :: readvan, readrrkj
  !
CONTAINS
  !---------------------------------------------------------------------
  subroutine readvan( is, iunps )
    !---------------------------------------------------------------------
    !
    !     Read Vanderbilt pseudopotential for species "is" from unit "iunps"
    !     Output parameters in module "uspp_param"
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
    !
    implicit none
    !
    !    First the arguments passed to the subroutine
    !
    integer                                                           &
         &      is,        &! The number of the pseudopotential
         &      iunps       ! The unit of the pseudo file
    !
    !   Local variables

    real(DP)                                                     &
         &       exfact,        &! index of the exchange and correlation used 
         &       etotpseu,      &! total pseudopotential energy
         &       ee(nchix),     &! the energy of the valence states
         &       eloc,          &! energy of the local potential
         &       dummy,         &! dummy real variable
         &       rc(nchix),     &! the cut-off radii of the pseudopotential
         &       eee(nbrx),     &! energies of the beta function
         &       ddd(nbrx,nbrx),&! the screened D_{\mu,\nu} parameters
         &       rcloc           ! the cut-off radius of the local potential 
    integer                                                           &
         &       idmy(3),       &! contains the date of creation of the pseudo
         &       ifpcor,        &! for core correction, 0 otherwise
         &       ios,           &! integer variable for I/O control
         &       i,             &! dummy counter 
         &       nnlz(nchix),   &! The nlm values of the valence states
         &       keyps,         &! the type of pseudopotential. Only US allowed
         &       irel,          &! it says if the pseudopotential is relativistic
         &       ifqopt(nsx),   &! level of Q optimization
         &       iptype(nbrx),  &! more recent parameters 
         &       npf,           &! as above
         &       nang,          &! number of angular momenta in pseudopotentials
         &       lloc,          &! angular momentum of the local part of PPs
         &       lmin, lmax,    &! min and max angular momentum in Q
         &       lp,            &! counter on Q angular momenta
         &       l,             &! counter on angular momenta
         &       iv, jv, ijv,   &! beta function counter
         &       ir              ! mesh points counter
    !
    character(len=20) :: title, dft_name
    character(len=60) fmt
    !
    !     We first check the input variables
    !
    if (is <= 0 .or. is >= nsx) &
         call errore('readvan','routine called with wrong 1st argument', 1)
    if (iunps <= 0 .or. iunps >= 100000) &
         call errore('readvan','routine called with wrong 2nd argument', 1)
    !
    read(iunps, *, err=100 ) &
         (iver(i,is),i=1,3), (idmy(i),i=1,3)
    !
    if ( iver(1,is) > 7 .or. iver(1,is) < 1 .or. &
         iver(2,is) > 9 .or. iver(2,is) < 0 .or. &
         iver(3,is) > 9 .or. iver(3,is) < 0 ) & 
         call errore('readvan','wrong file version read',1)
    !
    read( iunps, '(a20,3f15.9)', err=100, iostat=ios ) &
         title, zmesh(is), zv(is), exfact 
    !
    psd (is) = title(1:2)
    !
    if ( zmesh(is) < 1 .or. zmesh(is) > 100.d0) &
         call errore( 'readvan','wrong zmesh read', is )
    if ( zv(is) <= 0.d0 .or. zv(is) > 100.d0) &
         call errore('readvan','wrong zv read', is )
    if ( exfact < -6 .or. exfact > 6) &
         &     call errore('readvan','Wrong xc in pseudopotential',1)
    ! convert from "our" conventions to Vanderbilt conventions
    call dftname_cp (nint(exfact), dft_name)
    call set_dft_from_name( dft_name )
#if !defined (EXX)
    IF ( dft_is_hybrid() ) &
         CALL errore( 'readvan', 'HYBRID XC not implemented', 1 )
#endif
    IF ( dft_is_meta() ) &
         CALL errore( 'readvan ', 'META-GGA not implemented', 1 )
    !
    read( iunps, '(2i5,1pe19.11)', err=100, iostat=ios ) &
         nchi(is), mesh(is), etotpseu
    if ( nchi(is) < 0 .OR. nchi(is) > nchix ) &
         call errore( 'readvan', 'wrong or too large nchi read', nchi(is) )
    if ( mesh(is) > ndmx .or. mesh(is) < 0 ) &
         call errore( 'readvan','wrong mesh', is )
    !
    !     info on pseudo eigenstates - energies are not used
    !
    read( iunps, '(i5,2f15.9)', err=100, iostat=ios ) &
         ( nnlz(iv),  oc(iv,is), ee(iv), iv=1,nchi(is) )
    do iv = 1, nchi(is)
       i = nnlz(iv) / 100
       lchi(iv,is) = nnlz(iv)/10 - i * 10
    enddo
    read( iunps, '(2i5,f15.9)', err=100, iostat=ios ) &
         keyps, ifpcor, rinner(1,is)
    nlcc (is) = (ifpcor == 1)
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
    tvanp(is) = (keyps == 3)
    !
    !     Read information on the angular momenta, and on Q pseudization
    !     (version > 3.0)
    !
    if (iver(1,is) >= 3) then
       read( iunps, '(2i5,f9.5,2i5,f9.5)', err=100, iostat=ios )  &
            nang, lloc, eloc, ifqopt(is), nqf(is), dummy
!!! PWSCF: lmax(is)=nang, lloc(is)=lloc
       !
       !    NB: In the Vanderbilt atomic code the angular momentum goes 
       !        from 1 to nang
       !
       if ( nang > nchix + 1 .or. nang < 0 ) &
            call errore(' readvan', 'Wrong nang read', nang)
       if ( lloc == -1 ) lloc = nang+1
       if ( lloc > nang+1 .or. lloc < 0 ) &
            call errore( 'readvan', 'wrong lloc read', is )
       if ( nqf(is) > nqfx .or. nqf(is) < 0 ) &
            call errore(' readvan', 'Wrong nqf read', nqf(is))
       if ( ifqopt(is) < 0 ) &
            call errore( 'readvan', 'wrong ifqopt read', is )
    end if
    !
    !     Read and test the values of rinner (version > 5.1)
    !     rinner = radius at which to cut off partial core or q_ij
    !
    if (10*iver(1,is)+iver(2,is) >= 51) then
       !
       read( iunps, *, err=100, iostat=ios ) &
            (rinner(lp,is), lp=1,2*nang-1 )
       !
       do lp = 1, 2*nang-1
          if (rinner(lp,is) < 0.d0) &
               call errore('readvan','Wrong rinner read', is )
       enddo
    else if (iver(1,is) > 3) then
       do lp = 2, 2*nang-1
          rinner(lp,is)=rinner(1,is)
       end do
    end if
    !
    if (iver(1,is) >= 4) &
         read( iunps, '(i5)',err=100, iostat=ios ) irel
    !     
    !       set the number of angular momentum terms in q_ij to read in
    !
    if (iver(1,is) == 1) then
       oldvan(is) = .TRUE.
       ! old format: no distinction between nang and nchi
       nang = nchi(is)
       ! old format: no optimization of q_ij => 3-term taylor series
       nqf(is)=3
       nqlc(is)=5
    else if (iver(1,is) == 2) then
       nang = nchi(is)
       nqf(is)=3
       nqlc(is) = 2*nang - 1
    else
       nqlc(is) = 2*nang - 1
    end if
    !
    if ( nqlc(is) > lqmax .or. nqlc(is) < 0 ) &
         call errore(' readvan', 'Wrong  nqlc read', nqlc(is) )
    !
    read( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
         ( rc(l), l=1,nang )
    !
    !     reads the number of beta functions 
    !
    read( iunps, '(2i5)', err=100, iostat=ios ) &
         nbeta(is), kkbeta(is)
    !
    if( nbeta(is) > nbrx .or. nbeta(is) <0 ) &
         call errore( 'readvan','nbeta wrong or too large', is )
    if( kkbeta(is) > mesh(is) .or. kkbeta(is) < 0 ) &
         call errore( 'readvan','kkbeta wrong or too large', is )
    !
    !    Now reads the main Vanderbilt parameters
    !
    do iv=1,nbeta(is)
       read( iunps, '(i5)',err=100, iostat=ios ) &
            lll(iv,is)
       read( iunps, '(1p4e19.11)',err=100, iostat=ios ) &
            eee(iv), ( betar(ir,iv,is), ir=1,kkbeta(is) )
       if ( lll(iv,is) > lmaxx .or. lll(iv,is) < 0 ) &
            call errore( 'readvan', 'lll wrong or too large ', is )
       do jv=iv,nbeta(is)
          read( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
               dion(iv,jv,is), ddd(iv,jv), qqq(iv,jv,is),   &
               (qfunc(ir,iv,jv,is),ir=1,kkbeta(is)),        &
               ((qfcoef(i,lp,iv,jv,is),i=1,nqf(is)),lp=1,nqlc(is))
          !
          !     Use the symmetry of the coefficients
          !
          dion(jv,iv,is)=dion(iv,jv,is)
          qqq(jv,iv,is)=qqq(iv,jv,is)
          !
          do ir = 1, kkbeta(is)
             qfunc(ir,jv,iv,is)=qfunc(ir,iv,jv,is)
          enddo
          !
          do i = 1, nqf(is)
             do lp= 1, nqlc(is)
                qfcoef(i,lp,jv,iv,is)=qfcoef(i,lp,iv,jv,is)
             enddo
          enddo
       enddo
    enddo
    !
    !    for versions later than 7.2
    !
    if (10*iver(1,is)+iver(2,is) >= 72) then
       read( iunps, '(6i5)',err=100, iostat=ios ) &
            (iptype(iv), iv=1,nbeta(is))
       read( iunps, '(i5,f15.9)',err=100, iostat=ios ) &
            npf, dummy
    end if

    !
    !   read the local potential vloc_at
    !
    read( iunps, '(1p4e19.11)',err=100, iostat=ios ) &
         rcloc, ( vloc_at(ir,is), ir=1,mesh(is) )
    !
    !   If present reads the core charge rho_atc(r)=4*pi*r**2*rho_core(r)
    !
    if ( nlcc(is) ) then 
       if (iver(1,is) >= 7) &
            read( iunps, '(1p4e19.11)', err=100, iostat=ios ) dummy
       read( iunps, '(1p4e19.11)', err=100, iostat=ios )  &
            ( rho_atc(ir,is), ir=1,mesh(is) )
    endif
    !
    !     Read the screened local potential (not used)
    !
    read( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
         (rho_at(ir,is), ir=1,mesh(is))
    !
    !     Read the valence atomic charge
    !
    read( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
         (rho_at(ir,is), ir=1,mesh(is))
    !
    !     Read the logarithmic mesh (if version > 1)
    !
    if (iver(1,is) >1) then
       read( iunps, '(1p4e19.11)',err=100, iostat=ios ) &
            (r(ir,is),ir=1,mesh(is))
       read( iunps, '(1p4e19.11)',err=100, iostat=ios ) &
            (rab(ir,is),ir=1,mesh(is))
    else
       !
       !     generate herman-skillman mesh (if version = 1)
       !
       call herman_skillman_grid (mesh(is),zmesh(is),r(1,is),rab(1,is))
    end if
    !
    !     convert vloc_at to the conventions used in the rest of the code
    !     (as read from Vanderbilt's format it is r*v_loc(r))
    !
    do ir = 2, mesh (is)
       vloc_at (ir, is) = vloc_at (ir, is) / r (ir, is)
    enddo
    vloc_at (1, is) = vloc_at (2, is)
    !
    !     set rho_atc(r)=rho_core(r)  (without 4*pi*r^2 factor,
    !     for compatibility with rho_atc in the non-US case)
    !
    if (nlcc (is) ) then
       rho_atc(1,is) = 0.D0
       do ir=2,mesh(is)
          rho_atc(ir,is) = rho_atc(ir,is)/4.0/3.14159265/r(ir,is)**2
       enddo
    end if
    !
    !    Read the wavefunctions of the atom
    !      
    if (iver(1,is) >= 7) then
       read( iunps, *, err=100, iostat=ios ) i
       if (i /= nchi(is)) &
            call errore('readvan','unexpected or unimplemented case',1)
    end if
    !
    if (iver(1,is) >= 6) &
         read( iunps, *, err=100, iostat=ios ) &
         ((chi(ir,iv,is),ir=1,mesh(is)),iv=1,nchi(is))
    !
    if (iver(1,is) == 1) then
       !
       !   old version: read the q_l(r) and fit them with the Vanderbilt's form
       !
       call fit_qrl(iunps, is)
       !
    end if
    !
    !    Here we write on output information on the pseudopotential 
    !
    WRITE( stdout,200) is
200 format (/4x,60('=')/4x,'|  pseudopotential report',               &
         &        ' for atomic species:',i3,11x,'|')
    WRITE( stdout,300) 'pseudo potential version', &
         iver(1,is), iver(2,is), iver(3,is)
300 format (4x,'|  ',1a30,3i4,13x,' |' /4x,60('-'))
    WRITE( stdout,400) title, dft_name
400 format (4x,'|  ',2a20,' exchange-corr  |')
    WRITE( stdout,500) zmesh(is), is, zv(is), exfact
500 format (4x,'|  z =',f5.0,4x,'zv(',i2,') =',f5.0,4x,'exfact =',    &
         &     f10.5, 9x,'|')
    WRITE( stdout,600) ifpcor, etotpseu
600 format (4x,'|  ifpcor = ',i2,10x,' atomic energy =',f10.5,        &
         &     ' Ry',6x,'|')
    WRITE( stdout,700)
700 format(4x,'|  index    orbital      occupation    energy',14x,'|')
    WRITE( stdout,800) ( iv, nnlz(iv), oc(iv,is), ee(iv), iv=1,nchi(is) )
800 format(4x,'|',i5,i11,5x,f10.2,f12.2,15x,'|')
    if (iver(1,is) >= 3 .and. nang > 0) then
       write(fmt,900) 2*nang-1, 40-8*(2*nang-2)
900    format('(4x,''|  rinner ='',',i1,'f8.4,',i2,'x,''|'')')
       WRITE( stdout,fmt)  (rinner(lp,is),lp=1,2*nang-1)
    end if
    WRITE( stdout,1000)
1000 format(4x,'|    new generation scheme:',32x,'|')
    WRITE( stdout,1100) nbeta(is),kkbeta(is),rcloc
1100 format(4x,'|    nbeta = ',i2,5x,'kkbeta =',i5,5x,                 &
         &     'rcloc =',f10.4,4x,'|'/                                      &
         &     4x,'|    ibeta    l     epsilon   rcut',25x,'|')
    do iv = 1, nbeta(is)
       lp=lll(iv,is)+1
       WRITE( stdout,1200) iv,lll(iv,is),eee(iv),rc(lp)
1200   format(4x,'|',5x,i2,6x,i2,4x,2f7.2,25x,'|')
    enddo
    WRITE( stdout,1300)
1300 format (4x,60('='))
    !
    return
100 call errore('readvan','error reading pseudo file', abs(ios) )
  end subroutine readvan

  !-----------------------------------------------------------------------
  subroutine fit_qrl(iunps, is)
    !-----------------------------------------------------------------------
    !
    ! find coefficients qfcoef that fit the pseudized qrl in US PP
    ! these coefficients are written to file in newer versions of the 
    ! Vanderbilt PP generation code but not in some ancient versions
    !
    !
    implicit none
    integer, intent(in) :: iunps, is
    !
    real (kind=DP), allocatable :: qrl(:,:), a(:,:), ainv(:,:), b(:), x(:)
    real (kind=DP) :: deta
    integer :: iv, jv, ijv, lmin, lmax, l, ir, irinner, i,j
    !
    !
    allocate ( a(nqf(is),nqf(is)), ainv(nqf(is),nqf(is)) )
    allocate ( b(nqf(is)), x(nqf(is)) )
    ALLOCATE ( qrl(kkbeta(is), nqlc(is)) )
    !
    ijv = 0
    do iv=1,nbeta(is)
       do jv=iv,nbeta(is)
          ijv = ijv + 1
          !
          ! original version, assuming lll(jv) >= lll(iv) 
          !   lmin=lll(jv,is)-lll(iv,is)+1
          !   lmax=lmin+2*lll(iv,is)
          ! note that indices run from 1 to Lmax+1, not from 0 to Lmax
          !
          lmin = ABS( lll(jv,is) - lll(iv,is) ) + 1
          lmax =      lll(jv,is) + lll(iv,is)   + 1
          IF ( lmin < 1 .OR. lmax >  SIZE(qrl,2)) &
               CALL errore ('fit_qrl', 'bad 2rd dimension for array qrl', 1)
          !
          !  read q_l(r) for all l
          !
          read(iunps,*, err=100) &
                  ( (qrl(ir,l),ir=1,kkbeta(is)), l=lmin,lmax)
          !
          do l=lmin,lmax
             !
             ! reconstruct rinner
             !
             do ir=kkbeta(is),1,-1
                if ( abs(qrl(ir,l)-qfunc(ir,iv,jv,is)) > 1.0d-6) go to 10
             end do
10           irinner = ir+1
             rinner(l,is) = r(irinner,is)
             !
             ! least square minimization: find
             ! qrl = sum_i c_i r^{l+1}r^{2i-2} for r < rinner
             !
             a(:,:) = 0.d0
             b(:)   = 0.d0
             do i = 1, nqf(is)
                do ir=1,irinner
                   b(i) = b(i) + r(ir,is)**(2*i-2+l+1) * qrl(ir,l)
                end do
                do j = i, nqf(is)
                   do ir=1,irinner
                      a(i,j) = a(i,j) + r(ir,is)**(2*i-2+l+1) * &
                           r(ir,is)**(2*j-2+l+1) 
                   end do
                   if (j > i) a(j,i) = a(i,j) 
                end do
             end do
             !
             call invmat (nqf(is), a, ainv, deta)
             !
             do i = 1, nqf(is)
                qfcoef(i,l,iv,jv,is) = dot_product(ainv(i,:),b(:))
                if (iv /= jv) qfcoef(i,l,jv,iv,is) = qfcoef(i,l,iv,jv,is)
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
  !-----------------------------------------------------------------------
  SUBROUTINE herman_skillman_grid(mesh,z,r,rab)
    !-----------------------------------------------------------------------
    !
    !     generate Herman-Skillman radial grid (obsolescent)
    !     c    - 0.8853418/z**(1/3)
    !
    IMPLICIT NONE
    !
    INTEGER mesh
    REAL(DP) z, r(mesh), rab(mesh)
    !
    REAL(DP) deltax
    INTEGER nblock,i,j,k
    !
    nblock = mesh/40
    i=1
    r(i)=0.0
    deltax=0.0025*0.88534138/z**(1.d0/3.d0)
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
  subroutine readrrkj( is, iunps )
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
    implicit none
    !
    !    First the arguments passed to the subroutine
    !
    integer :: &
         is,        &! The index of the pseudopotential
         iunps       ! the unit from with pseudopotential is read
    !
    !    Local variables
    !
    integer::  iexch, icorr, igcx, igcc

    integer:: &
         nb,mb, nmb,&! counters on beta functions
         n,         &! counter on mesh points
         ir,        &! counters on mesh points
         pseudotype,&! the type of pseudopotential
         ios,       &! I/O control
         ndum,      &! dummy integer variable
         l,         &! counter on angular momentum
         ikk         ! the kkbeta for each beta
    real(DP):: &
         x,         &! auxiliary variable
         etotps,    &! total energy of the pseudoatom
         rdum        ! dummy real variable
    !
    logical :: & 
         rel      ! if true the atomic calculation is relativistic
    !
    character(len=75) :: &
         titleps    ! the title of the pseudo
    !
    real(DP) :: &
         rsatom(ndmx)    ! charge density of pseudoatom
    integer :: &
         lmax       ! max angular momentum
    !
    !     We first check the input variables
    !
    if (is <= 0 .or. is >= nsx) &
         call errore('readvan','routine called with wrong 1st argument', 1)
    if (iunps <= 0 .or. iunps >= 100000) &
         call errore('readvan','routine called with wrong 2nd argument', 1)
    !
    read( iunps, '(a75)', err=100, iostat=ios ) &
         titleps
    psd (is) = titleps(7:8)
    !
    read( iunps, '(i5)',err=100, iostat=ios ) &
         pseudotype
    tvanp(is) = (pseudotype == 3)
    
    if ( tvanp(is) .AND. meta_ionode ) then
       WRITE( stdout, &
              '(/,5X,"RRKJ3 Ultrasoft PP for ",a2,/)') titleps(7:8)
    else
       WRITE( stdout, &
              '(/,5X,"RRKJ3 norm-conserving PP for ",a2,/)') titleps(7:8)
    endif
    
    read( iunps, '(2l5)',err=100, iostat=ios ) &
         rel, nlcc(is)
    read( iunps, '(4i5)',err=100, iostat=ios ) &
         iexch, icorr, igcx,  igcc

    call set_dft_from_indices(iexch,icorr,igcx,igcc)

    read( iunps, '(2e17.11,i5)') &
         zv(is), etotps, lmax
    if ( zv(is) < 1 .or. zv(is) > 100 ) &
         call errore('readrrkj','wrong potential read',is)
    !
    read( iunps, '(4e17.11,i5)',err=100, iostat=ios ) &
         xmin (is), rdum, zmesh(is), dx (is), mesh(is)
    !
    if (mesh(is) > ndmx .or. mesh(is) < 0) &
         call errore('readrrkj', 'wrong mesh',is)
    !
    read( iunps, '(2i5)', err=100, iostat=ios ) &
         nchi(is), nbeta(is)
    !
    if (nbeta(is) > nbrx .or. nbeta(is) < 0) &
         call errore('readrrkj', 'wrong nbeta', is)
    if (nchi(is) > nchix .or. nchi(is) < 0) &
         call errore('readrrkj', 'wrong nchi', is)
    !
    read( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
         ( rdum, nb=1,nchi(is) )
    read( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
         ( rdum, nb=1,nchi(is) )
    !
    do nb=1,nchi(is)
       read(iunps,'(a2,2i3,f6.2)',err=100,iostat=ios) &
            rdum, ndum, lchi(nb,is), oc(nb,is)
       lll(nb,is)=lchi(nb,is)
       !
       ! oc < 0 distinguishes between bound states from unbound states
       !
       if (oc (nb, is) <= 0.d0) oc (nb, is) = -1.0
    enddo
    !
    kkbeta(is)=0
    do nb=1,nbeta(is)
       read ( iunps, '(i6)',err=100, iostat=ios ) ikk
       kkbeta(is)=max(kkbeta(is),ikk)
       read ( iunps, '(1p4e19.11)',err=100, iostat=ios ) &
            ( betar(ir,nb,is), ir=1,ikk)
       do ir=ikk+1,mesh(is)
          betar(ir,nb,is)=0.d0
       enddo
       do mb=1,nb
          read( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
               dion(nb,mb,is)
          dion(mb,nb,is)=dion(nb,mb,is)
          if (pseudotype.eq.3) then
             read(iunps,'(1p4e19.11)',err=100,iostat=ios) &
                  qqq(nb,mb,is)
             qqq(mb,nb,is)=qqq(nb,mb,is)
             read(iunps,'(1p4e19.11)',err=100,iostat=ios) &
                  (qfunc(n,nb,mb,is),n=1,mesh(is))
             do n=1,mesh(is)
                qfunc(n,mb,nb,is)=qfunc(n,nb,mb,is)
             enddo
          else
             qqq(nb,mb,is)=0.d0
             qqq(mb,nb,is)=0.d0
             do n=1,mesh(is)
                qfunc(n,nb,mb,is)=0.d0
                qfunc(n,mb,nb,is)=0.d0
             enddo
          endif
       enddo
    enddo
    !
    !   reads the local potential 
    !
    read( iunps, '(1p4e19.11)',err=100, iostat=ios ) &
         rdum, ( vloc_at(ir,is), ir=1,mesh(is) )
    !
    !   reads the atomic charge
    !
    read( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
         ( rho_at(ir, is), ir=1,mesh(is) )
    !
    !   if present reads the core charge
    !
    if ( nlcc(is) ) then 
       read( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
            ( rho_atc(ir,is), ir=1,mesh(is) )
    endif
    !
    !   read the pseudo wavefunctions of the atom
    !  
    read( iunps, '(1p4e19.11)', err=100, iostat=ios ) &
         ((chi(ir,nb,is),ir=1,mesh(is)),nb=1,nchi(is))
    !
    !    set several variables for compatibility with the rest of the code
    !
    nqlc(is)=2*lmax+1
    if ( nqlc(is) > lqmax .or. nqlc(is) < 0 ) &
         call errore(' readrrkj', 'Wrong  nqlc', nqlc(is) )
    do l=1,nqlc(is)
       rinner(l,is)=0.d0
    enddo
    !
    !    compute the radial mesh
    !
    do ir = 1, mesh(is)
       x = xmin(is) + DBLE(ir-1) * dx (is)
       r(ir,is) = exp(x) / zmesh(is)
       rab(ir,is) = dx(is) * r(ir,is)
    end do
    !
    !     set rho_atc(r)=rho_core(r)  (without 4*pi*r^2 factor)
    !
    if ( nlcc(is) ) then
       do ir=1,mesh(is)
          rho_atc(ir,is) = rho_atc(ir,is)/4.0/3.14159265d0/r(ir,is)**2
       enddo
    end if
    !
    return
100 call errore('readrrkj','Reading pseudo file',abs(ios))
    stop
  end subroutine readrrkj
end module read_uspp_module
