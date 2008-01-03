!---------------------------------------------------------------------
!
! Copyright (C) 2001-2002 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
module pseudo
  !
  ! All variables to be read from the UPF file
  ! (UPF = unified pseudopotential format)
  !
  integer ,parameter :: npsx = 2
  ! npsx  : maximum number of different pseudopotentials
  integer, parameter :: lmaxx  = 3, nchix  = 6, ndm = 2000
  ! lmaxx : maximum non local angular momentum in PP      
  ! nchix : maximum number of atomic wavefunctions per PP
  ! ndm   : maximum number of points in the radial mesh
  integer, parameter :: nbrx = 8, lqmax = 5, nqfx = 8
  ! nbrx  : maximum number of beta functions         
  ! lqmax : maximum number of angular momentum of Q  
  ! nqfx  : maximum number of coefficients in Q smoothing
  !
  ! pp_header
  character (len=80):: generated, date_author, comment
  character (len=2) :: psd(npsx), pseudotype
  character (len=20):: dft(npsx)
  integer :: lmax(npsx), mesh(npsx), nbeta(npsx), ntwfc(npsx)
  logical :: nlcc(npsx), isus(npsx)
  real(8) :: zp(npsx), ecutrho, ecutwfc, etotps
  real(8) :: oc(nchix,npsx)
  character(len=2) :: els(nchix,npsx)
  integer :: lchi(nchix,npsx)
  !
  ! pp_mesh
  real(8) :: r(ndm,npsx), rab(ndm,npsx)
  !   pp_nlcc
  real(8) :: rho_atc(ndm,npsx)
  !
  ! pp_local
  real(8) ::  vloc0(ndm,npsx)
  !
  ! pp_nonlocal
  ! pp_beta
  real(8) :: betar(ndm, nbrx, npsx)
  integer :: lll(nbrx,npsx), ikk2(nbrx,npsx)  
  ! pp_dij
  real(8) :: dion(nbrx,nbrx,npsx)
  ! pp_qij
  integer ::  nqf(npsx), nqlc(npsx)
  real(8) :: rinner(lqmax,npsx), qqq(nbrx,nbrx,npsx), &
       qfunc(ndm,nbrx,nbrx,npsx)
  ! pp_qfcoef
  real(8) :: qfcoef(nqfx,lqmax,nbrx,nbrx,npsx)
  !
  ! pp_pswfc
  real(8) :: chi(ndm,nchix,npsx)
  !
  ! pp_rhoatom
  real(8) :: rho_at(ndm,npsx)
end module pseudo
!
program virtual
  !---------------------------------------------------------------------
  !
  !  Read pseudopotentials in the Unified Pseudopotential Format (UPF)
  !
  implicit none
  integer :: is, ios, iunps = 4
  real (8) :: x
  character (len=256) :: filein(2), fileout
  print '('' '')'
  print '('' Generate the UPF pseudopotential for a virtual atom '')'
  print '('' combining two pseudopootentials in UPF format '')'
  print '('' '')'
  !
  do is=1,2
     print '(''  Input PP file # '',i2,'' in UPF format > '',$)', is
     read (5, '(a)', end = 20, err = 20) filein(is)
     open(unit=iunps,file=filein(is),status='old',form='formatted',iostat=ios)
     if (ios.ne.0) stop
     write (*,*) " IOS= ", ios, is, iunps
     call read_pseudo(is, iunps)
     close (unit=iunps)
     print '('' '')'
  end do
  print '('' New Pseudo = x '',a,'' + (1-x) '',a)', (trim(filein(is)), is=1,2)
10 continue
  print '('' mixing parameter x [0<x<1] = '',$)'
  read (5,*) x
  if (x<0.d0 .or. x>1)  go to 10

  call compute_virtual(x,filein)

  fileout='NewPseudo.UPF'
  print '(''Output PP file in UPF format :  '',a)', fileout

  open(unit=2,file=fileout,status='unknown',form='formatted')
  call write_upf(2)
  close (unit=2)

20 stop
end program virtual
!
!---------------------------------------------------------------------
subroutine compute_virtual(x,filein)
  use pseudo
  use upf, ONLY : &
           upf_rel => rel, upf_rcloc => rcloc, upf_nwfs => nwfs, &
           upf_oc => oc, upf_rcut => rcut, upf_rcutus => rcutus, &
           upf_epseu => epseu, upf_els => els, &
           upf_lchi => lchi, upf_nns => nns, &
           upf_generated => generated, upf_date_author => date_author, &
           upf_comment => comment, &
           upf_psd => psd, upf_pseudotype => pseudotype, &
           upf_iexch => iexch, &
           upf_icorr => icorr, &
           upf_igcx  => igcx, &
           upf_igcc => igcc, &
           upf_lmax => lmax, upf_mesh => mesh, &
           upf_nbeta => nbeta, upf_ntwfc => ntwfc, upf_nlcc => nlcc, &
           upf_zp => zp, upf_ecutrho => ecutrho, upf_ecutwfc => ecutwfc, &
           upf_etotps => etotps, upf_ocw => ocw, &
           upf_elsw => elsw, upf_lchiw =>lchiw, &
           upf_r => r, upf_rab => rab, &
           upf_rho_atc => rho_atc, &
           upf_vloc0   => vloc0, &
           upf_betar => betar, upf_lll => lll,  upf_ikk2 => ikk2, &
           upf_dion => dion, &
           upf_nqf => nqf, upf_nqlc => nqlc, &
           upf_rinner => rinner, upf_qqq => qqq, upf_qfunc => qfunc, &
           upf_qfcoef => qfcoef, &
           upf_chi => chi, &
           upf_rho_at  => rho_at
  use splinelib
  use funct, ONLY : set_dft_from_name, get_iexch, get_icorr, get_igcx, get_igcc
  implicit none
  integer :: i, j, ib
  character (len=256) :: filein(2)
  character (len=5) :: xlabel
  real (8) :: x, capel
  real (8), allocatable :: aux1(:,:), aux2(:,:)
  logical :: interpolate
  interpolate = .false.
  !
  !pp_info
  upf_rel = -1
  upf_rcloc = 0.d0
  !
  !pp_header
  upf_generated  = 'Generated using virtual.x code '
  upf_date_author= 'Author and generation date: unknown. '//&
                   'Refer to original pseudopotential files'
  write( xlabel, '(f5.3)' ) x
  upf_comment    = 'Pseudo = x '//trim(filein(1))//&
                   ' + (1-x) '//trim(filein(2))//', with x='//xlabel
  upf_psd = "Xx"
  upf_pseudotype = "NC"
  if (isus(1) .or. isus(2)) upf_pseudotype = "US"
  call set_dft_from_name(dft(1))
  upf_iexch = get_iexch()
  upf_icorr = get_icorr()
  upf_igcx  = get_igcx()
  upf_igcc  = get_igcc()
  call set_dft_from_name(dft(2))
  if (get_iexch().ne.upf_iexch .or. get_icorr().ne.upf_icorr .or. &
      get_igcx().ne.upf_igcx .or. get_igcc().ne.upf_igcc) &
      call errore ('virtual','conflicting DFT functionals',1)
  upf_lmax = max(lmax(1), lmax(2))
  if (mesh(1).ne.mesh(2) ) then
     write (*,*) " pseudopotentials have different mesh " 
     write (*,*) mesh(1),mesh(2)
     write (*,*) r(1,1), r(1,2)
     write (*,*) r(mesh(1),1),r(mesh(2),2)
     interpolate = .true.
  end if
  upf_mesh = mesh(1)
  upf_nbeta = nbeta(1)+nbeta(2)
  upf_ntwfc = ntwfc(1)
  upf_nlcc  = nlcc(1).or.nlcc(2)
  upf_ecutrho = ecutrho
  upf_ecutwfc = ecutwfc
  upf_etotps  = etotps
  allocate( upf_ocw(upf_ntwfc), upf_elsw(upf_ntwfc), upf_lchiw(upf_ntwfc) )
  upf_ocw(1:upf_ntwfc)  = oc(1:upf_ntwfc,1)
  upf_elsw(1:upf_ntwfc) = els(1:upf_ntwfc,1)
  upf_lchiw(1:upf_ntwfc) = lchi(1:upf_ntwfc,1)
  upf_zp    =  x * zp(1) + (1.d0-x) * zp(2)
  !
  !pp_mesh
  capel = 0.d0
  do i=1,upf_mesh
     capel = capel + abs(r(i,1)-r(i,2)) + abs(rab(i,1)-rab(i,2))
  end do
  if (capel.gt.1.d-6) then
     write (*,*) " pseudopotentials have different mesh " 
     interpolate = .true.
  end if
  write (*,*) "INTERPOLATE =", interpolate
  if (interpolate) call errore ("virtual", &
                  "grid interpolation is not working yet",1)

  if (interpolate) allocate ( aux1(1,mesh(1)), aux2(1,mesh(2)) )

  allocate( upf_r(upf_mesh), upf_rab(upf_mesh) )
  upf_r(1:upf_mesh)   = r(1:upf_mesh,1)
  upf_rab(1:upf_mesh) = rab(1:upf_mesh,1)
  !
  !pp_nlcc
  allocate( upf_rho_atc(upf_mesh) )
  if (interpolate) then 
     write (*,*) "interpolate rho_atc"
     aux2(1,1:mesh(2)) = rho_atc(1:mesh(2),2)
     call dosplineint( r(1:mesh(2),2), aux2, upf_r(1:upf_mesh), aux1 )
     rho_atc(1:upf_mesh,2) = aux1(1,1:upf_mesh)
     write (*,*) " done"
  end if
  upf_rho_atc(1:upf_mesh) =    x     * rho_atc(1:upf_mesh,1) + &
                            (1.d0-x) * rho_atc(1:upf_mesh,2)
  !
  !pp_local
  allocate( upf_vloc0(upf_mesh) )
  if (interpolate) then 
     write (*,*) " interpolate vloc0"
     aux2(1,1:mesh(2)) =  vloc0(1:mesh(2),2)
     write (*,*) " done a"
     call dosplineint( r(1:mesh(2),2), aux2, upf_r(1:upf_mesh), aux1 )
     write (*,*) " done b"
     vloc0(1:upf_mesh,2) = aux1(1,1:upf_mesh)
     write (*,*) " done"
  end if
  upf_vloc0(1:upf_mesh) =      x     * vloc0(1:upf_mesh,1) +  &
                            (1.d0-x) * vloc0(1:upf_mesh,2)
  !
  !pp_nonlocal
  !pp_beta
  allocate( upf_betar(upf_mesh,upf_nbeta), &
            upf_lll(upf_nbeta), upf_ikk2(upf_nbeta) )
  ib = 0
  do i=1,nbeta(1)
     ib  = ib + 1
     upf_betar(1:upf_mesh,ib) = betar(1:upf_mesh,i,1)
     upf_lll(ib)              = lll(i,1)
     upf_ikk2(ib)             = ikk2(i,1)
  end do
  do i=1,nbeta(2)
     ib  = ib + 1
     if (interpolate) then 
     write (*,*) " interpolate betar"
        aux2(1,1:mesh(2)) = betar(1:mesh(2),i,2)
        call dosplineint( r(1:mesh(2),2), aux2, upf_r(1:upf_mesh), aux1 )
        betar(1:upf_mesh,i,2) = aux1(1,1:upf_mesh)
     write (*,*) " done"
     end if
     upf_betar(1:upf_mesh,ib) = betar(1:upf_mesh,i,2)
     upf_lll(ib)              = lll(i,2)
     upf_ikk2(ib)             = ikk2(i,2)
  end do
  !
  !pp_dij
  allocate( upf_dion(upf_nbeta, upf_nbeta) )
  upf_dion(:,:) = 0.d0
  do i=1,nbeta(1)
     do j=1,nbeta(1)
        upf_dion(i,j) = x * dion(i,j,1)
     end do
  end do
  do i=1,nbeta(2)
     do j=1,nbeta(2)
        upf_dion(nbeta(1)+i,nbeta(1)+j) = (1.d0-x) * dion(i,j,2)
     end do
  end do
  !
  !pp_qij
  if (nqf(1).ne.nqf(2)) &
      call errore ("Virtual","different nqf are not implemented (yet)", 1)
  if (nqlc(1).ne.nqlc(2)) &
      call errore ("Virtual","different nqlc are not implemented (yet)", 1)
  upf_nqf = nqf(1)
  upf_nqlc = nqlc(1)
  allocate( upf_rinner(upf_nqlc), upf_qqq(upf_nbeta,upf_nbeta), &
            upf_qfunc(upf_mesh,upf_nbeta,upf_nbeta) )
  do i=1,upf_nqlc
    if(rinner(i,1).ne.rinner(i,2)) &
       call errore("Virtual","different rinner are not implemented (yet)",i)
  end do
  upf_rinner(1:upf_nqlc) = rinner(1:upf_nqlc,1)
  
  upf_qqq(:,:) = 0.d0
  upf_qfunc(:,:,:) = 0.d0
  do i=1,nbeta(1)
     do j=1,nbeta(1)
        upf_qqq(i,j) = x * qqq(i, j,1)
        upf_qfunc(1:upf_mesh,i,j) = x * qfunc(1:upf_mesh,i,j,1)
     end do
  end do
  do i=1,nbeta(2)
     do j=1,nbeta(2)
        upf_qqq(nbeta(1)+i,nbeta(1)+j) = (1.d0-x) * qqq(i, j, 2)
        if (interpolate) then 
     write (*,*) " interpolate qfunc"
           aux2(1,1:mesh(2) ) = qfunc(1:mesh(2),i,j,2)
           call dosplineint( r(1:mesh(2),2), aux2, upf_r(1:upf_mesh), aux1 )
           qfunc(1:upf_mesh,i,j,2) = aux1(1,1:upf_mesh)
     write (*,*) " done"
        end if
        upf_qfunc(1:upf_mesh,nbeta(1)+i,nbeta(1)+j) = (1.d0-x) * qfunc(1:upf_mesh,i,j,2)
     end do
  end do
  !
  !pp_qfcoef
  allocate( upf_qfcoef(upf_nqf,upf_nqlc,upf_nbeta,upf_nbeta) )
  upf_qfcoef(:,:,:,:) = 0.d0
  do i=1,nbeta(1)
     do j=1,nbeta(1)
        upf_qfcoef(1:upf_nqf,1:upf_nqlc,i,j) = &
            x * qfcoef(1:upf_nqf,1:upf_nqlc,i,j, 1)
     end do
  end do
  do i=1,nbeta(2)
     do j=1,nbeta(2)
        upf_qfcoef(1:upf_nqf,1:upf_nqlc,nbeta(1)+i,nbeta(1)+j) = &
            (1.d0-x) * qfcoef(1:upf_nqf,1:upf_nqlc,i,j, 2)
     end do
  end do
  !
  !pp_pswfc
  allocate (upf_chi(upf_mesh,upf_ntwfc) )
  upf_chi(1:upf_mesh,1:upf_ntwfc) = chi(1:upf_mesh,1:upf_ntwfc,1)
  !
  !pp_rhoatm
  allocate (upf_rho_at(upf_mesh) )
  if (interpolate) then 
     write (*,*) " interpolate rho_at"
     aux2(1,1:mesh(2)) = rho_at(1:mesh(2),2)
     call dosplineint( r(1:mesh(2),2), aux2, upf_r(1:upf_mesh), aux1 )
     rho_at(1:upf_mesh,2) = aux1(1,1:upf_mesh)
     write (*,*) " done"
  end if
  upf_rho_at(1:upf_mesh) =    x     * rho_at(1:upf_mesh,1) + &
                           (1.d0-x) * rho_at(1:upf_mesh,2) 

end subroutine compute_virtual
!
!---------------------------------------------------------------------
subroutine read_pseudo (is, iunps)  
  !---------------------------------------------------------------------
  !
  !  Read pseudopotential in the Unified Pseudopotential Format (UPF)
  !
  use pseudo
  implicit none
  !
  integer :: is, iunps  
  ! is   : index of this pseudopotential
  ! iunps: unit connected with pseudopotential file
  !
  if (is < 0 .or. is > npsx ) call errore ('read_pseudo', 'Wrong is number', 1)
  write ( *, * ) " Reading pseudopotential file in UPF format..."  
  !------->Search for Header
  call scan_begin (iunps, "HEADER", .true.)  
  call read_pseudo_header (is, iunps)  
  call scan_end (iunps, "HEADER")  

  !-------->Search for mesh information
  call scan_begin (iunps, "MESH", .true.)  
  call read_pseudo_mesh (is, iunps)  
  call scan_end (iunps, "MESH")  
  !-------->If  present, search for nlcc
  if (nlcc (is) ) then  
     call scan_begin (iunps, "NLCC", .true.)  
     call read_pseudo_nlcc (is, iunps)  
     call scan_end (iunps, "NLCC")  
  endif
  !-------->Search for Local potential
  call scan_begin (iunps, "LOCAL", .true.)  
  call read_pseudo_local (is, iunps)  
  call scan_end (iunps, "LOCAL")  
  !-------->Search for Nonlocal potential
  call scan_begin (iunps, "NONLOCAL", .true.)  
  call read_pseudo_nl (is, iunps)  
  call scan_end (iunps, "NONLOCAL")  
  !-------->Search for atomic wavefunctions
  call scan_begin (iunps, "PSWFC", .true.)  
  call read_pseudo_pswfc (is, iunps)  
  call scan_end (iunps, "PSWFC")  
  !-------->Search for atomic charge
  call scan_begin (iunps, "RHOATOM", .true.)  
  call read_pseudo_rhoatom (is, iunps)  
  call scan_end (iunps, "RHOATOM")  
  !
  write ( *, * ) " ...done"
  return
end subroutine read_pseudo
!---------------------------------------------------------------------

subroutine scan_begin (iunps, string, rew)  
  !---------------------------------------------------------------------
  !
  implicit none
  ! Unit of the input file
  integer :: iunps  
  ! Label to be matched
  character (len=*) :: string  
  logical :: rew  
  ! Flag: if .true. rewind the file
  character (len=80) :: rstring  
  ! String read from file
  integer :: ios
  logical, external :: matches 

  ios = 0
  if (rew) rewind (iunps)  
  do while (ios.eq.0)  
     read (iunps, *, iostat = ios, err = 300) rstring  
     if (matches ("<PP_"//string//">", rstring) ) return  
  enddo
300 call errore ('scan_begin', 'No '//string//' block', abs (ios) )  

end subroutine scan_begin
!---------------------------------------------------------------------

subroutine scan_end (iunps, string)  
  !---------------------------------------------------------------------
  implicit none
  ! Unit of the input file
  integer :: iunps
  ! Label to be matched
  character (len=*) :: string  
  ! String read from file
  character (len=80) :: rstring
  integer :: ios
  logical, external :: matches 

  read (iunps, '(a)', iostat = ios, err = 300) rstring  
  if (matches ("</PP_"//string//">", rstring) ) return  
300 call errore ('scan_end', &
       'No '//string//' block end statement, possibly corrupted file',  - 1)
end subroutine scan_end
!
!---------------------------------------------------------------------

subroutine read_pseudo_header (is, iunps)  
  !---------------------------------------------------------------------
  !
  use pseudo
  implicit none
  !
  integer :: is, iunps  
  !
  integer :: nv, ios, nw  
  character (len=75) :: dummy  
  logical, external :: matches 

  read (iunps, *, err = 100, iostat = ios) nv, dummy  
  read (iunps, *, err = 100, iostat = ios) psd (is), dummy  
  read (iunps, *, err = 100, iostat = ios) pseudotype
  if (matches (pseudotype, "US") ) isus (is) = .true.  
  read (iunps, *, err = 100, iostat = ios) nlcc (is), dummy  
  read (iunps, '(a20,t24,a)', err = 100, iostat = ios) dft(is), dummy
  read (iunps, * ) zp (is), dummy  
  read (iunps, * ) etotps, dummy  
  read (iunps, * ) ecutwfc, ecutrho
  read (iunps, * ) lmax (is), dummy  
  read (iunps, *, err = 100, iostat = ios) mesh (is), dummy  
  read (iunps, *, err = 100, iostat = ios) ntwfc(is), nbeta (is), dummy
  read (iunps, '(a)', err = 100, iostat = ios) dummy
  do nw = 1, ntwfc(is)
     read (iunps, * ) els (nw,is), lchi (nw, is), oc (nw, is)  
  enddo
  return  
100 call errore ('read_pseudo_header', 'Reading pseudo file', abs (ios))
end subroutine read_pseudo_header
!
!---------------------------------------------------------------------
subroutine read_pseudo_local (is, iunps)  
  !---------------------------------------------------------------------
  !
  use pseudo
  implicit none
  !
  integer :: is, iunps  
  !
  integer :: ir, ios  
  !
  read (iunps, *, err=100, iostat=ios) (vloc0(ir,is) , ir=1,mesh(is))

100 call errore ('read_pseudo_local','Reading pseudo file', abs(ios) )

  return  
end subroutine read_pseudo_local
!
!---------------------------------------------------------------------

subroutine read_pseudo_mesh (is, iunps)  
  !---------------------------------------------------------------------
  !
  use pseudo
  implicit none
  !
  integer :: is, iunps  
  !
  integer :: ir, ios
  !
  call scan_begin (iunps, "R", .false.)  
  read (iunps, *, err = 100, iostat = ios) (r(ir,is), ir=1,mesh(is) )
  call scan_end (iunps, "R")  
  call scan_begin (iunps, "RAB", .false.)  
  read (iunps, *, err = 100, iostat = ios) (rab(ir,is), ir=1,mesh(is) )
  call scan_end (iunps, "RAB")  

  return  

100 call errore ('read_pseudo_mesh', 'Reading pseudo file', abs (ios) )  
end subroutine read_pseudo_mesh
!
!---------------------------------------------------------------------

subroutine read_pseudo_nl (is, iunps)  
  !---------------------------------------------------------------------
  !
  use pseudo
  implicit none
  !
  integer :: is, iunps
  !
  integer :: nb, mb, n, ir, nd, ios, idum, ldum, icon, lp, i
  ! counters
  character (len=75) :: dummy  
  !
  do nb = 1, nbeta (is)  
     call scan_begin (iunps, "BETA", .false.)  
     read (iunps, *, err = 100, iostat = ios) idum, lll(nb,is), dummy
     read (iunps, '(i6)', err = 100, iostat = ios) ikk2(nb,is)  
     read (iunps, *, err = 100, iostat = ios) &
          (betar(ir,nb,is), ir=1,ikk2(nb,is))
     do ir = ikk2(nb,is) + 1, mesh (is)  
        betar (ir, nb, is) = 0.d0  
     enddo
     call scan_end (iunps, "BETA")  
  enddo

  call scan_begin (iunps, "DIJ", .false.)  
  read (iunps, *, err = 100, iostat = ios) nd, dummy  
  dion (:,:,is) = 0.d0
  do icon = 1, nd  
     read (iunps, *, err = 100, iostat = ios) nb, mb, dion(nb,mb,is)
     dion (mb,nb,is) = dion (nb,mb,is)  
  enddo
  call scan_end (iunps, "DIJ")  

  if (isus (is) ) then  
     call scan_begin (iunps, "QIJ", .false.)  
     read (iunps, *, err = 100, iostat = ios) nqf(is)
     nqlc (is)= 2 * lmax (is) + 1
     if (nqlc(is).gt.lqmax .or. nqlc(is).lt.0) &
          call errore (' read_pseudo_nl', 'Wrong  nqlc', nqlc (is) )
     if (nqf(is).ne.0) then
        call scan_begin (iunps, "RINNER", .false.)  
        read (iunps,*,err=100,iostat=ios) &
             (idum,rinner(i,is),i=1,nqlc(is))
        call scan_end (iunps, "RINNER")  
     end if
     do nb = 1, nbeta(is)  
        do mb = nb, nbeta(is)

           read (iunps,*,err=100,iostat=ios) idum, idum, ldum, dummy
           !"  i    j   (l)"
           if (ldum.ne.lll(mb,is) ) call errore ('read_pseudo_nl', &
                'inconsistent angular momentum for Q_ij', 1)

           read (iunps,*,err=100,iostat=ios) qqq(nb,mb,is), dummy
           ! "Q_int"
           qqq(mb,nb,is) = qqq(nb,mb,is)  

           read (iunps,*,err=100,iostat=ios) &
                        (qfunc(n,nb,mb,is), n=1,mesh(is))
           do n = 0, mesh (is)  
              qfunc(n,mb,nb,is) = qfunc(n,nb,mb,is)  
           enddo

           if (nqf(is).gt.0) then
              call scan_begin (iunps, "QFCOEF", .false.)  
              read (iunps,*,err=100,iostat=ios) &
                        ((qfcoef(i,lp,nb,mb,is),i=1,nqf(is)),lp=1,nqlc(is))
              call scan_end (iunps, "QFCOEF")  
           end if

        enddo
     enddo
     call scan_end (iunps, "QIJ")  
  else  
     qqq (:,:,is) = 0.d0
     qfunc(:,:,:,is) =0.d0
  endif

100 call errore ('read_pseudo_nl', 'Reading pseudo file', abs (ios) )  
  return  
end subroutine read_pseudo_nl
!
!---------------------------------------------------------------------
subroutine read_pseudo_nlcc (is, iunps)  
  !---------------------------------------------------------------------
  !
  use pseudo
  implicit none
  !
  integer :: is, iunps  
  !
  integer :: ir, ios  

  read (iunps, *, err = 100, iostat = ios) (rho_atc(ir,is), ir=1,mesh(is) )
  !
100 call errore ('read_pseudo_nlcc', 'Reading pseudo file', abs (ios) )  
  return  
end subroutine read_pseudo_nlcc
!
!---------------------------------------------------------------------
subroutine read_pseudo_pswfc (is, iunps)  
  !---------------------------------------------------------------------
  !
  use pseudo
  implicit none
  !
  integer :: is, iunps
  !
  character (len=75) :: dummy  
  integer :: nb, ir, ios  
  !
  do nb = 1, ntwfc(is)  
     read (iunps,*,err=100,iostat=ios) dummy  !Wavefunction labels
     read (iunps,*,err=100,iostat=ios) (chi(ir,nb,is), ir=1,mesh(is))
  enddo
100 call errore ('read_pseudo_pswfc', 'Reading pseudo file', abs(ios))
  return  

end subroutine read_pseudo_pswfc
!
!---------------------------------------------------------------------
subroutine read_pseudo_rhoatom (is, iunps)  
  !---------------------------------------------------------------------
  !
  use pseudo
  implicit none
  !
  integer :: is, iunps
  !
  integer :: ir, ios  

  read (iunps,*,err=100,iostat=ios) (rho_at(ir,is), ir=1,mesh(is))
  return  

100 call errore ('read_pseudo_rhoatom','Reading pseudo file',abs(ios))

end subroutine read_pseudo_rhoatom

