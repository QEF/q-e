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
  integer ,parameter :: npsx = 6
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
!---------------------------------------------------------------------
program read_ps
  !---------------------------------------------------------------------
  !
  !  Read pseudopotentials in the Unified Pseudopotential Format (UPF)
  !
  implicit none
  integer :: is, ios, iunps = 4
  character (len=256) :: filein
  !
  is = 0
10 print '(''  Input PP file # '',i2,'' in UPF format > '',$)', is+1
  read (5, '(a)', end = 20, err = 20) filein
  open(unit=iunps,file=filein,status='old',form='formatted',iostat=ios)
  if (ios.ne.0) stop
  is = is + 1
  call read_pseudo(is, iunps)
  close (unit=iunps)
  go to 10
20 stop
end program read_ps
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

