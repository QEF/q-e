
! Copyright (C) 2002-2003 PWSCF-FPMD-CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!  extracted from module "readpseudo" of FPMD

!=----------------------------------------------------------------------------=!
      MODULE read_pseudo_module
!=----------------------------------------------------------------------------=!

!  this module handles the reading of pseudopotential data

! ...   declare modules
        USE kinds, ONLY: dbl
        IMPLICIT NONE
        SAVE

      CONTAINS
!
!---------------------------------------------------------------------
subroutine read_pseudo_upf (iunps, upf, ierr)  
  !---------------------------------------------------------------------
  !
  !   read pseudopotential "upf" in the Unified Pseudopotential Format
  !   from unit "iunps" - return error code in "ierr" (success: ierr=0)
  !
  use pseudo_types
  !
  implicit none
  !
  integer :: iunps, ierr 
  TYPE (pseudo_upf), INTENT(INOUT) :: upf
  !
  !     Local variables
  !
  integer :: ios
  character (len=80) :: dummy  
  logical, external :: matches
  !
  !
  CALL nullify_pseudo_upf( upf )
  !
  ! First check if this pseudo-potential has spin-orbit information 
  !
  ierr = 1  
  ios = 0
  upf%has_so=.true.
  addinfo_loop: do while (ios == 0)  
     read (iunps, *, iostat = ios, err = 200) dummy  
     if (matches ("<PP_ADDINFO>", dummy) ) then  
        ierr = 0
        exit addinfo_loop
     endif
  enddo addinfo_loop
  if (ierr == 1) upf%has_so=.false. 

  !------->Search for Header
  !     This version doesn't use the new routine scan_begin
  !     because this search must set extra flags for
  !     compatibility with other pp format reading
  ierr = 1  
  ios = 0
  rewind(iunps)
  header_loop: do while (ios == 0)  
     read (iunps, *, iostat = ios, err = 200) dummy  
     if (matches ("<PP_HEADER>", dummy) ) then  
        ierr = 0
        call read_pseudo_header (upf, iunps)  
        exit header_loop
     endif
  enddo header_loop


  if (ierr .ne. 0) return
  
  call scan_end (iunps, "HEADER")  

  ! WRITE( stdout, * ) "Reading pseudopotential file in UPF format"  

  !-------->Search for mesh information
  call scan_begin (iunps, "MESH", .true.)  
  call read_pseudo_mesh (upf, iunps)  
  call scan_end (iunps, "MESH")  
  !-------->If  present, search for nlcc
  if ( upf%nlcc ) then  
     call scan_begin (iunps, "NLCC", .true.)  
     call read_pseudo_nlcc (upf, iunps)  
     call scan_end (iunps, "NLCC")  
  else
     allocate( upf%rho_atc( 0:upf%mesh ) )
     upf%rho_atc = 0.0d0
  endif
  !-------->Search for Local potential
  call scan_begin (iunps, "LOCAL", .true.)  
  call read_pseudo_local (upf, iunps)  
  call scan_end (iunps, "LOCAL")  
  !-------->Search for Nonlocal potential
  call scan_begin (iunps, "NONLOCAL", .true.)  
  call read_pseudo_nl (upf, iunps)  
  call scan_end (iunps, "NONLOCAL")  
  !-------->Search for atomic wavefunctions
  call scan_begin (iunps, "PSWFC", .true.)  
  call read_pseudo_pswfc (upf, iunps)  
  call scan_end (iunps, "PSWFC")  
  !-------->Search for atomic charge
  call scan_begin (iunps, "RHOATOM", .true.)  
  call read_pseudo_rhoatom (upf, iunps)  
  call scan_end (iunps, "RHOATOM")  
  !-------->Search for add_info
  if (upf%has_so) then
     call scan_begin (iunps, "ADDINFO", .true.)  
     call read_pseudo_addinfo (upf, iunps)  
     call scan_end (iunps, "ADDINFO")  
  endif

200 return  

end subroutine read_pseudo_upf
!---------------------------------------------------------------------


subroutine scan_begin (iunps, string, rew)  
  !---------------------------------------------------------------------
  !
  implicit none
  ! Unit of the input file
  integer :: iunps  
  ! Label to be matched
  character (len=*) :: string  
  ! String read from file
  character (len=75) :: rstring  
  ! Flag if .true. rewind the file
  logical, external :: matches
  logical :: rew  
  integer :: ios

  ios = 0
  if (rew) rewind (iunps)  
  do while (ios.eq.0)  
     read (iunps, *, iostat = ios, err = 300) rstring  
     if (matches ("<PP_"//string//">", rstring) ) return  
  enddo
  return
300 call errore ('scan_begin', 'No '//string//' block', abs (ios) )  
end subroutine scan_begin
!---------------------------------------------------------------------

subroutine scan_end (iunps, string)  
  !---------------------------------------------------------------------
  implicit none
  ! Unit of the input file
  integer :: iunps, ios
  ! Label to be matched
  character (len=*) :: string  
  ! String read from file
  character (len=75) :: rstring
  logical, external :: matches
  logical :: rew

  read (iunps, '(a)', iostat = ios, err = 300) rstring  
  if (matches ("</PP_"//string//">", rstring) ) return  
  return
300 call errore ('scan_end', &
       'No '//string//' block end statement, possibly corrupted file',  - 1)
end subroutine scan_end
!
!---------------------------------------------------------------------

subroutine read_pseudo_header (upf, iunps)  
  !---------------------------------------------------------------------
  !
  USE pseudo_types, ONLY: pseudo_upf
  USE kinds

  implicit none
  !
  TYPE (pseudo_upf), INTENT(INOUT) :: upf
  integer :: iunps  
  !
  integer :: is, ierr  
  integer :: nb, ios, nw  
  character (len=80) :: dummy  
  logical, external :: matches

  ! Version number (presently ignored)
  read (iunps, *, err = 100, iostat = ios) upf%nv , dummy  
  ! Element label
  read (iunps, *, err = 100, iostat = ios) upf%psd , dummy  
  ! Type of pseudo
  read (iunps, *, err = 100, iostat = ios) upf%typ  
  if (matches (upf%typ, "US") ) then
     upf%tvanp = .true.  
  else if (matches (upf%typ, "NC") ) then
     upf%tvanp = .false.  
  else
     call errore ('read_pseudo_header', 'unknown pseudo type', 1)
  endif

  read (iunps, *, err = 100, iostat = ios) upf%nlcc , dummy  

  read (iunps, '(a20,t24,a)', err = 100, iostat = ios) upf%dft, dummy  

  read (iunps, * ) upf%zp , dummy  
  read (iunps, * ) upf%etotps, dummy  
  read (iunps, * ) upf%ecutwfc, upf%ecutrho
  read (iunps, * ) upf%lmax , dummy  
  read (iunps, *, err = 100, iostat = ios) upf%mesh , dummy  
  read (iunps, *, err = 100, iostat = ios) upf%nwfc, upf%nbeta , dummy
  read (iunps, '(a)', err = 100, iostat = ios) dummy
  ALLOCATE( upf%els( upf%nwfc ), upf%lchi( upf%nwfc ), upf%oc( upf%nwfc ) )
  do nw = 1, upf%nwfc  
     read (iunps, * ) upf%els (nw), upf%lchi (nw), upf%oc (nw)  
  enddo

  return  

100  call errore ('read_pseudo_header', 'Reading pseudo file', abs (ios))
end subroutine read_pseudo_header

!---------------------------------------------------------------------

subroutine read_pseudo_mesh (upf, iunps)  
  !---------------------------------------------------------------------
  !
  USE kinds
  USE pseudo_types, ONLY: pseudo_upf

  implicit none
  !
  integer :: iunps  
  TYPE (pseudo_upf), INTENT(INOUT) :: upf
  !
  integer :: ir, ios
  character (len=75) :: dummy  

  ALLOCATE( upf%r( 0:upf%mesh ), upf%rab( 0:upf%mesh ) )
  upf%r   = 0.0d0
  upf%rab = 0.0d0

  call scan_begin (iunps, "R", .false.)  
  read (iunps, *, err = 100, iostat = ios) (upf%r(ir), ir=1,upf%mesh )
  call scan_end (iunps, "R")  
  call scan_begin (iunps, "RAB", .false.)  
  read (iunps, *, err = 100, iostat = ios) (upf%rab(ir), ir=1,upf%mesh )
  call scan_end (iunps, "RAB")  

  return  

100 call errore ('read_pseudo_mesh', 'Reading pseudo file', abs (ios) )  
end subroutine read_pseudo_mesh


!---------------------------------------------------------------------
subroutine read_pseudo_nlcc (upf, iunps)
  !---------------------------------------------------------------------
  !
  USE kinds
  USE pseudo_types, ONLY: pseudo_upf

  implicit none
  !
  integer :: iunps  
  TYPE (pseudo_upf), INTENT(INOUT) :: upf
  !
  integer :: ir, ios
  !
  ALLOCATE( upf%rho_atc( 0:upf%mesh ) )
  upf%rho_atc = 0.0d0

  read (iunps, *, err = 100, iostat = ios) (upf%rho_atc(ir), ir=1,upf%mesh )
  !
  return

100 call errore ('read_pseudo_nlcc', 'Reading pseudo file', abs (ios) )
  return
end subroutine read_pseudo_nlcc

!---------------------------------------------------------------------
subroutine read_pseudo_local (upf, iunps)
  !---------------------------------------------------------------------
  !
  USE kinds
  USE pseudo_types, ONLY: pseudo_upf

  implicit none
  !
  integer :: iunps  
  TYPE (pseudo_upf), INTENT(INOUT) :: upf
  !
  integer :: ir, ios
  character (len=75) :: dummy
  !
  ALLOCATE( upf%vloc( 0:upf%mesh ) )
  upf%vloc = 0.0d0

  read (iunps, *, err=100, iostat=ios) (upf%vloc(ir) , ir=1,upf%mesh )

  return

100 call errore ('read_pseudo_local','Reading pseudo file', abs(ios) )
  return
end subroutine read_pseudo_local

!---------------------------------------------------------------------

subroutine read_pseudo_nl (upf, iunps)  
  !---------------------------------------------------------------------
  !
  USE kinds
  USE pseudo_types, ONLY: pseudo_upf

  implicit none
  !
  integer :: iunps  
  TYPE (pseudo_upf), INTENT(INOUT) :: upf
  !
  integer :: nb, mb, n, ir, nd, ios, idum, ldum, icon, lp, i, ikk
  ! counters
  character (len=75) :: dummy  
  !
  if ( upf%nbeta == 0) then
     upf%nqf = 0
     upf%nqlc= 0
     ALLOCATE( upf%kkbeta( 1 ) )
     ALLOCATE( upf%lll( 1 ) )
     ALLOCATE( upf%beta( 0:upf%mesh, 1 ) )
     ALLOCATE( upf%dion( 1, 1 ) )
     ALLOCATE( upf%rinner( 1 ) )
     ALLOCATE( upf%qqq   ( 1, 1 ) )
     ALLOCATE( upf%qfunc ( 0:upf%mesh, 1, 1 ) )
     ALLOCATE( upf%qfcoef( 1, 1, 1, 1 ) )
     return
  end if
  ALLOCATE( upf%kkbeta( upf%nbeta ) )
  ALLOCATE( upf%lll( upf%nbeta ) )
  ALLOCATE( upf%beta( 0:upf%mesh, upf%nbeta ) )
  ALLOCATE( upf%dion( upf%nbeta, upf%nbeta ) )
  upf%kkbeta = 0  
  upf%lll    = 0  
  upf%beta   = 0.0d0
  upf%dion   = 0.0d0

  do nb = 1, upf%nbeta 
     call scan_begin (iunps, "BETA", .false.)  
     read (iunps, *, err = 100, iostat = ios) idum, upf%lll(nb), dummy
     read (iunps, '(i6)', err = 100, iostat = ios) ikk  
     upf%kkbeta(nb) = ikk
     read (iunps, *, err = 100, iostat = ios) (upf%beta(ir,nb), ir=1,ikk)
     call scan_end (iunps, "BETA")  
  enddo


  call scan_begin (iunps, "DIJ", .false.)  
  read (iunps, *, err = 100, iostat = ios) upf%nd, dummy  
  do icon = 1, upf%nd  
     read (iunps, *, err = 100, iostat = ios) nb, mb, upf%dion(nb,mb)
     upf%dion (mb,nb) = upf%dion (nb,mb)  
  enddo
  call scan_end (iunps, "DIJ")  


  if ( upf%tvanp ) then  
     call scan_begin (iunps, "QIJ", .false.)  
     read (iunps, *, err = 100, iostat = ios) upf%nqf
     upf%nqlc = 2 * upf%lmax  + 1
     ALLOCATE( upf%rinner( upf%nqlc ) )
     ALLOCATE( upf%qqq   ( upf%nbeta, upf%nbeta ) )
     ALLOCATE( upf%qfunc ( 0:upf%mesh, upf%nbeta, upf%nbeta ) )
     ALLOCATE( upf%qfcoef( MAX( upf%nqf,1 ), upf%nqlc, upf%nbeta, upf%nbeta ) )
     upf%rinner = 0.0d0
     upf%qqq    = 0.0d0
     upf%qfunc  = 0.0d0
     upf%qfcoef = 0.0d0
     if ( upf%nqf /= 0) then
        call scan_begin (iunps, "RINNER", .false.)  
        read (iunps,*,err=100,iostat=ios) ( idum, upf%rinner(i), i=1,upf%nqlc )
        call scan_end (iunps, "RINNER")  
     end if
     do nb = 1, upf%nbeta
        do mb = nb, upf%nbeta

           read (iunps,*,err=100,iostat=ios) idum, idum, ldum, dummy
           !"  i    j   (l)"
           if (ldum /= upf%lll(mb) ) then
             call errore ('read_pseudo_nl','inconsistent angular momentum for Q_ij', 1)
           end if

           read (iunps,*,err=100,iostat=ios) upf%qqq(nb,mb), dummy
           ! "Q_int"
           upf%qqq(mb,nb) = upf%qqq(nb,mb)  

           read (iunps, *, err=100, iostat=ios) (upf%qfunc(n,nb,mb), n=1,upf%mesh)
           do n = 0, upf%mesh 
              upf%qfunc(n,mb,nb) = upf%qfunc(n,nb,mb)  
           enddo

           if ( upf%nqf > 0 ) then
              call scan_begin (iunps, "QFCOEF", .false.)  
              read (iunps,*,err=100,iostat=ios) &
                        ( ( upf%qfcoef(i,lp,nb,mb), i=1,upf%nqf ), lp=1,upf%nqlc )
              call scan_end (iunps, "QFCOEF")  
           end if

        enddo
     enddo
     call scan_end (iunps, "QIJ")  
  else  
     upf%nqf  = 1
     upf%nqlc = 2 * upf%lmax  + 1
     ALLOCATE( upf%rinner( upf%nqlc ) )
     ALLOCATE( upf%qqq   ( upf%nbeta, upf%nbeta ) )
     ALLOCATE( upf%qfunc ( 0:upf%mesh, upf%nbeta, upf%nbeta ) )
     ALLOCATE( upf%qfcoef( upf%nqf, upf%nqlc, upf%nbeta, upf%nbeta ) )
     upf%rinner = 0.0d0
     upf%qqq    = 0.0d0
     upf%qfunc  = 0.0d0
     upf%qfcoef = 0.0d0
  endif


  return  

100 call errore ('read_pseudo_nl', 'Reading pseudo file', abs (ios) )  
end subroutine read_pseudo_nl


!---------------------------------------------------------------------
subroutine read_pseudo_pswfc (upf, iunps)  
  !---------------------------------------------------------------------
  !
  USE kinds  
  USE pseudo_types, ONLY: pseudo_upf
  !
  implicit none
  !
  integer :: iunps
  TYPE (pseudo_upf), INTENT(INOUT) :: upf
  !
  character (len=75) :: dummy  
  integer :: nb, ir, ios 

  ALLOCATE( upf%chi( 0:upf%mesh, MAX( upf%nwfc, 1 ) ) )
  upf%chi = 0.0d0
  do nb = 1, upf%nwfc  
     read (iunps, *, err=100, iostat=ios) dummy  !Wavefunction labels
     read (iunps, *, err=100, iostat=ios) ( upf%chi(ir,nb), ir=1,upf%mesh )
  enddo

  return  

100 call errore ('read_pseudo_pswfc', 'Reading pseudo file', abs(ios))
end subroutine read_pseudo_pswfc

!---------------------------------------------------------------------
subroutine read_pseudo_rhoatom (upf, iunps)  
  !---------------------------------------------------------------------
  !
  USE kinds 
  USE pseudo_types, ONLY: pseudo_upf
  !
  implicit none
  !
  integer :: iunps
  TYPE (pseudo_upf), INTENT(INOUT) :: upf
  !
  integer :: ir, ios  
  !
  ALLOCATE( upf%rho_at( 0:upf%mesh ) )
  upf%rho_at = 0.0d0
  read (iunps,*,err=100,iostat=ios) ( upf%rho_at(ir), ir=1,upf%mesh )
  !
  return  

100 call errore ('read_pseudo_rhoatom','Reading pseudo file',abs(ios))
end subroutine read_pseudo_rhoatom
!
!---------------------------------------------------------------------
  subroutine read_pseudo_addinfo (upf, iunps)
!---------------------------------------------------------------------
!
!     This routine reads from the new UPF file,
!     and the total angual momentum jjj of the beta and jchi of the
!     wave-functions.
!
USE pseudo_types, ONLY: pseudo_upf
USE kinds
implicit none
integer :: iunps

TYPE (pseudo_upf), INTENT(INOUT) :: upf
integer :: nb, ios

ALLOCATE( upf%nn(upf%nwfc), upf%rcut(upf%nwfc), upf%rcutus(upf%nwfc) )
ALLOCATE( upf%epseu(upf%nwfc), upf%jchi(upf%nwfc) )
ALLOCATE( upf%jjj(upf%nbeta) )

upf%nn=0
upf%rcut=0.d0
upf%rcutus=0.d0
upf%epseu=0.d0
upf%jchi=0.d0
do nb = 1, upf%nwfc
  read (iunps, '(a2,2i3,2f6.2)',err=100,iostat=ios) upf%els(nb),  &
       upf%nn(nb), upf%lchi(nb), upf%jchi(nb), upf%oc(nb)
enddo

upf%jjj=0.d0
do nb = 1, upf%nbeta
  read (iunps, '(i5,f6.2)', err=100,iostat=ios) upf%lll(nb), upf%jjj(nb)
enddo

read(iunps, '(4f15.8)') upf%xmin, upf%rmax, upf%zmesh, upf%dx

return
100 call errore ('read_pseudo_addinfo','Reading pseudo file', abs(ios))
end subroutine read_pseudo_addinfo

!=----------------------------------------------------------------------------=!
      END MODULE read_pseudo_module
!=----------------------------------------------------------------------------=!
