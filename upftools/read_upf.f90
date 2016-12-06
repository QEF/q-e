MODULE pseudo_mod
  !
  ! All variables to be read from the UPF file
  ! (UPF = unified pseudopotential format)
  !
  INTEGER ,PARAMETER :: npsx = 6
  ! npsx  : maximum number of different pseudopotentials
  INTEGER, PARAMETER :: lmaxx  = 3, nchix  = 6, ndm = 2000
  ! lmaxx : maximum non local angular momentum in PP
  ! nchix : maximum number of atomic wavefunctions per PP
  ! ndm   : maximum number of points in the radial mesh
  INTEGER, PARAMETER :: nbrx = 8, lqmax = 5, nqfx = 8
  ! nbrx  : maximum number of beta functions
  ! lqmax : maximum number of angular momentum of Q
  ! nqfx  : maximum number of coefficients in Q smoothing
  !
  ! pp_header
  CHARACTER (len=80):: generated, date_author, comment
  CHARACTER (len=2) :: psd(npsx), pseudotype
  CHARACTER (len=20):: dft(npsx)
  INTEGER :: lmax(npsx), mesh(npsx), nbeta(npsx), ntwfc(npsx)
  LOGICAL :: nlcc(npsx), isus(npsx)
  real(8) :: zp(npsx), ecutrho, ecutwfc, etotps
  real(8) :: oc(nchix,npsx)
  CHARACTER(len=2) :: els(nchix,npsx)
  INTEGER :: lchi(nchix,npsx)
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
  INTEGER :: lll(nbrx,npsx), ikk2(nbrx,npsx)
  ! pp_dij
  real(8) :: dion(nbrx,nbrx,npsx)
  ! pp_qij
  INTEGER ::  nqf(npsx), nqlc(npsx)
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
  !
  PRIVATE
  PUBLIC :: read_pseudo
  !
CONTAINS
!
!---------------------------------------------------------------------
SUBROUTINE read_pseudo (is, iunps)
  !---------------------------------------------------------------------
  !
  !  Read pseudopotential in the Unified Pseudopotential Format (UPF)
  !
  IMPLICIT NONE
  !
  INTEGER :: is, iunps
  ! is   : index of this pseudopotential
  ! iunps: unit connected with pseudopotential file
  !
  IF (is < 0 .or. is > npsx ) CALL errore ('read_pseudo', 'Wrong is number', 1)
  WRITE ( *, * ) " Reading pseudopotential file in UPF format..."
  !------->Search for Header
  CALL scan_begin (iunps, "HEADER", .true.)
  CALL read_pseudo_header (is, iunps)
  CALL scan_end (iunps, "HEADER")

  !-------->Search for mesh information
  CALL scan_begin (iunps, "MESH", .true.)
  CALL read_pseudo_mesh (is, iunps)
  CALL scan_end (iunps, "MESH")
  !-------->If  present, search for nlcc
  IF (nlcc (is) ) THEN
     CALL scan_begin (iunps, "NLCC", .true.)
     CALL read_pseudo_nlcc (is, iunps)
     CALL scan_end (iunps, "NLCC")
  ENDIF
  !-------->Search for Local potential
  CALL scan_begin (iunps, "LOCAL", .true.)
  CALL read_pseudo_local (is, iunps)
  CALL scan_end (iunps, "LOCAL")
  !-------->Search for Nonlocal potential
  CALL scan_begin (iunps, "NONLOCAL", .true.)
  CALL read_pseudo_nl (is, iunps)
  CALL scan_end (iunps, "NONLOCAL")
  !-------->Search for atomic wavefunctions
  CALL scan_begin (iunps, "PSWFC", .true.)
  CALL read_pseudo_pswfc (is, iunps)
  CALL scan_end (iunps, "PSWFC")
  !-------->Search for atomic charge
  CALL scan_begin (iunps, "RHOATOM", .true.)
  CALL read_pseudo_rhoatom (is, iunps)
  CALL scan_end (iunps, "RHOATOM")
  !
  WRITE ( *, * ) " ...done"
  RETURN
END SUBROUTINE read_pseudo
!---------------------------------------------------------------------

SUBROUTINE scan_begin (iunps, string, rew)
  !---------------------------------------------------------------------
  !
  IMPLICIT NONE
  ! Unit of the input file
  INTEGER :: iunps
  ! Label to be matched
  CHARACTER (len=*) :: string
  LOGICAL :: rew
  ! Flag: if .true. rewind the file
  CHARACTER (len=80) :: rstring
  ! String read from file
  INTEGER :: ios
  LOGICAL, EXTERNAL :: matches

  ios = 0
  IF (rew) REWIND (iunps)
  DO WHILE (ios==0)
     READ (iunps, *, iostat = ios, err = 300) rstring
     IF (matches ("<PP_"//string//">", rstring) ) RETURN
  ENDDO
300 CALL errore ('scan_begin', 'No '//string//' block', abs (ios) )

END SUBROUTINE scan_begin
!---------------------------------------------------------------------

SUBROUTINE scan_end (iunps, string)
  !---------------------------------------------------------------------
  IMPLICIT NONE
  ! Unit of the input file
  INTEGER :: iunps
  ! Label to be matched
  CHARACTER (len=*) :: string
  ! String read from file
  CHARACTER (len=80) :: rstring
  INTEGER :: ios
  LOGICAL, EXTERNAL :: matches

  READ (iunps, '(a)', iostat = ios, err = 300) rstring
  IF (matches ("</PP_"//string//">", rstring) ) RETURN
300 CALL errore ('scan_end', &
       'No '//string//' block end statement, possibly corrupted file',  - 1)
END SUBROUTINE scan_end
!
!---------------------------------------------------------------------

SUBROUTINE read_pseudo_header (is, iunps)
  !---------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  INTEGER :: is, iunps
  !
  INTEGER :: nv, ios, nw
  CHARACTER (len=75) :: dummy
  LOGICAL, EXTERNAL :: matches

  READ (iunps, *, err = 100, iostat = ios) nv, dummy
  READ (iunps, *, err = 100, iostat = ios) psd (is), dummy
  READ (iunps, *, err = 100, iostat = ios) pseudotype
  IF (matches (pseudotype, "US") ) isus (is) = .true.
  READ (iunps, *, err = 100, iostat = ios) nlcc (is), dummy
  READ (iunps, '(a20,t24,a)', err = 100, iostat = ios) dft(is), dummy
  READ (iunps, * ) zp (is), dummy
  READ (iunps, * ) etotps, dummy
  READ (iunps, * ) ecutwfc, ecutrho
  READ (iunps, * ) lmax (is), dummy
  READ (iunps, *, err = 100, iostat = ios) mesh (is), dummy
  READ (iunps, *, err = 100, iostat = ios) ntwfc(is), nbeta (is), dummy
  READ (iunps, '(a)', err = 100, iostat = ios) dummy
  DO nw = 1, ntwfc(is)
     READ (iunps, * ) els (nw,is), lchi (nw, is), oc (nw, is)
  ENDDO
  RETURN
100 CALL errore ('read_pseudo_header', 'Reading pseudo file', abs (ios))
END SUBROUTINE read_pseudo_header
!
!---------------------------------------------------------------------
SUBROUTINE read_pseudo_local (is, iunps)
  !---------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  INTEGER :: is, iunps
  !
  INTEGER :: ir, ios
  !
  READ (iunps, *, err=100, iostat=ios) (vloc0(ir,is) , ir=1,mesh(is))

100 CALL errore ('read_pseudo_local','Reading pseudo file', abs(ios) )

  RETURN
END SUBROUTINE read_pseudo_local
!
!---------------------------------------------------------------------

SUBROUTINE read_pseudo_mesh (is, iunps)
  !---------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  INTEGER :: is, iunps
  !
  INTEGER :: ir, ios
  !
  CALL scan_begin (iunps, "R", .false.)
  READ (iunps, *, err = 100, iostat = ios) (r(ir,is), ir=1,mesh(is) )
  CALL scan_end (iunps, "R")
  CALL scan_begin (iunps, "RAB", .false.)
  READ (iunps, *, err = 100, iostat = ios) (rab(ir,is), ir=1,mesh(is) )
  CALL scan_end (iunps, "RAB")

  RETURN

100 CALL errore ('read_pseudo_mesh', 'Reading pseudo file', abs (ios) )
END SUBROUTINE read_pseudo_mesh
!
!---------------------------------------------------------------------

SUBROUTINE read_pseudo_nl (is, iunps)
  !---------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  INTEGER :: is, iunps
  !
  INTEGER :: nb, mb, n, ir, nd, ios, idum, ldum, icon, lp, i
  ! counters
  CHARACTER (len=75) :: dummy
  !
  DO nb = 1, nbeta (is)
     CALL scan_begin (iunps, "BETA", .false.)
     READ (iunps, *, err = 100, iostat = ios) idum, lll(nb,is), dummy
     READ (iunps, '(i6)', err = 100, iostat = ios) ikk2(nb,is)
     READ (iunps, *, err = 100, iostat = ios) &
          (betar(ir,nb,is), ir=1,ikk2(nb,is))
     DO ir = ikk2(nb,is) + 1, mesh (is)
        betar (ir, nb, is) = 0.d0
     ENDDO
     CALL scan_end (iunps, "BETA")
  ENDDO

  CALL scan_begin (iunps, "DIJ", .false.)
  READ (iunps, *, err = 100, iostat = ios) nd, dummy
  dion (:,:,is) = 0.d0
  DO icon = 1, nd
     READ (iunps, *, err = 100, iostat = ios) nb, mb, dion(nb,mb,is)
     dion (mb,nb,is) = dion (nb,mb,is)
  ENDDO
  CALL scan_end (iunps, "DIJ")

  IF (isus (is) ) THEN
     CALL scan_begin (iunps, "QIJ", .false.)
     READ (iunps, *, err = 100, iostat = ios) nqf(is)
     nqlc (is)= 2 * lmax (is) + 1
     IF (nqlc(is)>lqmax .or. nqlc(is)<0) &
          CALL errore (' read_pseudo_nl', 'Wrong  nqlc', nqlc (is) )
     IF (nqf(is)/=0) THEN
        CALL scan_begin (iunps, "RINNER", .false.)
        READ (iunps,*,err=100,iostat=ios) &
             (idum,rinner(i,is),i=1,nqlc(is))
        CALL scan_end (iunps, "RINNER")
     ENDIF
     DO nb = 1, nbeta(is)
        DO mb = nb, nbeta(is)

           READ (iunps,*,err=100,iostat=ios) idum, idum, ldum, dummy
           !"  i    j   (l)"
           IF (ldum/=lll(mb,is) ) CALL errore ('read_pseudo_nl', &
                'inconsistent angular momentum for Q_ij', 1)

           READ (iunps,*,err=100,iostat=ios) qqq(nb,mb,is), dummy
           ! "Q_int"
           qqq(mb,nb,is) = qqq(nb,mb,is)

           READ (iunps,*,err=100,iostat=ios) &
                        (qfunc(n,nb,mb,is), n=1,mesh(is))
           DO n = 0, mesh (is)
              qfunc(n,mb,nb,is) = qfunc(n,nb,mb,is)
           ENDDO

           IF (nqf(is)>0) THEN
              CALL scan_begin (iunps, "QFCOEF", .false.)
              READ (iunps,*,err=100,iostat=ios) &
                        ((qfcoef(i,lp,nb,mb,is),i=1,nqf(is)),lp=1,nqlc(is))
              CALL scan_end (iunps, "QFCOEF")
           ENDIF

        ENDDO
     ENDDO
     CALL scan_end (iunps, "QIJ")
  ELSE
     qqq (:,:,is) = 0.d0
     qfunc(:,:,:,is) =0.d0
  ENDIF

100 CALL errore ('read_pseudo_nl', 'Reading pseudo file', abs (ios) )
  RETURN
END SUBROUTINE read_pseudo_nl
!
!---------------------------------------------------------------------
SUBROUTINE read_pseudo_nlcc (is, iunps)
  !---------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  INTEGER :: is, iunps
  !
  INTEGER :: ir, ios

  READ (iunps, *, err = 100, iostat = ios) (rho_atc(ir,is), ir=1,mesh(is) )
  !
100 CALL errore ('read_pseudo_nlcc', 'Reading pseudo file', abs (ios) )
  RETURN
END SUBROUTINE read_pseudo_nlcc
!
!---------------------------------------------------------------------
SUBROUTINE read_pseudo_pswfc (is, iunps)
  !---------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  INTEGER :: is, iunps
  !
  CHARACTER (len=75) :: dummy
  INTEGER :: nb, ir, ios
  !
  DO nb = 1, ntwfc(is)
     READ (iunps,*,err=100,iostat=ios) dummy  !Wavefunction labels
     READ (iunps,*,err=100,iostat=ios) (chi(ir,nb,is), ir=1,mesh(is))
  ENDDO
100 CALL errore ('read_pseudo_pswfc', 'Reading pseudo file', abs(ios))
  RETURN

END SUBROUTINE read_pseudo_pswfc
!
!---------------------------------------------------------------------
SUBROUTINE read_pseudo_rhoatom (is, iunps)
  !---------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  INTEGER :: is, iunps
  !
  INTEGER :: ir, ios

  READ (iunps,*,err=100,iostat=ios) (rho_at(ir,is), ir=1,mesh(is))
  RETURN

100 CALL errore ('read_pseudo_rhoatom','Reading pseudo file',abs(ios))

END SUBROUTINE read_pseudo_rhoatom

END MODULE pseudo_mod
