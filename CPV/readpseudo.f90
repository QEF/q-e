!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"

!  ----------------------------------------------
!  BEGIN manual

!=----------------------------------------------------------------------------=!
   MODULE read_pseudo_module_fpmd
!=----------------------------------------------------------------------------=!

!  this module handles the reading of pseudopotential data
!  ----------------------------------------------
!  routines in this module:
!  SUBROUTINE readpseudo( ap, upf, psdir, psfile, nsp, nspnl)
!  ----------------------------------------------
!
!  Description of the Native FPMD pseudopotential format
!
!  The format of the file must be as follows
!  (lowercase text and }'s are comments):
!
!  When POTTYP = 'ANALYTIC' the layout is:
!
!    TCC      TMIX                additional stuff on each line is ignored
!    POTTYP   LLOC LNL ( INDL(i), i = 1, LNL )
!    ( WGV(i), i = 1, LNL )       this line only if tmix(is) is true
!    ZV       IGAU                igau must be 1 or 3     }
!    WRC(1) RC(1) WRC(2) RC(2)    this line if igau = 3   }
!    RC(1)                        this one if igau = 1    }
!    RCL(1,1)    AL(1,1)    BL(1,1)         }             }  this
!     ...         ...        ...            }  l = 0      }  section
!    RCL(IGAU,1) AL(IGAU,1) BL(IGAU,1)      }             }  only if
!    RCL(1,2)    AL(1,2)    BL(1,2)      }                }  pottyp is
!     ...         ...        ...         }     l = 1      }  'ANALYTIC'
!    RCL(IGAU,2) AL(IGAU,2) BL(IGAU,2)   }                }
!    RCL(1,3)    AL(1,3)    BL(1,3)         }             }
!     ...         ...        ...            }  l = 2      }
!    RCL(IGAU,3) AL(IGAU,3) BL(IGAU,3)      }             }
!    NMESH NCHAN                                       }
!    RW( 1 )     ( RPS( 1, j ), j = 1, NCHAN )         }  pseudowave
!     ...         ...              ...                 }
!    RW( NMESH ) ( RPS( NMESH, j ), j = 1, NCHAN )     }
!
!  
!  When POTTYP = 'NUMERIC' the layout is:
!
!    TCC      TMIX             additional stuff on each line is ignored
!    POTTYP   LLOC LNL  ( INDL(i), i = 1, LNL )
!    ( WGV(i), i = 1, LNL )       this line only if tmix(is) is true
!    ZV                                             }
!    NMESH NCHAN                                    }    this if
!    RW( 1 )     ( VR( 1, j ), j = 1, NCHAN )       }    pottyp is
!     ...       ...             ...                 }    'NUMERIC'
!    RW( NMESH ) ( VR( NMESH, j ), j = 1, NCHAN )   }
!    NMESH NCHAN                                       }
!    RW( 1 )     ( RPS( 1, j ), j = 1, NCHAN )         }  pseudowave
!     ...         ...              ...                 }
!    RW( NMESH ) ( RPS( NMESH, j ), j = 1, NCHAN )     }
!
!  DETAILED DESCRIPTION OF INPUT PARAMETERS:
!
!    TCC      (logical)   True if Core Correction are required for this 
!                         pseudo
!
!    TMIX     (logical)   True if we want to mix nonlocal pseudopotential 
!                         components 
!
!    WGV(i)   (real)      wheight of the nonlocal components in the 
!                         pseudopotential mixing scheme 
!                         These parameters are present only if TMIX = .TRUE.
!                         1 <= i <= LNL
!                           
!    POTTYP   (character) pseudopotential type
!                         pottyp = 'ANALYTIC' : use an analytic expression
!                         pottyp = 'NUMERIC'  : read values from a table
!
!    ZV       (integer)   valence for each species
!
!    IGAU     (integer)   number of Gaussians in the pseudopotentials
!                         expression used only if pottyp='ANALYTIC'
!
!  parameters from Bachelet-Hamann-Schluter's table:
!
!    WRC(2)   (real)      c1, c2 (core)  parameters
!    RC(2)    (real)      alpha1, alpha2 parameters
!
!    RCL(i,3) (real)      alpha1, alpha2, alpha3 for each angular momentum
!                         1 <= i <= IGAU
!    AL(i,3)  (real)      parameters for each angular momentum
!                         1 <= i <= IGAU
!    BL(i,3)  (real)      parameters for each angular momentum
!                         1 <= i <= IGAU
!
!  nonlocality
!    IGAU     (integer)   number of Gaussians for analytic pseudopotentials
!    LLOC     (integer)   index of the angular momentum component added to 
!                          the local part  ( s = 1, p = 2, d = 3 )
!    LNL      (integer)   number of non local component
!    INDL(i)  (integer)   indices of non local components
!                         1 <= i <= LNL
!                         ( 1 3 means s and d taken as non local )
!
!  pseudo grids
!    NMESH    (integer)   number of points in the mesh mesh
!    NCHAN    (integer)   numbero of colums, radial components
!    RW(i)    (real)      distance from the core in A.U. (radial mesh)
!                         1 <= i <= NMESH
!    RPS(i,j) (real)      Atomic pseudo - wavefunctions
!                         1 <= i <= NMESH ; 1 <= j <= NCHAN
!    VP(i,j)  (real)      Atomic pseudo - potential
!                         1 <= i <= NMESH ; 1 <= j <= NCHAN
!
!  ----------------------------------------------
!  END manual


! ...   declare modules

        USE kinds, ONLY: dbl
        USE io_files, ONLY: pseudounit

        IMPLICIT NONE

        SAVE

        PRIVATE

        REAL(dbl) :: TOLMESH = 1.d-5

        PUBLIC :: read_pseudo_fpmd

!=----------------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------------=!

!  subroutines
!  ----------------------------------------------
!  ----------------------------------------------

   SUBROUTINE read_pseudo_fpmd( ap, upf, psdir, psfile, nsp, nspnl )

!  this subroutine reads pseudopotential parameters from file
!
!  Allowed format are:
!  Native FPMD  ( 'NUMERIC' and 'ANALYTIC' )
!  Native PWSCF ( 'GIANNOZ' )
!  UPF          ( 'UPF' )
!
! INPUT:
!  nsp       number of atomic species
!  psdir     directori containing pseudopotentials
!  psfile(i) pseudopotentials file name
!            1 <= i <= nsp
!
! OUTPUT:
!  ap(i)     atomic pseudopotential structure
!  upf(i)    atomic pseudopotential structure
!            1 <= i <= nsp
!  nspnl     number of atomic species  with non local pseudopotentials
!  
!  
!  ----------------------------------------------

      USE pseudo_types, ONLY: pseudo_ncpp, pseudo_upf
      USE mp, ONLY: mp_bcast
      USE mp_global, ONLY: mpime, root, group
      USE io_global, ONLY: stdout

      IMPLICIT NONE

! ... declare subroutine arguments
      TYPE (pseudo_ncpp), INTENT(OUT) :: ap(:)
      TYPE (pseudo_upf ), INTENT(OUT) :: upf(:)
      INTEGER, INTENT(IN) :: nsp
      INTEGER, INTENT(OUT) :: nspnl
      CHARACTER(LEN=*), INTENT(IN) :: psfile(:)
      CHARACTER(LEN=*), INTENT(IN) :: psdir

! ... declare other variables
      CHARACTER(LEN=20)  :: pottyp
      CHARACTER(LEN=80)  :: error_msg
      CHARACTER(LEN=256) :: filename
      INTEGER            :: is, ierr, info
      LOGICAL            :: rdlocal

!  end of declarations
!  ----------------------------------------------

      nspnl   = 0
      rdlocal = .FALSE.

      IF( nsp > SIZE(ap) ) THEN
        CALL errore(' READPOT ',' Wrong dimensions ', MAX(1,ABS(nsp)) )
      END IF

      ierr = 0
      info = 0
      error_msg = 'none'
     
      IF(mpime == root) THEN

        WRITE( stdout,4)
    4   FORMAT(//,3X,'Atomic Pseudopotentials Parameters',/, &
                  3X,'----------------------------------' )

        DO is = 1, nsp

          IF( psdir /= ' ' ) THEN
            filename = TRIM( ADJUSTL( psdir ) ) // '/' // TRIM( ADJUSTL( psfile(is) ) ) 
          ELSE
            filename = TRIM( ADJUSTL( psfile(is) ) ) 
          ENDIF

          WRITE( stdout,6) is
          WRITE( stdout,7) filename
    6     FORMAT( /,3X,'ATOMIC PSEUDOPOTENTIAL for SPECIE : ',I2)
    7     FORMAT(   3X,'Read from file ',A60)

          OPEN( UNIT = pseudounit, FILE = filename, STATUS = 'OLD' )
          REWIND( pseudounit )

          CALL check_file_type( pseudounit, info )
          CLOSE( pseudounit )

          OPEN( UNIT = pseudounit, FILE = filename, STATUS = 'OLD' )
          REWIND( pseudounit )

          IF( info == 1 ) THEN

!  ...      Pseudopotential form is UPF
            ap(is)%pottyp = 'UPF'

          ELSE

!  ...      Read pseudopotential header
            CALL read_head_pp( pseudounit, ap(is), error_msg, info)
            IF( info /= 0 ) GO TO 200

            IF ( ap(is)%lnl > 0) THEN
              !   A non-local pseudopotential is being read, increase nspnl.
              nspnl = nspnl + 1
              !   non-local species must come before local one,
              !   write an error message if a n-l pseudopotential follow a local one
              IF( rdlocal ) THEN
                info = 20
                error_msg = ' Local pseudopotentials should follow non local ones '
                GO TO 200
              END IF
            ELSE
              !   now follow local potentials
              rdlocal = .TRUE.
            END IF

          END IF

          pottyp = ap(is)%pottyp

          IF( pottyp == 'GIANNOZ' ) THEN

            CALL read_giannoz(pseudounit, ap(is), info)
            IF( info /= 0 ) GO TO 200

          ELSE IF( pottyp == 'UPF' ) THEN

            CALL read_pseudo(pseudounit, ap(is), upf(is), info)  
            IF( info /= 0 ) GO TO 200
            IF( ap(is)%lnl > 0 ) nspnl = nspnl + 1

          ELSE IF( pottyp == 'NUMERIC' ) THEN

            CALL read_numeric_pp( pseudounit, ap(is), error_msg, info)
            IF( info /= 0 ) GO TO 200

          ELSE IF( pottyp == 'ANALYTIC' ) THEN

            CALL read_analytic_pp( pseudounit, ap(is), error_msg, info)
            IF( info /= 0 ) GO TO 200

          ELSE 

            info = 1
            error_msg = ' Pseudopotential type '//TRIM(pottyp)//' not implemented '
            GO TO 200

          END IF

          CALL ap_info( ap(is) )

          CLOSE(pseudounit)

        END DO

      END IF

200   CONTINUE
      CALL mp_bcast(info, root, group)

      IF( info /= 0 ) THEN
        ierr = 1000 * is + info
        CALL errore(' readpseudo ', error_msg, ABS(ierr) )
      END IF

      DO is = 1, nsp

        CALL mp_bcast(ap(is)%psd,root,group)
        CALL mp_bcast(ap(is)%pottyp,root,group)
        CALL mp_bcast(ap(is)%tmix,root,group)
        CALL mp_bcast(ap(is)%tnlcc,root,group)
        CALL mp_bcast(ap(is)%igau,root,group)
        CALL mp_bcast(ap(is)%lloc,root,group)
        CALL mp_bcast(ap(is)%lnl,root,group)
        CALL mp_bcast(ap(is)%indl,root,group)
        CALL mp_bcast(ap(is)%nchan,root,group)
        CALL mp_bcast(ap(is)%mesh,root,group)
        CALL mp_bcast(ap(is)%zv,root,group)
        CALL mp_bcast(ap(is)%raggio,root,group)
        CALL mp_bcast(ap(is)%dx,root,group)
        CALL mp_bcast(ap(is)%rw,root,group)
        CALL mp_bcast(ap(is)%rab,root,group)
        CALL mp_bcast(ap(is)%vnl,root,group)
        CALL mp_bcast(ap(is)%vrps,root,group)
        CALL mp_bcast(ap(is)%vloc,root,group)
        CALL mp_bcast(ap(is)%wgv,root,group)
        CALL mp_bcast(ap(is)%rc,root,group)
        CALL mp_bcast(ap(is)%wrc,root,group)
        CALL mp_bcast(ap(is)%rcl,root,group)
        CALL mp_bcast(ap(is)%al,root,group)
        CALL mp_bcast(ap(is)%bl,root,group)

        CALL mp_bcast(ap(is)%nrps,root,group)
        CALL mp_bcast(ap(is)%lrps,root,group)
        CALL mp_bcast(ap(is)%oc,root,group)
        CALL mp_bcast(ap(is)%rps,root,group)
        CALL mp_bcast(ap(is)%rhoc,root,group)

      END DO

      CALL mp_bcast( nspnl, root, group )

      RETURN
      END SUBROUTINE read_pseudo_fpmd

!=----------------------------------------------------------------------------=!

      SUBROUTINE check_file_type( iunit, info )

! ...   This sub. check if a given fortran unit 'iunit' contains a UPF pseudopot.

        INTEGER, INTENT(IN) :: iunit
        INTEGER, INTENT(OUT) :: info
        CHARACTER(LEN=80) :: dummy
        LOGICAL, EXTERNAL :: matches
        INTEGER :: ios

        info = 0
        ios  = 0
        header_loop: do while (ios == 0)
          read (iunit, *, iostat = ios, err = 200) dummy  
          if (matches ("<PP_HEADER>", dummy) ) then
            info = 1
            exit header_loop
          endif
        enddo header_loop

200     continue
        RETURN
      END SUBROUTINE check_file_type

!=----------------------------------------------------------------------------=!

      SUBROUTINE analytic_to_numeric(ap)

!       This subroutine converts an Analytic pseudo into a numeric one

        USE constants, ONLY: pi
        USE pseudo_types, ONLY: pseudo_ncpp

        IMPLICIT NONE

        TYPE (pseudo_ncpp), INTENT(INOUT) :: ap
        INTEGER   :: ir, mesh, lmax, l, n, il, ib, ll
        REAL(dbl) :: xmin, zmesh, dx, x

! ...   declare external function
        REAL(dbl) :: erf, erfc
        EXTERNAL erf, erfc

        IF( ap%mesh == 0 ) THEN
! ...     Local pseudopotential, define a logaritmic grid
          mesh  =  400
          xmin  = -5.0d0
          zmesh =  6.0d0
          dx    =  0.025d0
          DO ir = 1, mesh
            x = xmin + REAL(ir-1) * dx
            ap%rw(ir)  = EXP(x) / zmesh
          END DO
          ap%mesh = mesh
          ap%dx   = dx
          ap%rab  = ap%dx * ap%rw
        END IF

        ap%vnl  = 0.0d0
        ap%vloc = 0.0d0
        ap%vrps = 0.0d0
        do l = 1, 3
          do ir = 1, ap%mesh
            ap%vnl(ir,l)= - ( ap%wrc(1) * erf( SQRT( ap%rc(1) ) * ap%rw(ir) ) + &
                              ap%wrc(2) * erf( SQRT( ap%rc(2) ) * ap%rw(ir) )   & 
                            ) * ap%zv / ap%rw(ir)
          end do
          do ir = 1, ap%mesh
            do n = 1, ap%igau
              ap%vnl(ir,l) = ap%vnl(ir,l) + &
                             ( ap%al(n,l) + ap%bl(n,l) * ap%rw(ir)**2 ) * &
                             EXP( - ap%rcl(n,l) * ap%rw(ir)**2 )
            end do
          end do
        end do

! ...   Copy local component to a separate array
        ap%vloc(:) = ap%vnl(:,ap%lloc)
        DO l = 1, ap%lnl
          ll=ap%indl(l)  ! find out the angular momentum (ll-1) of the component stored
                         ! in position l
          ap%vrps(:,l) = ( ap%vnl(:,ll) - ap%vloc(:) ) * ap%rps(:,ll)
        END DO

        RETURN
      END SUBROUTINE

!=----------------------------------------------------------------------------=!

      SUBROUTINE read_giannoz(uni, ap, ierr)
        USE constants, ONLY: pi
        USE pseudo_types, ONLY: pseudo_ncpp
        IMPLICIT NONE
        TYPE (pseudo_ncpp), INTENT(INOUT) :: ap
        INTEGER, INTENT(IN) :: uni
        INTEGER, INTENT(OUT) :: ierr
        REAL(dbl) :: chi( SIZE(ap%rps, 1), SIZE(ap%rps, 2) )
        REAL(dbl) :: vnl( SIZE(ap%vnl, 1), SIZE(ap%vnl, 2) )
        REAL(dbl) :: rho_core( SIZE(ap%rhoc, 1) )
        REAL(dbl) :: r, ra, rb, fac
        REAL(dbl) :: oc( SIZE(ap%rps, 2) )
        REAL(dbl) :: enl( SIZE(ap%rps, 2) )
        REAL(dbl) :: zmesh, xmin, dx, etot
        REAL(dbl) :: zval
        INTEGER   :: nn(SIZE(ap%rps, 2)), ll(SIZE(ap%rps, 2))
        INTEGER   :: nwf, mesh, i, j, in1, in2, in3, in4, m
        INTEGER   :: lmax, nlc, nnl, lloc, l, il
        LOGICAL   :: nlcc
        CHARACTER(len=80) :: dft
        CHARACTER(len=4)  :: atom
        CHARACTER(len=2)  :: el( SIZE(ap%rps, 2) )
        CHARACTER(len=80) :: ppinfo
        CHARACTER(len=80) :: strdum
        CHARACTER(len=2) :: sdum1, sdum2

!
        ierr = 0

        READ(uni,fmt='(a)') dft
        READ(uni,fmt='(a4,f5.1,3i2,a2,l1,a2,i2,a)') &
          atom, zval, lmax, nlc, nnl, sdum1, nlcc, sdum2, lloc, ppinfo

        ! WRITE( stdout,*) ' DEBUG ', atom, zval,lmax, nlc, nnl, nlcc, lloc, ppinfo

        IF( (lmax+1) > SIZE(ap%vnl, 2) ) THEN
          ierr = 1
          RETURN
        END IF
        IF( (nlcc .AND. .NOT.ap%tnlcc) .OR. (.NOT.nlcc .AND. ap%tnlcc) ) THEN
          ierr = 2
          RETURN
        END IF

        READ(uni,fmt='(f8.2,f8.4,f10.6,2i6)') zmesh, xmin, dx, mesh, nwf

        IF( mesh > SIZE(ap%rps, 1) ) THEN
          ierr = 3
          RETURN
        END IF
        IF( nwf > SIZE(ap%rps, 2) ) THEN
          ierr = 4
          RETURN
        END IF

        DO j = 0, lmax
           READ(uni,fmt="(A16,i1)") strdum, l
           READ(uni,'(4e16.8)') (vnl(i,j+1), i=1,mesh)
        END DO
        IF (nlcc) THEN
          READ(uni,fmt='(4e16.8)') (rho_core(i), i=1,mesh)
        END IF   
        DO j = 1, nwf
          READ(uni,fmt="(A16,a2)") strdum,el(j)
          READ(uni,fmt='(i5,f6.2)') ll(j),oc(j)
          READ(uni,fmt='(4e16.8)') (chi(i,j), i=1,mesh)
        END DO

        ap%zv = zval
        ap%nchan = lmax+1
        ap%mesh = mesh
        ap%rw = 0.0d0
        ap%vnl = 0.0d0
        ap%vrps = 0.0d0
        fac = 0.5d0

        ! WRITE( stdout,*) ' DEBUG ', ap%lloc, ap%numeric, ap%lnl, ap%raggio, ap%zv

        DO i = 1, mesh
          r = EXP(xmin+REAL(i-1)*dx)/zmesh
          ap%rw(i) = r
          DO j = 1, lmax+1
            ap%vnl(i,j) = vnl(i,j) * fac
          END DO
        END DO
        IF( MINVAL( ap%rw(1:mesh) ) <= 0.0d0 ) THEN
           ierr = 5
           RETURN
        END IF
        ap%dx  = dx
        ap%rab = ap%dx * ap%rw
        ap%vloc(:) = ap%vnl(:,ap%lloc)

        ap%lrps(1:nwf) = ll(1:nwf)
        ap%oc = 0.0d0
        ap%nrps = nwf
        ap%mesh = mesh
        ap%rps = 0.0d0
        fac = 1.0d0/SQRT(4.0d0*pi)
        fac = 1.0d0
        DO i = 1, mesh
          r = EXP(xmin+REAL(i-1)*dx)/zmesh
          DO j = 1, nwf
            ap%rps(i,j) = chi(i,j) * fac
          END DO
        END DO

        DO l = 1, ap%lnl
          il=ap%indl(l)  ! find out the angular momentum (il-1) of the component stored
                         ! in position l
          DO i = 1, mesh
            ap%vrps(i,l) = ( ap%vnl(i,il) - ap%vloc(i) ) * ap%rps(i,il)
          END DO
        END DO

        IF( nlcc ) THEN
          ap%rhoc = 0.0d0
          DO i = 1, mesh
            r = EXP(xmin+REAL(i-1)*dx)/zmesh
            ap%rhoc(i) = rho_core(i)
          END DO
        END IF

        RETURN
      END SUBROUTINE 

!=----------------------------------------------------------------------------=!


      SUBROUTINE ap_info( ap )
        USE pseudo_types, ONLY: pseudo_ncpp
        USE io_global, ONLY: stdout

        TYPE (pseudo_ncpp), INTENT(IN) :: ap
        INTEGER   :: in1, in2, in3, in4, m, il, ib, l, i

        IF (ap%lnl > 0) THEN
          WRITE( stdout,10) ap%pottyp
          IF (ap%tmix) THEN
            WRITE( stdout,107) 
            WRITE( stdout,106)  (ap%indl(l),l=1,ap%lnl)
            WRITE( stdout,105)  (ap%wgv(l),l=1,ap%lnl)
          ELSE
            WRITE( stdout,50) ap%lloc
          END IF
          WRITE( stdout,60) (ap%indl(l),l=1,ap%lnl)
        ELSE
! ...     A local pseudopotential has been read.
          WRITE( stdout,11) ap%pottyp
          WRITE( stdout,50) ap%lloc
        END IF

   10   FORMAT(   3X,'Type is ',A10,' and NONLOCAL. ')
  107   FORMAT(   3X,'Mixed reference potential:')
  106   FORMAT(   3X,'  L     :',3(9X,i1))
  105   FORMAT(   3X,'  Weight:',3(2X,F8.5))
   50   FORMAT(   3X,'Local component is ..... : ',I3)
   60   FORMAT(   3X,'Non local components are : ',4I3)
   11   FORMAT(   3X,'Type is ',A10,' and LOCAL. ')
   20   FORMAT(   3X,'Pseudo charge : ',F8.3,', pseudo radius : ',F8.3)

        WRITE( stdout,20) ap%zv, ap%raggio

        IF( ap%pottyp /= 'ANALYTIC' ) THEN

          WRITE( stdout,131) ap%nchan, ap%mesh, ap%dx
          in1=1
          in2=ap%mesh/4
          in3=ap%mesh/2
          in4=ap%mesh
          WRITE( stdout,132)
          WRITE( stdout,120) in1,ap%rw(in1),(ap%vnl(in1,m),m=1,ap%nchan)
          WRITE( stdout,120) in2,ap%rw(in2),(ap%vnl(in2,m),m=1,ap%nchan)
          WRITE( stdout,120) in3,ap%rw(in3),(ap%vnl(in3,m),m=1,ap%nchan)
          WRITE( stdout,120) in4,ap%rw(in4),(ap%vnl(in4,m),m=1,ap%nchan)
  131     FORMAT(/, 3X,'Pseudopotentials Grid    : Channels = ',I2,&
                   ', Mesh = ',I5,/,30X,'dx   = ',F16.14)
  132     FORMAT(   3X,'point      radius        pseudopotential')
  120     FORMAT(I8,F15.10,5F10.6)

        ELSE

          WRITE( stdout,25) ap%igau
          WRITE( stdout,30)
          WRITE( stdout,104) ap%wrc(1),ap%rc(1),ap%wrc(2),ap%rc(2)
   25     FORMAT(/, 3X,'Gaussians used : ',I2,'. Parameters are : ')
   30     FORMAT(   3X,'C (core), Alfa(core) : ')
  104     FORMAT(4(3X,F8.4))

          WRITE( stdout,40)
          DO il=1,3
            DO ib=1,ap%igau
              WRITE( stdout,103) ap%rcl(ib,il),ap%al(ib,il),ap%bl(ib,il)
            END DO
          END DO
   40     FORMAT(   3X,'Hsc radii and coeff. A and B :')
  103     FORMAT(3X,F8.4,2(3X,F15.7))


        END IF

        IF( ap%nrps > 0 .AND. ap%mesh > 0 ) THEN
          WRITE( stdout,141) ap%nrps, ap%mesh, ap%dx
          in1=1
          in2=ap%mesh/4
          in3=ap%mesh/2
          in4=ap%mesh
          WRITE( stdout,145) (ap%oc(i),i=1,ap%nrps)
          WRITE( stdout,142)
          WRITE( stdout,120) in1,ap%rw(in1),(ap%rps(in1,m),m=1,ap%nrps)
          WRITE( stdout,120) in2,ap%rw(in2),(ap%rps(in2,m),m=1,ap%nrps)
          WRITE( stdout,120) in3,ap%rw(in3),(ap%rps(in3,m),m=1,ap%nrps)
          WRITE( stdout,120) in4,ap%rw(in4),(ap%rps(in4,m),m=1,ap%nrps)
        END IF

  141   FORMAT(/, 3X,'Atomic wavefunction Grid : Channels = ',I2,&
                   ', Mesh = ',I5,/,30X,'dx   = ',F16.14)
  142   FORMAT(   3X,'point      radius        wavefunction')
  145   FORMAT(   3X,'Channels occupation number : ',5F10.4)

        IF( ap%tnlcc ) THEN
          WRITE( stdout,151) ap%mesh, ap%dx
          in1 = 1
          in2 = ap%mesh / 4
          in3 = ap%mesh / 2
          in4 = ap%mesh
          WRITE( stdout,152)
          WRITE( stdout,120) in1,ap%rw(in1),ap%rhoc(in1)
          WRITE( stdout,120) in2,ap%rw(in2),ap%rhoc(in2)
          WRITE( stdout,120) in3,ap%rw(in3),ap%rhoc(in3)
          WRITE( stdout,120) in4,ap%rw(in4),ap%rhoc(in4)
        END IF

  151   FORMAT(/, 3X,'Core correction Grid     : Mesh = ',I5, &
             ', dx   = ',F16.14)
  152   FORMAT(   3X,'point      radius        rho core')

        RETURN
      END SUBROUTINE 

!=----------------------------------------------------------------------------=!

      REAL(dbl) FUNCTION calculate_dx( a, m )
        REAL(dbl), INTENT(IN) :: a(:)
        INTEGER, INTENT(IN) :: m 
        INTEGER :: n
        REAL(dbl) :: ra, rb 
          n = MIN( SIZE( a ), m )
          ra = a(1)
          rb = a(n)
          calculate_dx = LOG( rb / ra ) / REAL( n - 1 )
        RETURN
      END FUNCTION 

!=----------------------------------------------------------------------------=!

SUBROUTINE read_atomic_wf( iunit, ap, err_msg, ierr)

  USE pseudo_types, ONLY: pseudo_ncpp
  USE parser, ONLY: field_count

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iunit
  TYPE (pseudo_ncpp), INTENT(INOUT) :: ap
  CHARACTER(LEN=*) :: err_msg
  INTEGER, INTENT(OUT) :: ierr
!
  CHARACTER(LEN=80) :: input_line
  INTEGER :: i, j, m, strlen, info, nf, mesh
  REAL(dbl) :: rdum

! ... read atomic wave functions
! ... nchan : indicate number of atomic wave functions ( s p d )

  ierr = 0
  err_msg = ' error while reading atomic wf '

  ap%rps  = 0.0_dbl
  ap%nrps = 0
  ap%oc   = 0.0d0
  ap%lrps = 0

  ! this is for local pseudopotentials
  IF( ap%lnl == 0 ) RETURN
              
  READ(iunit,'(A80)',end=100) input_line
  CALL field_count(nf, input_line)

  strlen = len_trim(input_line)

  IF( nf == 2 ) THEN
    READ(input_line(1:strlen),*,IOSTAT=ierr) mesh, ap%nrps
  ELSE
    READ(input_line(1:strlen),*,IOSTAT=ierr) mesh, ap%nrps, ( ap%oc(j), j=1, MIN(ap%nrps,SIZE(ap%oc)) )
  END IF
  IF( ap%nrps > SIZE(ap%rps,2) ) THEN
    ierr = 2   
    err_msg = ' NCHAN NOT PROGRAMMED '
    GO TO 110
  END IF
  IF( mesh > SIZE(ap%rw) .OR. mesh < 0) THEN
    ierr = 4
    err_msg = ' WAVMESH OUT OF RANGE '
    GO TO 110
  END IF

  DO j = 1, mesh
    READ(iunit,*,IOSTAT=ierr) rdum, (ap%rps(j,m),m=1,ap%nrps)
    IF( ap%mesh == 0 ) ap%rw(j) = rdum
    IF( ABS(rdum - ap%rw(j))/(rdum+ap%rw(j)) > TOLMESH ) THEN
      ierr = 5
      err_msg = ' radial meshes do not match '
      GO TO 110
    END IF
  END DO

  IF( ap%mesh == 0 ) THEN
    ap%mesh = mesh
    ap%dx = calculate_dx( ap%rw, ap%mesh )
    ap%rab  = ap%dx * ap%rw
  END IF

  GOTO 110
100 ierr = 1
110 CONTINUE
  
  RETURN
END SUBROUTINE

!=----------------------------------------------------------------------------=!

SUBROUTINE read_numeric_pp( iunit, ap, err_msg, ierr)
  USE pseudo_types, ONLY: pseudo_ncpp
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iunit
  TYPE (pseudo_ncpp), INTENT(INOUT) :: ap
  CHARACTER(LEN=*) :: err_msg
  INTEGER, INTENT(OUT) :: ierr
!
  CHARACTER(LEN=80) :: input_line
  INTEGER :: i, j, m, strlen, info, nf, l, ll

! ... read numeric atomic pseudopotential
! ... nchan : indicate number of atomic wave functions ( s p d )

  ierr = 0
  err_msg = ' error while reading atomic numeric pseudo '

  IF(ap%tmix) THEN
    READ(iunit,*) (ap%wgv(l),l=1,ap%lnl)
  END IF

  READ(iunit,*,IOSTAT=ierr) ap%zv
  READ(iunit,*,IOSTAT=ierr) ap%mesh, ap%nchan

  IF((ap%nchan > SIZE(ap%vnl,2) ) .OR. (ap%nchan < 1)) THEN
    ierr = 1
    err_msg = ' NCHAN NOT PROGRAMMED '
    GO TO 110
  END IF
  IF((ap%mesh > SIZE(ap%rw) ) .OR. (ap%mesh < 0)) THEN
    info = 2
    err_msg = ' NPOTMESH OUT OF RANGE '
    GO TO 110
  END IF

  ap%rw = 0.0d0
  ap%vnl = 0.0d0
  ap%vloc = 0.0d0
  ap%vrps = 0.0d0
  DO j = 1, ap%mesh
    READ(iunit,*,IOSTAT=ierr) ap%rw(j), (ap%vnl(j,l),l=1,ap%nchan)
  END DO

  IF( MINVAL( ap%rw(1:ap%mesh) ) <= 0.0d0 ) THEN
    info = 30
    err_msg = ' ap rw too small '
    GO TO 110
  END IF

! ...  mixed reference potential is in vr(lloc)
  IF(ap%tmix) THEN
    DO j=1,ap%mesh
      ap%vnl(j,ap%lloc)= 0.d0
      DO l=1,ap%nchan
        IF(l /= ap%lloc) THEN
          ap%vnl(j,ap%lloc)=  ap%vnl(j,ap%lloc) + ap%wgv(l) * ap%vnl(j,l)
        END IF
      END DO
    END DO
  END IF
  ap%vloc(:) = ap%vnl(:,ap%lloc)
  ap%dx = calculate_dx( ap%rw, ap%mesh )
  ap%rab  = ap%dx * ap%rw

  CALL read_atomic_wf( iunit, ap, err_msg, ierr)
  IF( ierr /= 0 ) GO TO 110

  DO l = 1, ap%lnl
    ll=ap%indl(l) 
    ap%vrps(:,l) = ( ap%vnl(:,ll) - ap%vloc(:) ) * ap%rps(:,ll)
  END DO

  IF(ap%tnlcc) THEN
    CALL read_atomic_cc( iunit, ap,  err_msg, ierr)
    IF( ierr /= 0 ) GO TO 110
  END IF

  GOTO 110
100 ierr = 1
110 CONTINUE
  
  RETURN
END SUBROUTINE

!=----------------------------------------------------------------------------=!

subroutine read_pseudo (iunps, ap, upf, ierr)  
  !
  !   read a pseudopotential in the Unified Pseudopotential Format
  !   from unit "iunps" - convert and copy to internal FPMD variables
  !   return error code in "ierr" (success: ierr=0)
  !
  use pseudo_types
  use read_pseudo_module, only: read_pseudo_upf
  use control_flags, only: tuspp
  !
  implicit none
  !
  integer :: is, iunps, ierr 
  TYPE (pseudo_ncpp), INTENT(INOUT) :: ap
  TYPE (pseudo_upf ), INTENT(INOUT) :: upf
  !
  !
  call read_pseudo_upf(iunps, upf, ierr)
  !
  if ( ierr /= 0 ) return
  !
  IF( upf%tvanp ) THEN
    CALL upf2uspp( upf, is )
    tuspp = .TRUE.
  ELSE
    CALL upf2ncpp( upf, ap )
  END IF
  !
  RETURN
  !
end subroutine read_pseudo

!=----------------------------------------------------------------------------=!

SUBROUTINE upf2uspp( upf, is )
  !
  !   convert and copy upf ultra-soft pseudo to internal modules
  !
  use pseudo_types
  use uspp_param, only: qfunc, qfcoef, rinner, qqq, vloc_at, &
                   lll, nbeta, kkbeta,  nqlc, nqf, betar, dion
  use atom, only: chi, lchi, nchi, rho_atc, r, rab, mesh, nlcc
  use ions_base, only: zv
  use cvan, only: ipp
  use funct, only: dft, which_dft

  INTEGER, INTENT(INOUT) :: is
  TYPE (pseudo_upf ), INTENT(INOUT) :: upf

  integer :: nb, exfact

  zv(is)  = upf%zp
  ! psd (is)= upf%psd
  ! tvanp(is)=upf%tvanp
  if (upf%tvanp) then
     ipp(is) = -2
  else
     ipp(is) = +4
  end if
  nlcc(is) = upf%nlcc
  !
  dft = upf%dft
  call which_dft (upf%dft)
  !
  mesh(is) = upf%mesh
  if (mesh(is) > ndmx ) call errore('read_pseudo','increase mmaxx',mesh(is))
  !
  nchi(is) = upf%nwfc
  lchi(1:upf%nwfc, is) = upf%lchi(1:upf%nwfc)
  ! oc(1:upf%nwfc, is) = upf%oc(1:upf%nwfc)
  chi(1:upf%mesh, 1:upf%nwfc, is) = upf%chi(1:upf%mesh, 1:upf%nwfc)

  !
  nbeta(is)= upf%nbeta
  kkbeta(is)=0
  do nb=1,upf%nbeta
     kkbeta(is)=max(upf%kkbeta(nb),kkbeta(is))
  end do
  betar(1:upf%mesh, 1:upf%nbeta, is) = upf%beta(1:upf%mesh, 1:upf%nbeta)
  dion(1:upf%nbeta, 1:upf%nbeta, is) = upf%dion(1:upf%nbeta, 1:upf%nbeta)
  !

  ! lmax(is) = upf%lmax
  nqlc(is) = upf%nqlc
  nqf (is) = upf%nqf
  lll(1:upf%nbeta,is) = upf%lll(1:upf%nbeta)
  rinner(1:upf%nqlc,is) = upf%rinner(1:upf%nqlc)
  qqq(1:upf%nbeta,1:upf%nbeta,is) = upf%qqq(1:upf%nbeta,1:upf%nbeta)
  qfunc (1:upf%mesh, 1:upf%nbeta, 1:upf%nbeta, is) = &
       upf%qfunc(1:upf%mesh,1:upf%nbeta,1:upf%nbeta)
  qfcoef(1:upf%nqf, 1:upf%nqlc, 1:upf%nbeta, 1:upf%nbeta, is ) = &
       upf%qfcoef( 1:upf%nqf, 1:upf%nqlc, 1:upf%nbeta, 1:upf%nbeta )

  !
  r  (1:upf%mesh, is) = upf%r  (1:upf%mesh)
  rab(1:upf%mesh, is) = upf%rab(1:upf%mesh)
  !
  if ( upf%nlcc) then
     rho_atc (1:upf%mesh, is) = upf%rho_atc(1:upf%mesh)
  else
     rho_atc (:,is) = 0.d0
  end if
  ! rsatom (1:upf%mesh, is) = upf%rho_at (1:upf%mesh)
  ! lloc(is) = 1

  !
  vloc_at (1:upf%mesh, is) = upf%vloc(1:upf%mesh)
  !

  ! compatibility with old Vanderbilt formats
  call fill_qrl(is)

  RETURN
END SUBROUTINE


!=----------------------------------------------------------------------------=!


SUBROUTINE upf2ncpp( upf, ap )

  !
  !   convert and copy upf norm conserving pseudo to internal FPMD variables
  !

  use pseudo_types

  TYPE (pseudo_ncpp), INTENT(INOUT) :: ap
  TYPE (pseudo_upf ), INTENT(INOUT) :: upf

  integer :: l, il

  ap%rw    = 0.0d0
  ap%vnl   = 0.0d0
  ap%vrps  = 0.0d0
  ap%rps   = 0.0d0
  ap%oc    = 0.0d0
  ap%tmix  = .FALSE.
  ap%tnlcc = upf%nlcc
  !
  ap%zv   = upf%zp
  ap%lnl  = upf%nbeta

  !  assume that lloc = lmax

  ap%lloc = upf%nbeta + 1 

  !  angular momentum indl: S = 1, P = 2, ecc ... while  lll: S = 0, P = 1, ecc ...

  ap%indl( 1:upf%nbeta ) = upf%lll( 1:upf%nbeta ) + 1
  ap%nchan = upf%nbeta + 1   ! projectors and local part
  ap%mesh  = upf%mesh
  ap%rw( 1:upf%mesh )     = upf%r( 1:upf%mesh )
  ap%vnl( 1:upf%mesh, 1 ) = upf%vloc( 1:upf%mesh ) / 2.0d0  ! Rydberg to Hartree atomic units
  ap%dx   = calculate_dx( ap%rw, ap%mesh )
  ap%rab  = ap%dx * ap%rw
  ap%vloc( 1:upf%mesh ) = upf%vloc( 1:upf%mesh ) / 2.0d0
  ap%nrps = upf%nwfc
  ap%lrps( 1:upf%nwfc ) = upf%lchi( 1:upf%nwfc )
  ap%rps( 1:upf%mesh, 1:upf%nwfc ) = upf%chi( 1:upf%mesh, 1:upf%nwfc )
  
  DO l = 1, ap%lnl

    !  find out the angular momentum (il-1) of the component stored in position l

    il = ap%indl(l)  
     
    !  vrps(i, il) = ( vnl(i, il) - vloc(i) ) * rps(i, il)
     
    ap%vrps( 1:upf%mesh, l ) = upf%beta( 1:upf%mesh, il ) / 2.0d0 

  END DO
 
  IF( ap%tnlcc ) THEN
    ap%rhoc = 0.0d0
    ap%rhoc(1:upf%mesh) = upf%rho_atc(1:upf%mesh)
  END IF

  RETURN
END SUBROUTINE

!=----------------------------------------------------------------------------=!

SUBROUTINE read_head_pp( iunit, ap, err_msg, ierr)
  USE pseudo_types, ONLY: pseudo_ncpp
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iunit
  TYPE (pseudo_ncpp), INTENT(INOUT) :: ap
  CHARACTER(LEN=*) :: err_msg
  INTEGER, INTENT(OUT) :: ierr
!
  INTEGER :: i, l

! ... read pseudo header

  ierr = 0
  err_msg = ' error while reading header pseudo '

  ap%indl = 0
  READ(iunit, *) ap%tnlcc, ap%tmix
  READ(iunit, *) ap%pottyp, ap%lloc, ap%lnl, (ap%indl(l), l = 1, MIN(ap%lnl, SIZE(ap%indl)) ) 

  IF( ap%lnl > SIZE(ap%indl) .OR. ap%lnl < 0 ) THEN
    ierr = 1
    err_msg = 'LNL out of range'
    GO TO 110
  END IF
  IF( ap%lloc < 0 .OR. ap%lloc > SIZE(ap%vnl,2) ) THEN
    ierr = 3
    err_msg = 'LLOC out of range'
    GO TO 110
  END IF
  IF( ap%tmix .AND. ap%pottyp /= 'NUMERIC' ) THEN
    ierr = 4
    err_msg = 'tmix not implemented for pseudo ' // ap%pottyp
    GO TO 110
  END IF
  DO l = 2, ap%lnl
    IF( ap%indl(l) <= ap%indl(l-1)) THEN
      ierr = 5
      err_msg =' NONLOCAL COMPONENTS MUST BE GIVEN IN ASCENDING ORDER'
      GO TO 110
    END IF
  END DO
  DO l = 1, ap%lnl
    IF( ap%indl(l) == ap%lloc) THEN
      ierr = 6
      err_msg = ' LLOC.EQ.L NON LOCAL!!' 
      GO TO 110
    END IF
  END DO

  GOTO 110
100 ierr = 1
110 CONTINUE
  
  RETURN
END SUBROUTINE

!=----------------------------------------------------------------------------=!

SUBROUTINE read_analytic_pp( iunit, ap, err_msg, ierr)
  USE pseudo_types, ONLY: pseudo_ncpp
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iunit
  TYPE (pseudo_ncpp), INTENT(INOUT) :: ap
  CHARACTER(LEN=*) :: err_msg
  INTEGER, INTENT(OUT) :: ierr
!
  INTEGER :: i, l

! ... read analytic pseudo gaussians

  ierr = 0
  err_msg = ' error while reading atomic analytic pseudo '

  READ(iunit,*,IOSTAT=ierr) ap%zv, ap%igau

  ap%mesh = 0 
  ap%nchan = 0 
  ap%dx = 0.0d0
  ap%rab  = 0.0d0
  ap%rw   = 0.0d0
  ap%vnl   = 0.0d0
  ap%vloc   = 0.0d0
  ap%vrps   = 0.0d0

  SELECT CASE (ap%igau)
    CASE ( 1 )
      READ(iunit,*,IOSTAT=ierr) ap%rc(1)
      ap%wrc(1) = 1.d0
      ap%wrc(2) = 0.d0
      ap%rc(2)  = 0.d0
    CASE ( 3 )
      READ(iunit,*,IOSTAT=ierr) ap%wrc(1), ap%rc(1), ap%wrc(2), ap%rc(2)
    CASE DEFAULT
      ierr = 1
      err_msg = ' IGAU NOT PROGRAMMED '
      GO TO 110
  END SELECT

  DO l=1,3
    DO i=1,ap%igau
      READ(iunit,*,IOSTAT=ierr) ap%rcl(i,l), ap%al(i,l), ap%bl(i,l)
    END DO
  END DO

  CALL read_atomic_wf( iunit, ap, err_msg, ierr)
  IF( ierr /= 0 ) GO TO 110

  IF(ap%tnlcc) THEN
    CALL read_atomic_cc( iunit, ap, err_msg, ierr)
    IF( ierr /= 0 ) GO TO 110
  END IF

! ... Analytic pseudo are not supported anymore, conversion
! ... to numeric form is forced
  CALL analytic_to_numeric( ap )

  GOTO 110
100 ierr = 1
110 CONTINUE
  
  RETURN
END SUBROUTINE

!=----------------------------------------------------------------------------=!

SUBROUTINE read_atomic_cc( iunit, ap, err_msg, ierr )

  !  this subroutine reads core correction charge mesh

  USE pseudo_types, ONLY: pseudo_ncpp

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iunit
  TYPE (pseudo_ncpp), INTENT(INOUT) :: ap
  CHARACTER(LEN=*) :: err_msg
  INTEGER, INTENT(OUT) :: ierr
!
  CHARACTER(LEN=80) :: input_line
  INTEGER :: j, mesh
  REAL(dbl) :: rdum

! ... read atomic core

  ierr = 0
  err_msg = ' error while reading atomic core pseudo '

  ap%rhoc = 0.0d0

  READ( iunit, *, IOSTAT = ierr ) mesh
  IF( mesh > SIZE( ap%rw ) .OR. mesh < 0 ) THEN
    ierr = 17
    err_msg = '  CORE CORRECTION MESH OUT OF RANGE '
    GO TO 110
  END IF
  DO j = 1, mesh
    READ( iunit, *, IOSTAT = ierr ) rdum, ap%rhoc(j)
    IF( ap%mesh == 0 ) ap%rw(j) = rdum
    IF( ABS( rdum - ap%rw(j) ) / ( rdum + ap%rw(j) ) > TOLMESH ) THEN
      ierr = 5
      err_msg = ' core cor. radial mesh does not match '
      GO TO 110
    END IF
  END DO

  IF( ap%mesh == 0 ) THEN
    ap%mesh = mesh
    ap%dx   = calculate_dx( ap%rw, ap%mesh )
    ap%rab  = ap%dx * ap%rw
  END IF

  GOTO 110
100 ierr = 1
110 CONTINUE
  
  RETURN
END SUBROUTINE

!=----------------------------------------------------------------------------=!
   END MODULE read_pseudo_module_fpmd
!=----------------------------------------------------------------------------=!

