!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"


MODULE diis


        USE kinds

        IMPLICIT NONE
        SAVE

        PRIVATE

! ...   declare module-scope variables

        LOGICAL :: o_diis, oqnr_diis, treset_diis
        LOGICAL :: tfortr_diis
        INTEGER :: max_diis, nreset, maxnstep
        REAL(DP) :: hthrs_diis, tol_diis, delt_diis

        REAL(DP) :: diis_fthr

        ! ...   declare module-scope variables
        INTEGER    :: srot0, srot, sroti, srotf
        REAL(DP)  :: tolene, tolrho, tolrhof, tolrhoi
        REAL(DP)  :: temp_elec


        REAL(DP) :: dethr

        COMPLEX(sgl), ALLOCATABLE :: parame(:,:,:,:)
        COMPLEX(sgl),  ALLOCATABLE :: grade(:,:,:,:)

        LOGICAL    :: updstrfac
        REAL(DP)  :: gemax

        PUBLIC :: diis_setup, diis_print_info, allocate_diis
        PUBLIC :: delt_diis, deallocate_diis
        PUBLIC :: simupd
        PUBLIC :: tolrhoi, sroti, tolrho, maxnstep, nreset
        PUBLIC :: srotf, tolene, srot0, srot, treset_diis
        PUBLIC :: tolrhof, temp_elec


CONTAINS

!  ----------------------------------------------

        SUBROUTINE diis_setup( diis_fthr_inp, oqnr_diis_inp, o_diis_inp, &
            max_diis_inp, hthrs_diis_inp, tol_diis_inp, &
            maxnstep_inp, reset_diis_inp, delt_diis_inp, & 
            temp_inp, srot0_inp, sroti_inp, srotf_inp, &
            tolrho_inp, tolrhoi_inp, tolrhof_inp, tolene_inp)

            LOGICAL, INTENT(IN) :: o_diis_inp, oqnr_diis_inp
            REAL(DP), INTENT(IN) :: diis_fthr_inp
            INTEGER, INTENT(IN) :: max_diis_inp, reset_diis_inp, maxnstep_inp
            REAL(DP), INTENT(IN) :: hthrs_diis_inp, tol_diis_inp, delt_diis_inp
            INTEGER, INTENT(IN) :: srot0_inp, sroti_inp, srotf_inp
            REAL(DP), INTENT(IN) ::  temp_inp
            REAL(DP), INTENT(IN) ::  tolene_inp, tolrho_inp, tolrhof_inp, tolrhoi_inp

            temp_elec = temp_inp
            srot0 = srot0_inp
            sroti = sroti_inp
            srotf = srotf_inp
            tolrho = tolrho_inp
            tolrhoi = tolrhoi_inp
            tolrhof = tolrhof_inp
            tolene = tolene_inp

            tfortr_diis = .FALSE.
            diis_fthr = diis_fthr_inp
            IF( diis_fthr > 0.0d0 ) THEN
              tfortr_diis = .TRUE.
              tol_diis = diis_fthr_inp
            ELSE
              tol_diis = tol_diis_inp
            END IF

            o_diis = o_diis_inp
            oqnr_diis = oqnr_diis_inp
            max_diis = max_diis_inp
            nreset = reset_diis_inp
            maxnstep = maxnstep_inp
            hthrs_diis = hthrs_diis_inp
            delt_diis = delt_diis_inp
            dethr = 1.d-6

          RETURN
        END SUBROUTINE diis_setup

!  ----------------------------------------------
!  ----------------------------------------------
        SUBROUTINE allocate_diis(ngwxm,nx,nk)
          INTEGER, INTENT(IN) :: ngwxm,nx,nk
          INTEGER :: ierr
          ALLOCATE( parame( ngwxm, nx, nk, max_diis ), STAT=ierr)
          IF( ierr/=0 ) CALL errore(' allocate_diis ', ' allocating parame ', ierr)
          ALLOCATE( grade( ngwxm, nx, nk, max_diis ), STAT=ierr)
          IF( ierr/=0 ) CALL errore(' allocate_diis ', ' allocating grade ', ierr)
          RETURN
        END SUBROUTINE allocate_diis

!  ----------------------------------------------
!  ----------------------------------------------
        SUBROUTINE deallocate_diis
          INTEGER :: ierr
          IF(ALLOCATED(parame)) THEN
            DEALLOCATE(parame, STAT=ierr)
            IF( ierr/=0 ) CALL errore(' allocate_diis ', ' deallocating parame ', ierr)
          END IF
          IF(ALLOCATED(grade))  THEN
            DEALLOCATE(grade, STAT=ierr)
            IF( ierr/=0 ) CALL errore(' allocate_diis ', ' deallocating grade ', ierr)
          END IF
          RETURN
        END SUBROUTINE deallocate_diis

!  ----------------------------------------------
!  ----------------------------------------------
        SUBROUTINE diis_print_info(unit)
          USE control_flags, ONLY: t_diis, t_diis_simple, t_diis_rot
          USE charge_mix, ONLY: charge_mix_print_info
          INTEGER, INTENT(IN) :: unit
            WRITE(unit,100)
            IF( t_diis_simple ) THEN
              WRITE(unit,200)
            ELSE
              WRITE(unit,300)
            END IF
            WRITE(unit,400) max_diis
            WRITE(unit,410) maxnstep
            WRITE(unit,420) delt_diis
            WRITE(unit,430) tol_diis
            WRITE(unit,450) temp_elec
            IF( .NOT. t_diis_simple ) THEN
              CALL charge_mix_print_info(unit)
              WRITE(unit,440) tolrho, tolrhoi, tolrhof
              WRITE(unit,460) srot0, sroti, srotf
              WRITE(unit,470) tolene
            END IF

 100      FORMAT(/,3X,'Using DIIS for electronic minimization')
 200      FORMAT(  3X,'DIIS without charge mixing and wave-function rotation')
 300      FORMAT(  3X,'DIIS with charge mixing and wave-function rotation')
 400      FORMAT(  3X,'order ............................ = ',1I6)
 410      FORMAT(  3X,'maximum number of iterations ..... = ',1I6)
 420      FORMAT(  3X,'time step used for a steepest step = ',1D10.4)
 430      FORMAT(  3X,'threshold for convergence ........ = ',1D10.4)
 450      FORMAT(  3X,'electronic temperature ........... = ',1D10.4)
 440      FORMAT(  3X,'tolerances on rho change ......... = ',3D10.4)
 460      FORMAT(  3X,'steps without new potentials ..... = ',3I5)
 470      FORMAT(  3X,'tolerances on energy change ...... = ',1D10.4)

          RETURN
        END SUBROUTINE diis_print_info

!  ----------------------------------------------

      SUBROUTINE converg( c, cdesc, gemax, cnorm, tprint )

!  this routine checks for convergence, by computing the norm of the
!  gradients of wavefunctions
!  version for the Gamma point
!  ----------------------------------------------

        USE wave_base, ONLY: dotp
        USE wave_types, ONLY: wave_descriptor
        USE mp_global, ONLY: intra_image_comm
        USE io_global, ONLY: stdout
        USE mp, ONLY: mp_max
        IMPLICIT NONE

! ...   declare subroutine arguments
        COMPLEX(DP), INTENT(IN) :: c(:,:)
        TYPE (wave_descriptor), INTENT(IN) :: cdesc
        LOGICAL, INTENT(IN) :: tprint
        REAL(DP), INTENT(OUT) :: gemax, cnorm

! ...   declare other variables
        INTEGER  :: iabs, i, nb, ngw
        INTEGER, EXTERNAL :: IZAMAX
        REAL(DP) :: gemax_l

! ...   end of declarations
!  ----------------------------------------------

        ngw     = cdesc%ngwl
        nb      = cdesc%nbl( 1 )
        gemax_l = 0.d0
        cnorm   = 0.d0
        DO i = 1, nb
          iabs = IZAMAX ( ngw, c(1,i), 1)
          IF( gemax_l < ABS( c(iabs,i) ) ) THEN
            gemax_l = ABS ( c(iabs,i) )
          END IF
          cnorm = cnorm + dotp( cdesc%gzero, ngw, c(:,i), c(:,i) )
        END DO
        CALL mp_max(gemax_l, intra_image_comm)
        gemax = gemax_l
        cnorm = SQRT( cnorm / ( cdesc%nbt( 1 ) * cdesc%ngwt ) )

        IF(tprint) WRITE( stdout,100) gemax, cnorm
 100    FORMAT(' Electrons: max. comp:',g16.10,'  rms norm :',g16.10)

        RETURN
      END SUBROUTINE converg

!  ----------------------------------------------
!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE hesele(svar2, vpp, eigr, sfac, vps)

!  (describe briefly what this routine does...)
!  ----------------------------------------------

! ...  declare modules
       USE control_flags, ONLY: gamma_only
       USE cell_base, ONLY: tpiba2
       USE pseudopotential, ONLY: nspnl
       USE ions_base, ONLY: nsp, na
       USE mp_global, ONLY: intra_image_comm
       USE mp, ONLY: mp_sum, mp_max
       USE reciprocal_vectors, ONLY: gstart, gzero, ggp
       USE reciprocal_space_mesh, ONLY: gkmask_l, gkcutz_l
       USE uspp_param, only: nh
       USE uspp, only: nhtol, indv


      IMPLICIT NONE

! ... declare subroutine arguments
      REAL(DP) :: svar2
      REAL(DP) :: vps(:,:)
      REAL(DP) :: vpp(:)
      COMPLEX(DP) :: sfac(:,:)
      COMPLEX(DP) :: eigr(:,:)

! ... declare other variables
      REAL(DP)  vp,ftpi,arg
      INTEGER l,ll,i,ig,igh,is,m,j, ih, iv

! ... end of declarations
!  ----------------------------------------------

! ... calculate H0 :The diagonal approximation to the 2nd derivative matrix

! ... local potential
      vp = 0.0d0
        IF(gzero) THEN
          DO i = 1, nsp
            vp = vp + DBLE( sfac(1,i) ) * vps(1,i)
          END DO
        END IF
        CALL mp_sum(vp, intra_image_comm)

      vpp = vp

! ... preconditioning for nonlocal potential  TO BE FIXED
!      DO is = 1, nspnl
!        DO ih = 1, nh( is )
!          ll = nhtol ( ih, is ) + 1
!          iv = indv  ( ih, is ) 
!          vpp(:) = vpp(:) + na(is) * dion( ih, ih, is ) * bec( ) ** 2
!        END DO
!      END DO

! ... kinetic energy
      ftpi = 0.5d0 * tpiba2

      vpp(:) = vpp(:) + ftpi * ggp( 1:SIZE(vpp) )

      WHERE ( ABS( vpp ) .LT. hthrs_diis ) vpp = hthrs_diis
      vpp = 1.0d0 / vpp

      RETURN
      END SUBROUTINE hesele

!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE fermi_diis( ent, occ, nb, nel, eig, wke, efermi, sume, temp)
        USE electrons_module, ONLY: fermi_energy
        REAL(DP)   :: occ(:,:,:)
        REAL(DP)   :: wke(:,:,:)
        REAL(DP)   :: eig(:,:,:), efermi, sume, ent, temp, entk, qtot
        INTEGER :: ispin, nel, nb

        qtot = DBLE( nel ) 
        CALL fermi_energy( eig, occ, wke, efermi, qtot, temp, sume)

! ...   compute the entropic correction
        ent = 0.0d0
        DO ispin = 1, SIZE( occ, 3 )
            CALL entropy( occ(:,1,ispin), temp, nb, entk )
            ent = ent + entk
        END DO
      RETURN
      END SUBROUTINE fermi_diis

!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE simupd(gemax,doions,c0,cgrad,cdesc,svar0,svarm,svar2,etot,f, &
                        eigr,sfac,vps,treset,istate,cnorm,eold,ndiis,nowv)

!  this routine performs a DIIS step
!  version for the Gamma point
!  ----------------------------------------------

! ... declare modules
      USE wave_types, ONLY: wave_descriptor
      USE mp_global, ONLY: intra_image_comm
      USE io_global, ONLY: stdout
      USE mp, ONLY: mp_sum, mp_max
      USE gvecw, ONLY: ngw
      USE reciprocal_vectors, ONLY: gstart, gzero

      IMPLICIT NONE

! ... declare subroutine arguments

! max_diis (integer ) max number of stored wave functions
! oqnr_diis   (logical ) if .true. use an appx. ham. as guess
! treset (logical ) if .true. reset DIIS
! tol_diis   (REAL(DP)  ) convergence tolerance
! hthrs_diis  (REAL(DP)  ) minimum value for a hessel matrix elm.
! etot   (REAL(DP)  ) Kohn-Sham energy

      COMPLEX(DP), INTENT(INOUT) :: c0(:,:,:), cgrad(:,:,:)
      TYPE (wave_descriptor), INTENT(IN) :: cdesc
      COMPLEX(DP) :: sfac(:,:) 
      COMPLEX(DP) :: eigr(:,:) 
      LOGICAL, INTENT(OUT) :: doions
      LOGICAL, INTENT(INOUT) :: treset
      REAL(DP) :: gemax, etot, svar0, svarm, svar2
      REAL(DP) :: f(:), vps(:,:)
      REAL(DP) :: eold
      INTEGER ::  ndiis, nowv

! ... declare other variables
      COMPLEX(DP), ALLOCATABLE :: cm(:,:) !  cm(ngw,c0(1)%nb_l)
      REAL(DP),    ALLOCATABLE :: bc(:,:) !  bc(max_diis+1,max_diis+1)
      REAL(DP),    ALLOCATABLE :: vc(:)   !  vc(max_diis+1)
      REAL(DP),    ALLOCATABLE :: fm1(:)
      REAL(DP),    ALLOCATABLE :: vpp(:)  !  vpp(ngw)
      REAL(DP)    cnorm
      REAL(DP)    var2,rc0rc0,ff,vvpp,fff,rri
      REAL(DP)    rrj,rii,rij,r1,r2
      LOGICAL   teinc
      INTEGER   istate, nsize, i, j, l, k, ierr
      INTEGER, SAVE :: prevv

! ... end of declarations
!  ----------------------------------------------

! c0 : current wavefunction, on output new estimated G.S. wave functions
! cgrad : current wavefunction gradient --> common fowf
! cm : scratch space; out --> old wavefunction

      IF( treset ) istate = 0

      ALLOCATE( fm1( SIZE(f) ), STAT=ierr )
      IF( ierr/=0 ) CALL errore(' simupd ', ' allocating fm1 ', ierr)
      ALLOCATE( cm( cdesc%ldg, cdesc%ldb ), STAT=ierr )
      IF( ierr/=0 ) CALL errore(' simupd ', ' allocating cm ', ierr)
      ALLOCATE( vpp( cdesc%ldg ), STAT=ierr )
      IF( ierr/=0 ) CALL errore(' simupd ', ' allocating vpp ', ierr)
      ALLOCATE( bc(max_diis+1, max_diis+1), STAT=ierr )
      IF( ierr/=0 ) CALL errore(' simupd ', ' allocating bc ', ierr)
      ALLOCATE( vc(max_diis+1), STAT=ierr )
      IF( ierr/=0 ) CALL errore(' simupd ', ' allocating vc ', ierr)

      WHERE ( ABS(f) .GT. 1.0D-10 )  
        fm1 = 1.0d0 / f
      ELSEWHERE
        fm1 = 1.0d-10
      END WHERE

! ... test for convergence
      CALL converg( cgrad(:,:,1), cdesc, gemax, cnorm, .FALSE.)
      doions = cnorm .LT. tol_diis
      IF( doions ) THEN
        GOTO 999
      END IF

! ... check for an increase in energy
! ... in a good step etot < eold and therefore (etot-eold) < 0 < dethr
      teinc = .FALSE. ;  IF( ( etot - eold ) > dethr ) teinc = .TRUE.

! ... statemachine
      SELECT CASE (istate)
      CASE (2)
! ...   if the energy has increased, reject the last step and then reset the DIIS 
        istate = 2 ; IF( teinc ) istate = 1
      CASE DEFAULT
! ...   reset DIIS
        istate = 2
        CALL updis(ndiis, nowv, prevv, nsize, max_diis, 0)
      END SELECT

      IF ( istate == 1 ) THEN

        treset = .TRUE.

        WRITE( stdout,*) 'Steepest descent step; DIIS reset!'

! ...   reject last move, if it wasn't a reset
        IF( ndiis .NE. 0 ) THEN
          c0(:,:,1)    = parame(:,:,1,nowv)
          cgrad(:,:,1) = grade(:,:,1,nowv)
          var2 = -svar2
        ELSE
          var2 = svar2
        END IF

! ...   do a steepest-descent step
        CALL diis_steepest(c0, cgrad, cdesc, var2)

      ELSE

        treset = .FALSE.

! ...   perform a electronic DIIS step
        CALL updis(ndiis, nowv, prevv, nsize, max_diis, 1)

! ...   update DIIS buffers
        CALL update_diis_buffers( c0, cgrad, cdesc, nowv)

! ...   calculate the new Hessian
        CALL hesele(svar2,vpp,eigr,sfac,vps)

! ...   set up DIIS matrix
        ff=1.0d0
        bc(1:nsize,1:nsize) = 0.0d0

        DO l=gstart,ngw
          vvpp=vpp(l)*vpp(l)
          DO k=1, cdesc%nbl( 1 )
            IF(oqnr_diis ) ff = fm1(k)
            fff=2.0d0*ff*ff*vvpp
            DO i=1,nsize-1
              rri=DBLE(grade(l,k,1,i))
              rii=AIMAG(grade(l,k,1,i))
              DO j=1,i
                rrj=DBLE(grade(l,k,1,j))
                rij=AIMAG(grade(l,k,1,j))
                bc(i,j)=bc(i,j)+(rri*rrj+rii*rij)*fff
              END DO
            END DO
          END DO
        END DO

        IF(gzero) THEN
          DO k = 1, cdesc%nbl( 1 )
            IF(oqnr_diis ) ff = fm1(k)
            fff=ff*ff*vpp(1)*vpp(1)
            DO i=1,nsize-1
              r1=DBLE(grade(1,k,1,i))
              DO j=1,i
                r2=DBLE(grade(1,k,1,j))
                bc(i,j)=bc(i,j)+r1*r2*fff
              END DO
            END DO
          END DO
        END IF

        DO i=1,nsize-1
          DO j=1,i
            bc(j,i)=bc(i,j)
          END DO
        END DO
        DO i=1,nsize-1
          CALL mp_sum( bc(1:nsize,i), intra_image_comm)
        END DO

        DO i=1,nsize-1
          vc(i)=0.d0
          bc(i,nsize)=-1.d0
          bc(nsize,i)=-1.d0
        END DO
        vc(nsize)=-1.d0
        bc(nsize,nsize)=0.d0

! ...   solve system of linear equations
        CALL solve(bc,max_diis+1,nsize,vc)

! ...   compute interpolated coefficient and gradient vectors
        cm           = CMPLX(0.0d0,0.0d0)
        c0(:,:,1)    = CMPLX(0.0d0,0.0d0)
        DO i = 1, nsize-1
          DO j = 1, cdesc%nbl( 1 )
            cm(:,j)   = cm(:,j)   + vc(i) *  grade(:,j,1,i)
            c0(:,j,1) = c0(:,j,1) + vc(i) * parame(:,j,1,i)
          END DO
        END DO

! ...   estimate new parameter vectors
        ff = 1.0d0
        DO i = 1, cdesc%nbl( 1 )
          IF(oqnr_diis ) ff = fm1(i)
          c0(:,i,1) = c0(:,i,1) - ff * vpp(:) * cm(:,i)
        END DO

        eold = etot

      END IF

 999  CONTINUE

      DEALLOCATE(fm1, cm, vpp, bc, vc, STAT=ierr)
      IF( ierr/=0 ) CALL errore(' simupd ', ' deallocating arrays ', ierr)

      RETURN
      END SUBROUTINE simupd

!  ----------------------------------------------
!  ----------------------------------------------
!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE updis(ndiis,nowv,prevv,nsize,max_diis,iact)

!  this routine sets up parameters for the next DIIS step
!  ----------------------------------------------

      IMPLICIT NONE
      INTEGER ndiis,nowv,nsize,max_diis,iact,prevv

      IF(iact.EQ.1) THEN
        prevv = nowv
        nowv = mod(ndiis, max_diis) + 1  ! index of current wavefunction
        ndiis = ndiis +  1               ! DIIS steps since last reset
        nsize = ndiis +  1               ! number of states stored
        IF( nsize .GT. max_diis ) nsize = max_diis + 1
      ELSE
! ...   reset DIIS
        ndiis=0
        nowv =1
        prevv =1
        nsize=0
      END IF

      RETURN
      END SUBROUTINE updis

!=----------------------------------------------------------------------------=!

        SUBROUTINE update_diis_buffers( c, cgrad, cdesc, nowv)
          USE wave_types, ONLY: wave_descriptor
          IMPLICIT NONE
          COMPLEX(DP), INTENT(IN) :: cgrad(:,:,:)
          COMPLEX(DP), INTENT(INOUT) :: c(:,:,:)
          TYPE (wave_descriptor), INTENT(IN) :: cdesc
          INTEGER, INTENT(IN) :: nowv
          INTEGER :: ib
            DO ib = 1, cdesc%nbl( 1 )
              CALL ZDSCAL( cdesc%ngwl, -1.0d0, cgrad(1,ib,1), 1 )
            END DO
            grade(:,:,1,nowv)  = cgrad(:,:,1)
            parame(:,:,1,nowv) = c(:,:,1)
          RETURN
        END SUBROUTINE update_diis_buffers

!=----------------------------------------------------------------------------=!

        SUBROUTINE diis_steepest(c, cgrad, cdesc, var2)
          USE wave_types, ONLY: wave_descriptor
          IMPLICIT NONE
          COMPLEX(DP), INTENT(IN) :: cgrad(:,:,:)
          COMPLEX(DP), INTENT(INOUT) :: c(:,:,:)
          TYPE (wave_descriptor), INTENT(IN) :: cdesc
          REAL(DP), INTENT(IN) :: var2
          c(:,:,1) = c(:,:,1) + var2 * cgrad(:,:,1)
          IF ( cdesc%gzero ) THEN
             c(1,:,1) = CMPLX( DBLE( c(1,:,1) ), 0.d0 )
          END IF
          RETURN
        END SUBROUTINE diis_steepest

!=----------------------------------------------------------------------------=!

   SUBROUTINE solve(b, ldb, ndim, v)

      USE mp_global, ONLY: root, intra_image_comm
      USE mp, ONLY: mp_bcast
      !
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ldb, ndim
      REAL(DP) :: b(:,:), v(:)

! ... declare other variables

      INTEGER  :: i, j, k, info
      REAL(DP) :: ap(ndim*(ndim+1)/2)
      INTEGER  :: ipiv(ndim)


      k = 0
      DO j = 1,ndim
        DO i = j,ndim
          k = k + 1
          ap(k) = b(i,j)
        END DO
      END DO

      CALL DSPTRF ( 'L', ndim, ap, ipiv, info )
      IF(info .NE. 0) THEN
        CALL errore(' solve ',' dsptrf has failed ', info)
      END IF
      CALL DSPTRS ( 'L', ndim, 1, ap, ipiv, v, ndim, info )
      IF(info .NE. 0) THEN
        CALL errore(' solve ',' dsptrs has failed ', info)
      END IF

      CALL mp_bcast(v, root, intra_image_comm)

      RETURN
      END SUBROUTINE solve



!=----------------------------------------------------------------------------=!
END MODULE diis
!=----------------------------------------------------------------------------=!

