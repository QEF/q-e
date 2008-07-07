!
! Copyright (C) 2008 Simon Binnie
! This file is distributed under the terms of the
! GNU General Public License. See the file 'License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#include "f_defs.h"
!---------------------------------------------------------------------
PROGRAM casino2upf
  !---------------------------------------------------------------------
  !
  !     Convert a pseudopotential written in CASINO tabulated
  !     format to unified pseudopotential format

  IMPLICIT NONE
  CHARACTER(len=256) filein, fileout
  CHARACTER(len=256), ALLOCATABLE:: wavefile(:)
  INTEGER nofiles, i

  PRINT*, 'CASINO2UPF Convertor'
  PRINT*, 'Enter CASINO pp.data filename:'

  CALL get_file ( filein )

  PRINT*, 'How many wavefunction *files* are you using?'
  READ(*,*) nofiles
  ALLOCATE(wavefile(nofiles))
  PRINT*, 'Enter wavefunction files, starting with the ground state:'
  DO i=1,nofiles
     CALL get_file ( wavefile(i) )
     OPEN(unit=i,file=wavefile(i),status='old',form='formatted')  
  ENDDO
  OPEN(unit=99,file=filein,status='old',form='formatted')  

  CALL read_casino(99,nofiles)
  CLOSE (unit=99)
  DO i=1,nofiles
     CLOSE (i)
  ENDDO

  ! convert variables read from CASINO format into those needed
  ! by the upf format - add missing quantities

  CALL convert_casino

  fileout=TRIM(filein)//'.UPF'
  PRINT '(''Output PP file in US format :  '',a)', fileout

  OPEN(unit=2,file=fileout,status='unknown',form='formatted')
  CALL write_upf(2)
  CLOSE (unit=2)
  STOP
20 CALL errore ('casino2upf', 'Reading pseudo file name ', 1)

END PROGRAM casino2upf

MODULE casino
  
  !
  ! All variables read from CASINO file format
  ! 
  ! trailing underscore means that a variable with the same name
  ! is used in module 'upf' containing variables to be written
  !


  CHARACTER(len=20) :: dft_
  CHARACTER(len=2)  :: psd_
  REAL(8) :: zp_
  INTEGER nlc, nnl, lmax_, lloc, nchi, rel_
  LOGICAL :: numeric, bhstype, nlcc_
  REAL(8) :: alpc(2), cc(2), alps(3,0:3), aps(6,0:3)
  REAL(8) :: a_nlcc, b_nlcc, alpha_nlcc

  REAL(8) :: zmesh, xmin, dx
  REAL(8), ALLOCATABLE::  r_(:), rab_(:)
  INTEGER :: mesh_

  REAL(8), ALLOCATABLE::  vnl(:,:), rho_atc_(:), rho_at_(:)
  INTEGER, ALLOCATABLE:: lchi_(:), nns_(:)
  REAL(8), ALLOCATABLE:: chi_(:,:),  oc_(:)

END MODULE casino
! 
!     ----------------------------------------------------------
SUBROUTINE read_casino(iunps,nofiles)
  !     ----------------------------------------------------------
  ! 
  USE casino
  USE upf , ONLY : els
  USE kinds
  IMPLICIT NONE
  TYPE :: wavfun_list
     INTEGER :: occ,eup,edwn, nquant, lquant 
     CHARACTER(len=2) :: label
     REAL*8, ALLOCATABLE :: wavefunc(:)
     TYPE (wavfun_list), POINTER :: p

  END TYPE wavfun_list

  TYPE (wavfun_list), POINTER :: mhead
  TYPE (wavfun_list), POINTER :: mptr
  TYPE (wavfun_list), POINTER :: mtail

  INTEGER :: iunps, nofiles
  !
  LOGICAL :: groundstate, found
  CHARACTER(len=1), DIMENSION(0:3) :: convel=(/'s','p','d','f'/)
  CHARACTER(len=2) :: label, rellab
  REAL (8) :: erf
  REAL(DP), parameter :: r_exp=20._dp/1500._dp


  INTEGER :: l, i, ir, nb, gsorbs, j,k,m,tmp, lquant, orbs, nquant

  INTEGER, ALLOCATABLE :: gs(:,:)
  EXTERNAL erf

  NULLIFY (  mhead, mptr, mtail )
  dft_ = 'HF'   !Hardcoded at the moment should eventually be HF anyway

  nlc = 0              !These two values are always 0 for numeric pps
  nnl = 0              !so lets just hard code them

  nlcc_ = .FALSE.       !Again these two are alwas false for CASINO pps
  bhstype = .FALSE.    


 
  READ(iunps,'(a2,35x,a2)') rellab, psd_  
  READ(iunps,*)
  IF ( rellab .EQ. 'DF' ) THEN
     rel_=1
  ELSE
     rel_=0
  ENDIF
 
  READ(iunps,*) zmesh,zp_  !Here we are reading zmesh (atomic #) and
  DO i=1,3                 !zp_ (pseudo charge)
     READ(iunps,*)
  ENDDO
  READ(iunps,*) lmax_               !reading in lmax
  IF ( zp_.LE.0d0 ) &
       CALL errore( 'read_casino','Wrong zp ',1 )
  IF ( lmax_.GT.3.OR.lmax_.LT.0 ) &
       CALL errore( 'read_casino','Wrong lmax ',1 )

  lloc=lmax_ !think lloc shoudl always = lmax for this case yes/no ??
 
  !
  !    compute the radial mesh
  !

  DO i=1,3
     READ(iunps,*)
  ENDDO
  READ(iunps,*) mesh_   !Reading in total no. of mesh points


  ALLOCATE(  r_(mesh_))
  ALLOCATE(rab_(mesh_))
  READ(iunps,*)
  DO i=1,mesh_
     READ(iunps,*) r_(i)
  ENDDO
  DO ir = 1, mesh_
     rab_(ir) = r_exp  * r_(ir) !hardcoded at the moment
  END DO


  ALLOCATE(vnl(mesh_,0:lmax_))

  DO l = 0, lmax_
     READ(iunps, '(a)', err=300)
     READ(iunps, *, err=300)  (vnl(ir,l),ir=1,mesh_)
  ENDDO

  DO l = 0, lmax_
     DO ir = 1, mesh_
        vnl(ir,l) = vnl(ir,l)/r_(ir) !Removing the factor of r CASINO has
     ENDDO
     vnl(1,l) = 0                   !correcting for the divide by zero 
  ENDDO

  ALLOCATE(rho_atc_(mesh_))
  IF(nlcc_) THEN
     READ(iunps, *, err=300) ( rho_atc_(ir), ir=1,mesh_ )
  ENDIF

  !
  ! subtract the local part
  !

  DO l = 0, lmax_
     IF ( l.NE.lloc ) vnl(:,l) = vnl(:,l) - vnl(:,lloc)
  ENDDO


  ALLOCATE(mhead)
  mtail => mhead

  mptr => mhead

  NULLIFY(mtail%p)
  groundstate=.TRUE.
  DO j=1,nofiles

     DO i=1,4

        READ(j,*)
     ENDDO

     READ(j,*) orbs

     IF ( groundstate ) THEN 

        ALLOCATE( gs(orbs,3) )

        gs = 0
        gsorbs = orbs
     END IF

     DO i=1,2
        READ(j,*)
     ENDDO

     READ(j,*) mtail%eup, mtail%edwn
     READ(j,*)

     DO i=1,mtail%eup+mtail%edwn
        READ(j,*) tmp, nquant, lquant

        IF ( groundstate ) THEN
           found = .TRUE.

           DO m=1,orbs

              IF ( (nquant==gs(m,1) .AND. lquant==gs(m,2)) ) THEN
                 gs(m,3) = gs(m,3) + 1
                 EXIT
              END IF

              found = .FALSE.   

           ENDDO

           IF (.NOT. found ) THEN

              DO m=1,orbs

                 IF ( gs(m,1) == 0 ) THEN
                    gs(m,1) = nquant
                    gs(m,2) = lquant
                    gs(m,3) = 1

                    EXIT
                 ENDIF

              ENDDO

           ENDIF

        ENDIF

     ENDDO

     READ(j,*)
     READ(j,*) 

     DO i=1,mesh_ 
        READ(j,*)
     ENDDO

     DO k=1,orbs
        READ(j,'(13x,a2)', err=300) label
        READ(j,*), tmp, nquant, lquant

        IF ( .NOT. groundstate ) THEN
           found = .FALSE.

           DO m = 1,gsorbs

              IF ( nquant == gs(m,1) .AND. lquant == gs(m,2) ) THEN
                 found = .TRUE.
                 EXIT
              END IF
           END DO
           mptr => mhead
           DO
              IF ( .NOT. ASSOCIATED(mptr) )EXIT
              IF ( nquant == mptr%nquant .AND. lquant == mptr%lquant ) found = .TRUE.        
              mptr =>mptr%p
           END DO
           IF ( found ) THEN
              DO i=1,mesh_
                 READ(j,*)
              ENDDO

              CYCLE
           END IF
        END IF
        IF ( ALLOCATED(mtail%wavefunc) ) THEN
           ALLOCATE(mtail%p)
           mtail=>mtail%p
           NULLIFY(mtail%p)
           ALLOCATE( mtail%wavefunc(mesh_) )
        ELSE
           ALLOCATE( mtail%wavefunc(mesh_) )
        END IF
        mtail%label = label
        mtail%nquant = nquant
        mtail%lquant = lquant


        READ(j, *, err=300) (mtail%wavefunc(ir),ir=1,mesh_)
     ENDDO
     groundstate = .FALSE.
  ENDDO

  nchi =0
  mptr => mhead
  DO
     IF ( .NOT. ASSOCIATED(mptr) )EXIT
     nchi=nchi+1

     mptr =>mptr%p
  END DO

  ALLOCATE(lchi_(nchi), els(nchi), nns_(nchi))
  ALLOCATE(oc_(nchi))
  ALLOCATE(chi_(mesh_,nchi))
  oc_ = 0

  !Sort out the occupation numbers
  DO i=1,gsorbs  
     oc_(i)=gs(i,3)
  ENDDO
  deallocate( gs )

  i=1
  mptr => mhead
  DO
     IF ( .NOT. ASSOCIATED(mptr) )EXIT
     nns_(i) = mptr%nquant 
     lchi_(i) = mptr%lquant
     els(i) = mptr%label

     DO ir=1,mesh_

        chi_(ir:,i) = mptr%wavefunc(ir)
     ENDDO
     deallocate( mptr%wavefunc )
     mptr =>mptr%p
     i=i+1
  END DO

  !Clean up the linked list (deallocate it)
  DO
     IF ( .NOT. ASSOCIATED(mhead) )EXIT
     mptr => mhead
     mhead => mhead%p
     deallocate( mptr )
  END DO


  !
  !    compute the atomic charges
  !
  ALLOCATE(rho_at_(mesh_))
  rho_at_(:)=0.d0
  DO nb = 1, nchi
     IF( oc_(nb).NE.0.d0) &
          &  rho_at_(:) = rho_at_(:) + oc_(nb)*chi_(:,nb)**2
  END DO
  !     ----------------------------------------------------------
  WRITE (6,'(a)') 'Pseudopotential successfully read'
  !     ----------------------------------------------------------
  RETURN

300 CALL errore('read_casino','pseudo file is empty or wrong',1)

END SUBROUTINE read_casino

  !     ----------------------------------------------------------
SUBROUTINE convert_casino
  !     ----------------------------------------------------------
  USE casino
  USE upf
  USE funct, ONLY : set_dft_from_name, get_iexch, get_icorr, get_igcx, get_igcc
  IMPLICIT NONE
  REAL(8), PARAMETER :: rmax = 10.0d0
  REAL(8), ALLOCATABLE :: aux(:)
  REAL(8) :: vll
  INTEGER :: kkbeta, l, iv, ir, i

  WRITE(generated, '("From a Trail & Needs tabulated PP for CASINO")')
  WRITE(date_author,'("Author: unknown    Generation date: as well")')
  comment = 'Info: automatically converted from CASINO Tabulated format'

  rel = rel_


  rcloc = 0.0d0
  nwfs  = nchi 
  ALLOCATE( oc(nwfs), epseu(nwfs))
  ALLOCATE(lchi(nwfs), nns(nwfs) )
  ALLOCATE(rcut (nwfs), rcutus (nwfs))
  DO i=1, nwfs
     nns (i)  = nns_(i)
     lchi(i)  = lchi_(i)
     rcut(i)  = 0.0d0
     rcutus(i)= 0.0d0
     oc (i)   = oc_(i)
     epseu(i) = 0.0d0
  END DO
  DEALLOCATE (lchi_, oc_, nns_)

  psd = psd_
  pseudotype = 'NC'
  nlcc = nlcc_
  zp = zp_
  etotps = 0.0d0
  ecutrho=0.0d0
  ecutwfc=0.0d0
  IF ( lmax_ == lloc) THEN
     lmax = lmax_-1
  ELSE
     lmax = lmax_
  END IF
  nbeta= lmax_ 
  mesh = mesh_
  ntwfc= nchi
  ALLOCATE( elsw(ntwfc), ocw(ntwfc), lchiw(ntwfc) )
  DO i=1, nchi
     lchiw(i) = lchi(i)
     ocw(i)   = oc(i)
     elsw(i)  = els(i)
  END DO
  CALL set_dft_from_name(dft_)
  iexch = get_iexch()
  icorr = get_icorr()
  igcx  = get_igcx()
  igcc  = get_igcc()

  ALLOCATE(rab(mesh))
  ALLOCATE(  r(mesh))
  rab = rab_
  r = r_

  ALLOCATE (rho_atc(mesh))
  rho_atc = rho_atc_
  DEALLOCATE (rho_atc_)

  ALLOCATE (vloc0(mesh))
  vloc0(:) = vnl(:,lloc)

  IF (nbeta > 0) THEN

     ALLOCATE(ikk2(nbeta), lll(nbeta))
     kkbeta=mesh
     DO ir = 1,mesh
        IF ( r(ir) > rmax ) THEN
           kkbeta=ir
           EXIT
        END IF
     END DO

     ! make sure kkbeta is odd as required for simpson
     IF(MOD(kkbeta,2) == 0) kkbeta=kkbeta-1 
     ikk2(:) = kkbeta
     ALLOCATE(aux(kkbeta))
     ALLOCATE(betar(mesh,nbeta))
     ALLOCATE(qfunc(mesh,nbeta,nbeta))
     ALLOCATE(dion(nbeta,nbeta))
     ALLOCATE(qqq (nbeta,nbeta))
     qfunc(:,:,:)=0.0d0
     dion(:,:) =0.d0
     qqq(:,:)  =0.d0
     iv=0
     DO i=1,nchi
        l=lchi(i)
        IF (l.NE.lloc) THEN
           iv=iv+1
           lll(iv)=l
           DO ir=1,kkbeta
              betar(ir,iv)=chi_(ir,i)*vnl(ir,l)
              aux(ir) = chi_(ir,i)**2*vnl(ir,l)

           END DO
           CALL simpson(kkbeta,aux,rab,vll)
           dion(iv,iv) = 1.0d0/vll
        END IF
        IF(iv >= nbeta) EXIT  ! skip additional pseudo wfns
     ENDDO


     DEALLOCATE (vnl, aux)

     !
     !   redetermine ikk2
     !
     DO iv=1,nbeta
        ikk2(iv)=kkbeta
        DO ir = kkbeta,1,-1
           IF ( ABS(betar(ir,iv)) > 1.d-12 ) THEN
              ikk2(iv)=ir
              EXIT
           END IF
        END DO
     END DO
  END IF
  ALLOCATE (rho_at(mesh))
  rho_at = rho_at_
  DEALLOCATE (rho_at_)

  ALLOCATE (chi(mesh,ntwfc))
  chi = chi_
  DEALLOCATE (chi_)

  RETURN
END SUBROUTINE convert_casino
