!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file 'License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
PROGRAM rrkj2upf
  !---------------------------------------------------------------------
  !
  !     Convert a pseudopotential written in "rrkj3" format
  !     (Rabe-Rappe-Kaxiras-Joannopoulos with 3 Bessel functions)
  !     to unified pseudopotential format
  !
  IMPLICIT NONE
  CHARACTER(len=256) filein, fileout
  !
  !
  IF ( trim(filein) == ' ') &
       CALL errore ('rrkj2upf', 'usage: rrkj2upf "file-to-be-converted"', 1)
  CALL get_file ( filein )
  OPEN (unit = 1, file = filein, status = 'old', form = 'formatted')
  CALL read_rrkj(1)
  CLOSE (1)

  ! convert variables read from rrkj3 format into those needed
  ! by the upf format - add missing quantities

  CALL convert_rrkj

  fileout=trim(filein)//'.UPF'
  PRINT '(''Output PP file in UPF format :  '',a)', fileout

  OPEN(unit=2,file=fileout,status='unknown',form='formatted')
  CALL write_upf_v1(2)
  CLOSE (unit=2)

STOP
20 WRITE (6,'("rrkj2upf: error reading pseudopotential file name")')
   STOP

END PROGRAM rrkj2upf

MODULE rrkj3
  !
  ! All variables read from RRKJ3 file format
  !
  ! trailing underscore means that a variable with the same name
  ! is used in module 'upf' containing variables to be written
  !
  CHARACTER(len=75):: titleps
  CHARACTER (len=2), ALLOCATABLE :: els_(:)
  INTEGER :: pseudotype_, iexch_, icorr_, igcx_, igcc_, mesh_, &
       nwfs_, nbeta_, lmax_
  LOGICAL :: rel_, nlcc_
  real (8) :: zp_, etotps_, xmin, rmax, zmesh, dx, rcloc_
  INTEGER, ALLOCATABLE:: lchi_(:), nns_(:), ikk2_(:)
  real (8), ALLOCATABLE :: rcut_(:), rcutus_(:), oc_(:), &
       beta(:,:), dion_(:,:), qqq_(:,:), ddd(:,:), qfunc_(:,:,:), &
       rho_atc_(:), rho_at_(:), chi_(:,:), vloc_(:)
END MODULE rrkj3
!
!     ----------------------------------------------------------
SUBROUTINE read_rrkj(iunps)
  !     ----------------------------------------------------------
  !
  USE rrkj3
  IMPLICIT NONE
  INTEGER :: iunps
  INTEGER :: nb, mb, n, ir, ios

  !---  > Start the header reading
  READ (iunps, '(a75)', err = 100) titleps
  READ (iunps, *, err = 100)  pseudotype_
  READ (iunps, *, err = 100) rel_, nlcc_
  READ (iunps, *, err=100) iexch_, icorr_, igcx_, igcc_
  READ (iunps, '(2e17.11,i5)') zp_, etotps_, lmax_
  READ (iunps, '(4e17.11,i5)', err=100) xmin, rmax, zmesh, dx, mesh_
  READ (iunps, *, err=100) nwfs_, nbeta_

  ALLOCATE(rcut_(nwfs_), rcutus_(nwfs_))
  READ (iunps, *, err=100) (rcut_(nb), nb=1,nwfs_)
  READ (iunps, *, err=100) (rcutus_(nb), nb=1,nwfs_)

  ALLOCATE(els_(nwfs_), nns_(nwfs_), lchi_(nwfs_), oc_(nwfs_))
  DO nb = 1, nwfs_
     READ (iunps, '(a2,2i3,f6.2)', err = 100) els_(nb), &
          nns_(nb), lchi_(nb) , oc_(nb)
  ENDDO

  ALLOCATE(ikk2_(nbeta_))
  ALLOCATE(beta( mesh_,nbeta_))
  ALLOCATE(dion_(nbeta_,nbeta_))
  ALLOCATE(ddd (nbeta_,nbeta_))
  ALLOCATE(qqq_(nbeta_,nbeta_))
  ALLOCATE(qfunc_(mesh_,nbeta_,nbeta_))

  DO nb = 1, nbeta_
     READ (iunps, *, err = 100) ikk2_(nb)
     READ (iunps, *, err = 100) (beta (ir, nb) , ir = 1,ikk2_(nb) )
     DO ir = ikk2_(nb) + 1, mesh_
        beta (ir, nb) = 0.d0
     ENDDO
     DO mb = 1, nb
        READ (iunps, *, err = 100) dion_(nb, mb)
        dion_(mb, nb) = dion_(nb, mb)
        IF (pseudotype_==3) THEN
           READ (iunps, *, err = 100) qqq_(nb, mb)
           qqq_(mb, nb) = qqq_(nb, mb)
           READ (iunps, *, err = 100) (qfunc_(n,nb, mb), n = 1, mesh_)
           DO n = 1, mesh_
              qfunc_(n, mb, nb) = qfunc_(n, nb, mb)
           ENDDO
        ELSE
           qqq_(nb, mb) = 0.d0
           qqq_(mb, nb) = 0.d0
           DO n = 1, mesh_
              qfunc_(n, nb, mb) = 0.d0
              qfunc_(n, mb, nb) = 0.d0
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  !
  !     read the local potential
  !
  ALLOCATE(vloc_(mesh_))
  READ (iunps, *, err = 100) rcloc_, (vloc_(ir ) , ir = 1, mesh_ )
  !
  !     read the atomic charge
  !
  ALLOCATE(rho_at_(mesh_))
  READ (iunps, *, err=100) (rho_at_(ir), ir=1,mesh_)
  !
  !     if present read the core charge
  !
  ALLOCATE(rho_atc_(mesh_))
  IF (nlcc_) THEN
     READ (iunps, *, err=100) (rho_atc_(ir), ir=1, mesh_)
  ENDIF
  !
  !     read the pseudo wavefunctions of the atom
  !
  ALLOCATE(chi_(mesh_,nwfs_))
  READ (iunps, *, err=100) ( (chi_(ir,nb), ir = 1,mesh_) , nb = 1, nwfs_)
  !
  !     ----------------------------------------------------------
  WRITE (6,'(a)') 'Pseudopotential successfully read'
  !     ----------------------------------------------------------
  !
  RETURN
100 WRITE (6,'("read_rrkj: error reading pseudopotential file")')
    STOP

END SUBROUTINE read_rrkj

SUBROUTINE convert_rrkj
  !     ----------------------------------------------------------
  !
  USE rrkj3
  USE upf
  USE constants, ONLY : fpi
  IMPLICIT NONE
  INTEGER i, n
  real(8) :: x


  WRITE(generated, '("Generated using Andrea Dal Corso code (rrkj3)")')
  WRITE(date_author,'("Author: Andrea Dal Corso   Generation date: unknown")')
  comment = 'Info:'//titleps
  IF (rel_) THEN
     rel = 1
  ELSE
     rel = 0
  ENDIF
  rcloc = rcloc_
  nwfs = nwfs_
  ALLOCATE( els(nwfs), oc(nwfs), epseu(nwfs))
  ALLOCATE(lchi(nwfs), nns(nwfs) )
  ALLOCATE(rcut (nwfs), rcutus (nwfs))
  DO i=1, nwfs
     nns (i)  = nns_(i)
     lchi(i)  = lchi_(i)
     rcut(i)  = rcut_(i)
     rcutus(i)= rcutus_(i)
     oc (i)   = oc_(i)
     els(i)   = els_(i)
     epseu(i) = 0.0d0
  ENDDO
  DEALLOCATE (els_, oc_, rcutus_, rcut_, nns_)

  psd  = titleps (7:8)
  IF (pseudotype_==3) THEN
     pseudotype = 'US'
  ELSE
     pseudotype = 'NC'
  ENDIF
  nlcc = nlcc_
  zp = zp_
  etotps = etotps_
  ecutrho=0.0d0
  ecutwfc=0.0d0
  lmax = lmax_
  mesh = mesh_
  nbeta = nbeta_
  ntwfc = 0
  DO i=1, nwfs
     IF (oc(i) > 1.0d-12) ntwfc = ntwfc + 1
  ENDDO
  ALLOCATE( elsw(ntwfc), ocw(ntwfc), lchiw(ntwfc) )
  n = 0
  DO i=1, nwfs
     IF (oc(i) > 1.0d-12) THEN
        n = n + 1
        elsw(n) = els(i)
        ocw (n) = oc (i)
        lchiw(n)=lchi(i)
     ENDIF
  ENDDO
  iexch = iexch_
  icorr = icorr_
  igcx  = igcx_
  igcc  = igcc_

  ALLOCATE(rab(mesh))
  ALLOCATE(  r(mesh))
  ! define logarithmic mesh
  DO i = 1, mesh
     x = xmin + dble(i-1) * dx
     r  (i) = exp(x) / zmesh
     rab(i) = dx * r(i)
  ENDDO

  ALLOCATE (rho_atc(mesh))
  ! rrkj rho_core(r) =  4pi*r^2*rho_core(r) UPF
  rho_atc (:) = rho_atc_(:) / fpi / r(:)**2
  DEALLOCATE (rho_atc_)

  ALLOCATE (vloc0(mesh))
  vloc0 = vloc_
  DEALLOCATE (vloc_)

  ALLOCATE(ikk2(nbeta), lll(nbeta))
  ikk2 = ikk2_
  lll  = lchi_
  DEALLOCATE (ikk2_, lchi_)
!  kkbeta  = 0
!  do nb=1,nbeta
!     kkbeta  = max (kkbeta , ikk2(nb) )
!  end do
  ALLOCATE(betar(mesh,nbeta))
  betar = 0.0d0
  DO i=1, nbeta
     betar(1:ikk2(i),i) = beta(1:ikk2(i),i)
  ENDDO
  DEALLOCATE (beta)

  ALLOCATE(dion(nbeta,nbeta))
  dion = dion_
  DEALLOCATE (dion_)

  ALLOCATE(qqq(nbeta,nbeta))
  qqq = qqq_
  DEALLOCATE (qqq_)

  ALLOCATE(qfunc(mesh,nbeta,nbeta))
  qfunc = qfunc_

  nqf = 0
  nqlc= 0

  ALLOCATE (rho_at(mesh))
  rho_at = rho_at_
  DEALLOCATE (rho_at_)

  ALLOCATE (chi(mesh,ntwfc))
  n = 0
  DO i=1, nwfs
     IF (oc(i) > 1.0d-12) THEN
        n = n + 1
        chi(:,n) = chi_(:,i)
     ENDIF
  ENDDO
  DEALLOCATE (chi_)

  RETURN
END SUBROUTINE convert_rrkj

