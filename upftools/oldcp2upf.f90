!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file 'License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
PROGRAM oldcp2upf
  !---------------------------------------------------------------------
  !
  !     Convert a pseudopotential written in the old CP90 format
  !     (without core correction) to unified pseudopotential format
  !
  IMPLICIT NONE
  CHARACTER(len=256) filein, fileout
  !
  !
  CALL get_file ( filein )
  OPEN (unit = 1, file = filein, status = 'old', form = 'formatted')
  CALL read_oldcp(1)
  CLOSE (1)

  ! convert variables read from old CP90 format into those needed
  ! by the upf format - add missing quantities

  CALL convert_oldcp

  fileout=trim(filein)//'.UPF'
  PRINT '("Output PP file in UPF format :  ",a)', fileout

  OPEN(unit=2,file=fileout,status='unknown',form='formatted')
  CALL write_upf_v1(2)
  CLOSE (unit=2)

STOP
20 CALL errore ('oldcp2upf', 'Reading pseudo file name ', 1)

END PROGRAM oldcp2upf

MODULE oldcp
  !
  ! All variables read from old CP90 file format
  !
  real(8) :: amesh, z, zv
  INTEGER :: exfact, lloc, nbeta_, mesh_
  real(8) :: wrc1, rc1, wrc2, rc2, rcl(3,3), al(3,3), bl(3,3)
  real(8), ALLOCATABLE :: r_(:), vnl(:,:), chi_(:,:)
  !
  !------------------------------

END MODULE oldcp
!
!     ----------------------------------------------------------
SUBROUTINE read_oldcp(iunps)
  !     ----------------------------------------------------------
  !
  USE oldcp
  IMPLICIT NONE
  INTEGER :: iunps
  !
  real(8), EXTERNAL :: qe_erf
  INTEGER :: i, l, j, jj
  !
  READ(iunps,*, end=10, err=10) z, zv, nbeta_, lloc, exfact
  IF (z < 1 .or. z > 100 .or. zv < 1 .or. zv > 25 ) &
       CALL errore ('read_oldcp','wrong potential read',1)
  READ(iunps,*, end=10, err=10) wrc1, rc1, wrc2, rc2
  READ(iunps,*, end=10, err=10) ( ( rcl(i,l), al(i,l), &
                     bl(i,l), i = 1, 3), l = 1, 3)
  READ(iunps,*, end=10, err=10) mesh_, amesh
  ALLOCATE(r_(mesh_))
  ALLOCATE (chi_(mesh_,nbeta_))
  DO l = 1, nbeta_
     IF (l > 1) READ(iunps,*, end=10, err=10) mesh_, amesh
     DO j = 1, mesh_
        READ(iunps,*, end=10, err=10) jj, r_(j), chi_(j,l)
     ENDDO
  ENDDO
  !
  !  convert analytic to numeric form
  !
  ALLOCATE (vnl(mesh_,0:nbeta_))
  DO l=0,nbeta_
     !
     !  DO NOT USE f90 ARRAY SYNTAX: qe_erf IS NOT AN INTRINSIC FUNCTION!!!
     !
     DO j=1, mesh_
        vnl(j,l)= - (wrc1*qe_erf(sqrt(rc1)*r_(j)) + &
                     wrc2*qe_erf(sqrt(rc2)*r_(j)) ) * zv/r_(j)
     ENDDO
     !
     DO i=1,3
        vnl(:,l)= vnl(:,l)+ (al(i,l+1)+ bl(i,l+1)*r_(:)**2) * &
             exp(-rcl(i,l+1)*r_(:)**2)
     ENDDO
  ENDDO

  RETURN
10 CALL errore('read_oldcp','error in reading file',1)

END SUBROUTINE read_oldcp

!     ----------------------------------------------------------
SUBROUTINE convert_oldcp
  !     ----------------------------------------------------------
  !
  USE oldcp
  USE upf
  IMPLICIT NONE
  real(8), PARAMETER :: rmax = 10.0d0
  real(8), ALLOCATABLE :: aux(:)
  real(8) :: vll
  CHARACTER (len=20):: dft
  CHARACTER (len=2), EXTERNAL :: atom_name
  INTEGER :: kkbeta
  INTEGER :: l, i, ir, iv
  !
  WRITE(generated, '("Generated using unknown code")')
  WRITE(date_author,'("Author: unknown    Generation date: as well")')
  comment = 'Info: automatically converted from old CP90 format'
  ! reasonable assumption
  IF (z > 18) THEN
     rel = 1
  ELSE
     rel = 0
  ENDIF
  rcloc = 0.0d0
  nwfs  = nbeta_
  ALLOCATE( els(nwfs), oc(nwfs), epseu(nwfs))
  ALLOCATE(lchi(nwfs), nns(nwfs) )
  ALLOCATE(rcut (nwfs), rcutus (nwfs))
  DO i=1, nwfs
     WRITE(*,'("Wavefunction # ",i1,": label, occupancy > ")', advance="NO") i
     READ (5,*) els(i), oc(i)
     nns (i)  = 0
     lchi(i)  = i-1
     rcut(i)  = 0.0d0
     rcutus(i)= 0.0d0
     epseu(i) = 0.0d0
  ENDDO
  psd   = atom_name (nint(z))
  pseudotype = 'NC'
  nlcc = .false.
  zp = nint(zv)
  etotps =0.0d0
  ecutrho=0.0d0
  ecutwfc=0.0d0
  lmax  = nbeta_ - 1
  nbeta = nbeta_
  mesh  = mesh_
  ntwfc = nwfs
  ALLOCATE( elsw(ntwfc), ocw(ntwfc), lchiw(ntwfc) )
  DO i=1, nwfs
     lchiw(i) = lchi(i)
     ocw(i)   = oc(i)
     elsw(i)  = els(i)
  ENDDO
  !
  IF ( exfact==0) THEN
     iexch=1; icorr=1; igcx=0; igcc=0 ! Perdew-Zunger
  ELSEIF ( exfact==1) THEN
     iexch=1; icorr=3; igcx=1; igcc=3 ! Becke-Lee-Yang-Parr
  ELSEIF ( exfact==2) THEN
     iexch=1; icorr=1; igcx=1; igcc=0 ! Becke88 exchange
  ELSEIF (exfact==-5.or.exfact==3) THEN
     iexch=1; icorr=1; igcx=1; igcc=1 ! Becke88-Perdew 86
  ELSEIF (exfact==-6.or.exfact==4) THEN
     iexch=1; icorr=4; igcx=2; igcc=2 ! Perdew-Wang 91
  ELSEIF (exfact== 5) THEN
     iexch=1; icorr=4; igcx=3; igcc=4 ! Perdew-Becke-Erkerhof
  ELSE
     CALL errore('convert','Wrong xc in pseudopotential',1)
  ENDIF

  ALLOCATE(rab(mesh))
  ALLOCATE(  r(mesh))
  r = r_
  rab = r * log( amesh )
  !
  !  convert analytic to numeric form
  !
  !
  ALLOCATE (vloc0(mesh))
  ! the factor 2 converts from Hartree to Rydberg
  vloc0(:) = vnl(:,lloc)*2.d0

  IF (nbeta > 0) THEN

     ALLOCATE(ikk2(nbeta), lll(nbeta))
     kkbeta=mesh
     DO ir = 1,mesh
        IF ( r(ir) > rmax ) THEN
           kkbeta=ir
           exit
        ENDIF
     ENDDO
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
     DO i=1,nwfs
        l=lchi(i)
        IF (l/=lloc) THEN
           iv=iv+1
           lll(iv)=l
           DO ir=1,kkbeta
              ! the factor 2 converts from Hartree to Rydberg
              betar(ir,iv) = 2.d0 * chi_(ir,l+1) * &
                   ( vnl(ir,l) - vnl(ir,lloc) )
              aux(ir) = chi_(ir,l+1) * betar(ir,iv)
           ENDDO
           CALL simpson(kkbeta,aux,rab,vll)
           dion(iv,iv) = 1.0d0/vll
        ENDIF
     ENDDO

  ENDIF

  ALLOCATE (rho_at(mesh))
  rho_at = 0.d0
  DO i=1,nwfs
     rho_at(:) = rho_at(:) + ocw(i) * chi_(:,i) ** 2
  ENDDO

  ALLOCATE (chi(mesh,ntwfc))
  chi = chi_

  !     ----------------------------------------------------------
  WRITE (6,'(a)') 'Pseudopotential successfully converted'
  !     ----------------------------------------------------------
  RETURN
END SUBROUTINE convert_oldcp
